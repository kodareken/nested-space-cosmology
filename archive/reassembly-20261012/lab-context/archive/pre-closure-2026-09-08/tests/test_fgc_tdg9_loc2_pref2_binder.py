from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
from types import SimpleNamespace
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from recursive_horizons.fgc.evolution import tdg9_loc2_pref2_binder as binder


ROOT = Path(__file__).resolve().parents[1]
CONFIG = (ROOT / binder.CONFIG_PATH).read_bytes()
RESULT = (ROOT / binder.RESULT_PATH).read_bytes()


class TDG9LOC2PREF2BinderTests(unittest.TestCase):
    def test_compact_result_is_raw_and_store_blind(self) -> None:
        forbidden = AssertionError("compact verification opened a live input")
        with (
            patch.object(binder, "_raw_snapshot", side_effect=forbidden),
            patch.object(binder, "_snapshot_store", side_effect=forbidden),
            patch.object(binder, "bind_raw_result", side_effect=forbidden),
        ):
            value = binder.validate_compact_result(CONFIG, RESULT)
        replay = value["artifact_payload"]["independent_replay"]
        self.assertEqual(replay["occurrence_count"], 10)
        self.assertEqual(replay["route_agreement_count"], 60)
        self.assertEqual(replay["both_components_independently_fail_count"], 10)

    def test_compact_mutations_fail_closed(self) -> None:
        value = json.loads(RESULT)
        mutations = (
            lambda item: item["artifact_payload"]["claims"].__setitem__(
                "temporal_instrument_remedy_selected", True
            ),
            lambda item: item["artifact_payload"]["independent_replay"].__setitem__(
                "route_agreement_count", 59
            ),
            lambda item: item["artifact_payload"]["independent_replay"][
                "occurrences"
            ][0].__setitem__("raw_occurrence_sha256", "0" * 64),
            lambda item: item["artifact_payload"]["conclusion"].__setitem__(
                "physics_inference_permitted", True
            ),
        )
        for mutation in mutations:
            changed = deepcopy(value)
            mutation(changed)
            with self.assertRaises(binder.TDG9LOC2PREF2Error):
                binder.validate_compact_result(
                    CONFIG,
                    binder.canonical_result(changed),
                )

    def test_json_rejects_duplicate_keys_and_noncanonical_bytes(self) -> None:
        with self.assertRaises(binder.TDG9LOC2PREF2Error):
            binder._json(b'{"a":1,"a":2}\n', "duplicate")
        with self.assertRaises(binder.TDG9LOC2PREF2Error):
            binder._json(b'{"a": 1}\n', "noncanonical")

    def test_raw_reader_rejects_symlink_and_extra_leaf(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / binder.RAW_NAMESPACE
            target.parent.mkdir(parents=True)
            real = root / "real"
            real.mkdir()
            target.symlink_to(real, target_is_directory=True)
            with self.assertRaises((binder.TDG9LOC2PREF2Error, OSError)):
                binder._raw_pair(root)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / binder.RAW_NAMESPACE
            target.mkdir(parents=True)
            for name in ("manifest.json", "terminal.json", "foreign"):
                (target / name).write_text("{}\n", encoding="ascii")
            with self.assertRaises(binder.TDG9LOC2PREF2Error):
                binder._raw_pair(root)

    def test_store_traversal_rejects_nested_symlink(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            store = root / "store"
            inner = store / "inner"
            inner.mkdir(parents=True)
            (inner / "leaf").write_text("x", encoding="ascii")
            (inner / "link").symlink_to(inner / "leaf")
            with patch.object(binder.ar1, "PREF2_STORE_PATH", "store"):
                with self.assertRaises(binder.TDG9LOC2PREF2Error):
                    binder._snapshot_store(root)

    def test_route_agreement_requires_digest_keys_and_intervals(self) -> None:
        candidate = SimpleNamespace(
            polynomial_ordinal=1,
            location="interior_stationary",
            location_ordinal=0,
            metadata=(("ordinal", 1),),
            parameter_lower=Fraction(1, 3),
            parameter_upper=Fraction(1, 3),
            absolute_lower=Fraction(2, 3),
            absolute_upper=Fraction(2, 3),
        )
        first = SimpleNamespace(
            polynomial_count=1,
            candidate_count=3,
            classification="unique_maximum",
            candidates=(candidate,),
            global_absolute_lower=Fraction(2, 3),
            global_absolute_upper=Fraction(2, 3),
        )
        second = SimpleNamespace(
            polynomial_count=1,
            candidate_count=3,
            classification="unique_maximum",
            stationary_count_stream_sha256="a" * 64,
            candidates=(candidate,),
            global_absolute_lower=Fraction(2, 3),
            global_absolute_upper=Fraction(2, 3),
        )
        facts = binder._route_agreement(first, second, "a" * 64, "good")
        self.assertTrue(all(facts.values()))
        second.stationary_count_stream_sha256 = "b" * 64
        with self.assertRaises(binder.TDG9LOC2PREF2Error):
            binder._route_agreement(first, second, "a" * 64, "bad")

    def test_contraction_and_component_ownership_are_independent(self) -> None:
        outer = SimpleNamespace(
            global_absolute_lower=Fraction(1),
            global_absolute_upper=Fraction(1),
        )
        finest = SimpleNamespace(
            global_absolute_lower=Fraction(1),
            global_absolute_upper=Fraction(1),
        )
        contraction = binder._contraction_assessment(outer, finest)
        self.assertEqual(
            contraction["classification"],
            "sufficient_contraction_failure",
        )
        components = {
            name: {
                "D01": outer,
                "D12": finest,
                "contraction": contraction,
            }
            for name in ("value_V", "slope_S", "complete_C")
        }
        ownership = binder._component_ownership(components)
        self.assertEqual(
            ownership["classification"],
            "both_components_independently_fail",
        )

    def test_binder_imports_no_loc2_runner(self) -> None:
        source = (
            ROOT
            / "src/recursive_horizons/fgc/evolution/tdg9_loc2_pref2_binder.py"
        ).read_text(encoding="utf-8")
        self.assertNotIn("scripts.run_fgc_tdg9_loc2", source)
        probe = subprocess.run(
            [
                sys.executable,
                "-I",
                "-B",
                "-c",
                (
                    "import sys; sys.path.insert(0, 'src'); "
                    "import recursive_horizons.fgc.evolution.tdg9_loc2_pref2_binder; "
                    "bad=[name for name in sys.modules if "
                    "'run_fgc_tdg9_loc2' in name]; assert not bad, bad"
                ),
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        self.assertEqual(probe.returncode, 0, probe.stdout + probe.stderr)


if __name__ == "__main__":
    unittest.main()
