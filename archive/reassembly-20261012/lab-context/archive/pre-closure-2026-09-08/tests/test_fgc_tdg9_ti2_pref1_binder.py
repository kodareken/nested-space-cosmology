from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
from types import SimpleNamespace
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from recursive_horizons.fgc.evolution import tdg9_ti2_pref1_binder as binder


ROOT = Path(__file__).resolve().parents[1]
CONFIG = (ROOT / binder.CONFIG_PATH).read_bytes()
RESULT_PATH = ROOT / binder.RESULT_PATH


class TDG9TI2PREF1BinderTests(unittest.TestCase):
    def test_compact_result_is_raw_store_and_shadow_blind(self) -> None:
        result = RESULT_PATH.read_bytes()
        forbidden = AssertionError("compact verification opened a live input")
        with (
            patch.object(binder, "_raw_snapshot", side_effect=forbidden),
            patch.object(binder, "_snapshot_store", side_effect=forbidden),
            patch.object(binder, "_restore_shadow", side_effect=forbidden),
            patch.object(binder, "bind_raw_result", side_effect=forbidden),
        ):
            value = binder.validate_compact_result(CONFIG, result)
        replay = value["artifact_payload"]["independent_replay"]
        self.assertEqual(replay["component_level_count"], 60)
        self.assertEqual(replay["shadow_complete_counts"], {
            "failure": 6,
            "inconclusive": 0,
            "pass": 4,
        })

    def test_compact_mutations_fail_closed(self) -> None:
        value = json.loads(RESULT_PATH.read_bytes())
        mutations = (
            lambda item: item["artifact_payload"]["raw"].__setitem__(
                "classification", "completed_all_ten_complete_failures_clear_under_tableau_shadow"
            ),
            lambda item: item["artifact_payload"]["authority_binding"].__setitem__(
                "commit", "0" * 40
            ),
            lambda item: item["artifact_payload"]["store_binding"][
                "binder_snapshot_after"
            ].__setitem__("sha256", "0" * 64),
            lambda item: item["artifact_payload"]["selection"].__setitem__(
                "tableau_selector", "RK4"
            ),
            lambda item: item["artifact_payload"]["selection"].__setitem__(
                "actual_spatial_operator", "SBP2"
            ),
            lambda item: item["artifact_payload"]["scope"].__setitem__(
                "production_SSPRK3_comparator", True
            ),
            lambda item: item["artifact_payload"]["claims"].__setitem__(
                "successor_remedy_selected", True
            ),
            lambda item: item["artifact_payload"]["independent_replay"].__setitem__(
                "component_candidate_counts", {"complete_C": 1, "slope_S": 1, "value_V": 1}
            ),
            lambda item: item["artifact_payload"]["independent_replay"].__setitem__(
                "ownership_transition_matrix", {"forged": 10}
            ),
            lambda item: item["artifact_payload"]["independent_replay"][
                "occurrences"
            ].reverse(),
            lambda item: item["artifact_payload"]["independent_replay"][
                "occurrences"
            ][2].__setitem__("shadow_ownership_class", "complete_not_failure"),
            lambda item: item["artifact_payload"]["conclusion"].__setitem__(
                "physics_inference_permitted", True
            ),
        )
        for mutation in mutations:
            changed = deepcopy(value)
            mutation(changed)
            changed_raw = binder.canonical_result(changed)
            with (
                patch.object(
                    binder,
                    "COMPACT_RESULT_SHA256",
                    sha256(changed_raw).hexdigest(),
                ),
                self.assertRaises(binder.TDG9TI2PREF1Error),
            ):
                binder.validate_compact_result(CONFIG, changed_raw)

    def test_json_rejects_duplicate_keys_and_noncanonical_bytes(self) -> None:
        with self.assertRaises(binder.TDG9TI2PREF1Error):
            binder._json(b'{"a":1,"a":2}\n', "duplicate")
        with self.assertRaises(binder.TDG9TI2PREF1Error):
            binder._json(b'{"a": 1}\n', "noncanonical")

    def test_raw_reader_rejects_symlink_extra_and_changed_hash(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / binder.RAW_NAMESPACE
            target.parent.mkdir(parents=True)
            real = root / "real"
            real.mkdir()
            target.symlink_to(real, target_is_directory=True)
            with self.assertRaises((binder.TDG9TI2PREF1Error, OSError)):
                binder._raw_pair(root)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / binder.RAW_NAMESPACE
            target.mkdir(parents=True)
            for name in ("manifest.json", "terminal.json", "foreign"):
                (target / name).write_text("{}\n", encoding="ascii")
            with self.assertRaises(binder.TDG9TI2PREF1Error):
                binder._raw_pair(root)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / binder.RAW_NAMESPACE
            target.mkdir(parents=True)
            for name in ("manifest.json", "terminal.json"):
                (target / name).write_text("{}\n", encoding="ascii")
            with self.assertRaises(binder.TDG9TI2PREF1Error):
                binder._raw_pair(root)

    def test_contraction_and_ownership_reduction_are_independent(self) -> None:
        outer = SimpleNamespace(
            global_absolute_lower=Fraction(1),
            global_absolute_upper=Fraction(1),
        )
        finest = SimpleNamespace(
            global_absolute_lower=Fraction(1),
            global_absolute_upper=Fraction(1),
        )
        contraction = binder.loc2_pref2._contraction_assessment(outer, finest)
        self.assertEqual(
            contraction["classification"], "sufficient_contraction_failure"
        )
        components = {
            name: {"D01": outer, "D12": finest, "contraction": contraction}
            for name in ("value_V", "slope_S", "complete_C")
        }
        ownership = binder.loc2_pref2._component_ownership(components)
        self.assertEqual(
            ownership["classification"], "both_components_independently_fail"
        )

    def test_binder_imports_no_ti2_runner(self) -> None:
        source = (
            ROOT
            / "src/recursive_horizons/fgc/evolution/tdg9_ti2_pref1_binder.py"
        ).read_text(encoding="utf-8")
        self.assertNotIn("scripts.run_fgc_tdg9_ti2", source)
        probe = subprocess.run(
            [
                sys.executable,
                "-I",
                "-B",
                "-c",
                (
                    "import sys; sys.path.insert(0, 'src'); "
                    "import recursive_horizons.fgc.evolution.tdg9_ti2_pref1_binder; "
                    "bad=[name for name in sys.modules if "
                    "'run_fgc_tdg9_ti2' in name]; assert not bad, bad"
                ),
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        self.assertEqual(probe.returncode, 0, probe.stdout + probe.stderr)


if __name__ == "__main__":
    unittest.main()
