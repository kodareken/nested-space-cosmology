from __future__ import annotations

import ast
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from recursive_horizons.fgc.evolution import tdg9_ac1_pref1_binder as binder
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (
    CertifiedMagnitudeInterval,
    classify_tdg6_channel,
)


ROOT = Path(__file__).resolve().parents[1]
CONFIG = (ROOT / binder.CONFIG_PATH).read_bytes()
RESULT_PATH = ROOT / binder.RESULT_PATH


class TDG9AC1PREF1BinderTests(unittest.TestCase):
    def test_compact_result_is_raw_store_and_runner_blind(self) -> None:
        result = RESULT_PATH.read_bytes()
        forbidden = AssertionError("compact verification opened a live input")
        with (
            patch.object(binder, "_raw_snapshot", side_effect=forbidden),
            patch.object(binder, "_snapshot_store", side_effect=forbidden),
            patch.object(binder, "bind_raw_result", side_effect=forbidden),
        ):
            value = binder.validate_compact_result(CONFIG, result)
        replay = value["artifact_payload"]["independent_replay"]
        self.assertEqual(replay["admitted_count"], 17)
        self.assertEqual(replay["failed_channels"], ["u:R"])
        self.assertEqual(
            replay["channels"][3]["classification"], "order_inconclusive"
        )

    def test_compact_mutations_fail_closed(self) -> None:
        value = json.loads(RESULT_PATH.read_bytes())
        mutations = (
            lambda item: item["artifact_payload"]["raw"].__setitem__(
                "classification", "completed_retry3_all_18_channels_pass"
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
                "failed_channels", ["u:alpha"]
            ),
            lambda item: item["artifact_payload"]["independent_replay"][
                "channels"
            ].reverse(),
            lambda item: item["artifact_payload"]["independent_replay"]["channels"][
                3
            ].__setitem__("classification", "resolved_order_pass"),
            lambda item: item["artifact_payload"]["independent_replay"][
                "replay_receipt"
            ].__setitem__("retry", 4),
            lambda item: item["artifact_payload"]["conclusion"].__setitem__(
                "physics_inference_permitted", True
            ),
            lambda item: item["artifact_payload"]["conclusion"].__setitem__(
                "production_method_earned", True
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
                self.assertRaises(binder.TDG9AC1PREF1Error),
            ):
                binder.validate_compact_result(CONFIG, changed_raw)

    def test_json_rejects_duplicate_keys_noncanonical_and_nonfinite(self) -> None:
        with self.assertRaises(binder.TDG9AC1PREF1Error):
            binder._json(b'{"a":1,"a":2}\n', "duplicate")
        with self.assertRaises(binder.TDG9AC1PREF1Error):
            binder._json(b'{"a": 1}\n', "noncanonical")
        with self.assertRaises(binder.TDG9AC1PREF1Error):
            binder._json(b'{"a":NaN}\n', "nan")
        with self.assertRaises(binder.TDG9AC1PREF1Error):
            binder._json(b'{"a":Infinity}\n', "infinity")

    def test_raw_reader_rejects_symlink_extra_and_changed_hash(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / binder.RAW_NAMESPACE
            target.parent.mkdir(parents=True)
            real = root / "real"
            real.mkdir()
            target.symlink_to(real, target_is_directory=True)
            with self.assertRaises((binder.TDG9AC1PREF1Error, OSError)):
                binder._raw_pair(root)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / binder.RAW_NAMESPACE
            target.mkdir(parents=True)
            for name in ("manifest.json", "terminal.json", "foreign"):
                (target / name).write_text("{}\n", encoding="ascii")
            with self.assertRaises(binder.TDG9AC1PREF1Error):
                binder._raw_pair(root)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / binder.RAW_NAMESPACE
            target.mkdir(parents=True)
            for name in ("manifest.json",):
                (target / name).write_text("{}\n", encoding="ascii")
            with self.assertRaises(binder.TDG9AC1PREF1Error):
                binder._raw_pair(root)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / binder.RAW_NAMESPACE
            target.mkdir(parents=True)
            for name in ("manifest.json", "terminal.json"):
                (target / name).write_text("{}\n", encoding="ascii")
            with self.assertRaises(binder.TDG9AC1PREF1Error):
                binder._raw_pair(root)

    def test_independent_u_R_reclassification_is_order_inconclusive(self) -> None:
        outer = CertifiedMagnitudeInterval(
            Fraction(0),
            Fraction.from_float(float.fromhex("0x1.2d198e246e459p-38")),
        )
        finest = CertifiedMagnitudeInterval(
            Fraction(0),
            Fraction.from_float(float.fromhex("0x1.2d1c31bb91376p-38")),
        )
        decision = classify_tdg6_channel(outer, finest)
        self.assertEqual(decision.classification, "order_inconclusive")
        self.assertFalse(decision.admission_passed)
        self.assertTrue(decision.temporal_retry_permitted)
        self.assertEqual(decision.finest_pair_debit, finest.upper)
        self.assertEqual(
            binder.reduce_retry3_terminal(
                complete_admission_passed=False,
                failed_channels=("u:R",),
            ),
            "completed_retry3_one_or_more_channels_fail",
        )
        self.assertEqual(
            binder.reduce_retry3_terminal(
                complete_admission_passed=True,
                failed_channels=(),
            ),
            "completed_retry3_all_18_channels_pass",
        )

    def test_binder_imports_no_ac1_runner_or_ti2_shadow(self) -> None:
        source = (
            ROOT
            / "src/recursive_horizons/fgc/evolution/tdg9_ac1_pref1_binder.py"
        ).read_text(encoding="utf-8")
        forbidden = (
            "tdg9_ac1_authority",
            "tdg9_ar1_authority",
            "tdg9_loc2_pref2_binder",
            "tdg9_ti2_pref1_binder",
            "tdg9_ti2_authority",
            "tdg6_temporal_admission_runtime",
            "hlt16_campaign",
            "proto19_gr0_static_factory",
            "run_fgc_tdg9_ac1",
            "run_fgc_tdg9_ti2",
            "serialize_continuous_admission",
            "validate_all_of_identity",
            "prepare_tdg6_gr0_compositor",
            "_restore_replay",
        )
        imported: set[str] = set()
        for node in ast.walk(ast.parse(source)):
            if isinstance(node, ast.Import):
                imported.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                if node.module:
                    imported.add(node.module)
                imported.update(alias.name for alias in node.names)
        self.assertFalse(
            {
                name
                for name in imported
                if any(needle in name for needle in forbidden)
            }
        )
        for line in source.splitlines():
            stripped = line.strip()
            if stripped.startswith(("import ", "from ")):
                for needle in forbidden:
                    self.assertNotIn(needle, stripped)
        probe = subprocess.run(
            [
                sys.executable,
                "-I",
                "-B",
                "-c",
                (
                    "import sys; sys.path.insert(0, 'src'); "
                    "import recursive_horizons.fgc.evolution.tdg9_ac1_pref1_binder; "
                    "needles=("
                    "'tdg9_ac1_authority',"
                    "'tdg9_ti2_pref1_binder',"
                    "'tdg9_ti2_authority',"
                    "'tdg9_loc2_pref2_binder',"
                    "'tdg6_temporal_admission_runtime',"
                    "'hlt16_campaign',"
                    "'proto19_gr0_static_factory',"
                    "'run_fgc_tdg9_ac1',"
                    "'run_fgc_tdg9_ti2',"
                    "); "
                    "bad=[name for name in sys.modules if any("
                    "needle in name for needle in needles)]; "
                    "assert not bad, bad"
                ),
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        self.assertEqual(probe.returncode, 0, probe.stdout + probe.stderr)


if __name__ == "__main__":
    unittest.main()
