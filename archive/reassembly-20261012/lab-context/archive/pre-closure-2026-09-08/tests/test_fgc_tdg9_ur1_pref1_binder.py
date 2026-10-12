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

from recursive_horizons.fgc.evolution import tdg9_ur1_pref1_binder as binder
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (
    CertifiedMagnitudeInterval,
    classify_tdg6_channel,
)


ROOT = Path(__file__).resolve().parents[1]
CONFIG = (ROOT / binder.CONFIG_PATH).read_bytes()
RESULT_PATH = ROOT / binder.RESULT_PATH


def _compact_bytes() -> bytes:
    return RESULT_PATH.read_bytes()


class TDG9UR1PREF1BinderTests(unittest.TestCase):
    def test_compact_result_is_raw_store_and_shadow_blind(self) -> None:
        result = _compact_bytes()
        forbidden = AssertionError("compact verification opened a live input")
        with (
            patch.object(binder, "_raw_snapshot", side_effect=forbidden),
            patch.object(binder, "_snapshot_store", side_effect=forbidden),
            patch.object(binder, "bind_raw_result", side_effect=forbidden),
            patch.object(binder, "_git", side_effect=forbidden),
        ):
            value = binder.validate_compact_result(CONFIG, result)
        replay = value["artifact_payload"]["independent_replay"]
        self.assertEqual(replay["both_owners"], binder.OWNER_NAME)
        self.assertEqual(replay["classifier"]["classification"], "order_inconclusive")
        self.assertTrue(replay["construction_alone_sufficient"])
        self.assertFalse(replay["arithmetic_alone_sufficient"])
        self.assertEqual(
            replay["row_identity"]["sha256"],
            "b4a48c7bf72e015f3f3d9a340ea550d608514fb8cd886562265158c33c0b3969",
        )

    def test_compact_mutations_fail_closed(self) -> None:
        value = json.loads(_compact_bytes())
        mutations = (
            lambda item: item["artifact_payload"]["raw"].__setitem__(
                "classification",
                "completed_u_R_rows_match_TI2_but_D01_D12_lower_owners_differ",
            ),
            lambda item: item["artifact_payload"]["authority_binding"].__setitem__(
                "commit", "0" * 40
            ),
            lambda item: item["artifact_payload"]["authority"].__setitem__(
                "blob_count", 4
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
                "both_owners", "raw_candidate_maximum_is_zero"
            ),
            lambda item: item["artifact_payload"]["independent_replay"].__setitem__(
                "construction_alone_sufficient", False
            ),
            lambda item: item["artifact_payload"]["independent_replay"].__setitem__(
                "arithmetic_alone_sufficient", True
            ),
            lambda item: item["artifact_payload"]["independent_replay"][
                "classifier"
            ].__setitem__("classification", "resolved_order_pass"),
            lambda item: item["artifact_payload"]["independent_replay"][
                "replay_receipt"
            ].__setitem__("retry", 4),
            lambda item: item["artifact_payload"]["independent_replay"][
                "row_identity"
            ].__setitem__("sha256", "0" * 64),
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
                self.assertRaises(binder.TDG9UR1PREF1Error),
            ):
                binder.validate_compact_result(CONFIG, changed_raw)

    def test_json_rejects_duplicate_keys_noncanonical_and_nonfinite(self) -> None:
        with self.assertRaises(binder.TDG9UR1PREF1Error):
            binder._json(b'{"a":1,"a":2}\n', "duplicate")
        with self.assertRaises(binder.TDG9UR1PREF1Error):
            binder._json(b'{"a": 1}\n', "noncanonical")
        with self.assertRaises(binder.TDG9UR1PREF1Error):
            binder._json(b'{"a":NaN}\n', "nan")
        with self.assertRaises(binder.TDG9UR1PREF1Error):
            binder._json(b'{"a":Infinity}\n', "infinity")

    def test_raw_reader_rejects_symlink_extra_and_changed_hash(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / binder.RAW_NAMESPACE
            target.parent.mkdir(parents=True)
            real = root / "real"
            real.mkdir()
            target.symlink_to(real, target_is_directory=True)
            with self.assertRaises((binder.TDG9UR1PREF1Error, OSError)):
                binder._raw_pair(root)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / binder.RAW_NAMESPACE
            target.mkdir(parents=True)
            real = root / "real.json"
            real.write_text("{}\n", encoding="ascii")
            (target / "manifest.json").symlink_to(real)
            (target / "terminal.json").write_text("{}\n", encoding="ascii")
            with self.assertRaises((binder.TDG9UR1PREF1Error, OSError)):
                binder._raw_pair(root)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / binder.RAW_NAMESPACE
            target.mkdir(parents=True)
            for name in ("manifest.json", "terminal.json", "foreign"):
                (target / name).write_text("{}\n", encoding="ascii")
            with self.assertRaises(binder.TDG9UR1PREF1Error):
                binder._raw_pair(root)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / binder.RAW_NAMESPACE
            target.mkdir(parents=True)
            for name in ("manifest.json",):
                (target / name).write_text("{}\n", encoding="ascii")
            with self.assertRaises(binder.TDG9UR1PREF1Error):
                binder._raw_pair(root)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / binder.RAW_NAMESPACE
            target.mkdir(parents=True)
            for name in ("manifest.json", "terminal.json"):
                (target / name).write_text("{}\n", encoding="ascii")
            with self.assertRaises(binder.TDG9UR1PREF1Error):
                binder._raw_pair(root)

    def test_exact_owner_formulas_exclude_bernstein_and_name_debit_clip(self) -> None:
        d01 = binder.reproduce_envelope_lower(
            raw=Fraction.from_float(float.fromhex("0x1.78661cc000000p-43")),
            construction=Fraction.from_float(
                float.fromhex("0x1.2969e2a3be458p-38")
            ),
            arithmetic=Fraction.from_float(float.fromhex("0x1.6b2daa9000016p-90")),
        )
        d12 = binder.reproduce_envelope_lower(
            raw=Fraction.from_float(float.fromhex("0x1.23c13aa500000p-43")),
            construction=Fraction.from_float(
                float.fromhex("0x1.2969e2a67fea0p-38")
            ),
            arithmetic=Fraction.from_float(float.fromhex("0x1.b9799e940001bp-91")),
        )
        self.assertEqual(d01["owner"], binder.OWNER_NAME)
        self.assertEqual(d12["owner"], binder.OWNER_NAME)
        self.assertEqual(d01["construction_over_raw"], binder.D01_CONSTRUCTION_OVER_RAW)
        self.assertEqual(d12["construction_over_raw"], binder.D12_CONSTRUCTION_OVER_RAW)
        self.assertEqual(
            d01["construction_over_raw"],
            Fraction(654019454794891, 25865933815808),
        )
        self.assertEqual(
            d12["construction_over_raw"],
            Fraction(163504863789045, 5012309316608),
        )
        self.assertTrue(d01["construction_alone_sufficient"])
        self.assertTrue(d12["construction_alone_sufficient"])
        self.assertFalse(d01["arithmetic_alone_sufficient"])
        self.assertFalse(d12["arithmetic_alone_sufficient"])
        self.assertEqual(d01["exact_clipped"], Fraction(0))
        self.assertEqual(d12["exact_clipped"], Fraction(0))
        self.assertLessEqual(d01["exact_unclipped"], 0)
        self.assertLessEqual(d12["exact_unclipped"], 0)
        self.assertEqual(d01["reproduced_stored_lower"], 0.0)
        self.assertEqual(d12["reproduced_stored_lower"], 0.0)

        unclipped_without_bernstein = (
            Fraction(1) - Fraction(2, 5) - Fraction(2, 5)
        )
        unclipped_with_bernstein = (
            Fraction(1) - Fraction(2, 5) - Fraction(2, 5) - Fraction(3, 10)
        )
        self.assertEqual(unclipped_without_bernstein, Fraction(1, 5))
        self.assertLess(unclipped_with_bernstein, 0)
        self.assertEqual(
            d01["exact_unclipped"],
            Fraction.from_float(float.fromhex("0x1.78661cc000000p-43"))
            - Fraction.from_float(float.fromhex("0x1.2969e2a3be458p-38"))
            - Fraction.from_float(float.fromhex("0x1.6b2daa9000016p-90")),
        )
        self.assertNotEqual(
            d01["exact_unclipped"],
            Fraction.from_float(float.fromhex("0x1.78661cc000000p-43"))
            - Fraction.from_float(float.fromhex("0x1.2969e2a3be458p-38"))
            - Fraction.from_float(float.fromhex("0x1.6b2daa9000016p-90"))
            - Fraction.from_float(float.fromhex("0x1.0p+0")),
        )
        self.assertEqual(
            binder.name_zero_lower_owner(
                raw=Fraction(0),
                unclipped=Fraction(0),
                reproduced=0.0,
            ),
            "raw_candidate_maximum_is_zero",
        )
        self.assertEqual(
            binder.name_zero_lower_owner(
                raw=Fraction(1),
                unclipped=Fraction(1, 10**400),
                reproduced=0.0,
            ),
            "downward_binary64_rounding_of_positive_exact_lower",
        )

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
        self.assertFalse(decision.order_threshold_resolved)
        self.assertIsNone(decision.order_threshold_passed)
        self.assertEqual(decision.finest_pair_debit, finest.upper)

    def test_closed_terminal_reduction_has_four_arms(self) -> None:
        self.assertEqual(
            binder.reduce_ur1_terminal(
                rows_matched=False,
                intervals_match_ac1=True,
                owners_equal=True,
            ),
            "row_stream_identity_mismatch",
        )
        self.assertEqual(
            binder.reduce_ur1_terminal(
                rows_matched=False,
                intervals_match_ac1=False,
                owners_equal=False,
            ),
            "row_stream_identity_mismatch",
        )
        self.assertEqual(
            binder.reduce_ur1_terminal(
                rows_matched=True,
                intervals_match_ac1=False,
                owners_equal=True,
            ),
            "replayed_production_intervals_differ_from_sealed_AC1",
        )
        self.assertEqual(
            binder.reduce_ur1_terminal(
                rows_matched=True,
                intervals_match_ac1=True,
                owners_equal=True,
            ),
            "completed_u_R_rows_match_TI2_and_zero_lowers_owned_by_named_envelope_quantity",
        )
        self.assertEqual(
            binder.reduce_ur1_terminal(
                rows_matched=True,
                intervals_match_ac1=True,
                owners_equal=False,
            ),
            "completed_u_R_rows_match_TI2_but_D01_D12_lower_owners_differ",
        )

    def test_binder_imports_no_ur1_ac1_ti2_runtime_or_loc1(self) -> None:
        source = (
            ROOT
            / "src/recursive_horizons/fgc/evolution/tdg9_ur1_pref1_binder.py"
        ).read_text(encoding="utf-8")
        forbidden = (
            "tdg9_ur1_authority",
            "tdg9_ac1_authority",
            "tdg9_ac1_pref1_binder",
            "tdg9_ti2_pref1_binder",
            "tdg9_ti2_authority",
            "tdg9_loc2_pref2_binder",
            "tdg6_temporal_admission_runtime",
            "tdg5_stage_complete_refinement_runtime",
            "hlt16_campaign",
            "proto19_gr0_static_factory",
            "run_fgc_tdg9_ur1",
            "run_fgc_tdg9_ac1",
            "run_fgc_tdg9_ti2",
            "run_fgc_tdg9_loc1",
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
        self.assertIn("CertifiedMagnitudeInterval", source)
        self.assertIn("classify_tdg6_channel", source)
        probe = subprocess.run(
            [
                sys.executable,
                "-I",
                "-B",
                "-c",
                (
                    "import sys; sys.path.insert(0, 'src'); "
                    "import recursive_horizons.fgc.evolution.tdg9_ur1_pref1_binder; "
                    "needles=("
                    "'tdg9_ur1_authority',"
                    "'tdg9_ac1_authority',"
                    "'tdg9_ac1_pref1_binder',"
                    "'tdg9_ti2_pref1_binder',"
                    "'tdg9_ti2_authority',"
                    "'tdg6_temporal_admission_runtime',"
                    "'tdg5_stage_complete_refinement_runtime',"
                    "'hlt16_campaign',"
                    "'proto19_gr0_static_factory',"
                    "'run_fgc_tdg9_ur1',"
                    "'run_fgc_tdg9_ac1',"
                    "'run_fgc_tdg9_ti2',"
                    "'run_fgc_tdg9_loc1',"
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
