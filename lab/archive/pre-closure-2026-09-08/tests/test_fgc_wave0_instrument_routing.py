"""Partial Wave-0 instruments must not become a frontier or execution gate."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import json
import sys
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from scripts.build_artifact_catalog import (  # noqa: E402
    CURRENT_FRONTIER,
    NEXT_DESIGN,
    PROSPECTIVE_INSTRUMENTS,
    build_catalog,
)
from scripts.repo_checks import catalog as checks  # noqa: E402


class WaveZeroInstrumentRoutingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = json.loads((ROOT / "configs/fgc/artifact-catalog.json").read_text())

    def test_partial_instruments_are_explicitly_nonexecuting_and_not_compact(self):
        by_id = {entry["artifact_id"]: entry for entry in self.catalog["artifacts"]}
        self.assertEqual(set(PROSPECTIVE_INSTRUMENTS), {
            "FGC-1-SGB1-CTL1", "FGC-1-DEF1-STAB1",
            "FGC-1-TDG11-C1R1", "FGC-1-HLT17-MON17",
            "FGC-1-SGB1-CTL1-TRAP-REFINEMENT2",
            "FGC-1-SGB1-CTL1-TRAP-TAYLOR",
            "FGC-1-SGB1-CTL1-TRAP-BARRIER",
            "FGC-1-SGB1-CTL1-SOL1",
        })
        for artifact_id, specification in PROSPECTIVE_INSTRUMENTS.items():
            entry = by_id[artifact_id]
            self.assertEqual(entry["status"], "prospective_authority")
            self.assertEqual(entry["classification"], specification["classification"])
            for key in ("source_paths", "owner_documents", "test_paths"):
                self.assertEqual(set(entry[key]), set(specification[key]))
                self.assertTrue(all((ROOT / path).is_file() for path in entry[key]))
            for key in ("config_paths", "result_paths", "raw_dependencies", "reproducer_paths"):
                self.assertEqual(entry[key], [])
            self.assertIs(entry["may_execute_state_change"], False)
            self.assertIsNone(entry["explicit_live_target_or_none"])
            self.assertIsNone(entry["compact_verify_target_or_none"])
        self.assertEqual(self.catalog["frontier"], {
            "current_artifact_id": CURRENT_FRONTIER,
            "next_design_artifact_id": NEXT_DESIGN,
            "state_advance_authorized": False,
        })
        self.assertEqual([entry["artifact_id"] for entry in by_id.values()
                          if entry["status"] == "current_frontier"], [CURRENT_FRONTIER])
        sol1_frz1 = by_id["FGC-1-SGB1-CTL1-SOL1-FRZ1"]
        self.assertEqual(sol1_frz1["status"], "compact_active")
        self.assertEqual(
            sol1_frz1["classification"],
            "prospective_sol1_frozen_nominal_interval_inconclusive_no_health",
        )
        self.assertIs(sol1_frz1["may_execute_state_change"], False)
        self.assertIsNone(sol1_frz1["explicit_live_target_or_none"])
        self.assertEqual(
            sol1_frz1["compact_verify_target_or_none"],
            "fgc-sgb1-ctl1-sol1-frz1",
        )
        self.assertEqual(sol1_frz1["raw_dependencies"], [])
        self.assertIn("FGC-1-SGB1-CTL1-SOL1", sol1_frz1["predecessor_ids"])
        self.assertIn("FGC-1-SGB1-CTL1-SOL1-PREF1", sol1_frz1["successor_ids"])
        sol1_pref1 = by_id["FGC-1-SGB1-CTL1-SOL1-PREF1"]
        self.assertEqual(sol1_pref1["status"], "compact_active")
        self.assertEqual(
            sol1_pref1["classification"],
            "independently_bound_sol1_nominal_interval_inconclusive_no_health",
        )
        self.assertIs(sol1_pref1["may_execute_state_change"], False)
        self.assertIsNone(sol1_pref1["explicit_live_target_or_none"])
        self.assertEqual(
            sol1_pref1["compact_verify_target_or_none"],
            "fgc-sgb1-ctl1-sol1-pref1",
        )
        self.assertEqual(sol1_pref1["raw_dependencies"], [])
        self.assertIn("FGC-1-SGB1-CTL1-SOL1-FRZ1", sol1_pref1["predecessor_ids"])
        self.assertNotEqual(sol1_pref1["status"], "current_frontier")
        frz1 = by_id["FGC-1-DEF1-STAB1-FRZ1"]
        self.assertEqual(frz1["status"], "compact_active")
        self.assertEqual(
            frz1["classification"],
            "candidate_blind_conversion_error_map_instrument_freeze_no_trajectory",
        )
        self.assertIs(frz1["may_execute_state_change"], False)
        self.assertIsNone(frz1["explicit_live_target_or_none"])
        self.assertEqual(frz1["compact_verify_target_or_none"], "fgc-def1-stab1-frz1")
        self.assertEqual(frz1["raw_dependencies"], [])
        self.assertIn("FGC-1-DEF1-STAB1", frz1["predecessor_ids"])
        self.assertIn("FGC-1-DEF1-STAB1-PREF1", frz1["successor_ids"])
        self.assertNotEqual(frz1["status"], "current_frontier")
        pref1 = by_id["FGC-1-DEF1-STAB1-PREF1"]
        self.assertEqual(pref1["status"], "compact_active")
        self.assertEqual(
            pref1["classification"],
            "independently_bound_candidate_blind_error_map_readiness_only_no_trajectory",
        )
        self.assertIs(pref1["may_execute_state_change"], False)
        self.assertIsNone(pref1["explicit_live_target_or_none"])
        self.assertEqual(pref1["compact_verify_target_or_none"], "fgc-def1-stab1-pref1")
        self.assertEqual(pref1["raw_dependencies"], [])
        self.assertIn("FGC-1-DEF1-STAB1-FRZ1", pref1["predecessor_ids"])
        self.assertNotEqual(pref1["status"], "current_frontier")

    def test_generated_catalog_and_validator_agree(self):
        self.assertEqual(build_catalog(), self.catalog)
        self.assertEqual(checks.check_artifact_catalog(), [])

    def test_promoted_claims_results_or_execution_are_rejected(self):
        original = checks._read_unique_regular_bytes
        for artifact_id in PROSPECTIVE_INSTRUMENTS:
            for key, value in (
                ("status", "compact_active"),
                ("status", "current_frontier"),
                ("classification", "completed_branch_or_error_map_gate"),
                ("may_execute_state_change", True),
                ("compact_verify_target_or_none", "verify-fgc-wave0-algebra"),
                ("raw_dependencies", ["runs/fgc-2-sf1/forbidden"]),
            ):
                with self.subTest(artifact_id=artifact_id, key=key, value=value):
                    altered = deepcopy(self.catalog)
                    entry = next(item for item in altered["artifacts"]
                                 if item["artifact_id"] == artifact_id)
                    entry[key] = value
                    raw = (json.dumps(altered, indent=2, sort_keys=True) + "\n").encode()

                    def read(relative, *args, **kwargs):
                        if relative == "configs/fgc/artifact-catalog.json":
                            return raw
                        return original(relative, *args, **kwargs)

                    with patch.object(checks, "_read_unique_regular_bytes", side_effect=read):
                        failures = checks.check_artifact_catalog()
                    self.assertIn(f"artifact catalog partial instrument scope differs: {artifact_id}", failures)

    def test_make_route_is_bounded_math_without_a_live_operation(self):
        source = (ROOT / "mk/current-foundation.mk").read_text()
        block = source.split("verify-fgc-wave0-algebra:\n", 1)[1].split("\n\n", 1)[0]
        for module in (
            "test_fgc_sgb1_ctl1_constraints", "test_fgc_sgb1_ctl1_time_affinity",
            "test_fgc_sgb1_ctl1_source", "test_fgc_def1_stab1",
            "test_fgc_sgb1_ctl1_family", "test_fgc_sgb1_ctl1_center",
            "test_fgc_sgb1_ctl1_principal",
            "test_fgc_sgb1_ctl1_interval_health",
            "test_fgc_sgb1_ctl1_runtime",
            "test_fgc_sgb1_ctl1_admission",
            "test_fgc_sgb1_ctl1_initial_health",
            "test_fgc_sgb1_ctl1_controls",
            "test_fgc_sgb1_ctl1_cone",
            "test_fgc_sgb1_ctl1_continuity",
            "test_fgc_sgb1_ctl1_trap_refinement",
            "test_fgc_sgb1_ctl1_trap_refinement2",
            "test_fgc_sgb1_ctl1_trap_taylor",
            "test_fgc_sgb1_ctl1_family_principal",
            "test_fgc_sgb1_ctl1_cell_admission",
            "test_fgc_sgb1_ctl1_local_symmetrizer",
            "test_fgc_sgb1_ctl1_trap_barrier",
            "test_fgc_sgb1_ctl1_sol1",
            "test_fgc_sgb1_ctl1_sol1_frz1_certificate",
            "test_check_repo_sgb1_ctl1_sol1_frz1",
            "test_fgc_sgb1_ctl1_sol1_pref1_binder",
            "test_check_repo_sgb1_ctl1_sol1_pref1",
            "test_fgc_def1_geometry_error", "test_fgc_wave0_instrument_routing",
            "test_fgc_def1_stab1_providers", "test_fgc_def1_stab1_qualification",
            "test_fgc_def1_stab1_freeze_contract",
            "test_fgc_def1_stab1_frz1_certificate",
            "test_check_repo_def1_stab1_frz1",
            "test_fgc_def1_stab1_pref1_binder",
            "test_check_repo_def1_stab1_pref1",
        ):
            self.assertIn(module, block)
        for forbidden in ("scripts/run_", "--live", "--write-result", "runs/"):
            self.assertNotIn(forbidden, block)
        current = next(line for line in source.splitlines()
                       if line.startswith("verify-current-development:"))
        self.assertIn("verify-fgc-wave0-algebra", current)
        self.assertIn("verify-fgc-hlt17-components", current)
        self.assertIn("verify-fgc-pro20-ev1-components", current)
        self.assertRegex(source, r"(?m)^fgc-def1-stab1-frz1:")
        self.assertRegex(source, r"(?m)^fgc-def1-stab1-pref1:")
        self.assertRegex(source, r"(?m)^fgc-sgb1-ctl1-sol1-pref1:")
        self.assertNotRegex(source, r"(?m)^run-fgc-def1-stab1-frz1:")
        self.assertNotRegex(source, r"(?m)^status-fgc-def1-stab1-frz1:")
        self.assertNotRegex(source, r"(?m)^run-fgc-def1-stab1-pref1:")
        self.assertNotRegex(source, r"(?m)^status-fgc-def1-stab1-pref1:")
        self.assertNotRegex(source, r"(?m)^run-fgc-sgb1-ctl1-sol1-pref1:")
        self.assertNotRegex(source, r"(?m)^status-fgc-sgb1-ctl1-sol1-pref1:")
        hlt17 = source.split("verify-fgc-hlt17-components:\n", 1)[1].split("\n\n", 1)[0]
        for module in (
            "test_fgc_tdg11_c1r1_ring", "test_fgc_tdg11_c1r1_enclosure",
            "test_fgc_tdg11_c1r1_runtime", "test_fgc_hlt17_c1r1_integration",
            "test_fgc_protocol_v19", "test_fgc_hlt17_campaign_store",
        ):
            self.assertIn(module, hlt17)
        for forbidden in ("scripts/run_", "--measure", "--live", "runs/"):
            self.assertNotIn(forbidden, hlt17)
        pro20 = source.split("verify-fgc-pro20-ev1-components:\n", 1)[1].split(
            "\n\n", 1
        )[0]
        self.assertIn("test_fgc_pro20_ev1_protocol", pro20)
        self.assertIn("test_fgc_pro20_ev1_store", pro20)
        self.assertIn("test_fgc_pro20_ev1_runtime", pro20)
        self.assertIn("test_fgc_pro20_ev1_authority", pro20)
        for forbidden in ("scripts/run_", "--measure", "--live", "runs/"):
            self.assertNotIn(forbidden, pro20)


if __name__ == "__main__":
    unittest.main()
