from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sys
import tomllib
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.fgc.scoped_run_authorization import (  # noqa: E402
    FUTURE_COMMON_NONCLAIMS,
    FUTURE_EVIDENCE_CLASSIFICATIONS,
    FUTURE_EVIDENCE_REQUIREMENTS,
    FUTURE_PAYLOAD_TRUE_FIELDS,
    RUN1_SYM1_ARTIFACT_ID,
    SF1_PROTOCOL_ARTIFACT_ID,
    SF1_PROTOCOL_V1_ARTIFACT_ID,
    SF1_PROTOCOL_V2_ARTIFACT_ID,
    scoped_classical_spherical_run_audit,
    validate_sf1_protocol,
)


def _load_json(name: str):
    return json.loads((REPOSITORY / "results" / name).read_text(encoding="utf-8"))


class FGCRUN1AuthorizationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.inputs = {
            "action": _load_json("fgc-1-action-gate.json"),
            "variation": _load_json("fgc-1-metric-variation.json"),
            "uhyp1": _load_json("fgc-1-hyp1-dom3-uhyp1.json"),
            "con2": _load_json("fgc-1-hyp1-con2-mprop1.json"),
            "con3": _load_json("fgc-1-hyp1-con3-cau1.json"),
            "bnd1": _load_json("fgc-1-hyp1-bnd1-md1.json"),
            "def0": _load_json("fgc-1-def0-obs1.json"),
            "eft1": _load_json("fgc-1-eft1-open1.json"),
        }
        cls.protocol = tomllib.loads(
            (REPOSITORY / "configs/fgc/fgc-2-sf1-protocol-v3.toml").read_text(
                encoding="utf-8"
            )
        )
        cls.protocol_certificate = validate_sf1_protocol(cls.protocol)

    def audit(self, *, future=None, protocol=None, **changes):
        values = deepcopy(self.inputs)
        values.update(changes)
        return scoped_classical_spherical_run_audit(
            **values,
            protocol=deepcopy(self.protocol if protocol is None else protocol),
            future_evidence=(
                {name: None for name in FUTURE_EVIDENCE_REQUIREMENTS}
                if future is None
                else future
            ),
        )

    def complete_future_evidence(self):
        digest = "0" * 64
        scope = {
            "run1_artifact_id": RUN1_SYM1_ARTIFACT_ID,
            "run1_config_sha256": digest,
            "protocol_artifact_id": SF1_PROTOCOL_ARTIFACT_ID,
            "protocol_config_sha256": digest,
            "protocol_semantic_holdout_contract_sha256": self.protocol_certificate[
                "semantic_holdout_contract_sha256"
            ],
            "action_artifact_id": "FGC-1-ACT1",
            "action_config_sha256": digest,
            "action_result_sha256": digest,
            "variation_artifact_id": "FGC-1-VAR1",
            "variation_config_sha256": digest,
            "variation_result_sha256": digest,
            "target_branch": "FGC-QR",
            "physical_equations": "unredefined_ACT1_VAR1",
            "declared_run_envelope_sha256": digest,
        }
        future = {
            name: {
                "schema_version": 1,
                "artifact_id": artifact_id,
                "project_version": "0.11.0",
                "classification": FUTURE_EVIDENCE_CLASSIFICATIONS[name],
                "generated_by": "scripts/synthetic_future_fixture.py",
                "derivation_document": "docs/synthetic-future-fixture.md",
                "derivation_document_sha256": digest,
                "source_config_sha256": {"configs/synthetic.toml": digest},
                "predecessor_sha256": {"results/synthetic-parent.json": digest},
                "implementation_sha256": {"scripts/synthetic_future_fixture.py": digest},
                "scope_bindings": dict(scope),
                "certificate_contract": {
                    "artifact_specific_payload_validated": True,
                    "canonical_reproduction_passed": True,
                    "scope_bindings_verified": True,
                    "outcome_data_not_used_to_select_certificate_contract": True,
                },
                "artifact_payload": {"synthetic_combinator_fixture": True},
                "gate_status": {gate: True},
                "nonclaims": dict(FUTURE_COMMON_NONCLAIMS),
            }
            for name, (artifact_id, gate) in FUTURE_EVIDENCE_REQUIREMENTS.items()
        }
        for name, true_fields in FUTURE_PAYLOAD_TRUE_FIELDS.items():
            future[name]["artifact_payload"] = {
                "run_envelope_semantic_sha256": digest,
                "quantitative_evidence": {"synthetic_combinator_fixture": True},
                **{field: True for field in true_fields},
            }
        candidates = self.protocol_certificate["amplitude_candidates"]
        future["protocol_holdout"]["artifact_payload"] = {
            "selected_GR0_amplitude": candidates[0],
            "calibration_candidate_ledger": [
                {
                    "amplitude": amplitude,
                    "config_sha256": digest,
                    "result_sha256": digest,
                    "eligible": index == 0,
                    "eligibility_reason": (
                        "first converged eligible GR-0 control"
                        if index == 0
                        else "not selected after the first eligible candidate"
                    ),
                }
                for index, amplitude in enumerate(candidates)
            ],
            "static_all_case_initial_premise_ledger": [
                {
                    "case_id": case_id,
                    "chi_width_factor": width_factor,
                    "phi_seed_factor": seed_factor,
                    "exact_center_buffer_over_L0": center_buffer,
                    "config_sha256": digest,
                    "initial_data_result_sha256": digest,
                    "constraint_compatible": True,
                    "regular_center": True,
                    "finite_mass": True,
                    "initial_compactness_in_protocol_window": True,
                    "no_initial_trapped_sphere": True,
                    "static_premises_passed": True,
                    "FGCQR_evolution_outcome_inspected": False,
                }
                for case_id, width_factor, seed_factor, center_buffer in (
                    ("FGCQR-CENTRAL", "1", "1", "5/2"),
                    ("FGCQR-SEED-HALF", "1", "1/2", "5/2"),
                    ("FGCQR-SEED-DOUBLE", "1", "2", "5/2"),
                    ("FGCQR-WIDTH-SEVEN-EIGHTHS", "7/8", "1", "41/16"),
                    ("FGCQR-WIDTH-NINE-EIGHTHS", "9/8", "1", "39/16"),
                )
            ],
            "ordered_held_out_case_ids": self.protocol_certificate[
                "held_out_case_ids"
            ],
            "held_out_case_config_sha256": {
                case: digest for case in self.protocol_certificate["held_out_case_ids"]
            },
            "repository_sequence_receipt": {
                "preexisting_FGCQR_output_paths": [],
                "runner_refuses_overwrite": True,
                "append_only_event_log_sha256": digest,
                "sealed_before_FGCQR_execution": True,
                "human_private_observation_not_machine_verifiable": True,
            },
        }
        return future

    def test_current_evidence_validates_protocol_but_stops_all_runs(self) -> None:
        output = self.audit()
        self.assertEqual(output["artifact_id"], RUN1_SYM1_ARTIFACT_ID)
        scoped = output["authorization_audits"]["classical_spherical_diagnostic"]
        self.assertEqual(scoped["passed_predicate_count"], 0)
        self.assertEqual(scoped["required_predicate_count"], 8)
        self.assertEqual(len(scoped["missing_predicate_ids"]), 8)
        self.assertTrue(output["protocol_validation"]["outcome_neutral_contract_validated"])
        self.assertEqual(output["protocol_validation"]["artifact_id"], SF1_PROTOCOL_ARTIFACT_ID)
        self.assertEqual(output["protocol_validation"]["protocol_version"], 3)
        self.assertEqual(
            output["protocol_validation"]["exact_widest_center_buffer_over_L0"],
            "39/16",
        )
        self.assertEqual(
            output["protocol_validation"]["widest_case_coarsest_grid_intervals"],
            78,
        )
        self.assertFalse(output["protocol_validation"]["resolved_holdout_manifest_present"])
        self.assertTrue(all(output["validated_predecessor_evidence"].values()))
        self.assertTrue(output["gate_status"]["authorization_contract_and_protocol_validation_passed"])
        for key in (
            "classical_spherical_diagnostic_authorized",
            "retained_EFT_evolution_authorized",
            "physical_transition_claim_authorized",
            "FGCQR_holdout_execution_authorized",
        ):
            self.assertFalse(output["gate_status"][key])

    def test_complete_scoped_evidence_does_not_promote_EFT_or_physical_claim(self) -> None:
        output = self.audit(future=self.complete_future_evidence())
        scoped = output["authorization_audits"]["classical_spherical_diagnostic"]
        self.assertEqual(scoped["passed_predicate_count"], 8)
        self.assertTrue(scoped["authorized"])
        self.assertTrue(output["protocol_validation"]["resolved_holdout_manifest_present"])
        self.assertTrue(output["gate_status"]["classical_spherical_diagnostic_authorized"])
        self.assertTrue(output["gate_status"]["FGCQR_holdout_execution_authorized"])
        self.assertFalse(output["gate_status"]["retained_EFT_evolution_authorized"])
        self.assertFalse(output["gate_status"]["physical_transition_claim_authorized"])
        self.assertFalse(output["stop_mask"]["report_scoped_classical_mechanism_result"])

    def test_each_scoped_predicate_is_all_of_its_constituent_evidence(self) -> None:
        future = self.complete_future_evidence()
        future["multidirectional_health"] = None
        output = self.audit(future=future)
        scoped = output["authorization_audits"]["classical_spherical_diagnostic"]
        self.assertEqual(scoped["passed_predicate_count"], 7)
        self.assertIn(
            "classical_spherical_domain_branch_and_multidirectional_early_kill_envelope",
            scoped["missing_predicate_ids"],
        )
        self.assertFalse(output["gate_status"]["classical_spherical_diagnostic_authorized"])

    def test_protocol_mutations_fail_closed(self) -> None:
        zero_seed = deepcopy(self.protocol)
        zero_seed["initial_data"]["phi_seed"]["amplitude"] = "0"
        with self.assertRaisesRegex(ValueError, "must be positive"):
            validate_sf1_protocol(zero_seed)

        inspected = deepcopy(self.protocol)
        inspected["partition"]["held_out_outcomes_inspected_before_freeze"] = True
        with self.assertRaisesRegex(ValueError, "must be false"):
            validate_sf1_protocol(inspected)

        weak_margin = deepcopy(self.protocol)
        weak_margin["analysis"]["positive_margin_over_error_factor"] = 1
        with self.assertRaisesRegex(ValueError, "must remain four"):
            validate_sf1_protocol(weak_margin)

        rescalable_affine = deepcopy(self.protocol)
        rescalable_affine["analysis"]["affine_initial_normalization"] = (
            "arbitrary_positive_scale"
        )
        with self.assertRaisesRegex(ValueError, "affine normalization"):
            validate_sf1_protocol(rescalable_affine)

        unverifiable_blindness = deepcopy(self.protocol)
        unverifiable_blindness["provenance"][
            "human_private_observation_not_machine_verifiable"
        ] = False
        with self.assertRaisesRegex(ValueError, "must be true"):
            validate_sf1_protocol(unverifiable_blindness)

        sample_only_robustness = deepcopy(self.protocol)
        sample_only_robustness["robustness"][
            "open_neighborhood_requires_interval_or_continuity_certificate"
        ] = False
        with self.assertRaisesRegex(ValueError, "must be true"):
            validate_sf1_protocol(sample_only_robustness)

        one_stop_rejection = deepcopy(self.protocol)
        one_stop_rejection["negative_scope"][
            "single_stop_or_solver_failure_cannot_reject_branch"
        ] = False
        with self.assertRaisesRegex(ValueError, "must be true"):
            validate_sf1_protocol(one_stop_rejection)

        missing_phi_momentum = deepcopy(self.protocol)
        del missing_phi_momentum["initial_data"]["phi_seed"]["unit_normal_momentum"]
        with self.assertRaisesRegex(ValueError, "keys differ"):
            validate_sf1_protocol(missing_phi_momentum)

        outcome_tainted_amendment = deepcopy(self.protocol)
        outcome_tainted_amendment["amendment"][
            "FGCQR_evolution_outcomes_inspected_before_revision"
        ] = True
        with self.assertRaisesRegex(ValueError, "immutable diagnosis provenance differs"):
            validate_sf1_protocol(outcome_tainted_amendment)

        changed_repair_box = deepcopy(self.protocol)
        changed_repair_box["initial_data"]["chi"]["amplitude_candidates"][0] = "3/2"
        with self.assertRaisesRegex(ValueError, "GR0 repair scan"):
            validate_sf1_protocol(changed_repair_box)

        changed_buffer = deepcopy(self.protocol)
        changed_buffer["initial_data"]["chi"]["minimum_center_buffer_over_L0"] = "5/2"
        with self.assertRaisesRegex(ValueError, "centre-buffer repair arithmetic"):
            validate_sf1_protocol(changed_buffer)

        changed_diagnosis_hash = deepcopy(self.protocol)
        changed_diagnosis_hash["amendment"]["diagnosis_result_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "immutable diagnosis provenance differs"):
            validate_sf1_protocol(changed_diagnosis_hash)

        changed_manifest = deepcopy(self.protocol)
        changed_manifest["partition"]["resolved_holdout_manifest"] = (
            "configs/fgc/fgc-1-pro2-hld1.toml"
        )
        with self.assertRaisesRegex(ValueError, "manifest namespace differs"):
            validate_sf1_protocol(changed_manifest)

        changed_output_root = deepcopy(self.protocol)
        changed_output_root["provenance"]["holdout_output_root"] = (
            "runs/fgc-2-sf1/proto2/holdout"
        )
        with self.assertRaisesRegex(ValueError, "output namespace differs"):
            validate_sf1_protocol(changed_output_root)

        removed_static_precheck = deepcopy(self.protocol)
        removed_static_precheck["amendment"][
            "static_all_case_initial_premises_required_before_GR0_dynamic_eligibility"
        ] = False
        with self.assertRaisesRegex(ValueError, "immutable diagnosis provenance differs"):
            validate_sf1_protocol(removed_static_precheck)

    def test_immutable_PROTO1_history_still_validates(self) -> None:
        protocol = tomllib.loads(
            (REPOSITORY / "configs/fgc/fgc-2-sf1-protocol.toml").read_text(
                encoding="utf-8"
            )
        )
        certificate = validate_sf1_protocol(protocol)
        self.assertEqual(certificate["artifact_id"], SF1_PROTOCOL_V1_ARTIFACT_ID)
        self.assertEqual(certificate["protocol_version"], 1)

    def test_immutable_PROTO2_history_still_validates(self) -> None:
        protocol = tomllib.loads(
            (REPOSITORY / "configs/fgc/fgc-2-sf1-protocol-v2.toml").read_text(
                encoding="utf-8"
            )
        )
        certificate = validate_sf1_protocol(protocol)
        self.assertEqual(certificate["artifact_id"], SF1_PROTOCOL_V2_ARTIFACT_ID)
        self.assertEqual(certificate["protocol_version"], 2)

    def test_predecessor_promotion_and_EFT_state_mutation_are_rejected(self) -> None:
        promoted = deepcopy(self.inputs["eft1"])
        promoted["nonclaims"]["evolution_authorized"] = True
        with self.assertRaisesRegex(ValueError, "promoted nonclaim"):
            self.audit(eft1=promoted)

        relaxed = deepcopy(self.inputs["eft1"])
        relaxed["gate_status"]["retained_EFT_open_run_envelope_passed"] = True
        with self.assertRaisesRegex(ValueError, "must be false"):
            self.audit(eft1=relaxed)

    def test_future_artifact_or_gate_substitution_is_rejected(self) -> None:
        future = self.complete_future_evidence()
        future["regular_center"]["artifact_id"] = "FGC-1-ID1-FAM1"
        with self.assertRaisesRegex(ValueError, "artifact identity differs"):
            self.audit(future=future)

        future = self.complete_future_evidence()
        future["protocol_holdout"]["artifact_payload"][
            "static_all_case_initial_premise_ledger"
        ][-1]["static_premises_passed"] = False
        with self.assertRaisesRegex(ValueError, "static initial premise"):
            self.audit(future=future)

        future = self.complete_future_evidence()
        required_gate = FUTURE_EVIDENCE_REQUIREMENTS["regular_center"][1]
        future["regular_center"]["gate_status"][required_gate] = False
        with self.assertRaisesRegex(ValueError, "must be true"):
            self.audit(future=future)

    def test_future_scope_and_common_nonclaim_mutations_are_rejected(self) -> None:
        future = self.complete_future_evidence()
        future["run_domain"]["scope_bindings"][
            "protocol_semantic_holdout_contract_sha256"
        ] = "1" * 64
        with self.assertRaisesRegex(ValueError, "different RUN1 scope"):
            self.audit(future=future)

        future = self.complete_future_evidence()
        del future["run_domain"]["nonclaims"]["dark_sector_mechanism_derived"]
        with self.assertRaisesRegex(ValueError, "common nonclaim"):
            self.audit(future=future)

        future = self.complete_future_evidence()
        future["run_domain"]["scope_bindings"]["declared_run_envelope_sha256"] = "1" * 64
        future["run_domain"]["artifact_payload"]["run_envelope_semantic_sha256"] = "1" * 64
        with self.assertRaisesRegex(ValueError, "different run envelopes"):
            self.audit(future=future)


if __name__ == "__main__":
    unittest.main()
