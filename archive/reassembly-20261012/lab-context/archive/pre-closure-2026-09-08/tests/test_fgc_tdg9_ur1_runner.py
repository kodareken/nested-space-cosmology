"""Focused synthetic controls for the bounded TDG9 UR1 runner."""

from __future__ import annotations

import ast
from copy import deepcopy
from fractions import Fraction
from pathlib import Path
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]

from scripts import run_fgc_tdg9_ur1 as runner  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_ur1_authority as authority  # noqa: E402


AUTHORITY_COMMIT = "a" * 40


def _envelope(
    *,
    subintervals: int,
    raw: float,
    construction: float,
    arithmetic: float,
    bernstein: float,
    certified: float,
) -> runner.tdg5.Binary64CubicEnvelope:
    return runner.tdg5.Binary64CubicEnvelope(
        polynomial_count=subintervals,
        endpoint_candidate_count=2 * subintervals,
        real_interior_root_count=0,
        ambiguous_discriminant_count=0,
        raw_candidate_maximum=raw,
        coefficient_construction_debit=construction,
        outward_arithmetic_debit=arithmetic,
        bernstein_certification_slack=bernstein,
        certified_continuous_upper_bound=certified,
    )


def _zero_lower_interval(
    upper_hex: str, *, subintervals: int
) -> runner.tdg6.TDG6Binary64MagnitudeInterval:
    upper = float.fromhex(upper_hex)
    return runner.tdg6.TDG6Binary64MagnitudeInterval(
        lower_bound=0.0,
        upper_bound=upper,
        subinterval_count=subintervals,
        owned_row_count=1,
        envelope=_envelope(
            subintervals=subintervals,
            raw=0.0,
            construction=0.0,
            arithmetic=0.0,
            bernstein=upper,
            certified=upper,
        ),
    )


def _synthetic_admission() -> runner.tdg6.TDG6ContinuousAdmissionEvidence:
    outer = {
        channel: _zero_lower_interval(
            authority.AC1_U_R_D01_UPPER_HEX, subintervals=2
        )
        for channel in authority.ac1.CHANNEL_ORDER
    }
    finest = {
        channel: _zero_lower_interval(
            authority.AC1_U_R_D12_UPPER_HEX, subintervals=4
        )
        for channel in authority.ac1.CHANNEL_ORDER
    }
    return runner.tdg6.classify_tdg6_runtime_intervals(
        method=runner.COMPARATOR_METHOD,
        outer_intervals=outer,
        finest_intervals=finest,
    )


def _matched_row_stream() -> dict[str, object]:
    return {
        "algorithm": "sha256",
        "domain": authority.ROW_HASH_DOMAIN,
        "row_count": authority.OWNED_ROW_COUNT,
        "observed_sha256": authority.TI2_RETRY3_U_R_ROW_SHA256,
        "expected_sha256": authority.TI2_RETRY3_U_R_ROW_SHA256,
        "matched": True,
        "related_TI2_retry3_u_R_radius_free_complete_C_class": (
            authority.TI2_RETRY3_U_R_RELATED_COMPLETE_C_CLASS
        ),
        "related_not_replacement_production_admission": True,
    }


class UR1RunnerTests(unittest.TestCase):
    def test_manifest_uses_the_bound_ac1_pref1_identity(self) -> None:
        manifest = runner._manifest(AUTHORITY_COMMIT)
        self.assertEqual(
            manifest,
            {
                "schema": "UR1-raw-v1",
                "artifact_id": authority.ARTIFACT_ID,
                "runner_id": "FGC-1-TDG9-UR1-RUN1",
                "authority_commit": AUTHORITY_COMMIT,
                "AC1_PREF1_result_sha256": authority.AC1_PREF1_RESULT_SHA256,
                "experiment_label": authority.SEMANTICS["experiment_label"],
                "tableau_selector": "SSPRK3",
                "actual_spatial_operator": "inherited_RK4_2049_SBP4",
                "retry": 3,
                "published_channel": "u:R",
                "production_SSPRK3_comparator": False,
                "output_leaves": ["manifest.json", "terminal.json"],
            },
        )

    def test_runner_has_no_temporal_admission_commit_or_state_write_path(self) -> None:
        source = (ROOT / authority.RUNNER_PATH).read_text()
        tree = ast.parse(source)
        imported = {
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, (ast.Import, ast.ImportFrom))
            for alias in node.names
        }
        self.assertFalse(any("candidate" in name.lower() for name in imported))
        for forbidden in (
            "require_tdg6_temporal_admission(",
            "commit_tdg6",
            "accept_step(",
            "acquire_writer",
            "publish_checkpoint",
            "write_campaign",
            "fine_path_committed = True",
        ):
            self.assertNotIn(forbidden, source)
        self.assertIn("_restore_replay", source)
        self.assertIn("_prepare_shadow", source)
        self.assertIn("loc1._surface", source)
        self.assertIn("loc1._rows", source)
        self.assertIn("loc1._row_hash", source)
        self.assertEqual(runner.RAW_SCHEMA, "UR1-raw-v1")

    def test_closed_zero_lower_owner_enum_covers_all_named_cases(self) -> None:
        raw_zero = _envelope(
            subintervals=2,
            raw=0.0,
            construction=0.0,
            arithmetic=0.0,
            bernstein=9.0,
            certified=9.0,
        )
        zero_evidence = runner.classify_zero_lower_owner(
            raw_zero, 0.0, label="raw-zero"
        )
        self.assertEqual(zero_evidence["owner"], "raw_candidate_maximum_is_zero")

        debit_clip = _envelope(
            subintervals=2,
            raw=1.0,
            construction=0.75,
            arithmetic=0.25,
            bernstein=128.0,
            certified=129.0,
        )
        debit_evidence = runner.classify_zero_lower_owner(
            debit_clip, 0.0, label="debit-clip"
        )
        self.assertEqual(
            debit_evidence["owner"],
            "coefficient_plus_arithmetic_debit_clips_positive_raw_maximum",
        )
        self.assertEqual(
            runner._fraction_from_rational(
                debit_evidence["exact_unclipped"], label="unclipped"
            ),
            Fraction(0),
        )
        self.assertTrue(
            debit_evidence[
                "bernstein_certification_slack_not_part_of_lower_clip"
            ]
        )
        self.assertEqual(
            runner.name_zero_lower_owner(
                raw=Fraction(1),
                unclipped=Fraction(1, 10**400),
                reproduced=0.0,
            ),
            "downward_binary64_rounding_of_positive_exact_lower",
        )
        self.assertEqual(set(authority.ZERO_LOWER_OWNERS), {
            zero_evidence["owner"],
            debit_evidence["owner"],
            "downward_binary64_rounding_of_positive_exact_lower",
        })

    def test_serialized_u_R_has_exact_ac1_bounds_public_envelopes_and_classifier(self) -> None:
        serialized = runner.serialize_u_R(_synthetic_admission())
        self.assertEqual(serialized["channel"], "u:R")
        self.assertEqual(serialized["classification"], "order_inconclusive")
        self.assertFalse(serialized["admission_passed"])
        self.assertTrue(serialized["matches_sealed_AC1_production_intervals"])
        self.assertTrue(serialized["D01_D12_lower_owners_equal"])
        self.assertEqual(
            serialized["D01"]["lower_bound"]["binary64_hex"], "0x0.0p+0"
        )
        self.assertEqual(
            serialized["D01"]["upper_bound"]["binary64_hex"],
            "0x1.2d198e246e459p-38",
        )
        self.assertEqual(
            serialized["D12"]["upper_bound"]["binary64_hex"],
            "0x1.2d1c31bb91376p-38",
        )
        self.assertEqual(
            set(serialized["D01"]["envelope"]), runner._ENVELOPE_KEYS
        )
        self.assertEqual(
            serialized["D01"]["zero_lower_owner"]["owner"],
            "raw_candidate_maximum_is_zero",
        )

    def test_interval_validator_recomputes_clip_and_binds_owner_to_envelope(self) -> None:
        serialized = runner.serialize_u_R(_synthetic_admission())
        forged = deepcopy(serialized["D01"])
        forged["zero_lower_owner"]["bernstein_certification_slack"] = (
            runner._exact_binary64(0.0, label="forged")
        )
        with self.assertRaisesRegex(runner.UR1RunnerError, "owner_envelope_identity"):
            runner.validate_magnitude_interval(forged, label="D01")

        forged = deepcopy(serialized)
        forged["D01"]["upper_bound"] = runner._exact_binary64(
            float.fromhex("0x1.0p-37"), label="forged.upper"
        )
        forged["D01"]["envelope"]["certified_continuous_upper_bound"] = (
            deepcopy(forged["D01"]["upper_bound"])
        )
        with self.assertRaises(runner.UR1RunnerError):
            runner.validate_u_R(forged)

    def test_row_hash_publishes_related_ti2_class_only_after_exact_match(self) -> None:
        rows = [object()] * authority.OWNED_ROW_COUNT
        with (
            patch.object(runner.loc1, "_surface", return_value="surface"),
            patch.object(runner.loc1, "_rows", return_value=iter(rows)),
            patch.object(
                runner.loc1,
                "_row_hash",
                return_value=authority.TI2_RETRY3_U_R_ROW_SHA256,
            ),
        ):
            matched = runner.hash_u_R_rows(object())
        self.assertTrue(matched["matched"])
        self.assertEqual(matched["row_count"], 2044)
        self.assertEqual(
            matched["related_TI2_retry3_u_R_radius_free_complete_C_class"],
            "sufficient_contraction_pass",
        )
        self.assertTrue(matched["related_not_replacement_production_admission"])

        with (
            patch.object(runner.loc1, "_surface", return_value="surface"),
            patch.object(runner.loc1, "_rows", return_value=iter(rows)),
            patch.object(runner.loc1, "_row_hash", return_value="0" * 64),
        ):
            mismatch = runner.hash_u_R_rows(object())
        self.assertFalse(mismatch["matched"])
        self.assertNotIn(
            "related_TI2_retry3_u_R_radius_free_complete_C_class", mismatch
        )
        self.assertNotIn("related_not_replacement_production_admission", mismatch)

    def test_terminal_reduction_is_closed_and_ordered(self) -> None:
        self.assertEqual(
            runner.reduce_ur1_terminal(
                rows_matched=False,
                intervals_match_ac1=True,
                owners_equal=True,
            ),
            "row_stream_identity_mismatch",
        )
        self.assertEqual(
            runner.reduce_ur1_terminal(
                rows_matched=True,
                intervals_match_ac1=False,
                owners_equal=True,
            ),
            "replayed_production_intervals_differ_from_sealed_AC1",
        )
        self.assertEqual(
            runner.reduce_ur1_terminal(
                rows_matched=True,
                intervals_match_ac1=True,
                owners_equal=False,
            ),
            "completed_u_R_rows_match_TI2_but_D01_D12_lower_owners_differ",
        )
        self.assertEqual(
            runner.reduce_ur1_terminal(
                rows_matched=True,
                intervals_match_ac1=True,
                owners_equal=True,
            ),
            "completed_u_R_rows_match_TI2_and_zero_lowers_owned_by_named_envelope_quantity",
        )

    def test_completed_terminal_authenticates_only_u_R_and_all_nonclaims(self) -> None:
        sealed = (
            authority.SEALED_STORE_LEAF_COUNT,
            authority.SEALED_STORE_SNAPSHOT_SHA256,
        )
        u_r = runner.serialize_u_R(_synthetic_admission())
        classification = runner.reduce_ur1_terminal(
            rows_matched=True,
            intervals_match_ac1=True,
            owners_equal=True,
        )
        terminal = {
            **runner._terminal_base(
                AUTHORITY_COMMIT, classification, sealed, sealed
            ),
            "replay_receipt": runner.expected_replay_receipt(),
            "shadow_path_count": 7,
            "shadow_proposal_count": 7,
            "SSPRK3_stage_and_endpoint_record_count": 28,
            "row_stream": _matched_row_stream(),
            "u_R": u_r,
        }
        validated = runner._validate_terminal(
            terminal, authority_commit=AUTHORITY_COMMIT
        )
        self.assertEqual(set(validated) - set(runner._terminal_base(
            AUTHORITY_COMMIT, classification, sealed, sealed
        )), {
            "replay_receipt",
            "shadow_path_count",
            "shadow_proposal_count",
            "SSPRK3_stage_and_endpoint_record_count",
            "row_stream",
            "u_R",
        })
        self.assertFalse(validated["PDE_state_committed"])
        self.assertFalse(validated["temporal_retry_admission_called"])
        self.assertFalse(validated["fine_path_committed"])
        self.assertFalse(validated["production_method_earned"])
        self.assertTrue(
            validated[
                "completed_result_licenses_only_later_prospectively_frozen_envelope_or_admission_design"
            ]
        )

    def test_status_accepts_existing_raw_while_run_requires_absence(self) -> None:
        receipt = authority.UR1Authority(authority_commit=AUTHORITY_COMMIT)
        sealed = (
            authority.SEALED_STORE_LEAF_COUNT,
            authority.SEALED_STORE_SNAPSHOT_SHA256,
        )
        with (
            patch.object(runner, "_execution_authority", return_value=receipt),
            patch.object(runner, "_snapshot_store", return_value=sealed),
            patch.object(
                runner,
                "_inspect_output",
                return_value={
                    "state": "terminal",
                    "hashes": {"manifest.json": "1" * 64, "terminal.json": "2" * 64},
                    "classification": "shadow_proposal_premise_stop",
                },
            ),
            patch.object(
                authority,
                "require_output_absent",
                side_effect=AssertionError("status required absence"),
            ) as absence,
        ):
            status = runner.status(ROOT, authority_commit=AUTHORITY_COMMIT)
        absence.assert_not_called()
        self.assertFalse(status["safe_to_run"])

        with (
            patch.object(runner, "_execution_authority", return_value=receipt),
            patch.object(
                authority,
                "require_output_absent",
                side_effect=authority.UR1AuthorityError("already exists"),
            ) as absence,
            patch.object(runner, "_snapshot_store") as snapshot,
            patch.object(runner, "_restore_replay") as restore,
        ):
            with self.assertRaises(authority.UR1AuthorityError):
                runner.run(ROOT, authority_commit=AUTHORITY_COMMIT)
        absence.assert_called_once()
        snapshot.assert_not_called()
        restore.assert_not_called()


if __name__ == "__main__":
    unittest.main()
