"""Compact IMP1 schema/publication dispatch controls, with no trajectories."""

from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons import evidence_io as io  # noqa: E402
from scripts import reproduce_fgc_tdg11_imp1 as certificate  # noqa: E402


def _report():
    """Validator unit fixture, never a measured or publishable certificate."""

    methods = []
    for index, method in enumerate(certificate.METHODS):
        stages = 5 if index == 0 else 4
        methods.append(
            {
                "method": method,
                "point_count": 9,
                "owned_row_count": 4,
                "channels": list(certificate.CHANNELS),
                "shadow_stage_records": 7 * stages,
                "next_family_shadow_stage_records": 7 * stages,
                "total_shadow_proposals": 14,
                "committed_proposals": 4,
                "committed_stage_records": 4 * stages,
                "independent_reference_channels": 18,
                "endpoint_sha256": "1" * 64,
                "assessment_sha256": "2" * 64,
                "checkpoint_sha256": "3" * 64,
                **{name: True for name in certificate.CONTROL_FLAGS},
            }
        )
    sizes = [
        {
            "method": method,
            "point_count": count,
            "owned_row_count": count - 5,
            "channel_count": 18,
            "admission_passed": True,
            "assessment_sha256": "4" * 64,
            "shadow_stage_records": 35 if method == certificate.METHODS[0] else 28,
            "physical_source_tested": False,
        }
        for method in certificate.METHODS
        for count in (129, 2049)
    ]
    return {
        "method_controls": methods,
        "saturation_control": {
            "method": certificate.METHODS[0],
            "channel": "u:phi",
            "initial_value_hex": (2.0**52).hex(),
            "width_hex": (0.25).hex(),
            "actual_endpoint_unchanged": True,
            "raw_D01": {"numerator": "1", "denominator": "36"},
            "raw_D12": {"numerator": "1", "denominator": "72"},
            "corrected_D01": {"numerator": "0", "denominator": "1"},
            "corrected_D12": {"numerator": "0", "denominator": "1"},
            "public_fine_debit": {"numerator": "1", "denominator": "4"},
            "raw_failure_not_relabeled": True,
        },
        "size_controls": sizes,
    }


def _payload(report):
    config = {
        "implementation": [],
        "inherited": [],
        "environment": {"synthetic": "unit_fixture"},
    }
    raw = b"unit-test configuration only"
    result = {
        "schema_version": 1,
        "artifact_id": certificate.ARTIFACT_ID,
        "classification": certificate.CLASSIFICATION,
        "evaluator_id": certificate.EVALUATOR_ID,
        "config_sha256": sha256(raw).hexdigest(),
        "selection_commit": certificate.PREF1_COMMIT,
        "selection_result_sha256": certificate.PREF1_SHA256,
        "implementation": [],
        "inherited": [],
        "limits": certificate.LIMITS,
        "scope": certificate.SCOPE,
        "environment": config["environment"],
        "qualification": report,
        "qualification_sha256": sha256(io.canonical_json_bytes(report)).hexdigest(),
    }
    return config, raw, result


class IMP1ReproductionTests(unittest.TestCase):
    def test_small_independent_reference_controls_execute_the_real_interfaces(self):
        from recursive_horizons.fgc.evolution.tdg11_imp1_qualification import (
            _method_control,
            _saturation_control,
        )

        for method in certificate.METHODS:
            control = _method_control(method)
            certificate._validate_method_control(control, method)
        self.assertEqual(_saturation_control(), _report()["saturation_control"])

    def test_report_schema_requires_both_methods_all_channels_and_sizes(self):
        report = _report()
        certificate.validate_qualification_report(report)
        mutations = []
        first = deepcopy(report)
        first["method_controls"].pop()
        mutations.append(first)
        first = deepcopy(report)
        first["method_controls"][0]["channels"].pop()
        mutations.append(first)
        first = deepcopy(report)
        first["method_controls"].reverse()
        mutations.append(first)
        first = deepcopy(report)
        first["size_controls"].pop()
        mutations.append(first)
        first = deepcopy(report)
        first["method_controls"][0]["committed_proposals"] = 1
        mutations.append(first)
        for value in mutations:
            with self.assertRaises(certificate.IMP1CertificateError):
                certificate.validate_qualification_report(value)

    def test_boolean_aliases_and_false_claims_are_rejected(self):
        for flag in certificate.CONTROL_FLAGS:
            for value in (False, 1, "true"):
                report = _report()
                report["method_controls"][0][flag] = value
                with self.subTest(flag=flag, value=value):
                    with self.assertRaises(certificate.IMP1CertificateError):
                        certificate.validate_qualification_report(report)

    def test_rounding_debit_cannot_be_erased_or_raw_order_relabelled(self):
        for key, value in (
            ("public_fine_debit", {"numerator": "0", "denominator": "1"}),
            ("raw_D01", {"numerator": "0", "denominator": "1"}),
            ("raw_failure_not_relabeled", False),
        ):
            report = _report()
            report["saturation_control"][key] = value
            with self.assertRaises(certificate.IMP1CertificateError):
                certificate.validate_qualification_report(report)

    def test_qualification_has_closed_scope_and_no_endpoint_payload(self):
        for field, value in (
            ("physical_source_tested", True),
            ("state", [0.0]),
            ("point_count", 1025),
        ):
            report = _report()
            report["size_controls"][0][field] = value
            with self.assertRaises(certificate.IMP1CertificateError):
                certificate.validate_qualification_report(report)

    def test_placeholder_cannot_emit_a_success_config(self):
        with patch.object(certificate, "QUALIFICATION_SHA256", None):
            with self.assertRaisesRegex(
                certificate.IMP1CertificateError, "not yet been bound"
            ):
                certificate.emit_config_bytes(Path("/no-live-input"))

    def test_result_commitment_rejects_other_wellformed_hashes(self):
        report = _report()
        config, config_raw, result = _payload(report)
        commitment = result["qualification_sha256"]
        with (
            patch.object(certificate, "validate_config", return_value=config),
            patch.object(certificate, "QUALIFICATION_SHA256", commitment),
        ):
            certificate.validate_result(config_raw, io.canonical_json_bytes(result))
            altered = deepcopy(result)
            altered["qualification"]["method_controls"][0]["endpoint_sha256"] = "9" * 64
            altered["qualification_sha256"] = sha256(
                io.canonical_json_bytes(altered["qualification"])
            ).hexdigest()
            with self.assertRaises(certificate.IMP1CertificateError):
                certificate.validate_result(
                    config_raw, io.canonical_json_bytes(altered)
                )

    def test_result_scope_or_schema_promotion_fails(self):
        report = _report()
        config, config_raw, result = _payload(report)
        with (
            patch.object(certificate, "validate_config", return_value=config),
            patch.object(
                certificate, "QUALIFICATION_SHA256", result["qualification_sha256"]
            ),
        ):
            for field in certificate.SCOPE:
                if certificate.SCOPE[field] is False:
                    changed = deepcopy(result)
                    changed["scope"][field] = True
                    with self.subTest(field=field):
                        with self.assertRaises(certificate.IMP1CertificateError):
                            certificate.validate_result(
                                config_raw, io.canonical_json_bytes(changed)
                            )
            result["schema_version"] = True
            with self.assertRaises(certificate.IMP1CertificateError):
                certificate.validate_result(config_raw, io.canonical_json_bytes(result))

    def test_compact_dispatch_reads_only_declared_compact_inputs(self):
        report = _report()
        config, config_raw, result = _payload(report)
        reads = []

        def reader(root, path):
            reads.append(path)
            if path == certificate.CONFIG_PATH:
                return config_raw
            if path == certificate.RESULT_PATH:
                return io.canonical_json_bytes(result)
            raise AssertionError(f"unexpected compact input {path}")

        with (
            patch.object(certificate, "validate_config", return_value=config),
            patch.object(
                certificate, "QUALIFICATION_SHA256", result["qualification_sha256"]
            ),
            patch.object(io, "read_regular_file", side_effect=reader),
            patch.object(
                certificate, "qualify", side_effect=AssertionError("synthetic run")
            ),
        ):
            certificate.verify_compact(Path("/no-live-input"))
        self.assertEqual(reads, [certificate.CONFIG_PATH, certificate.RESULT_PATH])

    def test_existing_destination_refuses_before_configuration_or_work(self):
        with (
            patch.object(certificate.os.path, "lexists", return_value=True),
            patch.object(io, "read_regular_file", side_effect=AssertionError("read")),
            patch.object(certificate, "qualify", side_effect=AssertionError("work")),
        ):
            with self.assertRaisesRegex(
                certificate.IMP1CertificateError, "already exists"
            ):
                certificate.main(["--qualify", "--write-result"])

    def test_write_flag_without_qualification_is_invalid(self):
        with patch("sys.stderr"):
            with self.assertRaises(SystemExit) as caught:
                certificate.main(["--check", "--write-result"])
        self.assertEqual(caught.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
