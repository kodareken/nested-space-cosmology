"""Outcome-blind controls for the independent PREF1 wire/reduction boundary."""

from copy import deepcopy
import ast
from dataclasses import dataclass
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

from recursive_horizons.fgc.evolution import tdg11_msel1_contract as frozen
from recursive_horizons.fgc.evolution import tdg11_msel1_pref1_protocol as p
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (
    CertifiedMagnitudeInterval,
    classify_tdg6_channel,
)


AUTHORITY = "a04265794fe0a3266e34b50fcf034b9dcb462263"
CONFIG = "4d468b3569b34eea77fe5613bfb733c1db05321118ae2dad266474317997b717"
FREEZE = "69acb67a760e0a8e43223e13e2fa6ed1c0234b7e26e9ea93d4f8a5366250a561"
ENVIRONMENT = {"test_fixture_only": "no real measurement"}


def evidence(bounds=(Q(4), Q(4), Q(1), Q(1))):
    l0, u0, l1, u1 = bounds

    def difference(lower, upper, n, ceiling, coefficient):
        return {
            "lower": lower,
            "upper": upper,
            "polynomial_count": n,
            "candidate_count": 2 * n,
            "survivor_count": 2 * n,
            "localization_classification": "nonunique_or_interval_inconclusive",
            "coefficient_stream_sha256": coefficient,
            "survivor_key_stream_sha256": "4" * 64,
            "primary_stationary_count_stream_sha256": "5" * 64,
            "independent_stationary_count_stream_sha256": "5" * 64,
            "primary_evaluator_id": "tdg9_loc1_derivative_monotone_bisection_v1",
            "independent_evaluator_id": "tdg9_loc2_endpoint_deflated_adaptive_discriminant_v2",
            "maximum_candidates": ceiling,
            "routes_agree": True,
        }

    d0 = difference(l0, u0, 2, 16352, "2" * 64)
    d1 = difference(l1, u1, 4, 32704, "3" * 64)
    combined = sha256(
        b"TDG11-RATIONAL-COMPLETE-C-COMBINED-v1\n"
        + ("2" * 64 + "\n" + "3" * 64 + "\n").encode()
    ).hexdigest()
    zero = u0 == 0 and u1 == 0
    passed = not zero and 8 * u1**2 <= l0**2
    failed = not zero and not passed and 8 * l1**2 > u0**2
    return p.wire(
        {
            "d01": d0,
            "d12": d1,
            "decision": classify_tdg6_channel(
                CertifiedMagnitudeInterval(l0, u0), CertifiedMagnitudeInterval(l1, u1)
            ),
            "row_count": 1,
            "expected_row_count": 1,
            "row_stream_sha256": "1" * 64,
            "combined_coefficient_stream_sha256": combined,
            "sufficient_pass_left": 8 * u1**2,
            "sufficient_pass_right": l0**2,
            "sufficient_contraction_pass": passed,
            "sufficient_contraction_failure": failed,
            "threshold_inconclusive": not zero and not passed and not failed,
            "maximum_candidates_D01": 16352,
            "maximum_candidates_D12": 32704,
            "refinement_depth": 160,
            "schema_version": 1,
            "evaluator_id": "tdg11_rational_complete_c_dual_route_v1",
            "absolute_tolerance_used": False,
            "physical_signal_used_for_normalization": False,
            "declared_cubic_is_exact_PDE_history": False,
            "production_authority": False,
        }
    )


def channel(candidate, name, *, outcome="pass"):
    magnitudes = (Q(16), Q(4), Q(1)) if outcome == "pass" else (Q(1), Q(1), Q(1))
    rounds = (Q(1, 7), Q(1, 8), Q(1, 9))
    if candidate == frozen.CANDIDATES[2]:
        e0, e1, e2 = magnitudes
        return p.wire(
            {
                "channel": name,
                "estimators": magnitudes,
                "roundoff_bounds": rounds,
                "contraction_01": classify_tdg6_channel(
                    CertifiedMagnitudeInterval(e0, e0),
                    CertifiedMagnitudeInterval(e1, e1),
                ),
                "contraction_12": classify_tdg6_channel(
                    CertifiedMagnitudeInterval(e1, e1),
                    CertifiedMagnitudeInterval(e2, e2),
                ),
                "public_fine_debit": e2 + rounds[2],
                "global_PDE_enclosure": False,
            }
        )
    bounds = {
        "pass": (Q(4), Q(4), Q(1), Q(1)),
        "nonpass": (Q(1), Q(1), Q(1), Q(1)),
        "inconclusive": (Q(1), Q(4), Q(1), Q(1)),
    }[outcome]
    result = {
        "channel": name,
        "exact_complete_C": evidence(bounds),
        "public_fine_debit": p.wire(bounds[3]),
    }
    if candidate == frozen.CANDIDATES[0]:
        result["roundoff_bounds"] = p.wire(rounds)
        result["public_fine_debit"] = p.wire(bounds[3] + rounds[2])
    return result


def terminal(statuses=("pass", "pass", "pass")):
    widths = []
    for spec in frozen.replay_specs():

        def group(candidate, outcome):
            return {
                "status": "complete",
                "stop": None,
                "channels": [
                    channel(candidate, name, outcome=outcome)
                    for name in frozen.TDG6_COMPLETE_STATE_CHANNELS
                ],
            }

        widths.append(
            {
                "retry": spec["retry"],
                "generation": spec["generation"],
                "width_hex": spec["width_hex"],
                "restored": {
                    "checkpoint_sha256": spec["checkpoint_sha256"],
                    "descriptor_sha256": frozen.DESCRIPTOR_SHA256,
                    "state_sha256": frozen.PHYSICAL_STATE_SHA256,
                    "fingerprint_sha256": "6" * 64,
                },
                "families": {
                    arithmetic: {
                        "status": "complete",
                        "family_sha256": "7" * 64,
                        "accepted_source_prechecks": 7,
                        "stage_and_endpoint_records": 35,
                        "rhs_calls": 42,
                    }
                    for arithmetic in (
                        frozen.ORIGINAL_ARITHMETIC_ID,
                        frozen.COMPENSATED_ARITHMETIC_ID,
                    )
                },
                "baseline": group("baseline", "nonpass"),
                "candidates": {
                    candidate: group(candidate, status)
                    for candidate, status in zip(
                        frozen.CANDIDATES, statuses, strict=True
                    )
                },
            }
        )
    return {
        "schema": p.RAW_SCHEMA,
        "artifact_id": frozen.ARTIFACT_ID,
        "runner_id": p.RUNNER_ID,
        "authority_commit": AUTHORITY,
        "config_sha256": CONFIG,
        "freeze_sha256": FREEZE,
        "environment": ENVIRONMENT.copy(),
        "store_snapshot_before": [115, frozen.STORE_SHA256],
        "store_snapshot_after": [115, frozen.STORE_SHA256],
        "widths": widths,
        "accounting": {
            "accepted_source_prechecks": 42,
            "stage_and_endpoint_records": 210,
            "rhs_calls": 252,
            "static_shells": 6,
        },
        "global_stop": None,
        "elapsed_seconds_hex": (1.0).hex(),
        **p.reduce_widths(widths, expected_rows=1),
        "nonclaims": frozen.nonclaims(),
    }


def validate(value):
    return p.validate_raw_terminal(
        value,
        authority_commit=AUTHORITY,
        config_sha256=CONFIG,
        freeze_sha256=FREEZE,
        environment=ENVIRONMENT,
        expected_rows=1,
    )


def groups_for(width, arithmetic):
    if arithmetic == frozen.ORIGINAL_ARITHMETIC_ID:
        return [
            width["baseline"],
            width["candidates"][frozen.CANDIDATES[0]],
            width["candidates"][frozen.CANDIDATES[2]],
        ]
    return [width["candidates"][frozen.CANDIDATES[1]]]


def refresh(value):
    for key in ("accepted_source_prechecks", "stage_and_endpoint_records", "rhs_calls"):
        value["accounting"][key] = sum(
            family[key]
            for width in value["widths"]
            for family in width["families"].values()
            if family is not None
        )
    value.update(p.reduce_widths(value["widths"], expected_rows=1))


def clear_family(value, index, arithmetic):
    width = value["widths"][index]
    width["families"][arithmetic] = None
    for group in groups_for(width, arithmetic):
        group.update(status="not_attempted", channels=[], stop=None)


class PREF1ProtocolTests(unittest.TestCase):
    def test_large_exact_integer_codec_keeps_global_limit(self):
        before = sys.get_int_max_str_digits()
        for value in (0, 1, -1, 10**100, -(2**32767), 2**32768 - 1):
            self.assertEqual(p.parse_integer_text(p.integer_text(value)), value)
        self.assertEqual(sys.get_int_max_str_digits(), before)
        for bad in ("-0", "+1", "01", "١", "1.0", "", 1, True):
            with self.subTest(bad=bad), self.assertRaises(p.PREF1ProtocolError):
                p.parse_integer_text(bad)
        with self.assertRaises(p.PREF1ProtocolError):
            p.integer_text(2**32768)

    def test_rational_canonicality_and_recursive_serialization(self):
        self.assertEqual(p.fraction(p.wire(Q(1, 3))), Q(1, 3))
        for bad in (
            {"numerator": "2", "denominator": "6"},
            {"numerator": "1", "denominator": "0"},
            {"numerator": "-1", "denominator": "3"},
            {"numerator": True, "denominator": "3"},
        ):
            with self.assertRaises(p.PREF1ProtocolError):
                p.fraction(bad)

        @dataclass(frozen=True)
        class Box:
            magnitude: Q

        self.assertEqual(p.wire(Box(Q(1, 3))), {"magnitude": p.wire(Q(1, 3))})
        self.assertEqual(p.wire(-0.0), {"binary64_hex": "-0x0.0p+0"})
        recursive = []
        recursive.append(recursive)
        for bad in (float("nan"), float("inf"), recursive, {1: "bad key"}):
            with self.assertRaises(p.PREF1ProtocolError):
                p.wire(bad)

    def test_independent_classifier_matches_declared_law_and_typed_flags(self):
        intervals = [
            (Q(0), Q(0)),
            (Q(0), Q(1)),
            (Q(1), Q(1)),
            (Q(1), Q(4)),
            (Q(4), Q(4)),
        ]
        for first in intervals:
            for second in intervals:
                expected = classify_tdg6_channel(
                    CertifiedMagnitudeInterval(*first),
                    CertifiedMagnitudeInterval(*second),
                )
                self.assertEqual(
                    p.wire(p.decision_from_bounds(*first, *second)), p.wire(expected)
                )
        bad = p.wire(p.decision_from_bounds(Q(0), Q(0), Q(0), Q(0)))
        bad["admission_passed"] = 1
        with self.assertRaises(p.PREF1ProtocolError):
            p.validate_decision(bad, (Q(0), Q(0), Q(0), Q(0)))

    def test_sufficient_failure_is_separate_from_threshold_inconclusive(self):
        candidate = frozen.CANDIDATES[0]
        for expected in ("pass", "nonpass", "inconclusive"):
            self.assertEqual(
                p.classify_channel(
                    channel(candidate, "u:alpha", outcome=expected),
                    candidate,
                    expected_rows=1,
                ),
                expected,
            )
        for bounds in ((Q(0), Q(0), Q(0), Q(0)), (Q(0), Q(4), Q(0), Q(2))):
            result = p.validate_complete_c(evidence(bounds), expected_rows=1)
            self.assertTrue(result["decision"]["admission_passed"])
            self.assertFalse(result["sufficient_failure"])

    def test_selection_precedence_and_nonpromotion(self):
        for statuses, selected in (
            (("pass", "pass", "pass"), frozen.CANDIDATES[0]),
            (("nonpass", "pass", "pass"), frozen.CANDIDATES[1]),
            (("nonpass", "nonpass", "pass"), frozen.CANDIDATES[2]),
            (("nonpass", "nonpass", "nonpass"), None),
        ):
            value = terminal(statuses)
            validate(value)
            self.assertEqual(value["selected_candidate"], selected)
            self.assertFalse(value["licenses_only_separate_TDG11_IMP1"])
            self.assertTrue(value["independent_PREF1_required_before_IMP1"])
            self.assertTrue(
                all(
                    item["classified_channels"] == 54
                    for item in value["candidate_summaries"]
                )
            )

    def test_partial_valid_counterexample_differs_from_partial_pass(self):
        for start, expected in (("nonpass", "nonpass"), ("pass", "inconclusive")):
            value = terminal((start, "pass", "pass"))
            group = value["widths"][0]["candidates"][frozen.CANDIDATES[0]]
            group["channels"] = group["channels"][:1]
            group["status"] = "resource_stop"
            group["stop"] = {
                "owner": "exact_localizer",
                "code": "candidate_ceiling_exhausted",
                "detail": "synthetic bounded-prefix control",
            }
            value.update(p.reduce_widths(value["widths"], expected_rows=1))
            validate(value)
            self.assertEqual(value["candidate_summaries"][0]["status"], expected)
            self.assertEqual(value["candidate_summaries"][0]["classified_channels"], 37)
            self.assertFalse(value["candidate_summaries"][0]["complete_all_widths"])
            self.assertEqual(value["selected_candidate"], frozen.CANDIDATES[1])

    def test_authentication_shape_mutations_fail_closed(self):
        baseline = terminal()
        mutations = [
            (("authority_commit",), "0" * 40),
            (("config_sha256",), "0" * 64),
            (("freeze_sha256",), "0" * 64),
            (("environment", "test_fixture_only"), "changed"),
            (("store_snapshot_after", 1), "0" * 64),
            (("licenses_only_separate_TDG11_IMP1",), True),
            (("independent_PREF1_required_before_IMP1",), False),
            (("widths", 0, "generation"), 10),
            (("widths", 1, "width_hex"), frozen.replay_specs()[0]["width_hex"]),
            (("widths", 2, "restored", "checkpoint_sha256"), "0" * 64),
            (("accounting", "rhs_calls"), 251),
            (("accounting", "static_shells"), True),
            (("elapsed_seconds_hex",), "nan"),
        ]
        for path, replacement in mutations:
            value = deepcopy(baseline)
            cursor = value
            for key in path[:-1]:
                cursor = cursor[key]
            cursor[path[-1]] = replacement
            with self.subTest(path=path), self.assertRaises(p.PREF1ProtocolError):
                validate(value)
        value = deepcopy(baseline)
        value["accepted_endpoint"] = {}
        with self.assertRaises(p.PREF1ProtocolError):
            validate(value)
        value = deepcopy(baseline)
        value["widths"][0], value["widths"][1] = value["widths"][1], value["widths"][0]
        with self.assertRaises(p.PREF1ProtocolError):
            validate(value)

    def test_changed_bounds_hashes_counts_and_debits_are_not_labels(self):
        original = channel(frozen.CANDIDATES[0], "u:alpha")
        for key, replacement in (
            ("schema_version", True),
            ("row_count", True),
            ("sufficient_contraction_failure", True),
            ("combined_coefficient_stream_sha256", "0" * 64),
            ("production_authority", True),
            ("absolute_tolerance_used", 0),
        ):
            value = deepcopy(original)
            value["exact_complete_C"][key] = replacement
            with self.subTest(key=key), self.assertRaises(p.PREF1ProtocolError):
                p.classify_channel(value, frozen.CANDIDATES[0], expected_rows=1)
        for key, replacement in (
            ("candidate_count", True),
            ("candidate_count", 1),
            ("survivor_count", 0),
            ("routes_agree", 1),
            ("primary_stationary_count_stream_sha256", "0" * 64),
            ("independent_evaluator_id", "other route"),
        ):
            value = deepcopy(original)
            value["exact_complete_C"]["d01"][key] = replacement
            with self.subTest(key=key), self.assertRaises(p.PREF1ProtocolError):
                p.classify_channel(value, frozen.CANDIDATES[0], expected_rows=1)
        value = deepcopy(original)
        value["roundoff_bounds"][2] = p.wire(Q(0))
        with self.assertRaises(p.PREF1ProtocolError):
            p.classify_channel(value, frozen.CANDIDATES[0], expected_rows=1)
        value = channel(frozen.CANDIDATES[2], "u:alpha")
        value["global_PDE_enclosure"] = True
        with self.assertRaises(p.PREF1ProtocolError):
            p.classify_channel(value, frozen.CANDIDATES[2], expected_rows=1)

    def test_compact_validation_is_io_and_shadow_blind(self):
        value = terminal()
        with (
            patch("builtins.open", side_effect=AssertionError("file open")),
            patch.object(Path, "read_bytes", side_effect=AssertionError("raw read")),
            patch("os.scandir", side_effect=AssertionError("store traversal")),
            patch("subprocess.Popen", side_effect=AssertionError("Git/process")),
        ):
            validate(value)
        tree = ast.parse(Path(p.__file__).read_text())
        forbidden = (
            "numpy",
            "tdg11_msel1_reconstruction",
            "tdg11_rational_complete_c",
            "tdg11_msel1_runtime",
            "tdg11_msel1_authority",
            "run_fgc_tdg11_msel1",
        )
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [node.module or "", *(alias.name for alias in node.names)]
            else:
                continue
            self.assertFalse(
                any(part in name for name in names for part in forbidden), names
            )

    def test_global_stop_forbids_later_work_but_localizer_stop_does_not(self):
        value = terminal()
        stop = {
            "owner": "resource",
            "code": "wall_time_ceiling",
            "detail": "wall_time_ceiling",
        }
        group = value["widths"][0]["candidates"][frozen.CANDIDATES[0]]
        group.update(status="resource_stop", channels=group["channels"][:1], stop=stop)
        value["global_stop"] = stop
        refresh(value)
        with self.assertRaisesRegex(p.PREF1ProtocolError, "after global stop"):
            validate(value)
        # The alarm may instead happen immediately after all groups completed.
        value = terminal()
        value["global_stop"] = stop
        value["elapsed_seconds_hex"] = (14401.0).hex()
        validate(value)

    def test_legitimate_global_prefix_including_eighteen_record_interruption(self):
        for retained in (0, 1, 18):
            value = terminal()
            stop = {
                "owner": "resource",
                "code": "wall_time_ceiling",
                "detail": "wall_time_ceiling",
            }
            value["global_stop"] = stop
            value["elapsed_seconds_hex"] = (14401.0).hex()
            for index in range(3):
                clear_family(value, index, frozen.COMPENSATED_ARITHMETIC_ID)
                if index:
                    clear_family(value, index, frozen.ORIGINAL_ARITHMETIC_ID)
            group = value["widths"][0]["candidates"][frozen.CANDIDATES[0]]
            group.update(
                status="resource_stop", channels=group["channels"][:retained], stop=stop
            )
            value["widths"][0]["candidates"][frozen.CANDIDATES[2]].update(
                status="not_attempted", channels=[], stop=None
            )
            refresh(value)
            validate(value)
            self.assertEqual(value["classification"], p.INCONCLUSIVE)

    def test_completed_family_cannot_have_premise_or_foreign_owned_group_stop(self):
        for status, owner, code in (
            ("premise_stop", "guarded_shadow", "ValueError"),
            ("resource_stop", "guarded_shadow", "ValueError"),
            ("resource_stop", "exact_localizer", "unknown"),
        ):
            value = terminal()
            group = value["widths"][0]["candidates"][frozen.CANDIDATES[0]]
            group.update(
                status=status,
                channels=group["channels"][:1],
                stop={"owner": owner, "code": code, "detail": "synthetic control"},
            )
            refresh(value)
            with (
                self.subTest(status=status, owner=owner),
                self.assertRaises(p.PREF1ProtocolError),
            ):
                validate(value)

    def test_stopped_family_prefix_counts_and_source_retry(self):
        cases = [
            (1, 0, 1, "CompensatedArithmeticStop", True),
            (1, 0, 1, "MSEL1ResourceStop", False),
            (1, 0, 6, "ValueError", True),
            (2, 10, 12, "source_retry", True),
            (7, 35, 0, "ValueError", False),
            (1, 1, 2, "ValueError", False),
            (3, 5, 8, "ValueError", False),
            (1, 5, 7, "source_retry", False),
            (0, 0, 0, "source_retry", False),
        ]
        for prechecks, records, calls, code, valid in cases:
            value = terminal()
            width = value["widths"][0]
            stop = {
                "owner": "guarded_shadow",
                "code": code,
                "detail": "synthetic prefix",
            }
            width["families"][frozen.COMPENSATED_ARITHMETIC_ID] = {
                "status": "premise_stop",
                "stop": stop,
                "accepted_source_prechecks": prechecks,
                "stage_and_endpoint_records": records,
                "rhs_calls": calls,
            }
            for group in groups_for(width, frozen.COMPENSATED_ARITHMETIC_ID):
                group.update(status="premise_stop", channels=[], stop=stop)
            refresh(value)
            with self.subTest(counts=(prechecks, records, calls), code=code):
                if valid:
                    validate(value)
                else:
                    with self.assertRaises(p.PREF1ProtocolError):
                        validate(value)

    def test_zero_rhs_attempt_still_requires_complete_preparation(self):
        value = terminal()
        for index in range(3):
            for arithmetic in (
                frozen.ORIGINAL_ARITHMETIC_ID,
                frozen.COMPENSATED_ARITHMETIC_ID,
            ):
                clear_family(value, index, arithmetic)
        stop = {
            "owner": "resource",
            "code": "MSEL1ResourceStop",
            "detail": "synthetic interruption",
        }
        first = value["widths"][0]
        first["families"][frozen.ORIGINAL_ARITHMETIC_ID] = {
            "status": "resource_stop",
            "stop": stop,
            "accepted_source_prechecks": 0,
            "stage_and_endpoint_records": 0,
            "rhs_calls": 0,
        }
        for group in groups_for(first, frozen.ORIGINAL_ARITHMETIC_ID):
            group.update(status="resource_stop", channels=[], stop=stop)
        value["global_stop"] = {
            "owner": "resource",
            "code": "wall_time_ceiling",
            "detail": "wall_time_ceiling",
        }
        value["elapsed_seconds_hex"] = (14401.0).hex()
        refresh(value)
        validate(value)
        value["accounting"]["static_shells"] = 0
        with self.assertRaises(p.PREF1ProtocolError):
            validate(value)

    def test_embedded_estimator_cannot_report_a_localizer_stop(self):
        value = terminal(("nonpass", "nonpass", "pass"))
        group = value["widths"][0]["candidates"][frozen.CANDIDATES[2]]
        group.update(
            status="resource_stop",
            channels=[],
            stop={
                "owner": "exact_localizer",
                "code": "candidate_ceiling_exhausted",
                "detail": "impossible embedded localizer",
            },
        )
        refresh(value)
        with self.assertRaisesRegex(
            p.PREF1ProtocolError, "never invokes an exact localizer"
        ):
            validate(value)
