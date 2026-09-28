#!/usr/bin/env python3
"""Check compact IMP1 by default; explicitly qualify synthetic states only."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import sys
import tomllib
from typing import Callable

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons import evidence_io as io  # noqa: E402


ARTIFACT_ID = "FGC-1-TDG11-IMP1"
CONFIG_PATH = "configs/fgc/fgc-1-tdg11-imp1.toml"
RESULT_PATH = "results/fgc-1-tdg11-imp1.json"
PREF1_PATH = "results/fgc-1-tdg11-msel1-pref1.json"
PREF1_CONFIG_PATH = "configs/fgc/fgc-1-tdg11-msel1-pref1.toml"
PREF1_CONFIG_SHA256 = "7f5dd59a0a6026a52a366f0965cf851fb2b0e0434c2d15bc5f585c4cb1e0ad5a"
PREF1_SHA256 = "542c195b77bf79b317dbebec88d08f2833d51cf171c99b87ea88bc3375ae9632"
PREF1_COMMIT = "c69b8777c3a5d5862c5ebccc0924bacc9c6ac669"
EVALUATOR_ID = "tdg11_imp1_exact_bernstein_with_dual_fallback_v1"
CLASSIFICATION = "synthetic_runtime_qualification_passed_no_campaign"
# Filled only from a completed, independently checked deterministic synthetic
# report before this new owner is sealed. A placeholder cannot publish a pass.
QUALIFICATION_SHA256: str | None = (
    "a3170dd4fcd97aba6f2787b48f0addbdd78155559d4e6a753793b44bab82944e"
)
SOURCE_PREFIX = "src/recursive_horizons/fgc/evolution/"
METHODS = (
    "fourth_order_diagonal_norm_SBP_with_second_order_boundary_closure_plus_RK4",
    "second_order_diagonal_norm_SBP_plus_SSPRK3",
)
CHANNELS = tuple(
    f"{block}:{field}"
    for block in ("u", "p", "q")
    for field in ("alpha", "v", "lambda", "R", "phi", "chi")
)
IMPLEMENTATION_PATHS = (
    "scripts/reproduce_fgc_tdg11_imp1.py",
    SOURCE_PREFIX + "tdg11_imp1_enclosure.py",
    SOURCE_PREFIX + "tdg11_imp1_ledger.py",
    SOURCE_PREFIX + "tdg11_imp1_runtime.py",
    SOURCE_PREFIX + "tdg11_imp1_qualification.py",
    "tests/test_fgc_tdg11_imp1_enclosure.py",
    "tests/test_fgc_tdg11_imp1_ledger.py",
    "tests/test_fgc_tdg11_imp1_runtime.py",
    "tests/test_fgc_tdg11_imp1_reproduction.py",
)
INHERITED_PATHS = (
    "src/recursive_horizons/evidence_io.py",
    SOURCE_PREFIX + "numerical_engine.py",
    SOURCE_PREFIX + "proto5_runtime.py",
    SOURCE_PREFIX + "proto6_runtime.py",
    SOURCE_PREFIX + "proto7_runtime.py",
    SOURCE_PREFIX + "proto19_gr0_static_factory.py",
    SOURCE_PREFIX + "boundary_domain.py",
    SOURCE_PREFIX + "tdg5_stage_complete_refinement_runtime.py",
    SOURCE_PREFIX + "tdg6_temporal_admission_design.py",
    SOURCE_PREFIX + "tdg6_temporal_admission_runtime.py",
    SOURCE_PREFIX + "tdg7_binary64_subdivision_lattice.py",
    SOURCE_PREFIX + "tdg9_exact_temporal_arithmetic.py",
    SOURCE_PREFIX + "tdg9_local_extrema.py",
    SOURCE_PREFIX + "tdg9_local_extrema_independent_v2.py",
    SOURCE_PREFIX + "tdg11_msel1_reconstruction.py",
    SOURCE_PREFIX + "tdg11_rational_complete_c.py",
    SOURCE_PREFIX + "tdg11_msel1_pref1_reconstruction.py",
    SOURCE_PREFIX + "tdg11_msel1_pref1_localization.py",
    SOURCE_PREFIX + "tdg11_msel1_authority.py",
)
LIMITS = {
    "maximum_owned_rows": 2044,
    "maximum_rational_bits": 32768,
    "maximum_subdivision_depth": 12,
    "maximum_subdivisions_per_pair": 32768,
    "fallback_maximum_candidates_D01": 16352,
    "fallback_maximum_candidates_D12": 32704,
    "fallback_refinement_depth": 160,
    "maximum_rejection_records": 65536,
    "maximum_checkpoint_bytes": 16 * 1024 * 1024,
}
SCOPE = {
    "synthetic_states_only": True,
    "both_original_integrators_tested": True,
    "full_accumulation_debit_retained": True,
    "binary64_fine_endpoint_only": True,
    "inherited_TDG6_snapshot_preserved": True,
    "strict_in_memory_checkpoint_codec": True,
    "corrected_endpoint_adopted": False,
    "raw_store_or_campaign_read": False,
    "physical_source_operator_binding_qualified": False,
    "external_durable_cursor_implemented": False,
    "HLT17_qualified": False,
    "PROTO19_frozen": False,
    "production_state_advanced": False,
    "GR0_common_event_completed": False,
    "GR0_case_eligible": False,
    "SGBL_execution_authorized": False,
    "FGCQR_holdout_execution_authorized": False,
    "global_PDE_error_theorem_proved": False,
    "physical_claim_authorized": False,
}
CONTROL_FLAGS = (
    "admission_passed",
    "inherited_snapshot_preserved",
    "exact_debit_added",
    "fine_endpoint_adopted",
    "checkpoint_roundtrip_identical",
    "next_macro_prepared",
    "independent_reference_agreement",
)


class IMP1CertificateError(ValueError):
    """A partial, altered or overclaimed implementation certificate."""


def _fail(message: str) -> None:
    raise IMP1CertificateError(message)


def _sha(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _digest(value: object, name: str) -> str:
    if (
        type(value) is not str
        or len(value) != 64
        or any(c not in "0123456789abcdef" for c in value)
    ):
        _fail(f"invalid {name} SHA-256")
    return value


def _exact(actual: object, expected: object, name: str) -> None:
    # Canonical byte comparison prevents bool/int and float/int aliases.
    if io.canonical_json_bytes(actual) != io.canonical_json_bytes(expected):
        _fail(f"{name} differs from the prospective contract")


def _pins(root: Path, paths: tuple[str, ...]) -> list[dict[str, str]]:
    return [
        {"path": path, "sha256": _sha(io.read_regular_file(root, path))}
        for path in paths
    ]


def _expected_config(root: Path) -> dict[str, object]:
    if QUALIFICATION_SHA256 is None:
        _fail("IMP1 synthetic qualification has not yet been bound")
    predecessor_config = io.read_regular_file(root, PREF1_CONFIG_PATH)
    if _sha(predecessor_config) != PREF1_CONFIG_SHA256:
        _fail("PREF1 configuration changed")
    environment = tomllib.loads(predecessor_config.decode("ascii"))["environment"]
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "selection_commit": PREF1_COMMIT,
        "selection_result_sha256": PREF1_SHA256,
        "evaluator_id": EVALUATOR_ID,
        "selected_candidate": "exact_accumulation_reconstruction_with_debit",
        "qualification_sha256": QUALIFICATION_SHA256,
        "order_squared_multiplier": 8,
        "channels": list(CHANNELS),
        "methods": list(METHODS),
        "reference_point_count": 9,
        "size_control_point_counts": [129, 2049],
        "limits": dict(LIMITS),
        "environment": environment,
        "scope": dict(SCOPE),
        "implementation": _pins(root, IMPLEMENTATION_PATHS),
        "inherited": _pins(root, INHERITED_PATHS),
    }


def emit_config_bytes(root: Path = ROOT) -> bytes:
    config = _expected_config(root)
    lines = []
    for key, value in config.items():
        if key in {"limits", "environment", "scope", "implementation", "inherited"}:
            continue
        lines.append(f"{key} = {json.dumps(value)}")
    for table in ("limits", "environment", "scope"):
        lines.append(f"\n[{table}]")
        for key, value in config[table].items():
            lines.append(f"{key} = {json.dumps(value)}")
    for table in ("implementation", "inherited"):
        for record in config[table]:
            lines.extend(
                (
                    f"\n[[{table}]]",
                    f"path = {json.dumps(record['path'])}",
                    f"sha256 = {json.dumps(record['sha256'])}",
                )
            )
    return ("\n".join(lines) + "\n").encode("ascii")


def validate_config(raw: bytes, root: Path = ROOT) -> dict[str, object]:
    try:
        config = tomllib.loads(raw.decode("ascii"))
    except (UnicodeError, tomllib.TOMLDecodeError) as error:
        raise IMP1CertificateError("malformed IMP1 config") from error
    _exact(config, _expected_config(root), "configuration/source pins")
    if _sha(io.read_regular_file(root, PREF1_PATH)) != PREF1_SHA256:
        _fail("independent numerical selection certificate changed")
    # The predecessor's ordinary compact verifier is itself raw/store/shadow/
    # Git-blind and authenticates its complete inherited source closure.
    from recursive_horizons.fgc.evolution.tdg11_msel1_pref1_binder import (
        verify_compact as verify_selection,
    )

    selected = verify_selection(root)
    if selected["selected_candidate"] != "exact_accumulation_reconstruction_with_debit":
        _fail("PREF1 does not license the implemented route")
    return config


def _validate_method_control(record: object, method: str) -> None:
    if type(record) is not dict:
        _fail("method control must be a mapping")
    required = {
        "method",
        "point_count",
        "owned_row_count",
        "channels",
        "shadow_stage_records",
        "next_family_shadow_stage_records",
        "total_shadow_proposals",
        "committed_proposals",
        "committed_stage_records",
        "independent_reference_channels",
        "endpoint_sha256",
        "assessment_sha256",
        "checkpoint_sha256",
        *CONTROL_FLAGS,
    }
    if set(record) != required:
        _fail("partial or broadened method control")
    for flag in CONTROL_FLAGS:
        if record[flag] is not True:
            _fail(f"synthetic method control did not prove {flag}")
    stages = 5 if method == METHODS[0] else 4
    for name, expected in (
        ("method", method),
        ("point_count", 9),
        ("owned_row_count", 4),
        ("channels", list(CHANNELS)),
        ("shadow_stage_records", 7 * stages),
        ("next_family_shadow_stage_records", 7 * stages),
        ("total_shadow_proposals", 14),
        ("committed_proposals", 4),
        ("committed_stage_records", 4 * stages),
        ("independent_reference_channels", 18),
    ):
        _exact(record[name], expected, name)
    for name in ("endpoint_sha256", "assessment_sha256", "checkpoint_sha256"):
        _digest(record[name], name)


def validate_qualification_report(report: object) -> None:
    if type(report) is not dict or set(report) != {
        "method_controls",
        "saturation_control",
        "size_controls",
    }:
        _fail("partial or broadened synthetic qualification report")
    controls = report["method_controls"]
    if type(controls) is not list or len(controls) != 2:
        _fail("both original integrators must be qualified")
    for record, method in zip(controls, METHODS, strict=True):
        _validate_method_control(record, method)
    saturation = report["saturation_control"]
    _exact(
        saturation,
        {
            "method": METHODS[0],
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
        "rounding saturation control",
    )
    sizes = report["size_controls"]
    expected_cases = [(method, count) for method in METHODS for count in (129, 2049)]
    if type(sizes) is not list or len(sizes) != len(expected_cases):
        _fail("synthetic size control coverage differs")
    for record, (method, count) in zip(sizes, expected_cases, strict=True):
        if type(record) is not dict or set(record) != {
            "method",
            "point_count",
            "owned_row_count",
            "channel_count",
            "admission_passed",
            "assessment_sha256",
            "shadow_stage_records",
            "physical_source_tested",
        }:
            _fail("synthetic size control is partial or broadened")
        for name, value in (
            ("method", method),
            ("point_count", count),
            ("owned_row_count", count - 5),
            ("channel_count", 18),
            ("admission_passed", True),
            ("shadow_stage_records", 35 if method == METHODS[0] else 28),
            ("physical_source_tested", False),
        ):
            _exact(record[name], value, name)
        _digest(record["assessment_sha256"], "size-control assessment")


def validate_result(
    config_raw: bytes, result_raw: bytes, root: Path = ROOT
) -> dict[str, object]:
    config = validate_config(config_raw, root)
    payload = result_raw[:-1] if result_raw.endswith(b"\n") else result_raw
    result = io.load_canonical_json(payload)
    if type(result) is not dict or set(result) != {
        "schema_version",
        "artifact_id",
        "classification",
        "evaluator_id",
        "config_sha256",
        "selection_commit",
        "selection_result_sha256",
        "implementation",
        "inherited",
        "limits",
        "scope",
        "environment",
        "qualification",
        "qualification_sha256",
    }:
        _fail("partial or broadened IMP1 result")
    for name, expected in (
        ("schema_version", 1),
        ("artifact_id", ARTIFACT_ID),
        ("classification", CLASSIFICATION),
        ("evaluator_id", EVALUATOR_ID),
        ("config_sha256", _sha(config_raw)),
        ("selection_commit", PREF1_COMMIT),
        ("selection_result_sha256", PREF1_SHA256),
        ("limits", LIMITS),
        ("scope", SCOPE),
        ("environment", config["environment"]),
        ("implementation", config["implementation"]),
        ("inherited", config["inherited"]),
    ):
        _exact(result[name], expected, name)
    validate_qualification_report(result["qualification"])
    _exact(
        result["qualification_sha256"], QUALIFICATION_SHA256, "qualification commitment"
    )
    _exact(
        _sha(io.canonical_json_bytes(result["qualification"])),
        QUALIFICATION_SHA256,
        "actual qualification report commitment",
    )
    return result


def verify_compact(root: Path = ROOT) -> dict[str, object]:
    """Read only tracked compact/config/source paths, never run a shadow."""

    return validate_result(
        io.read_regular_file(root, CONFIG_PATH),
        io.read_regular_file(root, RESULT_PATH),
        root,
    )


def qualify(
    config_raw: bytes,
    root: Path = ROOT,
    *,
    progress: Callable[[dict], None] | None = None,
) -> dict:
    config = validate_config(config_raw, root)
    from recursive_horizons.fgc.evolution.tdg11_msel1_authority import (
        observe_environment,
    )
    from recursive_horizons.fgc.evolution.tdg11_imp1_qualification import (
        run_qualification,
    )

    _exact(observe_environment(), config["environment"], "qualification environment")
    report = run_qualification(progress=progress)
    validate_qualification_report(report)
    _exact(
        _sha(io.canonical_json_bytes(report)),
        QUALIFICATION_SHA256,
        "independently bound synthetic report",
    )
    if io.read_regular_file(root, CONFIG_PATH) != config_raw:
        _fail("IMP1 config changed during synthetic qualification")
    validate_config(config_raw, root)
    _exact(
        observe_environment(), config["environment"], "environment after qualification"
    )
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "evaluator_id": EVALUATOR_ID,
        "config_sha256": _sha(config_raw),
        "selection_commit": PREF1_COMMIT,
        "selection_result_sha256": PREF1_SHA256,
        "implementation": config["implementation"],
        "inherited": config["inherited"],
        "limits": dict(LIMITS),
        "scope": dict(SCOPE),
        "environment": config["environment"],
        "qualification": report,
        "qualification_sha256": QUALIFICATION_SHA256,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--check", action="store_true", help="compact check (default)")
    action.add_argument(
        "--emit-config", action="store_true", help="print prospective config"
    )
    action.add_argument(
        "--qualify", action="store_true", help="explicit synthetic qualification"
    )
    parser.add_argument(
        "--write-result",
        action="store_true",
        help="exclusively publish the qualified compact result",
    )
    args = parser.parse_args(argv)
    if args.write_result and not args.qualify:
        parser.error("--write-result requires --qualify")
    if args.emit_config:
        sys.stdout.buffer.write(emit_config_bytes())
        return 0
    if args.qualify:
        if args.write_result and os.path.lexists(ROOT / RESULT_PATH):
            _fail(
                "IMP1 result already exists; no qualification or replacement was attempted"
            )
        config_raw = io.read_regular_file(ROOT, CONFIG_PATH)
        result = qualify(
            config_raw,
            progress=lambda event: print(
                json.dumps(event, sort_keys=True), file=sys.stderr, flush=True
            ),
        )
        raw = io.canonical_json_bytes(result) + b"\n"
        if not args.write_result:
            sys.stdout.buffer.write(raw)
            return 0
        validate_result(config_raw, raw)
        if io.read_regular_file(ROOT, CONFIG_PATH) != config_raw:
            _fail("IMP1 config changed before publication")
        io.publish_exclusive_file(ROOT, RESULT_PATH, raw)
        if io.read_regular_file(ROOT, CONFIG_PATH) != config_raw:
            _fail(
                "IMP1 config changed after publication; preserve the result as forensic evidence"
            )
        result = verify_compact()
    else:
        result = verify_compact()
    print(
        json.dumps(
            {
                "artifact_id": ARTIFACT_ID,
                "classification": result["classification"],
                "next_boundary": "separate_PRO20_PROTO19_HLT17_freeze",
                "production_state_advanced": False,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
