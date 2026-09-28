#!/usr/bin/env python3
"""Build, print, or verify the conditional FGC-1-EFT0-LED1 certificate."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys
import tomllib
from typing import Any, Mapping


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

DEFAULT_CONFIG = REPOSITORY / "configs" / "fgc" / "fgc-1-eft0-led1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results" / "fgc-1-eft0-led1.json"
OWNER_DOCUMENT = REPOSITORY / "docs" / "fgc-eft0-led1.md"
ACT1_CONFIG = REPOSITORY / "configs" / "fgc" / "fgc-1-action-gate.toml"
ACT1_RESULT = REPOSITORY / "results" / "fgc-1-action-gate.json"
COMP1_CONFIG = REPOSITORY / "configs" / "fgc" / "fgc-1-hyp1-con1-comp1.toml"
COMP1_RESULT = REPOSITORY / "results" / "fgc-1-hyp1-con1-comp1.json"

from recursive_horizons.fgc.eft_ledger import (  # noqa: E402
    EFT0_LED1_ARTIFACT_ID,
    EFTLedgerConfig,
    load_eft_ledger_config,
    local_eft_ledger_certificate,
)


IMPLEMENTATION_FILES = (
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "eft_ledger.py",
    Path(__file__).resolve(),
)


@dataclass(frozen=True, slots=True)
class EFT0LED1ReproductionConfig:
    """Validated config surface plus immutable provenance paths."""

    ledger: EFTLedgerConfig
    raw: dict[str, Any]
    config_path: Path
    source_sha256: str


def _mapping(name: str, value: object) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be an object/table")
    return value


def _reject_duplicate_json_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON object key: {key}")
        result[key] = value
    return result


def _canonical_json_text(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


def _file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _object_sha256(value: Any) -> str:
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")).hexdigest()


def _relative(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _contains_float(value: Any) -> bool:
    if isinstance(value, float):
        return True
    if isinstance(value, Mapping):
        return any(_contains_float(item) for item in value.values())
    if isinstance(value, list):
        return any(_contains_float(item) for item in value)
    return False


def _repository_file(name: str, path: Path) -> Path:
    resolved = path.resolve()
    try:
        resolved.relative_to(REPOSITORY)
    except ValueError as exc:
        raise ValueError(f"{name} must stay inside the repository") from exc
    if not resolved.is_file():
        raise ValueError(f"{name} must exist")
    return resolved


def _load_canonical_predecessor(
    path: Path, *, name: str, artifact_id: str
) -> dict[str, Any]:
    """Check one provenance predecessor without promoting its theorem."""

    source = _repository_file(f"{name} result", path).read_text(encoding="utf-8")
    try:
        payload = json.loads(source, object_pairs_hook=_reject_duplicate_json_pairs)
    except (json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"{name} predecessor must be valid unique-key JSON") from exc
    if not isinstance(payload, dict) or source != _canonical_json_text(payload):
        raise ValueError(f"{name} predecessor must use canonical JSON")
    if payload.get("artifact_id") != artifact_id:
        raise ValueError(f"{name} predecessor identity differs from {artifact_id}")
    return payload


def _bind_action_predecessor(ledger: EFTLedgerConfig) -> None:
    """Require the exact ledger mirror to agree with ACT1's frozen FGC-QR row."""

    action_path = _repository_file("ACT1 config", ACT1_CONFIG)
    raw = tomllib.loads(action_path.read_text(encoding="utf-8"))
    try:
        model = raw["models"]["FGC-QR"]
        mirrored = (
            Q(str(model["planck_mass"])),
            Q(str(model["scalar_mass"])),
            Q(str(model["quartic_coupling"])),
            Q(str(model["ricci_coupling"])),
            Q(str(model["quadratic_gb_coupling"])),
        )
    except (KeyError, TypeError, ValueError, ZeroDivisionError) as exc:
        raise ValueError("ACT1 FGC-QR action parameters are malformed") from exc
    declared = (
        ledger.planck_mass,
        ledger.scalar_mass,
        ledger.quartic_coupling,
        ledger.beta,
        ledger.eta,
    )
    if mirrored != declared:
        raise ValueError("EFT0 exact action ledger differs from ACT1 FGC-QR")
    action_result = _load_canonical_predecessor(
        ACT1_RESULT, name="ACT1", artifact_id="FGC-1-ACT1"
    )
    if (
        action_result.get("gate_status", {}).get(
            "declared_action_and_scalar_algebra_certificate_passed"
        )
        is not True
    ):
        raise ValueError("ACT1 declared-action predecessor gate is not passed")


def load_config(path: Path = DEFAULT_CONFIG) -> EFT0LED1ReproductionConfig:
    """Parse the strict core config and prove its source provenance is fixed."""

    config_path = _repository_file("EFT0 config", Path(path))
    ledger = load_eft_ledger_config(config_path)
    if ledger.artifact_id != EFT0_LED1_ARTIFACT_ID:
        raise ValueError("EFT0 artifact identity drifted")
    raw = tomllib.loads(config_path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict) or _contains_float(raw):
        raise ValueError("EFT0 config must be an exact non-floating TOML object")
    _bind_action_predecessor(ledger)
    _repository_file("COMP1 config", COMP1_CONFIG)
    _load_canonical_predecessor(
        COMP1_RESULT, name="COMP1", artifact_id="FGC-1-HYP1-CON1-COMP1"
    )
    return EFT0LED1ReproductionConfig(
        ledger=ledger,
        raw=raw,
        config_path=config_path,
        source_sha256=_file_sha256(config_path),
    )


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    """Load only unique-key, canonical, non-floating certificate JSON."""

    try:
        source = Path(path).read_text(encoding="utf-8")
        payload = json.loads(source, object_pairs_hook=_reject_duplicate_json_pairs)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError("EFT0 result must be valid unique-key JSON") from exc
    if not isinstance(payload, dict):
        raise ValueError("EFT0 result must be a JSON object")
    if source != _canonical_json_text(payload):
        raise ValueError("EFT0 result must use canonical sorted indented JSON encoding")
    if _contains_float(payload):
        raise ValueError("EFT0 result must not contain floating-point data")
    return payload


def _validate_conditional_certificate(certificate: Mapping[str, Any]) -> None:
    if certificate.get("artifact_id") != EFT0_LED1_ARTIFACT_ID:
        raise ValueError("EFT0 certificate identity drifted")
    if certificate.get("retained_eft_validity") is not False:
        raise ValueError("EFT0 must not promote retained EFT validity")
    frequency = _mapping("certificate.frequency_control", certificate.get("frequency_control"))
    if frequency.get("available") is not False:
        raise ValueError("EFT0 must not promote frequency control")
    box = _mapping("certificate.component_box", certificate.get("component_box"))
    if box.get("full_jet_box") is not False:
        raise ValueError("EFT0 must not promote a component box into a full jet box")
    operators = _mapping("certificate.operator_basis", certificate.get("operator_basis"))
    if (
        operators.get("scope") != "representative_only"
        or operators.get("complete_wilsonian_basis") is not False
        or operators.get("representative_omitted_remainders_locally_evaluated") is not False
    ):
        raise ValueError("EFT0 operator basis or omitted-remainder boundary drifted")
    controls = _mapping("certificate.local_controls", certificate.get("local_controls"))
    if not controls or any(_mapping(f"certificate.local_controls.{name}", item).get("passes") is not True for name, item in controls.items()):
        raise ValueError("EFT0 local controls are not all strictly passed")
    nonclaims = _mapping("certificate.nonclaims", certificate.get("nonclaims"))
    expected_nonclaims = {
        "uv_completion", "complete_wilsonian_basis", "wilson_coefficients_predicted",
        "radiative_stability", "frequency_control", "global_or_open_eft_validity",
        "initial_data", "evolution", "collapse", "defocusing", "singularity_resolution",
    }
    if set(nonclaims) != expected_nonclaims or any(value is not False for value in nonclaims.values()):
        raise ValueError("EFT0 nonclaim ledger has drifted")


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    """Produce the deterministic conditional local-ledger certificate."""

    config = load_config(config_path)
    certificate = local_eft_ledger_certificate(config.ledger)
    _validate_conditional_certificate(certificate)
    source_configs = {
        "eft0_ledger": _relative(config.config_path),
        "act1_action_predecessor": _relative(ACT1_CONFIG),
        "comp1_predecessor": _relative(COMP1_CONFIG),
    }
    source_results = {
        "act1_action_predecessor": _relative(ACT1_RESULT),
        "comp1_predecessor": _relative(COMP1_RESULT),
    }
    return {
        "schema_version": 1,
        "artifact_id": EFT0_LED1_ARTIFACT_ID,
        "project_version": config.ledger.project_version,
        "artifact_label": "FGC-1 conditional local retained-EFT declaration and component-box ledger",
        "classification": "conditional_declared_local_power_counting_ledger_not_an_eft_validity_certificate",
        "development_status": "exact_local_component_controls_passed_under_external_cutoff_and_representative_wilson_assumptions_global_eft_validity_open",
        "generated_by": _relative(Path(__file__)),
        "source_configs": source_configs,
        "source_config_sha256": {
            _relative(config.config_path): config.source_sha256,
            _relative(ACT1_CONFIG): _file_sha256(ACT1_CONFIG),
            _relative(COMP1_CONFIG): _file_sha256(COMP1_CONFIG),
        },
        "source_results": source_results,
        "source_results_sha256": {
            _relative(ACT1_RESULT): _file_sha256(ACT1_RESULT),
            _relative(COMP1_RESULT): _file_sha256(COMP1_RESULT),
        },
        "derivation_document": _relative(OWNER_DOCUMENT),
        "derivation_document_sha256": _file_sha256(OWNER_DOCUMENT),
        "implementation_sha256": {_relative(path): _file_sha256(path) for path in IMPLEMENTATION_FILES},
        "configuration": config.raw,
        "conditional_local_ledger_certificate": certificate,
        "conditional_local_ledger_certificate_sha256": _object_sha256(certificate),
        "gate_status": {
            "conditional_declared_local_component_controls_passed": True,
            "frequency_or_covariant_derivative_norm_control_passed": False,
            "representative_omitted_remainder_bounds_passed": False,
            "retained_eft_validity_envelope_passed": False,
            "global_or_open_eft_validity_passed": False,
            "initial_data_or_evolution_authorized": False,
        },
        "nonclaims": certificate["nonclaims"],
    }


def _check_result(path: Path, expected: dict[str, Any]) -> None:
    actual = load_canonical_result(path)
    if actual != expected:
        raise ValueError("EFT0 result differs from the deterministic recomputed certificate")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--mode", choices=("write", "check", "record"), default="write")
    arguments = parser.parse_args()
    payload = record(arguments.config)
    if arguments.mode == "record":
        print(_canonical_json_text(payload), end="")
    elif arguments.mode == "check":
        _check_result(arguments.output, payload)
        print(f"verified {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(_canonical_json_text(payload), encoding="utf-8")
        print(f"wrote {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
