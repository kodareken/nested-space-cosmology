"""Prospective candidate-blind contract for DEF1-STAB1 FRZ1.

This module aggregates the already-owned conversion instruments into one
deterministic, hashable contract.  It does not create the FRZ1 config, emit a
PREF1 result, read a trajectory, or set ``DEF1_error_map_passed``.  The later
freeze must bind these bytes and an independent PREF1 must re-execute the
controls before the conversion-instrument gate can be promoted.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from types import MappingProxyType
from typing import Mapping

from recursive_horizons.evidence_io import canonical_json_bytes

from .def1_geometry_error import INPUT_NAMES
from .def1_stab1 import DEF1_BOOLEAN_NAMES, ERROR_BUDGET_COMPONENTS
from .def1_stab1_qualification import (
    BASE_METRIC_TWO_JET_ORDER,
    CONVERSION_BASE_TO_ADM,
    CONVERSION_COVERAGE_MATRIX,
    CONVERSION_FO1_ADM_GEOMETRY,
    CONVERSION_IMP1_TO_Q,
    CONVERSION_RED1_MHG2,
    IMP1_18_CHANNEL_ORDER,
    PROVIDER_ROUTES,
    RED1_TO_MHG2_EQUATION_PAIRS,
    BaseToAdmTwoJetInverse,
    Def1Stab1QualificationError,
    QualificationCoverageMatrix,
    convert_base_to_adm_geometry_slots,
    convert_imp1_channels_to_q,
    conversion_identity_digest,
    execute_qualification_coverage_matrix,
    fo1_adm_geometry_slot_map,
    invert_base_metric_two_jets,
    owner_file_sha256,
    red1_mhg2_equation_identity_map,
)
from .spherical_reduction import Jet2


Q = Fraction
SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-DEF1-STAB1-FRZ1"
PREF1_ARTIFACT_ID = "FGC-1-DEF1-STAB1-PREF1"
PARTIAL_INSTRUMENT_ID = "FGC-1-DEF1-STAB1"

OWNER_PATHS = (
    "src/recursive_horizons/fgc/def1_stab1.py",
    "src/recursive_horizons/fgc/def1_geometry_error.py",
    "src/recursive_horizons/fgc/def1_stab1_providers.py",
    "src/recursive_horizons/fgc/def1_stab1_qualification.py",
)
CONVERSION_NAMES = (
    CONVERSION_RED1_MHG2,
    CONVERSION_FO1_ADM_GEOMETRY,
    CONVERSION_BASE_TO_ADM,
    CONVERSION_IMP1_TO_Q,
    CONVERSION_COVERAGE_MATRIX,
)
NONCLAIMS = (
    "FRZ1 contains no candidate or control trajectory values.",
    "FRZ1 does not establish a global nonlinear PDE error theorem.",
    "Richardson and IMP1 debits remain conditional premises.",
    "FRZ1 does not evaluate the nine DEF1 booleans or ROB1.",
    "FRZ1 does not set DEF1_error_map_passed or authorize holdout.",
)


def _rational_mapping(value: Fraction) -> dict[str, int]:
    return {"numerator": value.numerator, "denominator": value.denominator}


def _jet_mapping(jet: Jet2) -> dict[str, dict[str, int]]:
    return {
        name: _rational_mapping(getattr(jet, name))
        for name in ("value", "dt", "dr", "dtt", "dtr", "drr")
    }


def _base_inverse_control() -> BaseToAdmTwoJetInverse:
    alpha = Jet2(Q(2), Q(1, 3), -Q(1, 5), Q(1, 7), -Q(1, 11), Q(1, 13))
    shift = Jet2(Q(1, 4), -Q(1, 6), Q(1, 8), -Q(1, 10), Q(1, 12), -Q(1, 14))
    radial = Jet2(Q(3, 2), Q(1, 9), -Q(1, 15), Q(1, 17), Q(1, 19), -Q(1, 21))
    radius = Jet2(Q(4), -Q(1, 7), Q(5, 4), Q(1, 23), -Q(1, 25), Q(1, 27))
    radial_squared = radial * radial
    return invert_base_metric_two_jets(
        h_tt=-(alpha * alpha) + radial_squared * shift * shift,
        h_tr=radial_squared * shift,
        h_rr=radial_squared,
        areal_radius=radius,
        lapse_root=alpha.value,
        radial_scale_root=radial.value,
    )


def _canonical_contract_payload(
    *,
    owner_hashes: Mapping[str, str],
    conversion_hashes: Mapping[str, str],
    base_inverse: BaseToAdmTwoJetInverse,
    coverage: QualificationCoverageMatrix,
    imp1_additive_q: Fraction,
) -> dict[str, object]:
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "partial_instrument_id": PARTIAL_INSTRUMENT_ID,
        "pref1_artifact_id": PREF1_ARTIFACT_ID,
        "owner_hashes": dict(owner_hashes),
        "conversion_identity_hashes": dict(conversion_hashes),
        "red1_to_mhg2_pairs": [list(pair) for pair in RED1_TO_MHG2_EQUATION_PAIRS],
        "base_metric_two_jet_order": list(BASE_METRIC_TWO_JET_ORDER),
        "geometry_input_order": list(INPUT_NAMES),
        "imp1_channel_order": list(IMP1_18_CHANNEL_ORDER),
        "error_budget_component_order": list(ERROR_BUDGET_COMPONENTS),
        "provider_routes": [
            [route.component, route.route_id, route.conversion]
            for route in PROVIDER_ROUTES
        ],
        "def1_boolean_names": list(DEF1_BOOLEAN_NAMES),
        "base_inverse_control": {
            "h_tt": _jet_mapping(base_inverse.h_tt),
            "h_tr": _jet_mapping(base_inverse.h_tr),
            "h_rr": _jet_mapping(base_inverse.h_rr),
            "areal_radius": _jet_mapping(base_inverse.areal_radius),
            "alpha": _jet_mapping(base_inverse.alpha),
            "shift": _jet_mapping(base_inverse.shift),
            "lambda": _jet_mapping(base_inverse.lambda_jet),
            "roundtrip_holds": base_inverse.roundtrip_holds,
        },
        "coverage_control_count": len(coverage.routes),
        "coverage_component_count": len(coverage.components),
        "imp1_control_additive_q": {
            "numerator": imp1_additive_q.numerator,
            "denominator": imp1_additive_q.denominator,
        },
        "conversion_instrument_contract_complete": True,
        "candidate_or_control_trajectory_read": False,
        "global_pde_error_certified": False,
        "def1_booleans_evaluated": False,
        "def1_error_map_passed": False,
        "rob1_passed": False,
        "holdout_authorized": False,
        "physical_claimed": False,
        "nonclaims": list(NONCLAIMS),
    }


@dataclass(frozen=True, slots=True)
class Def1Stab1FreezeContract:
    """Re-derived prospective contract.  Cached fields cannot promote it."""

    payload: Mapping[str, object]
    payload_sha256: str
    conversion_instrument_contract_complete: bool
    def1_error_map_passed: bool
    candidate_or_control_trajectory_read: bool
    global_pde_error_certified: bool
    def1_booleans_evaluated: bool
    rob1_passed: bool
    holdout_authorized: bool
    physical_claimed: bool

    def __post_init__(self) -> None:
        expected = _derive_freeze_contract()
        if dict(self.payload) != expected["payload"]:
            raise Def1Stab1QualificationError(
                "DEF1-STAB1 freeze payload differs from live conversion owners"
            )
        if self.payload_sha256 != expected["payload_sha256"]:
            raise Def1Stab1QualificationError(
                "DEF1-STAB1 freeze payload hash differs"
            )
        for name in (
            "conversion_instrument_contract_complete",
            "def1_error_map_passed",
            "candidate_or_control_trajectory_read",
            "global_pde_error_certified",
            "def1_booleans_evaluated",
            "rob1_passed",
            "holdout_authorized",
            "physical_claimed",
        ):
            if getattr(self, name) is not expected[name]:
                raise Def1Stab1QualificationError(
                    f"DEF1-STAB1 freeze flag {name} differs"
                )
        object.__setattr__(self, "payload", MappingProxyType(dict(self.payload)))


def _derive_freeze_contract() -> dict[str, object]:
    red1_mhg2_equation_identity_map()
    fo1_adm_geometry_slot_map()
    base_inverse = _base_inverse_control()
    packed = convert_base_to_adm_geometry_slots(
        base_inverse,
        tangent=(("k.t", Q(5, 2)), ("k.r", -Q(1, 4))),
    )
    if not packed.complete or tuple(name for name, _ in packed.supplied) != INPUT_NAMES:
        raise Def1Stab1QualificationError(
            "BASE-to-ADM control did not produce all 26 geometry inputs"
        )
    channel_debits = tuple(
        (name, Q(index + 1, 1024))
        for index, name in enumerate(IMP1_18_CHANNEL_ORDER)
    )
    lipschitz = tuple((name, Q(1, 32)) for name in IMP1_18_CHANNEL_ORDER)
    imp1 = convert_imp1_channels_to_q(
        channel_debits=channel_debits,
        lipschitz=lipschitz,
        context="DEF1-STAB1 prospective conversion control",
        unit="declared complete-Q control units",
    )
    coverage = execute_qualification_coverage_matrix()
    owner_hashes = {path: owner_file_sha256(path) for path in OWNER_PATHS}
    conversion_hashes = {
        name: conversion_identity_digest(name) for name in CONVERSION_NAMES
    }
    payload = _canonical_contract_payload(
        owner_hashes=owner_hashes,
        conversion_hashes=conversion_hashes,
        base_inverse=base_inverse,
        coverage=coverage,
        imp1_additive_q=imp1.additive_q,
    )
    payload_sha256 = sha256(canonical_json_bytes(payload)).hexdigest()
    return {
        "payload": payload,
        "payload_sha256": payload_sha256,
        "conversion_instrument_contract_complete": True,
        "def1_error_map_passed": False,
        "candidate_or_control_trajectory_read": False,
        "global_pde_error_certified": False,
        "def1_booleans_evaluated": False,
        "rob1_passed": False,
        "holdout_authorized": False,
        "physical_claimed": False,
    }


def build_def1_stab1_freeze_contract() -> Def1Stab1FreezeContract:
    return Def1Stab1FreezeContract(**_derive_freeze_contract())  # type: ignore[arg-type]


__all__ = [
    "ARTIFACT_ID",
    "Def1Stab1FreezeContract",
    "NONCLAIMS",
    "OWNER_PATHS",
    "PREF1_ARTIFACT_ID",
    "build_def1_stab1_freeze_contract",
]
