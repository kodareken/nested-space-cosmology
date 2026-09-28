#!/usr/bin/env python3
"""Run one retry-3 SSPRK3-on-SBP4 u:R row-hash and envelope-owner diagnostic."""

from __future__ import annotations

import argparse
import ctypes
import errno
from fractions import Fraction
from hashlib import sha256
import json
from math import isfinite, nextafter
import os
from pathlib import Path
import secrets
import stat
import sys
from typing import Mapping, NoReturn, Sequence

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from scripts import run_fgc_tdg9_loc1 as loc1  # noqa: E402
from scripts import run_fgc_tdg9_ti2 as ti2_runner  # noqa: E402
from recursive_horizons.fgc.evolution import tdg5_stage_complete_refinement_runtime as tdg5  # noqa: E402
from recursive_horizons.fgc.evolution import tdg6_temporal_admission_runtime as tdg6  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_ar1_authority as ar1  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_ur1_authority as authority  # noqa: E402
from recursive_horizons.fgc.evolution.hlt16_campaign_store import (  # noqa: E402
    HLT16CampaignStore,
)
from recursive_horizons.fgc.evolution.numerical_engine import COMPARATOR_METHOD  # noqa: E402
from recursive_horizons.fgc.evolution.proto19_gr0_static_factory import (  # noqa: E402
    build_static_gr0_shells,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (  # noqa: E402
    CertifiedMagnitudeInterval,
    classify_tdg6_channel,
)


RUNNER_ID = "FGC-1-TDG9-UR1-RUN1"
RAW_SCHEMA = "UR1-raw-v1"
TERMINAL_CLASSES = frozenset(
    {
        "completed_u_R_rows_match_TI2_and_zero_lowers_owned_by_named_envelope_quantity",
        "completed_u_R_rows_match_TI2_but_D01_D12_lower_owners_differ",
        "row_stream_identity_mismatch",
        "replayed_production_intervals_differ_from_sealed_AC1",
        "shadow_proposal_premise_stop",
        "invalid_provenance_or_implementation",
    }
)
PUBLISHABLE_TERMINAL_CLASSES = TERMINAL_CLASSES - {
    "invalid_provenance_or_implementation"
}
COMPLETED_TERMINAL_CLASSES = frozenset(
    {
        "completed_u_R_rows_match_TI2_and_zero_lowers_owned_by_named_envelope_quantity",
        "completed_u_R_rows_match_TI2_but_D01_D12_lower_owners_differ",
        "row_stream_identity_mismatch",
        "replayed_production_intervals_differ_from_sealed_AC1",
    }
)
OWNER_COMPLETED_TERMINAL_CLASSES = frozenset(
    {
        "completed_u_R_rows_match_TI2_and_zero_lowers_owned_by_named_envelope_quantity",
        "completed_u_R_rows_match_TI2_but_D01_D12_lower_owners_differ",
    }
)
_OUTPUT_LEAF_LIMIT = 512 * 1024 * 1024
_HEX = frozenset("0123456789abcdef")
_EXACT_KEYS = frozenset({"binary64_hex", "numerator", "denominator"})
_RATIONAL_KEYS = frozenset({"numerator", "denominator"})
_ENVELOPE_KEYS = frozenset(
    {
        "ambiguous_discriminant_count",
        "bernstein_certification_slack",
        "certified_continuous_upper_bound",
        "coefficient_construction_debit",
        "endpoint_candidate_count",
        "outward_arithmetic_debit",
        "polynomial_count",
        "raw_candidate_maximum",
        "real_interior_root_count",
    }
)
_OWNER_KEYS = frozenset(
    {
        "bernstein_certification_slack",
        "bernstein_certification_slack_not_part_of_lower_clip",
        "coefficient_construction_debit",
        "exact_clipped",
        "exact_unclipped",
        "outward_arithmetic_debit",
        "owner",
        "raw_candidate_maximum",
        "reproduced_stored_lower",
        "source",
    }
)
_INTERVAL_KEYS = frozenset(
    {
        "envelope",
        "lower_bound",
        "owned_row_count",
        "subinterval_count",
        "upper_bound",
        "zero_lower_owner",
    }
)
_U_R_KEYS = frozenset(
    {
        "D01",
        "D01_D12_lower_owners_equal",
        "D12",
        "admission_passed",
        "channel",
        "classification",
        "matches_sealed_AC1_production_intervals",
        "order_threshold_passed",
        "order_threshold_resolved",
        "temporal_retry_permitted",
    }
)
_ROW_BASE_KEYS = frozenset(
    {
        "algorithm",
        "domain",
        "expected_sha256",
        "matched",
        "observed_sha256",
        "row_count",
    }
)
_ROW_MATCH_KEYS = _ROW_BASE_KEYS | {
    "related_TI2_retry3_u_R_radius_free_complete_C_class",
    "related_not_replacement_production_admission",
}


class UR1RunnerError(RuntimeError):
    """Typed fail-closed UR1 runner error."""

    def __init__(self, owner: str, code: str, detail: object) -> None:
        self.owner = str(owner)
        self.code = str(code)
        self.detail = " ".join(str(detail).replace(str(ROOT), "<repo>").split())[:640]
        super().__init__(f"{self.owner}/{self.code}: {self.detail}")


def _fail(owner: str, code: str, detail: object) -> NoReturn:
    raise UR1RunnerError(owner, code, detail)


def _pretty(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def _unique(items: list[tuple[str, object]]) -> dict[str, object]:
    answer: dict[str, object] = {}
    for key, value in items:
        if key in answer:
            _fail("serialization", "duplicate_JSON_key", key)
        answer[key] = value
    return answer


def _exact_binary64(value: object, *, label: str) -> dict[str, str]:
    if type(value) is not float:
        _fail("serialization", "nonexact_binary64", label)
    if not isfinite(value):
        _fail("serialization", "nonfinite_binary64", label)
    fraction = Fraction.from_float(value)
    encoded = {
        "binary64_hex": value.hex(),
        "numerator": str(fraction.numerator),
        "denominator": str(fraction.denominator),
    }
    rebuilt = float.fromhex(encoded["binary64_hex"])
    if rebuilt.hex() != value.hex() or Fraction.from_float(rebuilt) != fraction:
        _fail("serialization", "binary64_roundtrip", label)
    return encoded


def _exact_rational(value: Fraction, *, label: str) -> dict[str, str]:
    if not isinstance(value, Fraction) or isinstance(value, bool):
        _fail("serialization", "nonexact_rational", label)
    encoded = {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
    }
    rebuilt = Fraction(int(encoded["numerator"]), int(encoded["denominator"]))
    if (
        rebuilt != value
        or str(rebuilt.numerator) != encoded["numerator"]
        or str(rebuilt.denominator) != encoded["denominator"]
    ):
        _fail("serialization", "rational_roundtrip", label)
    return encoded


def _fraction_from_exact(value: object, *, label: str) -> Fraction:
    if not isinstance(value, Mapping):
        _fail("serialization", "exact_mapping", label)
    if set(value) != _EXACT_KEYS:
        _fail("serialization", "exact_keys", label)
    hex_value = value["binary64_hex"]
    numerator = value["numerator"]
    denominator = value["denominator"]
    if (
        not isinstance(hex_value, str)
        or not isinstance(numerator, str)
        or not isinstance(denominator, str)
    ):
        _fail("serialization", "exact_types", label)
    try:
        rebuilt = float.fromhex(hex_value)
        fraction = Fraction(int(numerator), int(denominator))
    except (ValueError, ZeroDivisionError) as exc:
        raise UR1RunnerError("serialization", "exact_parse", label) from exc
    if (
        not isfinite(rebuilt)
        or rebuilt.hex() != hex_value
        or str(fraction.numerator) != numerator
        or str(fraction.denominator) != denominator
        or Fraction.from_float(rebuilt) != fraction
    ):
        _fail("serialization", "exact_roundtrip", label)
    return fraction


def _fraction_from_rational(value: object, *, label: str) -> Fraction:
    if not isinstance(value, Mapping) or set(value) != _RATIONAL_KEYS:
        _fail("serialization", "rational_keys", label)
    numerator = value["numerator"]
    denominator = value["denominator"]
    if not isinstance(numerator, str) or not isinstance(denominator, str):
        _fail("serialization", "rational_types", label)
    try:
        fraction = Fraction(int(numerator), int(denominator))
    except (ValueError, ZeroDivisionError) as exc:
        raise UR1RunnerError("serialization", "rational_parse", label) from exc
    if (
        str(fraction.numerator) != numerator
        or str(fraction.denominator) != denominator
    ):
        _fail("serialization", "rational_roundtrip", label)
    return fraction


def downward_binary64_lower(exact: Fraction) -> float:
    """Reproduce TDG6 stored-lower downward rounding without calling it."""

    if not isinstance(exact, Fraction) or exact < 0:
        _fail("envelope", "nonnegative_exact_lower", exact)
    lower = float(exact)
    if not isfinite(lower):
        _fail("envelope", "nonfinite_stored_lower", lower)
    if Fraction.from_float(lower) > exact:
        lower = nextafter(lower, 0.0)
    if not isfinite(lower) or lower < 0.0:
        _fail("envelope", "rounded_lower_invalid", lower)
    return lower


def name_zero_lower_owner(
    *,
    raw: Fraction,
    unclipped: Fraction,
    reproduced: float,
) -> str:
    """Name the exact term that clips a stored lower bound to zero."""

    if reproduced != 0.0:
        _fail("envelope", "stored_lower_not_zero", reproduced)
    if raw == 0:
        return "raw_candidate_maximum_is_zero"
    if unclipped <= 0:
        return "coefficient_plus_arithmetic_debit_clips_positive_raw_maximum"
    if unclipped > 0 and reproduced == 0.0:
        return "downward_binary64_rounding_of_positive_exact_lower"
    _fail("envelope", "zero_lower_owner_unclassified", (str(raw), str(unclipped)))


def classify_zero_lower_owner(
    envelope: tdg5.Binary64CubicEnvelope,
    stored_lower: float,
    *,
    label: str,
) -> dict[str, object]:
    if not isinstance(envelope, tdg5.Binary64CubicEnvelope):
        _fail("envelope", "type", type(envelope).__name__)
    if type(stored_lower) is not float or not isfinite(stored_lower) or stored_lower < 0.0:
        _fail("envelope", "stored_lower", label)
    raw = envelope.raw_candidate_maximum
    construction = envelope.coefficient_construction_debit
    arithmetic = envelope.outward_arithmetic_debit
    slack = envelope.bernstein_certification_slack
    for name, value in (
        ("raw_candidate_maximum", raw),
        ("coefficient_construction_debit", construction),
        ("outward_arithmetic_debit", arithmetic),
        ("bernstein_certification_slack", slack),
    ):
        if type(value) is not float or not isfinite(value) or value < 0.0:
            _fail("envelope", "nonnegative_binary64", f"{label}.{name}")
    exact_unclipped = (
        Fraction.from_float(raw)
        - Fraction.from_float(construction)
        - Fraction.from_float(arithmetic)
    )
    exact_clipped = max(Fraction(0), exact_unclipped)
    reproduced = downward_binary64_lower(exact_clipped)
    if reproduced != stored_lower:
        _fail("envelope", "reproduced_lower_mismatch", label)
    owner = name_zero_lower_owner(
        raw=Fraction.from_float(raw),
        unclipped=exact_unclipped,
        reproduced=reproduced,
    )
    if owner == "raw_candidate_maximum_is_zero":
        source = {
            "term": "raw_candidate_maximum",
            "value": _exact_binary64(raw, label=f"{label}.source.raw"),
        }
    elif owner == "coefficient_plus_arithmetic_debit_clips_positive_raw_maximum":
        source = {
            "term": "coefficient_construction_debit_plus_outward_arithmetic_debit",
            "raw_candidate_maximum": _exact_binary64(
                raw, label=f"{label}.source.raw"
            ),
            "coefficient_construction_debit": _exact_binary64(
                construction, label=f"{label}.source.construction"
            ),
            "outward_arithmetic_debit": _exact_binary64(
                arithmetic, label=f"{label}.source.arithmetic"
            ),
        }
    else:
        source = {
            "term": "positive_exact_lower_before_downward_binary64_rounding",
            "exact_unclipped": _exact_rational(
                exact_unclipped, label=f"{label}.source.unclipped"
            ),
        }
    return {
        "owner": owner,
        "raw_candidate_maximum": _exact_binary64(raw, label=f"{label}.raw"),
        "coefficient_construction_debit": _exact_binary64(
            construction, label=f"{label}.construction"
        ),
        "outward_arithmetic_debit": _exact_binary64(
            arithmetic, label=f"{label}.arithmetic"
        ),
        "bernstein_certification_slack": _exact_binary64(
            slack, label=f"{label}.bernstein"
        ),
        "bernstein_certification_slack_not_part_of_lower_clip": True,
        "exact_unclipped": _exact_rational(
            exact_unclipped, label=f"{label}.unclipped"
        ),
        "exact_clipped": _exact_rational(exact_clipped, label=f"{label}.clipped"),
        "reproduced_stored_lower": _exact_binary64(
            reproduced, label=f"{label}.reproduced"
        ),
        "source": source,
    }


def serialize_envelope(
    envelope: tdg5.Binary64CubicEnvelope, *, label: str
) -> dict[str, object]:
    if not isinstance(envelope, tdg5.Binary64CubicEnvelope):
        _fail("envelope", "type", type(envelope).__name__)
    return {
        "polynomial_count": envelope.polynomial_count,
        "endpoint_candidate_count": envelope.endpoint_candidate_count,
        "real_interior_root_count": envelope.real_interior_root_count,
        "ambiguous_discriminant_count": envelope.ambiguous_discriminant_count,
        "raw_candidate_maximum": _exact_binary64(
            envelope.raw_candidate_maximum, label=f"{label}.raw"
        ),
        "coefficient_construction_debit": _exact_binary64(
            envelope.coefficient_construction_debit, label=f"{label}.construction"
        ),
        "outward_arithmetic_debit": _exact_binary64(
            envelope.outward_arithmetic_debit, label=f"{label}.arithmetic"
        ),
        "bernstein_certification_slack": _exact_binary64(
            envelope.bernstein_certification_slack, label=f"{label}.bernstein"
        ),
        "certified_continuous_upper_bound": _exact_binary64(
            envelope.certified_continuous_upper_bound, label=f"{label}.certified_upper"
        ),
    }


def serialize_magnitude_interval(
    interval: tdg6.TDG6Binary64MagnitudeInterval, *, label: str
) -> dict[str, object]:
    if not isinstance(interval, tdg6.TDG6Binary64MagnitudeInterval):
        _fail("interval", "type", type(interval).__name__)
    if interval.upper_bound != interval.envelope.certified_continuous_upper_bound:
        _fail("interval", "upper_identity", label)
    owner = classify_zero_lower_owner(
        interval.envelope, interval.lower_bound, label=f"{label}.lower"
    )
    serialized = {
        "lower_bound": _exact_binary64(interval.lower_bound, label=f"{label}.lower"),
        "upper_bound": _exact_binary64(interval.upper_bound, label=f"{label}.upper"),
        "subinterval_count": interval.subinterval_count,
        "owned_row_count": interval.owned_row_count,
        "envelope": serialize_envelope(interval.envelope, label=f"{label}.envelope"),
        "zero_lower_owner": owner,
    }
    validate_magnitude_interval(serialized, label=label)
    return serialized


def validate_magnitude_interval(value: object, *, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != _INTERVAL_KEYS:
        _fail("status", "interval_keys", label)
    envelope = value.get("envelope")
    if not isinstance(envelope, Mapping) or set(envelope) != _ENVELOPE_KEYS:
        _fail("status", "envelope_keys", label)
    for name in (
        "polynomial_count",
        "endpoint_candidate_count",
        "real_interior_root_count",
        "ambiguous_discriminant_count",
        "subinterval_count",
        "owned_row_count",
    ):
        observed = envelope[name] if name in envelope else value[name]
        if type(observed) is not int or isinstance(observed, bool) or observed < 0:
            _fail("status", "count_type", f"{label}.{name}")
    if envelope["endpoint_candidate_count"] != 2 * envelope["polynomial_count"]:
        _fail("status", "endpoint_candidate_identity", label)
    if envelope["polynomial_count"] != value["subinterval_count"] * value["owned_row_count"]:
        _fail("status", "polynomial_count_identity", label)
    lower = _fraction_from_exact(value["lower_bound"], label=f"{label}.lower")
    upper = _fraction_from_exact(value["upper_bound"], label=f"{label}.upper")
    raw_envelope = _fraction_from_exact(
        envelope["raw_candidate_maximum"], label=f"{label}.envelope.raw"
    )
    construction_envelope = _fraction_from_exact(
        envelope["coefficient_construction_debit"],
        label=f"{label}.envelope.construction",
    )
    arithmetic_envelope = _fraction_from_exact(
        envelope["outward_arithmetic_debit"],
        label=f"{label}.envelope.arithmetic",
    )
    bernstein_envelope = _fraction_from_exact(
        envelope["bernstein_certification_slack"],
        label=f"{label}.envelope.bernstein",
    )
    certified = _fraction_from_exact(
        envelope["certified_continuous_upper_bound"], label=f"{label}.certified_upper"
    )
    if any(
        item < 0
        for item in (
            lower,
            upper,
            raw_envelope,
            construction_envelope,
            arithmetic_envelope,
            bernstein_envelope,
            certified,
        )
    ):
        _fail("status", "nonnegative_interval_identity", label)
    if upper != certified or upper < lower or certified < raw_envelope:
        _fail("status", "certified_upper_identity", label)
    owner = value.get("zero_lower_owner")
    if not isinstance(owner, Mapping) or set(owner) != _OWNER_KEYS:
        _fail("status", "owner_keys", label)
    if owner.get("owner") not in authority.ZERO_LOWER_OWNERS:
        _fail("status", "owner_enum", owner.get("owner"))
    if owner.get("bernstein_certification_slack_not_part_of_lower_clip") is not True:
        _fail("status", "bernstein_not_clip_flag", label)
    raw = _fraction_from_exact(
        owner["raw_candidate_maximum"], label=f"{label}.owner.raw"
    )
    construction = _fraction_from_exact(
        owner["coefficient_construction_debit"], label=f"{label}.owner.construction"
    )
    arithmetic = _fraction_from_exact(
        owner["outward_arithmetic_debit"], label=f"{label}.owner.arithmetic"
    )
    bernstein = _fraction_from_exact(
        owner["bernstein_certification_slack"], label=f"{label}.owner.bernstein"
    )
    if (
        raw != raw_envelope
        or construction != construction_envelope
        or arithmetic != arithmetic_envelope
        or bernstein != bernstein_envelope
    ):
        _fail("status", "owner_envelope_identity", label)
    unclipped = _fraction_from_rational(
        owner["exact_unclipped"], label=f"{label}.unclipped"
    )
    clipped = _fraction_from_rational(owner["exact_clipped"], label=f"{label}.clipped")
    reproduced = _fraction_from_exact(
        owner["reproduced_stored_lower"], label=f"{label}.reproduced"
    )
    expected_unclipped = raw - construction - arithmetic
    expected_clipped = max(Fraction(0), expected_unclipped)
    expected_reproduced = Fraction.from_float(downward_binary64_lower(expected_clipped))
    if (
        unclipped != expected_unclipped
        or clipped != expected_clipped
        or reproduced != expected_reproduced
        or reproduced != lower
    ):
        _fail("status", "owner_recompute", label)
    expected_owner = name_zero_lower_owner(
        raw=raw,
        unclipped=unclipped,
        reproduced=float.fromhex(owner["reproduced_stored_lower"]["binary64_hex"]),
    )
    if owner["owner"] != expected_owner:
        _fail("status", "owner_identity", label)
    source = owner.get("source")
    if not isinstance(source, Mapping):
        _fail("status", "owner_source", label)
    if expected_owner == "raw_candidate_maximum_is_zero":
        if set(source) != {"term", "value"} or source.get("term") != (
            "raw_candidate_maximum"
        ):
            _fail("status", "owner_source", label)
        if _fraction_from_exact(
            source["value"], label=f"{label}.source.raw"
        ) != raw:
            _fail("status", "owner_source_identity", label)
    elif expected_owner == (
        "coefficient_plus_arithmetic_debit_clips_positive_raw_maximum"
    ):
        if set(source) != {
            "term",
            "raw_candidate_maximum",
            "coefficient_construction_debit",
            "outward_arithmetic_debit",
        } or source.get("term") != (
            "coefficient_construction_debit_plus_outward_arithmetic_debit"
        ):
            _fail("status", "owner_source", label)
        if (
            _fraction_from_exact(
                source["raw_candidate_maximum"], label=f"{label}.source.raw"
            )
            != raw
            or _fraction_from_exact(
                source["coefficient_construction_debit"],
                label=f"{label}.source.construction",
            )
            != construction
            or _fraction_from_exact(
                source["outward_arithmetic_debit"],
                label=f"{label}.source.arithmetic",
            )
            != arithmetic
        ):
            _fail("status", "owner_source_identity", label)
    else:
        if set(source) != {"term", "exact_unclipped"} or source.get("term") != (
            "positive_exact_lower_before_downward_binary64_rounding"
        ):
            _fail("status", "owner_source", label)
        if _fraction_from_rational(
            source["exact_unclipped"], label=f"{label}.source.unclipped"
        ) != unclipped:
            _fail("status", "owner_source_identity", label)
    return value


def locate_u_R(
    admission: tdg6.TDG6ContinuousAdmissionEvidence,
) -> tdg6.TDG6RuntimeChannelAdmission:
    if not isinstance(admission, tdg6.TDG6ContinuousAdmissionEvidence):
        _fail("admission", "type", type(admission).__name__)
    matches = [
        item
        for item in admission.channel_admissions
        if item.channel == authority.PUBLISHED_CHANNEL
    ]
    if len(matches) != 1:
        _fail("admission", "u_R_count", len(matches))
    return matches[0]


def serialize_u_R(
    admission: tdg6.TDG6ContinuousAdmissionEvidence,
) -> dict[str, object]:
    item = locate_u_R(admission)
    outer = serialize_magnitude_interval(item.outer_difference, label="D01")
    finest = serialize_magnitude_interval(item.finest_difference, label="D12")
    decision = classify_tdg6_channel(
        item.outer_difference.exact_interval(),
        item.finest_difference.exact_interval(),
    )
    if (
        decision.classification != item.classification
        or decision.admission_passed is not item.admission_passed
        or decision.temporal_retry_permitted is not item.temporal_retry_permitted
        or decision.order_threshold_resolved is not item.order_threshold_resolved
        or decision.order_threshold_passed is not item.order_threshold_passed
        or decision.finest_pair_debit != item.finest_difference.exact_interval().upper
        or item.classification != "order_inconclusive"
        or item.admission_passed is not False
    ):
        _fail("admission", "reclassification", item.channel)
    matches_ac1 = (
        item.outer_difference.lower_bound.hex() == authority.AC1_U_R_D01_LOWER_HEX
        and item.outer_difference.upper_bound.hex() == authority.AC1_U_R_D01_UPPER_HEX
        and item.finest_difference.lower_bound.hex() == authority.AC1_U_R_D12_LOWER_HEX
        and item.finest_difference.upper_bound.hex() == authority.AC1_U_R_D12_UPPER_HEX
        and item.classification == authority.AC1_U_R_CLASS
        and item.outer_difference.exact_interval()
        == CertifiedMagnitudeInterval(
            Fraction.from_float(float.fromhex(authority.AC1_U_R_D01_LOWER_HEX)),
            Fraction.from_float(float.fromhex(authority.AC1_U_R_D01_UPPER_HEX)),
        )
        and item.finest_difference.exact_interval()
        == CertifiedMagnitudeInterval(
            Fraction.from_float(float.fromhex(authority.AC1_U_R_D12_LOWER_HEX)),
            Fraction.from_float(float.fromhex(authority.AC1_U_R_D12_UPPER_HEX)),
        )
    )
    d01_owner = str(outer["zero_lower_owner"]["owner"])
    d12_owner = str(finest["zero_lower_owner"]["owner"])
    serialized = {
        "channel": authority.PUBLISHED_CHANNEL,
        "classification": item.classification,
        "admission_passed": item.admission_passed,
        "temporal_retry_permitted": item.temporal_retry_permitted,
        "order_threshold_resolved": item.order_threshold_resolved,
        "order_threshold_passed": item.order_threshold_passed,
        "D01": outer,
        "D12": finest,
        "matches_sealed_AC1_production_intervals": matches_ac1,
        "D01_D12_lower_owners_equal": d01_owner == d12_owner,
    }
    validate_u_R(serialized)
    return serialized


def validate_u_R(value: object) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != _U_R_KEYS:
        _fail("status", "u_R_keys", value)
    if value.get("channel") != authority.PUBLISHED_CHANNEL:
        _fail("status", "published_channel", value.get("channel"))
    validate_magnitude_interval(value.get("D01"), label="D01")
    validate_magnitude_interval(value.get("D12"), label="D12")
    outer = CertifiedMagnitudeInterval(
        _fraction_from_exact(value["D01"]["lower_bound"], label="D01.lower"),
        _fraction_from_exact(value["D01"]["upper_bound"], label="D01.upper"),
    )
    finest = CertifiedMagnitudeInterval(
        _fraction_from_exact(value["D12"]["lower_bound"], label="D12.lower"),
        _fraction_from_exact(value["D12"]["upper_bound"], label="D12.upper"),
    )
    decision = classify_tdg6_channel(outer, finest)
    if (
        decision.classification != value.get("classification")
        or decision.admission_passed is not value.get("admission_passed")
        or decision.temporal_retry_permitted is not value.get("temporal_retry_permitted")
        or decision.order_threshold_resolved is not value.get("order_threshold_resolved")
        or decision.order_threshold_passed is not value.get("order_threshold_passed")
        or decision.classification != "order_inconclusive"
        or value.get("admission_passed") is not False
    ):
        _fail("status", "u_R_reclassification", value.get("classification"))
    expected_match = (
        value["D01"]["lower_bound"]["binary64_hex"] == authority.AC1_U_R_D01_LOWER_HEX
        and value["D01"]["upper_bound"]["binary64_hex"] == authority.AC1_U_R_D01_UPPER_HEX
        and value["D12"]["lower_bound"]["binary64_hex"] == authority.AC1_U_R_D12_LOWER_HEX
        and value["D12"]["upper_bound"]["binary64_hex"] == authority.AC1_U_R_D12_UPPER_HEX
        and value.get("classification") == authority.AC1_U_R_CLASS
    )
    if value.get("matches_sealed_AC1_production_intervals") is not expected_match:
        _fail("status", "ac1_interval_identity", expected_match)
    owners_equal = (
        value["D01"]["zero_lower_owner"]["owner"]
        == value["D12"]["zero_lower_owner"]["owner"]
    )
    if value.get("D01_D12_lower_owners_equal") is not owners_equal:
        _fail("status", "owner_equality_identity", owners_equal)
    return value


def hash_u_R_rows(prepared: tdg6.TDG6PreparedGR0Compositor) -> dict[str, object]:
    surface = loc1._surface(prepared)
    rows = list(loc1._rows(surface, authority.PUBLISHED_CHANNEL))
    if len(rows) != authority.OWNED_ROW_COUNT:
        _fail("rows", "row_count", len(rows))
    observed = loc1._row_hash(rows)
    matched = observed == authority.TI2_RETRY3_U_R_ROW_SHA256
    payload: dict[str, object] = {
        "algorithm": "sha256",
        "domain": authority.ROW_HASH_DOMAIN,
        "row_count": authority.OWNED_ROW_COUNT,
        "observed_sha256": observed,
        "expected_sha256": authority.TI2_RETRY3_U_R_ROW_SHA256,
        "matched": matched,
    }
    if matched:
        payload["related_TI2_retry3_u_R_radius_free_complete_C_class"] = (
            authority.TI2_RETRY3_U_R_RELATED_COMPLETE_C_CLASS
        )
        payload["related_not_replacement_production_admission"] = True
    validate_row_stream(payload)
    return payload


def validate_row_stream(value: object) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        _fail("status", "row_stream_type", type(value).__name__)
    matched = value.get("matched") is True
    expected_keys = _ROW_MATCH_KEYS if matched else _ROW_BASE_KEYS
    if set(value) != expected_keys:
        _fail("status", "row_stream_keys", sorted(value))
    if (
        value.get("algorithm") != "sha256"
        or value.get("domain") != authority.ROW_HASH_DOMAIN
        or value.get("row_count") != authority.OWNED_ROW_COUNT
        or value.get("expected_sha256") != authority.TI2_RETRY3_U_R_ROW_SHA256
        or not isinstance(value.get("observed_sha256"), str)
        or len(str(value.get("observed_sha256"))) != 64
        or any(character not in _HEX for character in str(value.get("observed_sha256")))
        or matched is not (value.get("observed_sha256") == authority.TI2_RETRY3_U_R_ROW_SHA256)
    ):
        _fail("status", "row_stream_identity", value.get("observed_sha256"))
    if matched:
        if (
            value.get("related_TI2_retry3_u_R_radius_free_complete_C_class")
            != authority.TI2_RETRY3_U_R_RELATED_COMPLETE_C_CLASS
            or value.get("related_not_replacement_production_admission") is not True
        ):
            _fail("status", "related_TI2_class", value)
    return value


def reduce_ur1_terminal(
    *,
    rows_matched: bool,
    intervals_match_ac1: bool,
    owners_equal: bool,
) -> str:
    if not rows_matched:
        return "row_stream_identity_mismatch"
    if not intervals_match_ac1:
        return "replayed_production_intervals_differ_from_sealed_AC1"
    if owners_equal:
        return (
            "completed_u_R_rows_match_TI2_and_zero_lowers_owned_by_named_envelope_quantity"
        )
    return "completed_u_R_rows_match_TI2_but_D01_D12_lower_owners_differ"


def _execution_authority(root: Path, commit: str) -> authority.UR1Authority:
    try:
        return authority.authorize(
            root,
            authority.read_leaf(root, authority.CONFIG_PATH),
            authority.read_leaf(root, authority.RESULT_PATH),
            commit,
        )
    except Exception as exc:
        raise UR1RunnerError("authority", "rejected", exc) from exc


def _snapshot_store(root: Path) -> tuple[int, str]:
    try:
        return ti2_runner._snapshot_store(root)
    except Exception as exc:
        raise UR1RunnerError("provenance", "store_snapshot", exc) from exc


def _open_output_parent(root: Path) -> tuple[int, str]:
    relative = Path(authority.OUTPUT_NAMESPACE)
    descriptor = os.open(
        root,
        os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        for part in relative.parent.parts:
            try:
                child = os.open(
                    part,
                    os.O_RDONLY
                    | getattr(os, "O_DIRECTORY", 0)
                    | getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=descriptor,
                )
            except FileNotFoundError:
                os.mkdir(part, mode=0o755, dir_fd=descriptor)
                child = os.open(
                    part,
                    os.O_RDONLY
                    | getattr(os, "O_DIRECTORY", 0)
                    | getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=descriptor,
                )
            os.close(descriptor)
            descriptor = child
        return descriptor, relative.name
    except Exception:
        os.close(descriptor)
        raise


def _publish(
    root: Path, manifest: Mapping[str, object], terminal: Mapping[str, object]
) -> None:
    parent_fd, target_name = _open_output_parent(root)
    stage_name = f"{authority.STAGING_PREFIX}{os.getpid()}-{secrets.token_hex(8)}"
    stage_fd: int | None = None
    try:
        try:
            os.stat(target_name, dir_fd=parent_fd, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            _fail("output", "namespace_exists", authority.OUTPUT_NAMESPACE)
        os.mkdir(stage_name, mode=0o700, dir_fd=parent_fd)
        stage_fd = os.open(
            stage_name,
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=parent_fd,
        )
        for name, value in (("manifest.json", manifest), ("terminal.json", terminal)):
            leaf = os.open(
                name,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
                0o644,
                dir_fd=stage_fd,
            )
            try:
                pending = memoryview(_pretty(value))
                while pending:
                    written = os.write(leaf, pending)
                    if written <= 0:
                        _fail("output", "short_write", name)
                    pending = pending[written:]
                os.fsync(leaf)
            finally:
                os.close(leaf)
        os.fsync(stage_fd)
        if os.uname().sysname != "Darwin":
            _fail("output", "exclusive_adoption_unavailable", os.uname().sysname)
        rename = getattr(ctypes.CDLL(None, use_errno=True), "renameatx_np", None)
        if rename is None:
            _fail("output", "exclusive_adoption_unavailable", "renameatx_np")
        rename.argtypes = (
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint,
        )
        rename.restype = ctypes.c_int
        if rename(
            parent_fd,
            os.fsencode(stage_name),
            parent_fd,
            os.fsencode(target_name),
            0x00000004,
        ):
            error = ctypes.get_errno()
            if error in {errno.EEXIST, errno.ENOTEMPTY}:
                _fail("output", "namespace_arrived", authority.OUTPUT_NAMESPACE)
            _fail("output", "exclusive_adoption_failed", error)
        os.close(stage_fd)
        stage_fd = None
        os.fsync(parent_fd)
    except Exception:
        if stage_fd is not None:
            for name in ("manifest.json", "terminal.json"):
                try:
                    os.unlink(name, dir_fd=stage_fd)
                except FileNotFoundError:
                    pass
            os.close(stage_fd)
        try:
            os.rmdir(stage_name, dir_fd=parent_fd)
        except FileNotFoundError:
            pass
        raise
    finally:
        os.close(parent_fd)


def _publish_terminal(
    root: Path,
    manifest: Mapping[str, object],
    terminal: Mapping[str, object],
    *,
    store_before: tuple[int, str],
    store_after: tuple[int, str],
) -> None:
    if store_after != store_before:
        _fail("provenance", "campaign_store_mutated", (store_before, store_after))
    _publish(root, manifest, terminal)


def _identity(value: os.stat_result) -> tuple[int, int, int, int, int, int]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
        value.st_nlink,
    )


def _read_output_leaf(directory_fd: int, name: str, before: os.stat_result) -> bytes:
    if (
        stat.S_ISLNK(before.st_mode)
        or not stat.S_ISREG(before.st_mode)
        or before.st_nlink != 1
        or before.st_size < 0
        or before.st_size > _OUTPUT_LEAF_LIMIT
    ):
        _fail("status", "unsafe_output_leaf", name)
    descriptor = -1
    try:
        descriptor = os.open(
            name,
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=directory_fd,
        )
        active = os.fstat(descriptor)
        expected = _identity(before)
        if _identity(active) != expected:
            _fail("status", "output_leaf_substituted", name)
        chunks: list[bytes] = []
        remaining = before.st_size
        while remaining:
            block = os.read(descriptor, min(1 << 20, remaining))
            if not block:
                _fail("status", "output_leaf_short_read", name)
            chunks.append(block)
            remaining -= len(block)
        if os.read(descriptor, 1):
            _fail("status", "output_leaf_grew", name)
        if _identity(os.fstat(descriptor)) != expected:
            _fail("status", "output_leaf_changed", name)
        after = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
        if _identity(after) != expected:
            _fail("status", "output_leaf_substituted", name)
        return b"".join(chunks)
    except OSError as exc:
        raise UR1RunnerError("status", "unsafe_output_leaf", name) from exc
    finally:
        if descriptor != -1:
            os.close(descriptor)


def _snapshot_output(root: Path) -> tuple[dict[str, object], dict[str, str]] | None:
    flags = (
        os.O_RDONLY
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    root_before = root.lstat()
    if stat.S_ISLNK(root_before.st_mode) or not stat.S_ISDIR(root_before.st_mode):
        _fail("status", "unsafe_repository_root", root)
    descriptors: list[int] = []
    identities: list[tuple[int, int, int, int, int, int]] = []
    edge_names: list[str] = []
    leaf_identities: dict[str, tuple[int, int, int, int, int, int]] = {}
    try:
        current = os.open(root, flags)
        descriptors.append(current)
        identities.append(_identity(os.fstat(current)))
        if identities[0] != _identity(root_before):
            _fail("status", "repository_root_substituted", root)
        for part in Path(authority.OUTPUT_NAMESPACE).parts:
            try:
                before = os.stat(part, dir_fd=current, follow_symlinks=False)
            except FileNotFoundError:
                return None
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
                _fail("status", "unsafe_namespace", part)
            child = os.open(part, flags, dir_fd=current)
            if _identity(os.fstat(child)) != _identity(before):
                os.close(child)
                _fail("status", "output_directory_substituted", part)
            descriptors.append(child)
            identities.append(_identity(before))
            edge_names.append(part)
            current = child
        output_fd = descriptors[-1]
        names = tuple(sorted(os.listdir(output_fd)))
        expected_names = ("manifest.json", "terminal.json")
        if names != expected_names:
            _fail("status", "partial_or_foreign_namespace", names)
        values: dict[str, object] = {}
        hashes: dict[str, str] = {}
        for name in names:
            before = os.stat(name, dir_fd=output_fd, follow_symlinks=False)
            leaf_identities[name] = _identity(before)
            raw = _read_output_leaf(output_fd, name, before)
            value = _parse_output_json(raw, name=name)
            values[name] = value
            hashes[name] = sha256(raw).hexdigest()
        if tuple(sorted(os.listdir(output_fd))) != expected_names:
            _fail("status", "output_changed_during_snapshot", names)
        for name, expected in leaf_identities.items():
            after = os.stat(name, dir_fd=output_fd, follow_symlinks=False)
            if _identity(after) != expected:
                _fail("status", "output_leaf_substituted", name)
        for index, (descriptor, expected) in enumerate(
            zip(descriptors, identities, strict=True)
        ):
            if _identity(os.fstat(descriptor)) != expected:
                _fail("status", "output_directory_changed", index)
        for index, name in enumerate(edge_names, start=1):
            after = os.stat(name, dir_fd=descriptors[index - 1], follow_symlinks=False)
            if _identity(after) != identities[index]:
                _fail("status", "output_directory_substituted", name)
        if _identity(root.lstat()) != identities[0]:
            _fail("status", "repository_root_substituted", root)
        return values, hashes
    except OSError as exc:
        raise UR1RunnerError(
            "status", "output_snapshot_unreadable", authority.OUTPUT_NAMESPACE
        ) from exc
    finally:
        for descriptor in reversed(descriptors):
            os.close(descriptor)


def _parse_output_json(raw: bytes, *, name: str) -> dict[str, object]:
    def reject_constant(token: str) -> NoReturn:
        _fail("status", "nonfinite_json_constant", (name, token))

    try:
        value = json.loads(
            raw,
            object_pairs_hook=_unique,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise UR1RunnerError("status", "noncanonical_output", name) from exc
    if not isinstance(value, dict) or raw != _pretty(value):
        _fail("status", "noncanonical_output", name)
    return value


def _require_keys(
    value: object, expected: set[str], *, label: str
) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != expected:
        _fail("status", "terminal_schema", (label, sorted(expected)))
    return value


def expected_replay_receipt() -> dict[str, object]:
    return authority.expected_raw_replay_receipt()


def _validate_replay_receipt(value: object) -> dict[str, object]:
    expected = expected_replay_receipt()
    if not isinstance(value, Mapping) or set(value) != set(expected):
        _fail("status", "replay_receipt_schema", (sorted(expected), value))
    receipt = {key: value[key] for key in expected}
    if receipt != expected:
        _fail("status", "replay_receipt_identity", receipt)
    return receipt


def _terminal_base(
    authority_commit: str,
    classification: str,
    before: tuple[int, str],
    after: tuple[int, str],
) -> dict[str, object]:
    if classification not in TERMINAL_CLASSES:
        _fail("reduction", "terminal_class", classification)
    return {
        "schema": RAW_SCHEMA,
        "artifact_id": authority.ARTIFACT_ID,
        "runner_id": RUNNER_ID,
        "classification": classification,
        "authority_commit": authority_commit,
        "experiment_label": authority.SEMANTICS["experiment_label"],
        "tableau_selector": "SSPRK3",
        "tableau_runtime_selector": COMPARATOR_METHOD,
        "actual_spatial_operator": "inherited_RK4_2049_SBP4",
        "retry": authority.RETRY,
        "attempted_width_hex": authority.WIDTH_HEX,
        "published_channel": authority.PUBLISHED_CHANNEL,
        "production_SSPRK3_comparator": False,
        "independent_method_agreement": False,
        "production_method_earned": False,
        "store_snapshot_before": {"leaf_count": before[0], "sha256": before[1]},
        "store_snapshot_after": {"leaf_count": after[0], "sha256": after[1]},
        "store_unchanged": before == after,
        "PDE_state_committed": False,
        "temporal_retry_admission_called": False,
        "fine_path_committed": False,
        "successor_remedy_selected": False,
        "common_event_completed": False,
        "GR0_calibration_completed": False,
        "candidate_branch_opened": False,
        "mechanism_result_earned": False,
        "physical_result_earned": False,
        "completed_result_licenses_only_later_prospectively_frozen_envelope_or_admission_design": (
            classification in OWNER_COMPLETED_TERMINAL_CLASSES
        ),
        "related_TI2_complete_C_class_is_not_replacement_production_admission": True,
    }


def _manifest(authority_commit: str) -> dict[str, object]:
    return {
        "schema": RAW_SCHEMA,
        "artifact_id": authority.ARTIFACT_ID,
        "runner_id": RUNNER_ID,
        "authority_commit": authority_commit,
        "AC1_PREF1_result_sha256": authority.AC1_PREF1_RESULT_SHA256,
        "experiment_label": authority.SEMANTICS["experiment_label"],
        "tableau_selector": "SSPRK3",
        "actual_spatial_operator": "inherited_RK4_2049_SBP4",
        "retry": authority.RETRY,
        "published_channel": authority.PUBLISHED_CHANNEL,
        "production_SSPRK3_comparator": False,
        "output_leaves": ["manifest.json", "terminal.json"],
    }


def _validate_terminal(
    value: object, *, authority_commit: str
) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        _fail("status", "terminal_schema", "terminal is not a mapping")
    classification = value.get("classification")
    if classification not in PUBLISHABLE_TERMINAL_CLASSES:
        _fail("status", "terminal_class_not_publishable", classification)
    sealed = (
        authority.SEALED_STORE_LEAF_COUNT,
        authority.SEALED_STORE_SNAPSHOT_SHA256,
    )
    expected_base = _terminal_base(
        authority_commit, str(classification), sealed, sealed
    )
    if any(value.get(name) != expected for name, expected in expected_base.items()):
        _fail("status", "terminal_identity", classification)
    base_keys = set(expected_base)
    if classification == "shadow_proposal_premise_stop":
        _require_keys(
            value,
            base_keys | {"typed_stop", "replay_receipt"},
            label="premise stop terminal",
        )
        stop = _require_keys(
            value["typed_stop"], {"type", "detail"}, label="premise stop"
        )
        if (
            stop["type"] != "TDG6RefinementPathStop"
            or not isinstance(stop["detail"], str)
            or not stop["detail"]
            or len(stop["detail"]) > 640
        ):
            _fail("status", "premise_stop_terminal", stop)
        _validate_replay_receipt(value["replay_receipt"])
        return value
    _require_keys(
        value,
        base_keys
        | {
            "replay_receipt",
            "shadow_path_count",
            "shadow_proposal_count",
            "SSPRK3_stage_and_endpoint_record_count",
            "row_stream",
            "u_R",
        },
        label="completed terminal",
    )
    _validate_replay_receipt(value["replay_receipt"])
    row_stream = validate_row_stream(value.get("row_stream"))
    u_r = validate_u_R(value.get("u_R"))
    expected = reduce_ur1_terminal(
        rows_matched=bool(row_stream["matched"]),
        intervals_match_ac1=bool(u_r["matches_sealed_AC1_production_intervals"]),
        owners_equal=bool(u_r["D01_D12_lower_owners_equal"]),
    )
    if (
        classification != expected
        or value.get("shadow_path_count") != 7
        or value.get("shadow_proposal_count") != 7
        or value.get("SSPRK3_stage_and_endpoint_record_count") != 28
        or value.get("production_method_earned") is not False
        or value.get("published_channel") != authority.PUBLISHED_CHANNEL
    ):
        _fail("status", "completed_terminal_counts", classification)
    return value


def _inspect_output(root: Path, *, authority_commit: str) -> dict[str, object]:
    snapshot = _snapshot_output(root)
    if snapshot is None:
        return {"state": "absent"}
    values, hashes = snapshot
    manifest = values["manifest.json"]
    if manifest != _manifest(authority_commit):
        _fail("status", "manifest_identity", manifest)
    terminal = _validate_terminal(
        values["terminal.json"], authority_commit=authority_commit
    )
    return {
        "state": "terminal",
        "hashes": hashes,
        "classification": terminal["classification"],
    }


def status(root: Path, *, authority_commit: str) -> dict[str, object]:
    _execution_authority(root.resolve(), authority_commit)
    before = _snapshot_store(root.resolve())
    if before != (
        authority.SEALED_STORE_LEAF_COUNT,
        authority.SEALED_STORE_SNAPSHOT_SHA256,
    ):
        _fail("status", "sealed_store_identity", before)
    output = _inspect_output(root.resolve(), authority_commit=authority_commit)
    return {
        "artifact_id": authority.ARTIFACT_ID,
        "authority_commit": authority_commit,
        "output": output,
        "store": {"leaf_count": before[0], "sha256": before[1]},
        "safe_to_run": output["state"] == "absent",
        "state_advance_authorized": False,
    }


def _restore_replay(root: Path) -> ti2_runner.RestoredReplay:
    store = HLT16CampaignStore(root / ar1.PREF2_STORE_PATH)
    shells = build_static_gr0_shells(root)
    try:
        return ti2_runner._restore_replay(root, store, shells, authority.REPLAY)
    except ti2_runner.TI1RunnerError as exc:
        raise UR1RunnerError(exc.owner, exc.code, exc.detail) from exc


def _prepare_shadow(
    restored: ti2_runner.RestoredReplay,
) -> tdg6.TDG6PreparedGR0Compositor:
    try:
        return ti2_runner._prepare_shadow(restored)
    except tdg6.TDG6RefinementPathStop:
        raise
    except ti2_runner.TI1RunnerError as exc:
        raise UR1RunnerError(exc.owner, exc.code, exc.detail) from exc


def _replay_receipt(restored: ti2_runner.RestoredReplay) -> dict[str, object]:
    fingerprint = restored.fingerprint
    if not isinstance(fingerprint, Mapping):
        _fail("replay", "fingerprint_type", type(fingerprint).__name__)
    return _validate_replay_receipt(
        {
            "member_key": fingerprint.get("member_key"),
            "retry": authority.RETRY,
            "attempted_width_hex": authority.WIDTH_HEX,
            "accepted_time_hex": fingerprint.get("accepted_time_hex"),
            "state_sha256": fingerprint.get("state_sha256"),
            "descriptor_sha256": fingerprint.get("descriptor_sha256"),
            "transaction_sha256": fingerprint.get("transaction_sha256"),
            "historical_journal_sha256": restored.historical_journal_sha256,
        }
    )


def run(root: Path, *, authority_commit: str) -> dict[str, object]:
    repository = root.resolve()
    receipt = _execution_authority(repository, authority_commit)
    if (
        receipt.tableau_selector != "SSPRK3"
        or receipt.retry != authority.RETRY
        or receipt.published_channel != authority.PUBLISHED_CHANNEL
        or receipt.state_advance_authorized
    ):
        _fail("authority", "receipt_semantics", receipt)
    authority.require_output_absent(repository)
    before = _snapshot_store(repository)
    sealed = (
        authority.SEALED_STORE_LEAF_COUNT,
        authority.SEALED_STORE_SNAPSHOT_SHA256,
    )
    if before != sealed:
        _fail("provenance", "sealed_store_identity", (before, sealed))
    replay_receipt: dict[str, object] = {}
    try:
        restored = _restore_replay(repository)
        replay_receipt = _replay_receipt(restored)
        prepared = _prepare_shadow(restored)
    except tdg6.TDG6RefinementPathStop as exc:
        after = _snapshot_store(repository)
        if after != before:
            _fail("provenance", "campaign_store_mutated", (before, after))
        terminal = {
            **_terminal_base(
                authority_commit, "shadow_proposal_premise_stop", before, after
            ),
            "typed_stop": {
                "type": type(exc).__name__,
                "detail": " ".join(str(exc).split())[:640],
            },
            "replay_receipt": replay_receipt,
        }
        _publish_terminal(
            repository,
            _manifest(authority_commit),
            terminal,
            store_before=before,
            store_after=after,
        )
        return terminal

    row_stream = hash_u_R_rows(prepared)
    serialized = serialize_u_R(prepared.continuous_admission)
    classification = reduce_ur1_terminal(
        rows_matched=bool(row_stream["matched"]),
        intervals_match_ac1=bool(serialized["matches_sealed_AC1_production_intervals"]),
        owners_equal=bool(serialized["D01_D12_lower_owners_equal"]),
    )
    after = _snapshot_store(repository)
    if after != before:
        _fail("provenance", "campaign_store_mutated", (before, after))
    terminal = {
        **_terminal_base(authority_commit, classification, before, after),
        "replay_receipt": replay_receipt,
        "shadow_path_count": 7,
        "shadow_proposal_count": 7,
        "SSPRK3_stage_and_endpoint_record_count": 28,
        "row_stream": row_stream,
        "u_R": serialized,
    }
    _publish_terminal(
        repository,
        _manifest(authority_commit),
        terminal,
        store_before=before,
        store_after=after,
    )
    return terminal


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authority-commit", required=True)
    parser.add_argument("--status", action="store_true")
    arguments = parser.parse_args(argv)
    try:
        result = (
            status(ROOT, authority_commit=arguments.authority_commit)
            if arguments.status
            else run(ROOT, authority_commit=arguments.authority_commit)
        )
    except UR1RunnerError as exc:
        print(
            json.dumps(
                {
                    "artifact_id": authority.ARTIFACT_ID,
                    "classification": "invalid_provenance_or_implementation",
                    "owner": exc.owner,
                    "code": exc.code,
                    "detail": exc.detail,
                },
                sort_keys=True,
                separators=(",", ":"),
            )
        )
        return 2
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
