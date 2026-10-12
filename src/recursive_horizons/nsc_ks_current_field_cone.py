"""Whole-cone field accumulator over streamed difference-reconstruction cells.

This owner sums directed residual enclosures of history-minus-reference
cells. It restores an existing v4 local-Fourier payload by hash. It does
not run a family campaign, fill the physical rho=1 source error, or issue
EXISTENCE or NON_EXISTENCE.
"""
from copy import deepcopy
from collections import deque
from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import json
from types import MappingProxyType

from flint import arb, ctx
import numpy as np

from .nsc_ks_ball_operator import AnalyticRadiusFamily, time_operator_enclosure
from .nsc_ks_ball_trajectory import exact_upper, restored_upper
from .nsc_ks_current_history_bounds import radius_bounds, rational_record, value_integral_bounds
from .nsc_ks_difference_error import difference_matter_error, propagate_difference_error
from .nsc_ks_finite_matter_error import source_norm_upper
from .nsc_ks_reference_error import homogeneous_reference_norms
from .nsc_ks_residual_error import pure_radius_commutator_integrals
from .nsc_ks_difference_residual_polynomial import (
    DifferenceResidualPolynomial, ball_split, difference_operator_remainder_bounds,
    difference_polynomial_bounds, monomial_keys,
)
from .nsc_ks_fourier_residual_bound import profile_sup_bounds
from .nsc_ks_local_profile_fourier_bound import (
    LocalProfileFourierBound, profile_payload_digest, restore_local_profile_fourier,
)
from .nsc_ks_profile_identity import profile_identity
from .nsc_ks_radius_enclosure import reciprocal_radius_tail
from .nsc_ks_trajectory import TrajectorySegment
from .nsc_local_incoming_family import LocalIncomingFamily


SCHEMA = "NSC-KS-CURRENT-FIELD-CONE-v1"
CHECKPOINT_SCHEMA = "NSC-KS-CURRENT-FIELD-CONE-CHECKPOINT-v1"
V4_SCHEMA = "NSC-KS-CURRENT-FIELD-PILOT-v4"
V4_PROFILE_IDENTITY = (
    "0b0e4cedfdb695d342c7f7347740631e916bdaee3a977df56808c725b25dacb0")
RUNTIME_RECORD_KEYS = ("runtime", "memory")
NONSCIENTIFIC_RECORD_KEYS = RUNTIME_RECORD_KEYS + ("forecast", "scientific_digest")
CONTINUOUS_INPUTS_SCHEMA = "NSC-KS-WHOLE-CONE-CONTINUOUS-INPUTS-v1"
PROPAGATION_INPUT_NAMES = (
    "reference_initial",
    "reference_residual",
    "reference_offdiagonal_integral",
    "difference_initial",
    "difference_residual_integrals",
    "offdiagonal_integral",
    "Bz_integral",
    "M_integral",
    "Mz_integral",
    "maximum_absolute_energy",
)
MATTER_INPUT_NAMES = (
    "reference_norm",
    "difference_norm",
    "reference_axial_norm",
    "difference_axial_norm",
    "reference_error",
    "difference_error",
    "reference_axial_error",
    "difference_axial_error",
    "source_norm",
    "mass",
    "absolute_angular",
    "axial_lower",
    "radius_lower",
    "multiplicity",
)


def canonical_dumps(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      allow_nan=False)


def canonical_digest(value):
    return sha256(canonical_dumps(value).encode("utf-8")).hexdigest()


def scientific_record(value):
    """Drop runtime, memory, forecast, and the digest field from a scientific digest."""
    if not isinstance(value, dict):
        raise TypeError("mapping required")
    return {key: item for key, item in value.items()
            if key not in NONSCIENTIFIC_RECORD_KEYS}


def scientific_digest(value):
    return canonical_digest(scientific_record(value))


def _positive_int(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError("positive integer " + name + " required")
    return value


def _finite_float(value, name):
    value = float(value)
    if not np.isfinite(value):
        raise ValueError("finite " + name + " required")
    return value


def _sha256_text(value, name):
    if (not isinstance(value, str) or len(value) != 64
            or any(char not in "0123456789abcdef" for char in value)):
        raise ValueError("lowercase SHA-256 " + name + " required")
    return value


def _hex_floats(values, name):
    return [_finite_float(value, name).hex() for value in values]


def _arb_q(value):
    if isinstance(value, arb):
        return value
    if isinstance(value, Fraction):
        return arb(value.numerator) / arb(value.denominator)
    if hasattr(value, "numerator") and hasattr(value, "denominator"):
        return arb(int(value.numerator)) / arb(int(value.denominator))
    return arb(value)


def _upper_pair(values, name):
    if not isinstance(values, (tuple, list)) or len(values) != 2:
        raise ValueError(name + " pair required")
    return tuple(value.upper() if isinstance(value, arb) else restored_upper(value)
                 for value in values)


def _pack_pair(values):
    return [exact_upper(value) for value in values]


def _zero_pair(bits):
    with ctx.workprec(bits):
        return (arb(0), arb(0))


def packed_radius_tail(family, rho_up, angular, order, *, bits):
    """Exact-upper potential remainder of the owned reciprocal-radius tail."""
    if not isinstance(family, LocalIncomingFamily):
        raise TypeError("LocalIncomingFamily required")
    bounds = radius_bounds(family, rho_up)
    tail = reciprocal_radius_tail(
        bounds["w"]["profile_bounds"], bounds["U"]["profile_bounds"],
        normal_support=bounds["normal_support"],
        reference_radius_lower=bounds["reference_radius_lower"],
        axial_lower=bounds["axial_lower"],
        absolute_angular=abs(angular), order=order)
    with ctx.workprec(bits):
        return tuple(exact_upper(_arb_q(value).upper())
                     for value in tail["potential_derivative_tail_bounds"][:2])


def profile_rows_from_bounds(rows):
    """Convert local Fourier balls to residual-bound coefficient tables."""
    if not isinstance(rows, dict) or not rows:
        raise ValueError("local Fourier profile rows required")
    result = {}
    for key, row in rows.items():
        if not isinstance(row, LocalProfileFourierBound):
            raise TypeError("LocalProfileFourierBound required")
        result[tuple(key)] = {
            "coefficients": list(row.coefficients),
            "tail": list(row.uniform_tail_bounds),
        }
    return result


def restore_v4_profile_payload(record, *, expected_payload_sha256=None,
                               expected_file_sha256=None, file_bytes=None):
    """Hash-validate and restore the immutable v4 coefficient payload."""
    if expected_file_sha256 is not None:
        if not isinstance(file_bytes, (bytes, bytearray)):
            raise ValueError("mismatched profile payload")
        if sha256(file_bytes).hexdigest() != expected_file_sha256:
            raise ValueError("mismatched profile payload")
    if not isinstance(record, dict) or record.get("schema") != V4_SCHEMA:
        raise ValueError("unexpected field-v4 schema")
    payload = record.get("profile_coefficient_payload")
    digest = profile_payload_digest(payload)
    if expected_payload_sha256 is not None and digest != expected_payload_sha256:
        raise ValueError("mismatched profile payload")
    return restore_local_profile_fourier(payload, expected_sha256=digest), digest


def enclose_difference_segment(config, family, profile_rows, radius_tail,
                               segment):
    """One-cell directed residual: history-minus-reference before norms."""
    if not isinstance(config, WholeConeFieldConfig):
        raise TypeError("WholeConeFieldConfig required")
    if not isinstance(family, LocalIncomingFamily):
        raise TypeError("LocalIncomingFamily required")
    if not isinstance(segment, TrajectorySegment):
        raise TypeError("saved anchored TrajectorySegment required")
    keys = monomial_keys(config.reciprocal_order)
    if set(profile_rows) != set(keys):
        raise ValueError("profile keys differ from the reciprocal keys")
    with ctx.workprec(config.bits):
        length = (arb(config.period_length_numerator)
                  / arb(config.period_length_denominator))
        tail = tuple(restored_upper(value) for value in config.radius_tail)
        if radius_tail is not None:
            tail = tuple(value.upper() if isinstance(value, arb)
                         else restored_upper(value) for value in radius_tail)
        field_x, field_a, field_d = ball_split(
            segment, np.asarray(config.source_weights, float),
            config.spatial_count, length, bits=config.bits)
        polynomials, time_errors = time_operator_enclosure(
            AnalyticRadiusFamily(family), field_x.rho_start, field_x.rho_end,
            config.angular, degree=config.time_degree,
            reciprocal_order=config.reciprocal_order, bits=config.bits)
        potential = {key: value for key, value in polynomials.items()
                     if isinstance(key, tuple)}
        if set(potential) != set(keys):
            raise ValueError("time-series monomials differ from the reciprocal keys")
        residual = DifferenceResidualPolynomial(
            field_x, field_a, polynomials["inv_a2"], polynomials["inv_a"],
            polynomials["inv_ar"], potential, config.mass, config.angular,
            np.asarray(config.source_energies, float))
        polynomial = difference_polynomial_bounds(residual, profile_rows)
        remainder = difference_operator_remainder_bounds(
            field_d, field_a, time_errors,
            profile_sup_bounds(profile_rows, field_x.length),
            tail, config.mass, config.angular,
            np.asarray(config.source_energies, float))
        total = tuple((left + right).upper()
                      for left, right in zip(polynomial, remainder))
        width = (arb(segment.rho_start) - arb(segment.rho_end)).upper()
        # DifferenceResidualPolynomial bounds S = X_xi - h L X.  Since
        # d rho = h d xi, integral |R| d rho = integral |S| d xi.  The
        # time/radius remainder owner includes the same h already.  A second
        # multiplication by the physical cell width would understate the
        # Gronwall input by one step-size factor.
        integral = tuple(bound.upper() for bound in total)
        return {
            "polynomial": polynomial,
            "remainder": remainder,
            "total": total,
            "width": width,
            "integral": integral,
        }


def whole_cone_propagate_difference_error(
        residual_integrals, *, reference_initial, reference_residual,
        reference_offdiagonal_integral, difference_initial,
        offdiagonal_integral, Bz_integral, M_integral, Mz_integral,
        maximum_absolute_energy, bits=160):
    """Adapter: Gronwall majorant with physical rho=1 source error left open."""
    bound = dict(propagate_difference_error(
        reference_initial, reference_residual, reference_offdiagonal_integral,
        difference_initial, residual_integrals, offdiagonal_integral,
        Bz_integral, M_integral, Mz_integral, maximum_absolute_energy,
        bits=bits))
    bound["source_accuracy_included"] = False
    bound["tangent_accuracy_included"] = False
    bound["physical_rho1_source_error"] = None
    bound["physical_EXISTENCE_certificate"] = False
    bound["physical_NONEXISTENCE_certificate"] = False
    bound["physical_local_gate"] = (
        "OPEN: physical rho=1 source error unresolved")
    return bound


def whole_cone_difference_matter_error(
        reference_norm, difference_norm, reference_axial_norm,
        difference_axial_norm, reference_error, difference_error,
        reference_axial_error, difference_axial_error, source_norm, *,
        mass, absolute_angular, axial_lower, radius_lower, multiplicity,
        bits=160):
    """Adapter: difference contraction with source accuracy left open."""
    bound = dict(difference_matter_error(
        reference_norm, difference_norm, reference_axial_norm,
        difference_axial_norm, reference_error, difference_error,
        reference_axial_error, difference_axial_error, source_norm,
        mass=mass, absolute_angular=absolute_angular, axial_lower=axial_lower,
        radius_lower=radius_lower, multiplicity=multiplicity, bits=bits))
    bound["source_accuracy_included"] = False
    bound["physical_rho1_source_error"] = None
    bound["physical_EXISTENCE_certificate"] = False
    bound["physical_NONEXISTENCE_certificate"] = False
    bound["physical_local_gate"] = (
        "OPEN: physical rho=1 source error unresolved")
    return bound


def _input_slot(value, *, owner=None, reason=None):
    if value is None:
        if not isinstance(reason, str) or not reason:
            raise ValueError("named null requires an explicit reason")
        return {
            "value": None,
            "owner": None,
            "reason": reason,
        }
    if not isinstance(owner, str) or not owner:
        raise ValueError("owned input requires an analytic owner")
    return {
        "value": value,
        "owner": owner,
        "reason": None,
    }


def _pack_scalar(value, bits):
    with ctx.workprec(bits):
        return exact_upper(_arb_q(value).upper())


def _slot_complete(slot):
    return isinstance(slot, dict) and slot.get("value") is not None


def restore_continuous_slot(slot):
    """Restore a packed slot to an adapter argument, preserving named nulls."""
    if not isinstance(slot, dict) or "value" not in slot:
        raise TypeError("continuous-input slot required")
    value = slot["value"]
    if value is None:
        return None
    if isinstance(value, dict) and "mantissa" in value and "exponent" in value:
        return restored_upper(value)
    if isinstance(value, (list, tuple)):
        restored = []
        for item in value:
            if isinstance(item, dict) and "mantissa" in item and "exponent" in item:
                restored.append(restored_upper(item))
            else:
                restored.append(item)
        return tuple(restored)
    return value


def _channel_multiplicity(channel):
    if not isinstance(channel, Mapping):
        return None
    copies, degeneracy = channel.get("copy_count"), channel.get("degeneracy")
    angular = channel.get("angular_eigenvalue")
    if copies is None or degeneracy is None or angular is None:
        return None
    if any(isinstance(value, bool) or int(value) != value or value <= 0
           for value in (copies, degeneracy)):
        raise ValueError("positive integer copies/degeneracy required")
    signs = 1 if float(angular) == 0 else 2
    return int(copies) * int(degeneracy) / signs


def whole_cone_continuous_inputs(
        family, *, mass, angular, rho_up, source_energies,
        source_weights=None, source_covariance=None, channel=None,
        residual_integrals=None, reference_amplitudes=None, bits=160):
    """Map existing analytic owners onto Gronwall/contraction slots.

    Unowned quantities stay explicit named nulls. Numerical zero is not
    substituted for missing source accuracy, sampled defects, or an unowned
    M / M_z integral.
    """
    if not isinstance(family, LocalIncomingFamily):
        raise TypeError("LocalIncomingFamily required")
    bits = _positive_int(bits, "bits")
    energies = tuple(_finite_float(value, "source energy")
                     for value in source_energies)
    if not energies:
        raise ValueError("source energies required")
    bounds = radius_bounds(family, rho_up)
    integrals = value_integral_bounds(bounds, mass, angular)
    commutators = pure_radius_commutator_integrals(
        bounds["w"]["profile_bounds"], bounds["U"]["profile_bounds"],
        bounds["normal_support"], bounds["radius_lower"], abs(angular))
    k0 = _pack_scalar(integrals["K0"], bits)
    bz = _pack_scalar(commutators[0], bits)
    with ctx.workprec(bits):
        difference_initial = _pack_pair(_zero_pair(bits))
        energy = _pack_scalar(max(abs(value) for value in energies), bits)
        mass_packed = _pack_scalar(abs(mass), bits)
        angular_packed = _pack_scalar(abs(angular), bits)
        axial = _pack_scalar(bounds["axial_lower"], bits)
        radius = _pack_scalar(bounds["radius_lower"], bits)
    residual_slot = _input_slot(
        None if residual_integrals is None else list(residual_integrals),
        owner=("WholeConeAccumulator.residual_integrals"
               if residual_integrals is not None else None),
        reason=(None if residual_integrals is not None else
                "difference residual integrals exist only after enclosed cells"))
    source_slot = _input_slot(
        None if source_covariance is None else _pack_scalar(
            source_norm_upper(source_covariance, bits=bits), bits),
        owner=("nsc_ks_finite_matter_error.source_norm_upper"
               if source_covariance is not None else None),
        reason=(None if source_covariance is not None else
                "coherent source matrix was not supplied"))
    multiplicity = _channel_multiplicity(channel)
    multiplicity_slot = _input_slot(
        None if multiplicity is None else _pack_scalar(multiplicity, bits),
        owner=("channel copy_count*degeneracy/angular_signs"
               if multiplicity is not None else None),
        reason=(None if multiplicity is not None else
                "retained-ledger copies and degeneracy were not supplied"))
    if reference_amplitudes is None or source_weights is None:
        reference_norm = reference_axial = None
        reference_owner = None
        reference_reason = (
            "computed homogeneous A norms require the evolved reference "
            "amplitudes; grid samples of D are not a continuous owner")
    else:
        ref_norm, ref_axial = homogeneous_reference_norms(
            reference_amplitudes, source_weights, energies, bits=bits)
        reference_norm = _pack_scalar(ref_norm, bits)
        reference_axial = _pack_scalar(ref_axial, bits)
        reference_owner = "nsc_ks_reference_error.homogeneous_reference_norms"
        reference_reason = None
    missing_field = (
        "no owned continuous sup-norm of computed delta F; nodal samples "
        "are not promoted")
    missing_error = (
        "propagate_difference_error outputs require every Gronwall input "
        "to be owned")
    propagation = {
        "reference_initial": _input_slot(
            None, reason=(
                "physical rho=1 source column error unresolved; not replaced "
                "by numerical zero")),
        "reference_residual": _input_slot(
            None, reason=(
                "no owned continuous enclosure of A'-G_ref A; DOP853 "
                "tolerances are not a residual certificate")),
        "reference_offdiagonal_integral": _input_slot(
            k0, owner="nsc_ks_current_history_bounds.value_integral_bounds.K0"),
        "difference_initial": _input_slot(
            list(difference_initial),
            owner="evolve_ks_difference_envelope exact zero D initialization"),
        "difference_residual_integrals": residual_slot,
        "offdiagonal_integral": _input_slot(
            k0, owner="nsc_ks_current_history_bounds.value_integral_bounds.K0"),
        "Bz_integral": _input_slot(
            bz, owner="nsc_ks_residual_error.pure_radius_commutator_integrals"),
        "M_integral": _input_slot(
            None, reason=(
                "no owned integral of the A-to-D radius coupling M")),
        "Mz_integral": _input_slot(
            None, reason="no owned integral of M_z"),
        "maximum_absolute_energy": _input_slot(
            energy, owner="selected source energy labels"),
    }
    matter = {
        "reference_norm": _input_slot(
            reference_norm, owner=reference_owner, reason=reference_reason),
        "difference_norm": _input_slot(None, reason=missing_field),
        "reference_axial_norm": _input_slot(
            reference_axial, owner=reference_owner, reason=reference_reason),
        "difference_axial_norm": _input_slot(None, reason=missing_field),
        "reference_error": _input_slot(None, reason=missing_error),
        "difference_error": _input_slot(None, reason=missing_error),
        "reference_axial_error": _input_slot(None, reason=missing_error),
        "difference_axial_error": _input_slot(None, reason=missing_error),
        "source_norm": source_slot,
        "mass": _input_slot(mass_packed, owner="selected source channel mass"),
        "absolute_angular": _input_slot(
            angular_packed, owner="selected source channel |angular|"),
        "axial_lower": _input_slot(
            axial, owner="nsc_ks_current_history_bounds.radius_bounds"),
        "radius_lower": _input_slot(
            radius, owner="nsc_ks_current_history_bounds.radius_bounds"),
        "multiplicity": multiplicity_slot,
    }
    if set(propagation) != set(PROPAGATION_INPUT_NAMES):
        raise ValueError("propagation input inventory changed")
    if set(matter) != set(MATTER_INPUT_NAMES):
        raise ValueError("matter input inventory changed")
    return {
        "schema": CONTINUOUS_INPUTS_SCHEMA,
        "propagation": propagation,
        "matter": matter,
        "geometry": {
            "rho_up": rational_record(bounds["rho_up"]),
            "axial_lower": rational_record(bounds["axial_lower"]),
            "radius_lower": rational_record(bounds["radius_lower"]),
            "normal_support": rational_record(bounds["normal_support"]),
            "w_profile_bounds": rational_record(bounds["w"]["profile_bounds"]),
            "U_profile_bounds": rational_record(bounds["U"]["profile_bounds"]),
            "K0": rational_record(integrals["K0"]),
            "K1": rational_record(integrals["K1"]),
            "characteristic_length": rational_record(integrals["D"]),
            "history_tangent_bound": integrals["history_tangent_bound"],
            "Bz_integral": rational_record(commutators[0]),
            "Bzz_integral": rational_record(commutators[1]),
        },
        "physical_rho1_source_error": None,
        "source_accuracy_included": False,
        "complete_for_propagate_difference_error": all(
            _slot_complete(slot) for slot in propagation.values()),
        "complete_for_difference_matter_error": all(
            _slot_complete(slot) for slot in matter.values()),
    }


def with_residual_integrals(inputs, residual_integrals):
    """Attach enclosed residual integrals without filling any other null."""
    if not isinstance(inputs, Mapping) or inputs.get("schema") != CONTINUOUS_INPUTS_SCHEMA:
        raise ValueError("whole-cone continuous-input inventory required")
    if residual_integrals is None:
        raise ValueError("residual integrals required")
    updated = deepcopy(dict(inputs))
    updated["propagation"] = dict(updated["propagation"])
    updated["propagation"]["difference_residual_integrals"] = _input_slot(
        list(residual_integrals),
        owner="WholeConeAccumulator.residual_integrals")
    updated["complete_for_propagate_difference_error"] = all(
        _slot_complete(slot) for slot in updated["propagation"].values())
    updated["complete_for_difference_matter_error"] = all(
        _slot_complete(slot) for slot in updated["matter"].values())
    return updated


def propagation_adapter_kwargs(inputs):
    """Return adapter kwargs only when every Gronwall slot is owned."""
    if not isinstance(inputs, Mapping) or inputs.get("schema") != CONTINUOUS_INPUTS_SCHEMA:
        raise ValueError("whole-cone continuous-input inventory required")
    if not inputs.get("complete_for_propagate_difference_error"):
        return None
    values = {name: restore_continuous_slot(inputs["propagation"][name])
              for name in PROPAGATION_INPUT_NAMES}
    values["residual_integrals"] = values.pop("difference_residual_integrals")
    return values


def matter_adapter_kwargs(inputs):
    """Return contraction kwargs only when every matter slot is owned."""
    if not isinstance(inputs, Mapping) or inputs.get("schema") != CONTINUOUS_INPUTS_SCHEMA:
        raise ValueError("whole-cone continuous-input inventory required")
    if not inputs.get("complete_for_difference_matter_error"):
        return None
    return {name: restore_continuous_slot(inputs["matter"][name])
            for name in MATTER_INPUT_NAMES}


@dataclass(frozen=True)
class WholeConeFieldConfig:
    family: tuple
    profile_identity: str
    bits: int
    time_degree: int
    reciprocal_order: int
    spatial_count: int
    period_origin: float
    period_length_numerator: int
    period_length_denominator: int
    mass: float
    angular: float
    source_weights: tuple
    source_energies: tuple
    radius_tail: tuple
    profile_payload_sha256: str
    settings: object
    source_binding_sha256: str = ""
    settings_digest: str = ""

    def __post_init__(self):
        if (not isinstance(self.family, (tuple, list)) or len(self.family) != 2
                or any(isinstance(value, bool) or not isinstance(value, int)
                       for value in self.family)):
            raise ValueError("family pair required")
        family = tuple(self.family)
        if family[0] < 0 or family[1] not in (-1, 1):
            raise ValueError("family pair required")
        object.__setattr__(self, "family", family)
        identity = _sha256_text(self.profile_identity, "profile identity")
        object.__setattr__(self, "profile_identity", identity)
        object.__setattr__(self, "bits", _positive_int(self.bits, "bits"))
        if self.bits < 64:
            raise ValueError("at least 64 Arb bits required")
        object.__setattr__(self, "time_degree",
                           _positive_int(self.time_degree, "time degree"))
        object.__setattr__(self, "reciprocal_order", _positive_int(
            self.reciprocal_order, "reciprocal order"))
        object.__setattr__(self, "spatial_count", _positive_int(
            self.spatial_count, "spatial count"))
        if self.spatial_count < 8:
            raise ValueError("owned spatial mode count required")
        object.__setattr__(self, "period_origin",
                           _finite_float(self.period_origin, "period origin"))
        num = self.period_length_numerator
        den = self.period_length_denominator
        if (isinstance(num, bool) or isinstance(den, bool)
                or not isinstance(num, int) or not isinstance(den, int)
                or den < 1 or num < 1):
            raise ValueError("positive rational period required")
        object.__setattr__(self, "mass", _finite_float(self.mass, "mass"))
        object.__setattr__(self, "angular",
                           _finite_float(self.angular, "angular"))
        weights = tuple(_finite_float(value, "source weight")
                        for value in self.source_weights)
        energies = tuple(_finite_float(value, "source energy")
                         for value in self.source_energies)
        if not weights or len(weights) != len(energies) or any(w <= 0 for w in weights):
            raise ValueError("original positive source column weights required")
        object.__setattr__(self, "source_weights", weights)
        object.__setattr__(self, "source_energies", energies)
        tail = tuple(self.radius_tail)
        if len(tail) != 2:
            raise ValueError("radius tail pair required")
        if any(not restored_upper(value) >= 0 for value in tail):
            raise ValueError("nonnegative radius tail pair required")
        object.__setattr__(self, "radius_tail", tail)
        payload = _sha256_text(
            self.profile_payload_sha256, "profile payload digest")
        object.__setattr__(self, "profile_payload_sha256", payload)
        if not isinstance(self.settings, Mapping):
            raise TypeError("settings mapping required")
        settings = dict(self.settings)
        canonical_dumps(settings)
        object.__setattr__(self, "settings", MappingProxyType(settings))
        source_binding = canonical_digest({
            "family": list(family),
            "profile_identity": identity,
            "spatial_count": self.spatial_count,
            "period_origin": self.period_origin.hex(),
            "period_length": [num, den],
            "mass": self.mass.hex(),
            "angular": self.angular.hex(),
            "source_weights": _hex_floats(weights, "source weight"),
            "source_energies": _hex_floats(energies, "source energy"),
        })
        settings_digest = canonical_digest({
            "bits": self.bits,
            "time_degree": self.time_degree,
            "reciprocal_order": self.reciprocal_order,
            "settings": settings,
            "history_minus_reference_before_norms": True,
        })
        if self.source_binding_sha256 and self.source_binding_sha256 != source_binding:
            raise ValueError("mismatched source")
        if self.settings_digest and self.settings_digest != settings_digest:
            raise ValueError("mismatched settings")
        object.__setattr__(self, "source_binding_sha256", source_binding)
        object.__setattr__(self, "settings_digest", settings_digest)

    def identity_payload(self):
        return {
            "family": list(self.family),
            "profile_identity": self.profile_identity,
            "bits": self.bits,
            "time_degree": self.time_degree,
            "reciprocal_order": self.reciprocal_order,
            "spatial_count": self.spatial_count,
            "period_origin": self.period_origin.hex(),
            "period_length": [self.period_length_numerator,
                              self.period_length_denominator],
            "mass": self.mass.hex(),
            "angular": self.angular.hex(),
            "source_weights": _hex_floats(self.source_weights, "source weight"),
            "source_energies": _hex_floats(self.source_energies, "source energy"),
            "radius_tail": list(self.radius_tail),
            "profile_payload_sha256": self.profile_payload_sha256,
            "settings": dict(self.settings),
            "source_binding_sha256": self.source_binding_sha256,
            "settings_digest": self.settings_digest,
            "history_minus_reference_before_norms": True,
        }

    def digest(self):
        return canonical_digest(self.identity_payload())


@dataclass(frozen=True)
class WholeConeCheckpoint:
    schema: str
    config_digest: str
    profile_payload_sha256: str
    source_binding_sha256: str
    settings_digest: str
    prefix_digest: str
    completed_cells: int
    last_cell_index: int
    last_rho_end: str
    polynomial_sum: tuple
    remainder_sum: tuple
    total_sum: tuple
    residual_integrals: tuple
    cells: tuple
    history_minus_reference_before_norms: bool = True
    physical_rho1_source_error: object = None
    physical_EXISTENCE_certificate: bool = False
    physical_NONEXISTENCE_certificate: bool = False

    def __post_init__(self):
        if self.schema != CHECKPOINT_SCHEMA:
            raise ValueError("unexpected whole-cone checkpoint schema")
        if self.physical_EXISTENCE_certificate or self.physical_NONEXISTENCE_certificate:
            raise ValueError("checkpoint cannot issue a physical gate")
        if self.physical_rho1_source_error is not None:
            raise ValueError("physical rho=1 source error is unresolved")
        if not self.history_minus_reference_before_norms:
            raise ValueError("history-minus-reference must be formed before norms")
        object.__setattr__(self, "cells", tuple(self.cells))
        object.__setattr__(self, "polynomial_sum", tuple(self.polynomial_sum))
        object.__setattr__(self, "remainder_sum", tuple(self.remainder_sum))
        object.__setattr__(self, "total_sum", tuple(self.total_sum))
        object.__setattr__(self, "residual_integrals",
                           tuple(self.residual_integrals))
        if (isinstance(self.completed_cells, bool)
                or not isinstance(self.completed_cells, int)
                or isinstance(self.last_cell_index, bool)
                or not isinstance(self.last_cell_index, int)):
            raise ValueError("integer checkpoint cell counts required")
        for values, name in (
                (self.polynomial_sum, "polynomial sum"),
                (self.remainder_sum, "remainder sum"),
                (self.total_sum, "total sum"),
                (self.residual_integrals, "residual integrals")):
            restored = _upper_pair(values, name)
            if any(not item >= 0 for item in restored):
                raise ValueError("nonnegative checkpoint " + name + " required")
        if self.residual_integrals != self.total_sum:
            raise ValueError("normalized residual integrals must equal normalized totals")
        for value, name in (
                (self.config_digest, "config digest"),
                (self.profile_payload_sha256, "profile payload digest"),
                (self.source_binding_sha256, "source binding digest"),
                (self.settings_digest, "settings digest"),
                (self.prefix_digest, "prefix digest")):
            _sha256_text(value, name)
        recomputed = prefix_digest_from_cells(
            self.config_digest, self.cells,
            polynomial_sum=self.polynomial_sum,
            remainder_sum=self.remainder_sum,
            total_sum=self.total_sum,
            residual_integrals=self.residual_integrals)
        if recomputed != self.prefix_digest:
            raise ValueError("mismatched prefix digest")
        if self.completed_cells != len(self.cells) or self.completed_cells < 1:
            raise ValueError("checkpoint prefix is empty or inconsistent")
        previous = None
        cell_keys = {
            "cell_index", "rho_start", "rho_end", "segment_sha256",
            "polynomial", "remainder", "total", "width", "integral",
        }
        for offset, cell in enumerate(self.cells):
            if (not isinstance(cell, dict)
                    or set(cell) != cell_keys
                    or not isinstance(cell.get("cell_index"), int)
                    or isinstance(cell.get("cell_index"), bool)
                    or cell["cell_index"] < 0
                    or not isinstance(cell.get("segment_sha256"), str)
                    or len(cell["segment_sha256"]) != 64
                    or any(char not in "0123456789abcdef"
                           for char in cell["segment_sha256"])):
                raise ValueError("checkpoint cell identity is malformed")
            try:
                rho_start = float.fromhex(cell["rho_start"])
                rho_end = float.fromhex(cell["rho_end"])
                width = restored_upper(cell["width"])
                pairs = tuple(_upper_pair(cell[name], "cell " + name)
                              for name in ("polynomial", "remainder", "total",
                                           "integral"))
            except (TypeError, ValueError, KeyError) as exc:
                raise ValueError("checkpoint cell bounds are malformed") from exc
            if (not np.isfinite(rho_start) or not np.isfinite(rho_end)
                    or not rho_start > rho_end or not width > 0
                    or any(not item >= 0 for pair in pairs for item in pair)):
                raise ValueError("checkpoint cell bounds are malformed")
            if cell["integral"] != cell["total"]:
                raise ValueError(
                    "normalized residual integral must not receive a second width factor")
            if offset and cell["cell_index"] != self.cells[offset - 1]["cell_index"] + 1:
                raise ValueError("checkpoint cells are not a contiguous prefix")
            if previous is not None and cell.get("rho_start") != previous:
                raise ValueError("checkpoint cell endpoints are not contiguous")
            previous = cell.get("rho_end")
        if (self.last_cell_index != self.cells[-1]["cell_index"]
                or self.last_rho_end != self.cells[-1]["rho_end"]):
            raise ValueError("checkpoint terminal cell is inconsistent")

    def to_mapping(self):
        return {
            "schema": self.schema,
            "config_digest": self.config_digest,
            "profile_payload_sha256": self.profile_payload_sha256,
            "source_binding_sha256": self.source_binding_sha256,
            "settings_digest": self.settings_digest,
            "prefix_digest": self.prefix_digest,
            "completed_cells": self.completed_cells,
            "last_cell_index": self.last_cell_index,
            "last_rho_end": self.last_rho_end,
            "polynomial_sum": deepcopy(list(self.polynomial_sum)),
            "remainder_sum": deepcopy(list(self.remainder_sum)),
            "total_sum": deepcopy(list(self.total_sum)),
            "residual_integrals": deepcopy(list(self.residual_integrals)),
            "cells": deepcopy(list(self.cells)),
            "history_minus_reference_before_norms": True,
            "physical_rho1_source_error": None,
            "physical_EXISTENCE_certificate": False,
            "physical_NONEXISTENCE_certificate": False,
        }

    @classmethod
    def from_mapping(cls, value):
        if not isinstance(value, dict):
            raise TypeError("checkpoint mapping required")
        required = {
            "schema", "config_digest", "profile_payload_sha256",
            "source_binding_sha256", "settings_digest", "prefix_digest",
            "completed_cells", "last_cell_index", "last_rho_end",
            "polynomial_sum", "remainder_sum", "total_sum",
            "residual_integrals", "cells",
            "history_minus_reference_before_norms",
            "physical_rho1_source_error", "physical_EXISTENCE_certificate",
            "physical_NONEXISTENCE_certificate",
        }
        if set(value) != required:
            raise ValueError("malformed whole-cone checkpoint mapping")
        return cls(
            schema=value["schema"],
            config_digest=value["config_digest"],
            profile_payload_sha256=value["profile_payload_sha256"],
            source_binding_sha256=value["source_binding_sha256"],
            settings_digest=value["settings_digest"],
            prefix_digest=value["prefix_digest"],
            completed_cells=value["completed_cells"],
            last_cell_index=value["last_cell_index"],
            last_rho_end=value["last_rho_end"],
            polynomial_sum=tuple(value["polynomial_sum"]),
            remainder_sum=tuple(value["remainder_sum"]),
            total_sum=tuple(value["total_sum"]),
            residual_integrals=tuple(value["residual_integrals"]),
            cells=tuple(value["cells"]),
            history_minus_reference_before_norms=value[
                "history_minus_reference_before_norms"],
            physical_rho1_source_error=value["physical_rho1_source_error"],
            physical_EXISTENCE_certificate=value[
                "physical_EXISTENCE_certificate"],
            physical_NONEXISTENCE_certificate=value[
                "physical_NONEXISTENCE_certificate"],
        )

    def validate(self, config):
        if not isinstance(config, WholeConeFieldConfig):
            raise TypeError("WholeConeFieldConfig required")
        if self.profile_payload_sha256 != config.profile_payload_sha256:
            raise ValueError("mismatched profile payload")
        if self.source_binding_sha256 != config.source_binding_sha256:
            raise ValueError("mismatched source")
        if self.settings_digest != config.settings_digest:
            raise ValueError("mismatched settings")
        if self.config_digest != config.digest():
            raise ValueError("mismatched prefix digest")
        return True


def prefix_digest_from_cells(config_digest, cells, *, polynomial_sum,
                             remainder_sum, total_sum, residual_integrals):
    return canonical_digest({
        "config_digest": config_digest,
        "cells": list(cells),
        "polynomial_sum": list(polynomial_sum),
        "remainder_sum": list(remainder_sum),
        "total_sum": list(total_sum),
        "residual_integrals": list(residual_integrals),
    })


def _segment_sha256(segment):
    digest = sha256()
    digest.update(float(segment.rho_start).hex().encode("ascii"))
    digest.update(b"\0")
    digest.update(float(segment.rho_end).hex().encode("ascii"))
    for array in (segment.start, segment.end, segment.coefficients):
        value = np.ascontiguousarray(array, dtype="<c16")
        digest.update(str(value.shape).encode("ascii"))
        digest.update(b"\0")
        digest.update(value.tobytes())
    return digest.hexdigest()


def segment_sha256(segment):
    """Dense-segment identity used by checkpoint replay."""
    if not isinstance(segment, TrajectorySegment):
        raise TypeError("saved anchored TrajectorySegment required")
    return _segment_sha256(segment)


def segment_identity(segment, cell_index):
    if isinstance(cell_index, bool) or not isinstance(cell_index, int) or cell_index < 0:
        raise ValueError("nonnegative cell index required")
    return {
        "cell_index": int(cell_index),
        "rho_start": float(segment.rho_start).hex(),
        "rho_end": float(segment.rho_end).hex(),
        "segment_sha256": _segment_sha256(segment),
    }


def _cell_record(cell_index, segment, enclosed):
    return {
        "cell_index": int(cell_index),
        "rho_start": float(segment.rho_start).hex(),
        "rho_end": float(segment.rho_end).hex(),
        "segment_sha256": _segment_sha256(segment),
        "polynomial": _pack_pair(enclosed["polynomial"]),
        "remainder": _pack_pair(enclosed["remainder"]),
        "total": _pack_pair(enclosed["total"]),
        "width": exact_upper(enclosed["width"]),
        "integral": _pack_pair(enclosed["integral"]),
    }


@dataclass(frozen=True)
class WholeConeFieldResult:
    schema: str
    family: tuple
    profile_identity: str
    completed_cells: int
    bounds: object
    coverage: object
    status: str
    propagation: object = None
    matter: object = None
    source_accuracy_included: bool = False
    physical_rho1_source_error: object = None
    physical_EXISTENCE_certificate: bool = False
    physical_NONEXISTENCE_certificate: bool = False
    family_14_1_scientific_run: bool = False
    history_minus_reference_before_norms: bool = True
    prefix_digest: str = ""
    config_digest: str = ""

    def __post_init__(self):
        if self.schema != SCHEMA:
            raise ValueError("unexpected whole-cone result schema")
        if (self.source_accuracy_included
                or self.physical_rho1_source_error is not None
                or self.physical_EXISTENCE_certificate
                or self.physical_NONEXISTENCE_certificate
                or self.family_14_1_scientific_run):
            raise ValueError("whole-cone core v1 cannot issue physical evidence")
        if not self.history_minus_reference_before_norms:
            raise ValueError("history-minus-reference must be formed before norms")
        if (not isinstance(self.completed_cells, int)
                or isinstance(self.completed_cells, bool)
                or self.completed_cells < 1):
            raise ValueError("positive completed cell count required")
        if not isinstance(self.status, str) or not self.status.startswith("OPEN"):
            raise ValueError("whole-cone core result must remain OPEN")
        if not isinstance(self.bounds, Mapping) or not isinstance(self.coverage, Mapping):
            raise TypeError("whole-cone bounds and coverage mappings required")
        _sha256_text(self.profile_identity, "profile identity")
        _sha256_text(self.prefix_digest, "prefix digest")
        _sha256_text(self.config_digest, "config digest")

    def to_mapping(self):
        return {
            "schema": self.schema,
            "family": list(self.family),
            "profile_identity": self.profile_identity,
            "completed_cells": self.completed_cells,
            "bounds": dict(self.bounds),
            "coverage": dict(self.coverage),
            "status": self.status,
            "propagation": None if self.propagation is None else dict(self.propagation),
            "matter": None if self.matter is None else dict(self.matter),
            "source_accuracy_included": False,
            "physical_rho1_source_error": None,
            "physical_EXISTENCE_certificate": False,
            "physical_NONEXISTENCE_certificate": False,
            "family_14_1_scientific_run": False,
            "history_minus_reference_before_norms": True,
            "prefix_digest": self.prefix_digest,
            "config_digest": self.config_digest,
            "certificate_use": False,
        }


class WholeConeAccumulator:
    """Streamed directed totals; segments and field balls are not retained."""

    def __init__(self, config, profile_payload, family):
        if not isinstance(config, WholeConeFieldConfig):
            raise TypeError("WholeConeFieldConfig required")
        if not isinstance(family, LocalIncomingFamily):
            raise TypeError("LocalIncomingFamily required")
        identity = profile_identity(family, include_normal_window=True)
        if identity != config.profile_identity:
            raise ValueError("mismatched source")
        if (not isinstance(profile_payload, dict) or not profile_payload
                or not all(isinstance(value, dict)
                           for value in profile_payload.values())):
            raise TypeError("serialized local Fourier payload required")
        restored = restore_local_profile_fourier(
            profile_payload, expected_sha256=config.profile_payload_sha256)
        if any(len(row.uniform_tail_bounds) != 3 for row in restored.values()):
            raise ValueError("whole-cone profile requires derivative tails through order 2")
        rows = profile_rows_from_bounds(restored)
        keys = monomial_keys(config.reciprocal_order)
        if set(rows) != set(keys):
            raise ValueError("profile keys differ from the reciprocal keys")
        self._config = config
        self._family = family
        self._profile_rows = rows
        self._bits = config.bits
        self._cells = []
        zeros = _pack_pair(_zero_pair(self._bits))
        self._polynomial_sum = list(zeros)
        self._remainder_sum = list(zeros)
        self._total_sum = list(zeros)
        self._integral_sum = list(zeros)
        self._last_rho_end = None
        self._last_cell_index = None
        self._finalized = False
        self._resume_cells_pending = deque()

    def _add_packed(self, packed, values):
        with ctx.workprec(self._bits):
            return [
                exact_upper((restored_upper(left) + right).upper())
                for left, right in zip(packed, values)
            ]

    def _add_enclosed(self, enclosed):
        self._polynomial_sum = self._add_packed(
            self._polynomial_sum, enclosed["polynomial"])
        self._remainder_sum = self._add_packed(
            self._remainder_sum, enclosed["remainder"])
        self._total_sum = self._add_packed(
            self._total_sum, enclosed["total"])
        self._integral_sum = self._add_packed(
            self._integral_sum, enclosed["integral"])

    def accumulate_segment(self, segment, *, cell_index):
        if self._finalized:
            raise ValueError("cannot accumulate after finalize_family")
        if self._resume_cells_pending:
            raise ValueError("checkpoint prefix replay required")
        if not isinstance(segment, TrajectorySegment):
            raise TypeError("saved anchored TrajectorySegment required")
        if isinstance(cell_index, bool) or not isinstance(cell_index, int):
            raise ValueError("integer cell index required")
        if (self._last_cell_index is not None
                and cell_index != self._last_cell_index + 1):
            raise ValueError("segment index is not the next prefix cell")
        if (self._last_rho_end is not None
                and float(segment.rho_start).hex() != self._last_rho_end):
            raise ValueError("segment replay")
        enclosed = enclose_difference_segment(
            self._config, self._family, self._profile_rows,
            self._config.radius_tail, segment)
        # Packing outside the bound's precision can round an exact dyadic
        # endpoint down when the ambient Arb context is only binary64-sized.
        with ctx.workprec(self._bits):
            record = _cell_record(cell_index, segment, enclosed)
        self._add_enclosed(enclosed)
        self._cells.append(record)
        self._last_rho_end = record["rho_end"]
        self._last_cell_index = cell_index
        return record

    def replay_segment(self, segment, *, cell_index):
        """Authenticate the next replayed dense segment without rebounding it."""
        if not isinstance(segment, TrajectorySegment):
            raise TypeError("saved anchored TrajectorySegment required")
        if isinstance(cell_index, bool) or not isinstance(cell_index, int):
            raise ValueError("integer cell index required")
        if not self._resume_cells_pending:
            raise ValueError("segment replay")
        stored = self._resume_cells_pending[0]
        identity = {
            "cell_index": int(cell_index),
            "rho_start": float(segment.rho_start).hex(),
            "rho_end": float(segment.rho_end).hex(),
            "segment_sha256": _segment_sha256(segment),
        }
        if any(stored.get(name) != value for name, value in identity.items()):
            raise ValueError("segment replay")
        self._resume_cells_pending.popleft()
        return stored

    def checkpoint(self):
        if not self._cells:
            raise ValueError("checkpoint requires an accumulated prefix")
        config_digest = self._config.digest()
        mapping = {
            "schema": CHECKPOINT_SCHEMA,
            "config_digest": config_digest,
            "profile_payload_sha256": self._config.profile_payload_sha256,
            "source_binding_sha256": self._config.source_binding_sha256,
            "settings_digest": self._config.settings_digest,
            "prefix_digest": prefix_digest_from_cells(
                config_digest, self._cells,
                polynomial_sum=self._polynomial_sum,
                remainder_sum=self._remainder_sum,
                total_sum=self._total_sum,
                residual_integrals=self._integral_sum),
            "completed_cells": len(self._cells),
            "last_cell_index": self._last_cell_index,
            "last_rho_end": self._last_rho_end,
            "polynomial_sum": deepcopy(list(self._polynomial_sum)),
            "remainder_sum": deepcopy(list(self._remainder_sum)),
            "total_sum": deepcopy(list(self._total_sum)),
            "residual_integrals": deepcopy(list(self._integral_sum)),
            "cells": deepcopy(list(self._cells)),
            "history_minus_reference_before_norms": True,
            "physical_rho1_source_error": None,
            "physical_EXISTENCE_certificate": False,
            "physical_NONEXISTENCE_certificate": False,
        }
        return WholeConeCheckpoint.from_mapping(mapping)

    @classmethod
    def from_checkpoint(cls, checkpoint, config, profile_payload, family):
        if isinstance(checkpoint, dict):
            checkpoint = WholeConeCheckpoint.from_mapping(checkpoint)
        if not isinstance(checkpoint, WholeConeCheckpoint):
            raise TypeError("WholeConeCheckpoint required")
        checkpoint.validate(config)
        accumulator = cls(config, profile_payload, family)
        accumulator._cells = deepcopy(list(checkpoint.cells))
        accumulator._polynomial_sum = deepcopy(list(checkpoint.polynomial_sum))
        accumulator._remainder_sum = deepcopy(list(checkpoint.remainder_sum))
        accumulator._total_sum = deepcopy(list(checkpoint.total_sum))
        accumulator._integral_sum = deepcopy(list(checkpoint.residual_integrals))
        accumulator._last_rho_end = checkpoint.last_rho_end
        accumulator._last_cell_index = checkpoint.last_cell_index
        accumulator._resume_cells_pending = deque(deepcopy(list(checkpoint.cells)))
        return accumulator

    def finalize_family(self, *, propagation=None, matter=None,
                        source_error=None, source_accuracy_included=False):
        if self._finalized:
            raise ValueError("family already finalized")
        if not self._cells:
            raise ValueError("at least one accumulated segment required")
        if self._resume_cells_pending:
            raise ValueError("checkpoint prefix replay required")
        if source_accuracy_included or source_error is not None:
            raise ValueError("physical rho=1 source error is unresolved")
        self._finalized = True
        bounds = {
            "polynomial_sum": list(self._polynomial_sum),
            "remainder_sum": list(self._remainder_sum),
            "total_continuous_normalized_residual_sum": list(self._total_sum),
            "residual_integrals": list(self._integral_sum),
            "last_cell_total": self._cells[-1]["total"],
        }
        with ctx.workprec(self._bits):
            integrals = _upper_pair(self._integral_sum, "residual integrals")
        propagated = None
        if propagation is not None:
            propagated = whole_cone_propagate_difference_error(
                integrals, **dict(propagation))
        contracted = None
        if matter is not None:
            contracted = whole_cone_difference_matter_error(**dict(matter))
        return WholeConeFieldResult(
            schema=SCHEMA,
            family=self._config.family,
            profile_identity=self._config.profile_identity,
            completed_cells=len(self._cells),
            bounds=bounds,
            coverage={
                "completed_cells": len(self._cells),
                "all_history_cells": False,
                "all_source_families": False,
                "family_14_1_scientific_run": False,
                "source_columns": len(self._config.source_weights),
                "whole_spatial_period": True,
                "physical_rho1_source_error": None,
            },
            status=(
                "OPEN: whole-cone field core accumulated; physical rho=1 "
                "source error unresolved"),
            propagation=propagated,
            matter=contracted,
            prefix_digest=prefix_digest_from_cells(
                self._config.digest(), self._cells,
                polynomial_sum=self._polynomial_sum,
                remainder_sum=self._remainder_sum,
                total_sum=self._total_sum,
                residual_integrals=self._integral_sum),
            config_digest=self._config.digest(),
        )


def accumulate_segment(accumulator, segment, *, cell_index):
    if not isinstance(accumulator, WholeConeAccumulator):
        raise TypeError("WholeConeAccumulator required")
    return accumulator.accumulate_segment(segment, cell_index=cell_index)


def finalize_family(accumulator, **kwargs):
    if not isinstance(accumulator, WholeConeAccumulator):
        raise TypeError("WholeConeAccumulator required")
    return accumulator.finalize_family(**kwargs)
