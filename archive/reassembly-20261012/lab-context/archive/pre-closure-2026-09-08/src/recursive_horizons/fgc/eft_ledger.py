"""Fail-closed exact-rational core for the FGC-1 retained-EFT ledger.

This module intentionally evaluates only a *declared* local power-counting
ledger.  It cannot infer a cutoff, a Wilsonian operator basis, Wilson
coefficients, a frequency bound, or EFT validity from FGC-QR fixture units.
All numeric values accepted here are exact canonical rational strings.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
import tomllib
from typing import Any, Mapping


Q = Fraction
EFT0_LED1_ARTIFACT_ID = "FGC-1-EFT0-LED1"
FGC_PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = Path("configs/fgc/fgc-1-eft0-led1.toml")


def _mapping(name: str, value: object) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a table")
    return value


def _keys(name: str, value: Mapping[str, object], expected: set[str]) -> None:
    actual = set(value)
    if actual != expected:
        raise ValueError(f"{name} keys differ: expected {sorted(expected)}, got {sorted(actual)}")


def _text(name: str, value: object) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a nonempty string")
    return value


def _boolean(name: str, value: object, *, required: bool | None = None) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{name} must be boolean")
    if required is not None and value is not required:
        raise ValueError(f"{name} must be {str(required).lower()}")
    return value


def _rational(name: str, value: object, *, positive: bool = False, nonnegative: bool = False) -> Fraction:
    text = _text(name, value)
    if any(character in text.lower() for character in (".", "e")):
        raise ValueError(f"{name} must be an exact rational string, never a float")
    try:
        result = Q(text)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{name} must be an exact rational string") from exc
    if str(result) != text:
        raise ValueError(f"{name} must use canonical rational spelling")
    if positive and result <= 0:
        raise ValueError(f"{name} must be positive")
    if nonnegative and result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def _no_floats(value: object) -> None:
    if isinstance(value, float):
        raise ValueError("configuration contains a float; exact rational strings are required")
    if isinstance(value, Mapping):
        for nested in value.values():
            _no_floats(nested)
    elif isinstance(value, list):
        for nested in value:
            _no_floats(nested)


@dataclass(frozen=True, slots=True)
class EFTLedgerConfig:
    """Frozen inputs for a conditional, local retained-EFT diagnostic."""

    artifact_id: str
    project_version: str
    cutoff: Fraction
    cutoff_origin: str
    planck_mass: Fraction
    scalar_mass: Fraction
    quartic_coupling: Fraction
    beta: Fraction
    eta: Fraction
    phi_center: Fraction
    phi_half_width: Fraction
    phi_radial_derivative_center: Fraction
    phi_radial_derivative_half_width: Fraction
    retained_operators: tuple[Mapping[str, object], ...]
    omitted_operators: tuple[Mapping[str, object], ...]
    controls: Mapping[str, tuple[str, Fraction]]


_RETAINED_IDS = (
    "einstein_hilbert",
    "phi2_ricci",
    "phi_kinetic",
    "chi_kinetic",
    "phi_mass",
    "phi_quartic",
    "phi2_gauss_bonnet",
)
_OMITTED_IDS = ("x_phi_squared", "x_chi_squared", "x_phi_x_chi", "phi6")
_CONTROL_IDS = ("F_over_Mpl_squared_lower", "epsilon_F", "epsilon_eta_phi2", "epsilon_phi_amplitude")
_RETAINED_REFERENCES = {
    "einstein_hilbert": ("M_Pl",),
    "phi2_ricci": ("beta",),
    "phi_kinetic": (),
    "chi_kinetic": (),
    "phi_mass": ("mu",),
    "phi_quartic": ("g4",),
    "phi2_gauss_bonnet": ("eta",),
}


def _validate_operator_table(
    table: object, *, name: str, expected_ids: tuple[str, ...], omitted: bool
) -> tuple[Mapping[str, object], ...]:
    if not isinstance(table, list) or len(table) != len(expected_ids):
        raise ValueError(f"{name} must contain exactly the declared operator set")
    validated: list[Mapping[str, object]] = []
    for index, entry in enumerate(table):
        item = _mapping(f"{name}[{index}]", entry)
        expected = {"id", "operator", "canonical_dimension", "coefficient_dimension", "coefficient", "status"}
        if omitted:
            expected |= {"wilson_bound", "suppression_power", "basis_scope"}
        else:
            expected |= {"parameter_references"}
        _keys(f"{name}[{index}]", item, expected)
        if _text(f"{name}[{index}].id", item["id"]) != expected_ids[index]:
            raise ValueError(f"{name} has a noncanonical operator order")
        _text(f"{name}[{index}].operator", item["operator"])
        dimension = _rational(f"{name}[{index}].canonical_dimension", item["canonical_dimension"], nonnegative=True)
        coefficient_dimension = _rational(
            f"{name}[{index}].coefficient_dimension", item["coefficient_dimension"]
        )
        if dimension.denominator != 1 or coefficient_dimension != Q(4) - dimension:
            raise ValueError(f"{name}[{index}] has inconsistent canonical dimensions")
        expected_status = "omitted_representative_only" if omitted else "retained"
        if _text(f"{name}[{index}].status", item["status"]) != expected_status:
            raise ValueError(f"{name}[{index}].status is not fail-closed")
        _text(f"{name}[{index}].coefficient", item["coefficient"])
        if not omitted:
            references = item["parameter_references"]
            if not isinstance(references, list) or tuple(references) != _RETAINED_REFERENCES[expected_ids[index]]:
                raise ValueError(f"{name}[{index}].parameter_references are noncanonical")
        if omitted:
            _rational(f"{name}[{index}].wilson_bound", item["wilson_bound"], positive=True)
            suppression = _rational(f"{name}[{index}].suppression_power", item["suppression_power"], nonnegative=True)
            if suppression.denominator != 1 or suppression != dimension - Q(4):
                raise ValueError(f"{name}[{index}].suppression_power is dimensionally inconsistent")
            if _text(f"{name}[{index}].basis_scope", item["basis_scope"]) != "representative_only":
                raise ValueError(f"{name}[{index}].basis_scope must remain representative_only")
        validated.append(item)
    return tuple(validated)


def load_eft_ledger_config(path: Path = DEFAULT_CONFIG) -> EFTLedgerConfig:
    """Parse the strict EFT0 declaration without accepting inferred inputs."""

    source = Path(path)
    raw: object = tomllib.loads(source.read_text(encoding="utf-8"))
    _no_floats(raw)
    root = _mapping("root", raw)
    _keys(
        "root",
        root,
        {
            "schema_version", "artifact_id", "project_version", "units", "frame", "model",
            "cutoff", "operator_basis", "local_comp1_centered_component_box",
            "retained_operators", "omitted_operators", "controls", "claims",
        },
    )
    if isinstance(root["schema_version"], bool) or root["schema_version"] != 1:
        raise ValueError("schema_version must equal 1")
    if _text("artifact_id", root["artifact_id"]) != EFT0_LED1_ARTIFACT_ID:
        raise ValueError(f"artifact_id must equal {EFT0_LED1_ARTIFACT_ID}")
    if _text("project_version", root["project_version"]) != FGC_PROJECT_VERSION:
        raise ValueError(f"project_version must equal {FGC_PROJECT_VERSION}")
    if _text("units", root["units"]) != "c=hbar=1; four spacetime dimensions":
        raise ValueError("units must declare the frozen natural-unit convention")
    if _text("frame", root["frame"]) != "Jordan frame; physical metric g_ab; signature -+++":
        raise ValueError("frame must declare the frozen Jordan-frame convention")

    model = _mapping("model", root["model"])
    _keys("model", model, {"id", "planck_mass", "mu", "g4", "beta", "eta", "fixture_values_are_not_cutoff"})
    if _text("model.id", model["id"]) != "FGC-QR":
        raise ValueError("model.id must equal FGC-QR")
    planck_mass = _rational("model.planck_mass", model["planck_mass"], positive=True)
    scalar_mass = _rational("model.mu", model["mu"], positive=True)
    quartic_coupling = _rational("model.g4", model["g4"], positive=True)
    beta = _rational("model.beta", model["beta"])
    eta = _rational("model.eta", model["eta"], positive=True)
    _boolean("model.fixture_values_are_not_cutoff", model["fixture_values_are_not_cutoff"], required=True)
    if (planck_mass, scalar_mass, quartic_coupling, beta, eta) != (Q(2), Q(3), Q(1, 2), -Q(1, 4), Q(1, 2)):
        raise ValueError("model parameters must match the frozen FGC-QR action fixture")

    cutoff = _mapping("cutoff", root["cutoff"])
    _keys("cutoff", cutoff, {"Lambda", "origin", "external_assumption", "inferred_from_fixture"})
    cutoff_value = _rational("cutoff.Lambda", cutoff["Lambda"], positive=True)
    cutoff_origin = _text("cutoff.origin", cutoff["origin"])
    _boolean("cutoff.external_assumption", cutoff["external_assumption"], required=True)
    _boolean("cutoff.inferred_from_fixture", cutoff["inferred_from_fixture"], required=False)
    if cutoff_origin == "fixture" or cutoff_origin == "inferred_from_fixture":
        raise ValueError("cutoff must not be inferred from fixture values")

    basis = _mapping("operator_basis", root["operator_basis"])
    _keys("operator_basis", basis, {"basis_scope", "field_redefinition_convention", "phi_symmetry", "chi_symmetry", "complete_wilsonian_basis"})
    if _text("operator_basis.basis_scope", basis["basis_scope"]) != "representative_only":
        raise ValueError("operator_basis must remain representative_only")
    _text("operator_basis.field_redefinition_convention", basis["field_redefinition_convention"])
    _text("operator_basis.phi_symmetry", basis["phi_symmetry"])
    _text("operator_basis.chi_symmetry", basis["chi_symmetry"])
    _boolean("operator_basis.complete_wilsonian_basis", basis["complete_wilsonian_basis"], required=False)

    local_box = _mapping("local_comp1_centered_component_box", root["local_comp1_centered_component_box"])
    _keys(local_box_name := "local_comp1_centered_component_box", local_box, {"source_artifact", "qift1_source_artifact", "qift1_parameter_half_width", "coordinate_radius", "phi_center", "phi_half_width", "phi_radial_derivative_center", "phi_radial_derivative_half_width", "full_jet_box", "frequency_bound_available"})
    if _text(f"{local_box_name}.source_artifact", local_box["source_artifact"]) != "FGC-1-HYP1-CON1-COMP1":
        raise ValueError("local component box must be centered on COMP1")
    if _rational(f"{local_box_name}.coordinate_radius", local_box["coordinate_radius"], positive=True) != Q(4):
        raise ValueError("local component box coordinate radius must equal COMP1's 4")
    phi_center = _rational(f"{local_box_name}.phi_center", local_box["phi_center"])
    phi_half_width = _rational(f"{local_box_name}.phi_half_width", local_box["phi_half_width"], positive=True)
    phi_dr_center = _rational(f"{local_box_name}.phi_radial_derivative_center", local_box["phi_radial_derivative_center"])
    phi_dr_half_width = _rational(f"{local_box_name}.phi_radial_derivative_half_width", local_box["phi_radial_derivative_half_width"], positive=True)
    if phi_center != Q(1, 131072) or phi_dr_center != Q(1, 131072):
        raise ValueError("local component box must use COMP1's activated phi and radial derivative center")
    _boolean(f"{local_box_name}.full_jet_box", local_box["full_jet_box"], required=False)
    _boolean(f"{local_box_name}.frequency_bound_available", local_box["frequency_bound_available"], required=False)
    if _text(f"{local_box_name}.qift1_source_artifact", local_box["qift1_source_artifact"]) != "FGC-1-HYP1-DOM1-QIFT1":
        raise ValueError("local component box must declare QIFT1 predecessor domain")
    qift_half_width = _rational(f"{local_box_name}.qift1_parameter_half_width", local_box["qift1_parameter_half_width"], positive=True)
    if qift_half_width != Q(1, 65536):
        raise ValueError("local component box must use QIFT1's frozen parameter half width")
    if phi_center - phi_half_width <= 0 or phi_dr_center - phi_dr_half_width <= 0:
        raise ValueError("component box must keep phi and phi_r strictly positive")
    if abs(phi_center) + phi_half_width >= qift_half_width or abs(phi_dr_center) + phi_dr_half_width >= qift_half_width:
        raise ValueError("component box must be a strict subbox of QIFT1")

    retained = _validate_operator_table(root["retained_operators"], name="retained_operators", expected_ids=_RETAINED_IDS, omitted=False)
    omitted = _validate_operator_table(root["omitted_operators"], name="omitted_operators", expected_ids=_OMITTED_IDS, omitted=True)

    controls = root["controls"]
    if not isinstance(controls, list) or len(controls) != len(_CONTROL_IDS):
        raise ValueError("controls must contain exactly the declared local controls")
    controls_by_id: dict[str, tuple[str, Fraction]] = {}
    for index, raw_control in enumerate(controls):
        control = _mapping(f"controls[{index}]", raw_control)
        expected_control_keys = {"id", "definition", "comparison", "computable_on_component_box"}
        control_id = _text(f"controls[{index}].id", control["id"])
        if control_id != _CONTROL_IDS[index]:
            raise ValueError("controls have a noncanonical order")
        _text(f"controls[{index}].definition", control["definition"])
        comparison = _text(f"controls[{index}].comparison", control["comparison"])
        if control_id == "F_over_Mpl_squared_lower":
            expected_control_keys.add("strict_lower_bound")
            expected_comparison = "greater_than"
            bound_key = "strict_lower_bound"
        else:
            expected_control_keys.add("strict_upper_bound")
            expected_comparison = "less_than"
            bound_key = "strict_upper_bound"
        _keys(f"controls[{index}]", control, expected_control_keys)
        if comparison != expected_comparison:
            raise ValueError(f"controls[{index}].comparison must be {expected_comparison}")
        threshold = _rational(f"controls[{index}].{bound_key}", control[bound_key], positive=True)
        if control_id == "F_over_Mpl_squared_lower":
            if threshold >= 1:
                raise ValueError("normalized F lower bound must lie strictly below one")
        elif threshold >= 1:
            raise ValueError("epsilon control threshold must lie strictly below one")
        _boolean(f"controls[{index}].computable_on_component_box", control["computable_on_component_box"], required=True)
        controls_by_id[control_id] = (comparison, threshold)

    claims = _mapping("claims", root["claims"])
    _keys("claims", claims, {"retained_eft_validity", "complete_wilsonian_basis", "uv_completion", "wilson_coefficients_predicted", "radiative_stability", "frequency_control", "global_or_open_eft_validity", "initial_data", "evolution", "collapse", "defocusing", "singularity_resolution"})
    for key, value in claims.items():
        _boolean(f"claims.{key}", value, required=False)

    return EFTLedgerConfig(
        artifact_id=EFT0_LED1_ARTIFACT_ID,
        project_version=FGC_PROJECT_VERSION,
        cutoff=cutoff_value,
        cutoff_origin=cutoff_origin,
        planck_mass=planck_mass,
        scalar_mass=scalar_mass,
        quartic_coupling=quartic_coupling,
        beta=beta,
        eta=eta,
        phi_center=phi_center,
        phi_half_width=phi_half_width,
        phi_radial_derivative_center=phi_dr_center,
        phi_radial_derivative_half_width=phi_dr_half_width,
        retained_operators=retained,
        omitted_operators=omitted,
        controls=controls_by_id,
    )


def _fraction_text(value: Fraction) -> str:
    return str(value)


def local_eft_ledger_certificate(config: EFTLedgerConfig) -> dict[str, Any]:
    """Evaluate only controls determined by the declared scalar component box.

    The φ and ∂rφ intervals are coordinate-component inputs, not a
    covariant derivative norm or a frequency bound.  Hence the derivative is
    serialized for provenance but is deliberately not promoted to a control.
    """

    if not isinstance(config, EFTLedgerConfig):
        raise TypeError("config must be an EFTLedgerConfig")
    phi_abs_upper = abs(config.phi_center) + config.phi_half_width
    phi_dr_abs_upper = abs(config.phi_radial_derivative_center) + config.phi_radial_derivative_half_width
    F_lower = config.planck_mass**2 - abs(config.beta) * phi_abs_upper**2
    F_over_Mpl_squared_lower = F_lower / config.planck_mass**2
    epsilon_F = abs(config.beta) * phi_abs_upper**2 / config.planck_mass**2
    epsilon_eta_phi2 = config.eta * phi_abs_upper**2
    epsilon_phi_amplitude = phi_abs_upper / config.cutoff
    values = {
        "F_over_Mpl_squared_lower": F_over_Mpl_squared_lower,
        "epsilon_F": epsilon_F,
        "epsilon_eta_phi2": epsilon_eta_phi2,
        "epsilon_phi_amplitude": epsilon_phi_amplitude,
    }
    passes = {
        "F_over_Mpl_squared_lower": F_over_Mpl_squared_lower > config.controls["F_over_Mpl_squared_lower"][1],
        "epsilon_F": epsilon_F < config.controls["epsilon_F"][1],
        "epsilon_eta_phi2": epsilon_eta_phi2 < config.controls["epsilon_eta_phi2"][1],
        "epsilon_phi_amplitude": epsilon_phi_amplitude < config.controls["epsilon_phi_amplitude"][1],
    }
    if not all(passes.values()):
        raise ValueError("a declared local EFT control crosses its strict threshold")
    return {
        "artifact_id": config.artifact_id,
        "classification": "conditional_declared_local_retained_eft_ledger_only",
        "exact_arithmetic": True,
        "cutoff": {"Lambda": _fraction_text(config.cutoff), "origin": config.cutoff_origin, "external_assumption": True, "inferred_from_fixture": False},
        "component_box": {
            "source_artifact": "FGC-1-HYP1-CON1-COMP1",
            "phi_center": _fraction_text(config.phi_center),
            "phi_half_width": _fraction_text(config.phi_half_width),
            "phi_absolute_upper": _fraction_text(phi_abs_upper),
            "phi_radial_derivative_center": _fraction_text(config.phi_radial_derivative_center),
            "phi_radial_derivative_half_width": _fraction_text(config.phi_radial_derivative_half_width),
            "phi_radial_derivative_absolute_upper": _fraction_text(phi_dr_abs_upper),
            "F_lower": _fraction_text(F_lower),
            "F_lower_over_M_Pl_squared": _fraction_text(F_over_Mpl_squared_lower),
            "full_jet_box": False,
        },
        "local_controls": {
            key: {
                "value": _fraction_text(value),
                "comparison": config.controls[key][0],
                ("strict_lower_bound" if key == "F_over_Mpl_squared_lower" else "strict_upper_bound"): _fraction_text(config.controls[key][1]),
                "passes": passes[key],
            }
            for key, value in values.items()
        },
        "frequency_control": {"available": False, "reason": "a finite local component/jet box does not bound proper frequency or Fourier support"},
        "action_parameters": {"M_Pl": _fraction_text(config.planck_mass), "mu": _fraction_text(config.scalar_mass), "g4": _fraction_text(config.quartic_coupling), "beta": _fraction_text(config.beta), "eta": _fraction_text(config.eta)},
        "operator_basis": {"scope": "representative_only", "complete_wilsonian_basis": False, "retained_operators": [dict(item) for item in config.retained_operators], "representative_omitted_operators": [dict(item) for item in config.omitted_operators], "representative_omitted_remainders_locally_evaluated": False, "representative_omitted_remainders_reason": "the component-only box lacks the covariant derivative, curvature, and frequency data needed for omitted-operator remainder bounds"},
        "retained_eft_validity": False,
        "nonclaims": {
            "uv_completion": False,
            "complete_wilsonian_basis": False,
            "wilson_coefficients_predicted": False,
            "radiative_stability": False,
            "frequency_control": False,
            "global_or_open_eft_validity": False,
            "initial_data": False,
            "evolution": False,
            "collapse": False,
            "defocusing": False,
            "singularity_resolution": False,
        },
    }
