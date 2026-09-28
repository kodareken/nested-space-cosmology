"""Candidate-blind DEF1-STAB1 conversion-ownership contracts.

This layer makes the missing pre-FRZ1/PREF1 conversion ownership explicit. It
does not freeze FRZ1/PREF1, read a trajectory, set ``DEF1_error_map_passed``,
or authenticate future COL1 values. Richardson and IMP1 remain conditional.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from numbers import Integral
from pathlib import Path
from typing import Final, Mapping, Sequence

from recursive_horizons.evidence_io import (
    UnsafePathError,
    canonical_json_bytes,
    read_regular_file,
)

from .constraint_system import ACTIVE_SPHERICAL_GAUGE_COMPONENTS
from .def1_geometry_error import INPUT_NAMES, JET_COMPONENTS, METRIC_FIELDS
from .def1_stab1 import (
    DEF1_BOOLEAN_NAMES,
    ERROR_BUDGET_COMPONENTS,
    PREMISE_STATUS_CONDITIONAL,
    PREMISE_STATUS_PROVEN,
    SOURCE_EXACT_ALGEBRA,
    SOURCE_IMP1_ADMISSION_DEBIT,
    SOURCE_RICHARDSON,
    SOURCE_SUPPLIED_CERTIFIED,
    Def1Stab1Error,
    inverse_base_metric_from_adm,
    assess_mass_flux_ledger,
)
from .def1_stab1_providers import (
    FULL_RESIDUAL_LENGTH,
    METRIC_FIRST_DERIVATIVE_INPUTS,
    Def1Stab1ProviderError,
    ProviderRecord,
    affine_provider,
    affine_residual_times_step_refused,
    arithmetic_provider,
    boundary_provider,
    complete_zero_radii,
    conservation_provider,
    extraction_provider,
    gauge_constraint_provider,
    initial_data_provider,
    interpolation_provider,
    nonlinear_source_provider,
    nullness_provider,
    physical_constraint_provider,
    reduction_constraint_provider,
    spatial_temporal_provider,
    trajectory_alignment_provider,
)
from .evolution.tdg6_temporal_admission_design import (
    TDG6_COMPLETE_STATE_CHANNELS,
)
from .modified_harmonic import MHG_EQUATION_ORDER, MHG_GAUGE_CONSTRAINT_ORDER
from .modified_harmonic_constraints import PHYSICAL_PROJECTION_ORDER
from .modified_harmonic_first_order import (
    FO1_FIELD_ORDER,
    FO1_PHYSICAL_ARGUMENT_GROUPS,
    FO1_PHYSICAL_ARGUMENT_ORDER,
    FO1_PHYSICAL_EQUATION_ORDER,
    FO1_REDUCTION_CONSTRAINT_ORDER,
    FO1_STATE_ORDER,
)
from .modified_harmonic_implicit import IMPLICIT_EQUATION_ORDER
from .modified_harmonic_reference import MHG2_FULL_EQUATION_ORDER
from .regular_center import REGULAR_EQUATION_ORDER
from .sgb1_ctl1_source import SOURCE_EQUATION_ORDER, SOURCE_FIELD_ORDER
from .spherical_reduction import (
    ADM_FIELD_ORDER,
    BASE_FIELD_ORDER,
    INDEPENDENT_EQUATION_ORDER,
    Jet2,
    SECOND_DERIVATIVE_ORDER,
    state_from_generalized_adm_pg_fixture,
)


Q = Fraction
_REPO_ROOT: Final[Path] = Path(__file__).resolve().parents[3]
TANGENT_SLOTS: Final[tuple[str, ...]] = ("k.t", "k.r")
ADM_METRIC_TWO_JET_ORDER: Final[tuple[str, ...]] = tuple(
    name for name in INPUT_NAMES if name not in TANGENT_SLOTS
)
FO1_GROUP_TO_JET_SLOT: Final[tuple[tuple[str, str], ...]] = (
    ("u", "value"),
    ("p", "dt"),
    ("q", "dr"),
    ("p_t", "dtt"),
    ("p_r", "dtr"),
    ("q_r", "drr"),
)
RED1_TO_MHG2_EQUATION_PAIRS: Final[tuple[tuple[str, str], ...]] = (
    ("metric_tt", "metric_tt_mhg"),
    ("metric_tr", "metric_tr_mhg"),
    ("metric_rr", "metric_rr_mhg"),
    ("metric_theta_theta", "metric_theta_theta_mhg"),
    ("scalar_phi", "scalar_phi"),
    ("scalar_chi", "scalar_chi"),
)
IMP1_CHANNEL_BLOCKS: Final[tuple[str, ...]] = ("u", "p", "q")
IMP1_CHANNEL_FIELDS: Final[tuple[str, ...]] = (
    "alpha",
    "v",
    "lambda",
    "R",
    "phi",
    "chi",
)
IMP1_18_CHANNEL_ORDER: Final[tuple[str, ...]] = tuple(
    f"{block}:{field}"
    for block in IMP1_CHANNEL_BLOCKS
    for field in IMP1_CHANNEL_FIELDS
)
ADM_SOURCE_PHYSICAL_ARGUMENT_ORDER: Final[tuple[str, ...]] = tuple(
    f"{group}.{field}"
    for group in FO1_PHYSICAL_ARGUMENT_GROUPS
    for field in SOURCE_FIELD_ORDER
)
ADM_SOURCE_FIRST_ORDER_STATE_ORDER: Final[tuple[str, ...]] = tuple(
    f"{group}.{field}"
    for group in ("u", "p", "q")
    for field in SOURCE_FIELD_ORDER
)

VOCABULARY_RED1: Final[str] = "red1_independent"
VOCABULARY_MHG2: Final[str] = "mhg2_six_row"
SLOT_IDENTITY: Final[str] = "identity"
SLOT_NAME_IDENTITY: Final[str] = "documented_name_identity"
SLOT_EXACT_INVERSE: Final[str] = "exact_base_to_adm_two_jet_inverse"
SLOT_REFUSED_OWNER: Final[str] = "refused_owner_cannot_supply"
SLOT_NOT_GEOMETRY: Final[str] = "not_a_geometry_slot"

CONVERSION_RED1_MHG2: Final[str] = "red1_mhg2_equation_order_identity"
CONVERSION_FO1_ADM_GEOMETRY: Final[str] = "fo1_adm_26_input_geometry_slot"
CONVERSION_BASE_TO_ADM: Final[str] = "base_to_adm_two_jet_inverse"
CONVERSION_IMP1_TO_Q: Final[str] = "imp1_18_channel_to_q"
CONVERSION_COVERAGE_MATRIX: Final[str] = "def1_14_component_route_coverage"

INVENTORY_FO1_BASE_PHYSICAL: Final[str] = "fo1_base_physical_arguments"
INVENTORY_FO1_BASE_STATE: Final[str] = "fo1_base_first_order_state"
INVENTORY_FO1_BASE_TO_ADM: Final[str] = "fo1_base_to_adm_two_jet"
BASE_METRIC_TWO_JET_ORDER: Final[tuple[str, ...]] = (
    "h_tt",
    "h_tr",
    "h_rr",
    "areal_radius",
)
INVENTORY_ADM_SOURCE_PHYSICAL: Final[str] = "adm_source_physical_arguments"
INVENTORY_ADM_SOURCE_STATE: Final[str] = "adm_source_first_order_state"
INVENTORY_ADM_METRIC_TWO_JET: Final[str] = "adm_metric_two_jet"

OWNER_PATH_RED1: Final[str] = "src/recursive_horizons/fgc/spherical_reduction.py"
OWNER_PATH_MHG1: Final[str] = "src/recursive_horizons/fgc/modified_harmonic.py"
OWNER_PATH_MHG2: Final[str] = (
    "src/recursive_horizons/fgc/modified_harmonic_reference.py"
)
OWNER_PATH_MHG3_IMP1: Final[str] = (
    "src/recursive_horizons/fgc/modified_harmonic_implicit.py"
)
OWNER_PATH_FO1: Final[str] = (
    "src/recursive_horizons/fgc/modified_harmonic_first_order.py"
)
OWNER_PATH_SGB1_SOURCE: Final[str] = "src/recursive_horizons/fgc/sgb1_ctl1_source.py"
OWNER_PATH_CON4: Final[str] = "src/recursive_horizons/fgc/constraint_system.py"
OWNER_PATH_CON1: Final[str] = (
    "src/recursive_horizons/fgc/modified_harmonic_constraints.py"
)
OWNER_PATH_SRC1: Final[str] = (
    "src/recursive_horizons/fgc/evolution/nonlinear_source.py"
)
OWNER_PATH_GEOMETRY: Final[str] = "src/recursive_horizons/fgc/def1_geometry_error.py"
OWNER_PATH_PROVIDERS: Final[str] = "src/recursive_horizons/fgc/def1_stab1_providers.py"
OWNER_PATH_TDG11_IMP1: Final[str] = (
    "src/recursive_horizons/fgc/evolution/tdg11_imp1_ledger.py"
)
OWNER_PATH_REGULAR_CENTER: Final[str] = "src/recursive_horizons/fgc/regular_center.py"

EXACT_INVERSE_JUSTIFICATION: Final[str] = (
    "exact BASE-to-ADM two-jet inverse owned with Jet2.compose: "
    "lambda=sqrt(h_rr) at a supplied positive root with q'=1/(2 lambda) and "
    "q''=-1/(4 lambda^3); shift=h_tr/h_rr; alpha=sqrt(h_tr^2/h_rr-h_tt) with "
    "the same compose coefficients at the supplied positive lapse root. Roots "
    "must square exactly to the radicands. Roundtrip is "
    "state_from_generalized_adm_pg_fixture."
)
TANGENT_JUSTIFICATION: Final[str] = (
    "FO1 and ADM metric two-jets do not own the affine/null tangent (k.t, k.r)"
)
SECOND_JET_JUSTIFICATION: Final[str] = (
    "first-order state (u,p,q) does not include p_t, p_r, q_r; missing second-"
    "jet slots are refused rather than zeroed"
)
INCOMPLETE_TWO_JET_JUSTIFICATION: Final[str] = (
    "BASE-to-ADM two-jet inverse requires complete two-jets "
    "(value,dt,dr,dtt,dtr,drr); first-order (u,p,q) cannot fill ADM second-jet "
    "slots"
)
BASE_TO_ADM_INVERSE_IDENTITY: Final[dict[str, str]] = {
    "lambda": (
        "sqrt(h_rr) via Jet2.compose(root, 1/(2 root), -1/(4 root^3)) "
        "at a supplied positive radial-scale root that squares exactly to h_rr.value"
    ),
    "shift": "h_tr/h_rr by Jet2 division",
    "alpha": (
        "sqrt(h_tr^2/h_rr-h_tt) via Jet2.compose with the same q',q'' at a "
        "supplied positive lapse root that squares exactly to the radicand"
    ),
    "areal_radius": "identity onto geometry R",
    "phi_chi": "not geometry inputs",
    "roundtrip_owner": "state_from_generalized_adm_pg_fixture",
}
PHI_CHI_JUSTIFICATION: Final[str] = (
    "phi and chi are FO1/ADM/IMP1 fields and are not among the 26 geometry "
    "inputs; they are not filled as implicit zeros"
)

_QUALIFICATION_CONTEXT: Final[str] = (
    "synthetic qualification control; not COL1 and not a trajectory"
)
_QUALIFICATION_UNIT: Final[str] = "geometry-input units and complete-Q units"


class Def1Stab1QualificationError(Def1Stab1Error):
    """Fail-closed conversion-ownership, provenance, or coverage error."""


def _rational(name: str, value: object) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise TypeError(f"{name} must be an exact rational value")
    return value if isinstance(value, Fraction) else Q(int(value))


def _nonnegative_rational(name: str, value: object) -> Fraction:
    result = _rational(name, value)
    if result < 0:
        raise Def1Stab1QualificationError(f"{name} must be nonnegative")
    return result


def _boolean(name: str, value: object) -> bool:
    if not isinstance(value, bool):
        raise TypeError(f"{name} must be a bool")
    return value


def _text(name: str, value: object) -> str:
    if not isinstance(value, str) or not value.strip():
        raise Def1Stab1QualificationError(f"{name} must be a nonempty string")
    return value


def _jsonable(value: object) -> object:
    if isinstance(value, tuple):
        return [_jsonable(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    return value


def _canonical_digest(value: object) -> str:
    return sha256(canonical_json_bytes(_jsonable(value))).hexdigest()


def _owner_relative(relative: str) -> str:
    relative = _text("owner_path", relative)
    if relative.startswith("/") or "\\" in relative or ".." in Path(relative).parts:
        raise Def1Stab1QualificationError(
            f"owner_path {relative!r} must be a repo-relative POSIX file path"
        )
    return relative


def owner_file_sha256(relative: str) -> str:
    relative = _owner_relative(relative)
    try:
        raw = read_regular_file(_REPO_ROOT, relative)
    except UnsafePathError as error:
        raise Def1Stab1QualificationError(
            f"owner_path {relative!r} must be a unique regular non-symlink file"
        ) from error
    return sha256(raw).hexdigest()


def live_imp1_18_channel_order() -> tuple[str, ...]:
    if TDG6_COMPLETE_STATE_CHANNELS != IMP1_18_CHANNEL_ORDER:
        raise Def1Stab1QualificationError(
            "IMP1 18-channel order drifted from TDG6_COMPLETE_STATE_CHANNELS"
        )
    return IMP1_18_CHANNEL_ORDER


def _require_exact_order(
    value: object, expected: tuple[str, ...], *, label: str
) -> tuple[str, ...]:
    if isinstance(value, Mapping):
        raise Def1Stab1QualificationError(
            f"{label} refuses mapping aliases; supply the exact ordered sequence"
        )
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError(f"{label} must be a sequence of names")
    names = tuple(value)
    if names != expected:
        missing = [name for name in expected if name not in names]
        extra = [name for name in names if name not in expected]
        raise Def1Stab1QualificationError(
            f"{label} must be the exact owner order without omissions, aliases, "
            f"or reordering; missing={missing} extra={extra} got={names}"
        )
    return names


def _named_nonnegative_inventory(
    value: object, expected: tuple[str, ...], *, label: str
) -> tuple[tuple[str, Fraction], ...]:
    if isinstance(value, Mapping):
        missing = [name for name in expected if name not in value]
        extra = sorted(name for name in value if name not in expected)
        if missing or extra:
            raise Def1Stab1QualificationError(
                f"{label} must contain exactly {len(expected)} named slots; "
                f"missing={missing} extra={extra}"
            )
        return tuple(
            (name, _nonnegative_rational(f"{label}[{name}]", value[name]))
            for name in expected
        )
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError(f"{label} must be a mapping or a sequence")
    if len(value) != len(expected):
        raise Def1Stab1QualificationError(
            f"{label} must list exactly {len(expected)} entries in owner order"
        )
    if value and not isinstance(value[0], Sequence):
        return tuple(
            (
                expected[index],
                _nonnegative_rational(f"{label}[{expected[index]}]", item),
            )
            for index, item in enumerate(value)
        )
    names: list[str] = []
    ordered: list[tuple[str, Fraction]] = []
    for index, item in enumerate(value):
        if not isinstance(item, Sequence) or isinstance(item, (str, bytes)):
            raise TypeError(f"{label}[{index}] must be a (name, value) pair")
        if len(item) != 2:
            raise TypeError(f"{label}[{index}] must be a (name, value) pair")
        name, raw = item
        if not isinstance(name, str):
            raise TypeError(f"{label}[{index}] name must be a string")
        names.append(name)
        ordered.append((name, _nonnegative_rational(f"{label}[{name}]", raw)))
    if len(set(names)) != len(names):
        raise Def1Stab1QualificationError(
            f"{label} duplicate names are refused rather than collapsed"
        )
    if tuple(names) != expected:
        raise Def1Stab1QualificationError(
            f"{label} must follow owner order without aliases or reordering"
        )
    return tuple(ordered)


def conversion_identity_digest(conversion_name: str) -> str:
    name = _text("conversion_name", conversion_name)
    if name == CONVERSION_RED1_MHG2:
        payload = {
            "conversion": name,
            "pairs": RED1_TO_MHG2_EQUATION_PAIRS,
            "red1": INDEPENDENT_EQUATION_ORDER,
            "mhg2": MHG2_FULL_EQUATION_ORDER,
            "mhg1": MHG_EQUATION_ORDER,
            "fo1_physical": FO1_PHYSICAL_EQUATION_ORDER,
            "implicit": IMPLICIT_EQUATION_ORDER,
            "sgb1_source": SOURCE_EQUATION_ORDER,
        }
    elif name == CONVERSION_FO1_ADM_GEOMETRY:
        payload = {
            "conversion": name,
            "input_names": INPUT_NAMES,
            "jet_slots": FO1_GROUP_TO_JET_SLOT,
            "fo1_fields": FO1_FIELD_ORDER,
            "adm_fields": ADM_FIELD_ORDER,
            "source_fields": SOURCE_FIELD_ORDER,
            "fo1_physical": FO1_PHYSICAL_ARGUMENT_ORDER,
            "adm_physical": ADM_SOURCE_PHYSICAL_ARGUMENT_ORDER,
            "base_to_adm": BASE_TO_ADM_INVERSE_IDENTITY,
        }
    elif name == CONVERSION_BASE_TO_ADM:
        payload = {
            "conversion": name,
            "base_metric_order": BASE_METRIC_TWO_JET_ORDER,
            "adm_metric_fields": ("alpha", "shift", "lambda", "R"),
            "inverse": BASE_TO_ADM_INVERSE_IDENTITY,
            "roundtrip": "state_from_generalized_adm_pg_fixture",
        }
    elif name == CONVERSION_IMP1_TO_Q:
        payload = {
            "conversion": name,
            "channels": live_imp1_18_channel_order(),
            "geometry_metric_fields": METRIC_FIELDS,
        }
    elif name == CONVERSION_COVERAGE_MATRIX:
        payload = {
            "conversion": name,
            "components": ERROR_BUDGET_COMPONENTS,
            "routes": tuple(
                (route.component, route.route_id) for route in PROVIDER_ROUTES
            ),
        }
    else:
        raise Def1Stab1QualificationError(
            f"unknown conversion identity {conversion_name!r}"
        )
    return _canonical_digest(payload)


@dataclass(frozen=True, slots=True)
class QualificationProvenanceRecord:
    """Honest owner-byte provenance; does not authenticate future COL1 values."""

    artifact_id: str
    owner_path: str
    owner_sha256: str
    conversion_name: str
    conversion_identity: str
    claims_col1_values: bool
    authenticates_owner_bytes: bool
    authenticates_future_col1: bool

    def __post_init__(self) -> None:
        artifact = _text("artifact_id", self.artifact_id)
        path = _text("owner_path", self.owner_path)
        digest = _text("owner_sha256", self.owner_sha256)
        conversion = _text("conversion_name", self.conversion_name)
        identity = _text("conversion_identity", self.conversion_identity)
        claims = _boolean("claims_col1_values", self.claims_col1_values)
        owner_bytes = _boolean(
            "authenticates_owner_bytes", self.authenticates_owner_bytes
        )
        future = _boolean(
            "authenticates_future_col1", self.authenticates_future_col1
        )
        live = owner_file_sha256(path)
        if digest != live:
            raise Def1Stab1QualificationError(
                f"owner_sha256 for {path} does not match the live owner file"
            )
        expected_identity = conversion_identity_digest(conversion)
        if identity != expected_identity:
            raise Def1Stab1QualificationError(
                "conversion_identity does not match the live conversion digest"
            )
        if claims or future:
            raise Def1Stab1QualificationError(
                "provenance may bind owner artifact/path/hash/conversion "
                "identity but cannot claim future COL1 values"
            )
        if not owner_bytes:
            raise Def1Stab1QualificationError(
                "provenance that does not authenticate owner bytes is refused"
            )
        object.__setattr__(self, "artifact_id", artifact)
        object.__setattr__(self, "owner_path", path)
        object.__setattr__(self, "owner_sha256", digest)
        object.__setattr__(self, "conversion_name", conversion)
        object.__setattr__(self, "conversion_identity", identity)


def bind_qualification_provenance(
    *,
    artifact_id: str,
    owner_path: str,
    conversion_name: str,
) -> QualificationProvenanceRecord:
    return QualificationProvenanceRecord(
        artifact_id=_text("artifact_id", artifact_id),
        owner_path=_text("owner_path", owner_path),
        owner_sha256=owner_file_sha256(owner_path),
        conversion_name=_text("conversion_name", conversion_name),
        conversion_identity=conversion_identity_digest(conversion_name),
        claims_col1_values=False,
        authenticates_owner_bytes=True,
        authenticates_future_col1=False,
    )


@dataclass(frozen=True, slots=True)
class EquationOrderOwnerBinding:
    """One six-row residual owner, or a typed non-six-row refusal."""

    owner_id: str
    artifact_id: str
    owner_path: str
    vocabulary: str
    live_order: tuple[str, ...]
    maps_to_mhg2: bool
    refusal: str | None


def _six_row_binding(
    *,
    owner_id: str,
    artifact_id: str,
    owner_path: str,
    vocabulary: str,
    live_order: tuple[str, ...],
) -> EquationOrderOwnerBinding:
    expected = (
        INDEPENDENT_EQUATION_ORDER
        if vocabulary == VOCABULARY_RED1
        else MHG2_FULL_EQUATION_ORDER
    )
    if live_order != expected:
        raise Def1Stab1QualificationError(
            f"{owner_id} live equation order drifted from {vocabulary}"
        )
    if vocabulary == VOCABULARY_MHG2 and live_order != tuple(
        pair[1] for pair in RED1_TO_MHG2_EQUATION_PAIRS
    ):
        raise Def1Stab1QualificationError(
            f"{owner_id} MHG2 order is not the RED1 identity image"
        )
    if vocabulary == VOCABULARY_RED1 and live_order != tuple(
        pair[0] for pair in RED1_TO_MHG2_EQUATION_PAIRS
    ):
        raise Def1Stab1QualificationError(
            f"{owner_id} RED1 order is not the identity domain"
        )
    return EquationOrderOwnerBinding(
        owner_id=owner_id,
        artifact_id=artifact_id,
        owner_path=owner_path,
        vocabulary=vocabulary,
        live_order=live_order,
        maps_to_mhg2=True,
        refusal=None,
    )


def _refused_binding(
    *,
    owner_id: str,
    artifact_id: str,
    owner_path: str,
    live_order: tuple[str, ...],
    refusal: str,
) -> EquationOrderOwnerBinding:
    return EquationOrderOwnerBinding(
        owner_id=owner_id,
        artifact_id=artifact_id,
        owner_path=owner_path,
        vocabulary="not_mhg2_six_row",
        live_order=live_order,
        maps_to_mhg2=False,
        refusal=refusal,
    )


def _equation_order_owners() -> tuple[EquationOrderOwnerBinding, ...]:
    if FO1_FIELD_ORDER != BASE_FIELD_ORDER:
        raise Def1Stab1QualificationError("FO1 field order drifted from RED1 BASE")
    if SECOND_DERIVATIVE_ORDER != ("dtt", "dtr", "drr"):
        raise Def1Stab1QualificationError("RED1 second-jet order drifted")
    return (
        _six_row_binding(
            owner_id="RED1",
            artifact_id="FGC-1-HYP1-RED1",
            owner_path=OWNER_PATH_RED1,
            vocabulary=VOCABULARY_RED1,
            live_order=INDEPENDENT_EQUATION_ORDER,
        ),
        _six_row_binding(
            owner_id="MHG1",
            artifact_id="FGC-1-HYP1-MHG1",
            owner_path=OWNER_PATH_MHG1,
            vocabulary=VOCABULARY_MHG2,
            live_order=MHG_EQUATION_ORDER,
        ),
        _six_row_binding(
            owner_id="MHG2-REF1",
            artifact_id="FGC-1-HYP1-MHG2-REF1",
            owner_path=OWNER_PATH_MHG2,
            vocabulary=VOCABULARY_MHG2,
            live_order=MHG2_FULL_EQUATION_ORDER,
        ),
        _six_row_binding(
            owner_id="MHG3-IMP1",
            artifact_id="FGC-1-HYP1-MHG3-IMP1",
            owner_path=OWNER_PATH_MHG3_IMP1,
            vocabulary=VOCABULARY_MHG2,
            live_order=IMPLICIT_EQUATION_ORDER,
        ),
        _six_row_binding(
            owner_id="FO1-physical",
            artifact_id="FGC-1-HYP1-FO1-RC1",
            owner_path=OWNER_PATH_FO1,
            vocabulary=VOCABULARY_MHG2,
            live_order=FO1_PHYSICAL_EQUATION_ORDER,
        ),
        _six_row_binding(
            owner_id="SGB1-source",
            artifact_id="FGC-1-SGB1-CTL1",
            owner_path=OWNER_PATH_SGB1_SOURCE,
            vocabulary=VOCABULARY_MHG2,
            live_order=SOURCE_EQUATION_ORDER,
        ),
        _six_row_binding(
            owner_id="SRC1-NL1",
            artifact_id="FGC-1-SRC1-NL1",
            owner_path=OWNER_PATH_SRC1,
            vocabulary=VOCABULARY_MHG2,
            live_order=MHG2_FULL_EQUATION_ORDER,
        ),
        _six_row_binding(
            owner_id="CON4-unredefined-residual-names",
            artifact_id="FGC-1-CON4-PHY1",
            owner_path=OWNER_PATH_CON1,
            vocabulary=VOCABULARY_RED1,
            live_order=INDEPENDENT_EQUATION_ORDER,
        ),
        _refused_binding(
            owner_id="CON4-physical-projection",
            artifact_id="FGC-1-CON4-PHY1",
            owner_path=OWNER_PATH_CON1,
            live_order=PHYSICAL_PROJECTION_ORDER,
            refusal=(
                "Hamiltonian/momentum projections are not the MHG2 six-row residual"
            ),
        ),
        _refused_binding(
            owner_id="MHG1-gauge-constraint",
            artifact_id="FGC-1-HYP1-MHG1",
            owner_path=OWNER_PATH_MHG1,
            live_order=MHG_GAUGE_CONSTRAINT_ORDER,
            refusal="gauge constraints H^t, H^r are not the MHG2 six-row residual",
        ),
        _refused_binding(
            owner_id="CON4-gauge-vector",
            artifact_id="FGC-1-CON4-PHY1",
            owner_path=OWNER_PATH_CON4,
            live_order=ACTIVE_SPHERICAL_GAUGE_COMPONENTS,
            refusal="gauge vector C^t, C^r is not the MHG2 six-row residual",
        ),
        _refused_binding(
            owner_id="FO1-reduction-constraint",
            artifact_id="FGC-1-HYP1-FO1-RC1",
            owner_path=OWNER_PATH_FO1,
            live_order=FO1_REDUCTION_CONSTRAINT_ORDER,
            refusal=(
                "FO1 kinematic reduction constraints are not the MHG2 six-row residual"
            ),
        ),
        _refused_binding(
            owner_id="REG1-regular-center",
            artifact_id="FGC-1-HYP1-REG1",
            owner_path=OWNER_PATH_REGULAR_CENTER,
            live_order=REGULAR_EQUATION_ORDER,
            refusal=(
                "regular-centre equation names are not aliases of the MHG2 six-row order"
            ),
        ),
    )


def _derive_equation_order_identity() -> dict[str, object]:
    if FO1_GROUP_TO_JET_SLOT != tuple(
        zip(FO1_PHYSICAL_ARGUMENT_GROUPS, JET_COMPONENTS, strict=True)
    ):
        raise Def1Stab1QualificationError(
            "FO1 group-to-jet identity drifted from FO1 groups or geometry jets"
        )
    owners = _equation_order_owners()
    provenance = tuple(
        bind_qualification_provenance(
            artifact_id=owner.artifact_id,
            owner_path=owner.owner_path,
            conversion_name=CONVERSION_RED1_MHG2,
        )
        for owner in owners
        if owner.maps_to_mhg2
    )
    return {
        "pairs": RED1_TO_MHG2_EQUATION_PAIRS,
        "red1_order": INDEPENDENT_EQUATION_ORDER,
        "mhg2_order": MHG2_FULL_EQUATION_ORDER,
        "owners": owners,
        "provenance": provenance,
        "conversion_identity": conversion_identity_digest(CONVERSION_RED1_MHG2),
        "claims_col1_values": False,
        "global_pde_error_certified": False,
    }


@dataclass(frozen=True, slots=True)
class EquationOrderIdentityMap:
    """Exact RED1/MHG2 six-row identity plus live owner bindings."""

    pairs: tuple[tuple[str, str], ...]
    red1_order: tuple[str, ...]
    mhg2_order: tuple[str, ...]
    owners: tuple[EquationOrderOwnerBinding, ...]
    provenance: tuple[QualificationProvenanceRecord, ...]
    conversion_identity: str
    claims_col1_values: bool
    global_pde_error_certified: bool

    def __post_init__(self) -> None:
        if self.claims_col1_values or self.global_pde_error_certified:
            raise Def1Stab1QualificationError(
                "equation-order identity cannot claim COL1 or certify a global PDE error"
            )
        derived = _derive_equation_order_identity()
        for name, value in derived.items():
            if getattr(self, name) != value:
                raise Def1Stab1QualificationError(
                    f"equation-order identity {name} does not match the live owners"
                )


def red1_mhg2_equation_identity_map() -> EquationOrderIdentityMap:
    return EquationOrderIdentityMap(**_derive_equation_order_identity())  # type: ignore[arg-type]


def map_equation_order_to_mhg2(
    order: object, *, owner_id: str
) -> tuple[str, ...]:
    identity = red1_mhg2_equation_identity_map()
    owner = next((item for item in identity.owners if item.owner_id == owner_id), None)
    if owner is None:
        raise Def1Stab1QualificationError(f"unknown equation-order owner {owner_id!r}")
    if not owner.maps_to_mhg2:
        raise Def1Stab1QualificationError(
            f"{owner_id} cannot be mapped onto the MHG2 six-row order: {owner.refusal}"
        )
    _require_exact_order(order, owner.live_order, label=f"{owner_id}.equation_order")
    return identity.mhg2_order


@dataclass(frozen=True, slots=True)
class GeometrySlotAssignment:
    """One of the 26 geometry inputs classified against a source inventory."""

    geometry_slot: str
    source_inventory: str
    source_slot: str | None
    status: str
    justification: str


def _fo1_source_slot(field: str, jet: str, *, first_order_only: bool) -> str | None:
    group = next(group for group, slot in FO1_GROUP_TO_JET_SLOT if slot == jet)
    if first_order_only and group not in {"u", "p", "q"}:
        return None
    return f"{group}.{field}"


def _classify_geometry_slots(inventory: str) -> tuple[GeometrySlotAssignment, ...]:
    inventory = _text("source_inventory", inventory)
    first_order_only = inventory in {
        INVENTORY_FO1_BASE_STATE,
        INVENTORY_ADM_SOURCE_STATE,
    }
    fo1_base = inventory in {
        INVENTORY_FO1_BASE_PHYSICAL,
        INVENTORY_FO1_BASE_STATE,
    }
    adm_named = inventory in {
        INVENTORY_ADM_SOURCE_PHYSICAL,
        INVENTORY_ADM_SOURCE_STATE,
        INVENTORY_ADM_METRIC_TWO_JET,
    }
    if inventory not in {
        INVENTORY_FO1_BASE_PHYSICAL,
        INVENTORY_FO1_BASE_STATE,
        INVENTORY_ADM_SOURCE_PHYSICAL,
        INVENTORY_ADM_SOURCE_STATE,
        INVENTORY_ADM_METRIC_TWO_JET,
    }:
        raise Def1Stab1QualificationError(
            f"unknown FO1/ADM geometry inventory {inventory!r}"
        )
    assignments: list[GeometrySlotAssignment] = []
    for slot in INPUT_NAMES:
        if slot in TANGENT_SLOTS:
            assignments.append(
                GeometrySlotAssignment(
                    geometry_slot=slot,
                    source_inventory=inventory,
                    source_slot=None,
                    status=SLOT_REFUSED_OWNER,
                    justification=TANGENT_JUSTIFICATION,
                )
            )
            continue
        field, jet = slot.split(".", 1)
        if fo1_base and field in {"alpha", "shift", "lambda"}:
            if first_order_only:
                assignments.append(
                    GeometrySlotAssignment(
                        geometry_slot=slot,
                        source_inventory=inventory,
                        source_slot=None,
                        status=SLOT_REFUSED_OWNER,
                        justification=INCOMPLETE_TWO_JET_JUSTIFICATION,
                    )
                )
            else:
                assignments.append(
                    GeometrySlotAssignment(
                        geometry_slot=slot,
                        source_inventory=inventory,
                        source_slot=None,
                        status=SLOT_EXACT_INVERSE,
                        justification=EXACT_INVERSE_JUSTIFICATION,
                    )
                )
            continue
        if fo1_base and field == "R":
            source = _fo1_source_slot(
                "areal_radius", jet, first_order_only=first_order_only
            )
            if source is None:
                assignments.append(
                    GeometrySlotAssignment(
                        geometry_slot=slot,
                        source_inventory=inventory,
                        source_slot=None,
                        status=SLOT_REFUSED_OWNER,
                        justification=SECOND_JET_JUSTIFICATION,
                    )
                )
            else:
                assignments.append(
                    GeometrySlotAssignment(
                        geometry_slot=slot,
                        source_inventory=inventory,
                        source_slot=source,
                        status=SLOT_NAME_IDENTITY,
                        justification=(
                            "FO1/RED1 areal_radius is the same quantity as "
                            "geometry R; this is a documented name identity"
                        ),
                    )
                )
            continue
        if inventory == INVENTORY_ADM_METRIC_TWO_JET:
            assignments.append(
                GeometrySlotAssignment(
                    geometry_slot=slot,
                    source_inventory=inventory,
                    source_slot=slot,
                    status=SLOT_IDENTITY,
                    justification=(
                        "ADM metric two-jet slots identity-map onto geometry "
                        "INPUT_NAMES metric slots"
                    ),
                )
            )
            continue
        source_field = field
        source = _fo1_source_slot(
            source_field, jet, first_order_only=first_order_only
        )
        if source is None:
            assignments.append(
                GeometrySlotAssignment(
                    geometry_slot=slot,
                    source_inventory=inventory,
                    source_slot=None,
                    status=SLOT_REFUSED_OWNER,
                    justification=SECOND_JET_JUSTIFICATION,
                )
            )
            continue
        assignments.append(
            GeometrySlotAssignment(
                geometry_slot=slot,
                source_inventory=inventory,
                source_slot=source,
                status=SLOT_IDENTITY,
                justification=(
                    "FO1 group-to-jet identity on ADM/source field names "
                    f"({source} -> {slot})"
                ),
            )
        )
    if tuple(item.geometry_slot for item in assignments) != INPUT_NAMES:
        raise Def1Stab1QualificationError(
            "geometry-slot contract omitted or reordered INPUT_NAMES"
        )
    if adm_named and any(
        item.status == SLOT_EXACT_INVERSE for item in assignments
    ):
        raise Def1Stab1QualificationError(
            "ADM-named inventories must not use the BASE-to-ADM inverse slots"
        )
    return tuple(assignments)


def geometry_slot_zero_justification(slot: str) -> str | None:
    if slot not in INPUT_NAMES:
        raise Def1Stab1QualificationError(
            f"{slot!r} is not a 26-input geometry slot"
        )
    return None


@dataclass(frozen=True, slots=True)
class Fo1AdmGeometrySlotMap:
    """Complete 26-slot classification for every owned FO1/ADM inventory."""

    inventories: tuple[tuple[str, tuple[GeometrySlotAssignment, ...]], ...]
    field_identities: tuple[tuple[str, str, str | None, str], ...]
    jet_slot_identity: tuple[tuple[str, str], ...]
    provenance: tuple[QualificationProvenanceRecord, ...]
    conversion_identity: str
    claims_col1_values: bool
    global_pde_error_certified: bool

    def __post_init__(self) -> None:
        if self.claims_col1_values or self.global_pde_error_certified:
            raise Def1Stab1QualificationError(
                "FO1/ADM geometry map cannot claim COL1 or certify a global PDE error"
            )
        derived = _derive_fo1_adm_geometry_map()
        for name, value in derived.items():
            if getattr(self, name) != value:
                raise Def1Stab1QualificationError(
                    f"FO1/ADM geometry map {name} does not match live owners"
                )


def _derive_fo1_adm_geometry_map() -> dict[str, object]:
    if METRIC_FIELDS != ("alpha", "shift", "lambda", "R"):
        raise Def1Stab1QualificationError("geometry metric fields drifted")
    if ADM_FIELD_ORDER != (
        "alpha",
        "shift",
        "lambda",
        "areal_radius",
        "phi",
        "chi",
    ):
        raise Def1Stab1QualificationError("ADM_FIELD_ORDER drifted")
    if SOURCE_FIELD_ORDER != ("alpha", "shift", "lambda", "R", "phi", "chi"):
        raise Def1Stab1QualificationError("SOURCE_FIELD_ORDER drifted")
    inventories = tuple(
        (name, _classify_geometry_slots(name))
        for name in (
            INVENTORY_FO1_BASE_PHYSICAL,
            INVENTORY_FO1_BASE_STATE,
            INVENTORY_ADM_SOURCE_PHYSICAL,
            INVENTORY_ADM_SOURCE_STATE,
            INVENTORY_ADM_METRIC_TWO_JET,
        )
    )
    field_identities = (
        ("alpha", "alpha", "alpha", "alpha"),
        ("shift", "shift", "shift", "v"),
        ("lambda", "lambda", "lambda", "lambda"),
        ("areal_radius", "areal_radius", "R", "R"),
        ("phi", "phi", None, "phi"),
        ("chi", "chi", None, "chi"),
    )
    provenance = (
        bind_qualification_provenance(
            artifact_id="FGC-1-HYP1-FO1-RC1",
            owner_path=OWNER_PATH_FO1,
            conversion_name=CONVERSION_FO1_ADM_GEOMETRY,
        ),
        bind_qualification_provenance(
            artifact_id="FGC-1-HYP1-RED1",
            owner_path=OWNER_PATH_RED1,
            conversion_name=CONVERSION_FO1_ADM_GEOMETRY,
        ),
        bind_qualification_provenance(
            artifact_id="FGC-1-DEF1-STAB1",
            owner_path=OWNER_PATH_GEOMETRY,
            conversion_name=CONVERSION_FO1_ADM_GEOMETRY,
        ),
        bind_qualification_provenance(
            artifact_id="FGC-1-SGB1-CTL1",
            owner_path=OWNER_PATH_SGB1_SOURCE,
            conversion_name=CONVERSION_FO1_ADM_GEOMETRY,
        ),
    )
    return {
        "inventories": inventories,
        "field_identities": field_identities,
        "jet_slot_identity": FO1_GROUP_TO_JET_SLOT,
        "provenance": provenance,
        "conversion_identity": conversion_identity_digest(CONVERSION_FO1_ADM_GEOMETRY),
        "claims_col1_values": False,
        "global_pde_error_certified": False,
    }


def fo1_adm_geometry_slot_map() -> Fo1AdmGeometrySlotMap:
    return Fo1AdmGeometrySlotMap(**_derive_fo1_adm_geometry_map())  # type: ignore[arg-type]


def _inventory_expected_names(inventory: str) -> tuple[str, ...]:
    if inventory == INVENTORY_FO1_BASE_PHYSICAL:
        return FO1_PHYSICAL_ARGUMENT_ORDER
    if inventory == INVENTORY_FO1_BASE_STATE:
        return FO1_STATE_ORDER
    if inventory == INVENTORY_ADM_SOURCE_PHYSICAL:
        return ADM_SOURCE_PHYSICAL_ARGUMENT_ORDER
    if inventory == INVENTORY_ADM_SOURCE_STATE:
        return ADM_SOURCE_FIRST_ORDER_STATE_ORDER
    if inventory == INVENTORY_ADM_METRIC_TWO_JET:
        return ADM_METRIC_TWO_JET_ORDER
    raise Def1Stab1QualificationError(f"unknown FO1/ADM inventory {inventory!r}")


@dataclass(frozen=True, slots=True)
class Fo1GeometryConversion:
    """Partial or complete 26-input conversion; refused slots are not zeroed."""

    inventory: str
    assignments: tuple[GeometrySlotAssignment, ...]
    supplied: tuple[tuple[str, Fraction], ...]
    refused: tuple[str, ...]
    complete: bool
    claims_col1_values: bool
    global_pde_error_certified: bool
    used_implicit_zero: bool


def convert_fo1_to_geometry_slots(
    values: object,
    *,
    inventory: str,
    tangent: object | None = None,
    require_complete: bool = True,
) -> Fo1GeometryConversion:
    expected = _inventory_expected_names(inventory)
    pairs = _named_nonnegative_inventory(
        values, expected, label=f"{inventory}.values"
    )
    lookup = {name: value for name, value in pairs}
    assignments = _classify_geometry_slots(inventory)
    if tangent is not None:
        tangent_pairs = _named_nonnegative_inventory(
            tangent, TANGENT_SLOTS, label="tangent"
        )
        lookup.update(tangent_pairs)
        replaced: list[GeometrySlotAssignment] = []
        for item in assignments:
            if item.geometry_slot in TANGENT_SLOTS:
                replaced.append(
                    GeometrySlotAssignment(
                        geometry_slot=item.geometry_slot,
                        source_inventory=inventory,
                        source_slot=item.geometry_slot,
                        status=SLOT_IDENTITY,
                        justification=(
                            "caller-supplied affine/null tangent; FO1 still does "
                            "not own these slots"
                        ),
                    )
                )
            else:
                replaced.append(item)
        assignments = tuple(replaced)
    supplied: list[tuple[str, Fraction]] = []
    refused: list[str] = []
    for item in assignments:
        if item.status in {SLOT_IDENTITY, SLOT_NAME_IDENTITY} and item.source_slot:
            if item.source_slot not in lookup:
                refused.append(item.geometry_slot)
                continue
            supplied.append((item.geometry_slot, lookup[item.source_slot]))
        else:
            refused.append(item.geometry_slot)
    complete = tuple(name for name, _ in supplied) == INPUT_NAMES and not refused
    if require_complete and not complete:
        raise Def1Stab1QualificationError(
            f"{inventory} cannot supply complete 26-input geometry slots; "
            f"refused={refused}. Missing owners are not replaced by zero."
        )
    if any(geometry_slot_zero_justification(name) for name, _ in supplied):
        raise Def1Stab1QualificationError(
            "geometry conversion inserted an unjustified implicit zero"
        )
    return Fo1GeometryConversion(
        inventory=inventory,
        assignments=assignments,
        supplied=tuple(supplied),
        refused=tuple(refused),
        complete=complete,
        claims_col1_values=False,
        global_pde_error_certified=False,
        used_implicit_zero=False,
    )


def _complete_jet(name: str, value: object) -> Jet2:
    if isinstance(value, Jet2):
        return Jet2(
            value.value, value.dt, value.dr, value.dtt, value.dtr, value.drr
        )
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a Jet2 or a complete six-slot mapping")
    missing = [part for part in JET_COMPONENTS if part not in value]
    extra = sorted(part for part in value if part not in JET_COMPONENTS)
    if missing or extra:
        raise Def1Stab1QualificationError(
            f"{name} must be a complete two-jet; missing={missing} extra={extra}"
        )
    return Jet2(
        *(_rational(f"{name}.{part}", value[part]) for part in JET_COMPONENTS)
    )


def _ordered_complete_jets(
    value: object, expected: tuple[str, ...], *, label: str
) -> dict[str, Jet2]:
    if isinstance(value, Mapping):
        missing = [name for name in expected if name not in value]
        extra = sorted(name for name in value if name not in expected)
        if missing or extra:
            raise Def1Stab1QualificationError(
                f"{label} must contain exactly {list(expected)}; "
                f"missing={missing} extra={extra}"
            )
        return {
            name: _complete_jet(f"{label}.{name}", value[name]) for name in expected
        }
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError(
            f"{label} must be a mapping or a sequence of (name, jet) pairs"
        )
    if len(value) != len(expected):
        raise Def1Stab1QualificationError(
            f"{label} must list exactly {len(expected)} jets in owner order"
        )
    names: list[str] = []
    jets: dict[str, Jet2] = {}
    for index, item in enumerate(value):
        if not isinstance(item, Sequence) or isinstance(item, (str, bytes)):
            raise TypeError(f"{label}[{index}] must be a (name, jet) pair")
        if len(item) != 2:
            raise TypeError(f"{label}[{index}] must be a (name, jet) pair")
        name, raw = item
        if not isinstance(name, str):
            raise TypeError(f"{label}[{index}] name must be a string")
        names.append(name)
        jets[name] = _complete_jet(f"{label}.{name}", raw)
    if len(set(names)) != len(names):
        raise Def1Stab1QualificationError(
            f"{label} duplicate names are refused rather than collapsed"
        )
    if tuple(names) != expected:
        raise Def1Stab1QualificationError(
            f"{label} must follow owner order without aliases or reordering"
        )
    return {name: jets[name] for name in expected}


def _named_signed_inventory(
    value: object, expected: tuple[str, ...], *, label: str
) -> tuple[tuple[str, Fraction], ...]:
    if isinstance(value, Mapping):
        missing = [name for name in expected if name not in value]
        extra = sorted(name for name in value if name not in expected)
        if missing or extra:
            raise Def1Stab1QualificationError(
                f"{label} must contain exactly {len(expected)} named slots; "
                f"missing={missing} extra={extra}"
            )
        return tuple(
            (name, _rational(f"{label}[{name}]", value[name])) for name in expected
        )
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError(f"{label} must be a mapping or a sequence")
    if len(value) != len(expected):
        raise Def1Stab1QualificationError(
            f"{label} must list exactly {len(expected)} entries in owner order"
        )
    if value and not isinstance(value[0], Sequence):
        return tuple(
            (expected[index], _rational(f"{label}[{expected[index]}]", item))
            for index, item in enumerate(value)
        )
    names: list[str] = []
    ordered: list[tuple[str, Fraction]] = []
    for index, item in enumerate(value):
        if not isinstance(item, Sequence) or isinstance(item, (str, bytes)):
            raise TypeError(f"{label}[{index}] must be a (name, value) pair")
        if len(item) != 2:
            raise TypeError(f"{label}[{index}] must be a (name, value) pair")
        name, raw = item
        if not isinstance(name, str):
            raise TypeError(f"{label}[{index}] name must be a string")
        names.append(name)
        ordered.append((name, _rational(f"{label}[{name}]", raw)))
    if len(set(names)) != len(names):
        raise Def1Stab1QualificationError(
            f"{label} duplicate names are refused rather than collapsed"
        )
    if tuple(names) != expected:
        raise Def1Stab1QualificationError(
            f"{label} must follow owner order without aliases or reordering"
        )
    return tuple(ordered)


def _exact_positive_sqrt_jet(
    radicand: Jet2, root: object, *, label: str
) -> Jet2:
    root_value = _rational(f"{label}.root", root)
    if root_value <= 0:
        raise Def1Stab1QualificationError(
            f"{label} root must be the positive branch"
        )
    if root_value * root_value != radicand.value:
        raise Def1Stab1QualificationError(
            f"{label} root must square exactly to the radicand; "
            f"{root_value}^2 != {radicand.value}"
        )
    return radicand.compose(
        root_value,
        1 / (2 * root_value),
        -1 / (4 * root_value**3),
    )


def _derive_adm_from_base(
    h_tt: Jet2,
    h_tr: Jet2,
    h_rr: Jet2,
    areal_radius: Jet2,
    *,
    lapse_root: object,
    radial_scale_root: object,
) -> tuple[Jet2, Jet2, Jet2, Fraction, Fraction]:
    if areal_radius.value <= 0:
        raise Def1Stab1QualificationError("areal_radius.value must be positive")
    if h_rr.value <= 0:
        raise Def1Stab1QualificationError(
            "h_rr.value must be strictly positive for a positive radial scale"
        )
    determinant = h_tt.value * h_rr.value - h_tr.value * h_tr.value
    if determinant >= 0:
        raise Def1Stab1QualificationError(
            "BASE two-jet must be Lorentzian (h_tt h_rr - h_tr^2 < 0)"
        )
    try:
        lambda_jet = _exact_positive_sqrt_jet(
            h_rr, radial_scale_root, label="radial_scale"
        )
        shift = h_tr / h_rr
        lapse_radicand = (h_tr * h_tr) / h_rr - h_tt
        alpha = _exact_positive_sqrt_jet(lapse_radicand, lapse_root, label="lapse")
    except (ValueError, ZeroDivisionError) as error:
        raise Def1Stab1QualificationError(
            f"singular BASE-to-ADM inverse: {error}"
        ) from error
    if alpha.value <= 0 or lambda_jet.value <= 0:
        raise Def1Stab1QualificationError(
            "inverted ADM alpha and lambda must be positive"
        )
    return alpha, shift, lambda_jet, alpha.value, lambda_jet.value


def _jet_mapping(jet: Jet2) -> dict[str, Fraction]:
    return {part: getattr(jet, part) for part in JET_COMPONENTS}


def _roundtrip_base_metric(
    *,
    alpha: Jet2,
    shift: Jet2,
    lambda_jet: Jet2,
    areal_radius: Jet2,
    h_tt: Jet2,
    h_tr: Jet2,
    h_rr: Jet2,
) -> None:
    recovered = state_from_generalized_adm_pg_fixture(
        {
            "state": {
                "alpha": _jet_mapping(alpha),
                "shift": _jet_mapping(shift),
                "lambda": _jet_mapping(lambda_jet),
                "areal_radius": _jet_mapping(areal_radius),
                "phi": _jet_mapping(Jet2.constant(0)),
                "chi": _jet_mapping(Jet2.constant(0)),
            },
            "branch": "GR-0",
        }
    )
    if (
        recovered.h_tt != h_tt
        or recovered.h_tr != h_tr
        or recovered.h_rr != h_rr
        or recovered.areal_radius != areal_radius
    ):
        raise Def1Stab1QualificationError(
            "BASE-to-ADM inverse failed the exact roundtrip through "
            "state_from_generalized_adm_pg_fixture"
        )


def _fo1_field_jet(lookup: Mapping[str, Fraction], field: str) -> Jet2:
    return Jet2(
        *(lookup[f"{group}.{field}"] for group, _slot in FO1_GROUP_TO_JET_SLOT)
    )


@dataclass(frozen=True, slots=True)
class BaseToAdmTwoJetInverse:
    """Immutable exact BASE-to-ADM two-jet inverse with a checked roundtrip."""

    h_tt: Jet2
    h_tr: Jet2
    h_rr: Jet2
    areal_radius: Jet2
    alpha: Jet2
    shift: Jet2
    lambda_jet: Jet2
    lapse_root: Fraction
    radial_scale_root: Fraction
    roundtrip_holds: bool
    provenance: QualificationProvenanceRecord
    conversion_identity: str
    claims_col1_values: bool
    global_pde_error_certified: bool
    used_implicit_zero: bool

    def __post_init__(self) -> None:
        if self.claims_col1_values or self.global_pde_error_certified:
            raise Def1Stab1QualificationError(
                "BASE-to-ADM inverse cannot claim COL1 or certify a global PDE error"
            )
        if self.used_implicit_zero or not self.roundtrip_holds:
            raise Def1Stab1QualificationError(
                "BASE-to-ADM inverse cannot hide zeros or drop the roundtrip"
            )
        derived = _derive_adm_from_base(
            self.h_tt,
            self.h_tr,
            self.h_rr,
            self.areal_radius,
            lapse_root=self.lapse_root,
            radial_scale_root=self.radial_scale_root,
        )
        alpha, shift, lambda_jet, lapse_root, radial_root = derived
        if (
            self.alpha != alpha
            or self.shift != shift
            or self.lambda_jet != lambda_jet
            or self.lapse_root != lapse_root
            or self.radial_scale_root != radial_root
        ):
            raise Def1Stab1QualificationError(
                "BASE-to-ADM inverse fields do not match the live Jet2 inverse"
            )
        _roundtrip_base_metric(
            alpha=self.alpha,
            shift=self.shift,
            lambda_jet=self.lambda_jet,
            areal_radius=self.areal_radius,
            h_tt=self.h_tt,
            h_tr=self.h_tr,
            h_rr=self.h_rr,
        )
        if self.conversion_identity != conversion_identity_digest(CONVERSION_BASE_TO_ADM):
            raise Def1Stab1QualificationError(
                "BASE-to-ADM conversion_identity does not match the live digest"
            )
        if self.provenance.conversion_name != CONVERSION_BASE_TO_ADM:
            raise Def1Stab1QualificationError(
                "BASE-to-ADM provenance conversion_name mismatch"
            )


def invert_base_metric_two_jets(
    *,
    h_tt: object,
    h_tr: object,
    h_rr: object,
    areal_radius: object,
    lapse_root: object,
    radial_scale_root: object,
) -> BaseToAdmTwoJetInverse:
    base_tt = _complete_jet("h_tt", h_tt)
    base_tr = _complete_jet("h_tr", h_tr)
    base_rr = _complete_jet("h_rr", h_rr)
    radius = _complete_jet("areal_radius", areal_radius)
    alpha, shift, lambda_jet, lapse, radial = _derive_adm_from_base(
        base_tt,
        base_tr,
        base_rr,
        radius,
        lapse_root=lapse_root,
        radial_scale_root=radial_scale_root,
    )
    _roundtrip_base_metric(
        alpha=alpha,
        shift=shift,
        lambda_jet=lambda_jet,
        areal_radius=radius,
        h_tt=base_tt,
        h_tr=base_tr,
        h_rr=base_rr,
    )
    return BaseToAdmTwoJetInverse(
        h_tt=base_tt,
        h_tr=base_tr,
        h_rr=base_rr,
        areal_radius=radius,
        alpha=alpha,
        shift=shift,
        lambda_jet=lambda_jet,
        lapse_root=lapse,
        radial_scale_root=radial,
        roundtrip_holds=True,
        provenance=bind_qualification_provenance(
            artifact_id="FGC-1-HYP1-RED1",
            owner_path=OWNER_PATH_RED1,
            conversion_name=CONVERSION_BASE_TO_ADM,
        ),
        conversion_identity=conversion_identity_digest(CONVERSION_BASE_TO_ADM),
        claims_col1_values=False,
        global_pde_error_certified=False,
        used_implicit_zero=False,
    )


def invert_ordered_base_metric_two_jets(
    jets: object,
    *,
    lapse_root: object,
    radial_scale_root: object,
) -> BaseToAdmTwoJetInverse:
    ordered = _ordered_complete_jets(
        jets, BASE_METRIC_TWO_JET_ORDER, label="base.metric_two_jets"
    )
    return invert_base_metric_two_jets(
        h_tt=ordered["h_tt"],
        h_tr=ordered["h_tr"],
        h_rr=ordered["h_rr"],
        areal_radius=ordered["areal_radius"],
        lapse_root=lapse_root,
        radial_scale_root=radial_scale_root,
    )


def invert_fo1_physical_arguments_to_adm(
    values: object,
    *,
    lapse_root: object,
    radial_scale_root: object,
) -> BaseToAdmTwoJetInverse:
    pairs = _named_signed_inventory(
        values, FO1_PHYSICAL_ARGUMENT_ORDER, label="fo1.physical_arguments"
    )
    lookup = {name: value for name, value in pairs}
    return invert_base_metric_two_jets(
        h_tt=_fo1_field_jet(lookup, "h_tt"),
        h_tr=_fo1_field_jet(lookup, "h_tr"),
        h_rr=_fo1_field_jet(lookup, "h_rr"),
        areal_radius=_fo1_field_jet(lookup, "areal_radius"),
        lapse_root=lapse_root,
        radial_scale_root=radial_scale_root,
    )


def convert_base_to_adm_geometry_slots(
    inverse: BaseToAdmTwoJetInverse,
    *,
    tangent: object,
) -> Fo1GeometryConversion:
    if not isinstance(inverse, BaseToAdmTwoJetInverse):
        raise TypeError("inverse must be a BaseToAdmTwoJetInverse")
    tangent_pairs = _named_signed_inventory(
        tangent, TANGENT_SLOTS, label="tangent"
    )
    supplied: list[tuple[str, Fraction]] = []
    for field, jet in (
        ("alpha", inverse.alpha),
        ("shift", inverse.shift),
        ("lambda", inverse.lambda_jet),
        ("R", inverse.areal_radius),
    ):
        for part in JET_COMPONENTS:
            supplied.append((f"{field}.{part}", getattr(jet, part)))
    supplied.extend(tangent_pairs)
    if tuple(name for name, _ in supplied) != INPUT_NAMES:
        raise Def1Stab1QualificationError(
            "BASE-to-ADM geometry pack omitted or reordered INPUT_NAMES"
        )
    assignments: list[GeometrySlotAssignment] = []
    for slot in INPUT_NAMES:
        if slot in TANGENT_SLOTS:
            assignments.append(
                GeometrySlotAssignment(
                    geometry_slot=slot,
                    source_inventory=INVENTORY_FO1_BASE_TO_ADM,
                    source_slot=slot,
                    status=SLOT_IDENTITY,
                    justification=(
                        "caller-supplied affine/null tangent; FO1 still does "
                        "not own these slots"
                    ),
                )
            )
        elif slot.startswith("R."):
            assignments.append(
                GeometrySlotAssignment(
                    geometry_slot=slot,
                    source_inventory=INVENTORY_FO1_BASE_TO_ADM,
                    source_slot=f"areal_radius.{slot.split('.', 1)[1]}",
                    status=SLOT_NAME_IDENTITY,
                    justification=(
                        "FO1/RED1 areal_radius is the same quantity as geometry R"
                    ),
                )
            )
        else:
            assignments.append(
                GeometrySlotAssignment(
                    geometry_slot=slot,
                    source_inventory=INVENTORY_FO1_BASE_TO_ADM,
                    source_slot=None,
                    status=SLOT_EXACT_INVERSE,
                    justification=EXACT_INVERSE_JUSTIFICATION,
                )
            )
    return Fo1GeometryConversion(
        inventory=INVENTORY_FO1_BASE_TO_ADM,
        assignments=tuple(assignments),
        supplied=tuple(supplied),
        refused=(),
        complete=True,
        claims_col1_values=False,
        global_pde_error_certified=False,
        used_implicit_zero=False,
    )


def phi_chi_geometry_status() -> str:
    return SLOT_NOT_GEOMETRY + ": " + PHI_CHI_JUSTIFICATION


@dataclass(frozen=True, slots=True)
class Imp1ChannelGeometryAssignment:
    channel: str
    geometry_slot: str | None
    status: str
    justification: str


def imp1_channel_geometry_assignments() -> tuple[Imp1ChannelGeometryAssignment, ...]:
    live_imp1_18_channel_order()
    group_to_jet = dict(FO1_GROUP_TO_JET_SLOT)
    assignments: list[Imp1ChannelGeometryAssignment] = []
    for channel in IMP1_18_CHANNEL_ORDER:
        block, field = channel.split(":", 1)
        jet = group_to_jet[block]
        if field in {"phi", "chi"}:
            assignments.append(
                Imp1ChannelGeometryAssignment(
                    channel=channel,
                    geometry_slot=None,
                    status=SLOT_NOT_GEOMETRY,
                    justification=PHI_CHI_JUSTIFICATION,
                )
            )
            continue
        geometry_field = "shift" if field == "v" else field
        slot = f"{geometry_field}.{jet}"
        if slot not in INPUT_NAMES:
            raise Def1Stab1QualificationError(
                f"IMP1 channel {channel} produced a non-geometry slot"
            )
        status = SLOT_NAME_IDENTITY if field == "v" else SLOT_IDENTITY
        justification = (
            "IMP1 v is the same quantity as ADM/geometry shift; this is a "
            "documented name identity, not a silent alias in either inventory"
            if field == "v"
            else f"IMP1 {channel} identity-maps onto {slot}"
        )
        assignments.append(
            Imp1ChannelGeometryAssignment(
                channel=channel,
                geometry_slot=slot,
                status=status,
                justification=justification,
            )
        )
    if tuple(item.channel for item in assignments) != IMP1_18_CHANNEL_ORDER:
        raise Def1Stab1QualificationError("IMP1 channel geometry map reordered")
    return tuple(assignments)


@dataclass(frozen=True, slots=True)
class Imp1ToQConversion:
    """Conditional IMP1 18-channel debit; never a global PDE certificate."""

    channel_order: tuple[str, ...]
    channel_debits: tuple[tuple[str, Fraction], ...]
    lipschitz: tuple[tuple[str, Fraction], ...]
    additive_q: Fraction
    status: str
    source: str
    context: str
    unit: str
    provenance: QualificationProvenanceRecord
    conversion_identity: str
    used_measured_q: bool
    global_pde_error_certified: bool
    richardson_treated_as_global_pde_error: bool
    imp1_admission_debit_treated_as_global_pde_error: bool
    def1_error_map_passed: bool
    claims_col1_values: bool

    def __post_init__(self) -> None:
        channels = live_imp1_18_channel_order()
        if self.channel_order != channels:
            raise Def1Stab1QualificationError(
                "IMP1-to-Q channel_order drifted from TDG6_COMPLETE_STATE_CHANNELS"
            )
        if tuple(name for name, _ in self.channel_debits) != channels:
            raise Def1Stab1QualificationError(
                "IMP1-to-Q debits are not in the exact 18-channel owner order"
            )
        if tuple(name for name, _ in self.lipschitz) != channels:
            raise Def1Stab1QualificationError(
                "IMP1-to-Q Lipschitz evidence is not in the exact 18-channel order"
            )
        total = sum(
            (
                debit * lipschitz
                for (_, debit), (_, lipschitz) in zip(
                    self.channel_debits, self.lipschitz, strict=True
                )
            ),
            Q(0),
        )
        if self.additive_q != total:
            raise Def1Stab1QualificationError(
                "IMP1-to-Q additive_q must equal the exact Lipschitz sum"
            )
        if self.status != PREMISE_STATUS_CONDITIONAL:
            raise Def1Stab1QualificationError(
                "IMP1-to-Q conversion remains a conditional premise"
            )
        if self.source != SOURCE_IMP1_ADMISSION_DEBIT:
            raise Def1Stab1QualificationError(
                "IMP1-to-Q source must remain imp1_admission_debit"
            )
        if (
            self.used_measured_q
            or self.global_pde_error_certified
            or self.richardson_treated_as_global_pde_error
            or self.imp1_admission_debit_treated_as_global_pde_error
            or self.def1_error_map_passed
            or self.claims_col1_values
        ):
            raise Def1Stab1QualificationError(
                "IMP1-to-Q conversion cannot certify a global PDE error, set "
                "DEF1_error_map_passed, use measured Q, or claim COL1 values"
            )
        if self.conversion_identity != conversion_identity_digest(CONVERSION_IMP1_TO_Q):
            raise Def1Stab1QualificationError(
                "IMP1-to-Q conversion_identity does not match the live digest"
            )
        if self.provenance.conversion_name != CONVERSION_IMP1_TO_Q:
            raise Def1Stab1QualificationError(
                "IMP1-to-Q provenance conversion_name mismatch"
            )
        _text("context", self.context)
        _text("unit", self.unit)


def convert_imp1_channels_to_q(
    *,
    channel_debits: object,
    lipschitz: object,
    context: object,
    unit: object,
    provenance: QualificationProvenanceRecord | None = None,
) -> Imp1ToQConversion:
    if isinstance(lipschitz, (Fraction, Integral)) and not isinstance(lipschitz, bool):
        raise Def1Stab1QualificationError(
            "IMP1-to-Q refuses a universal numeric Lipschitz factor; supply "
            "nonnegative per-channel Lipschitz evidence in the 18-channel order"
        )
    channels = live_imp1_18_channel_order()
    debit_pairs = _named_nonnegative_inventory(
        channel_debits, channels, label="imp1.channel_debits"
    )
    lipschitz_pairs = _named_nonnegative_inventory(
        lipschitz, channels, label="imp1.lipschitz"
    )
    total = sum(
        (
            debit * bound
            for (_, debit), (_, bound) in zip(debit_pairs, lipschitz_pairs, strict=True)
        ),
        Q(0),
    )
    record = provenance or bind_qualification_provenance(
        artifact_id="FGC-1-TDG11-IMP1",
        owner_path=OWNER_PATH_TDG11_IMP1,
        conversion_name=CONVERSION_IMP1_TO_Q,
    )
    return Imp1ToQConversion(
        channel_order=channels,
        channel_debits=debit_pairs,
        lipschitz=lipschitz_pairs,
        additive_q=total,
        status=PREMISE_STATUS_CONDITIONAL,
        source=SOURCE_IMP1_ADMISSION_DEBIT,
        context=_text("context", context),
        unit=_text("unit", unit),
        provenance=record,
        conversion_identity=conversion_identity_digest(CONVERSION_IMP1_TO_Q),
        used_measured_q=False,
        global_pde_error_certified=False,
        richardson_treated_as_global_pde_error=False,
        imp1_admission_debit_treated_as_global_pde_error=False,
        def1_error_map_passed=False,
        claims_col1_values=False,
    )


@dataclass(frozen=True, slots=True)
class ProviderRouteSpec:
    component: str
    route_id: str
    conversion: str


PROVIDER_ROUTES: Final[tuple[ProviderRouteSpec, ...]] = (
    ProviderRouteSpec(
        "spatial_temporal", "declared_exact_algebra", "declared radii or debit"
    ),
    ProviderRouteSpec(
        "spatial_temporal", "richardson_conditional", "conditional Richardson"
    ),
    ProviderRouteSpec(
        "spatial_temporal", "imp1_conditional", "conditional IMP1 admission"
    ),
    ProviderRouteSpec(
        "physical_constraint", "complete_residual", "complete 4x4 residual"
    ),
    ProviderRouteSpec(
        "gauge_constraint", "extension_residual", "complete gauge extension"
    ),
    ProviderRouteSpec(
        "gauge_constraint", "c_to_jet_inverse", "certified C-to-jet inverse"
    ),
    ProviderRouteSpec(
        "reduction_constraint",
        "discrepancy_operator",
        "first-derivative discrepancies",
    ),
    ProviderRouteSpec("initial_data", "complete_radii", "complete 26-input radii"),
    ProviderRouteSpec(
        "nonlinear_source", "inverse_j_residual", "inverse-J times six-residual"
    ),
    ProviderRouteSpec("affine", "gronwall_transport", "Gronwall transport"),
    ProviderRouteSpec("nullness", "residual_tangent_bound", "null-to-tangent bound"),
    ProviderRouteSpec(
        "trajectory_alignment", "complete_radii", "complete 26-input radii"
    ),
    ProviderRouteSpec(
        "interpolation", "lipschitz_remainder", "first-derivative Lipschitz"
    ),
    ProviderRouteSpec(
        "interpolation", "second_order_remainder", "second-order remainder"
    ),
    ProviderRouteSpec("extraction", "unsampled_gap", "affine-gap remainder"),
    ProviderRouteSpec("boundary", "additional_debit", "supplied boundary debit"),
    ProviderRouteSpec("boundary", "zero_guarded", "guarded explicit zero"),
    ProviderRouteSpec("conservation", "additional_debit", "supplied conservation debit"),
    ProviderRouteSpec("conservation", "zero_guarded", "guarded explicit zero"),
    ProviderRouteSpec("arithmetic", "explicit_debit", "explicit arithmetic debit"),
    ProviderRouteSpec("arithmetic", "imp1_conditional", "conditional IMP1 admission"),
)


@dataclass(frozen=True, slots=True)
class ProviderRouteControl:
    component: str
    route_id: str
    conversion: str
    positive_control_passed: bool
    injected_failure_refused: bool
    trajectory_values_evaluated: bool
    def1_booleans_evaluated: bool
    def1_boolean_names: tuple[str, ...]


def _decl(component: str) -> dict[str, str]:
    return {
        "context": _QUALIFICATION_CONTEXT,
        "unit": _QUALIFICATION_UNIT,
        "provenance": f"{component}: synthetic qualification control",
    }


def _zero_tensor() -> tuple[tuple[Fraction, ...], ...]:
    zero = (Q(0), Q(0), Q(0), Q(0))
    return (zero, zero, zero, zero)


def _injected_tensor() -> tuple[tuple[Fraction, ...], ...]:
    rows = [list(row) for row in _zero_tensor()]
    rows[0][0] = Q(1)
    return tuple(tuple(row) for row in rows)


def _mass_flux(*, valid: bool = True):
    inverse = inverse_base_metric_from_adm(1, 0, 1)
    if valid:
        return assess_mass_flux_ledger(
            radius=4,
            inverse_metric=inverse,
            radius_derivatives=(Q(0), Q(1)),
            coupling_F=1,
            mixed_equation_residual=((0, 0), (0, 0)),
            delta_mass=0,
            integrated_flux=0,
            residual_enclosure=0,
        )
    return assess_mass_flux_ledger(
        radius=4,
        inverse_metric=inverse,
        radius_derivatives=(Q(0), Q(1)),
        coupling_F=1,
        mixed_equation_residual=((0, 0), (0, 0)),
        delta_mass=Q(1, 2),
        integrated_flux=0,
        residual_enclosure=Q(1, 4),
    )


def _expect_provider_error(operation, fragment: str) -> None:
    try:
        operation()
    except Def1Stab1ProviderError as error:
        if fragment not in str(error):
            raise Def1Stab1QualificationError(
                f"injected failure message {str(error)!r} missed {fragment!r}"
            ) from error
        return
    raise Def1Stab1QualificationError(
        f"injected failure was not refused (expected {fragment!r})"
    )


def _run_provider_route(route: ProviderRouteSpec) -> ProviderRouteControl:
    residual = _injected_tensor()
    flux = _mass_flux()
    if route.route_id == "declared_exact_algebra":
        record = spatial_temporal_provider(
            source=SOURCE_EXACT_ALGEBRA,
            status=PREMISE_STATUS_PROVEN,
            radii=complete_zero_radii(),
            additive_q=0,
            **_decl(route.component),
        )
        _expect_provider_error(
            lambda: spatial_temporal_provider(
                source=SOURCE_EXACT_ALGEBRA,
                status=PREMISE_STATUS_PROVEN,
                radii=complete_zero_radii(),
                hamiltonian=0,
                **_decl(route.component),
            ),
            "unexpected arguments",
        )
    elif route.route_id == "richardson_conditional" and route.component == "spatial_temporal":
        record = spatial_temporal_provider(
            source=SOURCE_RICHARDSON,
            status=PREMISE_STATUS_CONDITIONAL,
            additive_q=Q(1, 16),
            **_decl(route.component),
        )
        _expect_provider_error(
            lambda: spatial_temporal_provider(
                source=SOURCE_RICHARDSON,
                status=PREMISE_STATUS_PROVEN,
                additive_q=Q(1, 16),
                **_decl(route.component),
            ),
            "conditional",
        )
    elif route.route_id == "imp1_conditional" and route.component == "spatial_temporal":
        record = spatial_temporal_provider(
            source=SOURCE_IMP1_ADMISSION_DEBIT,
            status=PREMISE_STATUS_CONDITIONAL,
            additive_q=Q(1, 8),
            **_decl(route.component),
        )
        _expect_provider_error(
            lambda: spatial_temporal_provider(
                source=SOURCE_IMP1_ADMISSION_DEBIT,
                status=PREMISE_STATUS_CONDITIONAL,
                additive_q=0,
                **_decl(route.component),
            ),
            "zero-only",
        )
    elif route.route_id == "complete_residual":
        record = physical_constraint_provider(
            residual_E=residual,
            tangent=(Q(1), Q(1)),
            coupling_F=2,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_EXACT_ALGEBRA,
            **_decl(route.component),
        )
        _expect_provider_error(
            lambda: physical_constraint_provider(
                residual_E=residual,
                tangent=(Q(1), Q(1)),
                coupling_F=2,
                hamiltonian=0,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_EXACT_ALGEBRA,
                **_decl(route.component),
            ),
            "unexpected arguments",
        )
    elif route.route_id == "extension_residual":
        record = gauge_constraint_provider(
            extension_E=residual,
            tangent=(Q(1), Q(1)),
            coupling_F=2,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_EXACT_ALGEBRA,
            **_decl(route.component),
        )
        _expect_provider_error(
            lambda: gauge_constraint_provider(
                constraint_C=(Q(1), Q(0)),
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_EXACT_ALGEBRA,
                **_decl(route.component),
            ),
            "C-only",
        )
    elif route.route_id == "c_to_jet_inverse":
        record = gauge_constraint_provider(
            constraint_C=(Q(1, 8), Q(0)),
            c_to_jet_inverse_bound=2,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl(route.component),
        )
        _expect_provider_error(
            lambda: gauge_constraint_provider(
                constraint_C=(Q(1), Q(0)),
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                **_decl(route.component),
            ),
            "C-only",
        )
    elif route.route_id == "discrepancy_operator":
        discrepancies = {name: Q(0) for name in METRIC_FIRST_DERIVATIVE_INPUTS}
        discrepancies["alpha.dt"] = Q(1, 5)
        record = reduction_constraint_provider(
            first_derivative_discrepancies=discrepancies,
            derivative_operator_bound=3,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_EXACT_ALGEBRA,
            **_decl(route.component),
        )
        _expect_provider_error(
            lambda: reduction_constraint_provider(
                first_derivative_discrepancies={"alpha.dt": Q(1, 5)},
                derivative_operator_bound=3,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_EXACT_ALGEBRA,
                **_decl(route.component),
            ),
            "missing",
        )
    elif route.component == "initial_data":
        radii = complete_zero_radii()
        radii["R.dtt"] = Q(1, 32)
        record = initial_data_provider(
            radii=radii,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl(route.component),
        )
        _expect_provider_error(
            lambda: initial_data_provider(
                radii={"R.dtt": Q(1, 32)},
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                **_decl(route.component),
            ),
            "missing",
        )
    elif route.route_id == "inverse_j_residual":
        record = nonlinear_source_provider(
            full_residual=(Q(1, 8), 0, 0, 0, 0, 0),
            inverse_jacobian_bound=2,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_EXACT_ALGEBRA,
            **_decl(route.component),
        )
        _expect_provider_error(
            lambda: nonlinear_source_provider(
                full_residual=(Q(1, 8), 0, 0),
                inverse_jacobian_bound=2,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_EXACT_ALGEBRA,
                **_decl(route.component),
            ),
            str(FULL_RESIDUAL_LENGTH),
        )
    elif route.route_id == "gronwall_transport":
        record = affine_provider(
            transport_defect=Q(1, 10),
            affine_interval=Q(1, 8),
            gronwall_factor=3,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl(route.component),
        )
        _expect_provider_error(
            lambda: affine_residual_times_step_refused(
                transport_defect=Q(1, 10), affine_interval=Q(1, 8)
            ),
            "residual times step",
        )
    elif route.route_id == "residual_tangent_bound":
        record = nullness_provider(
            null_residual=Q(1, 6),
            null_to_tangent_bound=2,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl(route.component),
        )
        _expect_provider_error(
            lambda: nullness_provider(
                null_residual=Q(1, 6),
                null_to_tangent_bound=0,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                **_decl(route.component),
            ),
            "cannot default to zero",
        )
    elif route.component == "trajectory_alignment":
        radii = complete_zero_radii()
        radii["k.r"] = Q(1, 20)
        record = trajectory_alignment_provider(
            radii=radii,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl(route.component),
        )
        _expect_provider_error(
            lambda: trajectory_alignment_provider(
                radii={"k.r": Q(1, 20)},
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                **_decl(route.component),
            ),
            "missing",
        )
    elif route.route_id == "lipschitz_remainder":
        record = interpolation_provider(
            sample_spacing=Q(1, 8),
            first_derivative_enclosures=complete_zero_radii() | {"R.value": Q(2)},
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl(route.component),
        )
        _expect_provider_error(
            lambda: interpolation_provider(
                sample_spacing=Q(1, 8),
                first_derivative_enclosures={},
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                **_decl(route.component),
            ),
            "derivative bounds are missing",
        )
    elif route.route_id == "second_order_remainder":
        record = interpolation_provider(
            sample_spacing=Q(1, 4),
            first_derivative_enclosures=complete_zero_radii(),
            second_derivative_enclosures={"R.value": Q(2)},
            remainder_radii={
                name: Q(0) for name in INPUT_NAMES if name != "R.value"
            },
            claim_second_order_remainder=True,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl(route.component),
        )
        _expect_provider_error(
            lambda: interpolation_provider(
                sample_spacing=Q(1, 8),
                first_derivative_enclosures=complete_zero_radii(),
                claim_second_order_remainder=True,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                **_decl(route.component),
            ),
            "second-order",
        )
    elif route.route_id == "unsampled_gap":
        record = extraction_provider(
            affine_samples=(0, Q(1, 2), 1),
            derivative_bound=Q(1, 3),
            sample_q_values=(1, 1, 1),
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl(route.component),
        )
        _expect_provider_error(
            lambda: extraction_provider(
                affine_samples=(0, Q(1, 2), 1),
                derivative_bound=None,
                sample_q_values=(1, 1, 1),
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                **_decl(route.component),
            ),
            "continuous-Q",
        )
    elif route.route_id == "additional_debit" and route.component == "boundary":
        record = boundary_provider(
            additional_debit=Q(1, 7),
            physical_causality_passed=True,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl(route.component),
        )
        _expect_provider_error(
            lambda: boundary_provider(
                additional_debit=0,
                physical_causality_passed=True,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                **_decl(route.component),
            ),
            "physical causality alone",
        )
    elif route.route_id == "zero_guarded" and route.component == "boundary":
        record = boundary_provider(
            additional_debit=0,
            no_influence_premise=True,
            coverage_premise=True,
            independent_guard=True,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_EXACT_ALGEBRA,
            **_decl(route.component),
        )
        _expect_provider_error(
            lambda: boundary_provider(
                additional_debit=0,
                physical_causality_passed=True,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_EXACT_ALGEBRA,
                **_decl(route.component),
            ),
            "physical causality alone",
        )
    elif route.route_id == "additional_debit" and route.component == "conservation":
        record = conservation_provider(
            additional_debit=Q(1, 9),
            mass_flux=flux,
            coverage_premise=True,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl(route.component),
        )
        _expect_provider_error(
            lambda: conservation_provider(
                additional_debit=Q(1, 9),
                mass_flux=_mass_flux(valid=False),
                coverage_premise=True,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                **_decl(route.component),
            ),
            "mass-flux veto",
        )
    elif route.route_id == "zero_guarded" and route.component == "conservation":
        record = conservation_provider(
            additional_debit=0,
            mass_flux=flux,
            coverage_premise=True,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_EXACT_ALGEBRA,
            **_decl(route.component),
        )
        _expect_provider_error(
            lambda: conservation_provider(
                additional_debit=0,
                mass_flux=flux,
                coverage_premise=False,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_EXACT_ALGEBRA,
                **_decl(route.component),
            ),
            "coverage premise",
        )
    elif route.route_id == "explicit_debit":
        record = arithmetic_provider(
            additive_q=Q(1, 64),
            source=SOURCE_EXACT_ALGEBRA,
            status=PREMISE_STATUS_PROVEN,
            **_decl(route.component),
        )
        _expect_provider_error(
            lambda: arithmetic_provider(
                additive_q=Q(1, 64),
                source=SOURCE_EXACT_ALGEBRA,
                status=PREMISE_STATUS_PROVEN,
                measured_q=1,
                **_decl(route.component),
            ),
            "unexpected arguments",
        )
    elif route.route_id == "imp1_conditional" and route.component == "arithmetic":
        record = arithmetic_provider(
            additive_q=Q(1, 32),
            source=SOURCE_IMP1_ADMISSION_DEBIT,
            status=PREMISE_STATUS_CONDITIONAL,
            **_decl(route.component),
        )
        _expect_provider_error(
            lambda: arithmetic_provider(
                additive_q=0,
                source=SOURCE_IMP1_ADMISSION_DEBIT,
                status=PREMISE_STATUS_CONDITIONAL,
                **_decl(route.component),
            ),
            "zero stub",
        )
    else:
        raise Def1Stab1QualificationError(
            f"unregistered provider route {route.component}/{route.route_id}"
        )
    if not isinstance(record, ProviderRecord) or record.component != route.component:
        raise Def1Stab1QualificationError(
            f"{route.component}/{route.route_id} positive control did not return "
            "its provider record"
        )
    return ProviderRouteControl(
        component=route.component,
        route_id=route.route_id,
        conversion=route.conversion,
        positive_control_passed=True,
        injected_failure_refused=True,
        trajectory_values_evaluated=False,
        def1_booleans_evaluated=False,
        def1_boolean_names=DEF1_BOOLEAN_NAMES,
    )


@dataclass(frozen=True, slots=True)
class QualificationCoverageMatrix:
    """Fourteen-component route coverage; trajectory and DEF1 booleans stay open."""

    components: tuple[str, ...]
    routes: tuple[ProviderRouteControl, ...]
    provenance: QualificationProvenanceRecord
    conversion_identity: str
    trajectory_values_evaluated: bool
    def1_booleans_evaluated: bool
    def1_error_map_passed: bool
    global_pde_error_certified: bool
    claims_col1_values: bool
    used_measured_q: bool

    def __post_init__(self) -> None:
        if self.components != ERROR_BUDGET_COMPONENTS:
            raise Def1Stab1QualificationError(
                "coverage matrix must use the frozen fourteen-component order"
            )
        expected_ids = tuple((item.component, item.route_id) for item in PROVIDER_ROUTES)
        got_ids = tuple((item.component, item.route_id) for item in self.routes)
        if got_ids != expected_ids:
            raise Def1Stab1QualificationError(
                "coverage matrix omitted, aliased, or reordered provider routes"
            )
        seen = {item.component for item in self.routes}
        if seen != set(ERROR_BUDGET_COMPONENTS):
            missing = [name for name in ERROR_BUDGET_COMPONENTS if name not in seen]
            raise Def1Stab1QualificationError(
                f"coverage matrix missing components {missing}"
            )
        for item in self.routes:
            if not item.positive_control_passed or not item.injected_failure_refused:
                raise Def1Stab1QualificationError(
                    f"{item.component}/{item.route_id} lacks positive and "
                    "injected-failure controls"
                )
            if item.trajectory_values_evaluated or item.def1_booleans_evaluated:
                raise Def1Stab1QualificationError(
                    "coverage matrix cannot evaluate trajectory values or DEF1 booleans"
                )
        if (
            self.trajectory_values_evaluated
            or self.def1_booleans_evaluated
            or self.def1_error_map_passed
            or self.global_pde_error_certified
            or self.claims_col1_values
            or self.used_measured_q
        ):
            raise Def1Stab1QualificationError(
                "coverage matrix cannot promote DEF1, certify a PDE error, or "
                "claim COL1/trajectory values"
            )
        if self.conversion_identity != conversion_identity_digest(
            CONVERSION_COVERAGE_MATRIX
        ):
            raise Def1Stab1QualificationError(
                "coverage conversion_identity does not match the live digest"
            )


def execute_qualification_coverage_matrix() -> QualificationCoverageMatrix:
    routes = tuple(_run_provider_route(route) for route in PROVIDER_ROUTES)
    return QualificationCoverageMatrix(
        components=ERROR_BUDGET_COMPONENTS,
        routes=routes,
        provenance=bind_qualification_provenance(
            artifact_id="FGC-1-DEF1-STAB1",
            owner_path=OWNER_PATH_PROVIDERS,
            conversion_name=CONVERSION_COVERAGE_MATRIX,
        ),
        conversion_identity=conversion_identity_digest(CONVERSION_COVERAGE_MATRIX),
        trajectory_values_evaluated=False,
        def1_booleans_evaluated=False,
        def1_error_map_passed=False,
        global_pde_error_certified=False,
        claims_col1_values=False,
        used_measured_q=False,
    )


__all__ = [
    "ADM_METRIC_TWO_JET_ORDER",
    "ADM_SOURCE_FIRST_ORDER_STATE_ORDER",
    "ADM_SOURCE_PHYSICAL_ARGUMENT_ORDER",
    "BASE_METRIC_TWO_JET_ORDER",
    "BaseToAdmTwoJetInverse",
    "CONVERSION_BASE_TO_ADM",
    "CONVERSION_COVERAGE_MATRIX",
    "CONVERSION_FO1_ADM_GEOMETRY",
    "CONVERSION_IMP1_TO_Q",
    "CONVERSION_RED1_MHG2",
    "Def1Stab1QualificationError",
    "EquationOrderIdentityMap",
    "Fo1AdmGeometrySlotMap",
    "Fo1GeometryConversion",
    "IMP1_18_CHANNEL_ORDER",
    "Imp1ToQConversion",
    "INVENTORY_ADM_METRIC_TWO_JET",
    "INVENTORY_ADM_SOURCE_PHYSICAL",
    "INVENTORY_ADM_SOURCE_STATE",
    "INVENTORY_FO1_BASE_PHYSICAL",
    "INVENTORY_FO1_BASE_STATE",
    "INVENTORY_FO1_BASE_TO_ADM",
    "PROVIDER_ROUTES",
    "QualificationCoverageMatrix",
    "QualificationProvenanceRecord",
    "RED1_TO_MHG2_EQUATION_PAIRS",
    "SLOT_EXACT_INVERSE",
    "bind_qualification_provenance",
    "convert_base_to_adm_geometry_slots",
    "convert_fo1_to_geometry_slots",
    "convert_imp1_channels_to_q",
    "invert_base_metric_two_jets",
    "invert_fo1_physical_arguments_to_adm",
    "invert_ordered_base_metric_two_jets",
    "execute_qualification_coverage_matrix",
    "fo1_adm_geometry_slot_map",
    "geometry_slot_zero_justification",
    "imp1_channel_geometry_assignments",
    "live_imp1_18_channel_order",
    "map_equation_order_to_mhg2",
    "owner_file_sha256",
    "phi_chi_geometry_status",
    "red1_mhg2_equation_identity_map",
]
