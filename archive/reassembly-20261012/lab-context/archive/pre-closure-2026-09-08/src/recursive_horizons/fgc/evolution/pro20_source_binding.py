"""Configuration-bound GR-0 callbacks for the prospective PRO20 factory.

This module observes supplied objects; construction performs no source call,
initial-data solve, file I/O, or origin authentication. The actual source is
still Proto12/SRC4 and the projector is the preserved GR-0 factory closure.
An immutable byte image includes solver controls omitted by the inherited
HLT16 template projection. Callback use rechecks that image before and after
delegation. This local binding is not a source-image or execution authority.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
from math import isfinite
from numbers import Real
from types import CodeType, FunctionType
from typing import Any

import numpy as np

from recursive_horizons.evidence_io import canonical_json_bytes, load_canonical_json

from . import gr0_calibration
from .numerical_engine import (
    EvolutionRHS,
    EvolutionState,
    SBPFirstDerivative,
    UniformRadialGrid,
    array_content_sha256,
)
from .proto12_runtime import Proto12GR0EvolutionOperator


SCHEMA = "FGC-1-PRO20-source-configuration-binding-v1"
MAX_POINT_COUNT = 16385
OPERATOR_FIELDS = frozenset({
    "grid", "derivative", "ko_dissipation", "raw_tolerance",
    "kinetic_condition_maximum", "maximum_refinement_iterations", "point_batch_size",
})


class Pro20SourceBindingError(ValueError):
    """A supplied callback or its captured configuration changed or is unsupported."""


def _integer(value: object, label: str, *, minimum: int, maximum: int) -> int:
    if type(value) is not int or not minimum <= value <= maximum:
        raise Pro20SourceBindingError(f"{label} must be an integer in [{minimum}, {maximum}]")
    return value


def _number(value: object, label: str, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise Pro20SourceBindingError(f"{label} must be finite binary64")
    answer = float(value)
    if not isfinite(answer) or (positive and answer <= 0):
        raise Pro20SourceBindingError(f"{label} is outside its finite domain")
    return answer


def _code_atom(value: object, *, depth: int) -> object:
    if depth > 16:
        raise Pro20SourceBindingError("callback code identity is too deeply nested")
    if value is None or type(value) in (bool, str, int):
        return value
    if type(value) is float:
        return {"binary64_hex": _number(value, "code constant").hex()}
    if type(value) is bytes:
        return {"bytes_hex": value.hex()}
    if type(value) is tuple:
        return {"tuple": [_code_atom(item, depth=depth + 1) for item in value]}
    if type(value) is frozenset:
        encoded = [_code_atom(item, depth=depth + 1) for item in value]
        return {"frozenset": sorted(encoded, key=repr)}
    if type(value) is slice:
        return {"slice": [
            _code_atom(value.start, depth=depth + 1),
            _code_atom(value.stop, depth=depth + 1),
            _code_atom(value.step, depth=depth + 1),
        ]}
    if value is Ellipsis:
        return {"ellipsis": True}
    if type(value) is CodeType:
        return {"code": _code_identity(value, depth=depth + 1)}
    raise Pro20SourceBindingError(
        "callback code has an unsupported constant "
        f"{type(value).__module__}.{type(value).__qualname__}"
    )


def _code_identity(code: object, *, depth: int = 0) -> dict[str, object]:
    if type(code) is not CodeType:
        raise Pro20SourceBindingError("callback has no inspectable Python code")
    return {
        "bytecode_sha256": sha256(code.co_code).hexdigest(),
        "exception_table_hex": code.co_exceptiontable.hex(),
        "constants": [_code_atom(item, depth=depth + 1) for item in code.co_consts],
        "names": list(code.co_names),
        "variables": list(code.co_varnames),
        "free_variables": list(code.co_freevars),
        "cell_variables": list(code.co_cellvars),
        "argcount": code.co_argcount,
        "positional_only_argcount": code.co_posonlyargcount,
        "keyword_only_argcount": code.co_kwonlyargcount,
        "local_count": code.co_nlocals,
        "flags": code.co_flags,
    }


def _function_identity(function: object) -> dict[str, object]:
    if type(function) is not FunctionType:
        raise Pro20SourceBindingError("callback is not a Python function")
    if function.__defaults__ is not None or function.__kwdefaults__ is not None:
        raise Pro20SourceBindingError("callback defaults differ from the preserved kernel")
    return {
        "module": function.__module__,
        "qualname": function.__qualname__,
        "code": _code_identity(function.__code__),
    }


def _grid_material(grid: object) -> dict[str, object]:
    if type(grid) is not UniformRadialGrid:
        raise Pro20SourceBindingError("source grid must be the preserved UniformRadialGrid")
    count = _integer(grid.point_count, "point_count", minimum=9, maximum=MAX_POINT_COUNT)
    lower = _number(grid.minimum, "grid minimum")
    upper = _number(grid.maximum, "grid maximum")
    if lower != 0.0 or upper <= lower:
        raise Pro20SourceBindingError("source requires a positive-width center grid")
    spacing = _number(grid.spacing, "grid spacing", positive=True)
    coordinates = grid.coordinates
    if coordinates.shape != (count,) or not np.all(np.diff(coordinates) > 0):
        raise Pro20SourceBindingError("source coordinates differ from the ordered grid")
    return {
        "minimum_hex": lower.hex(),
        "maximum_hex": upper.hex(),
        "spacing_hex": spacing.hex(),
        "point_count": count,
        "coordinates_sha256": array_content_sha256(coordinates),
    }


def _operator_material(operator: object) -> dict[str, object]:
    if type(operator) is not Proto12GR0EvolutionOperator:
        raise Pro20SourceBindingError("source must be the preserved Proto12 GR-0 operator")
    if set(vars(operator)) != OPERATOR_FIELDS:
        raise Pro20SourceBindingError("source configuration has missing or extra fields")
    grid = _grid_material(operator.grid)
    derivative = operator.derivative
    if type(derivative) is not SBPFirstDerivative:
        raise Pro20SourceBindingError("source derivative type differs")
    if _grid_material(derivative.grid) != grid:
        raise Pro20SourceBindingError("source and derivative grids differ")
    if type(derivative.order) is not int or derivative.order not in (2, 4):
        raise Pro20SourceBindingError("source derivative order differs")
    ko = _number(operator.ko_dissipation, "KO dissipation")
    if ko < 0:
        raise Pro20SourceBindingError("KO dissipation must be nonnegative")
    call = type(operator).__call__
    if (
        type(call) is not FunctionType
        or call.__module__ != Proto12GR0EvolutionOperator.__module__
        or call.__qualname__ != "Proto12GR0EvolutionOperator.__call__"
        or call.__closure__ is not None
    ):
        raise Pro20SourceBindingError("source entrypoint differs from the preserved kernel")
    return {
        "type": f"{type(operator).__module__}.{type(operator).__qualname__}",
        "entrypoint": _function_identity(call),
        "grid": grid,
        "derivative_type": f"{type(derivative).__module__}.{type(derivative).__qualname__}",
        "spatial_order": derivative.order,
        "norm_weights_sha256": array_content_sha256(derivative.norm_weights),
        "stencil_reach_intervals": derivative.stencil_reach_intervals,
        "ko_dissipation_hex": ko.hex(),
        "raw_tolerance_hex": _number(operator.raw_tolerance, "source tolerance", positive=True).hex(),
        "kinetic_condition_maximum_hex": _number(
            operator.kinetic_condition_maximum, "kinetic condition maximum", positive=True
        ).hex(),
        "maximum_refinement_iterations": _integer(
            operator.maximum_refinement_iterations, "source refinement iterations", minimum=0, maximum=16
        ),
        "point_batch_size": _integer(
            operator.point_batch_size, "source point-batch size", minimum=1, maximum=MAX_POINT_COUNT
        ),
    }


def _projector_material(projector: object, point_count: int) -> dict[str, object]:
    if type(projector) is not FunctionType:
        raise Pro20SourceBindingError("projector must be the preserved factory function")
    expected = tuple(
        item for item in gr0_calibration.make_gr0_center_boundary_projector.__code__.co_consts
        if type(item) is CodeType and item.co_name == "projector"
    )
    if (
        len(expected) != 1
        or projector.__code__ is not expected[0]
        or projector.__globals__ is not gr0_calibration.__dict__
        or projector.__module__ != gr0_calibration.__name__
        or projector.__qualname__ != "make_gr0_center_boundary_projector.<locals>.projector"
    ):
        raise Pro20SourceBindingError("projector code/global owner differs")
    if projector.__closure__ is None or len(projector.__closure__) != 2:
        raise Pro20SourceBindingError("projector closure differs")
    try:
        cells = dict(zip(projector.__code__.co_freevars,
                         (cell.cell_contents for cell in projector.__closure__), strict=True))
    except (ValueError, TypeError) as error:
        raise Pro20SourceBindingError("projector closure is incomplete") from error
    if set(cells) != {"fixed_outer_rows", "reference"}:
        raise Pro20SourceBindingError("projector closure fields differ")
    rows = _integer(cells["fixed_outer_rows"], "fixed outer rows", minimum=1,
                    maximum=point_count // 2 - 1)
    reference = cells["reference"]
    if type(reference) is not EvolutionState or reference.shape != (point_count, 6):
        raise Pro20SourceBindingError("projector reference state shape/type differs")
    if gr0_calibration.EvolutionState is not EvolutionState:
        raise Pro20SourceBindingError("projector state-class global changed")
    for name in ("ADM_CENTER_PARITIES", "Q_CENTER_PARITIES"):
        parity = getattr(gr0_calibration, name)
        if not isinstance(parity, np.ndarray) or parity.shape != (6,) or not np.all(np.abs(parity) == 1):
            raise Pro20SourceBindingError("projector parity globals differ")
    return {
        "entrypoint": _function_identity(projector),
        "fixed_outer_rows": rows,
        "reference_shape": list(reference.shape),
        "reference_state_sha256": array_content_sha256(reference.u, reference.p, reference.q),
        "adm_center_parities_sha256": array_content_sha256(gr0_calibration.ADM_CENTER_PARITIES),
        "q_center_parities_sha256": array_content_sha256(gr0_calibration.Q_CENTER_PARITIES),
    }


def _capture_material(operator: object, projector: object) -> dict[str, object]:
    source = _operator_material(operator)
    return {
        "schema": SCHEMA,
        "operator": source,
        "projector": _projector_material(projector, source["grid"]["point_count"]),
        "physical_source_authenticated": False,
        "source_manifest_authenticated": False,
        "historical_origin_authenticated": False,
        "campaign_execution_authorized": False,
    }


def _evaluate_rhs(operator: Proto12GR0EvolutionOperator, time: float, state: EvolutionState) -> EvolutionRHS:
    """Single forwarding point; qualification calls the real preserved operator."""
    return operator(time, state)


def _state_hash(state: object, point_count: int) -> str:
    if type(state) is not EvolutionState or state.shape != (point_count, 6):
        raise Pro20SourceBindingError("callback input state shape/type differs")
    return array_content_sha256(state.u, state.p, state.q)


@dataclass(frozen=True, slots=True, eq=False)
class Pro20GR0SourceBinding:
    """Provided GR-0 callbacks bound to an immutable observed byte image.

    This constructor does not authenticate that image or its inputs. The
    future factory/qualification/authority must establish that separately.
    There is no caller-supplied success flag or caller-supplied image hash.
    """

    operator: Proto12GR0EvolutionOperator = field(repr=False)
    projector: FunctionType = field(repr=False)
    _captured_bytes: bytes = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if type(self) is not Pro20GR0SourceBinding:
            raise Pro20SourceBindingError("source binding must not be subclassed")
        object.__setattr__(self, "_captured_bytes", canonical_json_bytes(
            _capture_material(self.operator, self.projector)
        ))

    @property
    def captured_bytes(self) -> bytes:
        return self._captured_bytes

    @property
    def configuration_sha256(self) -> str:
        return sha256(self._captured_bytes).hexdigest()

    def as_mapping(self) -> dict[str, Any]:
        return load_canonical_json(self._captured_bytes)

    def validate(self) -> None:
        try:
            observed = canonical_json_bytes(_capture_material(self.operator, self.projector))
        except Pro20SourceBindingError as error:
            raise Pro20SourceBindingError(
                "provided GR-0 source/projector configuration changed"
            ) from error
        if observed != self._captured_bytes:
            raise Pro20SourceBindingError("provided GR-0 source/projector configuration changed")

    def hlt17_synthetic_callable_binding(self) -> dict[str, object]:
        """Use HLT17's local fingerprint seam without claiming physical authentication."""
        self.validate()
        return {
            "kind": "pro20_gr0_configuration_bound_callbacks",
            "configuration_schema": SCHEMA,
            "configuration_sha256": self.configuration_sha256,
            "physical_source_authenticated": False,
            "source_manifest_authenticated": False,
        }

    def _before(self, time: Real, state: EvolutionState) -> tuple[float, str]:
        self.validate()
        return _number(time, "callback time"), _state_hash(state, self.operator.grid.point_count)

    def _after(self, state: EvolutionState, before: str) -> None:
        self.validate()
        if _state_hash(state, self.operator.grid.point_count) != before:
            raise Pro20SourceBindingError("provided callback modified its input state")

    def rhs(self, time: Real, state: EvolutionState) -> EvolutionRHS:
        checked_time, before = self._before(time, state)
        try:
            result = _evaluate_rhs(self.operator, checked_time, state)
        finally:
            self._after(state, before)
        if type(result) is not EvolutionRHS or result.shape != state.shape:
            raise Pro20SourceBindingError("provided source returned an incompatible RHS")
        return result

    def project(self, time: Real, state: EvolutionState) -> EvolutionState:
        checked_time, before = self._before(time, state)
        try:
            result = self.projector(checked_time, state)
        finally:
            self._after(state, before)
        if type(result) is not EvolutionState or result.shape != state.shape:
            raise Pro20SourceBindingError("provided projector returned an incompatible state")
        return result


__all__ = ["SCHEMA", "Pro20GR0SourceBinding", "Pro20SourceBindingError"]
