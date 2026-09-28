"""Small, auditable method-of-lines engine for the scoped spherical study.

The engine owns only numerical mechanics: uniform radial grids, the two
predeclared SBP derivative operators, parity projection at the regular centre,
RK4/SSPRK3 stage construction, fail-closed step transactions, and canonical
restart checkpoints.  The physical right-hand side and scientific monitors
are injected explicitly.  This keeps the numerical validation controls from
silently defining a second set of FGC equations.

Arrays use ``(radial_point, field)`` order.  A second-order system is stored as
``(u, p, q)`` with ``p=partial_t u`` and ``q=partial_r u``.  The engine evolves
``partial_t u=p`` and leaves the supplied right-hand side responsible for
``partial_t p`` and ``partial_t q``.  Production FGC use sets
``partial_t q=D_r p`` and monitors ``q-D_r u`` as a reduction constraint.

This module is numerical infrastructure.  It is not an FGC trajectory,
collapse result, affine-defocusing result, or retained-EFT authorization.
"""

from __future__ import annotations

from dataclasses import dataclass
import base64
import hashlib
import json
from math import isfinite
from numbers import Real
from pathlib import Path
import struct
from typing import Callable, Mapping, Protocol, Sequence

import numpy as np


PRIMARY_METHOD = "fourth_order_diagonal_norm_SBP_with_second_order_boundary_closure_plus_RK4"
COMPARATOR_METHOD = "second_order_diagonal_norm_SBP_plus_SSPRK3"
METHODS = (PRIMARY_METHOD, COMPARATOR_METHOD)


def _finite(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    return answer


def _finite_array(name: str, value: object, *, ndim: int | None = None) -> np.ndarray:
    answer = np.asarray(value, dtype=np.float64)
    if ndim is not None and answer.ndim != ndim:
        raise ValueError(f"{name} must have {ndim} dimensions")
    if not np.all(np.isfinite(answer)):
        raise ValueError(f"{name} must contain only finite binary64 values")
    return np.ascontiguousarray(answer)


@dataclass(frozen=True, slots=True)
class UniformRadialGrid:
    minimum: float
    maximum: float
    point_count: int

    def __post_init__(self) -> None:
        lower = _finite("minimum", self.minimum)
        upper = _finite("maximum", self.maximum)
        count = self.point_count
        if isinstance(count, bool) or not isinstance(count, int) or count < 9:
            raise ValueError("point_count must be an integer at least nine")
        if upper <= lower:
            raise ValueError("maximum must exceed minimum")
        object.__setattr__(self, "minimum", lower)
        object.__setattr__(self, "maximum", upper)

    @property
    def spacing(self) -> float:
        return (self.maximum - self.minimum) / (self.point_count - 1)

    @property
    def coordinates(self) -> np.ndarray:
        return np.linspace(
            self.minimum, self.maximum, self.point_count, dtype=np.float64
        )


@dataclass(frozen=True, slots=True)
class SBPFirstDerivative:
    """The frozen diagonal-norm 4-2 or 2-1 first-derivative operator."""

    grid: UniformRadialGrid
    order: int

    def __post_init__(self) -> None:
        if self.order not in {2, 4}:
            raise ValueError("SBP order must be two or four")
        if self.order == 4 and self.grid.point_count < 9:
            raise ValueError("the 4-2 SBP operator requires at least nine points")

    @property
    def stencil_reach_intervals(self) -> int:
        # The 4-2 boundary closure reaches three points away; the 2-1
        # comparator reaches one.  BND2 consumes these exact integers.
        return 3 if self.order == 4 else 1

    @property
    def norm_weights(self) -> np.ndarray:
        count = self.grid.point_count
        if self.order == 2:
            weights = np.ones(count, dtype=np.float64)
            weights[[0, -1]] = 0.5
        else:
            weights = np.ones(count, dtype=np.float64)
            boundary = np.asarray((17 / 48, 59 / 48, 43 / 48, 49 / 48))
            weights[:4] = boundary
            weights[-4:] = boundary[::-1]
        return self.grid.spacing * weights

    def _differentiate_without_parity(self, values: np.ndarray) -> np.ndarray:
        h = self.grid.spacing
        output = np.empty_like(values)
        if self.order == 2:
            output[0] = (values[1] - values[0]) / h
            output[1:-1] = (values[2:] - values[:-2]) / (2 * h)
            output[-1] = (values[-1] - values[-2]) / h
            return output

        output[0] = (
            -24 / 17 * values[0]
            + 59 / 34 * values[1]
            - 4 / 17 * values[2]
            - 3 / 34 * values[3]
        ) / h
        output[1] = (-values[0] + values[2]) / (2 * h)
        output[2] = (
            4 / 43 * values[0]
            - 59 / 86 * values[1]
            + 59 / 86 * values[3]
            - 4 / 43 * values[4]
        ) / h
        output[3] = (
            3 / 98 * values[0]
            - 59 / 98 * values[2]
            + 32 / 49 * values[4]
            - 4 / 49 * values[5]
        ) / h
        output[4:-4] = (
            values[2:-6]
            - 8 * values[3:-5]
            + 8 * values[5:-3]
            - values[6:-2]
        ) / (12 * h)
        output[-4] = -(
            3 / 98 * values[-1]
            - 59 / 98 * values[-3]
            + 32 / 49 * values[-5]
            - 4 / 49 * values[-6]
        ) / h
        output[-3] = -(
            4 / 43 * values[-1]
            - 59 / 86 * values[-2]
            + 59 / 86 * values[-4]
            - 4 / 43 * values[-5]
        ) / h
        output[-2] = -( -values[-1] + values[-3]) / (2 * h)
        output[-1] = -(
            -24 / 17 * values[-1]
            + 59 / 34 * values[-2]
            - 4 / 17 * values[-3]
            - 3 / 34 * values[-4]
        ) / h
        return output

    def differentiate(
        self,
        values: object,
        *,
        center_parities: Sequence[int] | None = None,
    ) -> np.ndarray:
        """Differentiate along the radial axis.

        Without ``center_parities`` this is the exact diagonal-norm SBP
        operator at both ends.  On a grid beginning at ``r=0``, supplying one
        parity per field replaces only the inner coordinate-boundary rows by
        the smooth even/odd extension through the centre.  The outer rows stay
        the frozen SBP closure.  This is the regular-centre formulation used
        by the spherical engine; it is not described as a two-boundary SBP
        identity after that projection.
        """

        data = _finite_array("values", values)
        if data.ndim not in {1, 2} or data.shape[0] != self.grid.point_count:
            raise ValueError("values must have radial leading axis matching the grid")
        scalar = data.ndim == 1
        work = data[:, None] if scalar else data
        output = self._differentiate_without_parity(work)
        if center_parities is not None:
            if self.grid.minimum != 0.0:
                raise ValueError("centre parity requires a grid beginning exactly at zero")
            parity = np.asarray(center_parities, dtype=np.int64)
            if parity.shape != (work.shape[1],) or np.any(np.abs(parity) != 1):
                raise ValueError("center_parities must contain one +/-1 per field")
            if np.any((parity == -1) & (work[0] != 0.0)):
                raise ValueError("odd fields must vanish exactly at the centre")
            h = self.grid.spacing
            if self.order == 2:
                output[0] = (work[1] - parity * work[1]) / (2 * h)
            else:
                # Fourth-order centered stencils with negative-radius samples
                # supplied by f(-r)=parity*f(r).
                output[0] = (
                    parity * work[2]
                    - 8 * parity * work[1]
                    + 8 * work[1]
                    - work[2]
                ) / (12 * h)
                output[1] = (
                    parity * work[1]
                    - 8 * work[0]
                    + 8 * work[2]
                    - work[3]
                ) / (12 * h)
        return output[:, 0] if scalar else output

    def dense_matrix(self) -> np.ndarray:
        """Return the unprojected matrix for exact SBP-identity tests."""

        identity = np.eye(self.grid.point_count, dtype=np.float64)
        return self.differentiate(identity)

    def sbp_identity_residual(self) -> float:
        derivative = self.dense_matrix()
        norm = np.diag(self.norm_weights)
        boundary = np.zeros_like(derivative)
        boundary[0, 0] = -1.0
        boundary[-1, -1] = 1.0
        return float(np.max(np.abs(norm @ derivative + derivative.T @ norm - boundary)))

    def kreiss_oliger_dissipation(
        self,
        values: object,
        *,
        coefficient: Real,
        center_parities: Sequence[int] | None = None,
    ) -> np.ndarray:
        """Return the frozen high-order KO filter contribution to ``d/dt``.

        The primary uses the seven-point sixth-difference filter and the
        comparator the five-point fourth-difference filter.  At a regular
        centre, negative-radius samples are supplied by parity.  The outer
        closure is set to zero over the filter reach; the independently
        causal-isolated measurement region never consumes those closure
        points.
        """

        data = _finite_array("values", values)
        if data.ndim not in {1, 2} or data.shape[0] != self.grid.point_count:
            raise ValueError("values must have radial leading axis matching the grid")
        epsilon = _finite("coefficient", coefficient)
        if epsilon < 0.0:
            raise ValueError("dissipation coefficient must be nonnegative")
        scalar = data.ndim == 1
        work = data[:, None] if scalar else data
        output = np.zeros_like(work)
        if epsilon == 0.0:
            return output[:, 0] if scalar else output
        parity = None
        if center_parities is not None:
            if self.grid.minimum != 0.0:
                raise ValueError("centre parity requires a grid beginning at zero")
            parity = np.asarray(center_parities, dtype=np.int64)
            if parity.shape != (work.shape[1],) or np.any(np.abs(parity) != 1):
                raise ValueError("center_parities must contain one +/-1 per field")
            if np.any((parity == -1) & (work[0] != 0.0)):
                raise ValueError("odd fields must vanish exactly at the centre")

        def sample(index: int) -> np.ndarray:
            if index >= 0:
                return work[index]
            if parity is None:
                raise IndexError("unprojected KO stencil crossed the inner boundary")
            return parity * work[-index]

        h = self.grid.spacing
        if self.order == 4:
            weights = (1.0, -6.0, 15.0, -20.0, 15.0, -6.0, 1.0)
            reach = 3
            scale = epsilon / (64.0 * h)
        else:
            weights = (1.0, -4.0, 6.0, -4.0, 1.0)
            reach = 2
            scale = -epsilon / (16.0 * h)
        first = 0 if parity is not None else reach
        for index in range(first, self.grid.point_count - reach):
            output[index] = scale * sum(
                weight * sample(index + offset)
                for weight, offset in zip(weights, range(-reach, reach + 1), strict=True)
            )
        return output[:, 0] if scalar else output


def parity_project(values: object, parities: Sequence[int]) -> np.ndarray:
    """Return a copy satisfying the exact centre value parity."""

    data = _finite_array("values", values, ndim=2).copy()
    parity = np.asarray(parities, dtype=np.int64)
    if parity.shape != (data.shape[1],) or np.any(np.abs(parity) != 1):
        raise ValueError("parities must contain one +/-1 per field")
    data[0, parity == -1] = 0.0
    return data


@dataclass(frozen=True, slots=True)
class EvolutionState:
    u: np.ndarray
    p: np.ndarray
    q: np.ndarray

    def __post_init__(self) -> None:
        arrays = tuple(
            _finite_array(name, value, ndim=2)
            for name, value in (("u", self.u), ("p", self.p), ("q", self.q))
        )
        if arrays[0].shape != arrays[1].shape or arrays[0].shape != arrays[2].shape:
            raise ValueError("u, p, and q must have identical shapes")
        if arrays[0].shape[0] < 9 or arrays[0].shape[1] < 1:
            raise ValueError("evolution state must contain at least nine points and one field")
        for name, value in zip(("u", "p", "q"), arrays, strict=True):
            immutable = value.copy()
            immutable.setflags(write=False)
            object.__setattr__(self, name, immutable)

    @property
    def shape(self) -> tuple[int, int]:
        return self.u.shape

    def linear_combination(
        self,
        terms: Sequence[tuple[float, "EvolutionRHS"]],
    ) -> "EvolutionState":
        du = np.zeros_like(self.u)
        dp = np.zeros_like(self.p)
        dq = np.zeros_like(self.q)
        for coefficient, rhs in terms:
            scale = _finite("linear-combination coefficient", coefficient)
            if not isinstance(rhs, EvolutionRHS) or rhs.shape != self.shape:
                raise ValueError("linear-combination RHS shape differs")
            du += scale * rhs.du
            dp += scale * rhs.dp
            dq += scale * rhs.dq
        return EvolutionState(self.u + du, self.p + dp, self.q + dq)


@dataclass(frozen=True, slots=True)
class EvolutionRHS:
    du: np.ndarray
    dp: np.ndarray
    dq: np.ndarray
    diagnostics: Mapping[str, object]

    def __post_init__(self) -> None:
        arrays = tuple(
            _finite_array(name, value, ndim=2)
            for name, value in (("du", self.du), ("dp", self.dp), ("dq", self.dq))
        )
        if arrays[0].shape != arrays[1].shape or arrays[0].shape != arrays[2].shape:
            raise ValueError("du, dp, and dq must have identical shapes")
        if not isinstance(self.diagnostics, Mapping):
            raise TypeError("RHS diagnostics must be a mapping")
        for name, value in zip(("du", "dp", "dq"), arrays, strict=True):
            immutable = value.copy()
            immutable.setflags(write=False)
            object.__setattr__(self, name, immutable)
        object.__setattr__(self, "diagnostics", dict(self.diagnostics))

    @property
    def shape(self) -> tuple[int, int]:
        return self.du.shape


class RightHandSide(Protocol):
    def __call__(self, time: float, state: EvolutionState) -> EvolutionRHS: ...


StateProjector = Callable[[float, EvolutionState], EvolutionState]


@dataclass(frozen=True, slots=True)
class StageRecord:
    stage_name: str
    time: float
    state: EvolutionState
    rhs: EvolutionRHS


@dataclass(frozen=True, slots=True)
class StepProposal:
    method: str
    initial_time: float
    final_time: float
    initial_state: EvolutionState
    candidate_state: EvolutionState
    stages: tuple[StageRecord, ...]


class ProposalGuard(Protocol):
    def __call__(self, proposal: StepProposal) -> Mapping[str, object] | None: ...


@dataclass(frozen=True, slots=True)
class AcceptedStep:
    time: float
    step_index: int
    transaction_serial: int
    state: EvolutionState
    guard_receipts: tuple[Mapping[str, object], ...]
    stage_records: tuple[StageRecord, ...]


def _project(
    projector: StateProjector | None,
    time: float,
    state: EvolutionState,
) -> EvolutionState:
    if projector is None:
        return state
    projected = projector(time, state)
    if not isinstance(projected, EvolutionState) or projected.shape != state.shape:
        raise ValueError("state projector returned an incompatible state")
    return projected


def propose_step(
    *,
    method: str,
    time: Real,
    step_size: Real,
    state: EvolutionState,
    rhs: RightHandSide,
    projector: StateProjector | None = None,
) -> StepProposal:
    """Build every internal stage and an explicit candidate endpoint."""

    if method not in METHODS:
        raise ValueError("unknown evolution method")
    if not isinstance(state, EvolutionState):
        raise TypeError("state must be EvolutionState")
    start = _finite("time", time)
    dt = _finite("step_size", step_size)
    if dt <= 0.0:
        raise ValueError("step_size must be positive")

    records: list[StageRecord] = []

    def evaluate(name: str, stage_time: float, stage_state: EvolutionState) -> EvolutionRHS:
        projected = _project(projector, stage_time, stage_state)
        value = rhs(stage_time, projected)
        if not isinstance(value, EvolutionRHS) or value.shape != state.shape:
            raise ValueError("right-hand side returned an incompatible value")
        records.append(StageRecord(name, stage_time, projected, value))
        return value

    if method == PRIMARY_METHOD:
        k1 = evaluate("rk4_k1", start, state)
        y2 = _project(projector, start + dt / 2, state.linear_combination(((dt / 2, k1),)))
        k2 = evaluate("rk4_k2", start + dt / 2, y2)
        y3 = _project(projector, start + dt / 2, state.linear_combination(((dt / 2, k2),)))
        k3 = evaluate("rk4_k3", start + dt / 2, y3)
        y4 = _project(projector, start + dt, state.linear_combination(((dt, k3),)))
        k4 = evaluate("rk4_k4", start + dt, y4)
        candidate = state.linear_combination(
            (
                (dt / 6, k1),
                (dt / 3, k2),
                (dt / 3, k3),
                (dt / 6, k4),
            )
        )
    else:
        k1 = evaluate("ssprk3_s0", start, state)
        y1 = _project(projector, start + dt, state.linear_combination(((dt, k1),)))
        k2 = evaluate("ssprk3_s1", start + dt, y1)
        euler_y1 = y1.linear_combination(((dt, k2),))
        y2 = EvolutionState(
            0.75 * state.u + 0.25 * euler_y1.u,
            0.75 * state.p + 0.25 * euler_y1.p,
            0.75 * state.q + 0.25 * euler_y1.q,
        )
        y2 = _project(projector, start + dt / 2, y2)
        k3 = evaluate("ssprk3_s2", start + dt / 2, y2)
        euler_y2 = y2.linear_combination(((dt, k3),))
        candidate = EvolutionState(
            state.u / 3 + 2 * euler_y2.u / 3,
            state.p / 3 + 2 * euler_y2.p / 3,
            state.q / 3 + 2 * euler_y2.q / 3,
        )

    final_time = start + dt
    candidate = _project(projector, final_time, candidate)
    # The endpoint is always evaluated independently: RK4 k4 and SSPRK3 s2
    # are not the accepted Runge--Kutta combination.
    evaluate("candidate_endpoint", final_time, candidate)
    return StepProposal(
        method=method,
        initial_time=start,
        final_time=final_time,
        initial_state=state,
        candidate_state=candidate,
        stages=tuple(records),
    )


def accept_step(
    proposal: StepProposal,
    *,
    previous_step_index: int,
    previous_transaction_serial: int,
    guards: Sequence[ProposalGuard] = (),
) -> AcceptedStep:
    """Run every guard before exposing the proposed state as accepted."""

    if not isinstance(proposal, StepProposal):
        raise TypeError("proposal must be StepProposal")
    if previous_step_index < 0 or previous_transaction_serial < 0:
        raise ValueError("previous indices must be nonnegative")
    receipts: list[Mapping[str, object]] = []
    for guard in guards:
        receipt = guard(proposal)
        receipts.append({} if receipt is None else dict(receipt))
    return AcceptedStep(
        time=proposal.final_time,
        step_index=previous_step_index + 1,
        transaction_serial=previous_transaction_serial + len(proposal.stages),
        state=proposal.candidate_state,
        guard_receipts=tuple(receipts),
        stage_records=proposal.stages,
    )


def evolve_fixed_steps(
    *,
    method: str,
    initial_state: EvolutionState,
    initial_time: Real,
    step_size: Real,
    step_count: int,
    rhs: RightHandSide,
    projector: StateProjector | None = None,
    guard_factory: Callable[[int, int], Sequence[ProposalGuard]] | None = None,
    retain_every: int = 1,
) -> tuple[AcceptedStep, ...]:
    """Evolve a bounded number of steps and retain immutable endpoint records."""

    if isinstance(step_count, bool) or not isinstance(step_count, int) or step_count < 1:
        raise ValueError("step_count must be a positive integer")
    if isinstance(retain_every, bool) or not isinstance(retain_every, int) or retain_every < 1:
        raise ValueError("retain_every must be a positive integer")
    state = initial_state
    time = _finite("initial_time", initial_time)
    dt = _finite("step_size", step_size)
    step_index = 0
    serial = 0
    retained: list[AcceptedStep] = []
    for _ in range(step_count):
        proposal = propose_step(
            method=method,
            time=time,
            step_size=dt,
            state=state,
            rhs=rhs,
            projector=projector,
        )
        guards = () if guard_factory is None else tuple(guard_factory(step_index, serial))
        accepted = accept_step(
            proposal,
            previous_step_index=step_index,
            previous_transaction_serial=serial,
            guards=guards,
        )
        state = accepted.state
        time = accepted.time
        step_index = accepted.step_index
        serial = accepted.transaction_serial
        if step_index % retain_every == 0 or step_index == step_count:
            retained.append(accepted)
    return tuple(retained)


@dataclass(frozen=True, slots=True)
class EvolutionCheckpoint:
    """Canonical, checksum-protected restart state.

    The JSON representation uses little-endian binary64 bytes rather than
    decimal lists, so a decode/encode round trip is bitwise stable.  Scientific
    metadata must remain JSON scalar/list/mapping data; opaque Python objects
    are intentionally forbidden.
    """

    time: float
    step_index: int
    transaction_serial: int
    state: EvolutionState
    metadata: Mapping[str, object]

    def __post_init__(self) -> None:
        object.__setattr__(self, "time", _finite("checkpoint time", self.time))
        for name in ("step_index", "transaction_serial"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f"{name} must be a nonnegative integer")
        if not isinstance(self.state, EvolutionState):
            raise TypeError("checkpoint state must be EvolutionState")
        if not isinstance(self.metadata, Mapping):
            raise TypeError("checkpoint metadata must be a mapping")
        # Serialization is also the schema validator for nested metadata.
        json.dumps(self.metadata, sort_keys=True, separators=(",", ":"), allow_nan=False)
        object.__setattr__(self, "metadata", dict(self.metadata))

    @staticmethod
    def _encode_array(value: np.ndarray) -> dict[str, object]:
        data = np.asarray(value, dtype="<f8", order="C")
        raw = data.tobytes(order="C")
        return {
            "dtype": "<f8",
            "shape": list(data.shape),
            "sha256": hashlib.sha256(raw).hexdigest(),
            "base64": base64.b64encode(raw).decode("ascii"),
        }

    @staticmethod
    def _decode_array(name: str, payload: object) -> np.ndarray:
        if not isinstance(payload, Mapping) or set(payload) != {
            "dtype", "shape", "sha256", "base64"
        }:
            raise ValueError(f"checkpoint {name} array payload differs")
        if payload["dtype"] != "<f8":
            raise ValueError("checkpoint array dtype differs")
        shape = payload["shape"]
        if not isinstance(shape, list) or len(shape) != 2 or any(
            isinstance(item, bool) or not isinstance(item, int) or item < 1 for item in shape
        ):
            raise ValueError("checkpoint array shape differs")
        try:
            raw = base64.b64decode(str(payload["base64"]), validate=True)
        except ValueError as exc:
            raise ValueError("checkpoint array base64 is invalid") from exc
        digest = hashlib.sha256(raw).hexdigest()
        if digest != payload["sha256"]:
            raise ValueError("checkpoint array checksum differs")
        expected = int(np.prod(shape)) * 8
        if len(raw) != expected:
            raise ValueError("checkpoint array byte length differs")
        return np.frombuffer(raw, dtype="<f8").reshape(tuple(shape)).copy()

    def canonical_bytes(self) -> bytes:
        body = {
            "schema_version": 1,
            "time": self.time,
            "step_index": self.step_index,
            "transaction_serial": self.transaction_serial,
            "state": {
                "u": self._encode_array(self.state.u),
                "p": self._encode_array(self.state.p),
                "q": self._encode_array(self.state.q),
            },
            "metadata": self.metadata,
        }
        return json.dumps(
            body,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("utf-8") + b"\n"

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.canonical_bytes()).hexdigest()

    @classmethod
    def from_bytes(cls, value: bytes) -> "EvolutionCheckpoint":
        if not isinstance(value, bytes):
            raise TypeError("checkpoint input must be bytes")
        try:
            payload = json.loads(value.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise ValueError("checkpoint is not canonical JSON") from exc
        if not isinstance(payload, Mapping) or set(payload) != {
            "schema_version", "time", "step_index", "transaction_serial", "state", "metadata"
        }:
            raise ValueError("checkpoint top-level schema differs")
        if payload["schema_version"] != 1:
            raise ValueError("checkpoint schema version differs")
        state = payload["state"]
        if not isinstance(state, Mapping) or set(state) != {"u", "p", "q"}:
            raise ValueError("checkpoint state schema differs")
        checkpoint = cls(
            time=payload["time"],
            step_index=payload["step_index"],
            transaction_serial=payload["transaction_serial"],
            state=EvolutionState(
                cls._decode_array("u", state["u"]),
                cls._decode_array("p", state["p"]),
                cls._decode_array("q", state["q"]),
            ),
            metadata=payload["metadata"],
        )
        if checkpoint.canonical_bytes() != value:
            raise ValueError("checkpoint is not in canonical byte form")
        return checkpoint

    def write(self, path: Path) -> None:
        if not isinstance(path, Path):
            raise TypeError("checkpoint path must be pathlib.Path")
        path.write_bytes(self.canonical_bytes())

    @classmethod
    def read(cls, path: Path) -> "EvolutionCheckpoint":
        if not isinstance(path, Path):
            raise TypeError("checkpoint path must be pathlib.Path")
        return cls.from_bytes(path.read_bytes())


def array_content_sha256(*arrays: object) -> str:
    """Hash canonical shapes and little-endian contents for run receipts."""

    digest = hashlib.sha256()
    for index, value in enumerate(arrays):
        array = _finite_array(f"arrays[{index}]", value)
        normalized = np.asarray(array, dtype="<f8", order="C")
        digest.update(struct.pack("<I", normalized.ndim))
        digest.update(struct.pack(f"<{normalized.ndim}Q", *normalized.shape))
        digest.update(normalized.tobytes(order="C"))
    return digest.hexdigest()


__all__ = [
    "AcceptedStep",
    "COMPARATOR_METHOD",
    "EvolutionCheckpoint",
    "EvolutionRHS",
    "EvolutionState",
    "METHODS",
    "PRIMARY_METHOD",
    "ProposalGuard",
    "RightHandSide",
    "SBPFirstDerivative",
    "StageRecord",
    "StepProposal",
    "UniformRadialGrid",
    "accept_step",
    "array_content_sha256",
    "evolve_fixed_steps",
    "parity_project",
    "propose_step",
]
