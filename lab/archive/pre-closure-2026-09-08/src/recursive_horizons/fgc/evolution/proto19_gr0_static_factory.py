"""GR-0-only construction of the six persisted PROTO19 GR-0 shells.

The historical construction path lives in executable ``scripts`` modules.
That path was appropriate while the grids were being frozen, but importing it
at restart also imports authorization and study code that does not belong in a
post-persistence runtime.  This module reconstructs the same fixed objects
directly from source-library components and exactly four immutable inputs.

The returned members are fresh, time-zero :class:`Proto14RunMember` templates.
No persisted state is opened here.  HLT16 authenticates a descriptor and its
array payload before overlaying evolving state onto one of these templates.
Some imported analytic modules also contain dormant candidate-capable
definitions; this factory selects and evaluates only the frozen GR-0 path.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from fractions import Fraction
from hashlib import sha256
import json
from math import isfinite
import os
from pathlib import Path
import stat
import tomllib
from typing import Any, Mapping

import numpy as np

from .boundary_domain import BoundaryGeometry, CausalBudgetState
from .gr0_calibration import (
    construct_gr0_grid_initial_data,
    make_gr0_center_boundary_projector,
)
from .hlt16_member_codec import _digest, _template_identity
from .numerical_engine import EvolutionState, SBPFirstDerivative, array_content_sha256
from .proto11_runtime import project_gr0_reference_balanced_state
from .proto12_runtime import Proto12GR0EvolutionOperator
from .proto14_runtime import Proto14RunMember
from .proto5_runtime import (
    GR0RuntimeStageTransaction,
    GR0UniversalThresholds,
)
from .protocol_v17 import MEMBER_KEYS
from ..initial_data_preflight import PulseParameters


AMPLITUDE = "3"

STATIC_INPUT_PATHS = (
    "configs/fgc/fgc-1-cal9-run1.toml",
    "configs/fgc/fgc-1-rsp2-run1.toml",
    "results/fgc-1-hlt10-mon10.json",
    "results/fgc-1-rsp2-frz1.json",
)

_STATIC_INPUT_SHA256 = {
    "configs/fgc/fgc-1-cal9-run1.toml": (
        "8e6fafff638ffd49549efab8b6f969a01964c915116db9ce3698953201172176"
    ),
    "configs/fgc/fgc-1-rsp2-run1.toml": (
        "8213341f6cc5dd1d044a1ee9ae68b9ce7c7445ef9b933bcf8bcc944d30e4abd1"
    ),
    "results/fgc-1-hlt10-mon10.json": (
        "630d0d1843130180cb49dd237fe32e83e0336e7e44e7c684e4baf406b1883fe5"
    ),
    "results/fgc-1-rsp2-frz1.json": (
        "a6c147d7b369387e70318ee42eb56c97332702a227fd93d7ad20a1cfe2c1aa22"
    ),
}

# These are the HLT16 template identities referenced by the persisted
# generation-eight checkpoint.  They bind types and callables as well as
# scalar parameters, initial arrays, grids, projectors, and tracer identity.
_PERSISTED_TEMPLATE_SHA256 = {
    "RK4-2049": "c4edce6fbe30732b4b328f5a5dd1031f16a8ee2de1ce64bbe911a692d00e4213",
    "RK4-4097": "ee491381c58c76647f3bc68bbc75b5ced769ffded725a72fdaceeb46c41cf1ab",
    "RK4-8193": "1040629a09a92f36b8a056f3f748a6afe61be371b43926b164b0e5f3eda16e98",
    "SSPRK3-4097": "512ef13daf8bab8f4924af53a263abaeb67b6e166e97a31c757a26a6b0c2b9fb",
    "SSPRK3-8193": "7866cde0fcd0cf20b2867c9487c400b32c872d08e0a62e036f902d4b8a16d302",
    "SSPRK3-16385": "b46a78e67b7f5eba7210d3f53fbe0b22a9dad49b16aecf87d8ca995379f5b851",
}


class Proto19GR0StaticFactoryError(ValueError):
    """A static input or constructed shell differs from persisted GR-0."""


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[4]


def _relative_parts(relative: str) -> tuple[str, ...]:
    parts = tuple(relative.split("/"))
    if (
        not relative
        or relative.startswith("/")
        or not parts
        or any(part in {"", ".", ".."} for part in parts)
    ):
        raise Proto19GR0StaticFactoryError(
            f"static GR-0 input path escapes its root: {relative}"
        )
    return parts


def _read_regular_file_no_follow(root: Path, relative: str) -> bytes:
    """Read one root-relative regular file through no-follow descriptors."""

    parts = _relative_parts(relative)
    no_follow = getattr(os, "O_NOFOLLOW", None)
    directory_only = getattr(os, "O_DIRECTORY", None)
    if no_follow is None or directory_only is None:
        raise Proto19GR0StaticFactoryError(
            "static GR-0 input platform lacks no-follow directory opens"
        )
    close_on_exec = getattr(os, "O_CLOEXEC", 0)
    directory_flags = os.O_RDONLY | no_follow | directory_only | close_on_exec
    leaf_flags = os.O_RDONLY | no_follow | getattr(os, "O_NONBLOCK", 0) | close_on_exec
    root_path = os.path.abspath(os.fspath(root))
    directory_fd: int | None = None
    leaf_fd: int | None = None
    try:
        directory_fd = os.open(root_path, directory_flags)
        for component in parts[:-1]:
            next_fd = os.open(component, directory_flags, dir_fd=directory_fd)
            os.close(directory_fd)
            directory_fd = next_fd
        leaf_fd = os.open(parts[-1], leaf_flags, dir_fd=directory_fd)
        before = os.fstat(leaf_fd)
        if not stat.S_ISREG(before.st_mode) or before.st_nlink != 1:
            raise Proto19GR0StaticFactoryError(
                f"static GR-0 input is not a single-link regular file: {relative}"
            )
        chunks: list[bytes] = []
        while True:
            chunk = os.read(leaf_fd, 1024 * 1024)
            if not chunk:
                break
            chunks.append(chunk)
        payload = b"".join(chunks)
        after = os.fstat(leaf_fd)
        identity_before = (
            before.st_dev,
            before.st_ino,
            before.st_mode,
            before.st_nlink,
            before.st_size,
            before.st_mtime_ns,
            before.st_ctime_ns,
        )
        identity_after = (
            after.st_dev,
            after.st_ino,
            after.st_mode,
            after.st_nlink,
            after.st_size,
            after.st_mtime_ns,
            after.st_ctime_ns,
        )
        if (
            identity_after != identity_before
            or after.st_nlink != 1
            or len(payload) != after.st_size
        ):
            raise Proto19GR0StaticFactoryError(
                f"static GR-0 input changed while captured: {relative}"
            )
        return payload
    except Proto19GR0StaticFactoryError:
        raise
    except OSError as error:
        raise Proto19GR0StaticFactoryError(
            f"static GR-0 input is unavailable or unsafe: {relative}"
        ) from error
    finally:
        if leaf_fd is not None:
            os.close(leaf_fd)
        if directory_fd is not None:
            os.close(directory_fd)


def _validated_static_input_bytes(
    payloads: Mapping[str, bytes],
) -> dict[str, bytes]:
    try:
        keys = tuple(payloads)
    except TypeError as error:
        raise Proto19GR0StaticFactoryError(
            "captured static GR-0 input inventory differs"
        ) from error
    if len(keys) != len(STATIC_INPUT_PATHS) or set(keys) != set(STATIC_INPUT_PATHS):
        raise Proto19GR0StaticFactoryError(
            "captured static GR-0 input inventory differs"
        )
    answer: dict[str, bytes] = {}
    for relative in STATIC_INPUT_PATHS:
        try:
            raw = payloads[relative]
        except (KeyError, TypeError) as error:
            raise Proto19GR0StaticFactoryError(
                "captured static GR-0 input inventory differs"
            ) from error
        if not isinstance(raw, bytes):
            raise Proto19GR0StaticFactoryError(
                f"captured static GR-0 input is not immutable bytes: {relative}"
            )
        if sha256(raw).hexdigest() != _STATIC_INPUT_SHA256[relative]:
            raise Proto19GR0StaticFactoryError(
                f"static GR-0 input identity differs: {relative}"
            )
        answer[relative] = raw
    return answer


def _read_static_inputs(root: Path) -> dict[str, bytes]:
    return _validated_static_input_bytes(
        {
            relative: _read_regular_file_no_follow(root, relative)
            for relative in STATIC_INPUT_PATHS
        }
    )


def _decode_inputs(
    root: Path | None,
    static_input_bytes: Mapping[str, bytes] | None,
) -> tuple[dict[str, Any], ...]:
    if static_input_bytes is None:
        if root is None:
            raise Proto19GR0StaticFactoryError("static GR-0 input root is absent")
        raw = _read_static_inputs(root)
    else:
        if not isinstance(static_input_bytes, Mapping):
            raise Proto19GR0StaticFactoryError(
                "captured static GR-0 inputs are not a mapping"
            )
        raw = _validated_static_input_bytes(static_input_bytes)
    try:
        cal9 = tomllib.loads(raw[STATIC_INPUT_PATHS[0]].decode("utf-8"))
        rsp2 = tomllib.loads(raw[STATIC_INPUT_PATHS[1]].decode("utf-8"))
        hlt10 = json.loads(raw[STATIC_INPUT_PATHS[2]])
        rsp2_freeze = json.loads(raw[STATIC_INPUT_PATHS[3]])
    except (UnicodeDecodeError, tomllib.TOMLDecodeError, json.JSONDecodeError) as error:
        raise Proto19GR0StaticFactoryError(
            "static GR-0 input cannot be decoded"
        ) from error
    values = (cal9, rsp2, hlt10, rsp2_freeze)
    if any(not isinstance(value, dict) for value in values):
        raise Proto19GR0StaticFactoryError("static GR-0 input root is not a table")
    return values


def _fraction(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (str, int)):
        raise Proto19GR0StaticFactoryError(f"{label} is not frozen rational text")
    try:
        answer = float(Fraction(value))
    except (ValueError, ZeroDivisionError) as error:
        raise Proto19GR0StaticFactoryError(
            f"{label} is not frozen rational text"
        ) from error
    if not isfinite(answer):
        raise Proto19GR0StaticFactoryError(f"{label} is not finite")
    return answer


def _mapping(value: object, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise Proto19GR0StaticFactoryError(f"{label} is not a mapping")
    return value


def _interpolate_columns(
    coordinates: np.ndarray,
    values: np.ndarray,
    positions: np.ndarray,
) -> np.ndarray:
    if (
        coordinates.ndim != 1
        or values.ndim != 2
        or values.shape[0] != coordinates.size
        or positions.ndim != 1
    ):
        raise ValueError("tracer interpolation shapes differ")
    return np.column_stack(
        [
            np.interp(positions, coordinates, values[:, field])
            for field in range(values.shape[1])
        ]
    )


@dataclass
class NormalFlowTracers:
    """Local byte-compatible form of the historical GR-0 tracer ledger."""

    labels: np.ndarray
    positions: np.ndarray
    proper_times: np.ndarray
    event_proper_times: list[np.ndarray]
    event_fields: list[np.ndarray]
    cutoff: float
    outer_radius: float

    @classmethod
    def create(
        cls,
        *,
        minimum: float,
        maximum: float,
        spacing: float,
        state: EvolutionState,
        coordinates: np.ndarray,
        cutoff: float,
        outer_radius: float,
    ) -> "NormalFlowTracers":
        count_float = (maximum - minimum) / spacing
        count = int(round(count_float))
        if abs(count_float - count) > 1.0e-12:
            raise ValueError("tracer interval is not aligned")
        labels = np.linspace(minimum, maximum, count + 1, dtype=np.float64)
        answer = cls(
            labels=labels,
            positions=labels.copy(),
            proper_times=np.zeros(labels.size, dtype=np.float64),
            event_proper_times=[],
            event_fields=[],
            cutoff=cutoff,
            outer_radius=outer_radius,
        )
        answer.append_common_event(state, coordinates)
        return answer

    def _metric_samples(
        self,
        state: EvolutionState,
        coordinates: np.ndarray,
        positions: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        fields = _interpolate_columns(coordinates, state.u, positions)
        lapse = fields[:, 0]
        shift = fields[:, 1]
        if np.any(lapse <= 0.0):
            raise ValueError("normal-flow tracer encountered nonpositive lapse")
        return lapse, shift

    def preview_advance(
        self,
        *,
        old_state: EvolutionState,
        new_state: EvolutionState,
        coordinates: np.ndarray,
        step_size: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        old_lapse, old_shift = self._metric_samples(
            old_state, coordinates, self.positions
        )
        old_velocity = -old_shift
        predictor = self.positions + step_size * old_velocity
        if np.any(predictor <= 0.0) or np.any(predictor >= self.outer_radius):
            raise ValueError("normal-flow tracer predictor left the numerical domain")
        new_lapse, new_shift = self._metric_samples(new_state, coordinates, predictor)
        new_velocity = -new_shift
        positions = self.positions + 0.5 * step_size * (old_velocity + new_velocity)
        proper = self.proper_times + 0.5 * step_size * (old_lapse + new_lapse)
        if (
            np.any(positions <= 0.0)
            or np.any(positions >= self.outer_radius)
            or np.any(proper <= self.proper_times)
        ):
            raise ValueError("normal-flow tracer update lost its timelike domain")
        return positions, proper

    def commit_advance(
        self,
        positions: np.ndarray,
        proper_times: np.ndarray,
    ) -> None:
        if (
            positions.shape != self.positions.shape
            or proper_times.shape != self.proper_times.shape
        ):
            raise ValueError("normal-flow tracer commit shapes differ")
        self.positions = np.asarray(positions, dtype=np.float64).copy()
        self.proper_times = np.asarray(proper_times, dtype=np.float64).copy()

    def field_sample(
        self,
        state: EvolutionState,
        coordinates: np.ndarray,
    ) -> np.ndarray:
        fields = _interpolate_columns(coordinates, state.u, self.positions)
        radius = self.positions
        answer = np.column_stack(
            (
                fields[:, 0] - 1.0,
                fields[:, 1] / radius,
                fields[:, 2] - 1.0,
                fields[:, 3] / radius - 1.0,
                fields[:, 4] / self.cutoff,
                fields[:, 5] / self.cutoff,
            )
        )
        if not np.all(np.isfinite(answer)):
            raise ValueError("normal-flow tracer field sample became nonfinite")
        return answer

    def append_common_event(
        self,
        state: EvolutionState,
        coordinates: np.ndarray,
    ) -> None:
        self.event_proper_times.append(self.proper_times.copy())
        self.event_fields.append(self.field_sample(state, coordinates))


# HLT16 descriptors intentionally encode the historical type identity.  The
# compatibility name is metadata only; no module with this name is imported.
NormalFlowTracers.__module__ = "scripts.run_fgc_gr0_calibration"


def _pulse_parameters(physical: Mapping[str, Any]) -> PulseParameters:
    return PulseParameters(
        chi_amplitude=_fraction(AMPLITUDE, "chi amplitude"),
        center=_fraction(physical.get("chi_center"), "chi center"),
        half_width=_fraction(physical.get("chi_half_width"), "chi half width"),
        phi_amplitude=_fraction(
            physical.get("phi_seed_amplitude"), "phi seed amplitude"
        ),
        planck_mass=_fraction(physical.get("planck_mass"), "Planck mass"),
        scalar_mass=_fraction(physical.get("scalar_mass"), "scalar mass"),
        quartic_coupling=_fraction(
            physical.get("quartic_coupling"), "quartic coupling"
        ),
    )


def _build_member(
    *,
    physical: Mapping[str, Any],
    method: Mapping[str, Any],
    numerics: Mapping[str, Any],
    thresholds: Mapping[str, Any],
    frozen_input: Mapping[str, Any],
    point_count: int,
) -> Proto14RunMember:
    method_label = method.get("method_label")
    spatial_order = method.get("spatial_order")
    if (
        method_label not in {"RK4", "SSPRK3"}
        or isinstance(spatial_order, bool)
        or spatial_order not in {2, 4}
        or frozen_input.get("amplitude") != AMPLITUDE
        or frozen_input.get("method") != method_label
        or frozen_input.get("point_count") != point_count
        or frozen_input.get("spatial_order") != spatial_order
        or frozen_input.get("integrator_id") != method.get("integrator_id")
    ):
        raise Proto19GR0StaticFactoryError("frozen GR-0 member identity differs")

    initial = construct_gr0_grid_initial_data(
        _pulse_parameters(physical),
        point_count=point_count,
        outer_radius=_fraction(physical.get("outer_radius"), "outer radius"),
        constraint_method=str(method.get("constraint_solve_method")),
        diagnostic_spatial_order=int(spatial_order),
    )
    state = project_gr0_reference_balanced_state(
        initial, spatial_order=int(spatial_order)
    )
    state_sha256 = array_content_sha256(state.u, state.p, state.q)
    if state_sha256 != frozen_input.get("projected_state_sha256"):
        raise Proto19GR0StaticFactoryError(
            f"constructed GR-0 state differs: {method_label}-{point_count}"
        )
    projected_initial = replace(initial, state=state)

    operator = Proto12GR0EvolutionOperator(
        projected_initial.grid,
        spatial_order=int(spatial_order),
        ko_dissipation=_fraction(numerics.get("ko_dissipation"), "KO dissipation"),
        raw_tolerance=_fraction(
            thresholds.get("source_residual_infinity_max"), "source residual limit"
        ),
        kinetic_condition_maximum=_fraction(
            thresholds.get("kinetic_condition_number_max"),
            "kinetic condition limit",
        ),
        maximum_refinement_iterations=int(thresholds.get("source_iteration_max", -1)),
        point_batch_size=int(numerics.get("source_point_batch_size", -1)),
    )
    before = array_content_sha256(state.u, state.p, state.q)
    initial_rhs = operator(0.0, state)
    diagnostics = initial_rhs.diagnostics
    if (
        array_content_sha256(state.u, state.p, state.q) != before
        or diagnostics.get("source_raw_gate_passed") is not True
        or diagnostics.get("SRC4_tensor_contracted_reference_source") is not True
        or diagnostics.get("PROTO11_interior_q_reprojected") is not False
        or diagnostics.get("PROTO11_exact_reference_equilibrium_applied") is not False
    ):
        raise Proto19GR0StaticFactoryError(
            f"constructed GR-0 source preflight differs: {method_label}-{point_count}"
        )

    derivative = SBPFirstDerivative(projected_initial.grid, int(spatial_order))
    transaction = GR0RuntimeStageTransaction(
        thresholds=GR0UniversalThresholds(
            source_residual_maximum=_fraction(
                thresholds.get("source_residual_infinity_max"),
                "source residual limit",
            ),
            source_iteration_maximum=int(thresholds.get("source_iteration_max", -1)),
            kinetic_condition_maximum=_fraction(
                thresholds.get("kinetic_condition_number_max"),
                "kinetic condition limit",
            ),
        ),
        causal_state=CausalBudgetState(
            previous_speed_upper=float(diagnostics["coordinate_speed_upper"])
        ),
        boundary_geometry=BoundaryGeometry(
            _fraction(physical.get("outer_radius"), "outer radius"),
            _fraction(
                physical.get("measurement_radius_maximum"),
                "measurement radius",
            ),
            _fraction(
                thresholds.get("minimum_boundary_causal_buffer"),
                "minimum causal buffer",
            ),
            derivative.stencil_reach_intervals * projected_initial.grid.spacing,
        ),
        grid_spacing=projected_initial.grid.spacing,
        cfl_maximum=_fraction(numerics.get("cfl_maximum"), "CFL maximum"),
        hat_normal_factor=_fraction(
            thresholds.get("hat_normal_factor"), "hat normal factor"
        ),
    )
    tracers = NormalFlowTracers.create(
        minimum=_fraction(
            numerics.get("normal_flow_tracer_radius_minimum"),
            "tracer minimum",
        ),
        maximum=_fraction(
            numerics.get("normal_flow_tracer_radius_maximum"),
            "tracer maximum",
        ),
        spacing=_fraction(numerics.get("normal_flow_tracer_spacing"), "tracer spacing"),
        state=state,
        coordinates=projected_initial.grid.coordinates,
        cutoff=_fraction(physical.get("cutoff_Lambda"), "cutoff"),
        outer_radius=_fraction(physical.get("outer_radius"), "outer radius"),
    )
    member = Proto14RunMember(
        amplitude=AMPLITUDE,
        method_label=str(method_label),
        integrator_id=str(method.get("integrator_id")),
        spatial_order=int(spatial_order),
        point_count=point_count,
        input_hash=str(frozen_input.get("expanded_run_config_sha256")),
        initial=projected_initial,
        state=state,
        operator=operator,
        projector=make_gr0_center_boundary_projector(
            projected_initial,
            fixed_outer_rows=int(numerics.get("fixed_outer_rows", -1)),
        ),
        transaction=transaction,
        tracers=tracers,
    )
    if member.key != f"{method_label}-{point_count}":
        raise Proto19GR0StaticFactoryError("constructed GR-0 member key differs")
    return member


def _validated_sources(
    cal9: Mapping[str, Any],
    rsp2: Mapping[str, Any],
    hlt10: Mapping[str, Any],
    rsp2_freeze: Mapping[str, Any],
) -> tuple[
    Mapping[str, Any],
    Mapping[tuple[str, int], Mapping[str, Any]],
    Mapping[str, Any],
]:
    cal9_scope = _mapping(cal9.get("scope"), "CAL9 scope")
    cal9_numerics = _mapping(cal9.get("numerics"), "CAL9 numerics")
    hlt10_payload = _mapping(hlt10.get("artifact_payload"), "HLT10 payload")
    records = hlt10_payload.get("frozen_run_inputs")
    if (
        cal9.get("artifact_id") != "FGC-1-CAL9-RUN1-PLAN"
        or cal9_scope.get("branch") != "GR-0"
        or cal9_scope.get("target_protocol") != "FGC-2-SF1-PROTO12"
        or cal9_numerics.get("resolutions") != [2049, 4097, 8193]
        or hlt10.get("artifact_id") != "FGC-1-HLT10-MON10"
        or not isinstance(records, list)
        or len(records) != 12
    ):
        raise Proto19GR0StaticFactoryError("CAL9/HLT10 static scope differs")
    selected: dict[tuple[str, int], Mapping[str, Any]] = {}
    for item in records:
        record = _mapping(item, "HLT10 frozen input")
        if record.get("amplitude") == AMPLITUDE:
            key = (str(record.get("method")), int(record.get("point_count", -1)))
            if key in selected:
                raise Proto19GR0StaticFactoryError("HLT10 frozen input is duplicated")
            selected[key] = record
    if set(selected) != {
        ("RK4", 2049),
        ("RK4", 4097),
        ("RK4", 8193),
        ("SSPRK3", 2049),
        ("SSPRK3", 4097),
        ("SSPRK3", 8193),
    }:
        raise Proto19GR0StaticFactoryError("HLT10 amplitude-three ladder differs")

    rsp2_scope = _mapping(rsp2.get("scope"), "RSP2 scope")
    rsp2_numerics = _mapping(rsp2.get("numerics"), "RSP2 numerics")
    rsp2_payload = _mapping(rsp2_freeze.get("artifact_payload"), "RSP2 freeze payload")
    new_input = _mapping(rsp2_payload.get("new_run_input"), "RSP2 new input")
    if (
        rsp2.get("artifact_id") != "FGC-1-RSP2-RUN1-PLAN"
        or rsp2_scope.get("branch") != "GR-0"
        or rsp2_scope.get("amplitude") != AMPLITUDE
        or rsp2_scope.get("method") != "SSPRK3"
        or rsp2_numerics.get("new_resolution") != 16385
        or rsp2_freeze.get("artifact_id") != "FGC-1-RSP2-FRZ1"
        or new_input.get("amplitude") != AMPLITUDE
        or new_input.get("method") != "SSPRK3"
        or new_input.get("point_count") != 16385
    ):
        raise Proto19GR0StaticFactoryError("RSP2 static scope differs")
    return cal9, selected, new_input


def build_static_gr0_shells(
    repository_root: Path | None = None,
    *,
    static_input_bytes: Mapping[str, bytes] | None = None,
) -> Mapping[str, Proto14RunMember]:
    """Return six fresh GR-0 templates in the PROTO17 canonical order.

    Without ``static_input_bytes``, only :data:`STATIC_INPUT_PATHS` are opened
    through the no-follow regular-file reader.  A supplied capture must contain
    exactly those four immutable byte strings and completely replaces file
    access.  No historical checkpoint, raw trajectory, candidate branch,
    executable script, or persisted HLT16 state is imported or read.
    """

    root: Path | None = None
    if static_input_bytes is None:
        root = (
            _repository_root()
            if repository_root is None
            else Path(os.path.abspath(os.fspath(repository_root)))
        )
    cal9, rsp2, hlt10, rsp2_freeze = _decode_inputs(root, static_input_bytes)
    cal9, cal9_inputs, rsp2_input = _validated_sources(cal9, rsp2, hlt10, rsp2_freeze)
    cal9_physical = _mapping(cal9.get("physical_inputs"), "CAL9 physical inputs")
    cal9_numerics = _mapping(cal9.get("numerics"), "CAL9 numerics")
    cal9_thresholds = _mapping(
        cal9.get("universal_thresholds"), "CAL9 universal thresholds"
    )
    methods = {
        "RK4": _mapping(cal9.get("primary_method"), "CAL9 primary method"),
        "SSPRK3": _mapping(cal9.get("comparator_method"), "CAL9 comparator method"),
    }

    rsp2_physical = _mapping(rsp2.get("physical_inputs"), "RSP2 physical inputs")
    rsp2_method = _mapping(rsp2.get("method"), "RSP2 method")
    rsp2_numerics = _mapping(rsp2.get("numerics"), "RSP2 numerics")
    rsp2_thresholds = _mapping(
        rsp2.get("universal_thresholds"), "RSP2 universal thresholds"
    )

    answer: dict[str, Proto14RunMember] = {}
    for key in MEMBER_KEYS:
        method_label, point_text = key.split("-", 1)
        point_count = int(point_text)
        if key == "SSPRK3-16385":
            member = _build_member(
                physical=rsp2_physical,
                method=rsp2_method,
                numerics=rsp2_numerics,
                thresholds=rsp2_thresholds,
                frozen_input=rsp2_input,
                point_count=point_count,
            )
        else:
            member = _build_member(
                physical=cal9_physical,
                method=methods[method_label],
                numerics=cal9_numerics,
                thresholds=cal9_thresholds,
                frozen_input=cal9_inputs[(method_label, point_count)],
                point_count=point_count,
            )
        template_sha256 = _digest(_template_identity(member))
        if template_sha256 != _PERSISTED_TEMPLATE_SHA256[key]:
            raise Proto19GR0StaticFactoryError(
                f"constructed GR-0 template differs from persisted HLT16: {key}"
            )
        answer[key] = member
    if tuple(answer) != MEMBER_KEYS:
        raise Proto19GR0StaticFactoryError("constructed GR-0 shell order differs")
    return answer


__all__ = [
    "AMPLITUDE",
    "NormalFlowTracers",
    "Proto19GR0StaticFactoryError",
    "STATIC_INPUT_PATHS",
    "build_static_gr0_shells",
]
