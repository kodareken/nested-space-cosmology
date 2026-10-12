"""Lossless GR-0 member codec for the authenticated PROTO18 runtime.

The descriptor contains the predecessor cursor.  The successor cursor points
at this descriptor, so embedding the successor would create a circular hash.
GEN0 uses the sealed PROTO17 cursor; evolved descriptors use their exact
PROTO18 predecessor.  This module performs no I/O or trajectory work.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, is_dataclass
from hashlib import sha256
import json
import math
from numbers import Integral, Real
from typing import Any, Mapping

import numpy as np

from .boundary_domain import CausalBudgetState
from .numerical_engine import EvolutionState, array_content_sha256
from .proto14_runtime import (
    Proto14RunMember,
    proto14_checkpoint_extension,
    restore_proto14_checkpoint_extension,
)
from . import proto15_runtime as p15
from .proto17_pure_construction import (
    MEMBER_KEYS,
    Proto17ConstructionError,
    cursor as proto17_cursor,
)
from .proto5_runtime import GR0RuntimeMonitorState


RUNTIME_PROTOCOL = "FGC-2-SF1-PROTO18"
GEN0_PROTOCOL = "FGC-2-SF1-PROTO17"
ARRAY_NAMES = (
    "u",
    "p",
    "q",
    "grid_coordinates",
    "tracer_labels",
    "tracer_positions",
    "tracer_proper_times",
    "event_proper_times",
    "event_fields",
)


class HLT16MemberCodecError(ValueError):
    """A member snapshot differs from the frozen restart contract."""


@dataclass(frozen=True, slots=True)
class HLT16MemberSnapshot:
    arrays: Mapping[str, np.ndarray]
    metadata: Mapping[str, Any]


def _canonical(value: object) -> bytes:
    try:
        return json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("ascii")
    except (TypeError, ValueError) as error:
        raise HLT16MemberCodecError("metadata is not canonical JSON") from error


def _digest(value: object) -> str:
    return sha256(_canonical(value)).hexdigest()


def _integer(value: object, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
        raise HLT16MemberCodecError(f"{label} is not a nonnegative integer")
    return int(value)


def _finite_hex(value: object, label: str, *, positive: bool = False) -> str:
    if not isinstance(value, str):
        raise HLT16MemberCodecError(f"{label} is not canonical binary64 text")
    try:
        result = float.fromhex(value)
    except ValueError as error:
        raise HLT16MemberCodecError(f"{label} is not canonical binary64 text") from error
    if not math.isfinite(result) or result.hex() != value or (positive and result <= 0.0):
        raise HLT16MemberCodecError(f"{label} is not canonical binary64 text")
    return value


def _strict_array(value: object, name: str) -> np.ndarray:
    if not isinstance(value, np.ndarray):
        raise HLT16MemberCodecError(f"{name} is not a NumPy array")
    if (
        value.dtype != np.dtype("<f8")
        or not value.flags.c_contiguous
        or not value.shape
        or not np.isfinite(value).all()
    ):
        raise HLT16MemberCodecError(
            f"{name} is not finite little-endian C-order float64"
        )
    return value.copy(order="C")


def _normalized_json(value: object) -> object:
    """Encode static shell data without platform-dependent float text."""
    if is_dataclass(value):
        return _normalized_json(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _normalized_json(item) for key, item in sorted(value.items())}
    if isinstance(value, (tuple, list)):
        return [_normalized_json(item) for item in value]
    if isinstance(value, bool) or value is None or isinstance(value, str):
        return value
    if isinstance(value, Integral):
        return int(value)
    if isinstance(value, Real):
        result = float(value)
        if not math.isfinite(result):
            raise HLT16MemberCodecError("static identity contains a nonfinite value")
        return {"binary64_hex": result.hex()}
    raise HLT16MemberCodecError(
        f"unsupported static identity type: {type(value).__name__}"
    )


def _type_identity(value: object) -> str:
    if callable(value) and hasattr(value, "__module__") and hasattr(value, "__qualname__"):
        return f"{value.__module__}.{value.__qualname__}"
    cls = type(value)
    return f"{cls.__module__}.{cls.__qualname__}"


def _operator_identity(operator: object) -> Mapping[str, object]:
    fields: dict[str, object] = {"type": _type_identity(operator)}
    for name in (
        "ko_dissipation",
        "residual_tolerance",
        "raw_tolerance",
        "kinetic_condition_maximum",
    ):
        if hasattr(operator, name):
            fields[name] = _normalized_json(getattr(operator, name))
    derivative = getattr(operator, "derivative", None)
    if derivative is not None and hasattr(derivative, "order"):
        fields["derivative_order"] = _integer(derivative.order, "derivative order")
    return fields


def _template_identity(member: Proto14RunMember) -> dict[str, object]:
    initial = member.initial
    grid = initial.grid
    transaction = member.transaction
    return {
        "amplitude": member.amplitude,
        "method": member.method_label,
        "point_count": member.point_count,
        "integrator_id": member.integrator_id,
        "spatial_order": member.spatial_order,
        "input_hash": member.input_hash,
        "member_type": _type_identity(member),
        "initial_type": _type_identity(initial),
        "operator": _operator_identity(member.operator),
        "projector": {"callable": _type_identity(member.projector)},
        "transaction_type": _type_identity(transaction),
        "tracer_type": _type_identity(member.tracers),
        "grid": {
            "minimum_hex": float(grid.minimum).hex(),
            "maximum_hex": float(grid.maximum).hex(),
            "spacing_hex": float(grid.spacing).hex(),
            "point_count": grid.point_count,
            "coordinates_sha256": array_content_sha256(grid.coordinates),
        },
        "initial": {
            "constraint_method": initial.constraint_method,
            "parameters": _normalized_json(initial.parameters),
            "support_minimum_index": initial.support_minimum_index,
            "support_maximum_index": initial.support_maximum_index,
            "peak_compactness": _normalized_json(initial.peak_compactness),
            "peak_radius": _normalized_json(initial.peak_radius),
            "outer_mass": _normalized_json(initial.outer_mass),
            "vacuum_momentum_constant": _normalized_json(
                initial.vacuum_momentum_constant
            ),
            "reduction_constraint_infinity": _normalized_json(
                initial.reduction_constraint_infinity
            ),
            "no_initial_trapped_sphere": initial.no_initial_trapped_sphere,
            "state_sha256": array_content_sha256(
                initial.state.u, initial.state.p, initial.state.q
            ),
        },
        "transaction": {
            "thresholds": _normalized_json(transaction.thresholds),
            "boundary_geometry": _normalized_json(transaction.boundary_geometry),
            "grid_spacing": _normalized_json(transaction.grid_spacing),
            "cfl_maximum": _normalized_json(transaction.cfl_maximum),
            "hat_normal_factor": _normalized_json(transaction.hat_normal_factor),
        },
    }


def _validated_proto15_cursor(payload: object) -> p15.Proto15Cursor:
    try:
        # `_cursor` is a constructor and deliberately recomputes the chain
        # hash.  A persisted descriptor must instead verify the hash it
        # observed, otherwise a one-bit provenance mutation would be silently
        # normalized into a different cursor.
        result = p15.Proto15Cursor(dict(payload))  # type: ignore[arg-type]
        result.validate()
        return result
    except (TypeError, ValueError, p15.Proto15CursorContractError) as error:
        raise HLT16MemberCodecError("PROTO18 predecessor cursor differs") from error


def _validated_proto17_cursor(payload: object) -> Mapping[str, object]:
    try:
        return dict(proto17_cursor(payload))  # type: ignore[arg-type]
    except (TypeError, ValueError, Proto17ConstructionError) as error:
        raise HLT16MemberCodecError("PROTO17 GEN0 cursor differs") from error


def _validated_monitor(value: object) -> dict[str, object]:
    fields = {
        "accepted_stage_count",
        "last_transaction_serial",
        "last_accepted_time",
        "first_failed_premise",
        "first_failed_transaction_serial",
        "first_failed_time",
    }
    if not isinstance(value, Mapping) or set(value) != fields:
        raise HLT16MemberCodecError("monitor schema differs")
    result = dict(value)
    result["accepted_stage_count"] = _integer(
        value["accepted_stage_count"], "accepted stage count"
    )
    serial = value["last_transaction_serial"]
    if isinstance(serial, bool) or not isinstance(serial, Integral) or int(serial) < -1:
        raise HLT16MemberCodecError("monitor serial differs")
    result["last_transaction_serial"] = int(serial)
    for name in ("last_accepted_time", "first_failed_time"):
        item = value[name]
        if item is not None and (
            not isinstance(item, Real) or not math.isfinite(float(item))
        ):
            raise HLT16MemberCodecError(f"monitor {name} differs")
    failed_serial = value["first_failed_transaction_serial"]
    if failed_serial is not None:
        result["first_failed_transaction_serial"] = _integer(
            failed_serial, "first failed transaction serial"
        )
    failed = value["first_failed_premise"]
    if failed is not None and not isinstance(failed, str):
        raise HLT16MemberCodecError("monitor failed premise differs")
    if (failed is None) != (failed_serial is None) or (failed is None) != (
        value["first_failed_time"] is None
    ):
        raise HLT16MemberCodecError("monitor failure tuple differs")
    return result


def _validated_causal(value: object) -> dict[str, object]:
    if not isinstance(value, Mapping):
        raise HLT16MemberCodecError("causal state schema differs")
    try:
        state = CausalBudgetState(**dict(value))
    except (TypeError, ValueError) as error:
        raise HLT16MemberCodecError("causal state differs") from error
    return asdict(state)


def _validated_tdg6(value: object):
    if not isinstance(value, Mapping) or set(value) != {
        "extension",
        "accumulated_debit_hex",
    }:
        raise HLT16MemberCodecError("TDG6 schema differs")
    debit_hex = value["accumulated_debit_hex"]
    if not isinstance(debit_hex, list):
        raise HLT16MemberCodecError("TDG6 debit schema differs")
    try:
        debit = np.asarray(
            [float.fromhex(_finite_hex(item, "TDG6 debit")) for item in debit_hex],
            dtype=np.float64,
        )
        from .tdg6_temporal_admission_runtime import restore_tdg6_checkpoint_extension

        ledger = restore_tdg6_checkpoint_extension(value["extension"], debit)
    except (TypeError, ValueError) as error:
        raise HLT16MemberCodecError("TDG6 extension differs") from error
    return ledger, {
        "extension": dict(value["extension"]),
        "accumulated_debit_hex": [float(item).hex() for item in debit],
    }


def validate_codec_metadata(metadata: Mapping[str, Any]) -> dict[str, Any]:
    """Validate one complete, non-circular descriptor payload."""
    required = {
        "runtime_identity",
        "template",
        "template_sha256",
        "accepted_plan",
        "counters",
        "tdg6",
        "monitor",
        "causal",
        "tracer",
        "provenance_cursor",
    }
    if not isinstance(metadata, Mapping) or set(metadata) != required:
        raise HLT16MemberCodecError("codec metadata schema differs")
    runtime = metadata["runtime_identity"]
    if not isinstance(runtime, Mapping) or set(runtime) != {
        "protocol_artifact_id",
        "campaign_id",
        "branch",
        "amplitude",
        "member_key",
        "method",
        "point_count",
        "accepted_time",
        "generation",
    }:
        raise HLT16MemberCodecError("runtime identity schema differs")
    key = runtime["member_key"]
    generation = _integer(runtime["generation"], "generation")
    if (
        runtime["protocol_artifact_id"] != RUNTIME_PROTOCOL
        or runtime["branch"] != "GR-0"
        or runtime["amplitude"] != "3"
        or key not in MEMBER_KEYS
        or runtime["method"] != str(key).split("-", 1)[0]
        or runtime["point_count"] != int(str(key).split("-", 1)[1])
        or not isinstance(runtime["campaign_id"], str)
        or not runtime["campaign_id"]
    ):
        raise HLT16MemberCodecError("runtime identity differs")
    try:
        accepted_time = p15._time_identity(runtime["accepted_time"])
    except (TypeError, ValueError, p15.Proto15CursorContractError) as error:
        raise HLT16MemberCodecError("accepted time differs") from error

    template = metadata["template"]
    if not isinstance(template, Mapping) or _digest(template) != metadata["template_sha256"]:
        raise HLT16MemberCodecError("template identity differs")
    _canonical(template)
    if (
        template.get("amplitude") != "3"
        or template.get("method") != runtime["method"]
        or template.get("point_count") != runtime["point_count"]
    ):
        raise HLT16MemberCodecError("template/runtime identity differs")

    counters = metadata["counters"]
    if not isinstance(counters, Mapping) or set(counters) != {
        "step_index",
        "transaction_serial",
        "accepted_macro_steps",
        "source_retry_count",
        "CFL_retry_count",
    }:
        raise HLT16MemberCodecError("counter schema differs")
    checked_counters = {
        name: _integer(value, name) for name, value in counters.items()
    }
    ledger, tdg6 = _validated_tdg6(metadata["tdg6"])
    if (
        checked_counters["accepted_macro_steps"] != ledger.accepted_macro_step_count
        or ledger.last_accepted_time.hex() != accepted_time["binary64_hex"]
    ):
        raise HLT16MemberCodecError("TDG6 counter/time relation differs")

    monitor = _validated_monitor(metadata["monitor"])
    causal = _validated_causal(metadata["causal"])
    if (
        float(causal["accepted_time"]).hex() != accepted_time["binary64_hex"]
        or monitor["last_transaction_serial"] + 1
        != checked_counters["transaction_serial"]
    ):
        raise HLT16MemberCodecError("monitor/causal runtime relation differs")

    tracer = metadata["tracer"]
    if not isinstance(tracer, Mapping) or set(tracer) != {
        "label_count",
        "event_count",
        "cutoff_hex",
        "outer_radius_hex",
    }:
        raise HLT16MemberCodecError("tracer metadata schema differs")
    checked_tracer = {
        "label_count": _integer(tracer["label_count"], "tracer label count"),
        "event_count": _integer(tracer["event_count"], "tracer event count"),
        "cutoff_hex": _finite_hex(tracer["cutoff_hex"], "tracer cutoff", positive=True),
        "outer_radius_hex": _finite_hex(
            tracer["outer_radius_hex"], "tracer outer radius", positive=True
        ),
    }
    if checked_tracer["label_count"] < 1 or checked_tracer["event_count"] < 1:
        raise HLT16MemberCodecError("tracer counts differ")

    provenance = metadata["provenance_cursor"]
    if not isinstance(provenance, Mapping) or set(provenance) != {
        "source_protocol",
        "bridge",
        "payload",
    }:
        raise HLT16MemberCodecError("provenance cursor schema differs")
    source = provenance["source_protocol"]
    accepted_plan = metadata["accepted_plan"]
    if source == GEN0_PROTOCOL:
        bridge = provenance["bridge"]
        if (
            generation != 0
            or bridge not in {
                "authenticated_GEN0_predecessor",
                "authenticated_GEN0_campaign_rebase",
            }
            or accepted_plan is not None
        ):
            raise HLT16MemberCodecError("GEN0 provenance scope differs")
        checked_payload = _validated_proto17_cursor(provenance["payload"])
        if checked_payload["accepted_boundary_time"] != accepted_time:
            raise HLT16MemberCodecError("GEN0 cursor/runtime time differs")
        same_campaign = checked_payload["campaign_id"] == runtime["campaign_id"]
        if same_campaign != (bridge == "authenticated_GEN0_predecessor"):
            raise HLT16MemberCodecError("GEN0 campaign bridge relation differs")
    elif source == RUNTIME_PROTOCOL:
        if (
            generation < 1
            or provenance["bridge"] != "active_runtime_predecessor"
            or accepted_plan is None
        ):
            raise HLT16MemberCodecError("evolved provenance scope differs")
        cursor_object = _validated_proto15_cursor(provenance["payload"])
        try:
            accepted_plan = p15._validated_accepted_plan(
                cursor_object, accepted_plan, accepted_time
            )
        except (TypeError, ValueError, p15.Proto15CursorContractError) as error:
            raise HLT16MemberCodecError("accepted plan/provenance relation differs") from error
        checked_payload = dict(cursor_object.payload)
    else:
        raise HLT16MemberCodecError("provenance protocol is forbidden")
    if checked_payload["member_key"] != key:
        raise HLT16MemberCodecError("provenance/runtime member relation differs")
    if source == RUNTIME_PROTOCOL and checked_payload["campaign_id"] != runtime["campaign_id"]:
        raise HLT16MemberCodecError("provenance/runtime campaign relation differs")

    result = {
        "runtime_identity": dict(runtime),
        "template": dict(template),
        "template_sha256": str(metadata["template_sha256"]),
        "accepted_plan": accepted_plan,
        "counters": checked_counters,
        "tdg6": tdg6,
        "monitor": monitor,
        "causal": causal,
        "tracer": checked_tracer,
        "provenance_cursor": {
            "source_protocol": source,
            "bridge": provenance["bridge"],
            "payload": dict(checked_payload),
        },
    }
    _canonical(result)
    return result


def _validated_arrays(
    arrays: Mapping[str, np.ndarray], metadata: Mapping[str, Any]
) -> dict[str, np.ndarray]:
    if not isinstance(arrays, Mapping) or tuple(arrays) != ARRAY_NAMES:
        raise HLT16MemberCodecError("snapshot array inventory/order differs")
    checked = {name: _strict_array(arrays[name], name) for name in ARRAY_NAMES}
    points = int(metadata["runtime_identity"]["point_count"])
    labels = int(metadata["tracer"]["label_count"])
    events = int(metadata["tracer"]["event_count"])
    expected = {
        "u": (points, 6),
        "p": (points, 6),
        "q": (points, 6),
        "grid_coordinates": (points,),
        "tracer_labels": (labels,),
        "tracer_positions": (labels,),
        "tracer_proper_times": (labels,),
        "event_proper_times": (events, labels),
        "event_fields": (events, labels, 6),
    }
    if any(checked[name].shape != shape for name, shape in expected.items()):
        raise HLT16MemberCodecError("snapshot array shapes differ")
    template_grid = metadata["template"].get("grid")
    if (
        not isinstance(template_grid, Mapping)
        or array_content_sha256(checked["grid_coordinates"])
        != template_grid.get("coordinates_sha256")
        or np.any(np.diff(checked["tracer_labels"]) <= 0.0)
    ):
        raise HLT16MemberCodecError("grid or tracer labels differ from template")
    return checked


def encode_member(
    member: Proto14RunMember,
    cursor_payload: Mapping[str, object],
    *,
    protocol_artifact_id: str,
    campaign_id: str,
    generation: int,
    provenance_source_protocol: str = GEN0_PROTOCOL,
    executed_plan: Mapping[str, object] | object | None = None,
) -> HLT16MemberSnapshot:
    """Encode one accepted finite boundary and its predecessor provenance."""
    if (
        not isinstance(member, Proto14RunMember)
        or member.amplitude != "3"
        or member.key not in MEMBER_KEYS
        or protocol_artifact_id != RUNTIME_PROTOCOL
    ):
        raise HLT16MemberCodecError("only the frozen PROTO18 GR-0 member is admitted")
    generation = _integer(generation, "generation")
    if provenance_source_protocol == GEN0_PROTOCOL:
        cursor = _validated_proto17_cursor(cursor_payload)
        bridge = (
            "authenticated_GEN0_predecessor"
            if cursor["campaign_id"] == campaign_id
            else "authenticated_GEN0_campaign_rebase"
        )
        if generation != 0 or executed_plan is not None:
            raise HLT16MemberCodecError("GEN0 encode scope differs")
        plan = None
    elif provenance_source_protocol == RUNTIME_PROTOCOL:
        cursor_object = _validated_proto15_cursor(cursor_payload)
        cursor = dict(cursor_object.payload)
        bridge = "active_runtime_predecessor"
        if generation < 1 or executed_plan is None:
            raise HLT16MemberCodecError("evolved encode scope differs")
        try:
            plan = p15._validated_accepted_plan(
                cursor_object, executed_plan, p15._time_identity(member.time)
            )
        except (TypeError, ValueError, p15.Proto15CursorContractError) as error:
            raise HLT16MemberCodecError("executed plan differs") from error
    else:
        raise HLT16MemberCodecError("provenance protocol is forbidden")
    if cursor["member_key"] != member.key:
        raise HLT16MemberCodecError("cursor/member identity differs")
    if provenance_source_protocol == RUNTIME_PROTOCOL and cursor["campaign_id"] != campaign_id:
        raise HLT16MemberCodecError("cursor/campaign identity differs")

    extension, debit = proto14_checkpoint_extension(member)
    arrays = {
        "u": _strict_array(member.state.u, "u"),
        "p": _strict_array(member.state.p, "p"),
        "q": _strict_array(member.state.q, "q"),
        "grid_coordinates": _strict_array(
            member.initial.grid.coordinates, "grid coordinates"
        ),
        "tracer_labels": _strict_array(member.tracers.labels, "tracer labels"),
        "tracer_positions": _strict_array(member.tracers.positions, "tracer positions"),
        "tracer_proper_times": _strict_array(
            member.tracers.proper_times, "tracer proper times"
        ),
        "event_proper_times": _strict_array(
            np.stack(member.tracers.event_proper_times), "event proper history"
        ),
        "event_fields": _strict_array(
            np.stack(member.tracers.event_fields), "event field history"
        ),
    }
    template = _template_identity(member)
    metadata = {
        "runtime_identity": {
            "protocol_artifact_id": RUNTIME_PROTOCOL,
            "campaign_id": campaign_id,
            "branch": "GR-0",
            "amplitude": "3",
            "member_key": member.key,
            "method": member.method_label,
            "point_count": member.point_count,
            "accepted_time": p15._time_identity(member.time),
            "generation": generation,
        },
        "template": template,
        "template_sha256": _digest(template),
        "accepted_plan": plan,
        "counters": {
            "step_index": member.step_index,
            "transaction_serial": member.transaction_serial,
            "accepted_macro_steps": member.temporal_ledger.accepted_macro_step_count,
            "source_retry_count": member.source_retry_count,
            "CFL_retry_count": member.CFL_retry_count,
        },
        "tdg6": {
            "extension": extension,
            "accumulated_debit_hex": [float(item).hex() for item in debit],
        },
        "monitor": asdict(member.transaction.state),
        "causal": asdict(member.transaction.causal_state),
        "tracer": {
            "label_count": arrays["tracer_labels"].size,
            "event_count": arrays["event_fields"].shape[0],
            "cutoff_hex": float(member.tracers.cutoff).hex(),
            "outer_radius_hex": float(member.tracers.outer_radius).hex(),
        },
        "provenance_cursor": {
            "source_protocol": provenance_source_protocol,
            "bridge": bridge,
            "payload": cursor,
        },
    }
    checked_metadata = validate_codec_metadata(metadata)
    checked_arrays = _validated_arrays(arrays, checked_metadata)
    return HLT16MemberSnapshot(checked_arrays, checked_metadata)


def restore_member(member: Proto14RunMember, snapshot: HLT16MemberSnapshot) -> None:
    """Restore only after complete static and dynamic identity validation."""
    if not isinstance(member, Proto14RunMember) or not isinstance(
        snapshot, HLT16MemberSnapshot
    ):
        raise HLT16MemberCodecError("codec types differ")
    metadata = validate_codec_metadata(snapshot.metadata)
    arrays = _validated_arrays(snapshot.arrays, metadata)
    identity = metadata["runtime_identity"]
    if identity["member_key"] != member.key or identity["branch"] != "GR-0":
        raise HLT16MemberCodecError("restore branch/member differs")
    if metadata["template"] != _template_identity(member):
        raise HLT16MemberCodecError("restore shell template differs")

    member.state = EvolutionState(arrays["u"], arrays["p"], arrays["q"])
    member.time = float.fromhex(identity["accepted_time"]["binary64_hex"])
    counters = metadata["counters"]
    member.step_index = counters["step_index"]
    member.transaction_serial = counters["transaction_serial"]
    member.source_retry_count = counters["source_retry_count"]
    member.CFL_retry_count = counters["CFL_retry_count"]
    member.transaction.state = GR0RuntimeMonitorState(**metadata["monitor"])
    member.transaction.causal_state = CausalBudgetState(**metadata["causal"])
    member.tracers.labels = arrays["tracer_labels"].copy()
    member.tracers.positions = arrays["tracer_positions"].copy()
    member.tracers.proper_times = arrays["tracer_proper_times"].copy()
    member.tracers.event_proper_times = [
        row.copy() for row in arrays["event_proper_times"]
    ]
    member.tracers.event_fields = [row.copy() for row in arrays["event_fields"]]
    member.tracers.cutoff = float.fromhex(metadata["tracer"]["cutoff_hex"])
    member.tracers.outer_radius = float.fromhex(
        metadata["tracer"]["outer_radius_hex"]
    )
    debit = np.asarray(
        [float.fromhex(item) for item in metadata["tdg6"]["accumulated_debit_hex"]],
        dtype=np.float64,
    )
    restore_proto14_checkpoint_extension(
        member, metadata["tdg6"]["extension"], debit
    )
    member.__post_init__()


__all__ = [
    "ARRAY_NAMES",
    "GEN0_PROTOCOL",
    "HLT16MemberCodecError",
    "HLT16MemberSnapshot",
    "RUNTIME_PROTOCOL",
    "encode_member",
    "restore_member",
    "validate_codec_metadata",
]
