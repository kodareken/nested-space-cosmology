"""Production resumable whole-cone field campaign for representative family 14_1.

The driver re-evolves the declared current history from the unchanged
retained upstream source, consumes streamed reconstruction segments, and
encloses history-minus-reference cells. It reuses the authenticated v4
profile payload. Partial, interrupted, method-failed and resource-stopped
output remains OPEN and is not a feasibility verdict or physical certificate.
"""
from copy import deepcopy
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import shutil
import time

import numpy as np

from .nsc_ks_ball_geometry import SubdivisionNeeded
from .nsc_ks_current_field_cone import (
    CHECKPOINT_SCHEMA, V4_PROFILE_IDENTITY, V4_SCHEMA, WholeConeAccumulator,
    WholeConeCheckpoint, WholeConeFieldConfig, canonical_digest,
    matter_adapter_kwargs, packed_radius_tail, prefix_digest_from_cells,
    propagation_adapter_kwargs, restore_v4_profile_payload, scientific_digest,
    scientific_record, segment_identity, with_residual_integrals,
    whole_cone_continuous_inputs, whole_cone_difference_matter_error,
    whole_cone_propagate_difference_error,
)
from .nsc_ks_evaluation_binding import write_bytes_atomic, write_json_atomic
from .nsc_ks_profile_identity import profile_identity
from .nsc_ks_radius_coupling_bounds import with_radius_coupling_bounds
from .nsc_ks_reference_residual import reference_segment_defect
from .nsc_ks_endpoint_contraction import complete_endpoint_inputs, contract_signed_pair
from .nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from .nsc_ks_source_envelope import computational_z_grid, usual_axial_support
from .nsc_ks_streaming_trajectory import StreamSummary, stream_ks_trajectory
from .nsc_ks_trajectory import TrajectorySegment
from .nsc_local_incoming_family import LocalIncomingFamily


SCHEMA = "NSC-KS-WHOLE-CONE-FIELD-v1"
PRODUCTION_FAMILY = (14, 1)
PRODUCTION_PROFILE_IDENTITY = V4_PROFILE_IDENTITY
HISTORY_RELATIVE = "results/development/nsc-ks-gate-history-lm-broyden.json"
V4_RELATIVE = "results/development/nsc-ks-current-field-pilot-v4.json"
CHECKPOINT_INTERVAL = 8
GRID_NODES = 1024
GRID_LENGTH = 0.4
TARGET_NODES = 129
ENCLOSURE_BITS = 90
TIME_DEGREE = 8
RECIPROCAL_ORDER = 4
PRODUCTION_RTOL = 5e-14
PRODUCTION_ATOL = 5e-19
PRODUCTION_MAX_STEP = 1.0 / 8192
FREE_SPACE_FLOOR_BYTES = 60 * 1024 ** 3
STREAM_KIND_TRAJECTORY = "stream_ks_trajectory"
STREAM_KIND_INJECTED = "injected"


class DiagnosticStop(Exception):
    """Bounded diagnostic stop; not a scientific verdict."""


class ResourceStop(Exception):
    """Resource budget exhausted; a valid checkpoint may remain."""


class CampaignBindingError(ValueError):
    """A resume dependency, prefix identity or source binding failed."""


def _sha256_text(value, name):
    if (not isinstance(value, str) or len(value) != 64
            or any(char not in "0123456789abcdef" for char in value)):
        raise ValueError("lowercase SHA-256 " + name + " required")
    return value


def _positive_int(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError("positive integer " + name + " required")
    return value


def _repository_root(root):
    return Path(root).resolve()


def bind_authenticated_v4_payload(root):
    """Restore the immutable v4 coefficient payload without recomputing it."""
    path = _repository_root(root) / V4_RELATIVE
    raw = path.read_bytes()
    record = json.loads(raw)
    if record.get("schema") != V4_SCHEMA:
        raise CampaignBindingError("unexpected field-v4 schema")
    _profiles, digest = restore_v4_profile_payload(
        record, expected_file_sha256=sha256(raw).hexdigest(), file_bytes=raw)
    payload = record.get("profile_coefficient_payload")
    if not isinstance(payload, dict) or not payload:
        raise CampaignBindingError("v4 profile coefficient payload missing")
    return payload, digest, record


def bind_declared_history(root):
    path = _repository_root(root) / HISTORY_RELATIVE
    history = json.loads(path.read_text())
    family = LocalIncomingFamily(np.array(history["history"]["coefficients"]))
    identity = profile_identity(family, include_normal_window=True)
    if identity != PRODUCTION_PROFILE_IDENTITY:
        raise CampaignBindingError("declared history is not the current profile")
    if identity != history.get("profile_identity"):
        raise CampaignBindingError("history identity differs from its record")
    return family, identity, history


def bind_retained_source(root, select_source, family_key=PRODUCTION_FAMILY):
    """Unchanged RetainedUpstreamArchive rows selected by the v2 control."""
    if not callable(select_source):
        raise TypeError("derive_nsc_ks_source_control_v2.select_source required")
    archive = RetainedUpstreamArchive(root)
    source, negative, initial, batch, channel, selections = select_source(
        archive, family_key)
    return {
        "source": source,
        "negative": negative,
        "initial_columns": initial,
        "batch": batch,
        "channel": channel,
        "selections": selections,
        "archive_hashes": dict(archive.input_hashes),
        "source_via_select_source": True,
    }


def production_solver_settings(v4_settings):
    if not isinstance(v4_settings, dict):
        raise TypeError("authenticated v4 settings mapping required")
    return {
        "integrator": "dop853",
        "rtol": PRODUCTION_RTOL.hex(),
        "atol": PRODUCTION_ATOL.hex(),
        "max_step": PRODUCTION_MAX_STEP.hex(),
        "step_control": "joint",
        "tangents": "zero",
        "grid_nodes": GRID_NODES,
        "computational_length": GRID_LENGTH,
        "target_nodes": TARGET_NODES,
        "checkpoint_interval": CHECKPOINT_INTERVAL,
        "v4_profile_settings": dict(v4_settings),
        "history_minus_reference_before_norms": True,
        "profile_payload_recomputed": False,
    }


def computational_period(grid):
    grid = np.asarray(grid, float)
    if grid.ndim != 1 or len(grid) < 8:
        raise ValueError("owned computational envelope grid required")
    spacing = Fraction(grid[1]) - Fraction(grid[0])
    period = spacing * len(grid)
    if period <= 0:
        raise ValueError("positive computational period required")
    return period, float(grid[0])


def production_field_config(family, source, batch, payload_digest, v4_settings,
                            *, bits=ENCLOSURE_BITS):
    if not isinstance(family, LocalIncomingFamily):
        raise TypeError("LocalIncomingFamily required")
    identity = profile_identity(family, include_normal_window=True)
    if identity != PRODUCTION_PROFILE_IDENTITY:
        raise CampaignBindingError("production config requires the declared profile")
    grid = computational_z_grid(GRID_NODES, GRID_LENGTH)
    period, origin = computational_period(grid)
    weights = tuple(float(value) for value in source.column_weights)
    energies = tuple(float(value) for value in source.energies)
    angular = float(batch.angular)
    settings = production_solver_settings(v4_settings)
    settings["rho_up"] = float(batch.rho_up).hex()
    return WholeConeFieldConfig(
        family=PRODUCTION_FAMILY,
        profile_identity=identity,
        bits=bits,
        time_degree=TIME_DEGREE,
        reciprocal_order=RECIPROCAL_ORDER,
        spatial_count=GRID_NODES,
        period_origin=origin,
        period_length_numerator=period.numerator,
        period_length_denominator=period.denominator,
        mass=float(batch.mass),
        angular=angular,
        source_weights=weights,
        source_energies=energies,
        radius_tail=packed_radius_tail(
            family, float(batch.rho_up), angular, RECIPROCAL_ORDER, bits=bits),
        profile_payload_sha256=payload_digest,
        settings=settings,
    ), grid


def production_stream(source, initial, family, grid, target, mass, angular, rho_up):
    """Same DOP853 joint-control stream used by the production campaign."""
    def stream(on_segment):
        return stream_ks_trajectory(
            source, initial, family, grid, target, mass, angular, rho_up,
            axial_support=usual_axial_support(), on_segment=on_segment,
            rtol=PRODUCTION_RTOL, atol=PRODUCTION_ATOL,
            max_step=PRODUCTION_MAX_STEP, tangents="zero",
            step_control="joint", integrator="dop853")
    return stream


def iterable_segment_stream(segments):
    """Test/diagnostic injection: one callback per supplied segment."""
    parts = tuple(segments)
    if any(not isinstance(part, TrajectorySegment) for part in parts):
        raise TypeError("saved anchored TrajectorySegment required")

    def stream(on_segment):
        if not parts:
            raise ValueError("nonempty segment stream required")
        for part in parts:
            on_segment(part)
        return None, StreamSummary(
            len(parts), float(parts[0].rho_start), float(parts[-1].rho_end))
    return stream


def solver_matches_production(settings):
    if not isinstance(settings, dict):
        return False
    try:
        rtol = float.fromhex(settings["rtol"])
        atol = float.fromhex(settings["atol"])
        max_step = float.fromhex(settings["max_step"])
    except (KeyError, TypeError, ValueError):
        return False
    return (
        settings.get("integrator") == "dop853"
        and settings.get("step_control") == "joint"
        and settings.get("tangents") == "zero"
        and rtol == PRODUCTION_RTOL
        and atol == PRODUCTION_ATOL
        and max_step == PRODUCTION_MAX_STEP
        and settings.get("profile_payload_recomputed") is False
    )


def scientific_coverage(*, family, profile_identity, completed_cells,
                        accepted_segments, stream_reached_rho1, last_rho_end,
                        cell_indices, diagnostic, resource_stopped,
                        method_failed, interrupted, source_via_select_source,
                        stream_kind, settings):
    """Truth conditions for all_history_cells and family_14_1_scientific_run."""
    matching = (
        isinstance(completed_cells, int) and not isinstance(completed_cells, bool)
        and completed_cells >= 1
        and accepted_segments == completed_cells
        and stream_reached_rho1
        and float(last_rho_end) == 1.0
        and tuple(cell_indices) == tuple(range(completed_cells))
        and not resource_stopped and not method_failed and not interrupted
    )
    selected_source_run = (
        matching
        and not diagnostic
        and tuple(family) == PRODUCTION_FAMILY
        and profile_identity == PRODUCTION_PROFILE_IDENTITY
        and source_via_select_source
        and stream_kind == STREAM_KIND_TRAJECTORY
        and solver_matches_production(settings)
    )
    return {
        "all_history_cells": bool(matching),
        "selected_source_trajectory_complete": bool(selected_source_run),
        # select_source defaults to four energy rows (12 coherent columns).
        # Temporal completion does not authenticate full-family coverage.
        "full_family_source_coverage": False,
        "family_14_1_scientific_run": False,
        "stream_reached_rho1": bool(stream_reached_rho1),
        "matching_enclosed_cells": bool(
            accepted_segments == completed_cells and completed_cells >= 1
            and tuple(cell_indices) == tuple(range(completed_cells))),
        "accepted_segments": int(accepted_segments),
        "enclosed_cells": int(completed_cells),
    }


def runtime_forecast(*, completed_cells, rho_start, last_rho_end, cpu_seconds):
    """Non-scientific remaining-work estimate; excluded from the digest."""
    if (not isinstance(completed_cells, int) or isinstance(completed_cells, bool)
            or completed_cells < 1 or not np.isfinite(cpu_seconds)
            or cpu_seconds < 0 or not np.isfinite(rho_start)
            or not np.isfinite(last_rho_end)):
        return {
            "estimated_remaining_cells": None,
            "estimated_remaining_cpu_seconds": None,
            "scientific": False,
            "basis": "insufficient",
        }
    spanned = float(rho_start) - float(last_rho_end)
    remaining = float(last_rho_end) - 1.0
    if spanned <= 0 or remaining < 0:
        return {
            "estimated_remaining_cells": 0.0 if remaining == 0 else None,
            "estimated_remaining_cpu_seconds": 0.0 if remaining == 0 else None,
            "scientific": False,
            "basis": "terminal" if remaining == 0 else "nonpositive span",
        }
    cells_per_rho = completed_cells / spanned
    cpu_per_cell = cpu_seconds / completed_cells
    estimated_cells = remaining * cells_per_rho
    return {
        "estimated_remaining_cells": float(estimated_cells),
        "estimated_remaining_cpu_seconds": float(estimated_cells * cpu_per_cell),
        "cells_per_rho": float(cells_per_rho),
        "cpu_seconds_per_cell": float(cpu_per_cell),
        "remaining_rho": float(remaining),
        "scientific": False,
        "basis": "linear in enclosed cells versus remaining rho",
    }


def write_checkpoint_atomic(path, checkpoint):
    if isinstance(checkpoint, WholeConeCheckpoint):
        mapping = checkpoint.to_mapping()
    elif isinstance(checkpoint, dict):
        mapping = WholeConeCheckpoint.from_mapping(checkpoint).to_mapping()
    else:
        raise TypeError("WholeConeCheckpoint required")
    write_json_atomic(path, mapping, indent=2)
    return WholeConeCheckpoint.from_mapping(mapping)


def load_checkpoint(path, config):
    mapping = json.loads(Path(path).read_text())
    checkpoint = WholeConeCheckpoint.from_mapping(mapping)
    checkpoint.validate(config)
    return checkpoint


def _status_text(*, diagnostic, resource_stopped, method_failed, interrupted,
                 all_history_cells, family_run):
    if method_failed:
        return ("OPEN: whole-cone method failed; not a feasibility verdict or "
                "physical certificate; physical rho=1 source error unresolved")
    if interrupted:
        return ("OPEN: whole-cone stream interrupted; not a feasibility verdict "
                "or physical certificate; physical rho=1 source error unresolved")
    if resource_stopped:
        return ("OPEN: whole-cone resource stop; valid checkpoint retained; "
                "not a completed family record; physical rho=1 source error "
                "unresolved")
    if diagnostic:
        return ("OPEN: whole-cone diagnostic stop; valid checkpoint retained; "
                "not a completed family record; physical rho=1 source error "
                "unresolved")
    if family_run:
        return ("OPEN: family 14_1 scientific run enclosed every accepted "
                "cell to rho=1; physical rho=1 source error unresolved")
    if all_history_cells:
        return ("OPEN: every accepted history cell on this stream is enclosed; "
                "selected-source control only; full-family source coverage "
                "and physical rho=1 source error unresolved")
    return ("OPEN: whole-cone field campaign prefix; physical rho=1 source "
            "error unresolved")


@dataclass(frozen=True)
class WholeConeCampaignResult:
    schema: str
    family: tuple
    profile_identity: str
    status: str
    coverage: object
    bounds: object
    continuous_inputs: object
    propagation: object
    matter: object
    prefix_digest: str
    config_digest: str
    stop_reason: str
    source_accuracy_included: bool = False
    physical_rho1_source_error: object = None
    physical_EXISTENCE_certificate: bool = False
    physical_NONEXISTENCE_certificate: bool = False
    family_14_1_scientific_run: bool = False
    all_history_cells: bool = False
    history_minus_reference_before_norms: bool = True
    certificate_use: bool = False
    completed_family_record: bool = False

    def __post_init__(self):
        if self.schema != SCHEMA:
            raise ValueError("unexpected whole-cone campaign schema")
        if not isinstance(self.status, str) or not self.status.startswith("OPEN"):
            raise ValueError("campaign result must remain OPEN")
        if (self.source_accuracy_included
                or self.physical_rho1_source_error is not None
                or self.physical_EXISTENCE_certificate
                or self.physical_NONEXISTENCE_certificate
                or self.certificate_use):
            raise ValueError("campaign cannot issue a physical gate")
        if not self.history_minus_reference_before_norms:
            raise ValueError("history-minus-reference must be formed before norms")
        if not isinstance(self.coverage, dict) or not isinstance(self.bounds, dict):
            raise TypeError("campaign coverage and bounds mappings required")
        _sha256_text(self.profile_identity, "profile identity")
        _sha256_text(self.prefix_digest, "prefix digest")
        _sha256_text(self.config_digest, "config digest")
        coverage_run = bool(self.coverage.get("family_14_1_scientific_run"))
        coverage_all = bool(self.coverage.get("all_history_cells"))
        if self.family_14_1_scientific_run != coverage_run:
            raise ValueError("family_14_1_scientific_run disagrees with coverage")
        if self.all_history_cells != coverage_all:
            raise ValueError("all_history_cells disagrees with coverage")
        if self.family_14_1_scientific_run and not self.all_history_cells:
            raise ValueError("family run requires every accepted cell")
        if self.completed_family_record != self.family_14_1_scientific_run:
            raise ValueError("completed family record tracks the scientific run flag")
        if self.family_14_1_scientific_run and self.stop_reason != "complete":
            raise ValueError("family scientific run requires a complete stream")
        if (self.stop_reason in ("diagnostic", "resource", "method", "interrupted")
                and (self.family_14_1_scientific_run or self.completed_family_record)):
            raise ValueError("partial output cannot be a completed family record")

    def to_mapping(self):
        return {
            "schema": self.schema,
            "family": list(self.family),
            "profile_identity": self.profile_identity,
            "status": self.status,
            "coverage": dict(self.coverage),
            "bounds": dict(self.bounds),
            "continuous_inputs": deepcopy(self.continuous_inputs),
            "propagation": None if self.propagation is None else dict(self.propagation),
            "matter": None if self.matter is None else dict(self.matter),
            "prefix_digest": self.prefix_digest,
            "config_digest": self.config_digest,
            "stop_reason": self.stop_reason,
            "source_accuracy_included": False,
            "physical_rho1_source_error": None,
            "physical_EXISTENCE_certificate": False,
            "physical_NONEXISTENCE_certificate": False,
            "family_14_1_scientific_run": self.family_14_1_scientific_run,
            "all_history_cells": self.all_history_cells,
            "history_minus_reference_before_norms": True,
            "certificate_use": False,
            "completed_family_record": self.completed_family_record,
        }


def _empty_bounds():
    return {
        "polynomial_sum": None,
        "remainder_sum": None,
        "total_continuous_normalized_residual_sum": None,
        "residual_integrals": None,
        "completed_cells": 0,
    }


def _bounds_from_accumulator(accumulator):
    checkpoint = accumulator.checkpoint()
    mapping = checkpoint.to_mapping()
    return {
        "polynomial_sum": mapping["polynomial_sum"],
        "remainder_sum": mapping["remainder_sum"],
        "total_continuous_normalized_residual_sum": mapping["total_sum"],
        "residual_integrals": mapping["residual_integrals"],
        "completed_cells": mapping["completed_cells"],
        "last_cell_total": mapping["cells"][-1]["total"],
        "last_rho_end": mapping["last_rho_end"],
    }


def _free_bytes(path):
    return int(shutil.disk_usage(Path(path).parent if Path(path).suffix else path).free)


def run_whole_cone_campaign(
        config, profile_payload, family, stream, *,
        checkpoint_path, resume=False, max_new_cells=None, cpu_budget=None,
        free_space_floor=None, free_space=None,
        clock=time.process_time, checkpoint_interval=CHECKPOINT_INTERVAL,
        stream_kind=STREAM_KIND_INJECTED, source_via_select_source=False,
        continuous_source=None, rho_up=None):
    """Enclose streamed cells, checkpoint every eight new cells, and stay OPEN.

    `stream` is `callable(on_segment) -> (prepared, StreamSummary)`. Resume
    re-streams from the same start and authenticates every stored dense-segment
    identity before a new cell is enclosed. Diagnostic and resource stops write
    a valid checkpoint and a forecast; they cannot emit a completed family
    record.
    """
    if not isinstance(config, WholeConeFieldConfig):
        raise TypeError("WholeConeFieldConfig required")
    if not isinstance(family, LocalIncomingFamily):
        raise TypeError("LocalIncomingFamily required")
    if not callable(stream):
        raise TypeError("segment stream callable required")
    checkpoint_path = Path(checkpoint_path)
    interval = _positive_int(checkpoint_interval, "checkpoint interval")
    if max_new_cells is not None:
        max_new_cells = _positive_int(max_new_cells, "max new cells")
    if cpu_budget is not None:
        cpu_budget = float(cpu_budget)
        if not np.isfinite(cpu_budget) or cpu_budget <= 0:
            raise ValueError("positive cpu budget required")
    started = float(clock())
    started_rho = None
    diagnostic_mode = max_new_cells is not None
    resource_stopped = False
    method_failed = False
    interrupted = False
    stop_reason = "complete"
    prepared = None
    summary = None
    new_cells = 0
    accepted = 0
    last_rho = None
    endpoint_reference = None
    endpoint_difference = None
    # Recomputed on every replayed segment, including the saved prefix. This
    # keeps the existing immutable difference checkpoint format usable.
    reference_integrals = []
    if resume:
        checkpoint = load_checkpoint(checkpoint_path, config)
        accumulator = WholeConeAccumulator.from_checkpoint(
            checkpoint, config, profile_payload, family)
        pending_replay = checkpoint.completed_cells
    else:
        if checkpoint_path.exists():
            raise FileExistsError("checkpoint exists; pass resume=True")
        accumulator = WholeConeAccumulator(config, profile_payload, family)
        pending_replay = 0

    def remaining_replay():
        return len(accumulator._resume_cells_pending)

    def maybe_resource_stop():
        nonlocal resource_stopped, stop_reason
        if cpu_budget is not None and float(clock()) - started >= cpu_budget:
            resource_stopped = True
            stop_reason = "resource"
            raise ResourceStop("cpu budget exhausted")
        if free_space_floor:
            probe = free_space or (lambda: _free_bytes(checkpoint_path))
            if int(probe()) < int(free_space_floor):
                resource_stopped = True
                stop_reason = "resource"
                raise ResourceStop("free-space floor")

    def persist():
        if not accumulator._cells:
            return None
        return write_checkpoint_atomic(checkpoint_path, accumulator.checkpoint())

    def on_segment(segment):
        nonlocal accepted, new_cells, last_rho, started_rho, stop_reason
        nonlocal endpoint_reference, endpoint_difference
        if not isinstance(segment, TrajectorySegment):
            raise TypeError("saved anchored TrajectorySegment required")
        if started_rho is None:
            started_rho = float(segment.rho_start)
        identity = segment_identity(segment, accepted)
        if remaining_replay():
            stored = accumulator.replay_segment(segment, cell_index=accepted)
            if any(stored.get(name) != identity[name] for name in identity):
                raise CampaignBindingError("segment replay")
        else:
            maybe_resource_stop()
            if max_new_cells is not None and new_cells >= max_new_cells:
                stop_reason = "diagnostic"
                raise DiagnosticStop("max new cells")
            accumulator.accumulate_segment(segment, cell_index=accepted)
            new_cells += 1
            if new_cells % interval == 0:
                persist()
        last_rho = float(segment.rho_end)
        if segment.rho_end == 1.0:
            size = 2*len(config.source_weights)
            endpoint_reference = segment.end[:size].reshape(2,-1).copy()
            endpoint_difference = segment.end[size:size*(1+config.spatial_count)].reshape(
                2,len(config.source_weights),config.spatial_count).copy()
        reference_integrals.append(reference_segment_defect(
            segment, config.source_energies, config.source_weights,
            config.mass, config.angular, bits=config.bits
        )["weighted_row_residual_integral"])
        accepted += 1
        del segment

    try:
        prepared, summary = stream(on_segment)
    except DiagnosticStop:
        stop_reason = "diagnostic"
        persist()
    except ResourceStop:
        resource_stopped = True
        stop_reason = "resource"
        persist()
    except (ArithmeticError, SubdivisionNeeded):
        method_failed = True
        stop_reason = "method"
        persist()
    except Exception:
        interrupted = True
        stop_reason = "interrupted"
        persist()
        raise
    else:
        if remaining_replay():
            raise CampaignBindingError(
                "checkpoint prefix replay did not consume every stored cell")
        if summary is None or not isinstance(summary, StreamSummary):
            raise CampaignBindingError("stream summary required")
        if accepted != summary.accepted_steps:
            raise CampaignBindingError("accepted segment count disagrees")
        persist()

    cpu_seconds = max(0.0, float(clock()) - started)
    cells = list(accumulator._cells)
    completed = len(cells)
    stream_reached_rho1 = (
        stop_reason == "complete"
        and last_rho == 1.0
        and summary is not None
        and float(summary.rho_end) == 1.0
        and remaining_replay() == 0
    )
    coverage = scientific_coverage(
        family=config.family,
        profile_identity=config.profile_identity,
        completed_cells=completed,
        accepted_segments=accepted if stop_reason == "complete" else accepted,
        stream_reached_rho1=stream_reached_rho1,
        last_rho_end=1.0 if last_rho == 1.0 else (last_rho if last_rho is not None else 0.0),
        cell_indices=[cell["cell_index"] for cell in cells],
        diagnostic=diagnostic_mode or stop_reason == "diagnostic",
        resource_stopped=resource_stopped,
        method_failed=method_failed,
        interrupted=interrupted or stop_reason == "interrupted",
        source_via_select_source=source_via_select_source,
        stream_kind=stream_kind,
        settings=dict(config.settings),
    )
    if stop_reason != "complete":
        coverage["all_history_cells"] = False
        coverage["family_14_1_scientific_run"] = False
    coverage["diagnostic"] = stop_reason == "diagnostic"
    coverage["diagnostic_mode"] = diagnostic_mode
    coverage["resource_stopped"] = resource_stopped
    coverage["method_failed"] = method_failed
    coverage["interrupted"] = interrupted
    coverage["whole_spatial_period"] = True
    coverage["source_columns"] = len(config.source_weights)
    coverage["physical_rho1_source_error"] = None
    coverage["checkpoint_interval"] = interval
    coverage["new_cells_this_session"] = new_cells

    channel = None
    covariance = None
    weights = config.source_weights
    amplitudes = None
    if continuous_source is not None:
        channel = continuous_source.get("channel")
        source_obj = continuous_source.get("source")
        if source_obj is not None:
            covariance = source_obj.covariance
            weights = source_obj.column_weights
    if prepared is not None:
        amplitudes = getattr(prepared, "reference_amplitudes", None)
    if rho_up is not None:
        bound_rho = float(rho_up)
    elif isinstance(config.settings.get("rho_up"), str):
        bound_rho = float.fromhex(config.settings["rho_up"])
    else:
        bound_rho = 1.03
    inputs = whole_cone_continuous_inputs(
        family, mass=config.mass, angular=config.angular, rho_up=bound_rho,
        source_energies=config.source_energies, source_weights=weights,
        source_covariance=covariance, channel=channel,
        reference_amplitudes=amplitudes, bits=config.bits)
    inputs = with_radius_coupling_bounds(
        inputs, family, bound_rho, config.angular, bits=config.bits)
    inputs["reference_residual_enclosed_cells"] = len(reference_integrals)
    if (stream_reached_rho1 and len(reference_integrals) == completed
            and coverage["all_history_cells"]):
        from flint import arb, ctx
        from .nsc_ks_ball_trajectory import exact_upper, restored_upper
        with ctx.workprec(config.bits):
            total_reference = sum((restored_upper(v) for v in reference_integrals),arb(0))
            inputs["propagation"]["reference_residual"] = {
                "value": exact_upper(total_reference),
                "owner": "reference_segment_defect over every authenticated segment",
                "reason": None,
            }
            inputs["propagation"]["reference_initial"] = {
                "value": exact_upper(arb(0)),
                "owner": "same stored columns, exact only for numerical field subproblem",
                "reason": "physical preparation error is separate and remains unbounded",
            }
        inputs["reference_error_scope"] = (
            "numerical evolution from fixed supplied columns; excludes physical source error")
    bounds = _empty_bounds()
    prefix = ""
    propagation = None
    matter = None
    if completed:
        bounds = _bounds_from_accumulator(accumulator)
        inputs = with_residual_integrals(inputs, bounds["residual_integrals"])
        prefix = prefix_digest_from_cells(
            config.digest(), cells,
            polynomial_sum=bounds["polynomial_sum"],
            remainder_sum=bounds["remainder_sum"],
            total_sum=bounds["total_continuous_normalized_residual_sum"],
            residual_integrals=bounds["residual_integrals"])
        propagate_kwargs = propagation_adapter_kwargs(inputs)
        if propagate_kwargs is not None:
            propagation = whole_cone_propagate_difference_error(**propagate_kwargs)
        if propagation is not None and stream_reached_rho1:
            if endpoint_reference is None or endpoint_difference is None:
                raise CampaignBindingError('completed stream is missing its captured endpoint')
            inputs = complete_endpoint_inputs(
                inputs,endpoint_reference,endpoint_difference,config.source_weights,
                config.source_energies,Fraction(config.period_length_numerator,
                                               config.period_length_denominator),
                propagation,bits=config.bits)
        matter_kwargs = matter_adapter_kwargs(inputs)
        if matter_kwargs is not None:
            negative = None if continuous_source is None else continuous_source.get('negative')
            if negative is not None:
                if (not np.array_equal(negative.energies,-np.asarray(config.source_energies))
                        or not np.array_equal(negative.column_weights,config.source_weights)):
                    raise CampaignBindingError('negative source energy/weight binding differs')
            matter = contract_signed_pair(matter_kwargs,
                None if negative is None else negative.covariance,bits=config.bits)
    if not prefix:
        prefix = canonical_digest({
            "config_digest": config.digest(),
            "cells": cells,
            "stop_reason": stop_reason,
        })
    result = WholeConeCampaignResult(
        schema=SCHEMA,
        family=config.family,
        profile_identity=config.profile_identity,
        status=_status_text(
            diagnostic=stop_reason == "diagnostic",
            resource_stopped=resource_stopped,
            method_failed=method_failed,
            interrupted=interrupted,
            all_history_cells=coverage["all_history_cells"],
            family_run=coverage["family_14_1_scientific_run"]),
        coverage=coverage,
        bounds=bounds,
        continuous_inputs=inputs,
        propagation=propagation,
        matter=matter,
        prefix_digest=prefix,
        config_digest=config.digest(),
        stop_reason=stop_reason,
        family_14_1_scientific_run=coverage["family_14_1_scientific_run"],
        all_history_cells=coverage["all_history_cells"],
        completed_family_record=coverage["family_14_1_scientific_run"],
    )
    forecast = runtime_forecast(
        completed_cells=max(completed, 1) if completed else 0,
        rho_start=started_rho if started_rho is not None else float("nan"),
        last_rho_end=last_rho if last_rho is not None else float("nan"),
        cpu_seconds=cpu_seconds)
    runtime = {
        "CPU_seconds": cpu_seconds,
        "accepted_segments": accepted,
        "new_cells": new_cells,
        "replayed_cells": pending_replay if stop_reason != "complete" else pending_replay,
    }
    mapping = result.to_mapping()
    mapping["checkpoint_schema"] = CHECKPOINT_SCHEMA
    mapping["checkpoint_path"] = str(checkpoint_path)
    mapping["runtime"] = runtime
    mapping["forecast"] = forecast
    mapping["scientific_digest"] = scientific_digest(mapping)
    if scientific_record(mapping).get("runtime") is not None:
        raise ValueError("runtime leaked into the scientific record")
    if scientific_record(mapping).get("forecast") is not None:
        raise ValueError("forecast leaked into the scientific record")
    return mapping, result
