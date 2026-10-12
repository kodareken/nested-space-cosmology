#!/usr/bin/env python3
"""Run the bounded same-state TDG9 TI2 SSPRK3-tableau counterfactual."""

from __future__ import annotations

import argparse
import ctypes
from copy import deepcopy
from dataclasses import asdict, dataclass, is_dataclass
import errno
from fractions import Fraction
from hashlib import sha256
import json
from math import isfinite
import os
from pathlib import Path
import secrets
import stat
import sys
from typing import Mapping, NoReturn, Sequence

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from scripts import run_fgc_tdg9_ar1 as ar1_runner  # noqa: E402
from scripts import run_fgc_tdg9_loc1 as loc1  # noqa: E402
from scripts import run_fgc_tdg9_loc2 as loc2  # noqa: E402
from recursive_horizons.fgc.evolution import (  # noqa: E402
    hlt16_campaign_runtime as campaign_runtime,
)
from recursive_horizons.fgc.evolution import proto15_runtime as p15  # noqa: E402
from recursive_horizons.fgc.evolution import tdg6_temporal_admission_runtime as tdg6  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_ar1_authority as ar1  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_loc2_pref2_binder as loc2_pref2  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_ti2_authority as authority  # noqa: E402
from recursive_horizons.fgc.evolution.hlt16_campaign_store import (  # noqa: E402
    HLT16CampaignStore,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto19_gr0_static_factory import (  # noqa: E402
    build_static_gr0_shells,
)
from recursive_horizons.fgc.evolution import tdg9_local_extrema as primary  # noqa: E402
from recursive_horizons.fgc.evolution import (  # noqa: E402
    tdg9_local_extrema_independent_v2 as independent,
)


RUNNER_ID = "FGC-1-TDG9-TI2-RUN1"
RAW_SCHEMA = "FGC-1-TDG9-TI2-raw-v1"
TERMINAL_CLASSES = frozenset(
    {
        "completed_all_ten_complete_failures_clear_under_tableau_shadow",
        "completed_one_or_more_complete_failures_persist_under_tableau_shadow",
        "bounded_localization_or_resource_inconclusive",
        "shadow_proposal_premise_stop",
        "invalid_provenance_or_implementation",
    }
)
PUBLISHABLE_TERMINAL_CLASSES = TERMINAL_CLASSES - {
    "invalid_provenance_or_implementation"
}
_OUTPUT_LEAF_LIMIT = 512 * 1024 * 1024
_HEX = frozenset("0123456789abcdef")


class TI1RunnerError(RuntimeError):
    """Typed fail-closed TI1 runner error."""

    def __init__(self, owner: str, code: str, detail: object) -> None:
        self.owner = str(owner)
        self.code = str(code)
        self.detail = " ".join(str(detail).replace(str(ROOT), "<repo>").split())[:640]
        super().__init__(f"{self.owner}/{self.code}: {self.detail}")


class TI1BoundedInconclusive(TI1RunnerError):
    """The frozen exact localization budget did not decide every occurrence."""


def _fail(owner: str, code: str, detail: object) -> NoReturn:
    raise TI1RunnerError(owner, code, detail)


def _canonical(value: object) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()


def _pretty(value: object) -> bytes:
    return (
        json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n"
    ).encode()


def _unique(items: list[tuple[str, object]]) -> dict[str, object]:
    answer: dict[str, object] = {}
    for key, value in items:
        if key in answer:
            _fail("serialization", "duplicate_JSON_key", key)
        answer[key] = value
    return answer


def _json_exact(value: object) -> object:
    if isinstance(value, Fraction):
        return {
            "numerator": str(value.numerator),
            "denominator": str(value.denominator),
        }
    if is_dataclass(value) and not isinstance(value, type):
        return _json_exact(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _json_exact(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_exact(item) for item in value]
    if value is None or isinstance(value, (str, int, bool)):
        return value
    _fail("serialization", "nonexact_value", type(value).__name__)


def _fingerprint_exact(value: object) -> object:
    """Encode only replay-fingerprint binary64 values without decimal loss."""

    if isinstance(value, Fraction):
        return {
            "numerator": str(value.numerator),
            "denominator": str(value.denominator),
        }
    if type(value) is float:
        if not isfinite(value):
            _fail("serialization", "nonfinite_fingerprint_binary64", value)
        return {"binary64_hex": value.hex()}
    if is_dataclass(value) and not isinstance(value, type):
        return _fingerprint_exact(asdict(value))
    if isinstance(value, Mapping):
        return {
            str(key): _fingerprint_exact(item) for key, item in value.items()
        }
    if isinstance(value, (tuple, list)):
        return [_fingerprint_exact(item) for item in value]
    if value is None or type(value) in {str, int, bool}:
        return value
    _fail("serialization", "unsupported_fingerprint_value", type(value).__name__)


def _fingerprint_index(name: str, value: object) -> int:
    if type(value) is not int or value < 0:
        _fail("replay", "negative_or_noninteger_fingerprint_index", name)
    return value


def _execution_authority(root: Path, commit: str) -> authority.TI2Authority:
    try:
        return authority.authorize(
            root,
            authority.read_leaf(root, authority.CONFIG_PATH),
            authority.read_leaf(root, authority.RESULT_PATH),
            commit,
        )
    except Exception as exc:
        raise TI1RunnerError("authority", "rejected", exc) from exc


def _snapshot_store(root: Path) -> tuple[int, str]:
    try:
        return loc1._snapshot_store(root)
    except Exception as exc:
        raise TI1RunnerError("provenance", "store_snapshot", exc) from exc


def _original_occurrences(root: Path) -> dict[tuple[int, str], dict[str, object]]:
    config_raw = authority.read_leaf(root, authority.LOC2_PREF2_CONFIG_PATH)
    result_raw = authority.read_leaf(root, authority.LOC2_PREF2_RESULT_PATH)
    if (
        sha256(config_raw).hexdigest() != authority.LOC2_PREF2_CONFIG_SHA256
        or sha256(result_raw).hexdigest() != authority.LOC2_PREF2_RESULT_SHA256
    ):
        _fail("provenance", "LOC2_PREF2_hash", "compact predecessor differs")
    try:
        result = loc2_pref2.validate_compact_result(config_raw, result_raw)
    except Exception as exc:
        raise TI1RunnerError("provenance", "LOC2_PREF2_compact", exc) from exc
    replay = result["artifact_payload"]["independent_replay"]
    answer: dict[tuple[int, str], dict[str, object]] = {}
    for item in replay["occurrences"]:
        key = (int(item["retry"]), str(item["channel"]))
        complete = {
            str(component["contraction"])
            for component in item["components"]
            if component["component"] == "complete_C"
        }
        if complete != {"sufficient_contraction_failure"}:
            _fail("provenance", "LOC2_complete_failure", key)
        answer[key] = {
            "original_complete_classification": "sufficient_contraction_failure",
            "original_complete_failure": True,
            "original_ownership_class": item["component_ownership"]["classification"],
            "original_ownership_sha256": item["component_ownership_sha256"],
            "original_row_stream_sha256": item["PREF1_row_stream_sha256"],
        }
    if tuple(answer) != authority.FAILED_OCCURRENCES:
        _fail("provenance", "LOC2_occurrence_order", tuple(answer))
    return answer


def _class_name(value: object) -> str:
    kind = type(value)
    return f"{kind.__module__}.{kind.__qualname__}"


def _callable_name(value: object) -> str:
    return f"{getattr(value, '__module__', type(value).__module__)}.{getattr(value, '__qualname__', type(value).__qualname__)}"


def _tracer_history_hash(member: object) -> str:
    arrays = [
        np.asarray(member.tracers.labels, dtype=np.float64),
        np.asarray(member.tracers.positions, dtype=np.float64),
        np.asarray(member.tracers.proper_times, dtype=np.float64),
    ]
    arrays.extend(np.asarray(row, dtype=np.float64) for row in member.tracers.event_proper_times)
    arrays.extend(np.asarray(row, dtype=np.float64) for row in member.tracers.event_fields)
    return array_content_sha256(*arrays)


def _member_fingerprint(member: object, descriptor_sha256: str) -> dict[str, object]:
    ledger = member.temporal_ledger
    if ledger is None:
        _fail("replay", "temporal_ledger_absent", member.key)
    _fingerprint_index("point_count", member.point_count)
    _fingerprint_index("step_index", member.step_index)
    _fingerprint_index("transaction_serial", member.transaction_serial)
    _fingerprint_index(
        "monitor.accepted_stage_count", member.transaction.state.accepted_stage_count
    )
    _fingerprint_index(
        "monitor.last_transaction_serial",
        member.transaction.state.last_transaction_serial,
    )
    failed_serial = member.transaction.state.first_failed_transaction_serial
    if failed_serial is not None:
        _fingerprint_index("monitor.first_failed_transaction_serial", failed_serial)
    return {
        "member_key": member.key,
        "method_label": member.method_label,
        "source_integrator": member.integrator_id,
        "point_count": member.point_count,
        "accepted_time_hex": float(member.time).hex(),
        "state_sha256": array_content_sha256(
            member.state.u, member.state.p, member.state.q
        ),
        "coordinates_sha256": array_content_sha256(member.initial.grid.coordinates),
        "grid_spacing_hex": float(member.initial.grid.spacing).hex(),
        "outer_radius_hex": float(member.initial.grid.maximum).hex(),
        "descriptor_sha256": descriptor_sha256,
        "operator_class": _class_name(member.operator),
        "projector_callable": _callable_name(member.projector),
        "transaction_class": _class_name(member.transaction),
        "tracer_class": _class_name(member.tracers),
        "ledger_class": _class_name(ledger),
        "transaction_sha256": sha256(
            _canonical(
                _fingerprint_exact(
                    {
                        "monitor": asdict(member.transaction.state),
                        "causal": asdict(member.transaction.causal_state),
                        "ledger": asdict(ledger),
                        "step_index": member.step_index,
                        "transaction_serial": member.transaction_serial,
                    }
                )
            )
        ).hexdigest(),
        "tracer_history_sha256": _tracer_history_hash(member),
    }


def _validate_one_variable_fingerprint(
    fingerprint: Mapping[str, object], *, retry: object
) -> None:
    retry_index = _fingerprint_index("retry", retry)
    if (
        fingerprint.get("member_key") != authority.MEMBER_KEY
        or fingerprint.get("method_label") != "RK4"
        or fingerprint.get("source_integrator") != PRIMARY_METHOD
        or fingerprint.get("point_count") != authority.POINT_COUNT
        or fingerprint.get("accepted_time_hex") != authority.ACCEPTED_TIME_HEX
        or fingerprint.get("state_sha256") != authority.PHYSICAL_STATE_SHA256
        or fingerprint.get("descriptor_sha256")
        != authority.MEMBER_DESCRIPTOR_SHA256
        or fingerprint.get("transaction_sha256")
        != authority.TRANSACTION_SHA256_BY_RETRY.get(retry_index)
        or fingerprint.get("coordinates_sha256") != authority.COORDINATES_SHA256
        or fingerprint.get("grid_spacing_hex") != authority.GRID_SPACING_HEX
        or fingerprint.get("outer_radius_hex") != authority.OUTER_RADIUS_HEX
        or str(fingerprint.get("operator_class", "")).split(".")[-1]
        != authority.SEMANTICS["RHS_owner"]
        or not str(fingerprint.get("projector_callable", "")).endswith(
            str(authority.SEMANTICS["projector_owner"])
        )
        or str(fingerprint.get("transaction_class", "")).split(".")[-1]
        != authority.SEMANTICS["transaction_owner"]
        or str(fingerprint.get("tracer_class", "")).split(".")[-1]
        != authority.SEMANTICS["tracer_owner"]
        or str(fingerprint.get("ledger_class", "")).split(".")[-1]
        != authority.SEMANTICS["ledger_owner"]
    ):
        _fail("replay", "one_variable_precondition", retry)
    for key in (
        "coordinates_sha256",
        "descriptor_sha256",
        "transaction_sha256",
        "tracer_history_sha256",
    ):
        value = fingerprint.get(key)
        if (
            not isinstance(value, str)
            or len(value) != 64
            or value.lower() != value
            or any(character not in "0123456789abcdef" for character in value)
        ):
            _fail("replay", "one_variable_precondition", (retry, key))


@dataclass(frozen=True, slots=True)
class RestoredReplay:
    replay: Mapping[str, object]
    member: object
    fingerprint: Mapping[str, object]
    historical_journal_sha256: str


def _restore_replay(
    root: Path,
    store: HLT16CampaignStore,
    shells: Mapping[str, object],
    replay: Mapping[str, object],
) -> RestoredReplay:
    generation = int(replay["predecessor_generation"])
    try:
        checkpoint = store.authenticated_checkpoint_at_generation(generation)
    except Exception as exc:
        raise TI1RunnerError("replay", "checkpoint_unavailable", exc) from exc
    if checkpoint.sha256 != replay["checkpoint_sha256"]:
        _fail("replay", "checkpoint_content_hash", generation)
    relative = (
        f"{ar1.PREF2_STORE_PATH}/checkpoints/"
        f"{generation:020d}-{checkpoint.sha256}.json"
    )
    if sha256(authority.read_leaf(root, relative)).hexdigest() != replay[
        "checkpoint_raw_sha256"
    ]:
        _fail("replay", "checkpoint_raw_hash", generation)
    state = checkpoint.members.get(authority.MEMBER_KEY)
    if state is None:
        _fail("replay", "member_absent", authority.MEMBER_KEY)
    try:
        ar1_runner._require_checkpoint_selection(store, checkpoint, state)
        member = deepcopy(shells[authority.MEMBER_KEY])
        campaign_runtime.restore_member_with_overlay(
            store, checkpoint, member, key=authority.MEMBER_KEY
        )
        cursor = p15.Proto15Cursor(dict(state.cursor))
        cursor.validate()
    except Exception as exc:
        raise TI1RunnerError("replay", "member_restore", exc) from exc
    ledger = member.temporal_ledger
    if (
        ledger is None
        or member.key != authority.MEMBER_KEY
        or member.method_label != "RK4"
        or member.integrator_id != PRIMARY_METHOD
        or member.point_count != authority.POINT_COUNT
        or float(member.time).hex() != authority.ACCEPTED_TIME_HEX
        or ledger.current_macro_step_temporal_retry_count
        != replay["prior_retry_count"]
        or cursor.mode != "RETRY_PENDING"
        or state.pending_owner != "temporal"
    ):
        _fail("replay", "predecessor_identity", replay["retry"])
    fingerprint = _member_fingerprint(member, state.descriptor_sha256)
    _validate_one_variable_fingerprint(fingerprint, retry=replay["retry"])
    try:
        historical = ar1_runner._journal_evidence(root, replay)
    except Exception as exc:
        raise TI1RunnerError("replay", "journal_identity", exc) from exc
    return RestoredReplay(
        replay=dict(replay),
        member=member,
        fingerprint=fingerprint,
        historical_journal_sha256=sha256(_canonical(historical)).hexdigest(),
    )


def _prepare_shadow(restored: RestoredReplay) -> tdg6.TDG6PreparedGR0Compositor:
    member = restored.member
    replay = restored.replay
    ledger = member.temporal_ledger
    if ledger is None:
        _fail("replay", "temporal_ledger_absent", replay["retry"])
    before = _member_fingerprint(member, str(restored.fingerprint["descriptor_sha256"]))
    try:
        prepared = tdg6.prepare_tdg6_gr0_compositor(
            method=COMPARATOR_METHOD,
            time=member.time,
            step_size=float.fromhex(str(replay["attempted_width_hex"])),
            state=member.state,
            rhs=member.operator,
            projector=member.projector,
            transaction=member.transaction,
            tracers=member.tracers,
            coordinates=member.initial.grid.coordinates,
            temporal_ledger=ledger,
            previous_step_index=member.step_index,
            previous_transaction_serial=member.transaction_serial,
        )
    except tdg6.TDG6RefinementPathStop:
        raise
    except Exception as exc:
        raise TI1RunnerError("shadow", "proposal_construction", exc) from exc
    after = _member_fingerprint(member, str(restored.fingerprint["descriptor_sha256"]))
    paths = (prepared.outer, prepared.medium, prepared.fine)
    proposals = tuple(attempt.proposal for path in paths for attempt in path.attempts)
    if (
        before != restored.fingerprint
        or after != before
        or prepared.method != COMPARATOR_METHOD
        or prepared.initial_state_sha256 != authority.PHYSICAL_STATE_SHA256
        or tuple(len(path.attempts) for path in paths) != (1, 2, 4)
        or len(proposals) != 7
        or any(proposal.method != COMPARATOR_METHOD for proposal in proposals)
        or any(len(proposal.stages) != 4 for proposal in proposals)
        or any(
            tuple(record.stage_name for record in proposal.stages)
            != ("ssprk3_s0", "ssprk3_s1", "ssprk3_s2", "candidate_endpoint")
            for proposal in proposals
        )
    ):
        _fail("shadow", "one_variable_identity", replay["retry"])
    return prepared


def _localize_shadow_occurrence(
    rows: Sequence[object], retry: int, channel: str
) -> dict[str, object]:
    component_evidence: dict[str, dict[str, object]] = {}
    serialized: list[dict[str, object]] = []
    for component in ("value_V", "slope_S", "complete_C"):
        component_evidence[component] = {}
        for level, expected in (
            ("D01", 2 * authority.OWNED_ROW_COUNT),
            ("D12", 4 * authority.OWNED_ROW_COUNT),
        ):
            first_cubics, second_cubics = loc2._cubics(rows, level, component)
            if len(first_cubics) != expected or len(second_cubics) != expected:
                _fail("localization", "polynomial_count", (retry, channel, component, level))
            try:
                first = primary.localize_absolute_maximum(
                    first_cubics,
                    maximum_candidates=4 * expected,
                    refinement_depth=authority.PRIMARY_REFINEMENT_DEPTH,
                )
                second = independent.localize_absolute_maximum_independently_v2(
                    second_cubics,
                    maximum_candidates=4 * expected,
                )
            except independent.RootIsolationInconclusive as exc:
                raise TI1BoundedInconclusive(
                    "localization", exc.reason, (retry, channel, component, level)
                ) from exc
            except RuntimeError as exc:
                if "candidate_ceiling_exhausted" in str(exc):
                    raise TI1BoundedInconclusive(
                        "resource", "candidate_ceiling_exhausted", (retry, channel, component, level)
                    ) from exc
                raise
            digest = loc2._primary_stationary_count_digest(first_cubics)
            try:
                loc2._require_route_agreement(
                    first,
                    second,
                    digest,
                    (retry, channel, component, level),
                )
            except loc2.LOC2RunnerError as exc:
                raise TI1RunnerError(
                    "localization", "primary_v2_disagreement", exc.detail
                ) from exc
            component_evidence[component][level] = first
            serialized.append(
                {
                    "component": component,
                    "level": level,
                    "polynomial_count": expected,
                    "candidate_count": first.candidate_count,
                    "primary_stationary_count_stream_sha256": digest,
                    "independent_stationary_count_stream_sha256": (
                        second.stationary_count_stream_sha256
                    ),
                    "primary": _json_exact(first),
                    "independent_v2": _json_exact(second),
                    "survivor_keys_and_intervals_compared": True,
                }
            )
        contraction = loc1._contraction_assessment(
            component_evidence[component]["D01"],
            component_evidence[component]["D12"],
        )
        component_evidence[component]["contraction"] = contraction
        for item in serialized:
            if item["component"] == component:
                item["contraction"] = _json_exact(contraction)
    ownership = loc1._component_ownership(component_evidence)
    complete = str(
        component_evidence["complete_C"]["contraction"]["classification"]
    )
    if complete not in {
        "sufficient_contraction_failure",
        "sufficient_contraction_pass",
        "threshold_inconclusive",
        "exact_zero",
    }:
        _fail("localization", "complete_classification", complete)
    return {
        "retry": retry,
        "channel": channel,
        "row_stream_sha256": loc1._row_hash(rows),
        "components": serialized,
        "shadow_complete_classification": complete,
        "shadow_complete_failure": complete == "sufficient_contraction_failure",
        "shadow_complete_pass": complete
        in {"sufficient_contraction_pass", "exact_zero"},
        "shadow_complete_inconclusive": complete == "threshold_inconclusive",
        "shadow_ownership": _json_exact(ownership),
        "shadow_ownership_class": ownership["classification"],
    }


def _reduce_completed(
    occurrences: Sequence[Mapping[str, object]],
) -> tuple[str, dict[str, int]]:
    if len(occurrences) != 10:
        _fail("reduction", "occurrence_count", len(occurrences))
    transition: dict[str, int] = {}
    for item in occurrences:
        key = f"{item['original_ownership_class']}->{item['shadow_ownership_class']}"
        transition[key] = transition.get(key, 0) + 1
    if any(bool(item["shadow_complete_inconclusive"]) for item in occurrences):
        classification = "bounded_localization_or_resource_inconclusive"
    elif all(bool(item["complete_failure_cleared"]) for item in occurrences):
        classification = (
            "completed_all_ten_complete_failures_clear_under_tableau_shadow"
        )
    else:
        classification = (
            "completed_one_or_more_complete_failures_persist_under_tableau_shadow"
        )
    return classification, dict(sorted(transition.items()))


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


def _read_output_leaf(
    directory_fd: int, name: str, before: os.stat_result
) -> bytes:
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
        raise TI1RunnerError("status", "unsafe_output_leaf", name) from exc
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
            try:
                value = json.loads(raw, object_pairs_hook=_unique)
            except (UnicodeDecodeError, json.JSONDecodeError) as exc:
                raise TI1RunnerError(
                    "status", "noncanonical_output", name
                ) from exc
            if not isinstance(value, dict) or raw != _pretty(value):
                _fail("status", "noncanonical_output", name)
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
            after = os.stat(
                name, dir_fd=descriptors[index - 1], follow_symlinks=False
            )
            if _identity(after) != identities[index]:
                _fail("status", "output_directory_substituted", name)
        if _identity(root.lstat()) != identities[0]:
            _fail("status", "repository_root_substituted", root)
        return values, hashes
    except OSError as exc:
        raise TI1RunnerError(
            "status", "output_snapshot_unreadable", authority.OUTPUT_NAMESPACE
        ) from exc
    finally:
        for descriptor in reversed(descriptors):
            os.close(descriptor)


def _require_keys(
    value: object, expected: set[str], *, label: str
) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != expected:
        _fail("status", "terminal_schema", (label, sorted(expected)))
    return value


def _require_nonnegative_int(value: object, *, label: str) -> int:
    if type(value) is not int or value < 0:
        _fail("status", "terminal_count", label)
    return value


def _require_hex(value: object, *, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or value.lower() != value
        or any(character not in _HEX for character in value)
    ):
        _fail("status", "terminal_digest", label)
    return value


def _require_fraction(value: object, *, label: str) -> Fraction:
    mapping = _require_keys(value, {"numerator", "denominator"}, label=label)
    numerator = mapping["numerator"]
    denominator = mapping["denominator"]
    if not isinstance(numerator, str) or not isinstance(denominator, str):
        _fail("status", "terminal_fraction", label)
    try:
        n = int(numerator)
        d = int(denominator)
    except ValueError as exc:
        raise TI1RunnerError("status", "terminal_fraction", label) from exc
    if str(n) != numerator or str(d) != denominator or d <= 0:
        _fail("status", "terminal_fraction", label)
    return Fraction(n, d)


def _validate_replay_receipts(
    value: object, *, complete: bool
) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        _fail("status", "replay_receipts", "not a mapping")
    keys = tuple(value)
    expected_keys = tuple(str(item["retry"]) for item in authority.REPLAYS)
    if keys != (expected_keys if complete else expected_keys[: len(keys)]):
        _fail("status", "replay_receipts", keys)
    if not complete and len(keys) > 2:
        _fail("status", "replay_receipts", keys)
    fingerprint_keys = {
        "member_key",
        "method_label",
        "source_integrator",
        "point_count",
        "accepted_time_hex",
        "state_sha256",
        "coordinates_sha256",
        "grid_spacing_hex",
        "outer_radius_hex",
        "descriptor_sha256",
        "operator_class",
        "projector_callable",
        "transaction_class",
        "tracer_class",
        "ledger_class",
        "transaction_sha256",
        "tracer_history_sha256",
    }
    receipt_keys = {
        "predecessor",
        "historical_journal_sha256",
        "attempted_width_hex",
        "shadow_path_count",
        "shadow_proposal_count",
        "SSPRK3_stage_and_endpoint_record_count",
    }
    by_retry = {str(item["retry"]): item for item in authority.REPLAYS}
    for key, receipt_raw in value.items():
        receipt = _require_keys(receipt_raw, receipt_keys, label=f"replay {key}")
        predecessor = _require_keys(
            receipt["predecessor"], fingerprint_keys, label=f"predecessor {key}"
        )
        replay = by_retry[key]
        fixed = {
            "member_key": authority.MEMBER_KEY,
            "method_label": "RK4",
            "source_integrator": PRIMARY_METHOD,
            "point_count": authority.POINT_COUNT,
            "accepted_time_hex": authority.ACCEPTED_TIME_HEX,
            "state_sha256": authority.PHYSICAL_STATE_SHA256,
            "coordinates_sha256": authority.COORDINATES_SHA256,
            "grid_spacing_hex": authority.GRID_SPACING_HEX,
            "outer_radius_hex": authority.OUTER_RADIUS_HEX,
            "descriptor_sha256": authority.MEMBER_DESCRIPTOR_SHA256,
        }
        if any(predecessor[name] != expected for name, expected in fixed.items()):
            _fail("status", "replay_predecessor", key)
        if (
            str(predecessor["operator_class"]).split(".")[-1]
            != authority.SEMANTICS["RHS_owner"]
            or not str(predecessor["projector_callable"]).endswith(
                str(authority.SEMANTICS["projector_owner"])
            )
            or str(predecessor["transaction_class"]).split(".")[-1]
            != authority.SEMANTICS["transaction_owner"]
            or str(predecessor["tracer_class"]).split(".")[-1]
            != authority.SEMANTICS["tracer_owner"]
            or str(predecessor["ledger_class"]).split(".")[-1]
            != authority.SEMANTICS["ledger_owner"]
        ):
            _fail("status", "replay_predecessor", key)
        for digest_name in (
            "state_sha256",
            "coordinates_sha256",
            "descriptor_sha256",
            "transaction_sha256",
            "tracer_history_sha256",
        ):
            _require_hex(predecessor[digest_name], label=f"{key}:{digest_name}")
        _require_hex(
            receipt["historical_journal_sha256"],
            label=f"{key}:historical_journal_sha256",
        )
        if (
            receipt["attempted_width_hex"] != replay["attempted_width_hex"]
            or receipt["shadow_path_count"] != 7
            or receipt["shadow_proposal_count"] != 7
            or receipt["SSPRK3_stage_and_endpoint_record_count"] != 28
        ):
            _fail("status", "replay_budget", key)
    return value


def _validate_candidate(
    value: object, *, independent_route: bool, label: str
) -> dict[str, object]:
    keys = {
        "polynomial_ordinal",
        "location",
        "location_ordinal",
        "parameter_lower",
        "parameter_upper",
        "absolute_lower",
        "absolute_upper",
        "metadata",
    }
    if independent_route:
        keys |= {"root_method", "refinement_bits"}
    item = _require_keys(value, keys, label=label)
    polynomial = _require_nonnegative_int(
        item["polynomial_ordinal"], label=f"{label}:polynomial"
    )
    location_ordinal = _require_nonnegative_int(
        item["location_ordinal"], label=f"{label}:location ordinal"
    )
    if item["location"] not in {
        "left_endpoint",
        "right_endpoint",
        "interior_stationary",
    }:
        _fail("status", "terminal_candidate", label)
    lower = _require_fraction(item["parameter_lower"], label=f"{label}:parameter lower")
    upper = _require_fraction(item["parameter_upper"], label=f"{label}:parameter upper")
    absolute_lower = _require_fraction(
        item["absolute_lower"], label=f"{label}:absolute lower"
    )
    absolute_upper = _require_fraction(
        item["absolute_upper"], label=f"{label}:absolute upper"
    )
    if not (0 <= lower <= upper <= 1) or not (
        0 <= absolute_lower <= absolute_upper
    ):
        _fail("status", "terminal_candidate_interval", label)
    metadata = item["metadata"]
    if not isinstance(metadata, list) or any(
        not isinstance(pair, list)
        or len(pair) != 2
        or not isinstance(pair[0], str)
        for pair in metadata
    ):
        _fail("status", "terminal_candidate_metadata", label)
    metadata_keys = [pair[0] for pair in metadata]
    if len(set(metadata_keys)) != len(metadata_keys):
        _fail("status", "terminal_candidate_metadata", label)
    if independent_route and (
        not isinstance(item["root_method"], str)
        or not item["root_method"]
        or type(item["refinement_bits"]) is not int
        or item["refinement_bits"] < 0
    ):
        _fail("status", "terminal_candidate_route", label)
    return {
        "key": (
            polynomial,
            item["location"],
            location_ordinal,
            sha256(_canonical(metadata)).hexdigest(),
        ),
        "parameter_lower": lower,
        "parameter_upper": upper,
        "absolute_lower": absolute_lower,
        "absolute_upper": absolute_upper,
    }


def _validate_localization_pair(
    item: Mapping[str, object], *, expected_polynomials: int, label: str
) -> dict[str, object]:
    primary_keys = {
        "classification",
        "polynomial_count",
        "candidate_count",
        "candidates",
        "global_absolute_lower",
        "global_absolute_upper",
        "maximum_candidates",
        "refinement_depth",
        "evaluator_id",
        "tolerance_used",
    }
    independent_keys = {
        "classification",
        "polynomial_count",
        "candidate_count",
        "candidates",
        "global_absolute_lower",
        "global_absolute_upper",
        "maximum_candidates",
        "initial_refinement_bits",
        "refinement_schedule",
        "global_proof_bit_ceiling",
        "stationary_count_stream_sha256",
        "evaluator_id",
        "tolerance_used",
        "boundary_clipping_used",
    }
    first = _require_keys(item["primary"], primary_keys, label=f"{label}:primary")
    second = _require_keys(
        item["independent_v2"], independent_keys, label=f"{label}:independent"
    )
    candidates = _require_nonnegative_int(
        item["candidate_count"], label=f"{label}:candidate count"
    )
    maximum = 4 * expected_polynomials
    if (
        item["polynomial_count"] != expected_polynomials
        or candidates < expected_polynomials * 2
        or candidates > maximum
        or first["polynomial_count"] != expected_polynomials
        or second["polynomial_count"] != expected_polynomials
        or first["candidate_count"] != candidates
        or second["candidate_count"] != candidates
        or first["maximum_candidates"] != maximum
        or second["maximum_candidates"] != maximum
        or first["refinement_depth"] != authority.PRIMARY_REFINEMENT_DEPTH
        or first["evaluator_id"] != primary.EVALUATOR_ID
        or first["tolerance_used"] is not False
        or second["initial_refinement_bits"] != independent.INITIAL_REFINEMENT_BITS
        or second["refinement_schedule"] != independent.REFINEMENT_SCHEDULE
        or second["global_proof_bit_ceiling"]
        != authority.GLOBAL_PROOF_BIT_CEILING
        or second["evaluator_id"] != independent.EVALUATOR_ID
        or second["tolerance_used"] is not False
        or second["boundary_clipping_used"] is not False
        or item["survivor_keys_and_intervals_compared"] is not True
    ):
        _fail("status", "localization_route", label)
    digest = _require_hex(
        item["primary_stationary_count_stream_sha256"],
        label=f"{label}:primary stationary stream",
    )
    if (
        item["independent_stationary_count_stream_sha256"] != digest
        or second["stationary_count_stream_sha256"] != digest
        or first["classification"] != second["classification"]
        or first["classification"]
        not in {"unique_maximum", "nonunique_or_interval_inconclusive"}
    ):
        _fail("status", "localization_route", label)
    first_candidates_raw = first["candidates"]
    second_candidates_raw = second["candidates"]
    if not isinstance(first_candidates_raw, list) or not isinstance(
        second_candidates_raw, list
    ):
        _fail("status", "localization_survivors", label)
    first_candidates = [
        _validate_candidate(value, independent_route=False, label=label)
        for value in first_candidates_raw
    ]
    second_candidates = [
        _validate_candidate(value, independent_route=True, label=label)
        for value in second_candidates_raw
    ]
    if (
        not first_candidates
        or len(first_candidates) > candidates
        or len(second_candidates) != len(first_candidates)
        or (first["classification"] == "unique_maximum")
        != (len(first_candidates) == 1)
    ):
        _fail("status", "localization_survivors", label)
    first_by_key = {value["key"]: value for value in first_candidates}
    second_by_key = {value["key"]: value for value in second_candidates}
    if (
        len(first_by_key) != len(first_candidates)
        or len(second_by_key) != len(second_candidates)
        or set(first_by_key) != set(second_by_key)
    ):
        _fail("status", "localization_survivors", label)
    for key in first_by_key:
        left = first_by_key[key]
        right = second_by_key[key]
        if max(left["parameter_lower"], right["parameter_lower"]) > min(
            left["parameter_upper"], right["parameter_upper"]
        ) or max(left["absolute_lower"], right["absolute_lower"]) > min(
            left["absolute_upper"], right["absolute_upper"]
        ):
            _fail("status", "localization_survivor_interval", label)
    first_lower = _require_fraction(
        first["global_absolute_lower"], label=f"{label}:primary global lower"
    )
    first_upper = _require_fraction(
        first["global_absolute_upper"], label=f"{label}:primary global upper"
    )
    second_lower = _require_fraction(
        second["global_absolute_lower"], label=f"{label}:v2 global lower"
    )
    second_upper = _require_fraction(
        second["global_absolute_upper"], label=f"{label}:v2 global upper"
    )
    if (
        not (0 <= first_lower <= first_upper)
        or not (0 <= second_lower <= second_upper)
        or max(first_lower, second_lower) > min(first_upper, second_upper)
    ):
        _fail("status", "localization_global_interval", label)
    return {
        "candidate_count": candidates,
        "lower": first_lower,
        "upper": first_upper,
    }


def _contraction_from_levels(
    outer: Mapping[str, object], finest: Mapping[str, object]
) -> dict[str, object]:
    left_lower = outer["lower"] ** 2  # type: ignore[operator]
    left_upper = outer["upper"] ** 2  # type: ignore[operator]
    right_lower = 8 * finest["lower"] ** 2  # type: ignore[operator]
    right_upper = 8 * finest["upper"] ** 2  # type: ignore[operator]
    if outer["upper"] == 0 and finest["upper"] == 0:
        classification = "exact_zero"
    elif right_upper <= left_lower:
        classification = "sufficient_contraction_pass"
    elif right_lower > left_upper:
        classification = "sufficient_contraction_failure"
    else:
        classification = "threshold_inconclusive"
    return {
        "classification": classification,
        "sufficient_contraction_pass_certified": classification
        == "sufficient_contraction_pass",
        "sufficient_contraction_failure_certified": classification
        == "sufficient_contraction_failure",
        "threshold_inconclusive": classification == "threshold_inconclusive",
        "left_D01_squared": {
            "lower": _json_exact(left_lower),
            "upper": _json_exact(left_upper),
        },
        "right_8_D12_squared": {
            "lower": _json_exact(right_lower),
            "upper": _json_exact(right_upper),
        },
    }


def _validate_shadow_components(
    value: object, *, label: str
) -> tuple[dict[str, int], dict[str, str], dict[str, object]]:
    if not isinstance(value, list) or len(value) != 6:
        _fail("status", "shadow_components", label)
    component_keys = {
        "component",
        "level",
        "polynomial_count",
        "candidate_count",
        "primary_stationary_count_stream_sha256",
        "independent_stationary_count_stream_sha256",
        "primary",
        "independent_v2",
        "survivor_keys_and_intervals_compared",
        "contraction",
    }
    order = tuple(
        (component, level)
        for component in ("value_V", "slope_S", "complete_C")
        for level in ("D01", "D12")
    )
    levels: dict[str, dict[str, Mapping[str, object]]] = {
        component: {} for component in ("value_V", "slope_S", "complete_C")
    }
    contractions: dict[str, list[object]] = {
        component: [] for component in levels
    }
    candidate_counts = {component: 0 for component in levels}
    for raw, (component, level) in zip(value, order, strict=True):
        item = _require_keys(raw, component_keys, label=f"{label}:{component}:{level}")
        if item["component"] != component or item["level"] != level:
            _fail("status", "shadow_component_order", label)
        expected_polynomials = (2 if level == "D01" else 4) * authority.OWNED_ROW_COUNT
        evidence = _validate_localization_pair(
            item,
            expected_polynomials=expected_polynomials,
            label=f"{label}:{component}:{level}",
        )
        levels[component][level] = evidence
        contractions[component].append(item["contraction"])
        candidate_counts[component] += int(evidence["candidate_count"])
    classes: dict[str, str] = {}
    for component in levels:
        expected = _contraction_from_levels(
            levels[component]["D01"], levels[component]["D12"]
        )
        if contractions[component] != [expected, expected]:
            _fail("status", "shadow_contraction", f"{label}:{component}")
        classes[component] = str(expected["classification"])
    value_failure = classes["value_V"] == "sufficient_contraction_failure"
    slope_failure = classes["slope_S"] == "sufficient_contraction_failure"
    complete_failure = classes["complete_C"] == "sufficient_contraction_failure"
    value_pass = classes["value_V"] == "sufficient_contraction_pass"
    slope_pass = classes["slope_S"] == "sufficient_contraction_pass"
    value_clear = value_pass or classes["value_V"] == "exact_zero"
    slope_clear = slope_pass or classes["slope_S"] == "exact_zero"
    if complete_failure and value_failure and slope_clear:
        ownership_class = "endpoint_state_owned_failure"
    elif complete_failure and slope_failure and value_clear:
        ownership_class = "width_scaled_RHS_owned_failure"
    elif complete_failure and value_failure and slope_failure:
        ownership_class = "both_components_independently_fail"
    elif complete_failure and value_pass and slope_pass:
        ownership_class = "mixed_only_failure"
    elif complete_failure:
        ownership_class = "component_ownership_inconclusive"
    else:
        ownership_class = "complete_not_failure"
    cancellation: dict[str, bool] = {}
    reinforcement: dict[str, bool] = {}
    for level in ("D01", "D12"):
        value_level = levels["value_V"][level]
        slope_level = levels["slope_S"][level]
        complete_level = levels["complete_C"][level]
        cancellation[level] = complete_level["upper"] < max(  # type: ignore[operator]
            value_level["lower"], slope_level["lower"]
        )
        reinforcement[level] = complete_level["lower"] > max(  # type: ignore[operator]
            value_level["upper"], slope_level["upper"]
        )
    ownership = {
        "classification": ownership_class,
        "value_only_failure": complete_failure and value_failure and slope_clear,
        "slope_only_failure": complete_failure and slope_failure and value_clear,
        "both_components_independently_fail": value_failure and slope_failure,
        "mixed_only_failure": complete_failure and value_pass and slope_pass,
        "cancellation_certified_by_level": cancellation,
        "reinforcement_certified_by_level": reinforcement,
    }
    return candidate_counts, classes, ownership


def _validate_occurrences(
    value: object,
    *,
    originals: Mapping[tuple[int, str], Mapping[str, object]],
    complete: bool,
) -> tuple[list[Mapping[str, object]], dict[str, int], dict[str, int]]:
    if not isinstance(value, list) or (complete and len(value) != 10):
        _fail("status", "occurrences", "count")
    if not complete and len(value) > 9:
        _fail("status", "occurrences", "prefix count")
    occurrence_keys = {
        "original_complete_classification",
        "original_complete_failure",
        "original_ownership_class",
        "original_ownership_sha256",
        "original_row_stream_sha256",
        "retry",
        "channel",
        "row_stream_sha256",
        "components",
        "shadow_complete_classification",
        "shadow_complete_failure",
        "shadow_complete_pass",
        "shadow_complete_inconclusive",
        "shadow_ownership",
        "shadow_ownership_class",
        "complete_failure_cleared",
        "ownership_class_changed",
    }
    counts = {"value_V": 0, "slope_S": 0, "complete_C": 0}
    shadow_counts = {"failure": 0, "pass": 0, "inconclusive": 0}
    expected_prefix = authority.FAILED_OCCURRENCES[: len(value)]
    validated: list[Mapping[str, object]] = []
    for index, (raw, expected_key) in enumerate(
        zip(value, expected_prefix, strict=True)
    ):
        item = _require_keys(raw, occurrence_keys, label=f"occurrence {index}")
        key = (item["retry"], item["channel"])
        if key != expected_key:
            _fail("status", "occurrence_order", (index, key))
        original = originals[expected_key]
        if any(item[name] != expected for name, expected in original.items()):
            _fail("status", "original_occurrence", expected_key)
        _require_hex(item["row_stream_sha256"], label=f"occurrence {index}:row")
        component_counts, component_classes, ownership = _validate_shadow_components(
            item["components"], label=f"occurrence {index}"
        )
        for name in counts:
            counts[name] += component_counts[name]
        complete_class = component_classes["complete_C"]
        failed = complete_class == "sufficient_contraction_failure"
        passed = complete_class in {"sufficient_contraction_pass", "exact_zero"}
        inconclusive = complete_class == "threshold_inconclusive"
        if (
            item["shadow_complete_classification"] != complete_class
            or item["shadow_complete_failure"] is not failed
            or item["shadow_complete_pass"] is not passed
            or item["shadow_complete_inconclusive"] is not inconclusive
            or item["complete_failure_cleared"] is not passed
            or item["shadow_ownership"] != ownership
            or item["shadow_ownership_class"] != ownership["classification"]
            or item["ownership_class_changed"]
            is not (
                item["original_ownership_class"] != ownership["classification"]
            )
            or sum((failed, passed, inconclusive)) != 1
        ):
            _fail("status", "shadow_occurrence", expected_key)
        shadow_counts["failure"] += int(failed)
        shadow_counts["pass"] += int(passed)
        shadow_counts["inconclusive"] += int(inconclusive)
        validated.append(item)
    return validated, counts, shadow_counts


def _validate_completed_terminal(
    terminal: Mapping[str, object],
    *,
    originals: Mapping[tuple[int, str], Mapping[str, object]],
) -> None:
    base_keys = set(
        _terminal_base(
            str(terminal["authority_commit"]),
            str(terminal["classification"]),
            (authority.SEALED_STORE_LEAF_COUNT, authority.SEALED_STORE_SNAPSHOT_SHA256),
            (authority.SEALED_STORE_LEAF_COUNT, authority.SEALED_STORE_SNAPSHOT_SHA256),
        )
    )
    summary_keys = {
        "replay_receipts",
        "occurrence_count",
        "occurrences",
        "ownership_transition_matrix",
        "original_complete_failure_count",
        "shadow_complete_failure_count",
        "shadow_complete_pass_count",
        "shadow_complete_inconclusive_count",
        "component_candidate_counts",
        "total_VSC_candidate_count",
        "shadow_proposal_count",
        "SSPRK3_stage_and_endpoint_record_count",
        "primary_v2_route_agreement_count",
    }
    _require_keys(terminal, base_keys | summary_keys, label="completed terminal")
    _validate_replay_receipts(terminal["replay_receipts"], complete=True)
    occurrences, candidates, shadow = _validate_occurrences(
        terminal["occurrences"], originals=originals, complete=True
    )
    transition: dict[str, int] = {}
    for item in occurrences:
        key = f"{item['original_ownership_class']}->{item['shadow_ownership_class']}"
        transition[key] = transition.get(key, 0) + 1
    transition = dict(sorted(transition.items()))
    classification, reduced_transition = _reduce_completed(occurrences)
    if (
        terminal["classification"] != classification
        or terminal["occurrence_count"] != 10
        or terminal["ownership_transition_matrix"] != transition
        or reduced_transition != transition
        or terminal["original_complete_failure_count"] != 10
        or terminal["shadow_complete_failure_count"] != shadow["failure"]
        or terminal["shadow_complete_pass_count"] != shadow["pass"]
        or terminal["shadow_complete_inconclusive_count"]
        != shadow["inconclusive"]
        or terminal["component_candidate_counts"] != candidates
        or terminal["total_VSC_candidate_count"] != sum(candidates.values())
        or terminal["shadow_proposal_count"] != 21
        or terminal["SSPRK3_stage_and_endpoint_record_count"] != 84
        or terminal["primary_v2_route_agreement_count"] != 60
        or sum(shadow.values()) != 10
        or any(
            count > authority.FAMILY_CANDIDATE_CEILING
            for count in candidates.values()
        )
        or sum(candidates.values()) > authority.AGGREGATE_CANDIDATE_CEILING
    ):
        _fail("status", "completed_terminal_counts", terminal["classification"])


def _validate_terminal(
    root: Path, value: object, *, authority_commit: str
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
            base_keys | {"typed_stop", "replay_receipts", "occurrences"},
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
            or value["occurrences"] != []
        ):
            _fail("status", "premise_stop_terminal", stop)
        _validate_replay_receipts(value["replay_receipts"], complete=False)
        return value
    originals = _original_occurrences(root)
    if classification == "bounded_localization_or_resource_inconclusive" and (
        "typed_stop" in value
    ):
        _require_keys(
            value,
            base_keys
            | {
                "typed_stop",
                "replay_receipts",
                "occurrences",
                "component_candidate_counts",
            },
            label="bounded stop terminal",
        )
        _validate_replay_receipts(value["replay_receipts"], complete=True)
        stop = value["typed_stop"]
        if stop == {"owner": "resource", "code": "aggregate_candidate_ceiling"}:
            occurrences, candidates, _shadow = _validate_occurrences(
                value["occurrences"], originals=originals, complete=True
            )
            if len(occurrences) != 10 or not (
                any(
                    count > authority.FAMILY_CANDIDATE_CEILING
                    for count in candidates.values()
                )
                or sum(candidates.values()) > authority.AGGREGATE_CANDIDATE_CEILING
            ):
                _fail("status", "aggregate_ceiling_terminal", candidates)
        else:
            stop_value = _require_keys(
                stop, {"owner", "code", "detail"}, label="bounded stop"
            )
            allowed = {
                ("resource", "candidate_ceiling_exhausted"),
                ("localization", "root_proof_cap_exceeded"),
                ("localization", "root_location_inconclusive"),
                ("localization", "root_separation_inconclusive"),
                ("localization", "root_count_invariant_failure"),
            }
            if (
                (stop_value["owner"], stop_value["code"]) not in allowed
                or not isinstance(stop_value["detail"], str)
                or not stop_value["detail"]
                or len(stop_value["detail"]) > 640
            ):
                _fail("status", "bounded_stop_terminal", stop_value)
            occurrences, candidates, _shadow = _validate_occurrences(
                value["occurrences"], originals=originals, complete=False
            )
            if len(occurrences) > 9:
                _fail("status", "bounded_stop_terminal", len(occurrences))
        if value["component_candidate_counts"] != candidates:
            _fail("status", "bounded_candidate_counts", candidates)
        return value
    _validate_completed_terminal(value, originals=originals)
    return value


def _inspect_output(
    root: Path, *, authority_commit: str
) -> dict[str, object]:
    snapshot = _snapshot_output(root)
    if snapshot is None:
        return {"state": "absent"}
    values, hashes = snapshot
    manifest = values["manifest.json"]
    if manifest != _manifest(authority_commit):
        _fail("status", "manifest_identity", manifest)
    terminal = _validate_terminal(
        root, values["terminal.json"], authority_commit=authority_commit
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
        "production_SSPRK3_comparator": False,
        "independent_method_agreement": False,
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
    }


def run(root: Path, *, authority_commit: str) -> dict[str, object]:
    repository = root.resolve()
    receipt = _execution_authority(repository, authority_commit)
    if receipt.tableau_selector != "SSPRK3" or receipt.state_advance_authorized:
        _fail("authority", "receipt_semantics", receipt)
    authority.require_output_absent(repository)
    before = _snapshot_store(repository)
    sealed = (
        authority.SEALED_STORE_LEAF_COUNT,
        authority.SEALED_STORE_SNAPSHOT_SHA256,
    )
    if before != sealed:
        _fail("provenance", "sealed_store_identity", (before, sealed))
    originals = _original_occurrences(repository)
    store = HLT16CampaignStore(repository / ar1.PREF2_STORE_PATH)
    shells = build_static_gr0_shells(repository)
    prepared: dict[int, tdg6.TDG6PreparedGR0Compositor] = {}
    replay_receipts: dict[str, object] = {}
    try:
        for replay in authority.REPLAYS:
            restored = _restore_replay(repository, store, shells, replay)
            prepared[int(replay["retry"])] = _prepare_shadow(restored)
            replay_receipts[str(replay["retry"])] = {
                "predecessor": dict(restored.fingerprint),
                "historical_journal_sha256": restored.historical_journal_sha256,
                "attempted_width_hex": replay["attempted_width_hex"],
                "shadow_path_count": 7,
                "shadow_proposal_count": 7,
                "SSPRK3_stage_and_endpoint_record_count": 28,
            }
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
            "replay_receipts": replay_receipts,
            "occurrences": [],
        }
        manifest = _manifest(authority_commit)
        _publish_terminal(
            repository,
            manifest,
            terminal,
            store_before=before,
            store_after=after,
        )
        return terminal

    occurrences: list[dict[str, object]] = []
    candidate_counts = {"value_V": 0, "slope_S": 0, "complete_C": 0}
    try:
        for retry, channel in authority.FAILED_OCCURRENCES:
            surface = loc1._surface(prepared[retry])
            rows = tuple(loc1._rows(surface, channel))
            shadow = _localize_shadow_occurrence(rows, retry, channel)
            original = originals[(retry, channel)]
            occurrence = {
                **original,
                **shadow,
                "complete_failure_cleared": bool(shadow["shadow_complete_pass"]),
                "ownership_class_changed": (
                    original["original_ownership_class"]
                    != shadow["shadow_ownership_class"]
                ),
            }
            occurrences.append(occurrence)
            for component in shadow["components"]:
                candidate_counts[str(component["component"])] += int(
                    component["candidate_count"]
                )
    except TI1BoundedInconclusive as exc:
        after = _snapshot_store(repository)
        if after != before:
            _fail("provenance", "campaign_store_mutated", (before, after))
        terminal = {
            **_terminal_base(
                authority_commit,
                "bounded_localization_or_resource_inconclusive",
                before,
                after,
            ),
            "typed_stop": {"owner": exc.owner, "code": exc.code, "detail": exc.detail},
            "replay_receipts": replay_receipts,
            "occurrences": occurrences,
            "component_candidate_counts": candidate_counts,
        }
        _publish_terminal(
            repository,
            _manifest(authority_commit),
            terminal,
            store_before=before,
            store_after=after,
        )
        return terminal

    total_candidates = sum(candidate_counts.values())
    if (
        any(
            count > authority.FAMILY_CANDIDATE_CEILING
            for count in candidate_counts.values()
        )
        or total_candidates > authority.AGGREGATE_CANDIDATE_CEILING
    ):
        after = _snapshot_store(repository)
        terminal = {
            **_terminal_base(
                authority_commit,
                "bounded_localization_or_resource_inconclusive",
                before,
                after,
            ),
            "typed_stop": {"owner": "resource", "code": "aggregate_candidate_ceiling"},
            "replay_receipts": replay_receipts,
            "occurrences": occurrences,
            "component_candidate_counts": candidate_counts,
        }
        _publish_terminal(
            repository,
            _manifest(authority_commit),
            terminal,
            store_before=before,
            store_after=after,
        )
        return terminal
    classification, transition = _reduce_completed(occurrences)
    after = _snapshot_store(repository)
    if after != before:
        _fail("provenance", "campaign_store_mutated", (before, after))
    terminal = {
        **_terminal_base(authority_commit, classification, before, after),
        "replay_receipts": replay_receipts,
        "occurrence_count": len(occurrences),
        "occurrences": occurrences,
        "ownership_transition_matrix": transition,
        "original_complete_failure_count": sum(
            int(bool(item["original_complete_failure"])) for item in occurrences
        ),
        "shadow_complete_failure_count": sum(
            int(bool(item["shadow_complete_failure"])) for item in occurrences
        ),
        "shadow_complete_pass_count": sum(
            int(bool(item["shadow_complete_pass"])) for item in occurrences
        ),
        "shadow_complete_inconclusive_count": sum(
            int(bool(item["shadow_complete_inconclusive"])) for item in occurrences
        ),
        "component_candidate_counts": candidate_counts,
        "total_VSC_candidate_count": total_candidates,
        "shadow_proposal_count": 21,
        "SSPRK3_stage_and_endpoint_record_count": 84,
        "primary_v2_route_agreement_count": 60,
    }
    _publish_terminal(
        repository,
        _manifest(authority_commit),
        terminal,
        store_before=before,
        store_after=after,
    )
    return terminal


def _manifest(authority_commit: str) -> dict[str, object]:
    return {
        "schema": RAW_SCHEMA,
        "artifact_id": authority.ARTIFACT_ID,
        "runner_id": RUNNER_ID,
        "authority_commit": authority_commit,
        "LOC2_PREF2_result_sha256": authority.LOC2_PREF2_RESULT_SHA256,
        "experiment_label": authority.SEMANTICS["experiment_label"],
        "tableau_selector": "SSPRK3",
        "actual_spatial_operator": "inherited_RK4_2049_SBP4",
        "production_SSPRK3_comparator": False,
        "output_leaves": ["manifest.json", "terminal.json"],
    }


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
    except TI1RunnerError as exc:
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
