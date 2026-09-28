"""GR-0 member construction at the PROTO19 restart boundary.

The authenticated generation-zero import is deliberately separate from the
post-bridge shell factory.  The former is the *only* path permitted to reopen
the two historical raw containers.  Once HLT16 has published generation one,
restart code constructs static GR-0 shells and hydrates their evolving fields
only through a verified persisted descriptor/payload pair.

Neither path is an event runner: this module has no proposal, commit, or
checkpoint API.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
import sys
import tomllib
from typing import Any, Mapping

import numpy as np

from .proto18_auth1_inputs import reimport_matches_pref26_evidence, source_identities_from_pref26_result
from .proto18_pref27_binder import Pref27Evidence, verify_generation_zero_store
from .protocol_v17 import MEMBER_KEYS
from .numerical_engine import EvolutionState, array_content_sha256
from .proto14_runtime import Proto14RunMember
from .proto5_runtime import GR0RuntimeMonitorState
from .boundary_domain import CausalBudgetState


BRANCH = "GR-0"
AMPLITUDE = "3"
RESTART_TIME = 23.0 / 16.0
TARGET_TIME = 1.5


class Proto19ProgressionInputError(ValueError):
    """The sealed GEN0 boundary cannot lawfully seed GR-0 progression."""


@dataclass(frozen=True, slots=True)
class ReconstructedGR0MemberSet:
    """Immutable live GR-0 members at event 23, suitable for a later adapter."""

    evidence: Pref27Evidence
    members: Mapping[str, object]
    source_identity: Mapping[str, Mapping[str, Any]]
    restart_time: float = RESTART_TIME
    target_time: float = TARGET_TIME
    branch: str = BRANCH
    amplitude: str = AMPLITUDE

    def __post_init__(self) -> None:
        if self.branch != BRANCH or self.amplitude != AMPLITUDE:
            raise Proto19ProgressionInputError("candidate branch or amplitude is forbidden")
        if tuple(self.members) != MEMBER_KEYS:
            raise Proto19ProgressionInputError("live members are not in sealed canonical order")
        if not (np.float64(self.restart_time).tobytes() == np.float64(RESTART_TIME).tobytes()
                and np.float64(self.target_time).tobytes() == np.float64(TARGET_TIME).tobytes()):
            raise Proto19ProgressionInputError("event 23 to 24 time contract differs")


def _root() -> Path:
    return Path(__file__).resolve().parents[4]


def _load_pref26(root: Path) -> Mapping[str, Any]:
    raw = json.loads((root / "results/fgc-1-pro18-pref26.json").read_text("utf-8"))
    if raw.get("artifact_id") != "FGC-1-PRO18-PREF26":
        raise Proto19ProgressionInputError("PREF26 result identity differs")
    return raw


def _same_float(left: float, right: float) -> bool:
    return np.float64(left).tobytes() == np.float64(right).tobytes()


def _check_member(key: str, member: object, evidence: Pref27Evidence) -> None:
    checkpoint = evidence.checkpoint
    descriptor = next(
        item for item in evidence.genesis_spec["member_descriptors"] if item["member_key"] == key
    )
    cursor = checkpoint["cursors"][key]
    ledger = checkpoint["ledgers"][key]
    state_object = descriptor["state_object"]
    required = ("state", "time", "step_index", "transaction_serial", "CFL_retry_count",
                "source_retry_count", "transaction", "tracers", "initial", "operator", "projector",
                "method_label", "point_count", "amplitude")
    if any(not hasattr(member, name) for name in required):
        raise Proto19ProgressionInputError(f"{key} runtime shell omits a required field")
    state = member.state
    if not isinstance(state, EvolutionState) or not _same_float(member.time, RESTART_TIME):
        raise Proto19ProgressionInputError(f"{key} live state/time differs")
    if (member.method_label != descriptor["method"] or member.point_count != descriptor["point_count"]
            or member.amplitude != AMPLITUDE or member.step_index != descriptor["step_index"]
            or member.transaction_serial != descriptor["transaction_serial"]
            or member.CFL_retry_count != descriptor["state_object"]["CFL_retry_count"]
            or member.source_retry_count != descriptor["state_object"]["source_retry_count"]
            or array_content_sha256(state.u, state.p, state.q) != state_object["physical_state_sha256"]
            or cursor["accepted_state_sha256"] != descriptor["state_object_canonical_sha256"]
            or cursor["mode"] != "FRESH_READY" or cursor["committed_common_event_index"] != 23
            or ledger["accepted_macro_step_count"] != 0
            or ledger["cumulative_temporal_retry_count"] != 0
            or ledger["current_macro_step_temporal_retry_count"] != 0
            or any(float.fromhex(value) != 0.0 for value in ledger["accumulated_debit_vector_hex"])):
        raise Proto19ProgressionInputError(f"{key} sealed runtime/cursor/ledger identity differs")
    if (member.temporal_ledger.accepted_macro_step_count != 0
            or member.temporal_ledger.cumulative_temporal_retry_count != 0
            or member.temporal_ledger.current_macro_step_temporal_retry_count != 0
            or member.temporal_ledger.serialized_temporal_rejections
            or any(member.temporal_ledger.accumulated_debit_vector)):
        raise Proto19ProgressionInputError(f"{key} live TDG6 ledger differs from sealed generation zero")
    transaction = member.transaction
    if (transaction.state != type(transaction.state)(**state_object["runtime_monitor_state"])
            or transaction.causal_state != type(transaction.causal_state)(**state_object["causal_state"])
            or len(member.tracers.event_fields) != state_object["event_history_state"]["sample_count"]
            or state.shape != tuple(state_object["physical_arrays"][0]["shape"])):
        raise Proto19ProgressionInputError(f"{key} monitor, causal, tracer, or grid identity differs")
    before = array_content_sha256(state.u, state.p, state.q)
    rhs = member.operator(RESTART_TIME, state)
    if (rhs.diagnostics.get("source_raw_gate_passed") is not True
            or rhs.diagnostics.get("SRC4_tensor_contracted_reference_source") is not True
            or array_content_sha256(state.u, state.p, state.q) != before):
        raise Proto19ProgressionInputError(f"{key} GR-0 source preflight differs")


def build_static_gr0_shells(repository_root: Path | None = None) -> Mapping[str, Proto14RunMember]:
    """Construct only the fixed GR-0 runtime shell topology.

    This deliberately does not authenticate, open, or import a historical
    generation-zero raw container.  ``restore_member_with_overlay`` validates
    a persisted HLT16 descriptor against this static template before it copies
    any evolved field, monitor, tracer, or counter into the shell.
    """
    root = _root() if repository_root is None else Path(repository_root).resolve()
    sys.path[:0] = [str(root), str(root / "src")]
    from scripts import run_fgc_gr0_calibration_v13 as proto13  # GR-0-only static owner

    shells = proto13._build_frozen_member_shells()
    if tuple(shells) != MEMBER_KEYS:
        raise Proto19ProgressionInputError("static GR-0 shell order differs")
    ordered: dict[str, Proto14RunMember] = {}
    for key in MEMBER_KEYS:
        shell = Proto14RunMember.from_proto7(shells[key])
        if shell.key != key or shell.amplitude != AMPLITUDE:
            raise Proto19ProgressionInputError("static GR-0 shell identity differs")
        ordered[key] = shell
    return ordered


def reconstruct_gr0_members(repository_root: Path | None = None) -> ReconstructedGR0MemberSet:
    """Authenticate GEN0 and return six live GR-0 shells without mutation."""
    root = _root() if repository_root is None else Path(repository_root).resolve()
    evidence = verify_generation_zero_store(root)
    pref26 = _load_pref26(root)
    auth_evidence = json.loads(
        (root / "results/fgc-1-pro18-auth1.json").read_text("utf-8")
    )["artifact_payload"]["auth1_input_evidence"]
    identities = source_identities_from_pref26_result(pref26, auth_evidence)
    with (root / "configs/fgc/fgc-1-pro18-frz1.toml").open("rb") as handle:
        pro18 = tomllib.load(handle)
    corpus = reimport_matches_pref26_evidence(
        root, pro18, auth_evidence, expected_source_identities=identities
    )
    # Reuse the fixed PROTO13 shell/operator/projector/tracer construction;
    # this is distinct from the one authorized historical raw import above.
    shells = build_static_gr0_shells(root)
    ordered: dict[str, Proto14RunMember] = {}
    for key in MEMBER_KEYS:
        shell = shells[key]
        verified = corpus.members[key]
        metadata = verified.metadata
        arrays = verified.arrays
        shell.state = EvolutionState(arrays["u"], arrays["p"], arrays["q"])
        shell.time = float(metadata["time"])
        shell.step_index = int(metadata["step_index"])
        shell.transaction_serial = int(metadata["transaction_serial"])
        shell.CFL_retry_count = int(metadata["CFL_retry_count"])
        shell.source_retry_count = int(metadata["source_retry_count"])
        shell.transaction.state = GR0RuntimeMonitorState(**dict(metadata["runtime_monitor_state"]))
        shell.transaction.causal_state = CausalBudgetState(**dict(metadata["causal_state"]))
        shell.tracers.positions = np.asarray(arrays["tracer_positions"], dtype=np.float64).copy()
        shell.tracers.proper_times = np.asarray(arrays["tracer_proper_times"], dtype=np.float64).copy()
        shell.tracers.event_proper_times = [row.copy() for row in arrays["event_proper_times"]]
        shell.tracers.event_fields = [row.copy() for row in arrays["event_fields"]]
        # ``build_static_gr0_shells`` intentionally creates no evolved
        # temporal state.  The authenticated GEN0 bridge supplies its accepted
        # time, so derive a fresh zero ledger at that exact imported boundary.
        shell.temporal_ledger = None
        shell.__post_init__()
        ordered[key] = shell
    for key, member in ordered.items():
        _check_member(key, member, evidence)
    return ReconstructedGR0MemberSet(evidence=evidence, members=ordered, source_identity=identities)


__all__ = ["AMPLITUDE", "BRANCH", "RESTART_TIME", "TARGET_TIME", "Proto19ProgressionInputError", "ReconstructedGR0MemberSet", "build_static_gr0_shells", "reconstruct_gr0_members"]
