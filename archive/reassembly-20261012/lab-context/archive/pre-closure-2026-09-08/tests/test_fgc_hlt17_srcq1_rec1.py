"""Synthetic/mocked SRCQ1-REC1 tests. No physical source, live store, or production grid."""

from __future__ import annotations

import ast
from dataclasses import dataclass
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from recursive_horizons.fgc.evolution.hlt17_admission_runtime import (  # noqa: E402
    hlt17_implementation_identity,
)
from recursive_horizons.fgc.evolution.hlt17_imp1_cursor import (  # noqa: E402
    plan_imp1_lattice,
)
from recursive_horizons.fgc.evolution.hlt17_srcq1 import (  # noqa: E402
    STATIC_INPUT_PATHS,
    ResourceCeilings,
    derive_prospective_requested_caps,
    evaluate_bound_source_once,
)
from recursive_horizons.fgc.evolution.pro20_origin import (  # noqa: E402
    Pro20OriginError,
)
from recursive_horizons.fgc.evolution.pro20_source_factory import (  # noqa: E402
    Pro20SourceFactoryError,
)
from recursive_horizons.fgc.evolution.hlt17_srcq1_rec1 import (  # noqa: E402
    ARTIFACT_ID,
    ATTEMPT2_ARTIFACT_ID,
    ATTEMPT2_AUTHORITY_COMMIT,
    ATTEMPT2_CLASSIFICATION,
    ATTEMPT2_COMPACT_RELATIVE,
    ATTEMPT2_COMPACT_SHA256,
    ATTEMPT2_DERIVATION_DOCUMENT,
    ATTEMPT2_IMPLEMENTATION_COMMIT,
    ATTEMPT2_IMPLEMENTATION_SHA256,
    ATTEMPT2_METADATA_LINKED_COMPACT_SHA256,
    ATTEMPT2_RAW_DIRECTORY_PAYLOAD_SHA256,
    ATTEMPT2_RAW_LEAF_BYTES,
    ATTEMPT2_RAW_LEAF_RELATIVE,
    ATTEMPT2_RAW_LEAF_SHA256,
    ATTEMPT2_SEALED_COMPACT_SHA256,
    AUTHORITY_SEAM,
    FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER,
    HLT17SRCQ1REC1AuthoritySeamError,
    HLT17SRCQ1REC1CapError,
    HLT17SRCQ1REC1Error,
    HLT17SRCQ1REC1IdentityError,
    HLT17SRCQ1REC1InconclusiveError,
    HLT17SRCQ1REC1PublicationError,
    ORIGIN_CAPTURE_SHA256,
    PREDECESSOR_MEMBER_KEYS,
    RESUME_MEMBER_KEYS,
    SCHEMA,
    authenticate_attempt2_predecessor,
    implementation_identity,
    isolate_working_checkpoint,
    load_coordinator_authority_delta_validator,
    plan_status,
    predecessor_rk_evidence_references,
    publish_qualification_result,
    qualify_and_publish,
    qualify_srcq1_rec1,
    recommended_resource_ceilings,
    refuse_endpoint_payload,
    require_authority_delta,
    require_independent_cfl_ratio,
    retain_bridge_result,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.pro20_origin import (  # noqa: E402
    MEMBER_KEYS,
    ORIGIN_TIME_HEX,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (  # noqa: E402
    TDG6_COMPLETE_STATE_CHANNELS,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_runtime import (  # noqa: E402
    TDG6TemporalLedger,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_ledger import seed_imp1_ledger  # noqa: E402


CORE = ROOT / "src/recursive_horizons/fgc/evolution/hlt17_srcq1_rec1.py"
CLI = ROOT / "scripts/qualify_fgc_hlt17_srcq1_rec1.py"
REC1 = "recursive_horizons.fgc.evolution.hlt17_srcq1_rec1"
START = 23.0 / 16.0
TARGET = 3.0 / 2.0
CFL_MAXIMUM_HEX = "0x1.0000000000000p-3"
FROZEN_CAP_INPUTS = {
    "RK4-2049": ("0x1.0000000000000p-4", "0x1.33345e70b5188p+0"),
    "RK4-4097": ("0x1.0000000000000p-5", "0x1.3333333aecf3ep+0"),
    "RK4-8193": ("0x1.0000000000000p-6", "0x1.3333333333653p+0"),
    "SSPRK3-4097": ("0x1.0000000000000p-5", "0x1.3333333fe51c4p+0"),
    "SSPRK3-8193": ("0x1.0000000000000p-6", "0x1.33333333349e5p+0"),
    "SSPRK3-16385": ("0x1.0000000000000p-7", "0x1.3333333333346p+0"),
}
COUNTERS = {
    "RK4-2049": (342, 1710, 0, 217),
    "RK4-4097": (673, 3365, 0, 443),
    "RK4-8193": (976, 4880, 0, 161),
    "SSPRK3-4097": (631, 2524, 0, 358),
    "SSPRK3-8193": (987, 3948, 0, 187),
    "SSPRK3-16385": (1771, 7084, 0, 0),
}


def _digest(label: str) -> str:
    return sha256(label.encode("ascii")).hexdigest()


def _plan(cap: float):
    return plan_imp1_lattice(START, TARGET, cap)


def _frozen_plan(key: str):
    return _plan(float.fromhex(FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER[key]))


def _half_plan(plan):
    return plan_imp1_lattice(
        plan.current, plan.event_target, float(plan.macro_width) / 2.0
    )


def _diagnostics(*, raw_gate: bool = True) -> dict[str, object]:
    return {
        "source_residual_infinity": 0.0,
        "source_raw_gate_passed": raw_gate,
        "source_refinement_iterations": 1,
        "source_residual_decreased_monotonically": True,
        "kinetic_condition_infinity": 1.0,
        "reduction_constraint_infinity": 0.0,
        "coordinate_speed_upper": 1.0,
        "minimum_lapse": 1.0,
        "minimum_radial_metric": 1.0,
        "minimum_areal_radius_away_from_center": 0.125,
    }


def _state() -> EvolutionState:
    u = np.zeros((9, 6), dtype=np.float64)
    p = np.zeros_like(u)
    q = np.zeros_like(u)
    u[:, 0] = 1.0
    u[:, 2] = 1.0
    u[:, 3] = np.linspace(0.0, 1.0, 9)
    q[:, 3] = 1.0
    return EvolutionState(u, p, q)


@dataclass
class FakeMonitor:
    accepted_stage_count: int = 0

    def as_mapping(self) -> dict[str, object]:
        return {"accepted_stage_count": self.accepted_stage_count}


@dataclass
class FakeCausal:
    accepted_time: float = START
    previous_speed_upper: float = 1.2

    def as_mapping(self) -> dict[str, str]:
        return {
            "accepted_time_hex": float(self.accepted_time).hex(),
            "previous_speed_upper_hex": float(self.previous_speed_upper).hex(),
        }


@dataclass
class FakeTracers:
    positions: np.ndarray
    proper_times: np.ndarray


@dataclass
class FakeTransaction:
    state: FakeMonitor
    causal_state: FakeCausal
    grid_spacing: float = 0.5
    cfl_maximum: float = 0.125


@dataclass
class FakeIdentity:
    method_label: str
    point_count: int
    integrator_id: str
    spatial_order: int


class FakeBinding:
    def __init__(self, digest: str, *, mutate: bool = False, raw_gate: bool = True) -> None:
        self.configuration_sha256 = digest
        self.mutate = mutate
        self.raw_gate = raw_gate

    def validate(self) -> None:
        return None

    def rhs(self, time: float, state: EvolutionState) -> EvolutionRHS:
        del time
        if self.mutate:
            writable = np.array(state.u, copy=True)
            writable[0, 0] += 1.0
            state.u[:] = writable
        zeros = np.zeros_like(state.u)
        return EvolutionRHS(zeros, zeros, zeros, _diagnostics(raw_gate=self.raw_gate))


class FakeCaptured:
    def __init__(
        self, key: str, *, step: int, serial: int, source: int, cfl: int, speed: float
    ) -> None:
        self.member_key = key
        self.step_index = step
        self.transaction_serial = serial
        self.source_retry_count = source
        self.cfl_retry_count = cfl
        self.accepted_time = {"rational": "23/16", "binary64_hex": ORIGIN_TIME_HEX}
        self.original_cursor_bytes = canonical_json_bytes(
            {"mode": "FRESH_READY", "retry_successor_payload_or_none": None}
        )
        self.causal_bytes = canonical_json_bytes(
            {
                "accepted_time": START,
                "accumulated_characteristic_distance": 0.0,
                "previous_speed_upper": speed,
            }
        )


class FakeBridge:
    def __init__(
        self,
        disposition: str,
        *,
        evidence: dict[str, object] | None = None,
        attempted_plan: object | None = None,
        next_plan: object | None = None,
        prepared: object | None = None,
        cursor: object | None = None,
        accepted_state_advanced: bool = False,
        committed: object | None = None,
        sink_writes: int = 0,
    ) -> None:
        self.disposition = disposition
        self.evidence = evidence or {}
        self.attempted_plan = attempted_plan
        self.next_plan = next_plan
        self.prepared = object() if prepared is None and disposition.startswith("prepared") else prepared
        self.cursor = cursor or SimpleNamespace(
            source_current=0,
            source_total=0,
            cfl_current=0,
            cfl_total=0,
            overlay=None,
        )
        self.accepted_state_advanced = accepted_state_advanced
        self.committed = committed
        self.sink_writes = sink_writes
        self.replay = None


class FakeCheckpoint:
    def __init__(self, results: list[object] | None = None) -> None:
        self._results = list(results or [])
        self._index = 0
        self.live = True
        self.absorbed: list[object] = []
        self.requested_caps: list[str] = []
        self.source_retry_count = 0
        self.CFL_retry_count = 0

    def isolate(self) -> "FakeCheckpoint":
        cloned = FakeCheckpoint(self._results[self._index :])
        cloned.live = False
        cloned.source_retry_count = self.source_retry_count
        cloned.CFL_retry_count = self.CFL_retry_count
        return cloned

    def agree(self) -> None:
        return None

    def prepare(self, *, requested_cap: object) -> object:
        self.requested_caps.append(float(requested_cap).hex())
        if self._index >= len(self._results):
            raise AssertionError("unexpected extra prepare")
        result = self._results[self._index]
        self._index += 1
        return result

    def absorb_nonfine(self, result: object) -> object:
        if self.live:
            raise AssertionError("source/CFL overlay absorbed on the live checkpoint")
        self.absorbed.append(result)
        disposition = getattr(result, "disposition", "")
        if "cfl" in disposition:
            self.CFL_retry_count += 1
            result.cursor.cfl_current = self.CFL_retry_count
            result.cursor.cfl_total = int(getattr(result.cursor, "cfl_total", 0)) + 1
        if "source" in disposition:
            self.source_retry_count += 1
            result.cursor.source_current = self.source_retry_count
            result.cursor.source_total = int(getattr(result.cursor, "source_total", 0)) + 1
        return result


class FakeMember:
    def __init__(self, key: str, captured: FakeCaptured, binding: FakeBinding) -> None:
        label, points = key.split("-", 1)
        integrator = PRIMARY_METHOD if label == "RK4" else COMPARATOR_METHOD
        spacing_hex, speed_hex = FROZEN_CAP_INPUTS[key]
        self.identity = FakeIdentity(
            method_label=label,
            point_count=int(points),
            integrator_id=integrator,
            spatial_order=4 if label == "RK4" else 2,
        )
        self.state = _state()
        self.time = START
        self.step_index = captured.step_index
        self.transaction_serial = captured.transaction_serial
        self.source_retry_count = captured.source_retry_count
        self.CFL_retry_count = captured.cfl_retry_count
        self.transaction = FakeTransaction(
            FakeMonitor(),
            FakeCausal(previous_speed_upper=float.fromhex(speed_hex)),
            grid_spacing=float.fromhex(spacing_hex),
            cfl_maximum=float.fromhex(CFL_MAXIMUM_HEX),
        )
        self.tracers = FakeTracers(
            np.array([0.25, 0.75], dtype=np.float64),
            np.array([0.0, 0.0], dtype=np.float64),
        )
        self.temporal_ledger = seed_imp1_ledger(
            TDG6TemporalLedger.zero(initial_time=START),
            method=integrator,
            state_sha256=array_content_sha256(self.state.u, self.state.p, self.state.q),
            step_index=self.step_index,
            transaction_serial=self.transaction_serial,
            origin_receipt_sha256=_digest("origin-receipt"),
        )
        self.source_binding = binding


class FakeView:
    def __init__(
        self,
        key: str,
        captured: FakeCaptured,
        binding: FakeBinding,
        checkpoint: FakeCheckpoint,
    ) -> None:
        self.member_key = key
        self.source_binding = binding
        self.source_configuration_sha256 = binding.configuration_sha256
        self.member = FakeMember(key, captured, binding)
        self.checkpoint = checkpoint


class FakeOrigin:
    def __init__(self, captured: dict[str, FakeCaptured]) -> None:
        self.capture_sha256 = ORIGIN_CAPTURE_SHA256
        self.generation1_bridge_content_id = _digest("bridge")
        self.members = tuple(captured[key] for key in MEMBER_KEYS)


class Clock:
    def __init__(self, *, jump_after: int | None = None, jump: float = 1e9) -> None:
        self.t = 0.0
        self.n = 0
        self.jump_after = jump_after
        self.jump = jump

    def __call__(self) -> float:
        self.n += 1
        if self.jump_after is not None and self.n >= self.jump_after:
            self.t += self.jump
        else:
            self.t += 0.001
        return self.t


def _static_pin():
    pins = {path: sha256(b"x").hexdigest() for path in STATIC_INPUT_PATHS}
    return patch(
        "recursive_horizons.fgc.evolution.hlt17_srcq1.STATIC_INPUT_SHA256", pins
    )


def _eighteen(passed: bool = True, classification: str = "resolved_order_pass"):
    return [
        {
            "channel": name,
            "admission_passed": passed,
            "public_fine_debit": "0",
            "gate_debit": "0",
            "extra_enclosure_debit": "0",
            "corrected": {
                "decision": {
                    "classification": classification,
                    "admission_passed": passed,
                }
            },
        }
        for name in TDG6_COMPLETE_STATE_CHANNELS
    ]


def _passing_c1r1():
    return patch(
        f"{REC1}.c1r1_shadow_record",
        return_value={
            "eighteen_channel_decisions": _eighteen(),
            "admission_passed": True,
            "public_fine_debits": ["0"] * 18,
        },
    )


def _prepared(key: str = "SSPRK3-4097", disposition: str = "prepared_fresh") -> FakeBridge:
    plan = _frozen_plan(key)
    return FakeBridge(
        disposition,
        evidence={"kind": disposition, "requested_cap_hex": plan.requested_cap.hex()},
        attempted_plan=plan,
        next_plan=None,
        prepared=object(),
    )


def _cfl_evidence(plan, *, observed: float = 0.25, maximum: float = 0.125) -> dict[str, object]:
    return {
        "kind": "cfl_retry",
        "path": "cfl",
        "observed_ratio_hex": float(observed).hex(),
        "maximum_ratio_hex": float(maximum).hex(),
        "attempted_requested_cap_hex": plan.requested_cap.hex(),
        "attempted_macro_width_hex": plan.macro_width.hex(),
        "accepted_boundary_preserved": True,
        "physical_classification": False,
    }


def _cfl_bridge(
    key: str = "SSPRK3-4097",
    *,
    exhausted: bool = False,
    observed: float = 0.25,
    maximum: float = 0.125,
    plan=None,
    next_plan=None,
) -> FakeBridge:
    attempted = plan or _frozen_plan(key)
    successor = None if exhausted else (next_plan or _half_plan(attempted))
    overlay = SimpleNamespace(
        exhausted=exhausted,
        owner="cfl",
        next_plan=successor,
        attempted_plan=attempted,
        owner_current_count=1,
    )
    return FakeBridge(
        "cfl_retry_exhausted" if exhausted else "cfl_retry_required",
        evidence=_cfl_evidence(attempted, observed=observed, maximum=maximum),
        attempted_plan=attempted,
        next_plan=successor,
        prepared=None,
        cursor=SimpleNamespace(
            source_current=0,
            source_total=0,
            cfl_current=1,
            cfl_total=359,
            overlay=overlay,
        ),
    )


def _source_bridge(key: str = "SSPRK3-4097", *, exhausted: bool = False) -> FakeBridge:
    attempted = _frozen_plan(key)
    successor = None if exhausted else _half_plan(attempted)
    evidence = {
        "kind": "source_retry",
        "path": "source",
        "source_retry": {"wired": True},
        "attempted_requested_cap_hex": attempted.requested_cap.hex(),
        "attempted_macro_width_hex": attempted.macro_width.hex(),
        "accepted_boundary_preserved": True,
        "physical_classification": False,
    }
    overlay = SimpleNamespace(
        exhausted=exhausted,
        owner="source",
        next_plan=successor,
        attempted_plan=attempted,
        owner_current_count=1,
    )
    return FakeBridge(
        "source_retry_exhausted" if exhausted else "source_retry_required",
        evidence=evidence,
        attempted_plan=attempted,
        next_plan=successor,
        prepared=None,
        cursor=SimpleNamespace(
            source_current=1,
            source_total=1,
            cfl_current=0,
            cfl_total=358,
            overlay=overlay,
        ),
    )


def _attempt2_payloads() -> tuple[bytes, bytes]:
    raw = {
        "accepted_state_advanced": False,
        "artifact_id": "FGC-1-HLT17-SRCQ1",
        "campaign_execution_authorized": False,
        "endpoint_adopted": False,
        "halted_outcome": "cfl",
        "store_published": False,
        "members": [
            {
                "member_key": "RK4-2049",
                "outcome": "prepared",
                "prepare_disposition": "prepared_fresh",
                "c1r1": {"receipt_sha256": "a" * 64},
            },
            {
                "member_key": "RK4-4097",
                "outcome": "prepared",
                "prepare_disposition": "prepared_fresh",
                "c1r1": {"receipt_sha256": "b" * 64},
            },
            {
                "member_key": "RK4-8193",
                "outcome": "prepared",
                "prepare_disposition": "prepared_fresh",
                "c1r1": {"receipt_sha256": "c" * 64},
            },
            {
                "member_key": "SSPRK3-4097",
                "outcome": "cfl",
                "prepare_disposition": "cfl_retry_required",
                "c1r1": None,
            },
        ],
    }
    raw_bytes = json.dumps(raw, sort_keys=True, separators=(",", ":")).encode("utf-8")
    raw_digest = sha256(raw_bytes).hexdigest()
    compact = {
        "artifact_id": ATTEMPT2_ARTIFACT_ID,
        "authority_commit": ATTEMPT2_AUTHORITY_COMMIT,
        "classification": ATTEMPT2_CLASSIFICATION,
        "completed_member_count": 4,
        "implementation_commit": ATTEMPT2_IMPLEMENTATION_COMMIT,
        "implementation_sha256": ATTEMPT2_IMPLEMENTATION_SHA256,
        "member_outcomes": [
            {
                "member_key": "RK4-2049",
                "outcome": "prepared",
                "source_evaluation_recorded": True,
            },
            {
                "member_key": "RK4-4097",
                "outcome": "prepared",
                "source_evaluation_recorded": True,
            },
            {
                "member_key": "RK4-8193",
                "outcome": "prepared",
                "source_evaluation_recorded": True,
            },
            {
                "cfl_ratio_evidence_recorded": False,
                "member_key": "SSPRK3-4097",
                "outcome": "runner_reported_cfl",
                "source_evaluation_recorded": True,
            },
        ],
        "origin_capture_sha256": ORIGIN_CAPTURE_SHA256,
        "raw_terminal": {
            "directory_payload_sha256": raw_digest,
            "leaf_bytes": len(raw_bytes),
            "leaf_path": ATTEMPT2_RAW_LEAF_RELATIVE,
            "leaf_sha256": raw_digest,
        },
        "retry_authorized": False,
        "state_boundary": {
            "accepted_state_advanced": False,
            "campaign_store_written": False,
            "endpoint_adopted_or_serialized": False,
            "pro20_namespace_absent": True,
        },
    }
    compact_bytes = json.dumps(compact, indent=2, sort_keys=True).encode("utf-8")
    return compact_bytes, raw_bytes


def pin_attempt2(compact_bytes: bytes, raw_bytes: bytes):
    digest = sha256(compact_bytes).hexdigest()
    return patch.multiple(
        REC1,
        ATTEMPT2_COMPACT_SHA256=digest,
        ATTEMPT2_SEALED_COMPACT_SHA256=digest,
        ATTEMPT2_RAW_LEAF_SHA256=sha256(raw_bytes).hexdigest(),
        ATTEMPT2_RAW_DIRECTORY_PAYLOAD_SHA256=sha256(raw_bytes).hexdigest(),
        ATTEMPT2_RAW_LEAF_BYTES=len(raw_bytes),
    )


def _metadata_linked_compact_bytes(sealed: bytes) -> bytes:
    marker = (
        b'  "classification": "completed_terminal_reports_cfl_evidence_incomplete_no_state_advance",\n'
    )
    inserted = (
        marker + b'  "derivation_document": "docs/fgc-hlt17-srcq1-frz1.md",\n'
    )
    if marker not in sealed:
        raise AssertionError("sealed Attempt2 compact layout differs")
    return sealed.replace(marker, inserted, 1)


def _sealed_compact_bytes_from_live_checkout() -> bytes:
    live = (ROOT / ATTEMPT2_COMPACT_RELATIVE).read_bytes()
    digest = sha256(live).hexdigest()
    if digest == ATTEMPT2_SEALED_COMPACT_SHA256:
        return live
    if digest != ATTEMPT2_METADATA_LINKED_COMPACT_SHA256:
        raise AssertionError("live Attempt2 compact identity differs")
    marker = (
        b'  "derivation_document": "docs/fgc-hlt17-srcq1-frz1.md",\n'
    )
    if live.count(marker) != 1:
        raise AssertionError("metadata-linked Attempt2 compact layout differs")
    sealed = live.replace(marker, b"", 1)
    if sha256(sealed).hexdigest() != ATTEMPT2_SEALED_COMPACT_SHA256:
        raise AssertionError("derived sealed Attempt2 compact identity differs")
    return sealed


def _authority_receipt(
    output: Path,
    *,
    commit: str = "a" * 40,
    sha256_digest: str = "b" * 64,
    closure: str,
    **extra: object,
) -> dict[str, object]:
    receipt: dict[str, object] = {
        "authority_commit": commit,
        "implementation_sha256": sha256_digest,
        "source_closure_sha256": closure,
        "prospective_requested_cap_hex_by_member": dict(
            FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER
        ),
        "output_namespace": str(output),
    }
    receipt.update(extra)
    return receipt


def _world(scripts: dict[str, list[object]] | None = None):
    captured = {}
    views = {}
    configs = {}
    scripts = scripts or {}
    for key in MEMBER_KEYS:
        step, serial, source, cfl = COUNTERS[key]
        spacing_hex, speed_hex = FROZEN_CAP_INPUTS[key]
        item = FakeCaptured(
            key,
            step=step,
            serial=serial,
            source=source,
            cfl=cfl,
            speed=float.fromhex(speed_hex),
        )
        digest = _digest(f"config-{key}")
        binding = FakeBinding(digest)
        if key in PREDECESSOR_MEMBER_KEYS:
            checkpoint = FakeCheckpoint([])
        else:
            checkpoint = FakeCheckpoint(scripts.get(key, [_prepared(key)]))
        captured[key] = item
        views[key] = FakeView(key, item, binding, checkpoint)
        configs[key] = digest
    origin = FakeOrigin(captured)
    runtime = SimpleNamespace(
        members=views,
        captured_origin_sha256=origin.capture_sha256,
        source_closure_sha256=_digest("closure"),
        construction_sha256=_digest("construction"),
    )
    return origin, runtime, configs


def _qualify(origin, runtime, compact_bytes, raw_bytes, **kwargs):
    defaults = {
        "static_input_bytes": {path: b"x" for path in STATIC_INPUT_PATHS},
        "source_closure_sha256": runtime.source_closure_sha256,
        "environment": {"python_implementation": "CPython"},
        "attempt2_compact_bytes": compact_bytes,
        "attempt2_raw_bytes": raw_bytes,
        "runtime_origin": runtime,
        "runtime_origin_builder": lambda *a, **k: (_ for _ in ()).throw(
            AssertionError("factory")
        ),
        "rss_bytes": lambda: 1,
        "clock": Clock(),
        "expected_origin_sha256": ORIGIN_CAPTURE_SHA256,
    }
    defaults.update(kwargs)
    return qualify_srcq1_rec1(origin, **defaults)


def _load_cli():
    spec = importlib.util.spec_from_file_location("qualify_fgc_hlt17_srcq1_rec1_cli", CLI)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SRCQ1REC1ContractTests(unittest.TestCase):
    def test_plan_status_does_not_execute(self) -> None:
        status = plan_status()
        self.assertEqual(status["artifact_id"], ARTIFACT_ID)
        self.assertEqual(status["schema"], SCHEMA)
        self.assertIs(status["physical_source_qualification_executed"], False)
        self.assertIs(status["rk_members_remeasured"], False)
        self.assertIs(status["campaign_execution_authorized"], False)
        self.assertEqual(status["resume_member_keys"], list(RESUME_MEMBER_KEYS))
        self.assertEqual(
            status["frozen_prospective_requested_cap_hex_by_member"]["SSPRK3-4097"],
            FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER["SSPRK3-4097"],
        )
        self.assertEqual(status["attempt2"]["compact_sha256"], ATTEMPT2_COMPACT_SHA256)
        self.assertEqual(
            status["attempt2"]["sealed_compact_sha256"], ATTEMPT2_SEALED_COMPACT_SHA256
        )
        self.assertEqual(
            status["attempt2"]["metadata_linked_compact_sha256"],
            ATTEMPT2_METADATA_LINKED_COMPACT_SHA256,
        )
        self.assertIs(
            status["attempt2"]["compact_identities"]["sealed_historical_blob"][
                "derivation_document"
            ],
            None,
        )
        self.assertEqual(
            status["attempt2"]["compact_identities"]["metadata_linked_blob"][
                "derivation_document"
            ],
            ATTEMPT2_DERIVATION_DOCUMENT,
        )
        self.assertEqual(status["attempt2"]["raw_leaf_sha256"], ATTEMPT2_RAW_LEAF_SHA256)
        self.assertEqual(status["authority_delta_validation"], AUTHORITY_SEAM)
        self.assertFalse(status["resource_ceilings"]["scientific_threshold_fitted"])

    def test_module_has_no_old_runner_import(self) -> None:
        tree = ast.parse(CORE.read_text("utf-8"))
        imported: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imported.append(node.module or "")
                imported.extend(alias.name for alias in node.names)
        joined = " ".join(imported)
        self.assertNotIn("run_fgc_gr0_calibration", joined)
        self.assertNotIn("hlt17_srcq1_auth1", joined)
        self.assertNotIn("hlt17_campaign_store", joined)
        cli_tree = ast.parse(CLI.read_text("utf-8"))
        cli_imported: list[str] = []
        for node in ast.walk(cli_tree):
            if isinstance(node, ast.Import):
                cli_imported.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                cli_imported.append(node.module or "")
                cli_imported.extend(alias.name for alias in node.names)
        cli_joined = " ".join(cli_imported)
        self.assertNotIn("run_fgc_gr0_calibration", cli_joined)
        self.assertNotIn("hlt17_srcq1_auth1", cli_joined)
        self.assertNotIn("hlt17_campaign_store", cli_joined)

    def test_default_cli_is_status_only(self) -> None:
        cli = _load_cli()
        with patch.object(
            cli, "capture_pro20_historical_origin", side_effect=AssertionError("store")
        ), patch.object(
            cli, "observe_environment", side_effect=AssertionError("env")
        ), patch.object(
            cli, "qualify_and_publish", side_effect=AssertionError("qualify")
        ), patch.object(
            cli, "read_regular_file", side_effect=AssertionError("io")
        ):
            code = cli.main([])
        self.assertEqual(code, 0)

    def test_cli_status_with_authority_commit_does_not_build_source(self) -> None:
        cli = _load_cli()
        with patch.object(
            cli, "capture_pro20_historical_origin", side_effect=AssertionError("store")
        ), patch.object(
            cli, "observe_environment", side_effect=AssertionError("env")
        ), patch.object(
            cli, "qualify_and_publish", side_effect=AssertionError("qualify")
        ), patch.object(
            cli, "read_regular_file", side_effect=AssertionError("io")
        ), patch.object(
            cli, "read_static_factory_inputs", side_effect=AssertionError("factory")
        ):
            self.assertEqual(cli.main(["--status"]), 0)
            self.assertEqual(cli.main(["--status", "--authority-commit", "a" * 40]), 1)

    def test_cli_refuses_partial_real_run_args(self) -> None:
        cli = _load_cli()
        with patch.object(
            cli, "capture_pro20_historical_origin", side_effect=AssertionError("store")
        ):
            self.assertEqual(cli.main(["--authority-commit", "a" * 40]), 2)
            self.assertEqual(cli.main(["--output-directory", "/tmp/srcq1-rec1-out"]), 2)

    def test_cli_real_path_stops_at_authority_seam_before_store(self) -> None:
        cli = _load_cli()
        with TemporaryDirectory() as raw:
            target = Path(raw) / "ns"
            with patch.object(
                cli,
                "capture_pro20_historical_origin",
                side_effect=AssertionError("store"),
            ), patch.object(
                cli, "observe_environment", side_effect=AssertionError("env")
            ), patch.object(
                cli, "read_regular_file", side_effect=AssertionError("io")
            ):
                code = cli.main(
                    ["--authority-commit", "a" * 40, "--output-directory", str(target)]
                )
        self.assertEqual(code, 1)

    def test_authority_seam_refuses_boolean_and_missing_module(self) -> None:
        with patch(
            f"{REC1}.importlib.import_module",
            side_effect=ImportError("missing freeze module"),
        ), self.assertRaises(HLT17SRCQ1REC1AuthoritySeamError):
            load_coordinator_authority_delta_validator()
        with self.assertRaises(HLT17SRCQ1REC1AuthoritySeamError):
            require_authority_delta(
                repository_root=ROOT,
                authority_commit="a" * 40,
                implementation_sha256="b" * 64,
                validator=lambda **kwargs: True,
            )
        receipt = require_authority_delta(
            repository_root=ROOT,
            authority_commit="a" * 40,
            implementation_sha256="b" * 64,
            validator=lambda **kwargs: {
                "delta": "exact",
                **kwargs,
                "source_closure_sha256": "c" * 64,
            },
        )
        self.assertIs(receipt["boolean_flag_accepted"], False)
        with self.assertRaises(HLT17SRCQ1REC1AuthoritySeamError):
            require_authority_delta(
                repository_root=ROOT,
                authority_commit="a" * 40,
                implementation_sha256="b" * 64,
                validator=lambda **kwargs: {
                    "authority_commit": "d" * 40,
                    "implementation_sha256": kwargs["implementation_sha256"],
                    "source_closure_sha256": "c" * 64,
                },
            )
        with self.assertRaises(HLT17SRCQ1REC1AuthoritySeamError):
            require_authority_delta(
                repository_root=ROOT,
                authority_commit="a" * 40,
                implementation_sha256="b" * 64,
                validator=lambda **kwargs: {
                    "authority_commit": kwargs["authority_commit"],
                    "implementation_sha256": kwargs["implementation_sha256"],
                },
            )

    def test_attempt2_real_compact_and_state_flags(self) -> None:
        compact_bytes = (ROOT / ATTEMPT2_COMPACT_RELATIVE).read_bytes()
        digest = sha256(compact_bytes).hexdigest()
        self.assertIn(
            digest,
            {
                ATTEMPT2_SEALED_COMPACT_SHA256,
                ATTEMPT2_METADATA_LINKED_COMPACT_SHA256,
            },
        )
        record = authenticate_attempt2_predecessor(
            compact_bytes=compact_bytes, require_raw=False
        )
        expected_identity = (
            "sealed_historical_blob"
            if digest == ATTEMPT2_SEALED_COMPACT_SHA256
            else "metadata_linked_blob"
        )
        self.assertEqual(record["compact_identity"], expected_identity)
        self.assertEqual(record["compact_sha256"], digest)
        self.assertEqual(record["classification"], ATTEMPT2_CLASSIFICATION)
        self.assertIs(record["accepted_state_advanced"], False)
        self.assertIs(record["campaign_store_written"], False)
        self.assertIs(record["endpoint_adopted_or_serialized"], False)
        self.assertIs(record["cfl_terminal_independently_reproduced"], False)
        self.assertEqual(record["raw_leaf_sha256"], ATTEMPT2_RAW_LEAF_SHA256)
        self.assertEqual(
            record["raw_directory_payload_sha256"], ATTEMPT2_RAW_DIRECTORY_PAYLOAD_SHA256
        )
        self.assertEqual(record["raw_leaf_bytes"], ATTEMPT2_RAW_LEAF_BYTES)
        refs = predecessor_rk_evidence_references(record)
        self.assertEqual([item["member_key"] for item in refs], list(PREDECESSOR_MEMBER_KEYS))
        self.assertTrue(all(item["remeasured"] is False for item in refs))

    def test_attempt2_raw_and_compact_hash_validation(self) -> None:
        compact_bytes, raw_bytes = _attempt2_payloads()
        with pin_attempt2(compact_bytes, raw_bytes):
            record = authenticate_attempt2_predecessor(
                compact_bytes=compact_bytes, raw_bytes=raw_bytes
            )
        self.assertEqual(record["compact_sha256"], sha256(compact_bytes).hexdigest())
        self.assertEqual(record["raw_leaf_sha256"], sha256(raw_bytes).hexdigest())
        with self.assertRaises(HLT17SRCQ1REC1IdentityError):
            authenticate_attempt2_predecessor(
                compact_bytes=b'{"artifact_id":"nope"}', require_raw=False
            )
        with pin_attempt2(compact_bytes, raw_bytes), self.assertRaises(
            HLT17SRCQ1REC1IdentityError
        ):
            authenticate_attempt2_predecessor(
                compact_bytes=compact_bytes, raw_bytes=b"not-the-raw-terminal"
            )

    def test_independent_cfl_ratio_predicate(self) -> None:
        checked = require_independent_cfl_ratio(
            {
                "kind": "cfl_retry",
                "observed_ratio_hex": (0.25).hex(),
                "maximum_ratio_hex": (0.125).hex(),
            }
        )
        self.assertIs(checked["observed_exceeds_maximum"], True)
        self.assertIs(checked["runner_label_sufficient"], False)
        with self.assertRaises(HLT17SRCQ1REC1InconclusiveError):
            require_independent_cfl_ratio({"kind": "cfl_retry"})
        with self.assertRaises(HLT17SRCQ1REC1Error):
            require_independent_cfl_ratio(
                {
                    "kind": "cfl_retry",
                    "observed_ratio_hex": (0.1).hex(),
                    "maximum_ratio_hex": (0.125).hex(),
                }
            )
        with self.assertRaises(HLT17SRCQ1REC1Error):
            require_independent_cfl_ratio(
                {
                    "kind": "source_retry",
                    "observed_ratio_hex": (0.25).hex(),
                    "maximum_ratio_hex": (0.125).hex(),
                }
            )
        with self.assertRaises(HLT17SRCQ1REC1Error):
            require_independent_cfl_ratio(
                {
                    "kind": "cfl_retry",
                    "observed_ratio_hex": (0.0).hex(),
                    "maximum_ratio_hex": (0.125).hex(),
                }
            )

    def test_rk_members_are_skipped_and_ssp_source_is_evaluated_once(self) -> None:
        compact_bytes, raw_bytes = _attempt2_payloads()
        origin, runtime, configs = _world()
        called: list[str] = []

        def evaluate(view):
            called.append(view.member_key)
            return evaluate_bound_source_once(
                binding=view.source_binding,
                time=float(view.member.time),
                state=view.member.state,
                transaction=view.member.transaction,
                tracers=view.member.tracers,
            )

        with _static_pin(), _passing_c1r1(), pin_attempt2(compact_bytes, raw_bytes):
            result = _qualify(
                origin,
                runtime,
                compact_bytes,
                raw_bytes,
                source_evaluator=evaluate,
                expected_source_configurations=configs,
                expected_environment={"python_implementation": "CPython"},
                expected_implementation=hlt17_implementation_identity(),
            )
        self.assertEqual(called, list(RESUME_MEMBER_KEYS))
        self.assertIs(result["rk_members_remeasured"], False)
        self.assertIsNone(result["halted_outcome"])
        self.assertEqual(result["evaluated_member_keys"], list(RESUME_MEMBER_KEYS))
        for key in PREDECESSOR_MEMBER_KEYS:
            row = next(item for item in result["members"] if item["member_key"] == key)
            self.assertIs(row["remeasured"], False)
            self.assertEqual(row["role"], "predecessor_evidence")
        for key in RESUME_MEMBER_KEYS:
            row = next(item for item in result["members"] if item["member_key"] == key)
            self.assertEqual(row["outcome"], "prepared")
            self.assertEqual(row["source_evaluation"]["source_health"], "recorded")
            self.assertEqual(row["fingerprints_before"], row["fingerprints_after"])

    def test_first_ssp_cfl_requires_observed_max_evidence(self) -> None:
        compact_bytes, raw_bytes = _attempt2_payloads()
        first = _frozen_plan("SSPRK3-4097")
        first_cfl = _cfl_bridge("SSPRK3-4097", plan=first)
        first_cfl.next_plan = None
        first_cfl.cursor.overlay.next_plan = None
        first_cfl.cursor.overlay.exhausted = True
        origin, runtime, _configs = _world(
            {
                "SSPRK3-4097": [first_cfl],
            }
        )
        with _static_pin(), _passing_c1r1(), pin_attempt2(compact_bytes, raw_bytes):
            result = _qualify(origin, runtime, compact_bytes, raw_bytes)
        ssp = next(item for item in result["members"] if item["member_key"] == "SSPRK3-4097")
        self.assertEqual(result["halted_outcome"], "cfl")
        self.assertEqual(ssp["outcome"], "cfl")
        self.assertEqual(ssp["attempts"][0]["bridge"]["disposition"], "cfl_retry_required")
        independent = ssp["attempts"][0]["bridge"]["independent_cfl"]
        self.assertIs(independent["observed_exceeds_maximum"], True)
        self.assertIs(ssp["attempts"][0]["overlay_absorbed"], True)
        self.assertEqual(result["evaluated_member_keys"], ["SSPRK3-4097"])
        working = runtime.members["SSPRK3-4097"].checkpoint
        self.assertEqual(working.absorbed, [])
        self.assertIs(working.live, True)

        origin, runtime, _configs = _world(
            {
                "SSPRK3-4097": [
                    FakeBridge(
                        "cfl_retry_required",
                        evidence={"kind": "cfl_retry"},
                        attempted_plan=first,
                        next_plan=_half_plan(first),
                        prepared=None,
                    )
                ]
            }
        )
        with _static_pin(), _passing_c1r1(), pin_attempt2(compact_bytes, raw_bytes):
            missing = _qualify(origin, runtime, compact_bytes, raw_bytes)
        self.assertEqual(missing["halted_outcome"], "inconclusive")
        ssp_missing = next(
            item for item in missing["newly_owned_members"] if item["member_key"] == "SSPRK3-4097"
        )
        self.assertFalse(ssp_missing.get("attempts", [{}])[0].get("overlay_absorbed", False))

    def test_overlay_absorption_follows_next_plan_and_then_prepares(self) -> None:
        compact_bytes, raw_bytes = _attempt2_payloads()
        first = _frozen_plan("SSPRK3-4097")
        second = _half_plan(first)
        third = _half_plan(second)
        cfl1 = _cfl_bridge("SSPRK3-4097", plan=first, next_plan=second)
        cfl2 = _cfl_bridge("SSPRK3-4097", plan=second, next_plan=third)
        prepared = FakeBridge(
            "prepared_overlay",
            evidence={"kind": "prepared_overlay"},
            attempted_plan=third,
            next_plan=None,
            prepared=object(),
        )
        origin, runtime, _configs = _world(
            {
                "SSPRK3-4097": [cfl1, cfl2, prepared],
                "SSPRK3-8193": [_prepared("SSPRK3-8193")],
                "SSPRK3-16385": [_prepared("SSPRK3-16385")],
            }
        )
        with _static_pin(), _passing_c1r1(), pin_attempt2(compact_bytes, raw_bytes):
            result = _qualify(origin, runtime, compact_bytes, raw_bytes)
        ssp = next(item for item in result["members"] if item["member_key"] == "SSPRK3-4097")
        self.assertEqual(ssp["outcome"], "prepared")
        self.assertEqual(ssp["prepare_count"], 3)
        self.assertEqual(ssp["overlay_count"], 2)
        self.assertEqual(
            ssp["attempts"][1]["bridge"]["attempted_plan"]["requested_cap_hex"],
            second.requested_cap.hex(),
        )
        self.assertIsNone(result["halted_outcome"])
        self.assertEqual(result["evaluated_member_keys"], list(RESUME_MEMBER_KEYS))
        isolated_caps = runtime.members["SSPRK3-4097"].checkpoint.requested_caps
        self.assertEqual(isolated_caps, [])

    def test_source_overlay_and_retry_exhaustion(self) -> None:
        compact_bytes, raw_bytes = _attempt2_payloads()
        first = _frozen_plan("SSPRK3-4097")
        nxt = _half_plan(first)
        source = _source_bridge("SSPRK3-4097")
        prepared = FakeBridge(
            "prepared_overlay",
            evidence={"kind": "prepared_overlay"},
            attempted_plan=nxt,
            prepared=object(),
        )
        origin, runtime, _configs = _world({"SSPRK3-4097": [source, prepared]})
        with _static_pin(), _passing_c1r1(), pin_attempt2(compact_bytes, raw_bytes):
            recovered = _qualify(origin, runtime, compact_bytes, raw_bytes)
        ssp = next(
            item for item in recovered["members"] if item["member_key"] == "SSPRK3-4097"
        )
        self.assertEqual(ssp["outcome"], "prepared")
        self.assertEqual(ssp["attempts"][0]["bridge"]["disposition"], "source_retry_required")
        self.assertIn("independent_source", ssp["attempts"][0]["bridge"])

        origin, runtime, _configs = _world(
            {"SSPRK3-4097": [_cfl_bridge("SSPRK3-4097", exhausted=True)]}
        )
        with _static_pin(), _passing_c1r1(), pin_attempt2(compact_bytes, raw_bytes):
            exhausted = _qualify(origin, runtime, compact_bytes, raw_bytes)
        self.assertEqual(exhausted["halted_outcome"], "cfl")
        ssp_ex = next(
            item for item in exhausted["members"] if item["member_key"] == "SSPRK3-4097"
        )
        self.assertEqual(ssp_ex["prepare_disposition"], "cfl_retry_exhausted")
        self.assertIs(ssp_ex["attempts"][0]["overlay_absorbed"], True)

        origin, runtime, _configs = _world(
            {"SSPRK3-4097": [_source_bridge("SSPRK3-4097", exhausted=True)]}
        )
        with _static_pin(), _passing_c1r1(), pin_attempt2(compact_bytes, raw_bytes):
            source_ex = _qualify(origin, runtime, compact_bytes, raw_bytes)
        self.assertEqual(source_ex["halted_outcome"], "source")

    def test_counters_fingerprints_and_no_endpoints(self) -> None:
        compact_bytes, raw_bytes = _attempt2_payloads()
        first = _frozen_plan("SSPRK3-4097")
        nxt = _half_plan(first)
        origin, runtime, _configs = _world(
            {
                "SSPRK3-4097": [
                    _cfl_bridge("SSPRK3-4097", plan=first, next_plan=nxt),
                    FakeBridge(
                        "prepared_overlay",
                        evidence={"kind": "prepared_overlay"},
                        attempted_plan=nxt,
                        prepared=object(),
                    ),
                ]
            }
        )
        live = runtime.members["SSPRK3-4097"].member
        before_source = live.source_retry_count
        before_cfl = live.CFL_retry_count
        with _static_pin(), _passing_c1r1(), pin_attempt2(compact_bytes, raw_bytes):
            result = _qualify(origin, runtime, compact_bytes, raw_bytes)
        ssp = next(item for item in result["members"] if item["member_key"] == "SSPRK3-4097")
        self.assertEqual(live.source_retry_count, before_source)
        self.assertEqual(live.CFL_retry_count, before_cfl)
        self.assertEqual(ssp["fingerprints_before"], ssp["fingerprints_after"])
        self.assertGreaterEqual(ssp["attempts"][0]["counters"]["working_CFL_retry_count"], 1)
        self.assertIs(result["accepted_state_advanced"], False)
        self.assertIs(result["endpoint_adopted"], False)
        self.assertIs(result["store_published"], False)
        refuse_endpoint_payload(result)
        with self.assertRaises(Exception):
            refuse_endpoint_payload({"corrected_endpoint": {"u": [1.0]}})
        with self.assertRaises(HLT17SRCQ1REC1Error):
            retain_bridge_result(
                FakeBridge("prepared_fresh", committed=object(), prepared=object())
            )

    def test_temporal_resource_invalid_and_isolation(self) -> None:
        compact_bytes, raw_bytes = _attempt2_payloads()
        origin, runtime, _configs = _world(
            {
                "SSPRK3-4097": [
                    FakeBridge(
                        "temporal_retry_required",
                        evidence={"kind": "temporal_retry_required"},
                        attempted_plan=_frozen_plan("SSPRK3-4097"),
                        next_plan=_half_plan(_frozen_plan("SSPRK3-4097")),
                        prepared=object(),
                        sink_writes=0,
                    )
                ]
            }
        )
        with _static_pin(), _passing_c1r1(), pin_attempt2(compact_bytes, raw_bytes):
            temporal = _qualify(origin, runtime, compact_bytes, raw_bytes)
        self.assertEqual(temporal["halted_outcome"], "temporal")

        origin, runtime, _configs = _world()
        ceilings = ResourceCeilings(
            max_rss_bytes=1,
            max_wall_seconds_per_member={
                key: 900.0 if key != "SSPRK3-4097" else 0.0000001 for key in MEMBER_KEYS
            },
            max_total_wall_seconds=3600.0,
            coordinator_freeze_required=True,
        )
        with _static_pin(), _passing_c1r1(), pin_attempt2(compact_bytes, raw_bytes):
            resource = _qualify(
                origin,
                runtime,
                compact_bytes,
                raw_bytes,
                resource_ceilings=ceilings,
                clock=Clock(jump_after=3, jump=10.0),
            )
        self.assertEqual(resource["halted_outcome"], "resource")
        live = FakeCheckpoint([_prepared()])
        with self.assertRaises(HLT17SRCQ1REC1Error):
            isolate_working_checkpoint(object())
        isolated = isolate_working_checkpoint(live)
        self.assertIsNot(isolated, live)
        self.assertIs(isolated.live, False)
        ceilings_default = recommended_resource_ceilings()
        self.assertEqual(ceilings_default.max_rss_bytes, 4 * 1024 * 1024 * 1024)
        self.assertEqual(ceilings_default.max_total_wall_seconds, 3600.0)

    def test_frozen_caps_refuse_retune(self) -> None:
        compact_bytes, raw_bytes = _attempt2_payloads()
        origin, runtime, _configs = _world()
        swapped = dict(FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER)
        swapped["SSPRK3-4097"] = swapped["SSPRK3-8193"]
        with _static_pin(), pin_attempt2(compact_bytes, raw_bytes), self.assertRaises(
            HLT17SRCQ1REC1CapError
        ):
            _qualify(
                origin,
                runtime,
                compact_bytes,
                raw_bytes,
                prospective_requested_cap_hex_by_member=swapped,
            )
        caps, facts = derive_prospective_requested_caps(origin, runtime)
        for key in MEMBER_KEYS:
            self.assertEqual(
                caps[key].hex(), FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER[key]
            )
            self.assertEqual(
                facts[key]["formula"],
                "cfl_maximum*grid_spacing/inherited_previous_speed_upper",
            )

    def test_atomic_result_existing_absent_and_uncertain(self) -> None:
        result = {
            "schema": SCHEMA,
            "artifact_id": ARTIFACT_ID,
            "accepted_state_advanced": False,
            "endpoint_adopted": False,
        }
        with TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            existing = root / "ns"
            existing.mkdir()
            with self.assertRaises(HLT17SRCQ1REC1PublicationError) as failed:
                publish_qualification_result(existing, result)
            self.assertFalse(failed.exception.published)
            missing_parent = root / "missing" / "ns"
            with self.assertRaises(HLT17SRCQ1REC1PublicationError) as absent:
                publish_qualification_result(missing_parent, result)
            self.assertFalse(absent.exception.published)
            target = root / "fresh"
            receipt = publish_qualification_result(target, result)
            self.assertTrue(receipt["published"])
            self.assertTrue((target / "qualification.json").is_file())

            def fail_after(cut: str) -> None:
                if cut == "after_publish":
                    raise OSError("injected durability failure")

            uncertain = root / "uncertain"
            with self.assertRaises(HLT17SRCQ1REC1PublicationError) as posted:
                publish_qualification_result(uncertain, result, fault_hook=fail_after)
            self.assertTrue(posted.exception.published)
            self.assertEqual(posted.exception.outcome, "inconclusive")

            malformed = dict(result)
            malformed["physical_source_qualification_attempted"] = 1
            with self.assertRaisesRegex(
                HLT17SRCQ1REC1Error,
                "physical_source_qualification_attempted must be a built-in bool",
            ):
                publish_qualification_result(root / "malformed", malformed)

    def test_qualify_and_publish_uses_injected_authority_not_a_flag(self) -> None:
        compact_bytes, raw_bytes = _attempt2_payloads()
        origin, runtime, configs = _world()
        with TemporaryDirectory() as raw, _static_pin(), _passing_c1r1(), pin_attempt2(
            compact_bytes, raw_bytes
        ), patch(
            f"{REC1}.implementation_identity",
            return_value={**implementation_identity(), "sha256": "b" * 64},
        ):
            output = Path(raw).resolve() / "ns"
            published = qualify_and_publish(
                origin,
                repository_root=ROOT,
                output_directory=output,
                authority_commit="a" * 40,
                static_input_bytes={path: b"x" for path in STATIC_INPUT_PATHS},
                source_closure_sha256=runtime.source_closure_sha256,
                environment={"python_implementation": "CPython"},
                attempt2_compact_bytes=compact_bytes,
                attempt2_raw_bytes=raw_bytes,
                authority_delta_validator=lambda **kwargs: _authority_receipt(
                    output,
                    commit=kwargs["authority_commit"],
                    sha256_digest=kwargs["implementation_sha256"],
                    closure=runtime.source_closure_sha256,
                    delta="injected",
                ),
                runtime_origin=runtime,
                runtime_origin_builder=lambda *a, **k: (_ for _ in ()).throw(
                    AssertionError("factory")
                ),
                rss_bytes=lambda: 1,
                clock=Clock(),
                expected_source_configurations=configs,
            )
            self.assertEqual(published["authority_commit"], "a" * 40)
            self.assertTrue((output / "qualification.json").is_file())
            self.assertIs(published["rk_members_remeasured"], False)


class SRCQ1REC1AdversarialTests(unittest.TestCase):
    def test_attempt2_metadata_linked_compact_identity(self) -> None:
        sealed = _sealed_compact_bytes_from_live_checkout()
        self.assertEqual(sha256(sealed).hexdigest(), ATTEMPT2_SEALED_COMPACT_SHA256)
        linked = _metadata_linked_compact_bytes(sealed)
        self.assertEqual(sha256(linked).hexdigest(), ATTEMPT2_METADATA_LINKED_COMPACT_SHA256)
        record = authenticate_attempt2_predecessor(
            compact_bytes=linked, require_raw=False
        )
        self.assertEqual(record["compact_identity"], "metadata_linked_blob")
        self.assertEqual(record["derivation_document"], ATTEMPT2_DERIVATION_DOCUMENT)
        self.assertEqual(record["sealed_compact_sha256"], ATTEMPT2_SEALED_COMPACT_SHA256)
        self.assertEqual(
            record["metadata_linked_compact_sha256"],
            ATTEMPT2_METADATA_LINKED_COMPACT_SHA256,
        )
        sealed_record = authenticate_attempt2_predecessor(
            compact_bytes=sealed, require_raw=False
        )
        self.assertEqual(sealed_record["compact_identity"], "sealed_historical_blob")
        mutated_science = linked.replace(
            b'"completed_member_count": 4', b'"completed_member_count": 5', 1
        )
        with self.assertRaises(HLT17SRCQ1REC1IdentityError):
            authenticate_attempt2_predecessor(
                compact_bytes=mutated_science, require_raw=False
            )
        mutated_metadata = linked.replace(
            b"docs/fgc-hlt17-srcq1-frz1.md", b"docs/fgc-hlt17-srcq1-rec1-frz1.md", 1
        )
        with self.assertRaises(HLT17SRCQ1REC1IdentityError):
            authenticate_attempt2_predecessor(
                compact_bytes=mutated_metadata, require_raw=False
            )

    def test_mixed_independent_retry_owners_follow_next_plan(self) -> None:
        compact_bytes, raw_bytes = _attempt2_payloads()
        first = _frozen_plan("SSPRK3-4097")
        second = _half_plan(first)
        third = _half_plan(second)
        cfl = _cfl_bridge("SSPRK3-4097", plan=first, next_plan=second)
        source = FakeBridge(
            "source_retry_required",
            evidence={
                "kind": "source_retry",
                "path": "source",
                "source_retry": {"wired": True},
                "attempted_requested_cap_hex": second.requested_cap.hex(),
                "attempted_macro_width_hex": second.macro_width.hex(),
                "accepted_boundary_preserved": True,
                "physical_classification": False,
            },
            attempted_plan=second,
            next_plan=third,
            prepared=None,
            cursor=SimpleNamespace(
                source_current=1,
                source_total=1,
                cfl_current=1,
                cfl_total=359,
                overlay=SimpleNamespace(
                    exhausted=False,
                    owner="source",
                    next_plan=third,
                    attempted_plan=second,
                    owner_current_count=1,
                ),
            ),
        )
        prepared = FakeBridge(
            "prepared_overlay",
            evidence={"kind": "prepared_overlay"},
            attempted_plan=third,
            prepared=object(),
        )
        origin, runtime, _configs = _world({"SSPRK3-4097": [cfl, source, prepared]})
        with _static_pin(), _passing_c1r1(), pin_attempt2(compact_bytes, raw_bytes):
            result = _qualify(origin, runtime, compact_bytes, raw_bytes)
        ssp = next(item for item in result["members"] if item["member_key"] == "SSPRK3-4097")
        self.assertEqual(ssp["outcome"], "prepared")
        self.assertEqual(ssp["cfl_owner_count"], 1)
        self.assertEqual(ssp["source_owner_count"], 1)
        self.assertEqual(ssp["prepare_count"], 3)
        self.assertGreater(ssp["prepare_count"], 2)

    def test_count_33_owner_exhaustion_is_lawful(self) -> None:
        compact_bytes, raw_bytes = _attempt2_payloads()
        frozen = FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER["SSPRK3-4097"]
        plans = [{"requested_cap_hex": frozen, "macro_width_hex": frozen}]
        for index in range(33):
            plans.append(
                {
                    "requested_cap_hex": f"overlay-{index + 1}",
                    "macro_width_hex": f"overlay-{index + 1}",
                }
            )
        script: list[object] = []
        for index in range(33):
            exhausted = index == 32
            attempted = plans[index]
            nxt = None if exhausted else plans[index + 1]
            script.append(
                FakeBridge(
                    "cfl_retry_exhausted" if exhausted else "cfl_retry_required",
                    evidence={
                        "kind": "cfl_retry",
                        "path": "cfl",
                        "observed_ratio_hex": (0.25).hex(),
                        "maximum_ratio_hex": (0.125).hex(),
                        "attempted_requested_cap_hex": attempted["requested_cap_hex"],
                        "attempted_macro_width_hex": attempted["macro_width_hex"],
                        "accepted_boundary_preserved": True,
                        "physical_classification": False,
                    },
                    attempted_plan=attempted,
                    next_plan=nxt,
                    prepared=None,
                    cursor=SimpleNamespace(
                        source_current=0,
                        source_total=0,
                        cfl_current=index + 1,
                        cfl_total=358 + index + 1,
                        overlay=SimpleNamespace(
                            exhausted=exhausted,
                            owner="cfl",
                            next_plan=nxt,
                            attempted_plan=attempted,
                            owner_current_count=index + 1,
                        ),
                    ),
                )
            )
        origin, runtime, _configs = _world({"SSPRK3-4097": script})
        with _static_pin(), _passing_c1r1(), pin_attempt2(compact_bytes, raw_bytes):
            result = _qualify(origin, runtime, compact_bytes, raw_bytes)
        ssp = next(item for item in result["members"] if item["member_key"] == "SSPRK3-4097")
        self.assertEqual(ssp["halted_outcome"], "cfl")
        self.assertEqual(ssp["cfl_owner_count"], 33)
        self.assertEqual(ssp["prepare_count"], 33)
        self.assertEqual(ssp["attempts"][-1]["bridge"]["disposition"], "cfl_retry_exhausted")

    def test_nonzero_or_malformed_sink_writes_are_invalid(self) -> None:
        compact_bytes, raw_bytes = _attempt2_payloads()
        origin, runtime, _configs = _world(
            {
                "SSPRK3-4097": [
                    FakeBridge(
                        "prepared_fresh",
                        evidence={"kind": "prepared_fresh"},
                        attempted_plan=_frozen_plan("SSPRK3-4097"),
                        prepared=object(),
                        sink_writes=1,
                    )
                ]
            }
        )
        with _static_pin(), _passing_c1r1(), pin_attempt2(compact_bytes, raw_bytes):
            result = _qualify(origin, runtime, compact_bytes, raw_bytes)
        self.assertEqual(result["halted_outcome"], "invalid")
        origin, runtime, _configs = _world(
            {
                "SSPRK3-4097": [
                    FakeBridge(
                        "prepared_fresh",
                        evidence={"kind": "prepared_fresh"},
                        attempted_plan=_frozen_plan("SSPRK3-4097"),
                        prepared=object(),
                        sink_writes=True,  # type: ignore[arg-type]
                    )
                ]
            }
        )
        with _static_pin(), _passing_c1r1(), pin_attempt2(compact_bytes, raw_bytes):
            malformed = _qualify(origin, runtime, compact_bytes, raw_bytes)
        self.assertEqual(malformed["halted_outcome"], "invalid")

    def test_missing_source_closure_refuses_before_physical_work(self) -> None:
        compact_bytes, raw_bytes = _attempt2_payloads()
        origin, runtime, _configs = _world()
        with TemporaryDirectory() as raw, _static_pin(), pin_attempt2(
            compact_bytes, raw_bytes
        ), patch(
            f"{REC1}.implementation_identity",
            return_value={**implementation_identity(), "sha256": "b" * 64},
        ):
            output = Path(raw).resolve() / "ns"
            with self.assertRaises(HLT17SRCQ1REC1AuthoritySeamError) as raised:
                qualify_and_publish(
                    origin,
                    repository_root=ROOT,
                    output_directory=output,
                    authority_commit="a" * 40,
                    static_input_bytes={path: b"x" for path in STATIC_INPUT_PATHS},
                    source_closure_sha256=runtime.source_closure_sha256,
                    environment={"python_implementation": "CPython"},
                    attempt2_compact_bytes=compact_bytes,
                    attempt2_raw_bytes=raw_bytes,
                    authority_delta_validator=lambda **kwargs: {
                        "authority_commit": kwargs["authority_commit"],
                        "implementation_sha256": kwargs["implementation_sha256"],
                        "prospective_requested_cap_hex_by_member": dict(
                            FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER
                        ),
                        "output_namespace": str(output),
                    },
                    runtime_origin_builder=lambda *a, **k: (_ for _ in ()).throw(
                        AssertionError("physical factory")
                    ),
                )
            self.assertEqual(raised.exception.phase, "authority_seam")
            self.assertFalse(raised.exception.physical_work_started)
            self.assertFalse(output.exists())

    def test_early_typed_factory_failure_publishes_terminal(self) -> None:
        compact_bytes, raw_bytes = _attempt2_payloads()
        origin, runtime, _configs = _world()
        with TemporaryDirectory() as raw, _static_pin(), pin_attempt2(
            compact_bytes, raw_bytes
        ), patch(
            f"{REC1}.implementation_identity",
            return_value={**implementation_identity(), "sha256": "b" * 64},
        ):
            output = Path(raw).resolve() / "ns"
            published = qualify_and_publish(
                origin,
                repository_root=ROOT,
                output_directory=output,
                authority_commit="a" * 40,
                static_input_bytes={path: b"x" for path in STATIC_INPUT_PATHS},
                source_closure_sha256=runtime.source_closure_sha256,
                environment={"python_implementation": "CPython"},
                attempt2_compact_bytes=compact_bytes,
                attempt2_raw_bytes=raw_bytes,
                authority_delta_validator=lambda **kwargs: _authority_receipt(
                    output,
                    commit=kwargs["authority_commit"],
                    sha256_digest=kwargs["implementation_sha256"],
                    closure=runtime.source_closure_sha256,
                ),
                runtime_origin_builder=lambda *a, **k: (_ for _ in ()).throw(
                    Pro20SourceFactoryError("typed factory failure")
                ),
                rss_bytes=lambda: 1,
                clock=Clock(),
            )
            self.assertTrue((output / "qualification.json").is_file())
            self.assertEqual(published["halted_outcome"], "source")
            self.assertIs(published["physical_source_qualification_attempted"], True)
            self.assertEqual(published["phase"], "qualification")

    def test_origin_loader_failure_publishes_invalid_pre_physical_terminal(self) -> None:
        compact_bytes, raw_bytes = _attempt2_payloads()
        _origin, runtime, _configs = _world()
        with TemporaryDirectory() as raw, _static_pin(), pin_attempt2(
            compact_bytes, raw_bytes
        ), patch(
            f"{REC1}.implementation_identity",
            return_value={**implementation_identity(), "sha256": "b" * 64},
        ):
            output = Path(raw).resolve() / "ns"
            published = qualify_and_publish(
                None,
                repository_root=ROOT,
                output_directory=output,
                authority_commit="a" * 40,
                static_input_bytes={path: b"x" for path in STATIC_INPUT_PATHS},
                source_closure_sha256=runtime.source_closure_sha256,
                environment={"python_implementation": "CPython"},
                attempt2_compact_bytes=compact_bytes,
                attempt2_raw_bytes=raw_bytes,
                authority_delta_validator=lambda **kwargs: _authority_receipt(
                    output,
                    commit=kwargs["authority_commit"],
                    sha256_digest=kwargs["implementation_sha256"],
                    closure=runtime.source_closure_sha256,
                ),
                origin_loader=lambda: (_ for _ in ()).throw(
                    Pro20OriginError("typed origin capture failure")
                ),
                runtime_origin_builder=lambda *a, **k: (_ for _ in ()).throw(
                    AssertionError("physical factory")
                ),
                rss_bytes=lambda: 1,
                clock=Clock(),
            )
            self.assertTrue((output / "qualification.json").is_file())
            self.assertEqual(published["halted_outcome"], "invalid")
            self.assertIs(published["physical_source_qualification_attempted"], False)
            self.assertEqual(published["phase"], "repository_input")

    def test_pre_physical_publication_uncertainty_keeps_attempted_false(self) -> None:
        compact_bytes, raw_bytes = _attempt2_payloads()
        _origin, runtime, _configs = _world()
        with TemporaryDirectory() as raw, _static_pin(), pin_attempt2(
            compact_bytes, raw_bytes
        ), patch(
            f"{REC1}.implementation_identity",
            return_value={**implementation_identity(), "sha256": "b" * 64},
        ):
            output = Path(raw).resolve() / "ns"

            def fail_after(cut: str) -> None:
                if cut == "after_publish":
                    raise OSError("injected durability failure")

            with self.assertRaises(HLT17SRCQ1REC1PublicationError) as posted:
                qualify_and_publish(
                    None,
                    repository_root=ROOT,
                    output_directory=output,
                    authority_commit="a" * 40,
                    static_input_bytes={path: b"x" for path in STATIC_INPUT_PATHS},
                    source_closure_sha256=runtime.source_closure_sha256,
                    environment={"python_implementation": "CPython"},
                    attempt2_compact_bytes=compact_bytes,
                    attempt2_raw_bytes=raw_bytes,
                    authority_delta_validator=lambda **kwargs: _authority_receipt(
                        output,
                        commit=kwargs["authority_commit"],
                        sha256_digest=kwargs["implementation_sha256"],
                        closure=runtime.source_closure_sha256,
                    ),
                    origin_loader=lambda: (_ for _ in ()).throw(
                        Pro20OriginError("typed origin capture failure")
                    ),
                    runtime_origin_builder=lambda *a, **k: (_ for _ in ()).throw(
                        AssertionError("physical factory")
                    ),
                    rss_bytes=lambda: 1,
                    clock=Clock(),
                    fault_hook=fail_after,
                )
            self.assertTrue(posted.exception.published)
            self.assertFalse(posted.exception.physical_work_started)
            self.assertEqual(posted.exception.phase, "publication")
            self.assertEqual(posted.exception.outcome, "inconclusive")

    def test_programmer_exception_is_not_classified_as_science(self) -> None:
        compact_bytes, raw_bytes = _attempt2_payloads()
        origin, runtime, _configs = _world()
        with TemporaryDirectory() as raw, _static_pin(), pin_attempt2(
            compact_bytes, raw_bytes
        ), patch(
            f"{REC1}.implementation_identity",
            return_value={**implementation_identity(), "sha256": "b" * 64},
        ):
            output = Path(raw).resolve() / "ns"
            with self.assertRaises(AssertionError):
                qualify_and_publish(
                    origin,
                    repository_root=ROOT,
                    output_directory=output,
                    authority_commit="a" * 40,
                    static_input_bytes={path: b"x" for path in STATIC_INPUT_PATHS},
                    source_closure_sha256=runtime.source_closure_sha256,
                    environment={"python_implementation": "CPython"},
                    attempt2_compact_bytes=compact_bytes,
                    attempt2_raw_bytes=raw_bytes,
                    authority_delta_validator=lambda **kwargs: _authority_receipt(
                        output,
                        commit=kwargs["authority_commit"],
                        sha256_digest=kwargs["implementation_sha256"],
                        closure=runtime.source_closure_sha256,
                    ),
                    runtime_origin_builder=lambda *a, **k: (_ for _ in ()).throw(
                        AssertionError("programmer bug")
                    ),
                    rss_bytes=lambda: 1,
                    clock=Clock(),
                )
            self.assertFalse(output.exists())

    def test_publication_uncertainty_does_not_claim_unexecuted(self) -> None:
        compact_bytes, raw_bytes = _attempt2_payloads()
        origin, runtime, configs = _world()
        with TemporaryDirectory() as raw, _static_pin(), _passing_c1r1(), pin_attempt2(
            compact_bytes, raw_bytes
        ), patch(
            f"{REC1}.implementation_identity",
            return_value={**implementation_identity(), "sha256": "b" * 64},
        ):
            output = Path(raw).resolve() / "ns"

            def fail_after(cut: str) -> None:
                if cut == "after_publish":
                    raise OSError("injected durability failure")

            with self.assertRaises(HLT17SRCQ1REC1PublicationError) as posted:
                qualify_and_publish(
                    origin,
                    repository_root=ROOT,
                    output_directory=output,
                    authority_commit="a" * 40,
                    static_input_bytes={path: b"x" for path in STATIC_INPUT_PATHS},
                    source_closure_sha256=runtime.source_closure_sha256,
                    environment={"python_implementation": "CPython"},
                    attempt2_compact_bytes=compact_bytes,
                    attempt2_raw_bytes=raw_bytes,
                    authority_delta_validator=lambda **kwargs: _authority_receipt(
                        output,
                        commit=kwargs["authority_commit"],
                        sha256_digest=kwargs["implementation_sha256"],
                        closure=runtime.source_closure_sha256,
                    ),
                    runtime_origin=runtime,
                    runtime_origin_builder=lambda *a, **k: (_ for _ in ()).throw(
                        AssertionError("factory")
                    ),
                    rss_bytes=lambda: 1,
                    clock=Clock(),
                    expected_source_configurations=configs,
                    fault_hook=fail_after,
                )
            self.assertTrue(posted.exception.published)
            self.assertTrue(posted.exception.physical_work_started)
            self.assertEqual(posted.exception.phase, "publication")
            self.assertEqual(posted.exception.outcome, "inconclusive")

    def test_cfl_evidence_failure_retains_attempt_bridge(self) -> None:
        compact_bytes, raw_bytes = _attempt2_payloads()
        first = _frozen_plan("SSPRK3-4097")
        origin, runtime, _configs = _world(
            {
                "SSPRK3-4097": [
                    FakeBridge(
                        "cfl_retry_required",
                        evidence={"kind": "cfl_retry"},
                        attempted_plan=first,
                        next_plan=_half_plan(first),
                        prepared=None,
                    )
                ]
            }
        )
        with _static_pin(), _passing_c1r1(), pin_attempt2(compact_bytes, raw_bytes):
            result = _qualify(origin, runtime, compact_bytes, raw_bytes)
        ssp = next(
            item for item in result["newly_owned_members"] if item["member_key"] == "SSPRK3-4097"
        )
        self.assertEqual(ssp["outcome"], "inconclusive")
        self.assertEqual(len(ssp["attempts"]), 1)
        self.assertEqual(ssp["attempts"][0]["bridge"]["disposition"], "cfl_retry_required")
        self.assertIn("evidence", ssp["attempts"][0]["bridge"])


if __name__ == "__main__":
    unittest.main()
