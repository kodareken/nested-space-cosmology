"""Synthetic/adversarial PREF1 tests. No physical source or live store."""

from __future__ import annotations

import ast
from copy import deepcopy
from dataclasses import dataclass
from hashlib import sha256
import importlib.util
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from typing import Mapping
import unittest
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from recursive_horizons.fgc.evolution.hlt17_imp1_cursor import (  # noqa: E402
    plan_imp1_lattice,
)
from recursive_horizons.fgc.evolution.hlt17_srcq1 import (  # noqa: E402
    ResourceCeilings,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.pro20_origin import MEMBER_KEYS  # noqa: E402
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (  # noqa: E402
    TDG6_COMPLETE_STATE_CHANNELS,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_runtime import (  # noqa: E402
    TDG6TemporalLedger,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_ledger import (  # noqa: E402
    seed_imp1_ledger,
)
from recursive_horizons.fgc.evolution import (  # noqa: E402
    hlt17_srcq1_rec1_pref1_binder as binder,
)


CORE = ROOT / "src/recursive_horizons/fgc/evolution/hlt17_srcq1_rec1_pref1_binder.py"
SCRIPT = ROOT / "scripts/reproduce_fgc_hlt17_srcq1_rec1_pref1.py"
PREF1 = "recursive_horizons.fgc.evolution.hlt17_srcq1_rec1_pref1_binder"
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
    return _plan(float.fromhex(binder.FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER[key]))


def _half_plan(plan):
    return plan_imp1_lattice(plan.current, plan.event_target, float(plan.macro_width) / 2.0)


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
    def __init__(self, digest: str) -> None:
        self.configuration_sha256 = digest

    def validate(self) -> None:
        return None

    def rhs(self, time: float, state: EvolutionState) -> EvolutionRHS:
        del time
        zeros = np.zeros_like(state.u)
        return EvolutionRHS(zeros, zeros, zeros, _diagnostics())


class FakeCaptured:
    def __init__(
        self, key: str, *, step: int, serial: int, source: int, cfl: int, speed: float
    ) -> None:
        self.member_key = key
        self.step_index = step
        self.transaction_serial = serial
        self.source_retry_count = source
        self.cfl_retry_count = cfl
        self.accepted_time = {"rational": "23/16", "binary64_hex": "0x1.7000000000000p+0"}
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
        self.prepared = (
            object()
            if prepared is None and disposition.startswith("prepared")
            else prepared
        )
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
        self.capture_sha256 = binder.ORIGIN_CAPTURE_SHA256
        self.members = tuple(captured[key] for key in MEMBER_KEYS)


class FakeRuntime:
    def __init__(self, views: dict[str, FakeView]) -> None:
        self.members = {key: views[key] for key in MEMBER_KEYS}
        self.captured_origin_sha256 = binder.ORIGIN_CAPTURE_SHA256
        self.source_closure_sha256 = binder.SOURCE_CLOSURE_SHA256


def _eighteen() -> list[dict[str, object]]:
    return [
        {
            "channel": name,
            "admission_passed": True,
            "classification": "resolved_order_pass",
            "public_fine_debit": "0/1",
            "gate_debit": "0/1",
            "extra_enclosure_debit": "0/1",
        }
        for name in TDG6_COMPLETE_STATE_CHANNELS
    ]


def _c1r1(key: str) -> dict[str, object]:
    return {
        "eighteen_channel_decisions": _eighteen(),
        "admission_passed": True,
        "public_fine_debits": ["0/1"] * 18,
        "receipt_sha256": binder.EXPECTED_C1R1_RECEIPTS[key],
        "assessment_sha256": binder.EXPECTED_C1R1_ASSESSMENTS[key],
        "family_sha256": binder.EXPECTED_C1R1_FAMILIES[key],
        "prepared": True,
        "endpoint_adopted": False,
        "campaign_execution_authorized": False,
        "corrected_endpoint_committable": False,
        "outer_or_medium_committable": False,
    }


def _cfl_evidence(plan, *, observed: str, maximum: str) -> dict[str, object]:
    return {
        "kind": "cfl_retry",
        "path": "outer_0",
        "observed_ratio_hex": observed,
        "maximum_ratio_hex": maximum,
        "attempted_requested_cap_hex": plan.requested_cap.hex(),
        "attempted_macro_width_hex": plan.macro_width.hex(),
        "accepted_boundary_preserved": True,
        "physical_classification": False,
    }


def _cfl_bridge(key: str = "SSPRK3-4097") -> FakeBridge:
    attempted = _frozen_plan(key)
    successor = plan_imp1_lattice(
        attempted.current,
        attempted.event_target,
        float.fromhex(binder.EXPECTED_SSP4097_NEXT_PLAN_CAP_HEX),
    )
    overlay = SimpleNamespace(
        exhausted=False,
        owner="cfl",
        next_plan=successor,
        attempted_plan=attempted,
        owner_current_count=1,
    )
    return FakeBridge(
        "cfl_retry_required",
        evidence=_cfl_evidence(
            attempted,
            observed=binder.EXPECTED_SSP4097_CFL["observed_ratio_hex"],
            maximum=binder.EXPECTED_SSP4097_CFL["maximum_ratio_hex"],
        ),
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


def _prepared(key: str, disposition: str) -> FakeBridge:
    if disposition == "prepared_overlay":
        plan = plan_imp1_lattice(
            START,
            TARGET,
            float.fromhex(binder.EXPECTED_SSP4097_NEXT_PLAN_CAP_HEX),
        )
    else:
        plan = _frozen_plan(key)
    return FakeBridge(
        disposition,
        evidence={"kind": disposition, "requested_cap_hex": plan.requested_cap.hex()},
        attempted_plan=plan,
        next_plan=plan,
        prepared=object(),
    )


def _world(results_by_key: dict[str, list[object]]):
    captured: dict[str, FakeCaptured] = {}
    views: dict[str, FakeView] = {}
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
        captured[key] = item
        binding = FakeBinding(binder.SOURCE_CONFIGURATION_SHA256_BY_MEMBER[key])
        checkpoint = FakeCheckpoint(results_by_key.get(key, []))
        views[key] = FakeView(key, item, binding, checkpoint)
    origin = FakeOrigin(captured)
    runtime = FakeRuntime(views)
    return origin, runtime, views


def _raw_rk_members() -> list[dict[str, object]]:
    rows = []
    for key in binder.PREDECESSOR_MEMBER_KEYS:
        rows.append(
            {
                "member_key": key,
                "role": "predecessor_evidence",
                "remeasured": False,
                "outcome": "prepared",
                "attempt2_c1r1_receipt_sha256": binder.ATTEMPT2_RK_RECEIPTS[key],
                "attempt2_prepare_disposition": "prepared_fresh",
                "accepted_state_advanced": False,
                "endpoint_adopted": False,
                "source_evaluated_by_rec1": False,
                "source_evaluation_recorded": True,
            }
        )
    return rows


def _raw_ssp_member(record: Mapping) -> dict[str, object]:
    return {
        "member_key": record["member_key"],
        "prepare_disposition": record["reconstructed_prepare_disposition"],
        "prepare_count": record["prepare_count"],
        "overlay_count": record["overlay_count"],
        "c1r1": record["c1r1"],
        "outcome": "prepared",
        "admission_passed": True,
        "accepted_state_advanced": False,
        "endpoint_adopted": False,
    }


def _attempt2_raw_payload() -> dict[str, object]:
    members = []
    for key in binder.PREDECESSOR_MEMBER_KEYS:
        members.append(
            {
                "member_key": key,
                "prepare_disposition": "prepared_fresh",
                "c1r1": {"receipt_sha256": binder.ATTEMPT2_RK_RECEIPTS[key]},
            }
        )
    members.append(
        {
            "member_key": "SSPRK3-4097",
            "prepare_disposition": "cfl_retry_required",
            "c1r1": None,
        }
    )
    return {
        "accepted_state_advanced": False,
        "endpoint_adopted": False,
        "store_published": False,
        "campaign_execution_authorized": False,
        "members": members,
    }


def _source_eval(view: object) -> dict[str, object]:
    del view
    return {
        "diagnostics": {"source_raw_gate_passed": True},
        "source_health": "recorded",
        "raw_health": "recorded",
        "kinetic_health": "recorded",
        "constraint_health": "recorded",
        "mutated_input": False,
        "mutated_monitor": False,
        "mutated_causal": False,
        "mutated_tracer": False,
        "input_state_sha256": "a" * 64,
        "causal_sha256": "b" * 64,
        "monitor_sha256": "c" * 64,
        "tracer_sha256": "d" * 64,
    }


def _established() -> dict[str, object]:
    return {
        "authority_bound": True,
        "freeze_bound": True,
        "attempt2_bound": True,
        "rec1_raw_bound": True,
        "historical_runs_bound": True,
        "planck_bound": True,
        "pro20_absent": True,
        "caps_derived": True,
        "rk_receipts_bound": True,
        "ssp4097_cfl_then_prepared": True,
        "ssp8193_fresh": True,
        "ssp16385_fresh": True,
        "eighteen_channels_bound": True,
        "no_state_advance": True,
        "resources_within_ceilings": True,
        "endpoints_absent": True,
        "runner_label_sufficient": False,
    }


def _load_script():
    spec = importlib.util.spec_from_file_location("pref1_reproducer", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PREF1BinderTests(unittest.TestCase):
    def test_compact_check_is_absent_without_tracked_files(self) -> None:
        with TemporaryDirectory() as tmp:
            with self.assertRaises(binder.HLT17SRCQ1REC1PREF1Error) as raised:
                binder.verify_compact(Path(tmp))
        self.assertIn("absent", str(raised.exception))

    def test_validate_compact_accepts_pinned_certificate(self) -> None:
        config = binder.emit_config_bytes()
        result = binder.canonical_result(binder.expected_compact())
        parsed = binder.validate_compact_result(config, result)
        self.assertEqual(parsed["artifact_payload"]["classification"], binder.CLASSIFICATION)
        with patch.object(binder, "COMPACT_RESULT_SHA256", None):
            with self.assertRaisesRegex(
                binder.HLT17SRCQ1REC1PREF1Error, "not canonical"
            ):
                binder.validate_compact_result(config, b" " + result)

    def test_check_mode_never_accesses_raw_runs_source_or_git(self) -> None:
        config = binder.emit_config_bytes()
        result = binder.canonical_result(binder.expected_compact())
        with TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            (root / "configs/fgc").mkdir(parents=True)
            (root / "results").mkdir()
            (root / binder.CONFIG_PATH).write_bytes(config)
            (root / binder.RESULT_PATH).write_bytes(result)
            with (
                patch.object(binder, "git_read") as git_read,
                patch.object(binder, "inventory_historical_runs") as runs,
                patch.object(binder, "authenticate_rec1_raw") as raw,
                patch.object(binder, "capture_pro20_historical_origin", create=True),
                patch.object(binder, "build_pro20_runtime_origin", create=True),
                patch.object(binder, "independent_source_evaluation") as source,
                patch.object(binder, "reproduce_three_ssp_paths") as reconstruct,
            ):
                binder.verify_compact(root)
            git_read.assert_not_called()
            runs.assert_not_called()
            raw.assert_not_called()
            source.assert_not_called()
            reconstruct.assert_not_called()

    def test_default_script_is_check_blind_and_does_not_reconstruct(self) -> None:
        script = _load_script()
        with (
            patch.object(binder, "bind_live") as live,
            patch.object(binder, "git_read") as git_read,
            patch.object(binder, "inventory_historical_runs") as runs,
        ):
            code = script.main([])
        self.assertEqual(code, 0)
        live.assert_not_called()
        git_read.assert_not_called()
        runs.assert_not_called()

    def test_altered_compact_authority_is_rejected(self) -> None:
        compact = deepcopy(binder.expected_compact())
        compact["artifact_payload"]["authority"]["authority_commit"] = "0" * 40
        with self.assertRaises(binder.HLT17SRCQ1REC1PREF1Error):
            binder.validate_compact_result(
                binder.emit_config_bytes(),
                binder.canonical_result(compact),
            )

    def test_independent_cfl_requires_observed_greater_than_maximum(self) -> None:
        record = binder.independent_cfl_ratio(dict(binder.EXPECTED_SSP4097_CFL))
        self.assertTrue(record["observed_exceeds_maximum"])
        with self.assertRaises(binder.HLT17SRCQ1REC1PREF1Error):
            binder.independent_cfl_ratio({"kind": "cfl_retry"})
        with self.assertRaises(binder.HLT17SRCQ1REC1PREF1Error):
            binder.independent_cfl_ratio(
                {
                    "kind": "cfl_retry",
                    "observed_ratio_hex": "0x1.0000000000000p-4",
                    "maximum_ratio_hex": "0x1.0000000000000p-3",
                }
            )
        with self.assertRaises(binder.HLT17SRCQ1REC1PREF1Error):
            binder.independent_cfl_ratio(
                {
                    "kind": "source_retry",
                    "observed_ratio_hex": "0x1.000000003e80dp-3",
                    "maximum_ratio_hex": "0x1.0000000000000p-3",
                }
            )

    def test_classification_requires_every_independent_predicate(self) -> None:
        facts = _established()
        self.assertEqual(
            binder.classification_if_established(facts), binder.CLASSIFICATION
        )
        facts["ssp4097_cfl_then_prepared"] = False
        with self.assertRaises(binder.HLT17SRCQ1REC1PREF1Error):
            binder.classification_if_established(facts)
        facts = _established()
        facts["runner_label_sufficient"] = True
        with self.assertRaises(binder.HLT17SRCQ1REC1PREF1Error):
            binder.classification_if_established(facts)

    def test_historical_runs_digest_is_sha256sum_lines(self) -> None:
        digest = binder.historical_runs_digest(
            [("runs/a", 1, "a" * 64), ("runs/b", 2, "b" * 64)]
        )
        text = f"{'a' * 64}  runs/a\n{'b' * 64}  runs/b\n"
        self.assertEqual(digest, sha256(text.encode("ascii")).hexdigest())

    def test_ssp4097_cfl_overlay_then_prepared_on_isolated_checkpoint(self) -> None:
        origin, runtime, views = _world(
            {
                "SSPRK3-4097": [
                    _cfl_bridge(),
                    _prepared("SSPRK3-4097", "prepared_overlay"),
                ]
            }
        )
        view = views["SSPRK3-4097"]
        with patch(f"{PREF1}.independent_c1r1_record", return_value=_c1r1("SSPRK3-4097")):
            record = binder.reproduce_ssp_member(
                member_key="SSPRK3-4097",
                view=view,
                captured=origin.members[3],
                initial_cap=float.fromhex(
                    binder.FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER["SSPRK3-4097"]
                ),
                evaluate=_source_eval,
                expected_config=binder.SOURCE_CONFIGURATION_SHA256_BY_MEMBER[
                    "SSPRK3-4097"
                ],
            )
        self.assertEqual(record["reconstructed_prepare_disposition"], "prepared_overlay")
        self.assertEqual(record["overlay_count"], 1)
        self.assertEqual(record["prepare_count"], 2)
        self.assertEqual(
            record["independent_cfl"], dict(binder.EXPECTED_SSP4097_CFL)
        )
        self.assertEqual(
            record["next_plan_cap_hex"], binder.EXPECTED_SSP4097_NEXT_PLAN_CAP_HEX
        )
        self.assertTrue(view.checkpoint.live)
        self.assertEqual(view.checkpoint.absorbed, [])
        isolated = view.checkpoint.isolate()
        del isolated
        self.assertEqual(record["c1r1"]["receipt_sha256"], binder.EXPECTED_C1R1_RECEIPTS["SSPRK3-4097"])

    def test_ssp8193_and_16385_reconstruct_fresh(self) -> None:
        origin, _runtime, views = _world(
            {
                "SSPRK3-8193": [_prepared("SSPRK3-8193", "prepared_fresh")],
                "SSPRK3-16385": [_prepared("SSPRK3-16385", "prepared_fresh")],
            }
        )
        for key, index in (("SSPRK3-8193", 4), ("SSPRK3-16385", 5)):
            with patch(f"{PREF1}.independent_c1r1_record", return_value=_c1r1(key)):
                record = binder.reproduce_ssp_member(
                    member_key=key,
                    view=views[key],
                    captured=origin.members[index],
                    initial_cap=float.fromhex(
                        binder.FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER[key]
                    ),
                    evaluate=_source_eval,
                    expected_config=binder.SOURCE_CONFIGURATION_SHA256_BY_MEMBER[key],
                )
            self.assertEqual(record["reconstructed_prepare_disposition"], "prepared_fresh")
            self.assertEqual(record["overlay_count"], 0)
            self.assertEqual(record["prepare_count"], 1)
            self.assertIsNone(record["independent_cfl"])

    def test_swapped_cap_and_plan_are_rejected(self) -> None:
        origin, _runtime, views = _world(
            {"SSPRK3-4097": [_prepared("SSPRK3-8193", "prepared_fresh")]}
        )
        with patch(f"{PREF1}.independent_c1r1_record", return_value=_c1r1("SSPRK3-4097")):
            with self.assertRaises(binder.HLT17SRCQ1REC1PREF1Error):
                binder.reproduce_ssp_member(
                    member_key="SSPRK3-4097",
                    view=views["SSPRK3-4097"],
                    captured=origin.members[3],
                    initial_cap=float.fromhex(
                        binder.FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER["SSPRK3-8193"]
                    ),
                    evaluate=_source_eval,
                    expected_config=binder.SOURCE_CONFIGURATION_SHA256_BY_MEMBER[
                        "SSPRK3-4097"
                    ],
                )

    def test_false_or_absent_cfl_evidence_is_rejected(self) -> None:
        origin, _runtime, views = _world({})
        bad = _cfl_bridge()
        bad.evidence = {"kind": "cfl_retry"}
        views["SSPRK3-4097"].checkpoint = FakeCheckpoint(
            [bad, _prepared("SSPRK3-4097", "prepared_overlay")]
        )
        with patch(f"{PREF1}.independent_c1r1_record", return_value=_c1r1("SSPRK3-4097")):
            with self.assertRaises(binder.HLT17SRCQ1REC1PREF1Error):
                binder.reproduce_ssp_member(
                    member_key="SSPRK3-4097",
                    view=views["SSPRK3-4097"],
                    captured=origin.members[3],
                    initial_cap=float.fromhex(
                        binder.FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER["SSPRK3-4097"]
                    ),
                    evaluate=_source_eval,
                    expected_config=binder.SOURCE_CONFIGURATION_SHA256_BY_MEMBER[
                        "SSPRK3-4097"
                    ],
                )

    def test_changed_channel_debit_or_receipt_is_rejected(self) -> None:
        origin, runtime, views = _world(
            {
                "SSPRK3-4097": [
                    _cfl_bridge(),
                    _prepared("SSPRK3-4097", "prepared_overlay"),
                ],
                "SSPRK3-8193": [_prepared("SSPRK3-8193", "prepared_fresh")],
                "SSPRK3-16385": [_prepared("SSPRK3-16385", "prepared_fresh")],
            }
        )
        with patch(
            f"{PREF1}.independent_c1r1_record",
            side_effect=lambda _prepared: _c1r1("SSPRK3-4097"),
        ):
            reconstructed = binder.reproduce_ssp_member(
                member_key="SSPRK3-4097",
                view=views["SSPRK3-4097"],
                captured=origin.members[3],
                initial_cap=float.fromhex(
                    binder.FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER["SSPRK3-4097"]
                ),
                evaluate=_source_eval,
                expected_config=binder.SOURCE_CONFIGURATION_SHA256_BY_MEMBER[
                    "SSPRK3-4097"
                ],
            )
        raw = _raw_ssp_member(reconstructed)
        raw["c1r1"] = dict(raw["c1r1"])
        raw["c1r1"]["receipt_sha256"] = "0" * 64
        with self.assertRaises(binder.HLT17SRCQ1REC1PREF1Error):
            binder.compare_ssp_to_raw(reconstructed, raw)
        raw = _raw_ssp_member(reconstructed)
        raw["c1r1"] = dict(raw["c1r1"])
        decisions = list(raw["c1r1"]["eighteen_channel_decisions"])
        first = dict(decisions[0])
        first["public_fine_debit"] = "1/2"
        decisions[0] = first
        raw["c1r1"]["eighteen_channel_decisions"] = decisions
        with self.assertRaises(binder.HLT17SRCQ1REC1PREF1Error):
            binder.compare_ssp_to_raw(reconstructed, raw)

    def test_partial_terminal_and_endpoint_flags_are_rejected(self) -> None:
        attempt2 = _attempt2_raw_payload()
        rec1 = {
            "members": _raw_rk_members()[:2],
            "accepted_state_advanced": False,
            "endpoint_adopted": False,
            "store_published": False,
            "campaign_execution_authorized": False,
        }
        with self.assertRaises(binder.HLT17SRCQ1REC1PREF1Error):
            binder.bind_attempt2_rk_receipts(attempt2, rec1)
        with self.assertRaises(Exception):
            from recursive_horizons.fgc.evolution.hlt17_srcq1 import (
                refuse_endpoint_payload,
            )

            refuse_endpoint_payload({"corrected_endpoint": [1.0]})

    def test_resource_excess_is_rejected(self) -> None:
        ceilings = ResourceCeilings(
            max_rss_bytes=1024,
            max_wall_seconds_per_member={key: 900.0 for key in MEMBER_KEYS},
            max_total_wall_seconds=1.0,
        )
        with self.assertRaises(binder.HLT17SRCQ1REC1PREF1ResourceError):
            binder._check_resources(ceilings, 2.0, 2048)

    def test_live_synthetic_path_prints_without_write_or_publication(self) -> None:
        origin, runtime, _views = _world(
            {
                "SSPRK3-4097": [
                    _cfl_bridge(),
                    _prepared("SSPRK3-4097", "prepared_overlay"),
                ],
                "SSPRK3-8193": [_prepared("SSPRK3-8193", "prepared_fresh")],
                "SSPRK3-16385": [_prepared("SSPRK3-16385", "prepared_fresh")],
            }
        )
        rec1_members = _raw_rk_members()
        for key, disposition, overlay, prepares in (
            ("SSPRK3-4097", "prepared_overlay", 1, 2),
            ("SSPRK3-8193", "prepared_fresh", 0, 1),
            ("SSPRK3-16385", "prepared_fresh", 0, 1),
        ):
            rec1_members.append(
                _raw_ssp_member(
                    {
                        "member_key": key,
                        "reconstructed_prepare_disposition": disposition,
                        "prepare_count": prepares,
                        "overlay_count": overlay,
                        "c1r1": _c1r1(key),
                    }
                )
            )
        rec1_raw = {
            "accepted_state_advanced": False,
            "endpoint_adopted": False,
            "store_published": False,
            "campaign_execution_authorized": False,
            "halted_outcome": None,
            "rk_members_remeasured": False,
            "members": rec1_members,
            "resource_facts": {
                "total_wall_seconds": binder.RECORDED_RAW_WALL_SECONDS,
                "rss_bytes": binder.RECORDED_RAW_RSS_BYTES,
            },
        }
        calls: list[str] = []
        tracked_before = {
            binder.CONFIG_PATH: (ROOT / binder.CONFIG_PATH).read_bytes(),
            binder.RESULT_PATH: (ROOT / binder.RESULT_PATH).read_bytes(),
        }

        def _shadow(_prepared: object) -> dict[str, object]:
            key = ("SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385")[len(calls)]
            calls.append(key)
            return _c1r1(key)

        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            before = list(root.rglob("*"))
            with (
                patch(f"{PREF1}.independent_c1r1_record", side_effect=_shadow),
                patch(f"{PREF1}.static_input_identities", return_value={}),
            ):
                result = binder.bind_live(
                    ROOT,
                    origin=origin,
                    runtime=runtime,
                    source_evaluator=_source_eval,
                    static_input_bytes={"x": b"x"},
                    environment={key: "x" for key in binder.ENVIRONMENT_KEYS},
                    skip_git=True,
                    skip_raw=True,
                    rec1_raw_payload=rec1_raw,
                    attempt2_raw_payload=_attempt2_raw_payload(),
                    historical_runs={
                        "file_count": binder.HISTORICAL_RUNS_FILE_COUNT,
                        "byte_count": binder.HISTORICAL_RUNS_BYTE_COUNT,
                        "digest": binder.HISTORICAL_RUNS_DIGEST,
                        "excluded_namespace": binder.REC1_OUTPUT_NAMESPACE,
                    },
                    planck={
                        path: {
                            "present": False,
                            "byte_count": size,
                            "sha256": digest,
                        }
                        for path, (size, digest) in binder.PLANCK_FILES.items()
                    },
                    clock=lambda: 0.0,
                    rss_bytes=lambda: 1,
                )
            after = list(root.rglob("*"))
        self.assertEqual(before, after)
        self.assertFalse(result["written"])
        self.assertEqual(result["classification"], binder.CLASSIFICATION)
        self.assertEqual(result["result"], binder.expected_compact())
        for relative, payload in tracked_before.items():
            self.assertEqual((ROOT / relative).read_bytes(), payload)

    def test_live_refuses_endpoint_payloads(self) -> None:
        compact = deepcopy(binder.expected_compact())
        compact["artifact_payload"]["corrected_endpoint"] = [0.0]
        with self.assertRaises(binder.HLT17SRCQ1REC1PREF1Error):
            binder.canonical_result(compact)

    def test_altered_source_configuration_is_rejected(self) -> None:
        origin, _runtime, views = _world(
            {"SSPRK3-8193": [_prepared("SSPRK3-8193", "prepared_fresh")]}
        )
        with patch(f"{PREF1}.independent_c1r1_record", return_value=_c1r1("SSPRK3-8193")):
            with self.assertRaises(binder.HLT17SRCQ1REC1PREF1Error):
                binder.reproduce_ssp_member(
                    member_key="SSPRK3-8193",
                    view=views["SSPRK3-8193"],
                    captured=origin.members[4],
                    initial_cap=float.fromhex(
                        binder.FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER["SSPRK3-8193"]
                    ),
                    evaluate=_source_eval,
                    expected_config="0" * 64,
                )

    def test_store_or_state_flags_on_compact_are_rejected(self) -> None:
        compact = deepcopy(binder.expected_compact())
        compact["artifact_payload"]["claims"]["accepted_state_advanced"] = True
        with self.assertRaises(binder.HLT17SRCQ1REC1PREF1Error):
            binder.validate_compact_result(
                binder.emit_config_bytes(),
                binder.canonical_result(compact),
            )

    def test_binder_and_script_do_not_import_rec1_runner_or_authority(self) -> None:
        for path in (CORE, SCRIPT):
            source = path.read_text(encoding="utf-8")
            tree = ast.parse(source)
            imported: set[str] = set()
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    imported.update(alias.name for alias in node.names)
                elif isinstance(node, ast.ImportFrom) and node.module:
                    imported.add(node.module)
                    imported.update(alias.name for alias in node.names)
            self.assertFalse(
                {
                    name
                    for name in imported
                    if any(
                        needle in name
                        for needle in (
                            "hlt17_srcq1_rec1_auth1",
                            "qualify_fgc_hlt17_srcq1_rec1",
                        )
                    )
                    or name.endswith(".hlt17_srcq1_rec1")
                    or name == "hlt17_srcq1_rec1"
                }
            )
            self.assertTrue(
                {
                    "c1r1_shadow_record",
                    "derive_prospective_requested_caps",
                    "evaluate_bound_source_once",
                }.isdisjoint(imported)
            )
        probe = subprocess.run(
            [
                sys.executable,
                "-I",
                "-B",
                "-c",
                (
                    "import sys; sys.path.insert(0, 'src'); "
                    "import recursive_horizons.fgc.evolution.hlt17_srcq1_rec1_pref1_binder; "
                    "needles=("
                    "'hlt17_srcq1_rec1_auth1',"
                    "'qualify_fgc_hlt17_srcq1_rec1',"
                    "); "
                    "bad=[name for name in sys.modules if any("
                    "needle in name for needle in needles) or "
                    "name.endswith('.hlt17_srcq1_rec1')]; "
                    "assert not bad, bad"
                ),
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        self.assertEqual(probe.returncode, 0, probe.stdout + probe.stderr)

    def test_script_has_no_publication_or_write_paths(self) -> None:
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertNotIn("publish_exclusive_file", source)
        self.assertNotIn("publish_exclusive_directory", source)
        self.assertNotIn("write_result", source)
        self.assertNotIn("--write", source)
        self.assertIn("never publishes", source)

    def test_swapped_member_compare_is_rejected(self) -> None:
        record = {
            "member_key": "SSPRK3-8193",
            "reconstructed_prepare_disposition": "prepared_fresh",
            "prepare_count": 1,
            "overlay_count": 0,
            "independent_cfl": None,
            "c1r1": _c1r1("SSPRK3-8193"),
        }
        raw = _raw_ssp_member(record)
        raw["member_key"] = "SSPRK3-16385"
        with self.assertRaises(binder.HLT17SRCQ1REC1PREF1Error):
            binder.compare_ssp_to_raw(record, raw)

    def test_altered_environment_pin_in_compact_is_rejected(self) -> None:
        compact = deepcopy(binder.expected_compact())
        compact["artifact_payload"]["authority"]["source_closure_sha256"] = "1" * 64
        with self.assertRaises(binder.HLT17SRCQ1REC1PREF1Error):
            binder.validate_compact_result(
                binder.emit_config_bytes(),
                binder.canonical_result(compact),
            )

    def test_directory_and_history_pins_are_exact(self) -> None:
        compact = deepcopy(binder.expected_compact())
        compact["artifact_payload"]["raw_terminal"]["directory_payload_sha256"] = "2" * 64
        with self.assertRaises(binder.HLT17SRCQ1REC1PREF1Error):
            binder.validate_compact_result(
                binder.emit_config_bytes(),
                binder.canonical_result(compact),
            )
        compact = deepcopy(binder.expected_compact())
        compact["artifact_payload"]["historical_runs_baseline"]["digest"] = "3" * 64
        with self.assertRaises(binder.HLT17SRCQ1REC1PREF1Error):
            binder.validate_compact_result(
                binder.emit_config_bytes(),
                binder.canonical_result(compact),
            )


if __name__ == "__main__":
    unittest.main()
