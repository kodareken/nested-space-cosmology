"""Synthetic/mocked SRCQ1 tests. No physical source, live store, or production grid."""

from __future__ import annotations

import ast
from dataclasses import dataclass
from hashlib import sha256
import importlib.util
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
    encode_tdg7_plan,
    plan_imp1_lattice,
)
from recursive_horizons.fgc.evolution.hlt17_srcq1 import (  # noqa: E402
    ARTIFACT_ID,
    AUTHORITY_SEAM,
    HLT17SRCQ1AuthoritySeamError,
    HLT17SRCQ1CapError,
    HLT17SRCQ1IdentityError,
    HLT17SRCQ1PublicationError,
    HLT17SRCQ1SourceError,
    ResourceCeilings,
    SCHEMA,
    STATIC_INPUT_PATHS,
    c1r1_shadow_record,
    derive_prospective_requested_caps,
    evaluate_bound_source_once,
    implementation_identity,
    load_coordinator_authority_delta_validator,
    plan_status,
    publish_qualification_result,
    qualify_and_publish,
    qualify_physical_source_origin,
    recommended_resource_ceilings,
    refuse_endpoint_payload,
    require_authority_delta,
    static_input_identities,
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


CORE = ROOT / "src/recursive_horizons/fgc/evolution/hlt17_srcq1.py"
CLI = ROOT / "scripts/qualify_fgc_hlt17_srcq1.py"
START = 23.0 / 16.0
TARGET = 3.0 / 2.0
CAP = 1.0 / 16.0


def _digest(label: str) -> str:
    return sha256(label.encode("ascii")).hexdigest()


def _plan():
    return plan_imp1_lattice(START, TARGET, CAP)


def _plan_mapping() -> dict[str, object]:
    return encode_tdg7_plan(_plan())


def _cursor_bytes() -> bytes:
    return canonical_json_bytes({"executed_plan": _plan_mapping()})


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
    def __init__(self, key: str, *, step: int, serial: int, source: int, cfl: int) -> None:
        self.member_key = key
        self.step_index = step
        self.transaction_serial = serial
        self.source_retry_count = source
        self.cfl_retry_count = cfl
        self.accepted_time = {"rational": "23/16", "binary64_hex": ORIGIN_TIME_HEX}
        self.original_cursor_bytes = canonical_json_bytes({
            "mode": "FRESH_READY", "retry_successor_payload_or_none": None
        })
        self.causal_bytes = canonical_json_bytes({
            "accepted_time": START,
            "accumulated_characteristic_distance": 0.0,
            "previous_speed_upper": 1.2,
        })


class FakeCheckpoint:
    def __init__(self, disposition: str = "prepared_fresh", prepared: object | None = None) -> None:
        self.disposition = disposition
        self.prepared = object() if prepared is None and disposition == "prepared_fresh" else prepared
        self.accepted_state_advanced = False

    def agree(self) -> None:
        return None

    def prepare(self, *, requested_cap: object) -> SimpleNamespace:
        del requested_cap
        return SimpleNamespace(
            disposition=self.disposition,
            prepared=self.prepared,
            accepted_state_advanced=self.accepted_state_advanced,
        )


class FakeMember:
    def __init__(self, key: str, captured: FakeCaptured, binding: FakeBinding) -> None:
        label, points = key.split("-", 1)
        integrator = PRIMARY_METHOD if label == "RK4" else COMPARATOR_METHOD
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
            FakeCausal(),
            grid_spacing=128.0 / (self.identity.point_count - 1),
            cfl_maximum=1.0 / 8.0,
        )
        self.tracers = FakeTracers(
            np.array([0.25, 0.75], dtype=np.float64),
            np.array([0.0, 0.0], dtype=np.float64),
        )
        self.temporal_ledger = seed_imp1_ledger(
            TDG6TemporalLedger.zero(initial_time=START),
            method=integrator,
            state_sha256=array_content_sha256(
                self.state.u, self.state.p, self.state.q
            ),
            step_index=self.step_index,
            transaction_serial=self.transaction_serial,
            origin_receipt_sha256=_digest("origin-receipt"),
        )
        self.source_binding = binding


class FakeView:
    def __init__(self, key: str, captured: FakeCaptured, binding: FakeBinding) -> None:
        self.member_key = key
        self.source_binding = binding
        self.source_configuration_sha256 = binding.configuration_sha256
        self.member = FakeMember(key, captured, binding)
        self.checkpoint = FakeCheckpoint()


class FakeOrigin:
    def __init__(self, captured: dict[str, FakeCaptured]) -> None:
        self.capture_sha256 = _digest("origin")
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


COUNTERS = {
    "RK4-2049": (342, 1710, 0, 217),
    "RK4-4097": (673, 3365, 0, 443),
    "RK4-8193": (976, 4880, 0, 161),
    "SSPRK3-4097": (631, 2524, 0, 358),
    "SSPRK3-8193": (987, 3948, 0, 187),
    "SSPRK3-16385": (1771, 7084, 0, 0),
}


def _static_pin():
    pins = {path: sha256(b"x").hexdigest() for path in STATIC_INPUT_PATHS}
    return patch(
        "recursive_horizons.fgc.evolution.hlt17_srcq1.STATIC_INPUT_SHA256", pins
    )


def _world(*, mutate: str | None = None, raw_fail: str | None = None):
    captured = {}
    views = {}
    configs = {}
    for key in MEMBER_KEYS:
        step, serial, source, cfl = COUNTERS[key]
        item = FakeCaptured(key, step=step, serial=serial, source=source, cfl=cfl)
        digest = _digest(f"config-{key}")
        binding = FakeBinding(
            digest,
            mutate=mutate == key,
            raw_gate=raw_fail != key,
        )
        captured[key] = item
        views[key] = FakeView(key, item, binding)
        configs[key] = digest
    origin = FakeOrigin(captured)
    runtime = SimpleNamespace(
        members=views,
        captured_origin_sha256=origin.capture_sha256,
        source_closure_sha256=_digest("closure"),
        construction_sha256=_digest("construction"),
    )
    return origin, runtime, configs


def _cap_hex(origin, runtime) -> dict[str, str]:
    caps, _facts = derive_prospective_requested_caps(origin, runtime)
    return {key: caps[key].hex() for key in MEMBER_KEYS}


def _load_cli():
    spec = importlib.util.spec_from_file_location("qualify_fgc_hlt17_srcq1_cli", CLI)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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
        "recursive_horizons.fgc.evolution.hlt17_srcq1.c1r1_shadow_record",
        return_value={
            "eighteen_channel_decisions": _eighteen(),
            "admission_passed": True,
            "public_fine_debits": ["0"] * 18,
        },
    )


class SRCQ1ContractTests(unittest.TestCase):
    def test_plan_status_does_not_execute(self) -> None:
        status = plan_status()
        self.assertEqual(status["artifact_id"], ARTIFACT_ID)
        self.assertEqual(status["schema"], SCHEMA)
        self.assertIs(status["physical_six_member_qualification_executed"], False)
        self.assertIs(status["campaign_execution_authorized"], False)
        self.assertEqual(status["member_keys"], list(MEMBER_KEYS))
        self.assertIn("exact --authority-commit", status["one_real_run"]["requires"])
        self.assertEqual(status["authority_delta_validation"], AUTHORITY_SEAM)
        self.assertTrue(status["resource_ceilings"]["coordinator_freeze_required"])
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
        self.assertNotIn("run_fgc_gr0_calibration_v6", joined)
        self.assertNotIn("run_fgc_gr0_calibration_v14", joined)
        cli_tree = ast.parse(CLI.read_text("utf-8"))
        cli_imported: list[str] = []
        for node in ast.walk(cli_tree):
            if isinstance(node, ast.Import):
                cli_imported.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                cli_imported.append(node.module or "")
                cli_imported.extend(alias.name for alias in node.names)
        self.assertNotIn("run_fgc_gr0_calibration", " ".join(cli_imported))

    def test_default_cli_is_status_only(self) -> None:
        cli = _load_cli()
        with patch.object(
            cli, "capture_pro20_historical_origin", side_effect=AssertionError("store")
        ), patch.object(
            cli, "observe_environment", side_effect=AssertionError("env")
        ), patch.object(
            cli, "qualify_and_publish", side_effect=AssertionError("qualify")
        ):
            code = cli.main([])
        self.assertEqual(code, 0)

    def test_cli_refuses_partial_real_run_args(self) -> None:
        cli = _load_cli()
        with patch.object(
            cli, "capture_pro20_historical_origin", side_effect=AssertionError("store")
        ):
            self.assertEqual(cli.main(["--authority-commit", "a" * 40]), 2)
            self.assertEqual(cli.main(["--output-directory", "/tmp/srcq1-out"]), 2)

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
            ):
                code = cli.main(
                    ["--authority-commit", "a" * 40, "--output-directory", str(target)]
                )
        self.assertEqual(code, 1)

    def test_authority_seam_refuses_boolean_and_missing_module(self) -> None:
        with patch(
            "recursive_horizons.fgc.evolution.hlt17_srcq1.importlib.import_module",
            side_effect=ImportError("missing freeze module"),
        ), self.assertRaises(HLT17SRCQ1AuthoritySeamError):
            load_coordinator_authority_delta_validator()
        with self.assertRaises(HLT17SRCQ1AuthoritySeamError):
            require_authority_delta(
                repository_root=ROOT,
                authority_commit="a" * 40,
                implementation_sha256="b" * 64,
                validator=lambda **kwargs: True,
            )
        with self.assertRaises(HLT17SRCQ1AuthoritySeamError):
            require_authority_delta(
                repository_root=ROOT,
                authority_commit="HEAD",
                implementation_sha256="b" * 64,
                validator=lambda **kwargs: {"ok": True},
            )
        receipt = require_authority_delta(
            repository_root=ROOT,
            authority_commit="a" * 40,
            implementation_sha256="b" * 64,
            validator=lambda **kwargs: {"delta": "exact", **kwargs},
        )
        self.assertIs(receipt["boolean_flag_accepted"], False)
        self.assertEqual(receipt["authority_commit"], "a" * 40)

    def test_prospective_cap_formula_and_authority_refusal(self) -> None:
        origin, runtime, _configs = _world()
        caps, facts = derive_prospective_requested_caps(origin, runtime)
        expected = _cap_hex(origin, runtime)
        for key in MEMBER_KEYS:
            self.assertEqual(caps[key].hex(), expected[key])
            self.assertEqual(
                facts[key]["formula"],
                "cfl_maximum*grid_spacing/inherited_previous_speed_upper",
            )
            self.assertNotIn("executed_plan", facts[key])
        checked, _facts = derive_prospective_requested_caps(
            origin, runtime, expected_cap_hex_by_member=expected
        )
        self.assertEqual({key: checked[key].hex() for key in MEMBER_KEYS}, expected)
        malformed = dict(expected)
        malformed["RK4-2049"] = "not-a-hex"
        with self.assertRaises(HLT17SRCQ1CapError):
            derive_prospective_requested_caps(
                origin, runtime, expected_cap_hex_by_member=malformed
            )
        swapped = dict(expected)
        swapped["RK4-2049"], swapped["RK4-8193"] = (
            swapped["RK4-8193"], swapped["RK4-2049"]
        )
        with self.assertRaisesRegex(HLT17SRCQ1CapError, "authority requested cap"):
            derive_prospective_requested_caps(
                origin, runtime, expected_cap_hex_by_member=swapped
            )
        runtime.members["RK4-2049"].member.transaction.cfl_maximum = 0.25
        with self.assertRaisesRegex(HLT17SRCQ1CapError, "authority requested cap"):
            derive_prospective_requested_caps(
                origin, runtime, expected_cap_hex_by_member=expected
            )

    def test_six_member_order_and_counter_preservation(self) -> None:
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

        def builder(*args, **kwargs):
            raise AssertionError("physical factory")

        with _static_pin(), _passing_c1r1():
            result = qualify_physical_source_origin(
                origin,
                static_input_bytes={path: b"x" for path in STATIC_INPUT_PATHS},
                source_closure_sha256=runtime.source_closure_sha256,
                environment={"python_implementation": "CPython"},
                runtime_origin=runtime,
                runtime_origin_builder=builder,
                source_evaluator=evaluate,
                expected_origin_sha256=origin.capture_sha256,
                expected_implementation=hlt17_implementation_identity(),
                expected_source_configurations=configs,
                expected_environment={"python_implementation": "CPython"},
                rss_bytes=lambda: 1,
                clock=Clock(),
            )
        self.assertEqual(called, list(MEMBER_KEYS))
        self.assertEqual(result["member_order"], list(MEMBER_KEYS))
        self.assertEqual(result["completed_member_count"], 6)
        self.assertIsNone(result["halted_outcome"])
        self.assertIs(result["accepted_state_advanced"], False)
        for record, key in zip(result["members"], MEMBER_KEYS, strict=True):
            self.assertEqual(record["member_key"], key)
            self.assertEqual(record["fingerprints_before"], record["fingerprints_after"])
            self.assertEqual(record["owned"]["member_key"], key)

        reversed_members = {
            key: runtime.members[key] for key in reversed(MEMBER_KEYS)
        }
        runtime.members = reversed_members
        with _static_pin(), self.assertRaises(HLT17SRCQ1IdentityError):
            qualify_physical_source_origin(
                origin,
                static_input_bytes={path: b"x" for path in STATIC_INPUT_PATHS},
                source_closure_sha256=runtime.source_closure_sha256,
                environment={"python_implementation": "CPython"},
                runtime_origin=runtime,
                runtime_origin_builder=builder,
                rss_bytes=lambda: 1,
            )

        origin, runtime, configs = _world()
        runtime.members["RK4-2049"].member.step_index += 1
        with _static_pin(), self.assertRaises(HLT17SRCQ1IdentityError):
            qualify_physical_source_origin(
                origin,
                static_input_bytes={path: b"x" for path in STATIC_INPUT_PATHS},
                source_closure_sha256=runtime.source_closure_sha256,
                environment={"python_implementation": "CPython"},
                runtime_origin=runtime,
                runtime_origin_builder=builder,
                source_evaluator=lambda view: (_ for _ in ()).throw(
                    AssertionError("source")
                ),
                rss_bytes=lambda: 1,
            )

    def test_source_mutation_and_diagnostic_failure(self) -> None:
        origin, runtime, configs = _world()
        member = runtime.members["RK4-2049"].member

        class MutatingTracer(FakeBinding):
            def rhs(self, time: float, state: EvolutionState) -> EvolutionRHS:
                member.tracers.positions[0] += 1.0
                return super().rhs(time, state)

        with self.assertRaises(HLT17SRCQ1SourceError):
            evaluate_bound_source_once(
                binding=MutatingTracer(_digest("tracer")),
                time=START,
                state=member.state,
                transaction=member.transaction,
                tracers=member.tracers,
            )
        origin, runtime, configs = _world(raw_fail="RK4-2049")
        with _static_pin():
            result = qualify_physical_source_origin(
                origin,
                static_input_bytes={path: b"x" for path in STATIC_INPUT_PATHS},
                source_closure_sha256=runtime.source_closure_sha256,
                environment={"python_implementation": "CPython"},
                runtime_origin=runtime,
                runtime_origin_builder=lambda *a, **k: (_ for _ in ()).throw(
                    AssertionError("factory")
                ),
                rss_bytes=lambda: 1,
                clock=Clock(),
            )
        self.assertEqual(result["halted_outcome"], "source")
        self.assertEqual(result["completed_member_count"], 1)
        self.assertEqual(result["members"][0]["outcome"], "source")

        monitor = FakeMonitor()
        causal = FakeCausal()
        transaction = FakeTransaction(monitor, causal)
        tracers = FakeTracers(
            np.array([0.25], dtype=np.float64),
            np.array([0.0], dtype=np.float64),
        )

        class MutatingMonitor(FakeBinding):
            def rhs(self, time: float, state: EvolutionState) -> EvolutionRHS:
                transaction.state.accepted_stage_count += 1
                return super().rhs(time, state)

        with self.assertRaises(HLT17SRCQ1SourceError):
            evaluate_bound_source_once(
                binding=MutatingMonitor(_digest("mon")),
                time=START,
                state=_state(),
                transaction=transaction,
                tracers=tracers,
            )

    def test_foreign_identities_are_refused(self) -> None:
        origin, runtime, configs = _world()

        def builder(*args, **kwargs):
            raise AssertionError("factory")

        with _static_pin(), self.assertRaises(HLT17SRCQ1IdentityError):
            qualify_physical_source_origin(
                origin,
                static_input_bytes={path: b"x" for path in STATIC_INPUT_PATHS},
                source_closure_sha256=runtime.source_closure_sha256,
                environment={"python_implementation": "CPython"},
                runtime_origin=runtime,
                runtime_origin_builder=builder,
                expected_origin_sha256="a" * 64,
                rss_bytes=lambda: 1,
            )
        foreign = dict(hlt17_implementation_identity())
        foreign["implementation_id"] = "foreign_c1r1"
        with _static_pin(), self.assertRaises((HLT17SRCQ1IdentityError, ValueError)):
            qualify_physical_source_origin(
                origin,
                static_input_bytes={path: b"x" for path in STATIC_INPUT_PATHS},
                source_closure_sha256=runtime.source_closure_sha256,
                environment={"python_implementation": "CPython"},
                runtime_origin=runtime,
                runtime_origin_builder=builder,
                expected_implementation=foreign,
                rss_bytes=lambda: 1,
            )
        with _static_pin(), self.assertRaises(HLT17SRCQ1IdentityError):
            qualify_physical_source_origin(
                origin,
                static_input_bytes={path: b"x" for path in STATIC_INPUT_PATHS},
                source_closure_sha256=runtime.source_closure_sha256,
                environment={"python_implementation": "CPython"},
                runtime_origin=runtime,
                runtime_origin_builder=builder,
                expected_environment={"python_implementation": "PyPy"},
                rss_bytes=lambda: 1,
            )
        wrong_configs = dict(configs)
        wrong_configs["RK4-2049"] = "a" * 64
        with _static_pin(), self.assertRaises(HLT17SRCQ1IdentityError):
            qualify_physical_source_origin(
                origin,
                static_input_bytes={path: b"x" for path in STATIC_INPUT_PATHS},
                source_closure_sha256=runtime.source_closure_sha256,
                environment={"python_implementation": "CPython"},
                runtime_origin=runtime,
                runtime_origin_builder=builder,
                expected_source_configurations=wrong_configs,
                rss_bytes=lambda: 1,
            )
        runtime.source_closure_sha256 = "a" * 64
        with _static_pin(), self.assertRaises(HLT17SRCQ1IdentityError):
            qualify_physical_source_origin(
                origin,
                static_input_bytes={path: b"x" for path in STATIC_INPUT_PATHS},
                source_closure_sha256=_digest("closure"),
                environment={"python_implementation": "CPython"},
                runtime_origin=runtime,
                runtime_origin_builder=builder,
                rss_bytes=lambda: 1,
            )

    def test_endpoint_serialization_is_refused(self) -> None:
        refuse_endpoint_payload({"corrected_endpoint_committable": False})
        with self.assertRaises(Exception):
            refuse_endpoint_payload({"corrected_endpoint": {"u": [1.0]}})
        with self.assertRaises(Exception):
            refuse_endpoint_payload({"shadow_endpoint": True})
        with self.assertRaises(Exception):
            refuse_endpoint_payload({"fine_endpoint_state": {"q": 1}})
        with self.assertRaises(Exception):
            refuse_endpoint_payload(np.zeros(3))
        prepared = SimpleNamespace(
            implementation_id="foreign",
            reference_wire_evaluator_id="x",
            plan=_plan(),
            receipt_sha256="a" * 64,
            assessment=SimpleNamespace(
                as_mapping=lambda: {
                    "channels": _eighteen(),
                    "public_fine_debits": ["0"] * 18,
                    "admission_passed": True,
                    "family_sha256": "b" * 64,
                    "assessment_sha256": "c" * 64,
                },
                admission_passed=True,
                assessment_sha256="c" * 64,
                family_sha256="b" * 64,
            ),
        )
        with patch(
            "recursive_horizons.fgc.evolution.hlt17_srcq1.require_hlt17_prepared",
            return_value=prepared,
        ):
            with self.assertRaises(HLT17SRCQ1IdentityError):
                c1r1_shadow_record(prepared)

    def test_resource_and_typed_stops(self) -> None:
        origin, runtime, configs = _world()
        ceilings = ResourceCeilings(
            max_rss_bytes=1,
            max_wall_seconds_per_member={
                key: 900.0 if key != "RK4-2049" else 0.0000001 for key in MEMBER_KEYS
            },
            max_total_wall_seconds=3600.0,
            coordinator_freeze_required=True,
        )
        with _static_pin(), _passing_c1r1():
            result = qualify_physical_source_origin(
                origin,
                static_input_bytes={path: b"x" for path in STATIC_INPUT_PATHS},
                source_closure_sha256=runtime.source_closure_sha256,
                environment={"python_implementation": "CPython"},
                resource_ceilings=ceilings,
                runtime_origin=runtime,
                runtime_origin_builder=lambda *a, **k: (_ for _ in ()).throw(
                    AssertionError("factory")
                ),
                rss_bytes=lambda: 1,
                clock=Clock(jump_after=3, jump=10.0),
            )
        self.assertEqual(result["halted_outcome"], "resource")
        self.assertEqual(result["members"][0]["outcome"], "resource")

        origin, runtime, configs = _world()
        runtime.members["RK4-2049"].checkpoint.disposition = "cfl_retry_required"
        runtime.members["RK4-2049"].checkpoint.prepared = None
        with _static_pin(), _passing_c1r1():
            result = qualify_physical_source_origin(
                origin,
                static_input_bytes={path: b"x" for path in STATIC_INPUT_PATHS},
                source_closure_sha256=runtime.source_closure_sha256,
                environment={"python_implementation": "CPython"},
                runtime_origin=runtime,
                runtime_origin_builder=lambda *a, **k: (_ for _ in ()).throw(
                    AssertionError("factory")
                ),
                rss_bytes=lambda: 1,
                clock=Clock(),
            )
        self.assertEqual(result["halted_outcome"], "cfl")
        runtime.members["RK4-2049"].checkpoint.disposition = "temporal_retry_required"
        runtime.members["RK4-2049"].checkpoint.prepared = None
        with _static_pin(), _passing_c1r1():
            result = qualify_physical_source_origin(
                origin,
                static_input_bytes={path: b"x" for path in STATIC_INPUT_PATHS},
                source_closure_sha256=runtime.source_closure_sha256,
                environment={"python_implementation": "CPython"},
                runtime_origin=runtime,
                runtime_origin_builder=lambda *a, **k: (_ for _ in ()).throw(
                    AssertionError("factory")
                ),
                rss_bytes=lambda: 1,
                clock=Clock(),
            )
        self.assertEqual(result["halted_outcome"], "temporal")

        ceilings = recommended_resource_ceilings()
        self.assertTrue(ceilings.coordinator_freeze_required)
        self.assertEqual(ceilings.max_rss_bytes, 4 * 1024 * 1024 * 1024)
        self.assertEqual(ceilings.max_total_wall_seconds, 3600.0)
        self.assertTrue(all(value == 900.0 for value in ceilings.max_wall_seconds_per_member.values()))

    def test_static_input_hashes_and_inventory(self) -> None:
        payloads = {path: f"payload-{path}".encode() for path in STATIC_INPUT_PATHS}
        with self.assertRaises(HLT17SRCQ1IdentityError):
            static_input_identities(payloads)
        pins = {path: sha256(payloads[path]).hexdigest() for path in STATIC_INPUT_PATHS}
        with patch(
            "recursive_horizons.fgc.evolution.hlt17_srcq1.STATIC_INPUT_SHA256", pins
        ):
            identities = static_input_identities(payloads)
        self.assertEqual(tuple(identities), STATIC_INPUT_PATHS)
        broken = dict(payloads)
        broken[STATIC_INPUT_PATHS[0]] = b"other"
        with patch(
            "recursive_horizons.fgc.evolution.hlt17_srcq1.STATIC_INPUT_SHA256", pins
        ):
            with self.assertRaises(HLT17SRCQ1IdentityError):
                static_input_identities(broken)

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
            with self.assertRaises(HLT17SRCQ1PublicationError) as failed:
                publish_qualification_result(existing, result)
            self.assertFalse(failed.exception.published)
            missing_parent = root / "missing" / "ns"
            with self.assertRaises(HLT17SRCQ1PublicationError) as absent:
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
            with self.assertRaises(HLT17SRCQ1PublicationError) as posted:
                publish_qualification_result(uncertain, result, fault_hook=fail_after)
            self.assertTrue(posted.exception.published)
            self.assertEqual(posted.exception.outcome, "inconclusive")

    def test_qualify_and_publish_uses_injected_authority_not_a_flag(self) -> None:
        origin, runtime, configs = _world()
        with TemporaryDirectory() as raw, _static_pin(), _passing_c1r1(), patch(
            "recursive_horizons.fgc.evolution.hlt17_srcq1.implementation_identity",
            return_value={
                **implementation_identity(),
                "sha256": "b" * 64,
            },
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
                authority_delta_validator=lambda **kwargs: {
                    "delta": "injected",
                    "authority_commit": kwargs["authority_commit"],
                    "prospective_requested_cap_hex_by_member": _cap_hex(origin, runtime),
                    "output_namespace": str(output),
                },
                runtime_origin=runtime,
                runtime_origin_builder=lambda *a, **k: (_ for _ in ()).throw(
                    AssertionError("factory")
                ),
                rss_bytes=lambda: 1,
                clock=Clock(),
            )
            self.assertEqual(published["authority_commit"], "a" * 40)
            self.assertEqual(published["authority_delta"]["delta"], "injected")
            self.assertTrue((output / "qualification.json").is_file())
            self.assertEqual(
                published["c1r1_actual_versus_imp1_reference_wire"]["distinct"], True
            )

    def test_factory_failure_publishes_a_typed_terminal(self) -> None:
        origin, runtime, _configs = _world()
        with TemporaryDirectory() as raw, _static_pin(), patch(
            "recursive_horizons.fgc.evolution.hlt17_srcq1.implementation_identity",
            return_value={**implementation_identity(), "sha256": "b" * 64},
        ):
            output = Path(raw).resolve() / "failed"
            result = qualify_and_publish(
                origin,
                repository_root=ROOT,
                output_directory=output,
                authority_commit="a" * 40,
                static_input_bytes={path: b"x" for path in STATIC_INPUT_PATHS},
                source_closure_sha256=runtime.source_closure_sha256,
                environment={"python_implementation": "CPython"},
                authority_delta_validator=lambda **kwargs: {
                    "authority_commit": kwargs["authority_commit"],
                    "prospective_requested_cap_hex_by_member": _cap_hex(origin, runtime),
                    "output_namespace": str(output),
                },
                runtime_origin_builder=lambda *args, **kwargs: (_ for _ in ()).throw(
                    RuntimeError("injected factory failure")
                ),
                rss_bytes=lambda: 1,
                clock=Clock(),
            )
            self.assertEqual(result["halted_outcome"], "source")
            self.assertTrue(result["physical_source_qualification_attempted"])
            self.assertFalse(result["accepted_state_advanced"])
            self.assertTrue((output / "qualification.json").is_file())


if __name__ == "__main__":
    unittest.main()
