"""Pure PRO20-EV1 production protocol controls. No store I/O and no live namespace."""

from __future__ import annotations

from dataclasses import dataclass, replace
from hashlib import sha256
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.evidence_io import canonical_json_bytes, load_canonical_json  # noqa: E402
from recursive_horizons.fgc.evolution.boundary_domain import (  # noqa: E402
    BoundaryGeometry,
    CausalBudgetState,
)
from recursive_horizons.fgc.evolution.hlt17_admission_runtime import (  # noqa: E402
    hlt17_implementation_identity,
)
from recursive_horizons.fgc.evolution.hlt17_imp1_bridge import (  # noqa: E402
    prepare_hlt17_next,
    replay_hlt17_rejection,
)
from recursive_horizons.fgc.evolution.hlt17_imp1_cursor import (  # noqa: E402
    HLT17_MEMBER_KEYS,
    HLT17_SYNTHETIC_MEMBER_KEYS,
    encode_hlt17_cursor,
    restore_hlt17_cursor,
)
from recursive_horizons.fgc.evolution.hlt17_member_checkpoint import (  # noqa: E402
    ACCEPTED_FINE,
    HLT17InMemoryCheckpoint,
    origin_predecessor_identity,
)
from recursive_horizons.fgc.evolution.hlt17_runtime_member import (  # noqa: E402
    HLT17ExternalReference,
    HLT17MemberIdentity,
    HLT17RuntimeMember,
    IDENTITY_SCOPE_LIVE,
    IDENTITY_SCOPE_SYNTHETIC,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto19_gr0_static_factory import (  # noqa: E402
    NormalFlowTracers,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    GR0RuntimeStageTransaction,
    GR0UniversalThresholds,
)
from recursive_horizons.fgc.evolution.pro20_rsrc1_protocol import (  # noqa: E402
    EVENT_ORIGIN_HEX,
    EVENT_TARGET_HEX,
    EXACT_DELTA_AUTHORITY_SEAM,
    FREEZE_ARTIFACT_ID,
    INDEPENDENT_BINDER_SEAM,
    NONCLAIM_FLAG_NAMES,
    PINNED_PROTOCOL_V19_PRIVATE_DEPENDENCIES,
    PRODUCTION_NAMESPACE,
    PROTOCOL_ARTIFACT_ID,
    PRO20RSRC1AuthorityReceipt,
    PRO20RSRC1ProtocolError,
    RESOURCE_GATE_SEAM,
    RUNNER_SEAM,
    STORE_KIND_PRODUCTION,
    build_authority_receipt,
    build_root_manifest,
    build_seed_generation,
    build_successor_generation,
    encode_session_close,
    encode_session_open,
    encode_terminal_lock,
    live_member_keys,
    next_required_member_key,
    parse_validated_generation,
    production_cohort_complete,
    validate_authority_receipt,
    validate_root_manifest,
    validate_session_close,
    validate_session_open,
    validate_terminal_lock,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (  # noqa: E402
    TDG6_MINIMUM_MACRO_STEP,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_ledger import (  # noqa: E402
    IMP1_COMPLETE_STATE_CHANNELS,
    IMP1_REJECTION_EVENT_TYPE,
    append_imp1_rejection,
    imp1_checkpoint_extension,
    seed_imp1_ledger,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_qualification import (  # noqa: E402
    SYNTHETIC_START as START,
    SYNTHETIC_TARGET as TARGET,
    make_synthetic_fixture,
    synthetic_diagnostics,
    synthetic_exponential_rhs,
)
from recursive_horizons.fgc.evolution import tdg6_temporal_admission_runtime as tdg6  # noqa: E402


REQUESTED_CAP = 1.0 / 32.0
ARTIFACT_ID = "FGC-1-PRO20-EV1-RSRC1"
CAMPAIGN_ID = "PRO20-RSRC1-PROTOCOL-SYNTHETIC"


@dataclass
class ProvidedGrid:
    coordinates: np.ndarray

    @property
    def point_count(self) -> int:
        return int(self.coordinates.size)

    @property
    def minimum(self) -> float:
        return float(self.coordinates[0])

    @property
    def maximum(self) -> float:
        return float(self.coordinates[-1])

    @property
    def spacing(self) -> float:
        return (self.maximum - self.minimum) / (self.point_count - 1)


@dataclass
class ProvidedInitial:
    grid: ProvidedGrid


class ControllableRHS:
    def __init__(self) -> None:
        self.mode = "ok"
        self.impulse_from = START

    def hlt17_synthetic_callable_binding(self) -> dict[str, object]:
        code = type(self).__call__.__code__
        return {
            "kind": "synthetic_failure_injection_double",
            "type": f"{type(self).__module__}.{type(self).__qualname__}",
            "call_code_sha256": sha256(code.co_code).hexdigest(),
            "mutable_injection_mode_excluded": True,
            "physical_source_authenticated": False,
        }

    def __call__(self, time: float, current) -> EvolutionRHS:
        diagnostics = synthetic_diagnostics()
        zeros = np.zeros_like(current.u)
        if self.mode == "source":
            if time > self.impulse_from:
                diagnostics["source_residual_infinity"] = 1.0
            return EvolutionRHS(zeros, zeros, zeros, diagnostics)
        if self.mode == "cfl":
            if time > self.impulse_from:
                diagnostics["coordinate_speed_upper"] = 1.0e9
            return EvolutionRHS(zeros, zeros, zeros, diagnostics)
        if self.mode == "impulse":
            du = np.zeros_like(current.u)
            if time == self.impulse_from + 1.0 / 128.0:
                du[:, 4] = 1.0
            return EvolutionRHS(du, np.zeros_like(du), np.zeros_like(du), diagnostics)
        return synthetic_exponential_rhs(time, current)


def _production_fixture(method: str, point_count: int):
    coordinates = np.linspace(0.0, 1.0, point_count)
    u = np.zeros((point_count, 6))
    u[:, 0], u[:, 1], u[:, 2], u[:, 3] = 1.0, 0.025, 1.0, coordinates
    u[:, 4] = 0.125
    u[:, 5] = 0.0625
    p, q = np.zeros_like(u), np.zeros_like(u)
    q[:, 3] = 1.0
    state = EvolutionState(u, p, q)
    transaction = GR0RuntimeStageTransaction(
        thresholds=GR0UniversalThresholds(),
        causal_state=CausalBudgetState(
            accepted_time=START, previous_speed_upper=1.0
        ),
        boundary_geometry=BoundaryGeometry(128.0, 24.0, 16.0, 0.375),
        grid_spacing=1.0,
        cfl_maximum=1.0,
    )
    tracers = NormalFlowTracers.create(
        minimum=0.25,
        maximum=0.75,
        spacing=0.25,
        state=state,
        coordinates=coordinates,
        cutoff=1.0,
        outer_radius=128.0,
    )
    inherited = replace(
        tdg6.TDG6TemporalLedger.zero(initial_time=START),
        accumulated_debit_vector=tuple(float(index + 1) / 1024 for index in range(18)),
    )
    ledger = seed_imp1_ledger(
        inherited,
        method=method,
        state_sha256=array_content_sha256(state.u, state.p, state.q),
        step_index=0,
        transaction_serial=0,
        origin_receipt_sha256="1" * 64,
    )
    return state, transaction, tracers, ledger, coordinates


def _live_identity(method: str, point_count: int) -> HLT17MemberIdentity:
    label = "RK4" if method == PRIMARY_METHOD else "SSPRK3"
    return HLT17MemberIdentity(
        scope=IDENTITY_SCOPE_LIVE,
        method_label=label,
        point_count=point_count,
        integrator_id=method,
        spatial_order=4 if label == "RK4" else 2,
        campaign_id=CAMPAIGN_ID,
        amplitude="3",
    )


def _make_live_member(
    method: str,
    point_count: int,
    *,
    rhs=synthetic_exponential_rhs,
    inherited_source: int = 0,
    inherited_cfl: int = 0,
) -> HLT17RuntimeMember:
    state, transaction, tracers, ledger, coordinates = _production_fixture(
        method, point_count
    )
    member = HLT17RuntimeMember(
        identity=_live_identity(method, point_count),
        initial=ProvidedInitial(ProvidedGrid(np.array(coordinates, dtype="<f8", copy=True))),
        operator=rhs,
        projector=None,
        transaction=transaction,
        tracers=tracers,
        state=state,
        temporal_ledger=ledger,
        source_operator_closure=HLT17ExternalReference(
            kind="source_operator_closure", sha256="a" * 64
        ),
        origin_reference=HLT17ExternalReference(
            kind="origin_receipt", sha256=ledger.origin_receipt_sha256
        ),
        time=float(ledger.last_accepted_time),
        step_index=int(ledger.current_step_index),
        transaction_serial=int(ledger.current_transaction_serial),
    )
    member.source_retry_count = inherited_source
    member.CFL_retry_count = inherited_cfl
    member.__post_init__()
    return member


def _receipt(**overrides: str) -> PRO20RSRC1AuthorityReceipt:
    kwargs = {
        "campaign_id": CAMPAIGN_ID,
        "authority_sha256": "11" * 32,
        "implementation_sha256": "22" * 32,
        "config_sha256": "33" * 32,
        "source_sha256": "44" * 32,
        "origin_sha256": "55" * 32,
        "environment_sha256": "66" * 32,
    }
    kwargs.update(overrides)
    return build_authority_receipt(**kwargs)


def _root(receipt: PRO20RSRC1AuthorityReceipt | None = None):
    return build_root_manifest(receipt or _receipt())


def _seed_live(inherited_source: int = 0, inherited_cfl: int = 0, *, controllable: bool = False):
    members = {}
    checkpoints = {}
    specs = (
        (PRIMARY_METHOD, 2049),
        (PRIMARY_METHOD, 4097),
        (PRIMARY_METHOD, 8193),
        (COMPARATOR_METHOD, 4097),
        (COMPARATOR_METHOD, 8193),
        (COMPARATOR_METHOD, 16385),
    )
    for method, point_count in specs:
        rhs = ControllableRHS() if controllable else synthetic_exponential_rhs
        member = _make_live_member(
            method,
            point_count,
            rhs=rhs,
            inherited_source=inherited_source,
            inherited_cfl=inherited_cfl,
        )
        origin = origin_predecessor_identity(
            member,
            origin_cursor_sha256=sha256(member.key.encode("ascii")).hexdigest(),
        )
        checkpoint = HLT17InMemoryCheckpoint.seed(
            member, event_target=TARGET, origin_predecessor=origin
        )
        checkpoints[member.key] = checkpoint
        members[member.key] = checkpoint.generation_bundle()
    root = _root()
    generation = build_seed_generation(root=root, bundles=members)
    return root, generation, checkpoints


def _minimum_step_exhaustion_evidence(bundle) -> dict[str, object]:
    cursor = restore_hlt17_cursor(load_canonical_json(bundle.cursor_bytes))
    current = cursor.ledger
    attempted = float(TDG6_MINIMUM_MACRO_STEP) * 1.5
    record = {
        "event_type": IMP1_REJECTION_EVENT_TYPE,
        "method": current.method,
        "time_hex": current.last_accepted_time.hex(),
        "event_target_hex": cursor.event_target.hex(),
        "attempted_width_hex": attempted.hex(),
        "next_cap_hex": (attempted / 2.0).hex(),
        "accepted_state_sha256": current.current_state_sha256,
        "accepted_step_index": current.current_step_index,
        "accepted_transaction_serial": current.current_transaction_serial,
        "assessment_sha256": "11" * 32,
        "preparation_sha256": "22" * 32,
        "channel_order": list(IMP1_COMPLETE_STATE_CHANNELS),
        "failed_channels": [IMP1_COMPLETE_STATE_CHANNELS[0]],
        "retry_count_for_current_macro_step": current.current_macro_step_temporal_retry_count + 1,
        "cumulative_temporal_retry_count": current.cumulative_temporal_retry_count + 1,
    }
    updated = append_imp1_rejection(current, record)
    mapping = imp1_checkpoint_extension(updated)
    return {
        "kind": "temporal_retry_exhausted",
        "reason": "minimum_macro_step",
        "updated_imp1_ledger": mapping,
        "updated_ledger_sha256": sha256(canonical_json_bytes(mapping)).hexdigest(),
        "executable_retry_cursor": False,
    }


class PRO20RSRC1ProtocolTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.root, cls.seed, cls.checkpoints = _seed_live(
            inherited_source=3,
            inherited_cfl=1,
            controllable=True,
        )

    def test_receipt_and_root_bind_identities_without_issuing_authority(self) -> None:
        receipt = _receipt()
        self.assertEqual(receipt.mapping["protocol_artifact_id"], PROTOCOL_ARTIFACT_ID)
        self.assertEqual(receipt.mapping["freeze_artifact_id"], FREEZE_ARTIFACT_ID)
        self.assertEqual(receipt.mapping["namespace"], PRODUCTION_NAMESPACE)
        self.assertEqual(tuple(receipt.mapping["member_keys"]), HLT17_MEMBER_KEYS)
        self.assertEqual(receipt.mapping["event_origin_time_hex"], EVENT_ORIGIN_HEX)
        self.assertEqual(receipt.mapping["event_target_time_hex"], EVENT_TARGET_HEX)
        self.assertEqual(receipt.mapping["implementation_identity"], hlt17_implementation_identity())
        for name in NONCLAIM_FLAG_NAMES:
            self.assertIs(receipt.mapping[name], False)
        self.assertEqual(receipt.mapping["runner_seam"], RUNNER_SEAM)
        self.assertEqual(receipt.mapping["exact_delta_authority_seam"], EXACT_DELTA_AUTHORITY_SEAM)
        self.assertEqual(receipt.mapping["resource_gate_seam"], RESOURCE_GATE_SEAM)
        self.assertEqual(receipt.mapping["independent_binder_seam"], INDEPENDENT_BINDER_SEAM)
        root = build_root_manifest(receipt)
        self.assertEqual(root["store_kind"], STORE_KIND_PRODUCTION)
        self.assertEqual(root["authority_receipt_sha256"], receipt.sha256)
        self.assertEqual(root["origin_sha256"], "55" * 32)
        self.assertNotIn("runner_sha256", root)
        self.assertFalse(root["campaign_execution_authorized"])
        self.assertFalse(root["production_write_authorized"])
        self.assertFalse(root["live_runner_authorized"])
        self.assertEqual(tuple(live_member_keys(root["member_keys"])), tuple(root["member_keys"]))
        rebuilt = validate_authority_receipt(dict(receipt.mapping))
        self.assertEqual(rebuilt["authority_sha256"], receipt.mapping["authority_sha256"])
        self.assertEqual(
            PINNED_PROTOCOL_V19_PRIVATE_DEPENDENCIES,
            (
                "protocol_v19._attempt_payload",
                "protocol_v19._blob_path",
                "protocol_v19._cohort_complete",
                "protocol_v19._cursor_from_bundle",
                "protocol_v19._descriptor_mapping",
                "protocol_v19._new_blobs_for_cohort",
                "protocol_v19._require_fine_transition",
                "protocol_v19._require_retry_transition",
                "protocol_v19._require_unchanged_attempt",
            ),
        )

    def test_session_records_round_trip_and_refuse_claim_promotion(self) -> None:
        opened = encode_session_open(
            campaign_id=CAMPAIGN_ID,
            session_index=0,
            owner_token="owner",
            host="synthetic-host",
            pid=17,
            thread_id=19,
            exclusion_lock_sha256="77" * 32,
            authority_receipt_sha256="88" * 32,
        )
        self.assertEqual(validate_session_open(opened), opened)
        promoted_open = dict(opened)
        promoted_open["production_write_authorized"] = True
        with self.assertRaisesRegex(PRO20RSRC1ProtocolError, "session open"):
            validate_session_open(promoted_open)
        closed = encode_session_close(
            campaign_id=CAMPAIGN_ID,
            session_index=0,
            session_open_sha256="99" * 32,
            owner_token="owner",
            authority_receipt_sha256="88" * 32,
        )
        self.assertEqual(validate_session_close(closed), closed)
        promoted_close = dict(closed)
        promoted_close["freeze_artifact_id"] = "FGC-1-PRO20-EV1-FOREIGN"
        with self.assertRaisesRegex(PRO20RSRC1ProtocolError, "session close"):
            validate_session_close(promoted_close)

    def test_generation_manifest_nonclaims_are_authenticated(self) -> None:
        for name in (
            "production_write_authorized",
            "campaign_execution_authorized",
            "calibration_claimed",
            "eligibility_claimed",
            "physical_classification",
        ):
            files = dict(self.seed.files)
            manifest = load_canonical_json(files["generation_manifest.json"])
            manifest[name] = True
            files["generation_manifest.json"] = canonical_json_bytes(manifest)
            with self.subTest(flag=name), self.assertRaisesRegex(
                PRO20RSRC1ProtocolError, name
            ):
                parse_validated_generation(
                    files,
                    root=self.root,
                    predecessor=None,
                    predecessor_bundles=None,
                )

    def test_altered_identities_and_claim_promotion_are_refused(self) -> None:
        receipt = _receipt()
        mutated = dict(receipt.mapping)
        mutated["origin_sha256"] = "77" * 32
        with self.assertRaisesRegex(PRO20RSRC1ProtocolError, "authority_receipt_sha256"):
            validate_root_manifest(
                {
                    **build_root_manifest(receipt),
                    "origin_sha256": "77" * 32,
                }
            )
        for name in ("authority_sha256", "source_sha256", "environment_sha256", "config_sha256"):
            with self.subTest(hash_name=name):
                with self.assertRaisesRegex(PRO20RSRC1ProtocolError, "distinct|receipt"):
                    mapping = dict(build_root_manifest(receipt))
                    mapping[name] = "aa" * 32
                    validate_root_manifest(mapping)
        wrong_identity = dict(hlt17_implementation_identity())
        wrong_identity["implementation_id"] = "not-c1r1"
        with self.assertRaisesRegex(PRO20RSRC1ProtocolError, "implementation_identity"):
            build_authority_receipt(
                campaign_id=CAMPAIGN_ID,
                authority_sha256="11" * 32,
                implementation_sha256="22" * 32,
                config_sha256="33" * 32,
                source_sha256="44" * 32,
                origin_sha256="55" * 32,
                environment_sha256="66" * 32,
                implementation_identity=wrong_identity,
            )
        for name in NONCLAIM_FLAG_NAMES:
            promoted = dict(receipt.mapping)
            promoted[name] = True
            with self.subTest(flag=name):
                with self.assertRaisesRegex(PRO20RSRC1ProtocolError, name):
                    validate_authority_receipt(promoted)
        mutated_ns = dict(receipt.mapping)
        mutated_ns["namespace"] = "runs/fgc-2-sf1/proto19/calibration"
        with self.assertRaisesRegex(PRO20RSRC1ProtocolError, "namespace"):
            validate_authority_receipt(mutated_ns)
        mutated_order = dict(receipt.mapping)
        mutated_order["member_keys"] = list(reversed(HLT17_MEMBER_KEYS))
        with self.assertRaisesRegex(PRO20RSRC1ProtocolError, "canonical order"):
            validate_authority_receipt(mutated_order)

    def test_seed_requires_exact_six_live_c1r1_members(self) -> None:
        root, generation, checkpoints = self.root, self.seed, self.checkpoints
        self.assertEqual(generation.kind, "seed")
        self.assertEqual(tuple(generation.bundles), HLT17_MEMBER_KEYS)
        self.assertEqual(next_required_member_key(generation), "RK4-2049")
        self.assertFalse(generation.terminal)
        self.assertFalse(generation.first_event_complete)
        checkpoint = generation.checkpoint_mapping
        self.assertFalse(checkpoint["calibration_claimed"])
        self.assertFalse(checkpoint["sgbl_claimed"])
        self.assertFalse(checkpoint["def1_claimed"])
        self.assertFalse(checkpoint["candidate_claimed"])
        self.assertFalse(checkpoint["mechanism_claimed"])
        self.assertFalse(checkpoint["physics_claimed"])
        self.assertFalse(checkpoint["calibration_eligible"])
        for key in HLT17_MEMBER_KEYS:
            state = checkpoint["members"][key]
            self.assertEqual(state["source_total"], 3)
            self.assertEqual(state["cfl_total"], 1)
            self.assertEqual(state["source_current"], 0)
            self.assertEqual(state["accepted_generation"], 0)
            self.assertEqual(state["accepted_time_hex"], START.hex())
            self.assertEqual(state["implementation_id"], hlt17_implementation_identity()["implementation_id"])
            self.assertEqual(checkpoints[key].member.identity.scope, IDENTITY_SCOPE_LIVE)
        synthetic = {}
        for method in (PRIMARY_METHOD, COMPARATOR_METHOD):
            state, transaction, tracers, ledger, coordinates = make_synthetic_fixture(method)
            label = "RK4" if method == PRIMARY_METHOD else "SSPRK3"
            member = HLT17RuntimeMember(
                identity=HLT17MemberIdentity(
                    scope=IDENTITY_SCOPE_SYNTHETIC,
                    method_label=label,
                    point_count=9,
                    integrator_id=method,
                    spatial_order=4 if label == "RK4" else 2,
                    campaign_id=CAMPAIGN_ID,
                    amplitude="synthetic",
                ),
                initial=ProvidedInitial(ProvidedGrid(np.array(coordinates, dtype="<f8", copy=True))),
                operator=synthetic_exponential_rhs,
                projector=None,
                transaction=transaction,
                tracers=tracers,
                state=state,
                temporal_ledger=ledger,
                source_operator_closure=HLT17ExternalReference(
                    kind="source_operator_closure", sha256="a" * 64
                ),
                origin_reference=HLT17ExternalReference(
                    kind="origin_receipt", sha256=ledger.origin_receipt_sha256
                ),
                time=float(ledger.last_accepted_time),
            )
            checkpoint_obj = HLT17InMemoryCheckpoint.seed(member, event_target=TARGET)
            synthetic[member.key] = checkpoint_obj.generation_bundle()
        with self.assertRaisesRegex(PRO20RSRC1ProtocolError, "six live|synthetic|cohort"):
            build_seed_generation(root=root, bundles=synthetic)
        with self.assertRaisesRegex(PRO20RSRC1ProtocolError, "canonical order|six live"):
            live_member_keys(list(HLT17_SYNTHETIC_MEMBER_KEYS))

    def test_canonical_order_and_pending_retry_ownership(self) -> None:
        root, seed, checkpoints = self.root, self.seed, self.checkpoints
        first = checkpoints["RK4-2049"].advance_once(requested_cap=REQUESTED_CAP)
        self.assertEqual(first.disposition, ACCEPTED_FINE)
        with self.assertRaisesRegex(PRO20RSRC1ProtocolError, "canonical next live member"):
            build_successor_generation(
                root=root,
                predecessor=seed,
                predecessor_bundles=seed.bundles,
                kind="accepted_fine",
                member_key="RK4-4097",
                successor_bundle=seed.bundles["RK4-4097"],
            )
        accepted = build_successor_generation(
            root=root,
            predecessor=seed,
            predecessor_bundles=seed.bundles,
            kind="accepted_fine",
            member_key="RK4-2049",
            successor_bundle=first.bundle,
        )
        self.assertEqual(accepted.changed_member_key, "RK4-2049")
        self.assertFalse(accepted.checkpoint_mapping["campaign_execution_authorized"])
        rk4 = checkpoints["RK4-2049"]
        rhs = rk4.member.operator
        assert isinstance(rhs, ControllableRHS)
        rhs.impulse_from = rk4.member.time
        rhs.mode = "impulse"
        temporal = rk4.advance_once(requested_cap=REQUESTED_CAP)
        self.assertEqual(temporal.disposition, "temporal_retry_required")
        pending = build_successor_generation(
            root=root,
            predecessor=accepted,
            predecessor_bundles=accepted.bundles,
            kind="temporal_retry",
            member_key="RK4-2049",
            successor_bundle=temporal.bundle,
        )
        self.assertEqual(next_required_member_key(pending), "RK4-2049")
        with self.assertRaisesRegex(PRO20RSRC1ProtocolError, "canonical next live member"):
            build_successor_generation(
                root=root,
                predecessor=pending,
                predecessor_bundles=pending.bundles,
                kind="accepted_fine",
                member_key="RK4-4097",
                successor_bundle=seed.bundles["RK4-4097"],
            )
        rk4.restore(rk4.accepted_encoding, encode_hlt17_cursor(rk4.cursor))
        rhs.impulse_from = rk4.member.time
        rhs.mode = "impulse"
        replay = replay_hlt17_rejection(rk4.cursor, rk4.live_boundary())
        rhs.mode = "source"
        source = rk4.absorb_nonfine(
            prepare_hlt17_next(rk4.cursor, rk4.live_boundary(), replay=replay)
        )
        self.assertEqual(source.disposition, "source_retry_required")
        sourced = build_successor_generation(
            root=root,
            predecessor=pending,
            predecessor_bundles=pending.bundles,
            kind="source_retry",
            member_key="RK4-2049",
            successor_bundle=source.bundle,
        )
        state = sourced.checkpoint_mapping["members"]["RK4-2049"]
        self.assertEqual(state["source_current"], 1)
        self.assertEqual(state["source_total"], 4)
        self.assertEqual(
            state["physical_state_sha256"],
            accepted.checkpoint_mapping["members"]["RK4-2049"]["physical_state_sha256"],
        )
        self.assertEqual(
            state["tdg6_snapshot_sha256"],
            accepted.checkpoint_mapping["members"]["RK4-2049"]["tdg6_snapshot_sha256"],
        )
        self.assertTrue(sourced.journal_mapping["payload"]["accepted_state_preserved"])

    def test_typed_terminal_and_first_event_completion_stay_nonclaims(self) -> None:
        root, seed = self.root, self.seed
        evidence = _minimum_step_exhaustion_evidence(seed.bundles["RK4-2049"])
        terminal = build_successor_generation(
            root=root,
            predecessor=seed,
            predecessor_bundles=seed.bundles,
            kind="temporal_retry_exhausted",
            member_key="RK4-2049",
            successor_bundle=seed.bundles["RK4-2049"],
            terminal_evidence=evidence,
        )
        self.assertTrue(terminal.terminal)
        self.assertEqual(terminal.disposition, "temporal_retry_exhausted")
        self.assertFalse(terminal.first_event_complete)
        self.assertFalse(terminal.checkpoint_mapping["physics_claimed"])
        with self.assertRaisesRegex(PRO20RSRC1ProtocolError, "forbids further appends"):
            build_successor_generation(
                root=root,
                predecessor=terminal,
                predecessor_bundles=terminal.bundles,
                kind="accepted_fine",
                member_key="RK4-2049",
                successor_bundle=seed.bundles["RK4-2049"],
            )
        states = {
            key: dict(seed.checkpoint_mapping["members"][key])
            for key in HLT17_MEMBER_KEYS
        }
        self.assertFalse(production_cohort_complete(states))
        for state in states.values():
            state["accepted_time_hex"] = EVENT_TARGET_HEX
            state["mode"] = "FRESH_READY"
            state["latest_owner"] = None
            state["source_current"] = 0
            state["cfl_current"] = 0
            state["parked"] = True
            state["public_fresh"] = True
        self.assertTrue(production_cohort_complete(states))
        with self.assertRaisesRegex(PRO20RSRC1ProtocolError, "canonical six"):
            production_cohort_complete(dict(reversed(tuple(states.items()))))
        lock = encode_terminal_lock(
            campaign_id=CAMPAIGN_ID,
            checkpoint_sha256="77" * 32,
            journal_record_sha256="88" * 32,
            store_generation=7,
            disposition="first_event_complete",
            changed_member_key="SSPRK3-16385",
            first_event_complete=True,
        )
        validated_lock = validate_terminal_lock(lock)
        self.assertTrue(validated_lock["first_event_complete"])
        self.assertFalse(validated_lock["calibration_eligible"])
        self.assertFalse(validated_lock["physics_claimed"])
        promoted_lock = dict(lock)
        promoted_lock["physics_claimed"] = True
        with self.assertRaisesRegex(PRO20RSRC1ProtocolError, "terminal lock"):
            validate_terminal_lock(promoted_lock)

    def test_wrong_live_method_and_foreign_origin_are_refused(self) -> None:
        root, seed, checkpoints = self.root, self.seed, self.checkpoints
        with self.assertRaisesRegex(PRO20RSRC1ProtocolError, "canonical next|outside"):
            build_successor_generation(
                root=root,
                predecessor=seed,
                predecessor_bundles=seed.bundles,
                kind="accepted_fine",
                member_key="RK4-9",
                successor_bundle=seed.bundles["RK4-2049"],
            )
        descriptor = load_canonical_json(seed.bundles["RK4-2049"].descriptor)
        predecessor = dict(descriptor["predecessor_cursor_identity"])
        predecessor["protocol_artifact_id"] = "-".join(
            ("FGC", "2", "SF1", "PROTO18")
        )
        descriptor["predecessor_cursor_identity"] = predecessor
        with self.assertRaises(Exception):
            from recursive_horizons.fgc.evolution.hlt17_member_checkpoint import (
                HLT17GenerationBundle,
            )

            HLT17GenerationBundle(
                canonical_json_bytes(descriptor),
                seed.bundles["RK4-2049"].payload,
                seed.bundles["RK4-2049"].cursor_bytes,
            )
        del checkpoints


if __name__ == "__main__":
    unittest.main()
