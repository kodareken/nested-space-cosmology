"""PRO20-EV1 production store publication, exclusion, and adversarial controls."""

from __future__ import annotations

from dataclasses import dataclass, replace
from hashlib import sha256
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.evidence_io import (  # noqa: E402
    canonical_json_bytes,
    load_canonical_json,
)
from recursive_horizons.fgc.evolution.boundary_domain import (  # noqa: E402
    BoundaryGeometry,
    CausalBudgetState,
)
from recursive_horizons.fgc.evolution.hlt17_imp1_cursor import HLT17_MEMBER_KEYS  # noqa: E402
from recursive_horizons.fgc.evolution.hlt17_member_checkpoint import (  # noqa: E402
    ACCEPTED_FINE,
    HLT17InMemoryCheckpoint,
    origin_predecessor_identity,
)
from recursive_horizons.fgc.evolution.hlt17_member_codec import HLT17MemberEncoding  # noqa: E402
from recursive_horizons.fgc.evolution.hlt17_runtime_member import (  # noqa: E402
    HLT17ExternalReference,
    HLT17MemberIdentity,
    HLT17RuntimeMember,
    IDENTITY_SCOPE_LIVE,
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
from recursive_horizons.fgc.evolution.pro20_ev1_protocol import (  # noqa: E402
    FREEZE_ARTIFACT_ID,
    PRODUCTION_NAMESPACE,
    PROTOCOL_ARTIFACT_ID,
    PRO20EV1ProtocolError,
    build_authority_receipt,
    generation_directory_name,
)
from recursive_horizons.fgc.evolution.pro20_ev1_store import (  # noqa: E402
    PRO20EV1CampaignStore,
    PRO20EV1CampaignStoreError,
    PRO20EV1PostpublicationUncertainty,
    PRO20EV1TerminalStore,
    PRO20EV1UnclosedSession,
    PRO20EV1WriterConflict,
    ROOT_MANIFEST_NAME,
    TERMINAL_LOCK_NAME,
    WRITER_EXCLUSION_NAME,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_ledger import seed_imp1_ledger  # noqa: E402
from recursive_horizons.fgc.evolution.tdg11_imp1_qualification import (  # noqa: E402
    SYNTHETIC_START as START,
    SYNTHETIC_TARGET as TARGET,
    synthetic_diagnostics,
    synthetic_exponential_rhs,
)
from recursive_horizons.fgc.evolution import tdg6_temporal_admission_runtime as tdg6  # noqa: E402


REQUESTED_CAP = 1.0 / 32.0
ARTIFACT_ID = "FGC-1-PRO20-EV1-FRZ1"
CAMPAIGN_ID = "PRO20-EV1-STORE-SYNTHETIC"
REAL_NAMESPACE = ROOT / PRODUCTION_NAMESPACE


def _real_namespace_identity() -> tuple[tuple[object, ...], ...] | None:
    """Snapshot canonical namespace metadata without assuming pre- or post-run state."""

    if not REAL_NAMESPACE.exists():
        return None
    paths = [REAL_NAMESPACE, *REAL_NAMESPACE.rglob("*")]
    return tuple(
        (
            "." if path == REAL_NAMESPACE else path.relative_to(REAL_NAMESPACE).as_posix(),
            metadata.st_dev,
            metadata.st_ino,
            metadata.st_mode,
            metadata.st_size,
            metadata.st_mtime_ns,
        )
        for path in sorted(paths, key=lambda item: item.as_posix())
        for metadata in (path.lstat(),)
    )


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
        self.calls = 0

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
        self.calls += 1
        diagnostics = synthetic_diagnostics()
        zeros = np.zeros_like(current.u)
        if self.mode == "source":
            if time > self.impulse_from:
                diagnostics["source_residual_infinity"] = 1.0
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
        causal_state=CausalBudgetState(accepted_time=START, previous_speed_upper=1.0),
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


def _make_live_member(method: str, point_count: int, *, rhs=None, inherited_source: int = 0, inherited_cfl: int = 0):
    operator = synthetic_exponential_rhs if rhs is None else rhs
    state, transaction, tracers, ledger, coordinates = _production_fixture(method, point_count)
    label = "RK4" if method == PRIMARY_METHOD else "SSPRK3"
    member = HLT17RuntimeMember(
        identity=HLT17MemberIdentity(
            scope=IDENTITY_SCOPE_LIVE,
            method_label=label,
            point_count=point_count,
            integrator_id=method,
            spatial_order=4 if label == "RK4" else 2,
            campaign_id=CAMPAIGN_ID,
            amplitude="3",
        ),
        initial=ProvidedInitial(ProvidedGrid(np.array(coordinates, dtype="<f8", copy=True))),
        operator=operator,
        projector=None,
        transaction=transaction,
        tracers=tracers,
        state=state,
        temporal_ledger=ledger,
        source_operator_closure=HLT17ExternalReference(
            kind="source_operator_closure", sha256="b" * 64
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


def _live_cohort(*, controllable: bool = False, inherited_source: int = 0, inherited_cfl: int = 0):
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
        rhs = ControllableRHS() if controllable else None
        member = _make_live_member(
            method,
            point_count,
            rhs=rhs,
            inherited_source=inherited_source,
            inherited_cfl=inherited_cfl,
        )
        origin = origin_predecessor_identity(
            member, origin_cursor_sha256=sha256(member.key.encode("ascii")).hexdigest()
        )
        checkpoints[member.key] = HLT17InMemoryCheckpoint.seed(
            member, event_target=TARGET, origin_predecessor=origin
        )
    bundles = {key: item.generation_bundle() for key, item in checkpoints.items()}
    return checkpoints, bundles


def _receipt():
    return build_authority_receipt(
        campaign_id=CAMPAIGN_ID,
        authority_sha256="11" * 32,
        implementation_sha256="22" * 32,
        config_sha256="33" * 32,
        source_sha256="44" * 32,
        origin_sha256="55" * 32,
        environment_sha256="66" * 32,
    )


def _temp_repo():
    raw = tempfile.TemporaryDirectory()
    repo = Path(raw.name).resolve()
    (repo / "runs" / "fgc-2-sf1").mkdir(parents=True)
    return raw, repo


CHILD = r"""
import sys
from pathlib import Path
sys.path[:0] = sys.argv[1].split(":")
from recursive_horizons.fgc.evolution.pro20_ev1_store import PRO20EV1CampaignStore
root = Path(sys.argv[2])
action = sys.argv[3]
store = PRO20EV1CampaignStore.open(root)
try:
    capability = store.acquire_writer(owner_token="child", host="synthetic-host")
except Exception as exc:
    sys.stdout.write(type(exc).__name__ + ":" + str(exc))
    raise SystemExit(0)
sys.stdout.write("ACQUIRED")
sys.stdout.flush()
if action == "hold":
    sys.stdin.read()
    store.close_writer(capability)
"""


class PRO20EV1CampaignStoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        checkpoints, cls.bundles = _live_cohort(
            controllable=True,
            inherited_source=2,
            inherited_cfl=1,
        )
        accepted = checkpoints["RK4-2049"].advance_once(
            requested_cap=REQUESTED_CAP
        )
        if accepted.disposition != ACCEPTED_FINE:
            raise AssertionError("bounded RK4-2049 fixture did not accept one fine step")
        cls.accepted_rk4_bundle = accepted.bundle

    def test_atomic_seed_namespace_and_read_only_authentication(self) -> None:
        real_before = _real_namespace_identity()
        holder, repo = _temp_repo()
        with holder:
            store = PRO20EV1CampaignStore.publish_seed_store(
                repo, receipt=_receipt(), bundles=self.bundles
            )
            self.assertEqual(store.root, repo / PRODUCTION_NAMESPACE)
            self.assertTrue((store.root / ROOT_MANIFEST_NAME).is_file())
            self.assertTrue((store.root / WRITER_EXCLUSION_NAME).is_file())
            self.assertTrue((store.root / generation_directory_name(0)).is_dir())
            self.assertEqual(_real_namespace_identity(), real_before)
            view = PRO20EV1CampaignStore.authenticate_read_only(repo)
            self.assertEqual(len(view.generations), 1)
            self.assertEqual(view.generations[0].kind, "seed")
            self.assertEqual(tuple(view.bundles), HLT17_MEMBER_KEYS)
            self.assertFalse(view.terminal)
            self.assertFalse(view.unclosed_session)
            self.assertEqual(view.root_manifest["protocol_artifact_id"], PROTOCOL_ARTIFACT_ID)
            self.assertEqual(view.root_manifest["freeze_artifact_id"], FREEZE_ARTIFACT_ID)
            self.assertFalse(view.root_manifest["campaign_execution_authorized"])
            self.assertIn("runner", store.remaining_seams)
            for key in HLT17_MEMBER_KEYS:
                state = view.generations[0].checkpoint_mapping["members"][key]
                self.assertEqual(state["source_total"], 2)
                self.assertEqual(state["cfl_total"], 1)
            capability = store.acquire_writer(owner_token="owner", host="synthetic-host")
            latest = store.publish_attempt(
                capability,
                kind="accepted_fine",
                member_key="RK4-2049",
                successor_bundle=self.accepted_rk4_bundle,
            )
            self.assertEqual(latest.changed_member_key, "RK4-2049")
            with self.assertRaisesRegex(PRO20EV1CampaignStoreError, "canonical next"):
                store.publish_attempt(
                    capability,
                    kind="accepted_fine",
                    member_key="RK4-4097",
                    successor_bundle=self.bundles["RK4-4097"],
                )
            store.close_writer(capability)
            restored = PRO20EV1CampaignStore.authenticate_read_only(repo)
            encoding = HLT17MemberEncoding(
                restored.bundles["RK4-2049"].descriptor,
                restored.bundles["RK4-2049"].payload,
            )
            self.assertEqual(
                encoding.decoded.metadata["runtime_identity"]["member_key"], "RK4-2049"
            )
        self.assertEqual(_real_namespace_identity(), real_before)

    def test_thread_process_exclusion_and_refused_session_takeover(self) -> None:
        holder, repo = _temp_repo()
        with holder:
            store = PRO20EV1CampaignStore.publish_seed_store(
                repo, receipt=_receipt(), bundles=self.bundles
            )
            capability = store.acquire_writer(owner_token="owner", host="synthetic-host")
            errors: list[BaseException] = []

            def rival() -> None:
                try:
                    store.acquire_writer(owner_token="rival", host="synthetic-host")
                except BaseException as exc:
                    errors.append(exc)

            thread = threading.Thread(target=rival)
            thread.start()
            thread.join()
            self.assertTrue(errors)
            self.assertIsInstance(errors[0], PRO20EV1WriterConflict)
            env = os.environ.copy()
            proc = subprocess.run(
                [
                    sys.executable,
                    "-c",
                    CHILD,
                    str(ROOT) + ":" + str(ROOT / "src"),
                    str(repo),
                    "exit",
                ],
                check=False,
                capture_output=True,
                text=True,
                env=env,
            )
            self.assertIn("PRO20EV1WriterConflict", proc.stdout)
            store.close_writer(capability)
            child = subprocess.Popen(
                [
                    sys.executable,
                    "-c",
                    CHILD,
                    str(ROOT) + ":" + str(ROOT / "src"),
                    str(repo),
                    "hold",
                ],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                env=env,
            )
            assert child.stdout is not None
            status = child.stdout.read(8)
            self.assertEqual(status, "ACQUIRED")
            with self.assertRaises(PRO20EV1WriterConflict):
                store.acquire_writer(owner_token="parent", host="synthetic-host")
            child.kill()
            child.wait()
            if child.stdin is not None:
                child.stdin.close()
            if child.stdout is not None:
                child.stdout.close()
            if child.stderr is not None:
                child.stderr.close()
            with self.assertRaises(PRO20EV1UnclosedSession):
                store.acquire_writer(owner_token="parent", host="synthetic-host")

    def test_publication_faults_symlink_hardlink_partial_and_fork(self) -> None:
        holder, repo = _temp_repo()
        with holder:
            def before(cut: str) -> None:
                if cut == "before_publish":
                    raise OSError("injected before_publish")

            with self.assertRaises(PRO20EV1CampaignStoreError):
                PRO20EV1CampaignStore.publish_seed_store(
                    repo,
                    receipt=_receipt(),
                    bundles=self.bundles,
                    _fault_hook=before,
                )
            self.assertFalse((repo / PRODUCTION_NAMESPACE).exists())

            def after(cut: str) -> None:
                if cut == "after_publish":
                    raise OSError("injected after_publish")

            with self.assertRaises(PRO20EV1PostpublicationUncertainty) as failed:
                PRO20EV1CampaignStore.publish_seed_store(
                    repo,
                    receipt=_receipt(),
                    bundles=self.bundles,
                    _fault_hook=after,
                )
            self.assertTrue(failed.exception.published)
            self.assertTrue((repo / PRODUCTION_NAMESPACE).is_dir())
            view = PRO20EV1CampaignStore.authenticate_read_only(repo)
            self.assertEqual(view.generations[0].kind, "seed")
            store = PRO20EV1CampaignStore.open(repo)
            capability = store.acquire_writer(owner_token="owner", host="synthetic-host")

            def after_generation(cut: str) -> None:
                if cut == "after_publish":
                    raise OSError("injected generation after_publish")

            with self.assertRaises(PRO20EV1PostpublicationUncertainty) as gen_failed:
                store.publish_attempt(
                    capability,
                    kind="resource_exhausted",
                    member_key="RK4-2049",
                    successor_bundle=self.bundles["RK4-2049"],
                    terminal_evidence={
                        "kind": "resource_exhausted",
                        "reason": "typed terminal after seed",
                        "executable_retry_cursor": False,
                    },
                    _fault_hook=after_generation,
                )
            self.assertTrue(gen_failed.exception.published)
            self.assertTrue((store.root / generation_directory_name(1)).exists())
            with self.assertRaises(PRO20EV1PostpublicationUncertainty):
                store.close_writer(capability)
            store.abandon_writer(capability)
            with self.assertRaises(PRO20EV1UnclosedSession):
                store.acquire_writer(owner_token="later", host="synthetic-host")

        holder, repo = _temp_repo()
        with holder:
            store = PRO20EV1CampaignStore.publish_seed_store(
                repo, receipt=_receipt(), bundles=self.bundles
            )
            capability = store.acquire_writer(owner_token="owner", host="synthetic-host")
            store.publish_attempt(
                capability,
                kind="accepted_fine",
                member_key="RK4-2049",
                successor_bundle=self.accepted_rk4_bundle,
            )
            store.close_writer(capability)
            generation0 = store.root / generation_directory_name(0)
            (generation0 / "extra.json").write_bytes(b"{}")
            with self.assertRaises(PRO20EV1CampaignStoreError):
                PRO20EV1CampaignStore.authenticate_read_only(repo)
            (generation0 / "extra.json").unlink()
            hole = store.root / generation_directory_name(2)
            hole.mkdir()
            (hole / "journal.json").write_bytes((generation0 / "journal.json").read_bytes())
            (hole / "checkpoint.json").write_bytes((generation0 / "checkpoint.json").read_bytes())
            (hole / "generation_manifest.json").write_bytes(
                (generation0 / "generation_manifest.json").read_bytes()
            )
            with self.assertRaises(PRO20EV1CampaignStoreError):
                PRO20EV1CampaignStore.authenticate_read_only(repo)
            hole.joinpath("journal.json").unlink()
            hole.joinpath("checkpoint.json").unlink()
            hole.joinpath("generation_manifest.json").unlink()
            hole.rmdir()
            target = generation0 / "journal.json"
            alias = store.root / "link-journal.json"
            os.symlink(target, alias)
            with self.assertRaises(PRO20EV1CampaignStoreError):
                PRO20EV1CampaignStore.authenticate_read_only(repo)
            alias.unlink()
            linked = generation0 / "hard-journal.json"
            os.link(target, linked)
            with self.assertRaises(PRO20EV1CampaignStoreError):
                PRO20EV1CampaignStore.authenticate_read_only(repo)
            linked.unlink()

    def test_typed_terminal_locks_and_refuses_repair(self) -> None:
        holder, repo = _temp_repo()
        with holder:
            store = PRO20EV1CampaignStore.publish_seed_store(
                repo, receipt=_receipt(), bundles=self.bundles
            )
            capability = store.acquire_writer(owner_token="owner", host="synthetic-host")
            seed_view = store.authenticate()
            terminal = store.publish_attempt(
                capability,
                kind="resource_exhausted",
                member_key="RK4-2049",
                successor_bundle=seed_view.bundles["RK4-2049"],
                terminal_evidence={
                    "kind": "resource_exhausted",
                    "reason": "synthetic typed terminal",
                    "executable_retry_cursor": False,
                },
            )
            self.assertTrue(terminal.terminal)
            self.assertTrue((store.root / TERMINAL_LOCK_NAME).is_file())
            self.assertFalse(terminal.checkpoint_mapping["calibration_eligible"])
            with self.assertRaises(PRO20EV1TerminalStore):
                store.publish_attempt(
                    capability,
                    kind="accepted_fine",
                    member_key="RK4-2049",
                    successor_bundle=self.bundles["RK4-2049"],
                )
            store.close_writer(capability)
            with self.assertRaises(PRO20EV1TerminalStore):
                store.acquire_writer(owner_token="later", host="synthetic-host")
            view = PRO20EV1CampaignStore.authenticate_read_only(repo)
            self.assertTrue(view.terminal)
            self.assertFalse(view.first_event_complete)

    def test_namespace_symlink_and_altered_root_hashes(self) -> None:
        real_before = _real_namespace_identity()
        holder, repo = _temp_repo()
        with holder:
            with self.assertRaises(PRO20EV1CampaignStoreError):
                PRO20EV1CampaignStore.publish_seed_store(
                    repo / "missing", receipt=_receipt(), bundles=self.bundles
                )
            nested = repo / "nested"
            nested.mkdir()
            os.symlink(repo / "runs", nested / "runs")
            with self.assertRaises((PRO20EV1CampaignStoreError, PRO20EV1ProtocolError, Exception)):
                PRO20EV1CampaignStore.publish_seed_store(
                    nested, receipt=_receipt(), bundles=self.bundles
                )
            store = PRO20EV1CampaignStore.publish_seed_store(
                repo, receipt=_receipt(), bundles=self.bundles
            )
            mutated = load_canonical_json((store.root / ROOT_MANIFEST_NAME).read_bytes())
            mutated["origin_sha256"] = "77" * 32
            (store.root / ROOT_MANIFEST_NAME).write_bytes(canonical_json_bytes(mutated))
            with self.assertRaises(PRO20EV1CampaignStoreError):
                PRO20EV1CampaignStore.authenticate_read_only(repo)
            with self.assertRaises(PRO20EV1CampaignStoreError):
                store.authenticate()
        self.assertEqual(_real_namespace_identity(), real_before)
        self.assertNotEqual(Path.cwd(), REAL_NAMESPACE)


if __name__ == "__main__":
    unittest.main()
