"""Synthetic HLT17/PROTO19 campaign-store publication, fault, and exclusion controls."""

from __future__ import annotations

from dataclasses import dataclass, replace
from hashlib import sha256
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
from unittest.mock import patch
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.evidence_io import (  # noqa: E402
    PrepublicationError,
    canonical_json_bytes,
    load_canonical_json,
    publish_exclusive_file,
    read_regular_file,
)
from recursive_horizons.fgc.evolution.hlt17_campaign_store import (  # noqa: E402
    HLT17CampaignStore,
    HLT17CampaignStoreError,
    HLT17PostpublicationUncertainty,
    HLT17TerminalStore,
    HLT17UnclosedSession,
    HLT17WriterConflict,
    PRODUCTION_EVENT_PATH,
    ROOT_MANIFEST_NAME,
    TERMINAL_LOCK_NAME,
    WRITER_EXCLUSION_NAME,
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
)
from recursive_horizons.fgc.evolution.hlt17_member_codec import HLT17MemberEncoding  # noqa: E402
from recursive_horizons.fgc.evolution.hlt17_runtime_member import (  # noqa: E402
    HLT17ExternalReference,
    HLT17MemberIdentity,
    HLT17RuntimeMember,
    IDENTITY_SCOPE_SYNTHETIC,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
)
from recursive_horizons.fgc.evolution.protocol_v19 import (  # noqa: E402
    generation_directory_name,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (  # noqa: E402
    TDG6_MINIMUM_MACRO_STEP,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_ledger import (  # noqa: E402
    IMP1_COMPLETE_STATE_CHANNELS,
    IMP1_REJECTION_EVENT_TYPE,
    append_imp1_rejection,
    imp1_checkpoint_extension,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_qualification import (  # noqa: E402
    SYNTHETIC_START as START,
    SYNTHETIC_TARGET as TARGET,
    make_synthetic_fixture,
    synthetic_diagnostics,
    synthetic_exponential_rhs,
)


REQUESTED_CAP = 1.0 / 32.0
CAMPAIGN_ID = "HLT17-STORE-SYNTHETIC"


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
        "assessment_sha256": "33" * 32,
        "preparation_sha256": "44" * 32,
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


def _identity(method: str) -> HLT17MemberIdentity:
    label = "RK4" if method == PRIMARY_METHOD else "SSPRK3"
    return HLT17MemberIdentity(
        scope=IDENTITY_SCOPE_SYNTHETIC,
        method_label=label,
        point_count=9,
        integrator_id=method,
        spatial_order=4 if label == "RK4" else 2,
        campaign_id=CAMPAIGN_ID,
        amplitude="synthetic",
    )


def _make_member(method: str = PRIMARY_METHOD, *, rhs=None) -> HLT17RuntimeMember:
    operator = synthetic_exponential_rhs if rhs is None else rhs
    state, transaction, tracers, ledger, coordinates = make_synthetic_fixture(method)
    return HLT17RuntimeMember(
        identity=_identity(method),
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


def _live_cohort(*, controllable: bool = False, inherited_source: int = 0, inherited_cfl: int = 0):
    checkpoints = {}
    for method in (PRIMARY_METHOD, COMPARATOR_METHOD):
        rhs = ControllableRHS() if controllable else synthetic_exponential_rhs
        member = _make_member(method, rhs=rhs)
        member.source_retry_count = inherited_source
        member.CFL_retry_count = inherited_cfl
        member.__post_init__()
        checkpoints[member.key] = HLT17InMemoryCheckpoint.seed(member, event_target=TARGET)
    bundles = {key: item.generation_bundle() for key, item in checkpoints.items()}
    return checkpoints, bundles


CHILD = r"""
import sys
from pathlib import Path
sys.path[:0] = sys.argv[1].split(":")
from recursive_horizons.fgc.evolution.hlt17_campaign_store import HLT17CampaignStore
root = Path(sys.argv[2])
action = sys.argv[3]
store = HLT17CampaignStore.open(root)
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


class HLT17CampaignStoreTests(unittest.TestCase):
    def test_seed_fine_retry_restore_both_methods_without_source_on_read(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            checkpoints, bundles = _live_cohort(controllable=True, inherited_source=2, inherited_cfl=1)
            store = HLT17CampaignStore.create(
                root, campaign_id=CAMPAIGN_ID, member_keys=HLT17_SYNTHETIC_MEMBER_KEYS
            )
            self.assertIn("physical_origin", store.remaining_seams)
            capability = store.acquire_writer(owner_token="owner", host="synthetic-host")
            seed = store.publish_seed(capability, bundles)
            self.assertEqual(seed.kind, "seed")
            for key in HLT17_SYNTHETIC_MEMBER_KEYS:
                state = seed.checkpoint_mapping["members"][key]
                self.assertEqual(state["source_total"], 2)
                self.assertEqual(state["cfl_total"], 1)
            latest = seed
            for key in ("RK4-9", "SSPRK3-9"):
                result = checkpoints[key].advance_once(requested_cap=REQUESTED_CAP)
                self.assertEqual(result.disposition, ACCEPTED_FINE)
                latest = store.publish_attempt(
                    capability,
                    kind="accepted_fine",
                    member_key=key,
                    successor_bundle=result.bundle,
                )
                self.assertEqual(latest.changed_member_key, key)
            rk4 = checkpoints["RK4-9"]
            rhs = rk4.member.operator
            assert isinstance(rhs, ControllableRHS)
            calls_before_retry = rhs.calls
            rhs.impulse_from = rk4.member.time
            rhs.mode = "impulse"
            temporal = rk4.advance_once(requested_cap=REQUESTED_CAP)
            self.assertEqual(temporal.disposition, "temporal_retry_required")
            latest = store.publish_attempt(
                capability,
                kind="temporal_retry",
                member_key="RK4-9",
                successor_bundle=temporal.bundle,
            )
            self.assertTrue(latest.journal_mapping["payload"]["accepted_state_preserved"])
            rk4.restore(rk4.accepted_encoding, encode_hlt17_cursor(rk4.cursor))
            rhs.impulse_from = rk4.member.time
            rhs.mode = "impulse"
            replay = replay_hlt17_rejection(rk4.cursor, rk4.live_boundary())
            rhs.mode = "source"
            source = rk4.absorb_nonfine(
                prepare_hlt17_next(rk4.cursor, rk4.live_boundary(), replay=replay)
            )
            self.assertEqual(source.disposition, "source_retry_required")
            latest = store.publish_attempt(
                capability,
                kind="source_retry",
                member_key="RK4-9",
                successor_bundle=source.bundle,
            )
            self.assertEqual(latest.checkpoint_mapping["members"]["RK4-9"]["source_total"], 3)
            store.close_writer(capability)
            calls_after_write = rhs.calls
            self.assertGreater(calls_after_write, calls_before_retry)
            view = HLT17CampaignStore.authenticate_read_only(root)
            self.assertEqual(view.journal_sha256, latest.journal_sha256)
            self.assertEqual(len(view.generations), 5)
            self.assertFalse(view.terminal)
            self.assertEqual(rhs.calls, calls_after_write)
            for method, key in ((PRIMARY_METHOD, "RK4-9"), (COMPARATOR_METHOD, "SSPRK3-9")):
                restored_member = _make_member(method, rhs=ControllableRHS())
                restored = HLT17InMemoryCheckpoint.seed(restored_member, event_target=TARGET)
                bundle = view.bundles[key]
                encoding = HLT17MemberEncoding(bundle.descriptor, bundle.payload)
                restored.restore(encoding, load_canonical_json(bundle.cursor_bytes))
                self.assertEqual(restored.member.time.hex(), encoding.decoded.metadata["accepted_time_hex"])
                self.assertEqual(restored.cursor.source_total, view.generations[-1].checkpoint_mapping["members"][key]["source_total"])
                self.assertEqual(restored.member.temporal_ledger.inherited_snapshot, rk4.member.temporal_ledger.inherited_snapshot if key == "RK4-9" else restored.member.temporal_ledger.inherited_snapshot)
            names = {path.name for path in root.iterdir()}
            self.assertNotIn("latest", names)
            self.assertNotIn("journal_tip.json", names)
            self.assertEqual(rhs.calls, calls_after_write)

    def test_thread_and_process_exclusion_and_unclosed_session(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            _checkpoints, bundles = _live_cohort()
            store = HLT17CampaignStore.create(
                root, campaign_id=CAMPAIGN_ID, member_keys=HLT17_SYNTHETIC_MEMBER_KEYS
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
            self.assertIsInstance(errors[0], HLT17WriterConflict)
            env = os.environ.copy()
            proc = subprocess.run(
                [
                    sys.executable,
                    "-c",
                    CHILD,
                    str(ROOT) + ":" + str(ROOT / "src"),
                    str(root),
                    "exit",
                ],
                check=False,
                capture_output=True,
                text=True,
                env=env,
            )
            self.assertIn("HLT17WriterConflict", proc.stdout)
            store.publish_seed(capability, bundles)
            store.close_writer(capability)
            child = subprocess.Popen(
                [
                    sys.executable,
                    "-c",
                    CHILD,
                    str(ROOT) + ":" + str(ROOT / "src"),
                    str(root),
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
            with self.assertRaises(HLT17WriterConflict):
                store.acquire_writer(owner_token="parent", host="synthetic-host")
            child.kill()
            child.wait()
            if child.stdin is not None:
                child.stdin.close()
            if child.stdout is not None:
                child.stdout.close()
            if child.stderr is not None:
                child.stderr.close()
            with self.assertRaises(HLT17UnclosedSession):
                store.acquire_writer(owner_token="parent", host="synthetic-host")

    def test_publication_faults_are_complete_or_absent(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            actual = publish_exclusive_file

            def refuse_lock(directory, relative, payload, **kwargs):
                if relative == WRITER_EXCLUSION_NAME:
                    raise PrepublicationError("injected lock refusal")
                return actual(directory, relative, payload, **kwargs)

            with patch(
                "recursive_horizons.fgc.evolution.hlt17_campaign_store.publish_exclusive_file",
                side_effect=refuse_lock,
            ), self.assertRaises(HLT17PostpublicationUncertainty) as partial:
                HLT17CampaignStore.create(
                    root,
                    campaign_id=CAMPAIGN_ID,
                    member_keys=HLT17_SYNTHETIC_MEMBER_KEYS,
                )
            self.assertTrue(partial.exception.published)
            self.assertTrue((root / ROOT_MANIFEST_NAME).is_file())
            self.assertFalse((root / WRITER_EXCLUSION_NAME).exists())
            with self.assertRaises(HLT17CampaignStoreError):
                HLT17CampaignStore.open(root)

        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            _checkpoints, bundles = _live_cohort()
            store = HLT17CampaignStore.create(
                root, campaign_id=CAMPAIGN_ID, member_keys=HLT17_SYNTHETIC_MEMBER_KEYS
            )
            capability = store.acquire_writer(owner_token="owner", host="synthetic-host")

            def before(cut: str) -> None:
                if cut == "before_publish":
                    raise OSError("injected before_publish")

            with self.assertRaises(HLT17CampaignStoreError):
                store.publish_seed(capability, bundles, _fault_hook=before)
            self.assertFalse((root / generation_directory_name(0)).exists())
            leftovers = list(root.glob(".generation-*.evidence-io-stage-*"))
            self.assertTrue(not leftovers or all(item.is_dir() or item.is_file() for item in leftovers))
            with self.assertRaises(HLT17CampaignStoreError):
                HLT17CampaignStore.authenticate_read_only(root)

            def after(cut: str) -> None:
                if cut == "after_publish":
                    raise OSError("injected after_publish")

            with self.assertRaises(HLT17PostpublicationUncertainty) as failed:
                store.publish_seed(capability, bundles, _fault_hook=after)
            self.assertTrue(failed.exception.published)
            self.assertTrue((root / generation_directory_name(0)).exists())
            view = HLT17CampaignStore.authenticate_read_only(root)
            self.assertEqual(view.generations[0].kind, "seed")
            with self.assertRaises(HLT17PostpublicationUncertainty):
                store.publish_attempt(
                    capability,
                    kind="resource_exhausted",
                    member_key="RK4-9",
                    successor_bundle=bundles["RK4-9"],
                    terminal_evidence={
                        "kind": "resource_exhausted",
                        "reason": "must not append after uncertain seed",
                        "executable_retry_cursor": False,
                    },
                )
            self.assertFalse((root / generation_directory_name(1)).exists())
            with self.assertRaises(HLT17PostpublicationUncertainty):
                store.close_writer(capability)
            store.abandon_writer(capability)
            with self.assertRaises(HLT17UnclosedSession):
                store.acquire_writer(owner_token="later", host="synthetic-host")
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            _checkpoints, bundles = _live_cohort()
            store = HLT17CampaignStore.create(
                root, campaign_id=CAMPAIGN_ID, member_keys=HLT17_SYNTHETIC_MEMBER_KEYS
            )
            capability = store.acquire_writer(owner_token="owner", host="synthetic-host")

            def after_fsync(cut: str) -> None:
                if cut == "after_parent_fsync":
                    raise OSError("injected after_parent_fsync")

            with self.assertRaises(HLT17PostpublicationUncertainty):
                store.publish_seed(capability, bundles, _fault_hook=after_fsync)
            self.assertTrue((root / generation_directory_name(0)).exists())
            view = HLT17CampaignStore.authenticate_read_only(root)
            self.assertEqual(view.generations[0].kind, "seed")
            with self.assertRaises(HLT17PostpublicationUncertainty):
                store.close_writer(capability)
            store.abandon_writer(capability)

    def test_symlink_hardlink_extra_partial_hole_and_rehash_controls(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            checkpoints, bundles = _live_cohort()
            store = HLT17CampaignStore.create(
                root, campaign_id=CAMPAIGN_ID, member_keys=HLT17_SYNTHETIC_MEMBER_KEYS
            )
            capability = store.acquire_writer(owner_token="owner", host="synthetic-host")
            seed = store.publish_seed(capability, bundles)
            result = checkpoints["RK4-9"].advance_once(requested_cap=REQUESTED_CAP)
            store.publish_attempt(
                capability,
                kind="accepted_fine",
                member_key="RK4-9",
                successor_bundle=result.bundle,
            )
            store.close_writer(capability)
            generation0 = root / generation_directory_name(0)
            (generation0 / "extra.json").write_bytes(b"{}")
            with self.assertRaises(HLT17CampaignStoreError):
                HLT17CampaignStore.authenticate_read_only(root)
            (generation0 / "extra.json").unlink()
            hole = root / generation_directory_name(2)
            hole.mkdir()
            (hole / "journal.json").write_bytes(seed.journal_raw)
            (hole / "checkpoint.json").write_bytes(seed.checkpoint_raw)
            (hole / "generation_manifest.json").write_bytes(seed.manifest_raw)
            with self.assertRaises(HLT17CampaignStoreError):
                HLT17CampaignStore.authenticate_read_only(root)
            hole.joinpath("journal.json").unlink()
            hole.joinpath("checkpoint.json").unlink()
            hole.joinpath("generation_manifest.json").unlink()
            hole.rmdir()
            target = generation0 / "journal.json"
            alias = root / "link-journal.json"
            os.symlink(target, alias)
            with self.assertRaises(HLT17CampaignStoreError):
                HLT17CampaignStore.authenticate_read_only(root)
            alias.unlink()
            linked = generation0 / "hard-journal.json"
            os.link(target, linked)
            with self.assertRaises(HLT17CampaignStoreError):
                HLT17CampaignStore.authenticate_read_only(root)
            linked.unlink()
            mutated = load_canonical_json((generation0 / "checkpoint.json").read_bytes())
            mutated["kind"] = "accepted_fine"
            (generation0 / "checkpoint.json").write_bytes(canonical_json_bytes(mutated))
            manifest = load_canonical_json((generation0 / "generation_manifest.json").read_bytes())
            manifest["checkpoint_sha256"] = sha256(
                (generation0 / "checkpoint.json").read_bytes()
            ).hexdigest()
            manifest["leaves"]["checkpoint.json"] = manifest["checkpoint_sha256"]
            (generation0 / "generation_manifest.json").write_bytes(canonical_json_bytes(manifest))
            with self.assertRaises(HLT17CampaignStoreError):
                HLT17CampaignStore.authenticate_read_only(root)

    def test_terminal_lock_closed_session_and_production_namespace(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            checkpoints, bundles = _live_cohort()
            store = HLT17CampaignStore.create(
                root, campaign_id=CAMPAIGN_ID, member_keys=HLT17_SYNTHETIC_MEMBER_KEYS
            )
            capability = store.acquire_writer(owner_token="owner", host="synthetic-host")
            seed = store.publish_seed(capability, bundles)
            terminal = store.publish_attempt(
                capability,
                kind="resource_exhausted",
                member_key="RK4-9",
                successor_bundle=seed.bundles["RK4-9"],
                terminal_evidence={
                    "kind": "resource_exhausted",
                    "reason": "synthetic typed terminal",
                    "executable_retry_cursor": False,
                },
            )
            self.assertTrue(terminal.terminal)
            self.assertTrue((root / TERMINAL_LOCK_NAME).is_file())
            with self.assertRaises(HLT17TerminalStore):
                store.publish_attempt(
                    capability,
                    kind="accepted_fine",
                    member_key="RK4-9",
                    successor_bundle=checkpoints["RK4-9"].generation_bundle(),
                )
            store.close_writer(capability)
            with self.assertRaises(HLT17TerminalStore):
                store.acquire_writer(owner_token="later", host="synthetic-host")
            view = HLT17CampaignStore.authenticate_read_only(root)
            self.assertTrue(view.terminal)
            self.assertFalse(view.unclosed_session)
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            nested = root / PRODUCTION_EVENT_PATH.replace("/", "_")
            nested.mkdir()
            store = HLT17CampaignStore.create(
                nested, campaign_id=CAMPAIGN_ID, member_keys=HLT17_SYNTHETIC_MEMBER_KEYS
            )
            self.assertTrue(store.root.exists())
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            forbidden = root / "runs" / "fgc-2-sf1" / "pro20-event1" / "calibration"
            forbidden.mkdir(parents=True)
            with self.assertRaises(HLT17CampaignStoreError):
                HLT17CampaignStore.create(
                    forbidden, campaign_id=CAMPAIGN_ID, member_keys=HLT17_SYNTHETIC_MEMBER_KEYS
                )
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            store = HLT17CampaignStore.create(
                root, campaign_id=CAMPAIGN_ID, member_keys=HLT17_SYNTHETIC_MEMBER_KEYS
            )
            first = store.acquire_writer(owner_token="one", host="synthetic-host")
            store.close_writer(first)
            second = store.acquire_writer(owner_token="two", host="synthetic-host")
            self.assertEqual(second.session_index, 1)
            store.close_writer(second)
            with self.assertRaises(HLT17CampaignStoreError):
                HLT17CampaignStore.create(
                    root, campaign_id=CAMPAIGN_ID, member_keys=HLT17_MEMBER_KEYS
                )

    def test_open_store_rereads_root_manifest_and_held_lock_inode(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            _checkpoints, bundles = _live_cohort()
            store = HLT17CampaignStore.create(
                root, campaign_id=CAMPAIGN_ID, member_keys=HLT17_SYNTHETIC_MEMBER_KEYS
            )
            capability = store.acquire_writer(owner_token="owner", host="synthetic-host")
            store.publish_seed(capability, bundles)
            original = read_regular_file

            def fake_root(root_path, relative, **kwargs):
                if relative == ROOT_MANIFEST_NAME:
                    return b"{}"
                return original(root_path, relative, **kwargs)

            with patch(
                "recursive_horizons.fgc.evolution.hlt17_campaign_store.read_regular_file",
                side_effect=fake_root,
            ):
                with self.assertRaises(HLT17CampaignStoreError):
                    store.authenticate()
                with self.assertRaises(HLT17CampaignStoreError):
                    HLT17CampaignStore.authenticate_read_only(root)
            lock_path = root / WRITER_EXCLUSION_NAME
            previous = lock_path.read_bytes()
            os.unlink(lock_path)
            lock_path.write_bytes(previous)
            with self.assertRaisesRegex(HLT17CampaignStoreError, "inode|lock"):
                store.authenticate()
            with self.assertRaisesRegex(HLT17CampaignStoreError, "inode|lock"):
                store.publish_attempt(
                    capability,
                    kind="resource_exhausted",
                    member_key="RK4-9",
                    successor_bundle=bundles["RK4-9"],
                    terminal_evidence={
                        "kind": "resource_exhausted",
                        "reason": "replaced lock",
                        "executable_retry_cursor": False,
                    },
                )
            store.abandon_writer(capability)

    def test_terminal_lock_prepublication_is_postpublication_uncertainty(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            _checkpoints, bundles = _live_cohort()
            store = HLT17CampaignStore.create(
                root, campaign_id=CAMPAIGN_ID, member_keys=HLT17_SYNTHETIC_MEMBER_KEYS
            )
            capability = store.acquire_writer(owner_token="owner", host="synthetic-host")
            seed = store.publish_seed(capability, bundles)
            from recursive_horizons.fgc.evolution import hlt17_campaign_store as store_mod

            real_publish = store_mod.publish_exclusive_file

            def fake_publish(target_root, relative, payload, **kwargs):
                if relative == TERMINAL_LOCK_NAME:
                    raise PrepublicationError("injected terminal lock refusal")
                return real_publish(target_root, relative, payload, **kwargs)

            with patch.object(store_mod, "publish_exclusive_file", side_effect=fake_publish):
                with self.assertRaises(HLT17PostpublicationUncertainty) as failed:
                    store.publish_attempt(
                        capability,
                        kind="resource_exhausted",
                        member_key="RK4-9",
                        successor_bundle=seed.bundles["RK4-9"],
                        terminal_evidence={
                            "kind": "resource_exhausted",
                            "reason": "terminal lock fault",
                            "executable_retry_cursor": False,
                        },
                    )
            self.assertTrue(failed.exception.published)
            self.assertTrue((root / generation_directory_name(1)).exists())
            self.assertFalse((root / TERMINAL_LOCK_NAME).exists())
            with self.assertRaises(HLT17PostpublicationUncertainty):
                store.publish_attempt(
                    capability,
                    kind="invalid_premise",
                    member_key="RK4-9",
                    successor_bundle=seed.bundles["RK4-9"],
                    terminal_evidence={
                        "kind": "invalid_premise",
                        "reason": "no terminal-lock bypass",
                        "spends_retry": False,
                        "executable_retry_cursor": False,
                    },
                )
            self.assertFalse((root / generation_directory_name(2)).exists())
            view = HLT17CampaignStore.authenticate_read_only(root)
            self.assertTrue(view.terminal)
            self.assertEqual(len(view.generations), 2)
            store.abandon_writer(capability)
            with self.assertRaises(HLT17UnclosedSession):
                store.acquire_writer(owner_token="later", host="synthetic-host")

    def test_forged_generation_flags_cannot_resume_a_terminal_predecessor(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            _checkpoints, bundles = _live_cohort()
            store = HLT17CampaignStore.create(
                root, campaign_id=CAMPAIGN_ID, member_keys=HLT17_SYNTHETIC_MEMBER_KEYS
            )
            capability = store.acquire_writer(owner_token="owner", host="synthetic-host")
            seed = store.publish_seed(capability, bundles)
            terminal = store.publish_attempt(
                capability,
                kind="resource_exhausted",
                member_key="RK4-9",
                successor_bundle=seed.bundles["RK4-9"],
                terminal_evidence={
                    "kind": "resource_exhausted",
                    "reason": "typed stop",
                    "executable_retry_cursor": False,
                },
            )
            with self.assertRaisesRegex(Exception, "disagrees with the immutable"):
                replace(terminal, terminal=False)
            store.close_writer(capability)
            view = HLT17CampaignStore.authenticate_read_only(root)
            with self.assertRaisesRegex(Exception, "disagrees with the immutable"):
                replace(view.generations[-1], terminal=False, first_event_complete=True)


if __name__ == "__main__":
    unittest.main()
