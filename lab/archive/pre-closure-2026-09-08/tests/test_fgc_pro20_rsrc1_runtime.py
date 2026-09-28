"""Prospective parent-only RSRC1 coordinator controls."""

from __future__ import annotations

import ast
from dataclasses import replace
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from recursive_horizons.fgc.evolution import pro20_rsrc1_attempt as attempt
from recursive_horizons.fgc.evolution import pro20_rsrc1_isolation as isolation
from recursive_horizons.fgc.evolution import pro20_rsrc1_member as member_owner
from recursive_horizons.fgc.evolution import pro20_rsrc1_protocol as protocol
from recursive_horizons.fgc.evolution import pro20_rsrc1_runtime as runtime
from recursive_horizons.fgc.evolution import pro20_rsrc1_seed as seed
from recursive_horizons.fgc.evolution import pro20_rsrc1_store as store
from recursive_horizons.fgc.evolution.hlt17_member_checkpoint import (
    HLT17InMemoryCheckpoint,
    origin_predecessor_identity,
)
from recursive_horizons.fgc.evolution.numerical_engine import (
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
)
from recursive_horizons.fgc.evolution.pro20_origin import MEMBER_KEYS
from tests.test_fgc_pro20_ev1_protocol import _make_live_member
from tests.test_fgc_pro20_rsrc1_seed import _bundle


ROOT = Path(__file__).resolve().parents[1]
IMPLEMENTATION_SHA256 = "2" * 64
AUTHORITY_COMMIT = "a" * 40


def _authority(**changes: object) -> dict[str, object]:
    mapping: dict[str, object] = {
        "authority_commit": AUTHORITY_COMMIT,
        "implementation_commit": "b" * 40,
        "campaign_id": member_owner.CAMPAIGN_ID,
        "output_namespace": protocol.PRODUCTION_NAMESPACE,
        "receipt_identities": {
            "authority_sha256": "1" * 64,
            "implementation_sha256": IMPLEMENTATION_SHA256,
            "config_sha256": "3" * 64,
            "source_sha256": "4" * 64,
            "origin_sha256": "5" * 64,
            "environment_sha256": "6" * 64,
        },
        "source_configuration": {
            key: f"{index:x}" * 64 for index, key in enumerate(MEMBER_KEYS, start=1)
        },
        "historical_source_tree": {
            "generation1_bridge_content_id": "7" * 64,
        },
    }
    mapping.update(changes)
    return mapping


def _checkpoint(member_key: str) -> HLT17InMemoryCheckpoint:
    method = PRIMARY_METHOD if member_key.startswith("RK4") else COMPARATOR_METHOD
    point_count = int(member_key.rsplit("-", 1)[1])
    live = _make_live_member(method, point_count)
    live.identity = replace(live.identity, campaign_id=member_owner.CAMPAIGN_ID)
    live.__post_init__()
    predecessor = origin_predecessor_identity(
        live,
        origin_cursor_sha256=sha256(member_key.encode("ascii")).hexdigest(),
    )
    return HLT17InMemoryCheckpoint.seed(
        live,
        event_target=member_owner.EVENT_TARGET,
        origin_predecessor=predecessor,
    )


class _Sequence:
    def __init__(self, values):
        self.values = list(values)
        self.last = self.values[-1]

    def __call__(self):
        if self.values:
            self.last = self.values.pop(0)
        return self.last


class RSRC1ParentRuntimeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.bundles = {key: _bundle(key) for key in MEMBER_KEYS}

    def _seed_child(self, _root: Path, request: seed.IsolatedSeedRequest):
        index = MEMBER_KEYS.index(request.member_key) + 1
        raw = seed.encode_seed_response(
            request,
            bundle=self.bundles[request.member_key],
            construction_sha256=f"{index:x}" * 64,
            peak_rss_bytes=1,
            wall_seconds=1.0,
        )
        return isolation.CompletedChild(0, raw, b"", 1.0)

    def _resource_attempt(
        self,
        _root: Path,
        request: isolation.IsolatedAttemptRequest,
    ) -> isolation.CompletedChild:
        classified = attempt.classify_attempt(
            member_key=request.member_key,
            predecessor_bundle=request.predecessor_bundle,
            successor_bundle=request.predecessor_bundle,
            accepted_state_advanced=False,
            resource_stop={"reason": "synthetic child resource"},
        )
        raw = isolation.encode_child_response(
            request,
            classified,
            peak_rss_bytes=1,
            wall_seconds=1.0,
        )
        return isolation.CompletedChild(0, raw, b"", 1.0)

    def _hooks(self, **changes: object) -> runtime.RuntimeHooks:
        values: dict[str, object] = {
            "authority_validator": lambda **_kwargs: _authority(),
            "seed_child": self._seed_child,
            "attempt_child": self._resource_attempt,
            "now": lambda: 0.0,
            "parent_rss": lambda: 1,
            "available_memory": lambda: isolation.MIN_AVAILABLE_MEMORY_BYTES,
            "disk_free": lambda _root: runtime.MIN_FREE_DISK_BYTES + 1,
            "namespace_bytes": lambda _root: 0,
            "implementation_identity_fn": lambda _root: {
                "implementation_sha256": IMPLEMENTATION_SHA256
            },
        }
        values.update(changes)
        return runtime.RuntimeHooks(**values)  # type: ignore[arg-type]

    @staticmethod
    def _root(temporary: str) -> Path:
        root = Path(temporary).resolve()
        (root / "runs" / "fgc-2-sf1").mkdir(parents=True)
        return root

    def test_status_has_no_authority_or_live_target(self) -> None:
        status = runtime.plan_status()
        self.assertEqual(status["campaign_id"], member_owner.CAMPAIGN_ID)
        self.assertEqual(status["namespace"], protocol.PRODUCTION_NAMESPACE)
        self.assertFalse(status["authority_implemented"])
        self.assertFalse(status["live_target_present"])
        self.assertFalse(status["campaign_execution_authorized"])

    def test_current_rss_and_available_memory_observers_are_strict(self) -> None:
        with patch.object(
            runtime.subprocess,
            "run",
            return_value=SimpleNamespace(returncode=0, stdout="123\n", stderr=""),
        ):
            self.assertEqual(runtime.current_parent_rss_bytes(), 123 * 1024)
        vm = """Mach Virtual Memory Statistics: (page size of 4096 bytes)\nPages free: 1.\nPages active: 99.\nPages inactive: 2.\nPages speculative: 3.\nPages purgeable: 4.\n"""
        with (
            patch.object(runtime.sys, "platform", "darwin"),
            patch.object(
                runtime.subprocess,
                "run",
                return_value=SimpleNamespace(returncode=0, stdout=vm, stderr=""),
            ),
        ):
            self.assertEqual(runtime.available_memory_bytes(), 10 * 4096)
        with patch.object(
            runtime.subprocess,
            "run",
            return_value=SimpleNamespace(returncode=0, stdout="ambiguous", stderr=""),
        ), self.assertRaises(runtime.PRO20RSRC1ResourceError):
            runtime.current_parent_rss_bytes()

    def test_absent_or_wrong_authority_refuses_before_seed_and_namespace(self) -> None:
        with TemporaryDirectory() as temporary:
            root = self._root(temporary)
            hooks = self._hooks(seed_child=lambda *_args: self.fail("seed child ran"))
            with self.assertRaises(runtime.PRO20RSRC1AuthorityError):
                runtime.run_first_event(root, authority_commit=None, hooks=hooks)
            self.assertFalse((root / store.PRODUCTION_CONTAINER).exists())
        wrong = _authority()
        wrong["receipt_identities"] = {
            **wrong["receipt_identities"],
            "implementation_sha256": "9" * 64,
        }
        with TemporaryDirectory() as temporary:
            root = self._root(temporary)
            hooks = self._hooks(
                authority_validator=lambda **_kwargs: wrong,
                seed_child=lambda *_args: self.fail("seed child ran"),
            )
            with self.assertRaises(runtime.PRO20RSRC1AuthorityError):
                runtime.run_first_event(
                    root,
                    authority_commit=AUTHORITY_COMMIT,
                    hooks=hooks,
                )
            self.assertFalse((root / store.PRODUCTION_CONTAINER).exists())

    def test_seed_resource_stop_leaves_no_namespace(self) -> None:
        def resource_seed(_root: Path, request: seed.IsolatedSeedRequest):
            raw = seed.encode_seed_memory_stop(request, peak_rss_bytes=1)
            return isolation.CompletedChild(0, raw, b"", 1.0)

        with TemporaryDirectory() as temporary:
            root = self._root(temporary)
            hooks = self._hooks(seed_child=resource_seed)
            with self.assertRaises(runtime.PRO20RSRC1SeedResourceError):
                runtime.run_first_event(
                    root,
                    authority_commit=AUTHORITY_COMMIT,
                    hooks=hooks,
                )
            self.assertFalse((root / store.PRODUCTION_CONTAINER).exists())
            self.assertEqual(hooks.seed_requests, [MEMBER_KEYS[0]])

    def test_clean_child_resource_publishes_terminal_and_closes_writer(self) -> None:
        with TemporaryDirectory() as temporary:
            root = self._root(temporary)
            hooks = self._hooks()
            result = runtime.run_first_event(
                root,
                authority_commit=AUTHORITY_COMMIT,
                hooks=hooks,
            )
            self.assertTrue(result["terminal"])
            self.assertFalse(result["first_event_complete"])
            self.assertEqual(result["disposition"], "resource_exhausted")
            self.assertEqual(hooks.seed_requests, list(MEMBER_KEYS))
            self.assertEqual(hooks.attempt_requests, [MEMBER_KEYS[0]])
            view = store.PRO20RSRC1CampaignStore.authenticate_read_only(root)
            self.assertFalse(view.unclosed_session)

    def test_accepted_child_bytes_publish_before_a_later_resource_terminal(self) -> None:
        calls = 0

        def accept_then_stop(
            _root: Path,
            request: isolation.IsolatedAttemptRequest,
        ) -> isolation.CompletedChild:
            nonlocal calls
            calls += 1
            if calls == 1:
                checkpoint = _checkpoint(request.member_key)
                attempt.restore_checkpoint_from_bundle(
                    checkpoint,
                    request.predecessor_bundle,
                )
                classified = attempt.execute_scheduled_attempt(checkpoint)
            else:
                classified = attempt.classify_attempt(
                    member_key=request.member_key,
                    predecessor_bundle=request.predecessor_bundle,
                    successor_bundle=request.predecessor_bundle,
                    accepted_state_advanced=False,
                    resource_stop={"reason": "after one accepted fine"},
                )
            raw = isolation.encode_child_response(
                request,
                classified,
                peak_rss_bytes=1,
                wall_seconds=1.0,
            )
            return isolation.CompletedChild(0, raw, b"", 1.0)

        with TemporaryDirectory() as temporary:
            root = self._root(temporary)
            result = runtime.run_first_event(
                root,
                authority_commit=AUTHORITY_COMMIT,
                hooks=self._hooks(attempt_child=accept_then_stop),
            )
            self.assertEqual(result["published_attempts"], 2)
            view = store.PRO20RSRC1CampaignStore.authenticate_read_only(root)
            self.assertEqual(view.generations[-2].kind, "accepted_fine")
            self.assertEqual(view.generations[-1].kind, "resource_exhausted")
            self.assertGreater(
                view.generations[-1].bundles[MEMBER_KEYS[0]].accepted_generation,
                0,
            )

    def test_abnormal_attempt_abandons_writer_and_is_not_a_resource_terminal(self) -> None:
        def fail_child(_root: Path, _request: isolation.IsolatedAttemptRequest):
            return isolation.CompletedChild(9, b"", b"", 1.0)

        with TemporaryDirectory() as temporary:
            root = self._root(temporary)
            hooks = self._hooks(attempt_child=fail_child)
            with self.assertRaises(runtime.PRO20RSRC1ForensicError):
                runtime.run_first_event(
                    root,
                    authority_commit=AUTHORITY_COMMIT,
                    hooks=hooks,
                )
            view = store.PRO20RSRC1CampaignStore.authenticate_read_only(root)
            self.assertTrue(view.unclosed_session)
            self.assertFalse(view.terminal)

    def test_parent_resource_stop_after_seed_is_typed_and_closed(self) -> None:
        rss = _Sequence(
            [1] * 7 + [isolation.PARENT_MAX_CURRENT_RSS_BYTES + 1]
        )
        with TemporaryDirectory() as temporary:
            root = self._root(temporary)
            hooks = self._hooks(
                parent_rss=rss,
                attempt_child=lambda *_args: self.fail("attempt child ran"),
            )
            result = runtime.run_first_event(
                root,
                authority_commit=AUTHORITY_COMMIT,
                hooks=hooks,
            )
            self.assertTrue(result["terminal"])
            self.assertEqual(result["disposition"], "resource_exhausted")
            self.assertEqual(hooks.attempt_requests, [])
            self.assertFalse(
                store.PRO20RSRC1CampaignStore.authenticate_read_only(root).unclosed_session
            )

    def test_closed_old_container_is_untouched_by_new_terminal(self) -> None:
        with TemporaryDirectory() as temporary:
            root = self._root(temporary)
            old = root / "runs" / "fgc-2-sf1" / "pro20-event1"
            old.mkdir()
            sentinel = old / "immutable"
            sentinel.write_bytes(b"old-store\n")
            result = runtime.run_first_event(
                root,
                authority_commit=AUTHORITY_COMMIT,
                hooks=self._hooks(),
            )
            self.assertTrue(result["terminal"])
            self.assertEqual(sentinel.read_bytes(), b"old-store\n")

    def test_runtime_imports_only_the_new_store_and_protocol(self) -> None:
        path = ROOT / "src/recursive_horizons/fgc/evolution/pro20_rsrc1_runtime.py"
        tree = ast.parse(path.read_text(encoding="utf-8"))
        imports: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                imports.append(node.module or "")
            elif isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
        joined = " ".join(imports)
        self.assertNotIn("pro20_ev1_store", joined)
        self.assertNotIn("pro20_ev1_protocol", joined)
        self.assertNotIn("pro20_ev1_runtime", joined)


if __name__ == "__main__":
    unittest.main()
