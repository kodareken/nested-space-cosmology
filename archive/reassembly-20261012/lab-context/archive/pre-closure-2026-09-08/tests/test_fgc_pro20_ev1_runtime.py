"""Bounded PRO20-EV1 runner controls. No physical source or live namespace."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
import inspect
import tempfile
import unittest
from unittest.mock import patch

from recursive_horizons.fgc.evolution.hlt17_imp1_bridge import HLT17BridgeResult
from recursive_horizons.fgc.evolution.hlt17_member_checkpoint import (
    HLT17InMemoryCheckpoint,
)
from recursive_horizons.fgc.evolution.numerical_engine import (
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
)
from recursive_horizons.fgc.evolution.pro20_ev1_authority import (
    GENERATION1_BRIDGE_CONTENT_ID,
    ORIGIN_CAPTURE_SHA256,
    SOURCE_CONFIGURATION_SHA256_BY_MEMBER,
)
from recursive_horizons.fgc.evolution.pro20_ev1_protocol import (
    PRODUCTION_NAMESPACE,
)
from recursive_horizons.fgc.evolution.pro20_ev1_runtime import (
    MAX_NAMESPACE_BYTES,
    MAX_PUBLISHED_ATTEMPTS,
    MAX_RSS_BYTES,
    MAX_TOTAL_WALL_SECONDS,
    MIN_FREE_DISK_BYTES,
    PLANNING_ACCEPTED_STEP_ESTIMATE,
    PRO20EV1AuthoritySeamError,
    PRO20EV1PremiseError,
    PRO20EV1ResourceError,
    PRO20EV1UnexpectedError,
    RuntimeHooks,
    check_resource_ceilings,
    classify_attempt,
    execute_scheduled_attempt,
    fresh_requested_cap,
    pending_plan_from_cursor,
    require_authority_delta,
    run_first_event,
    schedule_next_member,
)
from recursive_horizons.fgc.evolution.pro20_ev1_store import (
    PRO20EV1CampaignStore,
)
from recursive_horizons.fgc.evolution.pro20_origin import MEMBER_KEYS
from recursive_horizons.fgc.evolution.pro20_source_factory import CAMPAIGN_ID

from tests.test_fgc_pro20_ev1_protocol import (
    CAMPAIGN_ID as FIXTURE_CAMPAIGN_ID,
    ControllableRHS,
    _make_live_member,
    _seed_live,
)


ARTIFACT_ID = "FGC-1-PRO20-EV1-FRZ1"
AUTHORITY_COMMIT = "a" * 40
IMPLEMENTATION_COMMIT = "b" * 40
IMPLEMENTATION_SHA256 = "2" * 64
SOURCE_SHA256 = "4" * 64


def _checkpoint(member_key: str = "RK4-2049") -> HLT17InMemoryCheckpoint:
    method = PRIMARY_METHOD if member_key.startswith("RK4") else COMPARATOR_METHOD
    point_count = int(member_key.rsplit("-", 1)[1])
    rhs = ControllableRHS()
    member = _make_live_member(
        method,
        point_count,
        rhs=rhs,
        inherited_source=3,
        inherited_cfl=1,
    )
    from hashlib import sha256

    from recursive_horizons.fgc.evolution.hlt17_member_checkpoint import (
        origin_predecessor_identity,
    )

    predecessor = origin_predecessor_identity(
        member,
        origin_cursor_sha256=sha256(member.key.encode("ascii")).hexdigest(),
    )
    return HLT17InMemoryCheckpoint.seed(
        member,
        event_target=1.5,
        origin_predecessor=predecessor,
    )


@dataclass
class _Clock:
    values: list[float]

    def __post_init__(self) -> None:
        self.last = self.values[-1]

    def __call__(self) -> float:
        if self.values:
            self.last = self.values.pop(0)
        return self.last


def _authority() -> dict[str, object]:
    return {
        "authority_commit": AUTHORITY_COMMIT,
        "implementation_commit": IMPLEMENTATION_COMMIT,
        "receipt_identities": {
            "authority_sha256": "1" * 64,
            "implementation_sha256": IMPLEMENTATION_SHA256,
            "config_sha256": "3" * 64,
            "source_sha256": SOURCE_SHA256,
            "origin_sha256": ORIGIN_CAPTURE_SHA256,
            "environment_sha256": "5" * 64,
        },
        "source_configuration": dict(SOURCE_CONFIGURATION_SHA256_BY_MEMBER),
        "historical_source_tree": {
            "generation1_bridge_content_id": GENERATION1_BRIDGE_CONTENT_ID,
        },
        "campaign_id": CAMPAIGN_ID,
    }


class PRO20EV1RuntimeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        _root, cls.seed, checkpoints = _seed_live(
            inherited_source=3,
            inherited_cfl=1,
            controllable=True,
        )
        members = {}
        for key in MEMBER_KEYS:
            members[key] = SimpleNamespace(
                checkpoint=checkpoints[key],
                bundle=cls.seed.bundles[key],
                source_configuration_sha256=SOURCE_CONFIGURATION_SHA256_BY_MEMBER[
                    key
                ],
                captured_origin_sha256=ORIGIN_CAPTURE_SHA256,
                source_closure_sha256=SOURCE_SHA256,
            )
        cls.origin = SimpleNamespace(
            capture_sha256=ORIGIN_CAPTURE_SHA256,
            generation1_bridge_content_id=GENERATION1_BRIDGE_CONTENT_ID,
        )
        cls.runtime_origin = SimpleNamespace(
            captured_origin_sha256=ORIGIN_CAPTURE_SHA256,
            source_closure_sha256=SOURCE_SHA256,
            members=members,
            construction_sha256="6" * 64,
        )

    def _hooks(self, **changes: object) -> RuntimeHooks:
        values: dict[str, object] = {
            "origin_loader": lambda _root: self.origin,
            "runtime_origin_builder": lambda *_args, **_kwargs: self.runtime_origin,
            "static_input_loader": lambda _root: {},
            "authority_validator": lambda **_kwargs: _authority(),
            "implementation_identity_fn": lambda _root: {
                "sha256": IMPLEMENTATION_SHA256
            },
            "rss": lambda: 1,
            "disk_free": lambda _root: MIN_FREE_DISK_BYTES + 1,
            "namespace_bytes": lambda _root: 0,
        }
        values.update(changes)
        return RuntimeHooks(**values)  # type: ignore[arg-type]

    def test_fresh_cap_schedule_and_pending_contract(self) -> None:
        checkpoint = _checkpoint()
        self.assertEqual(fresh_requested_cap(checkpoint).hex(), (1.0 / 16.0).hex())
        self.assertEqual(schedule_next_member(self.seed), "RK4-2049")
        self.assertIsNone(pending_plan_from_cursor(checkpoint.cursor))
        with self.assertRaises(TypeError):
            fresh_requested_cap(SimpleNamespace())  # type: ignore[arg-type]

    def test_execute_accepts_one_fine_endpoint(self) -> None:
        result = execute_scheduled_attempt(_checkpoint())
        self.assertEqual(result.kind, "accepted_fine")
        self.assertTrue(result.accepted_state_advanced)
        self.assertFalse(result.physical_state_preserved)

    def test_source_cfl_and_temporal_retries_are_owned(self) -> None:
        for mode, expected in (
            ("source", "source_retry"),
            ("cfl", "cfl_retry"),
            ("impulse", "temporal_retry"),
        ):
            checkpoint = _checkpoint()
            rhs = checkpoint.member.operator
            self.assertIsInstance(rhs, ControllableRHS)
            rhs.impulse_from = checkpoint.member.time
            rhs.mode = mode
            result = execute_scheduled_attempt(checkpoint)
            with self.subTest(mode=mode):
                self.assertEqual(result.kind, expected)
                self.assertTrue(result.physical_state_preserved)
                self.assertTrue(result.executable_retry_cursor)

    def test_prepare_time_temporal_exhaustion_keeps_bundle_unchanged(self) -> None:
        checkpoint = _checkpoint()
        exhausted = HLT17BridgeResult(
            disposition="temporal_retry_exhausted",
            cursor=checkpoint.cursor,
            evidence={
                "kind": "temporal_retry_exhausted",
                "reason": "maximum_temporal_retries",
                "updated_imp1_ledger": {"schema": "synthetic"},
                "executable_retry_cursor": False,
            },
        )
        with patch.object(
            HLT17InMemoryCheckpoint,
            "prepare",
            return_value=exhausted,
        ):
            result = execute_scheduled_attempt(checkpoint)
        self.assertEqual(result.kind, "temporal_retry_exhausted")
        self.assertFalse(result.executable_retry_cursor)
        self.assertEqual(
            result.successor_bundle.cursor_bytes,
            checkpoint.generation_bundle().cursor_bytes,
        )

    def test_resource_evidence_cannot_promote_terminal(self) -> None:
        bundle = self.seed.bundles["RK4-2049"]
        result = classify_attempt(
            member_key="RK4-2049",
            predecessor_bundle=bundle,
            successor_bundle=bundle,
            accepted_state_advanced=False,
            resource_stop={
                "reason": "rss",
                "kind": "accepted_fine",
                "executable_retry_cursor": True,
            },
        )
        self.assertEqual(result.kind, "resource_exhausted")
        self.assertFalse(result.terminal_evidence["executable_retry_cursor"])
        self.assertFalse(result.terminal_evidence["spends_retry"])

    def test_each_resource_ceiling_is_typed(self) -> None:
        base = {
            "repository_root": Path.cwd(),
            "started": 0.0,
            "now": lambda: 0.0,
            "rss": lambda: 1,
            "disk_free": lambda _root: MIN_FREE_DISK_BYTES + 1,
            "namespace_bytes": lambda _root: 0,
            "published_attempts": 0,
        }
        mutations = (
            {"now": lambda: MAX_TOTAL_WALL_SECONDS + 1},
            {"rss": lambda: MAX_RSS_BYTES + 1},
            {"disk_free": lambda _root: MIN_FREE_DISK_BYTES - 1},
            {"namespace_bytes": lambda _root: MAX_NAMESPACE_BYTES + 1},
            {"published_attempts": MAX_PUBLISHED_ATTEMPTS + 1},
        )
        for mutation in mutations:
            with self.subTest(mutation=tuple(mutation)), self.assertRaises(
                PRO20EV1ResourceError
            ):
                check_resource_ceilings(**{**base, **mutation})

    def test_post_seed_resource_stop_publishes_and_closes(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            (root / "runs" / "fgc-2-sf1").mkdir(parents=True)
            hooks = self._hooks(now=_Clock([0.0, 0.0, MAX_TOTAL_WALL_SECONDS + 1]))
            result = run_first_event(
                root,
                authority_commit=AUTHORITY_COMMIT,
                hooks=hooks,
            )
            self.assertEqual(result["disposition"], "resource_exhausted")
            self.assertTrue(result["terminal"])
            self.assertFalse(result["unclosed_session"])
            self.assertFalse(result["calibration_eligible"])
            view = PRO20EV1CampaignStore.authenticate_read_only(root)
            self.assertTrue(view.terminal)

    def test_unexpected_step_leaves_unclosed_forensic_session(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            (root / "runs" / "fgc-2-sf1").mkdir(parents=True)

            def fail(_checkpoint: HLT17InMemoryCheckpoint):
                raise RuntimeError("injected programmer failure")

            hooks = self._hooks(now=lambda: 0.0, step_fn=fail)
            with self.assertRaises(PRO20EV1UnexpectedError):
                run_first_event(
                    root,
                    authority_commit=AUTHORITY_COMMIT,
                    hooks=hooks,
                )
            view = PRO20EV1CampaignStore.authenticate_read_only(root)
            self.assertTrue(view.unclosed_session)

    def test_authority_errors_are_typed_and_precede_source(self) -> None:
        def reject(**_kwargs: object):
            raise ValueError("not the freeze")

        with self.assertRaises(PRO20EV1AuthoritySeamError):
            require_authority_delta(
                repository_root=Path.cwd(),
                authority_commit=AUTHORITY_COMMIT,
                implementation_sha256=IMPLEMENTATION_SHA256,
                validator=reject,
            )
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            (root / "runs" / "fgc-2-sf1").mkdir(parents=True)
            hooks = self._hooks(authority_validator=reject)
            with self.assertRaises(PRO20EV1AuthoritySeamError):
                run_first_event(
                    root,
                    authority_commit=AUTHORITY_COMMIT,
                    hooks=hooks,
                )
            self.assertEqual(hooks.source_constructions, [])
            self.assertFalse((root / PRODUCTION_NAMESPACE).exists())

    def test_runtime_origin_identity_drift_refuses_before_seed(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            (root / "runs" / "fgc-2-sf1").mkdir(parents=True)
            drifted = SimpleNamespace(
                captured_origin_sha256="7" * 64,
                source_closure_sha256=SOURCE_SHA256,
                members=self.runtime_origin.members,
                construction_sha256="6" * 64,
            )
            hooks = self._hooks(
                runtime_origin_builder=lambda *_args, **_kwargs: drifted
            )
            with self.assertRaises(PRO20EV1PremiseError) as caught:
                run_first_event(
                    root,
                    authority_commit=AUTHORITY_COMMIT,
                    hooks=hooks,
                )
            self.assertTrue(caught.exception.physical_source_constructed)
            self.assertFalse((root / "runs" / "fgc-2-sf1" / "pro20-event1").exists())

    def test_fixture_campaign_is_not_live_authority(self) -> None:
        self.assertNotEqual(FIXTURE_CAMPAIGN_ID, CAMPAIGN_ID)
        self.assertEqual(PLANNING_ACCEPTED_STEP_ESTIMATE, 224)
        self.assertNotIn("receipt", inspect.signature(run_first_event).parameters)


if __name__ == "__main__":
    unittest.main()
