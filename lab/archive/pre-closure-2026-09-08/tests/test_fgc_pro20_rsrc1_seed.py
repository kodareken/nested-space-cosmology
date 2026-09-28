"""Store-blind RSRC1 per-member seed and complete-cohort controls."""

from __future__ import annotations

import ast
from copy import deepcopy
from dataclasses import replace
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from recursive_horizons.evidence_io import canonical_json_bytes
from recursive_horizons.fgc.evolution.hlt17_member_checkpoint import (
    HLT17InMemoryCheckpoint,
    origin_predecessor_identity,
)
from recursive_horizons.fgc.evolution.numerical_engine import (
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
)
from recursive_horizons.fgc.evolution.pro20_origin import MEMBER_KEYS
from recursive_horizons.fgc.evolution import pro20_rsrc1_isolation as isolation
from recursive_horizons.fgc.evolution import pro20_rsrc1_member as member_owner
from recursive_horizons.fgc.evolution import pro20_rsrc1_seed as seed
from tests.test_fgc_pro20_ev1_protocol import _make_live_member


ROOT = Path(__file__).resolve().parents[1]


def _request(member_key: str) -> seed.IsolatedSeedRequest:
    index = MEMBER_KEYS.index(member_key) + 1
    return seed.IsolatedSeedRequest(
        member_key=member_key,
        authority_receipt_sha256="1" * 64,
        origin_capture_sha256="2" * 64,
        generation1_bridge_content_id="3" * 64,
        source_closure_sha256="4" * 64,
        source_configuration_sha256=f"{index:x}" * 64,
    )


def _bundle(member_key: str):
    method = PRIMARY_METHOD if member_key.startswith("RK4") else COMPARATOR_METHOD
    point_count = int(member_key.rsplit("-", 1)[1])
    runtime = _make_live_member(method, point_count)
    runtime.identity = replace(
        runtime.identity,
        campaign_id=member_owner.CAMPAIGN_ID,
    )
    runtime.__post_init__()
    predecessor = origin_predecessor_identity(
        runtime,
        origin_cursor_sha256=sha256(member_key.encode("ascii")).hexdigest(),
    )
    checkpoint = HLT17InMemoryCheckpoint.seed(
        runtime,
        event_target=member_owner.EVENT_TARGET,
        origin_predecessor=predecessor,
    )
    return checkpoint.generation_bundle()


class RSRC1SeedTests(unittest.TestCase):
    def test_seed_request_and_success_response_roundtrip(self) -> None:
        request = _request("RK4-2049")
        self.assertEqual(
            seed.decode_seed_request(seed.encode_seed_request(request)),
            request,
        )
        response = seed.encode_seed_response(
            request,
            bundle=_bundle(request.member_key),
            construction_sha256="a" * 64,
            peak_rss_bytes=1,
            wall_seconds=1.0,
        )
        outcome = seed.decode_seed_response(request, response)
        self.assertTrue(outcome.successful)
        self.assertEqual(outcome.member_key, request.member_key)
        self.assertEqual(outcome.construction_sha256, "a" * 64)

    def test_peak_overrun_and_memory_stop_return_no_bundle(self) -> None:
        request = _request("RK4-2049")
        overrun = seed.encode_seed_response(
            request,
            bundle=_bundle(request.member_key),
            construction_sha256="a" * 64,
            peak_rss_bytes=isolation.CHILD_MAX_PEAK_RSS_BYTES + 1,
            wall_seconds=1.0,
        )
        for raw in (
            overrun,
            seed.encode_seed_memory_stop(request, peak_rss_bytes=1),
        ):
            outcome = seed.decode_seed_response(request, raw)
            with self.subTest(raw=raw[:32]):
                self.assertFalse(outcome.successful)
                self.assertIsNone(outcome.bundle)
                self.assertIsNone(outcome.construction_sha256)
                self.assertEqual(
                    outcome.resource_evidence["limit_bytes"],
                    isolation.CHILD_MAX_PEAK_RSS_BYTES,
                )

    def test_response_hash_member_endpoint_and_resource_attacks_fail(self) -> None:
        request = _request("RK4-2049")
        response = json.loads(
            seed.encode_seed_response(
                request,
                bundle=_bundle(request.member_key),
                construction_sha256="a" * 64,
                peak_rss_bytes=1,
                wall_seconds=1.0,
            )
        )
        mutations = (
            ("request_sha256", "0" * 64),
            ("member_key", "RK4-4097"),
            ("outcome", "resource_stop"),
            ("child_peak_rss_bytes", isolation.CHILD_MAX_PEAK_RSS_BYTES + 1),
        )
        for key, value in mutations:
            changed = deepcopy(response)
            changed[key] = value
            with self.subTest(key=key), self.assertRaises(seed.PRO20RSRC1SeedError):
                seed.decode_seed_response(request, canonical_json_bytes(changed))

    def test_abnormal_completed_seed_is_forensic_not_a_resource_stop(self) -> None:
        request = _request("RK4-2049")
        for completed in (
            isolation.CompletedChild(9, b"", b"", 1.0),
            isolation.CompletedChild(0, b"{}", b"warning", 1.0),
            isolation.CompletedChild(0, b"{}", b"", 1.0),
        ):
            with self.subTest(completed=completed), self.assertRaises(
                isolation.PRO20RSRC1ForensicUncertainty
            ) as caught:
                seed.reduce_completed_seed_child(request, completed)
            self.assertFalse(caught.exception.publish_typed_terminal)

    def test_child_constructs_only_one_seed_and_writes_nothing(self) -> None:
        request = _request("RK4-2049")
        calls: list[str] = []
        bundle = _bundle(request.member_key)

        def builder(_origin, *, member_key, static_input_bytes, source_closure_sha256):
            calls.append(member_key)
            self.assertEqual(static_input_bytes, {})
            self.assertEqual(source_closure_sha256, request.source_closure_sha256)
            return SimpleNamespace(
                captured_origin_sha256=request.origin_capture_sha256,
                generation1_bridge_content_id=request.generation1_bridge_content_id,
                source_closure_sha256=request.source_closure_sha256,
                source_configuration_sha256=request.source_configuration_sha256,
                bundle=bundle,
                construction_sha256="a" * 64,
            )

        with TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            raw = seed.execute_seed_child_request(
                root,
                seed.encode_seed_request(request),
                origin_loader=lambda _root: object(),
                static_loader=lambda _root: {},
                member_builder=builder,
                peak_rss=lambda: 1,
                now=iter((0.0, 1.0)).__next__,
            )
            self.assertTrue(seed.decode_seed_response(request, raw).successful)
            self.assertEqual(calls, [request.member_key])
            self.assertEqual(list(root.iterdir()), [])

    def test_complete_six_member_cohort_is_required_before_future_publish(self) -> None:
        outcomes = {}
        for index, key in enumerate(MEMBER_KEYS, start=1):
            request = _request(key)
            raw = seed.encode_seed_response(
                request,
                bundle=_bundle(key),
                construction_sha256=f"{index:x}" * 64,
                peak_rss_bytes=1,
                wall_seconds=1.0,
            )
            outcomes[key] = seed.decode_seed_response(request, raw)
        cohort = seed.bind_seed_cohort(outcomes)
        self.assertEqual(tuple(cohort.bundles), MEMBER_KEYS)
        self.assertEqual(len(cohort.cohort_sha256), 64)
        incomplete = dict(outcomes)
        incomplete.pop(MEMBER_KEYS[-1])
        with self.assertRaises(seed.PRO20RSRC1SeedError):
            seed.bind_seed_cohort(incomplete)
        resource = dict(outcomes)
        request = _request(MEMBER_KEYS[-1])
        resource[MEMBER_KEYS[-1]] = seed.decode_seed_response(
            request,
            seed.encode_seed_memory_stop(request, peak_rss_bytes=1),
        )
        with self.assertRaises(seed.PRO20RSRC1SeedError):
            seed.bind_seed_cohort(resource)

    def test_seed_module_has_no_production_store_runner_or_write_path(self) -> None:
        path = ROOT / "src/recursive_horizons/fgc/evolution/pro20_rsrc1_seed.py"
        tree = ast.parse(path.read_text(encoding="utf-8"))
        imports: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                imports.append(node.module or "")
            elif isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
        joined = " ".join(imports)
        self.assertNotIn("pro20_ev1_store", joined)
        self.assertNotIn("pro20_ev1_runtime", joined)
        self.assertNotIn("campaign_store", joined)
        source = path.read_text(encoding="utf-8")
        self.assertNotIn("publish_", source)
        self.assertNotIn("write_bytes", source)


if __name__ == "__main__":
    unittest.main()
