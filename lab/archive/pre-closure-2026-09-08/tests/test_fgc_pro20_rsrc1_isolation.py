"""RSRC1 request/response, resource, and forensic-boundary controls."""

from __future__ import annotations

import ast
from copy import deepcopy
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from recursive_horizons.evidence_io import canonical_json_bytes
from recursive_horizons.fgc.evolution import pro20_rsrc1_attempt as attempt
from recursive_horizons.fgc.evolution import pro20_rsrc1_isolation as isolation
from tests.test_fgc_pro20_ev1_runtime import _checkpoint


ROOT = Path(__file__).resolve().parents[1]


def _request():
    return isolation.IsolatedAttemptRequest(
        member_key="RK4-2049",
        predecessor_bundle=_checkpoint().generation_bundle(),
        authority_receipt_sha256="1" * 64,
        origin_capture_sha256="2" * 64,
        generation1_bridge_content_id="3" * 64,
        source_closure_sha256="4" * 64,
        source_configuration_sha256="5" * 64,
    )


class RSRC1IsolationWireTests(unittest.TestCase):
    def test_request_and_accepted_response_roundtrip(self) -> None:
        request = _request()
        self.assertEqual(
            isolation.decode_attempt_request(isolation.encode_attempt_request(request)),
            request,
        )
        result = attempt.execute_scheduled_attempt(_checkpoint())
        raw = isolation.encode_child_response(
            request,
            result,
            peak_rss_bytes=1024,
            wall_seconds=1.0,
        )
        mapping = json.loads(raw)
        self.assertNotIn("kind", mapping)
        self.assertNotIn("structural_class", mapping)
        rebound, metrics = isolation.classify_child_response(request, raw)
        self.assertEqual(rebound.kind, "accepted_fine")
        self.assertTrue(rebound.accepted_state_advanced)
        self.assertEqual(metrics.peak_rss_bytes, 1024)

    def test_peak_overrun_discards_child_endpoint_as_typed_resource(self) -> None:
        request = _request()
        accepted = attempt.execute_scheduled_attempt(_checkpoint())
        raw = isolation.encode_child_response(
            request,
            accepted,
            peak_rss_bytes=isolation.CHILD_MAX_PEAK_RSS_BYTES + 1,
            wall_seconds=1.0,
        )
        rebound, metrics = isolation.classify_child_response(request, raw)
        self.assertEqual(rebound.kind, "resource_exhausted")
        self.assertFalse(rebound.accepted_state_advanced)
        self.assertEqual(
            rebound.successor_bundle.payload,
            request.predecessor_bundle.payload,
        )
        self.assertEqual(
            rebound.terminal_evidence["observed_peak_rss_bytes"],
            metrics.peak_rss_bytes,
        )

    def test_hash_canonicality_and_promotion_attacks_fail(self) -> None:
        request = _request()
        raw = isolation.encode_attempt_request(request)
        with self.assertRaises(isolation.PRO20RSRC1IsolationError):
            isolation.decode_attempt_request(raw + b"\n")
        changed = json.loads(raw)
        changed["request_sha256"] = "0" * 64
        with self.assertRaises(isolation.PRO20RSRC1IsolationError):
            isolation.decode_attempt_request(canonical_json_bytes(changed))
        result = attempt.execute_scheduled_attempt(_checkpoint())
        response = json.loads(
            isolation.encode_child_response(
                request,
                result,
                peak_rss_bytes=1,
                wall_seconds=1.0,
            )
        )
        for key, value in (
            ("member_key", "RK4-4097"),
            ("request_sha256", "0" * 64),
            ("stop_kind", "resource"),
            (
                "child_peak_rss_bytes",
                isolation.CHILD_MAX_PEAK_RSS_BYTES + 1,
            ),
        ):
            altered = deepcopy(response)
            altered[key] = value
            with self.subTest(key=key), self.assertRaises(
                isolation.PRO20RSRC1IsolationError
            ):
                isolation.classify_child_response(
                    request, canonical_json_bytes(altered)
                )
        altered = deepcopy(response)
        altered["stop_kind"] = "resource"
        altered["stop_evidence"] = {"reason": "rss"}
        altered["accepted_state_advanced"] = True
        with self.assertRaises(isolation.PRO20RSRC1IsolationError):
            isolation.classify_child_response(request, canonical_json_bytes(altered))

    def test_nonzero_stderr_and_malformed_children_are_forensic_not_resource(self) -> None:
        request = _request()
        cases = (
            isolation.CompletedChild(9, b"", b"", 1.0),
            isolation.CompletedChild(0, b"{}", b"warning", 1.0),
            isolation.CompletedChild(0, b"{}", b"", 1.0),
        )
        for completed in cases:
            with self.subTest(completed=completed), self.assertRaises(
                isolation.PRO20RSRC1ForensicUncertainty
            ) as caught:
                isolation.reduce_completed_child(request, completed)
            self.assertFalse(caught.exception.publish_typed_terminal)
            self.assertFalse(caught.exception.close_writer)

    def test_parent_resource_premises_are_exact_operational_limits(self) -> None:
        isolation.check_parent_resource_premises(
            parent_current_rss_bytes=isolation.PARENT_MAX_CURRENT_RSS_BYTES,
            available_memory_bytes=isolation.MIN_AVAILABLE_MEMORY_BYTES,
        )
        for current, available in (
            (
                isolation.PARENT_MAX_CURRENT_RSS_BYTES + 1,
                isolation.MIN_AVAILABLE_MEMORY_BYTES,
            ),
            (0, isolation.MIN_AVAILABLE_MEMORY_BYTES - 1),
        ):
            with self.subTest(current=current, available=available), self.assertRaises(
                isolation.PRO20RSRC1IsolationError
            ):
                isolation.check_parent_resource_premises(
                    parent_current_rss_bytes=current,
                    available_memory_bytes=available,
                )

    def test_live_child_rss_monitor_kills_overrun_as_forensic_uncertainty(self) -> None:
        class Process:
            pid = 123
            returncode = 0

            def __init__(self):
                self.calls = 0
                self.terminated = False

            def communicate(self, *, input, timeout):
                self.calls += 1
                if self.calls == 1:
                    raise isolation.subprocess.TimeoutExpired(("child",), timeout)
                return b"response", b""

            def poll(self):
                return None

            def terminate(self):
                self.terminated = True

            def wait(self, timeout=None):
                return 0

            def kill(self):
                self.terminated = True

        with TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            process = Process()
            with patch.object(isolation.subprocess, "Popen", return_value=process), self.assertRaises(
                isolation.PRO20RSRC1ForensicUncertainty
            ) as caught:
                isolation.launch_child_process(
                    ("child",),
                    repository_root=root,
                    request_raw=b"{}",
                    timeout_seconds=1.0,
                    child_rss=lambda _pid: isolation.CHILD_MAX_PEAK_RSS_BYTES + 1,
                    now=lambda: 0.0,
                )
            self.assertTrue(process.terminated)
            self.assertIn("current-RSS", str(caught.exception))
            self.assertFalse(caught.exception.publish_typed_terminal)

    def test_live_child_rss_monitor_allows_bounded_completion(self) -> None:
        class Process:
            pid = 123
            returncode = 0

            def __init__(self):
                self.calls = 0

            def communicate(self, *, input, timeout):
                self.calls += 1
                if self.calls == 1:
                    raise isolation.subprocess.TimeoutExpired(("child",), timeout)
                return b"response", b""

            def poll(self):
                return None

        clock = iter((0.0, 0.0, 0.1, 0.2)).__next__
        with TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            process = Process()
            with patch.object(isolation.subprocess, "Popen", return_value=process):
                completed = isolation.launch_child_process(
                    ("child",),
                    repository_root=root,
                    request_raw=b"{}",
                    timeout_seconds=1.0,
                    child_rss=lambda _pid: 1,
                    now=clock,
                )
            self.assertEqual(completed.stdout, b"response")
            self.assertEqual(completed.returncode, 0)

    def test_child_current_rss_observer_is_strict(self) -> None:
        with patch.object(
            isolation.subprocess,
            "run",
            return_value=SimpleNamespace(returncode=0, stdout="321\n", stderr=""),
        ):
            self.assertEqual(isolation.child_process_rss_bytes(123), 321 * 1024)
        with patch.object(
            isolation.subprocess,
            "run",
            return_value=SimpleNamespace(returncode=0, stdout="ambiguous", stderr=""),
        ), self.assertRaises(isolation.PRO20RSRC1IsolationError):
            isolation.child_process_rss_bytes(123)

    def test_child_seam_builds_only_the_requested_member_and_writes_nothing(self) -> None:
        request = _request()
        calls: list[str] = []

        def builder(_origin, *, member_key, static_input_bytes, source_closure_sha256):
            calls.append(member_key)
            self.assertEqual(static_input_bytes, {})
            self.assertEqual(source_closure_sha256, request.source_closure_sha256)
            return SimpleNamespace(
                captured_origin_sha256=request.origin_capture_sha256,
                generation1_bridge_content_id=request.generation1_bridge_content_id,
                source_closure_sha256=request.source_closure_sha256,
                source_configuration_sha256=request.source_configuration_sha256,
                checkpoint=_checkpoint(),
            )

        with TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            raw = isolation.execute_child_request(
                root,
                isolation.encode_attempt_request(request),
                origin_loader=lambda _root: object(),
                static_loader=lambda _root: {},
                member_builder=builder,
                peak_rss=lambda: 1,
                now=iter((0.0, 1.0)).__next__,
            )
            rebound, _metrics = isolation.classify_child_response(request, raw)
            self.assertEqual(rebound.kind, "accepted_fine")
            self.assertEqual(calls, [request.member_key])
            self.assertEqual(list(root.iterdir()), [])

    def test_child_module_does_not_import_the_production_store_or_runner(self) -> None:
        path = ROOT / "src/recursive_horizons/fgc/evolution/pro20_rsrc1_isolation.py"
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
        source = path.read_text(encoding="utf-8")
        self.assertNotIn("publish_attempt", source)
        self.assertNotIn("publish_seed", source)
        self.assertNotIn("write_bytes", source)

    def test_child_cli_has_no_status_live_store_or_synthetic_switch(self) -> None:
        path = ROOT / "scripts/run_fgc_pro20_rsrc1_child.py"
        source = path.read_text(encoding="utf-8")
        tree = ast.parse(source)
        imports: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                imports.append(node.module or "")
            elif isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
        joined = " ".join(imports)
        self.assertNotIn("pro20_ev1_store", joined)
        self.assertNotIn("campaign_store", joined)
        for forbidden in (
            "argparse",
            "--run",
            "--status",
            "--live",
            "synthetic",
            "publish_",
            "write_bytes",
        ):
            self.assertNotIn(forbidden, source)


if __name__ == "__main__":
    unittest.main()
