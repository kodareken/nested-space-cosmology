from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sys
import tomllib
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]
from recursive_horizons.fgc.evolution.proto18_historical_replay import (
    Proto18HistoricalReplayError, replay_proto12_history, replay_rsp2_history,
)


def source(name: str) -> tuple[bytes, bytes, dict]:
    root = ROOT / "runs/fgc-2-sf1" / name
    external = (root / "events.jsonl").read_bytes()
    with np.load(root / "latest-checkpoint.npz", allow_pickle=False) as archive:
        embedded = bytes(archive["event_log_utf8"])
        metadata = json.loads(bytes(archive["metadata_utf8"]))
    return external, embedded, metadata


def selector_hashes() -> dict[str, str]:
    with (ROOT / "configs/fgc/fgc-1-pro18-frz1.toml").open("rb") as handle:
        config = tomllib.load(handle)
    return {
        item["key"]: item["physical_state_sha256"]
        for item in config["selectors"]["member"]
    }


class Proto18HistoricalReplayTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.proto12 = source("proto12/calibration")
        cls.rsp2 = source("rsp2/amplitude3-ssprk3-momentum")
        cls.hashes = selector_hashes()

    def replay_proto12(self, *args):
        return replay_proto12_history(
            *args,
            {key: value for key, value in self.hashes.items() if key != "SSPRK3-16385"},
        )

    def replay_rsp2(self, *args):
        return replay_rsp2_history(*args, self.hashes["SSPRK3-16385"])

    def test_frozen_histories_replay_exactly(self) -> None:
        p = self.replay_proto12(*self.proto12)
        r = self.replay_rsp2(*self.rsp2)
        self.assertEqual(p["selected_coordinate_time"], "23/16")
        self.assertEqual(r["endpoint_coordinate_time"], "23/16")

    def test_bytes_duplicates_and_embedding_fail_closed(self) -> None:
        raw, embedded, metadata = self.proto12
        with self.assertRaises(Proto18HistoricalReplayError):
            self.replay_proto12(raw, embedded[:-1], metadata)
        # Duplicate keys must fail at parsing rather than silently choosing the last.
        line = b'{"amplitude":"3","amplitude":"3"}\n'
        with self.assertRaises(Proto18HistoricalReplayError):
            from recursive_horizons.fgc.evolution.proto18_historical_replay import parse_jsonl_exact
            parse_jsonl_exact(line, label="test", expected_bytes=len(line), expected_sha256=__import__("hashlib").sha256(line).hexdigest(), expected_lines=1)
        nested = b'{"outer":{"x":1,"x":2}}\n'
        with self.assertRaises(Proto18HistoricalReplayError):
            parse_jsonl_exact(nested, label="test", expected_bytes=len(nested), expected_sha256=__import__("hashlib").sha256(nested).hexdigest(), expected_lines=1)
        nonfinite = b'{"x":NaN}\n'
        with self.assertRaises(Proto18HistoricalReplayError):
            parse_jsonl_exact(nonfinite, label="test", expected_bytes=len(nonfinite), expected_sha256=__import__("hashlib").sha256(nonfinite).hexdigest(), expected_lines=1)
        overflow = b'{"x":1e999}\n'
        with self.assertRaises(Proto18HistoricalReplayError):
            parse_jsonl_exact(overflow, label="test", expected_bytes=len(overflow), expected_sha256=__import__("hashlib").sha256(overflow).hexdigest(), expected_lines=1)
        with self.assertRaises(Proto18HistoricalReplayError):
            self.replay_proto12(raw[:-1] + b" \n", raw[:-1] + b" \n", metadata)

    def test_proto12_selected_record_and_metadata_links_fail_closed(self) -> None:
        raw, embedded, metadata = self.proto12
        altered = deepcopy(metadata)
        altered["completed_amplitude_records"][1]["last_completed_common_event_index"] = 22
        with self.assertRaises(Proto18HistoricalReplayError): self.replay_proto12(raw, embedded, altered)
        altered = deepcopy(metadata); altered["members"]["RK4-2049"]["time"] = 1.5
        with self.assertRaises(Proto18HistoricalReplayError): self.replay_proto12(raw, embedded, altered)
        altered = deepcopy(metadata)
        altered["completed_amplitude_records"][1]["member_final_state_hashes"]["RK4-2049"] = "0" * 64
        with self.assertRaises(Proto18HistoricalReplayError): self.replay_proto12(raw, embedded, altered)

    def test_rsp2_order_terminal_and_metadata_links_fail_closed(self) -> None:
        raw, embedded, metadata = self.rsp2
        lines = raw.splitlines(); first, second = lines[1], lines[2]; lines[1], lines[2] = second, first
        reordered = b"\n".join(lines) + b"\n"
        with self.assertRaises(Proto18HistoricalReplayError): self.replay_rsp2(reordered, reordered, metadata)
        altered = deepcopy(metadata); altered["terminal_result"]["classification"] = "invalid"
        with self.assertRaises(Proto18HistoricalReplayError): self.replay_rsp2(raw, embedded, altered)
        altered = deepcopy(metadata); altered["terminal_result"]["endpoint_assessment"] = {}
        with self.assertRaises(Proto18HistoricalReplayError): self.replay_rsp2(raw, embedded, altered)
        altered = deepcopy(metadata); altered["members"]["SSPRK3-16385"]["state_sha256"] = "0" * 64
        with self.assertRaises(Proto18HistoricalReplayError): self.replay_rsp2(raw, embedded, altered)


if __name__ == "__main__":
    unittest.main()
