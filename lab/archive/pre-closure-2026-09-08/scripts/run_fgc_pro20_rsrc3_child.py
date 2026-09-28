#!/usr/bin/env python3
"""One isolated RSRC3 seed or scheduled-attempt child."""

from __future__ import annotations

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SOURCE = str(ROOT / "src")
if SOURCE not in sys.path:
    sys.path.insert(0, SOURCE)

from recursive_horizons.fgc.evolution.pro20_rsrc3_profile import (  # noqa: E402
    activate_rsrc3,
)


activate_rsrc3()

from recursive_horizons.fgc.evolution.pro20_rsrc1_attempt import (  # noqa: E402
    classify_attempt,
)
from recursive_horizons.fgc.evolution.pro20_rsrc1_isolation import (  # noqa: E402
    MAX_HANDOFF_BYTES,
    current_child_peak_rss_bytes,
    decode_attempt_request,
    encode_child_response,
    execute_child_request,
    parse_canonical_handoff,
)
from recursive_horizons.fgc.evolution.pro20_rsrc1_seed import (  # noqa: E402
    SEED_REQUEST_SCHEMA,
    decode_seed_request,
    encode_seed_memory_stop,
    execute_seed_child_request,
)


def _memory_response(request_raw: bytes) -> bytes:
    schema = parse_canonical_handoff(request_raw, "child request").get("schema")
    if schema == SEED_REQUEST_SCHEMA:
        return encode_seed_memory_stop(
            decode_seed_request(request_raw),
            peak_rss_bytes=current_child_peak_rss_bytes(),
        )
    request = decode_attempt_request(request_raw)
    classification = classify_attempt(
        member_key=request.member_key,
        predecessor_bundle=request.predecessor_bundle,
        successor_bundle=request.predecessor_bundle,
        accepted_state_advanced=False,
        resource_stop={
            "reason": "isolated RSRC3 child memory allocation failed",
            "scope": "one_scheduled_attempt",
        },
    )
    return encode_child_response(
        request,
        classification,
        peak_rss_bytes=current_child_peak_rss_bytes(),
        wall_seconds=0.0,
    )


def main() -> int:
    request_raw = sys.stdin.buffer.read(MAX_HANDOFF_BYTES + 1)
    if len(request_raw) > MAX_HANDOFF_BYTES:
        return 2
    try:
        schema = parse_canonical_handoff(request_raw, "child request").get("schema")
        if schema == SEED_REQUEST_SCHEMA:
            response = execute_seed_child_request(
                ROOT,
                request_raw,
                peak_rss=current_child_peak_rss_bytes,
            )
        else:
            response = execute_child_request(ROOT, request_raw)
    except MemoryError:
        try:
            response = _memory_response(request_raw)
        except Exception:
            return 1
    except Exception:
        return 1
    sys.stdout.buffer.write(response)
    sys.stdout.buffer.flush()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
