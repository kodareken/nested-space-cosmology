#!/usr/bin/env python3
"""Fail closed: the historical TI1 authority produced an invalid run."""

from __future__ import annotations

import json


def main() -> int:
    print(
        json.dumps(
            {
                "artifact_id": "FGC-1-TDG9-TI1-FRZ1",
                "classification": "historical_invalid_authority_retired",
                "detail": (
                    "TI1 is immutable evidence; use only a prospectively "
                    "authorized successor"
                ),
                "state_advance_authorized": False,
            },
            sort_keys=True,
            separators=(",", ":"),
        )
    )
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
