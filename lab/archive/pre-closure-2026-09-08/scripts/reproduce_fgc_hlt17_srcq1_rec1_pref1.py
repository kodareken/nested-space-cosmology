#!/usr/bin/env python3
"""Verify compact PREF1 by default; reconstruct independently only on --live.

Default/--check is raw-, runs-, source-, Git-, and shadow-blind. It validates
tracked compact bytes only and never reconstructs or writes. Explicit --live
performs the independent no-write reproduction and prints a small compact
result/config. It never publishes, serializes an endpoint, or writes state.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution import (  # noqa: E402
    hlt17_srcq1_rec1_pref1_binder as binder,
)


FORBIDDEN_LIVE_IMPORTS = (
    "hlt17_srcq1_rec1_auth1",
    "qualify_fgc_hlt17_srcq1_rec1",
    "publish_qualification_result",
    "qualify_and_publish",
    "publish_exclusive",
)


def _check() -> int:
    try:
        result = binder.verify_compact(ROOT)
    except binder.HLT17SRCQ1REC1PREF1Error as error:
        sys.stderr.write(f"{error}\n")
        return 2
    print(
        "FGC-1-HLT17-SRCQ1-REC1-PREF1 compact result verified; "
        f"{result['artifact_payload']['classification']}; "
        "no state, endpoint, store, or physics claim was opened"
    )
    return 0


def _live() -> int:
    imported = set(sys.modules)
    result = binder.bind_live(ROOT)
    leaked = [
        name
        for name in sys.modules
        if name not in imported
        and any(needle in name for needle in FORBIDDEN_LIVE_IMPORTS)
    ]
    if leaked:
        raise binder.HLT17SRCQ1REC1PREF1Error(
            "live PREF1 imported REC1 runner/authority or publication code: "
            + ",".join(sorted(leaked))
        )
    if result.get("written") is not False:
        raise binder.HLT17SRCQ1REC1PREF1Error("live PREF1 attempted a write")
    sys.stdout.buffer.write(binder.render_config(result["config"]))
    sys.stdout.buffer.write(binder.canonical_result(result["result"]))
    print(
        "FGC-1-HLT17-SRCQ1-REC1-PREF1 live binder printed compact result/config; "
        "nothing was written; PRO20 first-event authority remains a separate freeze",
        file=sys.stderr,
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument(
        "--check",
        action="store_true",
        help="compact, raw/runs/source/Git/shadow-blind check (default)",
    )
    action.add_argument(
        "--live",
        action="store_true",
        help="independent no-write reproduction; print compact result/config",
    )
    arguments = parser.parse_args(argv)
    if arguments.live:
        return _live()
    return _check()


if __name__ == "__main__":
    raise SystemExit(main())
