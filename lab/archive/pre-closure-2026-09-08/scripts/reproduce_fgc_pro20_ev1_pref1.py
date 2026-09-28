#!/usr/bin/env python3
"""Verify compact PREF1 by default; reconstruct independently only on --live.

Default/--check is Git-, raw-, store-, source-, Planck-, and replay-blind. It
validates tracked compact bytes only and never reconstructs or writes.
Explicit --live performs the independent no-write reproduction and prints a
small compact result/config. The binder and certificate never publish,
serialize an additional endpoint, or write state. The production campaign
store is the event's closed terminal.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution import (  # noqa: E402
    pro20_ev1_pref1_certificate as certificate,
)


FORBIDDEN_LIVE_IMPORTS = (
    "pro20_ev1_runtime",
    "pro20_ev1_protocol",
    "pro20_ev1_store",
    "pro20_ev1_authority",
    "run_fgc_pro20_ev1",
)


def _check() -> int:
    try:
        result = certificate.verify_compact(ROOT)
    except certificate.PRO20EV1PREF1CertificateError as error:
        sys.stderr.write(f"{error}\n")
        return 2
    print(
        "FGC-1-PRO20-EV1-PREF1 compact result verified; "
        f"{result['artifact_payload']['classification']}; "
        "event-written production store with 146 accepted fines; "
        "binder/certificate wrote nothing; "
        "no common event, calibration, PRO21/Wave2, or physics claim was opened"
    )
    return 0


def _live() -> int:
    imported = set(sys.modules)
    result = certificate.bind_live(ROOT)
    leaked = [
        name
        for name in sys.modules
        if name not in imported
        and any(needle in name for needle in FORBIDDEN_LIVE_IMPORTS)
    ]
    if leaked:
        raise certificate.PRO20EV1PREF1CertificateError(
            "live PREF1 imported PRO20 runner/protocol/store/authority code: "
            + ",".join(sorted(leaked))
        )
    if result.get("written") is not False:
        raise certificate.PRO20EV1PREF1CertificateError(
            "live PREF1 attempted a write"
        )
    sys.stdout.buffer.write(certificate.render_config(result["config"]))
    sys.stdout.buffer.write(certificate.canonical_result(result["result"]))
    print(
        "FGC-1-PRO20-EV1-PREF1 live certificate printed compact result/config; "
        "binder/certificate wrote nothing; production campaign store remains "
        "the event terminal; RSRC1 remains a planning owner only",
        file=sys.stderr,
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument(
        "--check",
        action="store_true",
        help="compact, Git/raw/store/source/Planck/replay-blind check (default)",
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
