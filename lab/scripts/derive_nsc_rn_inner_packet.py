"""Reduce the inner-packet march to the fate of a shell released inside the collar.

Reads the stations written by ``run_nsc_rn_inner_packet.py``, locates the
mass shell in each, and fits its inward motion. No time step is taken and
nothing is re-run.

The point of the measurement is the part the recorded neutral family could
not reach. There the packet is released at 12 r_m and by the end of every
run 61% to 99% of the accreted mass has crossed below the grid's inner
edge at r = 1.04, so the accumulation is outside the domain. Here the same
packet is released inside the collar, the shell forms where it can be seen,
and the profile it leaves behind is resolved rather than truncated.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_inner_packet.py
"""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_rn_mass_shell as shell  # noqa: E402
from recursive_horizons import nsc_rn_source as source  # noqa: E402

SOURCE = LAB / "results/development/nsc-rn-inner-packet-v1"
RECORD = LAB / "results/development/nsc-rn-inner-packet-shell-v2.json"
STATION = re.compile(r"station-([0-9]+(?:\.[0-9]+)?)\.npz$")


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_head() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=LAB,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return None


def _clock(path: Path) -> float:
    match = STATION.search(path.name)
    return float(match.group(1))


def _stations(directory: Path) -> List[Path]:
    found = [p for p in directory.glob("station-*.npz")
             if STATION.search(p.name)]
    if not found:
        raise SystemExit(f"no station payloads found in {directory}")
    return sorted(found, key=_clock)


def main() -> int:
    if RECORD.exists():
        print(f"refusing to overwrite sealed record: {RECORD}", file=sys.stderr)
        return 1
    run_summary_path = SOURCE / "summary.json"
    if not run_summary_path.is_file():
        print(f"missing march summary: {run_summary_path}", file=sys.stderr)
        return 1
    run_summary = json.loads(run_summary_path.read_text(encoding="utf-8"))

    base_mass = float(run_summary["mu"])
    rows: List[Dict[str, Any]] = []
    for station in _stations(SOURCE):
        with np.load(station) as payload:
            state = shell.shell_state(payload, base_mass=base_mass)
            # The station name carries t/r_m at the recorded 0.1 resolution.
            # It is the clock, not coupled_t scaled by the mass ratio: r_m is
            # the magnetic radius and is not mu.
            state["t_over_rm"] = _clock(station)
            state["clock_raw"] = float(payload["coupled_t"][0])
        rows.append(state)

    motion = shell.shell_motion(rows)

    # What happens after the shell reaches the inner edge. In the recorded
    # family this never happened inside the domain, so it is the part of the
    # picture that reduction could not deliver.
    arrival_t = motion["inner_edge_arrival_t"]
    post = [row for row in rows if arrival_t is not None and row["t_over_rm"] >= arrival_t]
    accumulation = None
    if post:
        amps = [row["shell_amplitude"] for row in post]
        best = max(range(len(amps)), key=lambda i: amps[i])
        accumulation = {
            "arrival_t": arrival_t,
            "amplitude_at_arrival": amps[0],
            "peak_amplitude": amps[best],
            "peak_t": post[best]["t_over_rm"],
            "amplitude_gain": amps[best] / amps[0] if amps[0] else None,
            "final_amplitude": amps[-1],
            "amplitude_after_peak_fraction": amps[-1] / amps[best] if amps[best] else None,
            "stations_after_arrival": len(post),
        }

    # The quantity the recorded family could not deliver: how much mass is
    # actually inside the collar, and what profile it leaves behind.
    total = rows[-1]["total_gain"]
    final = rows[-1]
    record = {
        "schema": "NSC-RN-INNER-PACKET-SHELL-v1",
        "api_version": shell.API_VERSION,
        "creation_only": True,
        "producing_commit": _git_head(),
        "march_summary_sha256": _sha256_file(run_summary_path),
        "construction": {
            "mu": run_summary["mu"],
            "eps": run_summary["eps"],
            "packet_center_over_rm": run_summary["packet_center_over_rm"],
            "packet_width_over_rm": run_summary["packet_width_over_rm"],
            "recorded_family_packet_center_over_rm":
                run_summary["recorded_family_packet_center_over_rm"],
            "r_inner": run_summary["r_inner"],
            "r_outer": run_summary["r_outer"],
            "r_minus": run_summary["r_minus"],
            "points": run_summary["points"],
            "target_over_rm": run_summary["target_over_rm"],
            "stop_reason": run_summary["stop_reason"],
        },
        "what_this_measures": (
            "The fate of the same packet that the recorded family drops from 12 r_m, "
            "released instead inside the collar so that the shell forms and falls where "
            "the grid can resolve it. The point is the accumulation profile, which the "
            "recorded family could not reach."
        ),
        "scope_note": (
            "A construction check on one finite marched collar. It does not measure the "
            "actual universe. It is a different construction from the recorded family: "
            "same packet width, same spectrum, same chart, different release radius, so "
            "its numbers are not comparable with the recorded family's except through "
            "quantities that depend only on the construction."
        ),
        "method": {
            "shell_density": "centred dM/dr of the recorded coupled_mass",
            "smooth_window": shell.SMOOTH_WINDOW,
            "search_limit": shell.SEARCH_LIMIT,
            "feature_continuity": shell.FEATURE_CONTINUITY,
            "coherence_test": (
                "position against the inner edge plus amplitude continuity. Deliberately "
                "not an amplitude test: a shell piling up at the boundary has a sharper "
                "density peak than the same shell in flight, so an amplitude test accepts "
                "the pinned state and rejects the moving one."
            ),
        },
        "mass_budget": {
            "total_gain": total,
            "gain_below_inner_edge": final["gain_below_inner_edge"],
            "fraction_below_inner_edge": (
                final["gain_below_inner_edge"] / total if total else None),
        },
        "motion": motion,
        "accumulation_after_arrival": accumulation,
        "stations": rows,
    }

    RECORD.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {RECORD.relative_to(LAB)}")
    print(f"  stations         : {len(rows)}  (falling {motion['stations_falling']})")
    print(f"  shell start r    : {motion.get('radius_first', float('nan')):.4f}")
    print(f"  shell end r      : {motion.get('radius_last', float('nan')):.4f}")
    print(f"  infall rate      : {motion.get('infall_rate_mean', float('nan')):.4f}")
    print(f"  r^2              : {motion.get('r_squared', float('nan')):.5f}")
    print(f"  edge arrival t   : {motion['inner_edge_arrival_t']}")
    if accumulation:
        print(f"  accumulation     : amplitude {accumulation['amplitude_at_arrival']:.3e}"
              f" -> peak {accumulation['peak_amplitude']:.3e}"
              f" at t={accumulation['peak_t']}"
              f"  (x{accumulation['amplitude_gain']:.3f})")
        print(f"                     then {accumulation['amplitude_after_peak_fraction']:.3f}"
              f" of the peak at the final station")
    print(f"  fraction of gain below the inner edge at the final station: "
          f"{record['mass_budget']['fraction_below_inner_edge']:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())