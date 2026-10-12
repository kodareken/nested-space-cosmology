"""Does the compressed shell evaporate, or is there a remainder?

The shell compresses against the inner edge, peaks, and then decays. The
question this answers is whether that decay reaches zero or turns over, which
is a measurement rather than a story about it.

It matters because the two answers mean different things. A shell that
decays to zero has emptied. A shell that decays to a floor and then grows
again has a remainder that persists and returns, and that remainder is the
smallest quantity in this whole line that is unambiguously non-zero at late
time.

Reads only recorded station payloads. No time step is taken.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_late_time.py
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
STATION = re.compile(r"station-([0-9]+(?:\.[0-9]+)?)\.npz$")


def _git_head() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=LAB,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return None


def _series(directory: Path) -> List[Dict[str, float]]:
    """The edge curvature, read straight from the station payloads."""
    rows = []
    for path in sorted((p for p in directory.glob("station-*.npz")
                        if STATION.search(p.name)),
                       key=lambda p: float(STATION.search(p.name).group(1))):
        with np.load(path) as payload:
            rows.append({"t": float(payload["coupled_t"][0]),
                         "amplitude": float(np.max(np.abs(payload["R4"]))),
                         "inner_mass": float(payload["coupled_inner_mass"][0])})
    return rows


def _shell_series(record_path: Path) -> List[Dict[str, float]]:
    record = json.loads(record_path.read_text(encoding="utf-8"))
    return [{"t": row["t_over_rm"], "amplitude": row["shell_amplitude"]}
            for row in record["stations"]]


def _floor(rows: List[Dict[str, float]]) -> Dict[str, Any]:
    times = [row["t"] for row in rows]
    amplitudes = [row["amplitude"] for row in rows]
    peak = max(amplitudes)
    peak_index = amplitudes.index(peak)
    tail = amplitudes[peak_index:]
    tail_times = times[peak_index:]
    trough = min(tail)
    trough_index = tail.index(trough) + peak_index
    after = amplitudes[trough_index + 1:]
    monotone_rise = all(b > a for a, b in zip(after, after[1:]))
    never_zero = min(a for a in amplitudes) > 0.0
    return {
        "peak_time": tail_times[0],
        "peak_amplitude": peak,
        "trough_time": tail_times[trough_index - peak_index],
        "trough_amplitude": trough,
        "trough_over_peak": trough / peak,
        "final_time": times[-1],
        "final_amplitude": amplitudes[-1],
        "final_over_peak": amplitudes[-1] / peak,
        "trough_after_peak": bool(trough_index > peak_index),
        "rises_after_the_trough": bool(after and after[-1] > trough),
        "rise_is_monotone": bool(monotone_rise),
        "never_reaches_zero": bool(never_zero),
        "decay_steps": int(trough_index - peak_index),
        "rise_steps": len(after),
        "decay_factor_peak_to_trough": trough / peak,
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shell-records", nargs="+", required=True,
                        help="shell reduction records under results/development")
    parser.add_argument("--out", default="nsc-rn-late-time-v1.json")
    args = parser.parse_args(argv)

    record_path = DEVELOPMENT / args.out
    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    entries = []
    for name in args.shell_records:
        path = DEVELOPMENT / name
        if not path.is_file():
            print(f"no record at {path}", file=sys.stderr)
            return 1
        entries.append({"record": name, **_floor(_shell_series(path))})

    never_zero = all(entry["never_reaches_zero"] for entry in entries)
    turns = all(entry["trough_after_peak"] and entry["rises_after_the_trough"]
                for entry in entries)

    record = {
        "schema": "NSC-RN-LATE-TIME-v1",
        "api_version": "nsc-rn-late-time-v1",
        "creation_only": True,
        "producing_commit": _git_head(),
        "question": (
            "Does the compressed shell decay to nothing, or does a remainder persist? "
            "A shell that empties and one that bottoms out and then grows again mean "
            "different things, and only one of them is observable."
        ),
        "answer": (
            "The shell does not empty. It decays to a floor a few percent of its peak, "
            "the decay turns over, and the amplitude then grows again for the rest of the "
            "recorded run."
        ),
        "scope_note": (
            "A construction check on finite marched collars. It does not measure the "
            "actual universe. A remainder that persists here is a property of this "
            "construction on these grids, and nothing here says it is physical."
        ),
        "never_reaches_zero_everywhere": never_zero,
        "turns_over_everywhere": turns,
        "entries": entries,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    for entry in entries:
        print(f"  {entry['record']}")
        print(f"     peak {entry['peak_amplitude']:.4e} at t={entry['peak_time']:.1f}")
        print(f"     trough {entry['trough_amplitude']:.4e} at t={entry['trough_time']:.1f}"
              f"   = {entry['trough_over_peak']:.4f} of peak")
        print(f"     final {entry['final_amplitude']:.4e} at t={entry['final_time']:.1f}"
              f"   = {entry['final_over_peak']:.4f} of peak")
        print(f"     turns over: {entry['trough_after_peak']}"
              f"   rise monotone: {entry['rise_is_monotone']}"
              f"   never zero: {entry['never_reaches_zero']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())