"""Parent-to-child response in the RN chart, matched against the nested pair.

The claim ledger records, for the responsive parent construction, that a
+5 percent exterior occupation change moves original child-column content by
1.114 percent of its baseline, with held forecasts matching to 0.36-0.37
percent. That is a parent-to-child transfer coefficient of 0.2228.

This measures the same quantity on the RN radial chart, whose child interval
is nsc_rn_reference's areal interval near the outer horizon and whose parent
is the large-radius exterior. Two marches, one baseline and one with the
exterior occupation raised by 5 percent, read at matched times.

What is compared is an occupation ratio, not a wave amplitude. The two
constructions use different measures (projector occupation in one, areal
integral in the other), so agreement is judged at the level of magnitude and
sign, not to three figures.

Does not measure the actual universe and does not transfer a number from one
construction to the other.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from recursive_horizons import nsc_rn_mass_shell as mass_shell
from recursive_horizons import nsc_rn_reference as reference


def load(run_dir: Path):
    stations = sorted(run_dir.glob("station-*.npz"))
    out = {}
    for path in stations:
        with np.load(path, allow_pickle=False) as data:
            t = float(data["t"])
            out[t] = {
                "file": path,
                "mass": np.array(data["mass_profile"], dtype=float),
                "r": np.array(data["radius"], dtype=float),
            }
    return out


def child_content(t, child_lo: float, child_hi: float) -> float:
    """Occupation in the child interval, with the canonical half-density measure."""
    r = t["r"]
    m = t["mass"]
    band = (r >= child_lo) & (r <= child_hi)
    if not np.any(band):
        return 0.0
    return float(np.sum(m[band]))


def exterior_content(t, child_hi: float, r_out: float) -> float:
    r = t["r"]
    m = t["mass"]
    band = (r > child_hi) & (r <= r_out)
    return float(np.sum(m[band]))


def matched(baseline: dict, perturbed: dict) -> list[float]:
    common = sorted(set(baseline) & set(perturbed))
    return common


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    out_path = root / "results" / "development" / "nsc-rn-parent-child-response-v1.json"
    if out_path.exists():
        raise SystemExit(f"refusing to overwrite {out_path}")

    baseline_dir = Path(r"C:/Users/swedi/AppData/Local/Temp/nsc-rn-pcr-base")
    perturbed_dir = Path(r"C:/Users/swedi/AppData/Local/Temp/nsc-rn-pcr-plus5")
    if not (baseline_dir.is_dir() and perturbed_dir.is_dir()):
        raise SystemExit("run directories missing; see --prepare")

    geom = reference.case_geometry(mu=1.04)
    child_lo = float(geom["child"]["start"])
    child_hi = float(geom["child"]["end"])
    r_out = float(geom["outer_edge"])

    base = load(baseline_dir)
    pert = load(perturbed_dir)
    times = matched(base, pert)
    if not times:
        raise SystemExit("no matched stations")

    rows = []
    for t in times:
        cb = child_content(base[t], child_lo, child_hi)
        cp = child_content(pert[t], child_lo, child_hi)
        eb = exterior_content(base[t], child_hi, r_out)
        ep = exterior_content(pert[t], child_hi, r_out)
        if cb <= 0.0 or eb <= 0.0:
            continue
        rows.append(
            {
                "t": t,
                "child_baseline": cb,
                "child_perturbed": cp,
                "exterior_baseline": eb,
                "exterior_perturbed": ep,
                "child_fraction": (cp - cb) / cb,
                "exterior_fraction": (ep - eb) / eb,
            }
        )
    if not rows:
        raise SystemExit("no usable matched stations")

    # Late rows carry the response; early rows are still the release transient.
    tail = rows[len(rows) // 2 :]
    ratios = [
        r["child_fraction"] / r["exterior_fraction"]
        for r in tail
        if abs(r["exterior_fraction"]) > 1e-12
    ]
    mean_ratio = float(np.mean(ratios)) if ratios else float("nan")
    spread = float(np.std(ratios)) if len(ratios) > 1 else float("nan")

    ledger_ratio = 1.114 / 5.0

    record = {
        "record_type": "rn_parent_child_response",
        "provenance": {
            "source": str(Path(__file__).name),
            "baseline": str(baseline_dir),
            "perturbed": str(perturbed_dir),
            "child_interval": [child_lo, child_hi],
            "outer_edge": r_out,
            "matched_times": len(rows),
        },
        "child_interval": [child_lo, child_hi],
        "rows": rows,
        "tail_row_count": len(tail),
        "rn_transfer_ratio": mean_ratio,
        "rn_transfer_spread": spread,
        "ledger_transfer_ratio": ledger_ratio,
        "ledger_source": "claim-ledger.md responsive parent: +5% parent occupation -> 1.114% child content",
        "comparison_note": (
            "Different measures in the two constructions: projector occupation "
            "in the nested pair, areal integral here. Compared at magnitude and sign."
        ),
    }
    record["sha256"] = hashlib.sha256(
        json.dumps(record, sort_keys=True).encode("utf-8")
    ).hexdigest()

    out_path.write_text(json.dumps(record, indent=2, sort_keys=True))
    print(f"wrote {out_path}")
    print(f"matched times          {len(rows)}")
    print(f"RN transfer ratio      {mean_ratio:+.4f}  spread {spread:.4f}")
    print(f"ledger transfer ratio  {ledger_ratio:+.4f}")
    print(f"ratio RN/ledger        {mean_ratio / ledger_ratio:+.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())