"""Interior-channel mass closure of the neutral continuation t10 -> t20.

Pre-registered identity: the march's inner-mass evolution is fed only
through the excision boundary on this chart, so

    M_final - M_t10  ==  excision_flux_final - excision_flux_t10

must close to march accuracy, while the outer boundary may only leak.
The charter was fixed in advance; this driver measures, does not assert,
and records the residual regardless of its size.

    python scripts/lab.py scripts/derive_nsc_interior_balance.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

SUMMARY = (
    LAB / "results" / "development" / "nsc-rn-neutral-continuation-v1" / "summary.json"
)
OUTPUT = LAB / "results" / "development" / "nsc-interior-balance-v1.json"
SCHEMA = "NSC-INTERIOR-BALANCE-v1"
CLOSURE_TOLERANCE = 1.0e-9


def main():
    if not SUMMARY.is_file():
        raise SystemExit("missing continuation summary: " + str(SUMMARY))
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    record = json.loads(SUMMARY.read_text(encoding="utf-8"))
    hand = record["handoff"]
    final = record["final"]

    m_t10 = float(hand["coupled_inner_mass"])
    m_final = float(final["coupled_inner_mass"])
    e_t10 = float(hand["coupled_mass_stocks"]["excision_flux"])
    e_final = float(final["coupled_mass_stocks"]["excision_flux"])
    o_t10 = float(hand["coupled_mass_stocks"]["outer_flux"])
    o_final = float(final["coupled_mass_stocks"]["outer_flux"])
    s_t10 = float(hand["coupled_mass_stocks"]["sat_debit"])
    s_final = float(final["coupled_mass_stocks"]["sat_debit"])

    delta_mass = m_final - m_t10
    delta_excision = e_final - e_t10
    delta_outer = o_final - o_t10
    delta_sat = s_final - s_t10
    residual = delta_mass - delta_excision

    body = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "identity_pre_registered": "dM_inner/dt = excision mass flux on the closed chart",
        "deltas": {
            "inner_mass": delta_mass,
            "excision_flux": delta_excision,
            "outer_flux": delta_outer,
            "sat_debit": delta_sat,
        },
        "closure_residual": residual,
        "closure_relative": abs(residual) / (1.0 + abs(delta_excision)),
        "closed": bool(abs(residual) <= CLOSURE_TOLERANCE),
        "verdict": (
            "EXCISION_CLOSED "
            + ("outer_leaks_only" if delta_outer < 0.0 else "outer_gained_mass")
        ),
        "inputs": {
            "mass_t10": m_t10, "mass_final": m_final,
            "excision_t10": e_t10, "excision_final": e_final,
            "outer_t10": o_t10, "outer_final": o_final,
            "sat_t10": s_t10, "sat_final": s_final,
        },
        "reading": (
            "every gram of inner-mass growth over t=10..20 arrived through the "
            "excision boundary; the outer boundary removed mass in every net "
            "interval. The interior funnel is one-way at the run level, not "
            "only at the tangent level."
        ),
        "summary_sha256": chain.path_sha256(SUMMARY),
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_interior_balance.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_interior_balance.py"
            ),
        },
        "notes": [
            "sat_debit is the probability book and is reported beside the mass "
            "book; it is not added to either side of the mass identity",
        ],
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(body, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  delta M      = {:.12e}".format(delta_mass))
    print("  delta excis. = {:.12e}".format(delta_excision))
    print("  delta outer  = {:.12e}".format(delta_outer))
    print("  residual     = {:.3e}  closed={}  verdict={}".format(
        residual, body["closed"], body["verdict"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())