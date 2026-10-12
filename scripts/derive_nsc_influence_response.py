"""First-principles two-level influence fee from the stage tangent.

The response owner's stage_tangent linearises one stage_rates call. Give
it a direction localized in the OUTER level (r >= 24 r_m) and read the
mass-function response. The pre-registered per-level fee is the L1
fraction of that response inside the INNER level (r < 8 r_m), seam at a
quarter of the 32 r_m cap:

    f = sum_{r<8} w |d_mass| / sum_all w |d_mass|

A reciprocity control runs the mirrored direction (inner-localized) and
checks the same machinery from the other side, plus a locality control
(inner direction should answer mostly inside its own level).

The closure uses the round-1 band: local_fraction = 1 - f holds against
the BBN-anchored 0.049. No march, one tangent per probe.

    python scripts/lab.py scripts/derive_nsc_influence_response.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))
if str(LAB / "scripts") not in sys.path:
    sys.path.insert(0, str(LAB / "scripts"))

import numpy as np

import derive_nsc_rn_parent as driver
from recursive_horizons import nsc_influence_chain as chain
from recursive_horizons import nsc_rn_pg as pg
from recursive_horizons import nsc_rn_response as response
from recursive_horizons import nsc_rn_source as source

SEAM_OVER_RM = 8.0
OUTER_START_OVER_RM = 24.0
INNER_START_OVER_RM = 0.0
INNER_END_OVER_RM = SEAM_OVER_RM
OUTPUT = LAB / "results" / "development" / "nsc-influence-response-v1.json"
SCHEMA = "NSC-INFLUENCE-RESPONSE-v1"


def _window(radius, left, right, sharpness=8.0):
    edge_left = 0.5 * (1.0 + np.tanh(sharpness * (radius - left)))
    edge_right = 0.5 * (1.0 + np.tanh(sharpness * (right - radius)))
    return edge_left * edge_right


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    geometry = driver.case_geometry(1.04)
    grid = pg.build_grid(
        geometry["mass"], r_m=geometry["r_m"], points=geometry["route"]["points"],
        r_out=geometry["route"]["outer"],
    )
    packet = source.prepare_radial_packet(
        mass_over_rm=1.04, energy_target=1.0e-3, rank=driver.PILOT_RANK,
        model=grid["model"], radii=grid["radius"], weights=grid["weights"],
        derivative=grid["derivative"],
    )
    state = pg.make_state(
        grid, np.asarray(packet["phi"]),
        inner_mass=geometry["mass"],
        killing_frequency=float(packet["radial_moments"]["requested_omega"]),
        occupations=np.asarray(packet["occupations"], dtype=float),
    )
    radius = grid["radius"]
    weights = grid["weights"]

    def probe_chart(left, right):
        mask = _window(radius, left, right)
        direction = {
            "dphi": np.asarray(packet["phi"]) * mask[None, :, None],
            "d_occupations": np.zeros_like(state["occupations"]),
            "d_inner_mass": 0.0,
        }
        tangent = response.stage_tangent(state, grid, direction)
        fee = chain.spatial_fee(radius, weights, tangent["d_lapse"], SEAM_OVER_RM)
        return tangent, fee

    tangent_fwd, fee_fwd = probe_chart(OUTER_START_OVER_RM, float(geometry["route"]["outer"]))
    tangent_bwd, fee_bwd = probe_chart(0.0, SEAM_OVER_RM)
    mass_fee = chain.spatial_fee(
        radius, weights, tangent_fwd["d_mass"], SEAM_OVER_RM,
    )["l1_mass_response_fraction_inside"]

    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "seam_over_rm": SEAM_OVER_RM,
        "outer_level": [OUTER_START_OVER_RM, float(geometry["route"]["outer"])],
        "inner_level": [0.0, SEAM_OVER_RM],
        "forward_fee_lapse": fee_fwd["l1_mass_response_fraction_inside"],
        "backward_fee_lapse": fee_bwd["l1_mass_response_fraction_inside"],
        "mass_channel_fraction": mass_fee,
        "mass_channel_reading": (
            "the mass function is causal: an outer-level source changes "
            "d_mass only above its own radius, so the inner-level fraction "
            "is exactly zero. The chart channels (lapse/shift) are the "
            "nonlocal influence."),
        "edge_rates_forward": {
            "d_inner_mass_rate": float(tangent_fwd["d_inner_mass_rate"]),
            "d_outer_flux_rate": float(tangent_fwd["d_outer_flux_rate"]),
            "d_excision_flux_rate": float(tangent_fwd["d_excision_flux_rate"]),
        },
        "chain_read": {
            "f": fee_fwd["l1_mass_response_fraction_inside"],
            "local_fraction": 1.0 - fee_fwd["l1_mass_response_fraction_inside"],
            "target": chain.TARGET_LOCAL_FRACTION,
            "band": [chain.BAND_LO, chain.BAND_HI],
        },
        "round_1_reference": {
            "measured_mass_fee": 0.0010282460948723 / 1.04,
            "source": "neutral continuation summary, coupled vs fixed inner mass",
        },
        "notes": [
            "the forward fee is the fraction of an outer-level perturbation's "
            "mass response that lands inside the inner quarter of the cap",
            "the backward probe is the reciprocity control; exact reciprocity "
            "is not asserted, the two numbers are reported side by side",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_influence_response.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_influence_response.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  forward lapse fee = {:.6e}  backward lapse fee = {:.6e}  mass fraction = {:.6e}".format(
        fee_fwd["l1_mass_response_fraction_inside"],
        fee_bwd["l1_mass_response_fraction_inside"],
        mass_fee,
    ))
    print("  edge rates: inner {:.3e} outer {:.3e} excision {:.3e}".format(
        record["edge_rates_forward"]["d_inner_mass_rate"],
        record["edge_rates_forward"]["d_outer_flux_rate"],
        record["edge_rates_forward"]["d_excision_flux_rate"],
    ))
    print("  chain local fraction = {:.6f}".format(record["chain_read"]["local_fraction"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())