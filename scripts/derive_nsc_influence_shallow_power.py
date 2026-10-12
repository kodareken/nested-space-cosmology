"""Shallow-power and boundary-channel closure from recorded NSC bytes.

Rows of the second channel of the nested-influence question:

    1. What power nu does 1 - R^-nu = target require at the recorded
       route ratio R = 32 (the outer parent cap)?
    2. What power does the charged lowest-Landau screen actually give at
       the matched coupling g4 = sqrt(g_squared) from the scale binding
       branch (mu^2 = |q| g4^2 / (4 pi^2 r^2) -> r^-nu, nu = g4 sqrt|q|/2pi)?
    3. The measured boundary-channel terms from the neutral continuation
       (outer lapse rate, outer mass flux, pressure and lapse integrals).

Saved-data consumer; creation-only output.

    python scripts/lab.py scripts/derive_nsc_influence_shallow_power.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

SCALE_BINDING = LAB / "results" / "development" / "scale-binding.json"
CONTINUATION = (
    LAB / "results" / "development" / "nsc-rn-neutral-continuation-v1" / "summary.json"
)
OUTPUT = LAB / "results" / "development" / "nsc-influence-shallow-power-v1.json"
SCHEMA = "NSC-INFLUENCE-SHALLOW-POWER-v1"


def _load_json(path):
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def main():
    for path in (SCALE_BINDING, CONTINUATION):
        if not path.is_file():
            raise SystemExit("missing input record: " + str(path))
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))

    scale_record = _load_json(SCALE_BINDING)
    continuation_record = _load_json(CONTINUATION)

    branch = scale_record["development_branch"]
    q = float(branch["magnetic_flux"])
    g_squared = float(branch["mmp_map"]["g_squared"])
    g4 = g_squared ** 0.5

    screen = chain.lll_screen_power(g4, q)
    required = chain.required_power(chain.ROUTE_RATIO)
    band_nu = chain.band_required_power(chain.ROUTE_RATIO)
    required_coupling = chain.coupling_for_power(required["required_nu"], q)
    actual_closure = chain.chain_fraction(screen["nu"], chain.ROUTE_RATIO)
    needed_ratio = {
        "target_local": chain.TARGET_LOCAL_FRACTION,
        "per_level_decay": required["per_level_decay"],
        "ratio_with_matched_nu": required["per_level_decay"] ** (-1.0 / screen["nu"]),
        "reading": (
            "the layer ratio R such that R^-nu(matched) = f_target; R > 1 means "
            "nearly degenerate nesting is required at the matched coupling"
        ),
    }

    final = continuation_record["final"]
    outer_flux = float(final["coupled_mass_stocks"]["outer_flux"])
    inner_mass_final = float(final["coupled_inner_mass"])
    boundary_fee = abs(outer_flux) / inner_mass_final
    boundary = {
        "final_outer_lapse_rate": float(final["outer_lapse_rate"]),
        "final_outer_flux": outer_flux,
        "final_inner_mass": inner_mass_final,
        "per_level_boundary_mass_leak": boundary_fee,
        "leak_over_target_decay": boundary_fee / required["per_level_decay"],
        "structural": (
            "N_t = 0 at the absorbing outer boundary in this run: the lapse "
            "channel cannot carry a uniform term into the vacuum chart"
        ),
        "final_pressure_integral": float(final["pressure_integral"]),
        "final_lapse_integral": float(final["lapse_integral"]),
        "final_metric_power": float(final["metric_power"]),
    }

    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "band": [chain.BAND_LO, chain.BAND_HI],
        "target_local_fraction": chain.TARGET_LOCAL_FRACTION,
        "route_ratio": chain.ROUTE_RATIO,
        "route_ratio_source": "ROUTE_OUTER_OVER_RM = 32.0 in lab/scripts/derive_nsc_rn_parent.py",
        "required_power_at_route_ratio": required,
        "band_required_power": band_nu,
        "required_g4_for_target": required_coupling,
        "recorded_scale_binding": {
            "magnetic_flux": q,
            "g_squared": g_squared,
            "g4": g4,
            "omega": branch["omega"],
            "zeta": branch["zeta"],
        },
        "charged_screen": screen,
        "actual_closure_at_matched_coupling": actual_closure,
        "needed_ratio": needed_ratio,
        "boundary_channel": boundary,
        "input_sha256": {
            str(SCALE_BINDING.relative_to(LAB)): chain.path_sha256(SCALE_BINDING),
            str(CONTINUATION.relative_to(LAB)): chain.path_sha256(CONTINUATION),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_influence_shallow_power.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_influence_shallow_power.py"
            ),
        },
        "notes": [
            "The LLL screen power is derived from mu^2 = |q| g4^2 / (4 pi^2 r^2), "
            "the bosonization scale already in the papers; nothing is fitted.",
            "All numbers are read from frozen records except the required-power "
            "inversion, which is pure algebra closed on recorded constants.",
        ],
    }
    encoded = json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n"
    OUTPUT.write_text(encoded, encoding="utf-8")
    print("wrote", OUTPUT.relative_to(LAB))
    print(
        "  required nu @R=32 for target: {:.6f}  band [{:.6f}, {:.6f}]".format(
            required["required_nu"], band_nu["nu_min"], band_nu["nu_max"]
        )
    )
    print(
        "  matched g4 = {:.4f} (q={:g}) -> nu = {:.6f} -> f = {:.6e} local = {:.6f} {}".format(
            g4, q, screen["nu"], actual_closure["f"],
            actual_closure["local_fraction"], actual_closure["verdict"],
        )
    )
    print(
        "  g4 demanded by target at q={:g}: {:.5f}  (matched/required = {:.1f}x)".format(
            q, required_coupling["g4"], g4 / required_coupling["g4"],
        )
    )
    print(
        "  boundary: outer_lapse_rate = {:.3g}, per-level mass leak = {:.3e}, "
        "leak/f_target = {:.3e}".format(
            boundary["final_outer_lapse_rate"],
            boundary["per_level_boundary_mass_leak"],
            boundary["leak_over_target_decay"],
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())