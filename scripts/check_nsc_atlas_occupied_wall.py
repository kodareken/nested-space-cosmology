"""Finite state-to-stress-to-observer calculation; no field grid."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import mpmath as mp

from recursive_horizons.nsc_atlas_occupied_wall import OccupiedWall, filled_shell

ROOT = Path(__file__).resolve().parents[1]


def calculation():
    started, wall_started = time.process_time(), time.perf_counter()
    w = OccupiedWall()
    states = [w.static_state(q) for q in (999, 1000, 1001)]
    signals = [w.central_parent_signal(s) for s in states]
    max_junction = max(abs(s[k]) for s in states for k in
                       ("angular_junction_residual", "temporal_junction_residual",
                        "potential_residual", "mass_stationarity_residual"))
    if max_junction > 1e-12:
        raise RuntimeError("signed equilibrium check failed")
    h = 1e-4
    clock_slope = (w.central_parent_signal(states[1], emission_time=h)["reception_time"]
                   -w.central_parent_signal(states[1], emission_time=-h)["reception_time"])/(2*h)
    with mp.workdps(50):
        h2, t = mp.mpf(49)/625, mp.mpf(1)/225
        gh = (mp.mpf(7)/45)/filled_shell(1000)["energy_coefficient"]
        d = gh*filled_shell(1001)["energy_coefficient"]
        def mass(r):
            ms = t*r*r+d/r
            return h2*r**3/2+mp.sqrt(1-h2*r*r)*ms-ms*ms/(2*r)
        r = mp.findroot(lambda x: mp.diff(mass, x), (.99, 1.02))
        root_gap = states[2]["radius"]-float(r)
        mass_gap = states[2]["mass"]-float(mass(r))
    clock_gap = clock_slope-signals[1]["one_plus_redshift"]
    if abs(root_gap) > 1e-12 or abs(clock_gap) > 2e-9:
        raise RuntimeError("independent root or clock comparison failed")
    # For an exactly copied normalized sphere spectrum, N and A stay fixed.
    # At fixed G and hbar the dimensionless gravitational state changes as R^-2.
    separation = 12.
    clearance = separation-signals[0]["initial_join_radius"]-signals[2]["initial_join_radius"]
    if clearance <= 0 or min(s["V_second"] for s in states) <= 0:
        raise RuntimeError("two-room geometric or radial-stability condition failed")
    if time.process_time()-started > 120:
        raise RuntimeError("finite calculation exceeded its CPU budget")
    paths = (
        "src/recursive_horizons/nsc_atlas_occupied_wall.py",
        "src/recursive_horizons/nsc_atlas_parent_observables.py",
        "scripts/check_nsc_atlas_occupied_wall.py",
        "scripts/plot_nsc_atlas_occupied_wall.py",
        "tests/test_nsc_atlas_occupied_wall.py",
    )
    return dict(
        schema="nsc-atlas-occupied-wall-v1",
        candidate="Einstein de Sitter core / occupied worldvolume Dirac membrane / Schwarzschild / compensated flat dust parent",
        action_data=dict(bare_k_Rref=w.bare_k, H_squared_Rref_squared=w.H_squared,
                         planck_length_squared_over_Rref_squared=w.planck_length_squared,
                         species=w.species, parent_H0_Rref=w.parent_H0),
        preparation="Two parity-paired species, positive levels1..q fully occupied. No occupation fit to sky data. Normal-ordered excitation stress; fixed renormalized surface tension, no extra surface curvature terms.",
        reference=dict(level=1000, beta_minus="24/25", beta_plus="4/5",
                       gas_compactness="7/45", bare_tension_compactness="1/225",
                       radial_V_second_Rref_squared="35159/7200",
                       radial_omega_squared_Rref_squared="35159/14400",
                       Rref_over_planck_length=1/math.sqrt(w.planck_length_squared),
                       Lambda_cap_Rref_squared=3*w.H_squared),
        states=states, central_to_cosmological_observer=signals,
        two_room_parent=dict(centres=[[-6.,0.,0.],[6.,0.,0.]], levels=[999,1001],
                             initial_comoving_clearance=clearance,
                             observer_locations=[[-6.,8.,0.],[6.,8.,0.]],
                             common_density=signals[0]["parent_density"],
                             statement="Disjoint comoving vacuoles share one exact dust exterior for tau>=0. Ordinary horizonless orientation; not black-hole offspring formation."),
        checks=dict(max_signed_junction_residual=max_junction,
                    neighboring_radius_50_digit_gap=root_gap,
                    neighboring_mass_50_digit_gap=mass_gap,
                    proper_clock_derivative_gap=clock_gap,
                    max_two_stability_formula_gap=max(abs(s["V_second"]-s["V_second_mass_check"]) for s in states)),
        inheritance=dict(conserved="N and its spectral coefficient along the conformal wall history",
                         unselected="Initial occupation, surface coupling and cap vacuum; branching transfer is not constructed.",
                         homothety="At fixed constants, a copy with R->aR and the same occupation has bare_k*R multiplied by a and ell_P^2*A/R^2 divided by a^2."),
        limits="Exact thin-surface mean-stress model and linear radial stability. No claim of nonspherical/fragmentation stability, microscopic localization, formation, generated cap Lambda, cosmic fractions or baryon asymmetry.",
        numerical_uncertainty="Independent arithmetic/root and clock differences only; no microscopic-model error bound or observational fit.",
        evaluated_spectral_states=3, field_grids=0,
        process_cpu_seconds=time.process_time()-started,
        elapsed_seconds=time.perf_counter()-wall_started,
        source_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        reproduce="python scripts/run_python.py scripts/check_nsc_atlas_occupied_wall.py",
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    encoded = json.dumps(calculation(), indent=2, allow_nan=False)+"\n"
    if len(encoded.encode()) > 18000:
        raise RuntimeError("compact record exceeds 18KB")
    if args.output:
        with args.output.open("x", encoding="utf-8", newline="\n") as f:
            f.write(encoded)
    print(encoded)


if __name__ == "__main__":
    main()
