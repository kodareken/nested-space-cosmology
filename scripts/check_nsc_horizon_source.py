#!/usr/bin/env python3
"""Apply the known Dirac/CFT source to the unwrapped NSC horizon benchmark."""
import argparse
import json
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from recursive_horizons.nsc_horizon_source import (
    conformal_stress, pg_source, reconstructed_neck,
    state_constants, static_reduced_equations,
)
from check_nsc_compact_boundary_action import authenticated_record
from check_nsc_vacuum_charge_matching import compare, hashes

OUTPUT = ROOT / "results/development/horizon-source.json"
SOURCES = ("scripts/check_nsc_horizon_source.py",
           "src/recursive_horizons/nsc_horizon_source.py",
           "docs/nsc-horizon-source.md",
           "src/recursive_horizons/nsc_lorentzian.py",
           "docs/nsc-lorentzian-transport.md",
           "scripts/check_nsc_compact_boundary_action.py",
           "scripts/check_nsc_vacuum_charge_matching.py")
INPUTS = ("results/nsc-5-clock-horizon.json",
          "results/development/charged-sector.json",
          "results/development/compact-casimir.json")


def schwarzian(mapping, variable):
    first = sp.diff(mapping, variable)
    return sp.diff(mapping, variable, 3) / first - sp.Rational(3, 2) * (sp.diff(mapping, variable, 2) / first) ** 2


def calculate():
    records = [authenticated_record(path) for path in INPUTS]
    clock = records[0]
    residuals = {}

    def zero(name, value):
        simplified = sp.simplify(value)
        assert simplified == 0, (name, simplified)
        residuals[name] = "0"

    x = sp.Symbol("rho", real=True)
    central, kappa = sp.symbols("c_eff kappa_h", positive=True)
    A = sp.Function("A")(x)
    t_u, t_v = sp.symbols("t_u t_v", real=True)
    uu, uv, vv = conformal_stress(A, sp.diff(A, x), sp.diff(A, x, 2), t_u, t_v, central, sp.pi)
    Ttt, Trr = (uu + 2 * uv + vv) / A, -(uu - 2 * uv + vv) / A
    zero("twoD_trace_anomaly", Ttt + Trr - central * sp.diff(A, x, 2) / (24 * sp.pi))
    zero("twoD_radial_Ward", sp.diff(Trr, x) + sp.diff(A, x) / (2 * A) * (Trr - Ttt))
    beta, fp, fpp, radius = sp.symbols("beta A_prime A_second radius", positive=True)
    f = 1 - beta ** 2
    components = conformal_stress(f, fp, fpp, t_u, t_v, central, sp.pi)
    pg = pg_source(f, beta, radius, *components, sp.pi)
    zero("PG_Killing_current", pg["parent_Killing_power"] - (t_u - t_v))

    theta, q = sp.symbols("theta q_mag", real=True)
    flux = 2 * sp.pi * sp.integrate(q * sp.sin(theta) / 2, (theta, 0, sp.pi))
    zero("magnetic_flux_normalization", flux - 2 * sp.pi * q)
    zero("source_free_angular_Maxwell_density", sp.diff(q / (2 * radius ** 2), theta))

    u = sp.Symbol("u", real=True)
    scale = sp.Symbol("scale", positive=True)
    offset = sp.Symbol("offset", real=True)
    affine = scale * u + offset
    ray = -sp.exp(-kappa * u) / kappa
    zero("affine_state_map_has_zero_Schwarzian", schwarzian(affine, u))
    zero("exponential_peeling_Schwarzian", schwarzian(ray, u) + kappa ** 2 / 2)
    entropy = -central * sp.log(sp.diff(ray, u)) / 12
    flux_ray = -central * schwarzian(ray, u) / (24 * sp.pi)
    zero("energy_entropy_rate", flux_ray - (6 * sp.diff(entropy, u) ** 2 / central
                                          + sp.diff(entropy, u, 2)) / (2 * sp.pi))
    zero("entropy_rate", sp.diff(entropy, u) - central * kappa / 12)
    # Standard composition law, applied without equating it to a physical
    # parent-child clock map that has not been derived.
    z = sp.Symbol("z", real=True)
    outer, inner = sp.Function("U")(z), sp.Function("V")(u)
    zero("state_map_composition",
         schwarzian(outer.subs(z, inner), u)
         - sp.diff(inner, u) ** 2 * schwarzian(outer, z).subs(z, inner)
         - schwarzian(inner, u))

    d, a1, a2, a3 = sp.symbols("d a1 a2 a3", real=True, nonzero=True)
    series_A = a1 * d + a2 * d ** 2 / 2 + a3 * d ** 3 / 6
    outgoing = conformal_stress(series_A, sp.diff(series_A, d), sp.diff(series_A, d, 2),
                                central * a1 ** 2 / (192 * sp.pi), 0, central, sp.pi)[0]
    horizon_limit = sp.limit(outgoing / series_A ** 2, d, 0)
    zero("future_horizon_regular_outgoing_limit", horizon_limit - central * a3 / (192 * sp.pi * a1))

    # Reuse the spherical Einstein equations; check their assembly with this
    # source and the missing-potential interface, not another curvature tensor.
    f0, f1, r, r1, xi, U0, U1, U2, ks = sp.symbols(
        "A0 A1 r r1 xi U0 U1 U2 k_state", nonzero=True, real=True)
    f2, r2, constraint = static_reduced_equations(f0, f1, r, r1, xi, U0, U1, ks)
    derivative = (sp.diff(constraint, f0) * f1 + sp.diff(constraint, f1) * f2
                  + sp.diff(constraint, r) * r1 + sp.diff(constraint, r1) * r2
                  + sp.diff(constraint, U0) * U1 * r1)
    zero("static_constraint_propagates", derivative)
    zero("static_angular_equation", 2 * f0 * r * r2 + 2 * f1 * r * r1 + r ** 2 * f2 + r * U1)
    zero("static_null_difference_equation", -4 * f0 ** 2 * r * r2
         - xi * (2 * f0 * f2 - f1 ** 2 + 4 * ks ** 2))
    # Check against the standard independent mixed Einstein components
    # and the conformal tensor, not just the solved ODE expressions.
    st_u, st_v = state_constants(ks, "Hartle-Hawking", central, sp.pi)
    s_uu, s_uv, s_vv = conformal_stress(f0, f1, f2, st_u, st_v, central, sp.pi)
    G_Newton = 12 * sp.pi * xi / central
    source_tt = U0 / r ** 2 + 2 * G_Newton * (s_uu + 2 * s_uv + s_vv) / (r ** 2 * f0)
    source_rr = U0 / r ** 2 - 2 * G_Newton * (s_uu - 2 * s_uv + s_vv) / (r ** 2 * f0)
    geometric_tt = (1 - f0 * r1 ** 2 - 2 * f0 * r * r2 - f1 * r * r1) / r ** 2
    geometric_rr = (1 - f0 * r1 ** 2 - f1 * r * r1) / r ** 2
    geometric_angle = -f0 * r2 / r - f1 * r1 / r - f2 / 2
    zero("independent_tt_source_assembly", r ** 2 * (geometric_tt - source_tt) - constraint)
    zero("independent_rr_source_assembly", r ** 2 * (geometric_rr - source_rr) - constraint)
    zero("independent_angular_source_assembly", geometric_angle - U1 / (2 * r))
    GM, rQ2, Lambda = sp.symbols("GM rQ2 Lambda4", real=True)
    RN = 1 - 2 * GM / r + rQ2 / r ** 2 - Lambda * r ** 2 / 3
    potential = Lambda * r ** 2 + rQ2 / r ** 2
    RN2, RNr2, RNC = static_reduced_equations(
        RN, sp.diff(RN, r), r, 1, 0, potential, sp.diff(potential, r), ks)
    zero("classical_RN_dS_A_equation", RN2 - sp.diff(RN, r, 2))
    zero("classical_RN_dS_radius_equation", RNr2)
    zero("classical_RN_dS_constraint", RNC)
    neck = reconstructed_neck(r, f1, xi, U0, U1, ks)
    zero("neck_constraint", neck["constraint"])
    alpha = xi / r ** 2
    deficit = 1 - U0
    zero("neck_curvature_reconstruction", neck["A_second"] - (-r * U1 - 2 * deficit) / (r ** 2 - xi))
    zero("neck_flare_reconstruction", r * neck["radius_second"]
         - (2 * deficit + alpha * r * U1) / (2 * neck["A"] * (1 - alpha)))

    kh = float(clock["horizon_affine_parameter"]["kappa_h"])
    rh = float(clock["observer_map"]["horizon_rho"])
    radius_h = float(clock["observer_map"]["horizon_areal_radius"])
    power = kh ** 2 / (48 * float(sp.pi))
    temperature = kh / (2 * float(sp.pi))
    A0 = 1 - 3 * sp.pi / 2
    beta0 = sp.sqrt(1 - A0)
    original_neck = []
    for state in ("Boulware", "Unruh", "Hartle-Hawking"):
        tu, tv = state_constants(sp.Float(kh, 40), state, 1, sp.pi)
        source = conformal_stress(A0, 6, -3 * sp.pi, tu, tv, 1, sp.pi)
        result = pg_source(A0, beta0, 1, *source, sp.pi)
        numerical = {key: float(sp.N(sp.simplify(value), 30)) for key, value in result.items()}
        assert numerical["fourD_null_plus"] > 0 and numerical["fourD_null_minus"] > 0
        original_neck.append({"state": state, **numerical})
    geometric_nulls = [float(sp.N(-2 * (1 - beta0) ** 2, 30)),
                       float(sp.N(-2 * (-1 - beta0) ** 2, 30))]
    assert max(geometric_nulls) < 0
    vacuum_neck = (9 * sp.pi ** 2 - 6 * sp.pi - 36) / (192 * sp.pi)
    assert vacuum_neck > 0
    return {
        "schema": "NSC-HORIZON-SOURCE-v1",
        "status": "domain/state match and conditional magnetic Dirac source; full closure open",
        "source_hashes": hashes(SOURCES), "input_hashes": hashes(INPUTS),
        "exact_residuals": residuals,
        "domain": {
            "actual_patch": "R_rho x S2 x I_Y on the stated PG Cauchy slice; no radial endpoint identification",
            "flat_axial_U1_Wilson_modulus_present": False,
            "development_cell": "S1_axial x S2 x I_Y; its AP holonomy result is a different global domain",
            "magnetic_S2_flux_allowed": "H2 includes Z; F=(q_mag/2) sin(theta) dtheta wedge dphi",
            "quantization": "integral F=2*pi*q_mag for unit electric charge",
            "q_mag_selected": False,
            "neutral_branch": "q_mag=0 has no angular zero mode; neutral greybody factors remain required",
            "magnetic_zero_sector": "|q_mag| complex 2D Dirac fields per one massless 4D Dirac",
            "two_bulk_copies_counted_again": False,
            "transparency": "free massless angular-zero sector has transmission one on the unwrapped radial channel",
        },
        "conventions": {"units": "hbar=c_light=kB=L_throat=1; numerical source is per |q_mag|",
                        "signature": "+---; reduced radial metric ds2=A du dv",
                        "curvature": "R2=A_second",
                        "normal_nulls": "l_plus/minus=partial_tau+(-beta+/-1)*partial_rho",
                        "projected_fourD_source": "T_ab^LLL=t_ab/(4*pi*r^2) in the spherical reduced sector",
                        "scope": "free massless magnetic LLL; no link mass, gauge-dynamical gap or higher-mode contribution included"},
        "state_source": {
            "T_uu": "c_eff*(2*A*A_second-A_prime^2)/(192*pi)+t_u",
            "T_vv": "c_eff*(2*A*A_second-A_prime^2)/(192*pi)+t_v",
            "T_uv": "c_eff*A*A_second/(96*pi)",
            "Unruh": "t_u=c_eff*kappa_h^2/(48*pi), t_v=0; future-horizon regular outgoing sector and no incoming parent particles",
            "Hartle_Hawking": "t_u=t_v=c_eff*kappa_h^2/(48*pi); equilibrium incoming bath",
            "Boulware": "t_u=t_v=0; no net flux, but singular on the fixed horizon",
            "full_Killing_power": "4*pi*r^2*T^rho_tau=t_u-t_v",
            "local_matching_term": "a term proportional to g2 is retained as an unresolved potential contribution; it does not change these null or Killing-flux projections",
            "future_Hadamard_theorem_for_full_NSC": False,
        },
        "clock_flux_map": {
            "state_ray_map": "physical propagation from incoming affine vacuum coordinate U to outgoing parent retarded time u",
            "flux": "-c_eff*{U,u}/(24*pi)",
            "entropy": "S_rad=-c_eff*log(U_prime)/12 up to the specified reference constant",
            "peeling": "kappa(u)=-U_second/U_prime; flux=c_eff*(kappa^2+2*kappa_prime)/(48*pi)",
            "composition": "{U(V(u)),u}=V_prime^2*{U,V}+{V,u}",
            "passive_coordinate_change_creates_flux": False,
            "actual_interroom_clock_map_derived": False,
            "global_entropy_or_Page_curve_derived": False,
        },
        "benchmark_Unruh_per_abs_q": {
            "reused_horizon_rho": rh, "reused_horizon_areal_radius": radius_h,
            "reused_surface_gravity": kh, "Hawking_temperature": temperature,
            "parent_Killing_power": power, "outgoing_channel_entropy_rate": kh / 12,
            "future_horizon_null_stress": -power / (4 * float(sp.pi) * radius_h ** 2),
            "power_formula": "|q_mag|*kappa_h^2/(48*pi)",
            "power_in_ordinary_units": "|q_mag|*hbar*c_light^2*kappa_h^2/(48*pi), with kappa_h in inverse length",
            "cosmological_density_rate_Q_derived": False,
        },
        "original_neck_comparison": {
            "per_abs_q": original_neck, "required_Einstein_nulls": geometric_nulls,
            "geometric_null_definition": "G(l_plus/minus,l_plus/minus), before division by 8*pi*G",
            "geometric_CFT_null_base": str(vacuum_neck),
            "magnetic_and_potential_radial_nulls": "0 for static radial Lorentz-invariant potential terms",
            "fixed_benchmark_sourced_by_this_minimal_sector": False,
            "nonzero_Unruh_flux_allows_exact_static_solution_without_other_flux": False,
        },
        "reduced_closure_interface": {
            "xi": "G*c_eff/(12*pi)",
            "U": "8*pi*G*r^2*rho_pot(r); an unmatched source output, not a fitted dark function",
            "example_potential_structure": "Lambda4*r^2+rQ^2/r^2+u0; u0 and other finite matching terms are NOT set to zero",
            "A_second": str(f2), "radius_second": str(r2), "constraint": str(constraint),
            "neck_A": str(sp.factor(neck["A"])),
            "neck_A_second": str(sp.factor(neck["A_second"])),
            "neck_radius_second": str(sp.factor(neck["radius_second"])),
            "conditional_bounce_test": "at r_prime=0, A0<0, 0<alpha=xi/r0^2<1: r_second>0 iff 2*(1-U0)+alpha*r0*U1<0",
            "regular_local_ODE_domain": "A0!=0, r0>0, r0^2!=xi; global horizon and state matching still required",
            "potential_or_coefficients_assigned_to_obtain_bounce": False,
            "global_self_sourced_solution_derived": False,
        },
        "neutral_flux_normalization": {
            "one_complex_4D_Dirac": "sum_{kappa>0} [8*kappa/(2*pi)] integral omega*Gamma_kappa/(exp(omega/T_H)+1) d omega",
            "signs_of_kappa_are_antiparticles": False,
            "magnetic_abs_q_multiplies_neutral_kappa_sum": False,
            "exact_neutral_transmission_computed_here": False,
        },
        "scope": {
            "magnetic_flux_generation_derived": False,
            "interacting_LLL_masslessness_proved": False,
            "LLL_saturates_complete_stress": False,
            "full_finite_cutoff_action_matched": False,
            "Dirac_self_sourced_black_universe_completed": False,
            "new_universal_Hawking_or_CFT_law_claimed": False,
        },
        "comparison": {"fields": "all", "float_atol": 3e-13, "float_rtol": 3e-13,
                       "exact": "keys, types, strings and hashes", "exceptions": []},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        raise FileExistsError("refusing to overwrite recorded evidence")
    expected = json.loads(OUTPUT.read_text()) if args.check else None
    if args.check:
        compare(expected["source_hashes"], hashes(SOURCES), "$/source_hashes")
        compare(expected["input_hashes"], hashes(INPUTS), "$/input_hashes")
    result = calculate()
    if args.check:
        compare(expected, result)
        print("horizon source: every field and dependency reproduced; prior generators were not run")
    elif args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
        print(args.output)
    else:
        print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
