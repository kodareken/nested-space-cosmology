#!/usr/bin/env python3
"""Apply known mixed heat coefficients to the declared compact chiral domain.

Only the new adjoint/domain and boundary-coefficient map is calculated.
Previously recorded bulk coefficients and warp integrals are authenticated
and reused. The determinant phase and transmitting throat remain separate.
"""
import argparse
import json
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from recursive_horizons.nsc_spinor_bridge import weyl_matrices
from check_nsc_vacuum_charge_matching import compare, hashes

OUTPUT = ROOT / "results/development/compact-boundary-action.json"
SOURCES = (
    "scripts/check_nsc_compact_boundary_action.py",
    "docs/nsc-compact-boundary-action.md",
    "src/recursive_horizons/nsc_spinor_bridge.py",
    "scripts/check_nsc_vacuum_charge_matching.py",
    "docs/nsc-compact-mass-map.md",
)
INPUTS = ("results/development/vacuum-charge-matching.json",
          "results/development/charged-sector.json")


def authenticated_record(path, seen=None):
    seen = set() if seen is None else seen
    if path in seen:
        return json.loads((ROOT / path).read_text())
    seen.add(path)
    value = json.loads((ROOT / path).read_text())
    for field in ("source_hashes", "input_hashes"):
        declared = value.get(field, {})
        compare(declared, hashes(tuple(declared)), f"{path}/{field}")
    for dependency in value.get("input_hashes", {}):
        if dependency.endswith(".json"):
            authenticated_record(dependency, seen)
    return value


def calculate():
    previous = authenticated_record(INPUTS[0])
    authenticated_record(INPUTS[1])
    frame = weyl_matrices()
    g5, beta = frame["gamma5"], frame["gamma"][0]
    tangent = [beta, *(sp.I * g for g in frame["gamma"][1:])]
    kron = sp.kronecker_product
    normal = kron(sp.eye(2), g5)
    current = kron(sp.eye(2), sp.I * beta * g5)
    tangent = [kron(sp.eye(2), g) for g in tangent]
    sheet = kron(sp.diag(1, -1), sp.eye(4))
    rank, identity = 8, sp.eye(8)
    k_entries = sp.symbols("k00 k01 k02 k03 k11 k12 k13 k22 k23 k33", real=True)
    curvature = sp.zeros(4)
    index = 0
    for a in range(4):
        for b in range(a, 4):
            curvature[a, b] = curvature[b, a] = k_entries[index]
            index += 1
    trace_k = sp.trace(curvature)
    norm_k = sp.trace(curvature * curvature)
    residuals, domains = {}, {}

    def zero(name, value):
        if isinstance(value, sp.MatrixBase):
            assert value.applyfunc(sp.simplify) == sp.zeros(*value.shape), name
        else:
            assert sp.simplify(value) == 0, name
        residuals[name] = "0"

    for orientation in (-1, 1):
        prefix = f"orientation_{orientation}"
        allowed = (identity - orientation * sheet * normal) / 2
        adjoint = identity - allowed
        chi = 2 * allowed - identity
        zero(prefix + "/Lorentzian_endpoint_current", allowed * current * allowed)
        zero(prefix + "/adjoint_boundary_pairing", adjoint * normal * allowed)
        green_rank = (allowed * normal * allowed).rank()
        assert green_rank == 4  # Nonzero: this is NOT a self-adjoint Euclidean D.
        for a, gamma in enumerate(tangent):
            zero(prefix + f"/projected_tangent_{a}", allowed * gamma * allowed)
        s = -trace_k * allowed / 2  # inward Robin convention; K is outward trace
        a1_trace = sp.trace(chi)
        a2_boundary = sp.trace(2 * trace_k * identity + 12 * s) / 6
        # The Euclidean boundary projector follows the normal spin frame.
        chi_derivatives = [orientation * sheet * sum(
            (curvature[a, b] * tangent[b] for b in range(4)), sp.zeros(8))
            for a in range(4)]
        derivative_trace = sum(sp.trace(x * x) for x in chi_derivatives)
        zero(prefix + "/projector_derivative_trace", derivative_trace - rank * norm_k)
        a3_numerator = sp.trace(
            (13 * allowed - 7 * adjoint) * trace_k ** 2
            + (2 * allowed + 10 * adjoint) * norm_k
            + 96 * s * trace_k + 192 * s * s) - 12 * derivative_trace
        zero(prefix + "/a1", a1_trace)
        zero(prefix + "/a2_boundary", a2_boundary + rank * trace_k / 6)
        zero(prefix + "/a3_boundary", a3_numerator - 3 * rank * (trace_k ** 2 - 2 * norm_k))
        domains[prefix] = {"allowed_rank": int(allowed.rank()),
                           "Euclidean_same_domain_Green_form_rank": int(green_rank),
                           "Euclidean_D_self_adjoint_on_this_domain": False,
                           "adjoint_uses_complementary_values": True}

    # Reuse the conformal Dirac law; only the normal/Robin conversion is new.
    y = sp.Symbol("Y", real=True)
    length = sp.Symbol("L_star", positive=True)
    sigma = -sp.Rational(18, 1015) * (y / length) ** 2
    field = sp.Function("chi")(y)
    psi = sp.exp(-2 * sigma) * field
    endpoint_residuals = {}
    for outward in (-1, 1):
        k_out = 4 * outward * sp.exp(-sigma) * sp.diff(sigma, y)
        robin = -outward * sp.exp(-sigma) * sp.diff(psi, y) - k_out * psi / 2
        zero(f"normal_{outward}/Robin_conformal_map",
             robin + outward * sp.exp(-3 * sigma) * sp.diff(field, y))
        # Apply the existing Dirichlet metric momentum to this endpoint;
        # do not infer free-boundary stationarity from the GHY coefficient.
        coefficient = sp.simplify((-3 * k_out * length / 4).subs(y, outward * length))
        target = sp.Rational(108, 1015) * sp.exp(sp.Rational(18, 1015))
        zero(f"normal_{outward}/metric_boundary_residual", coefficient - target)
        endpoint_residuals[str(outward)] = str(coefficient)

    a5, boundary, gradient = sp.symbols("A5 B J", real=True)
    zero("bulk_boundary_warp_cancellation",
         a5 * (8 * boundary - 12 * gradient) - 8 * a5 * boundary + 12 * a5 * gradient)
    k_ghy = previous["specified_warp_check"]["K_with_standard_GHY"]
    # No old warp integration or spectrum is executed again.
    correction = 5 * k_ghy / 3
    return {
        "schema": "NSC-COMPACT-BOUNDARY-ACTION-v1",
        "status": "compact determinant-magnitude boundary matching; full source open",
        "source_hashes": hashes(SOURCES), "input_hashes": hashes(INPUTS),
        "domain": {"geometry": "Euclidean conformal collar e^(2 sigma)*(g4+dY^2), interval [-L_star,L_star]",
                   "sigma": "-18*(Y/L_star)^2/1015",
                   "allowed_projector": "P=(I8-eta*tau3 tensor gamma5)/2, eta=+1 or -1",
                   "fields": "two torsionless Dirac copies; geometric terms at zero link and gauge backgrounds",
                   "Hilbert_measure": "integral sqrt(g5) Psi_dagger Psi",
                   "Lorentzian_Cauchy_current": "I2 tensor (i beta gamma5)",
                   "Euclidean_normal_gamma": "outward_sign * I2 tensor gamma5",
                   "positive_operator": "D_E^dagger D_E, not a self-adjoint D_E squared on the same chiral domain",
                   "mixed_domain": "(I-P)Psi=0; (nabla_in - K_out/2)P Psi=0",
                   "K_convention": "trace of extrinsic curvature with outward normal; Vassilevich L_aa",
                   "boundary_type": "reflecting compact endpoints, not the parent-child transmitting throat"},
        "domains": domains, "exact_residuals": residuals,
        "imported_heat_formulas": {
            "source": "https://arxiv.org/html/hep-th/0306138v3, equations 5.30-5.32",
            "rank": rank,
            "a1": "0",
            "a2": "-rank/[12*(4*pi)^(5/2)] * (integral_M R + 2*integral_boundary K)",
            "a3_geometric": "rank/[128*(4*pi)^2] * integral_boundary (K^2-2*K_ab*K_ab)",
            "a3_on_this_umbilic_collar": "rank/[256*(4*pi)^2] * integral_boundary K^2 >= 0",
            "assumptions": "unit heat smearing; E=-R/4; local mixed conditions; smooth Euclidean metric; no torsion or link potential"},
        "proper_time_boundary_action": {
            "a2_window_weight": "(Lambda^3-nu_match^3)/3",
            "a3_window_weight": "(Lambda^2-nu_match^2)/2",
            "Einstein_boundary_to_bulk_coefficient_ratio": "2",
            "order_Lambda4_boundary_volume_term": "0",
            "boundary_a3_vacuum_density": "N*(Lambda^2-nu_match^2)/(8*(4*pi)^2) * sum_endpoints exp(2*sigma)*(sigma_prime)^2",
            "boundary_a3_sign_on_static_umbilic_collar": "nonnegative; N=2 here",
            "V_warp_with_derived_a2_boundary": "-12*A5*J",
            "warp_ratio_from_authenticated_input": "-(f3/f5)*K_with_standard_GHY/zeta",
            "reused_K_with_standard_GHY": k_ghy,
            "maximum_negative_relative_warp_factor": correction,
            "retained_Xi_lower_bound_for_zeta_ge_1_in_units_q_squared": 3 * (1 - correction)},
        "metric_boundary_matching": {
            "source": "https://arxiv.org/html/1605.01603v1, equations 4-5",
            "Dirichlet_momentum": "pi^ab=-A5*sqrt(abs(h))*epsilon*(K^ab-K*h^ab)",
            "leading_source_free_boundary_condition_if_metric_varied": "K_ab-K*h_ab=0",
            "dimensionless_geometric_residual_per_endpoint": endpoint_residuals,
            "residual_definition": "L_star*(K_ab-K*h_ab)/h_ab on the umbilic collar, no component sum",
            "common_interface_equation": "pi_parent^ab+pi_child^ab+delta Gamma_rest/delta h_ab=0",
            "Gamma_rest_excludes_already_included_geometric_terms": True,
            "current_Gaussian_endpoints_stationary_with_no_other_boundary_terms": False,
            "higher_order_or_transmission_metric_variation_solved": False},
        "scope": {"compact_geometric_a1_a2_a3_derived_for_declared_modulus_domain": True,
                  "Euclidean_D_dagger_D_is_positive": True,
                  "full_compact_boundary_action_derived": False,
                  "Dirac_determinant_phase_derived": False,
                  "transmitting_throat_boundary_action_derived": False,
                  "quantum_state_or_relative_scale_selected": False,
                  "full_renormalized_vacuum_fixed": False,
                  "self_sourced_solution_derived": False,
                  "higher_order_curvature_boundary_coefficients_included": False},
        "comparison": {"fields": "all", "exact": "identities, strings, types, keys and hashes",
                       "float_atol": 3e-13, "float_rtol": 3e-13, "exceptions": []},
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
        print("compact boundary action: all fields and dependencies reproduced; prior generators were not run")
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
