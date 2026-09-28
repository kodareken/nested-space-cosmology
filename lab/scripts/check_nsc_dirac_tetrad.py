#!/usr/bin/env python3
"""Derive the fixed-background PG Dirac connection, current and boundary data."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
import sympy as sp

from check_nsc_scale_closure import compare

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/nsc-4-dirac-tetrad.json"


def zero_matrix(matrix):
    """Assert an exact identity, returning every simplified entry for the record."""
    reduced = matrix.applyfunc(sp.simplify)
    assert reduced == sp.zeros(*reduced.shape), reduced
    return [[str(entry) for entry in row] for row in reduced.tolist()]


def wedge(left, right):
    return left * right.T - right * left.T


def exterior_derivative(form, coordinates):
    return sp.Matrix(4, 4, lambda a, b:
                     sp.diff(form[b], coordinates[a])
                     - sp.diff(form[a], coordinates[b]))


def calculate():
    tau, rho, theta, phi = sp.symbols("tau rho theta phi", real=True)
    coordinates = (tau, rho, theta, phi)
    radius = sp.Function("r")(rho)
    beta = sp.Function("beta")(rho)
    q = sp.diff(radius, rho) / radius
    beta_prime = sp.diff(beta, rho)
    coframe = sp.Matrix([[1, 0, 0, 0], [beta, 1, 0, 0],
                         [0, 0, radius, 0], [0, 0, 0, radius * sp.sin(theta)]])
    forms = [coframe.row(a).T for a in range(4)]
    eta = sp.diag(1, -1, -1, -1)
    inverse = coframe.inv()

    # Lowered Lorentz indices. Raising the first index is required in Cartan's
    # equation, whereas the spin connection uses these lowered forms directly.
    omega = [[sp.zeros(4, 1) for _ in range(4)] for _ in range(4)]
    entries = {(0, 1): -beta_prime * forms[1],
               (0, 2): -beta * q * forms[2],
               (0, 3): -beta * q * forms[3],
               (1, 2): q * forms[2],
               (1, 3): q * forms[3],
               (2, 3): sp.cot(theta) / radius * forms[3]}
    for (a, b), form in entries.items():
        omega[a][b], omega[b][a] = form, -form
    torsion = []
    for a in range(4):
        residual = exterior_derivative(forms[a], coordinates)
        for b in range(4):
            residual += wedge(eta[a, a] * omega[a][b], forms[b])
        torsion.append(zero_matrix(residual))
    for a in range(4):
        for b in range(4):
            zero_matrix(omega[a][b] + omega[b][a])
    metric = sp.simplify(coframe.T * eta * coframe)
    expected_metric = sp.Matrix([[1-beta**2, -beta, 0, 0], [-beta, -1, 0, 0],
                                 [0, 0, -radius**2, 0],
                                 [0, 0, 0, -radius**2*sp.sin(theta)**2]])
    metric_residual = zero_matrix(metric - expected_metric)

    sigma1 = sp.Matrix([[0, 1], [1, 0]])
    sigma2 = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sigma3 = sp.diag(1, -1)
    zero2 = sp.zeros(2)
    gamma = [sp.diag(1, 1, -1, -1)]
    for sigma in (sigma1, sigma2, sigma3):
        gamma.append(zero2.row_join(sigma).col_join((-sigma).row_join(zero2)))
    for a in range(4):
        for b in range(4):
            zero_matrix(gamma[a]*gamma[b]+gamma[b]*gamma[a]-2*eta[a, b]*sp.eye(4))

    # Explicit gamma-matrix contraction of all six connection forms. This
    # verifies the spin connection independently of radial Hermiticity.
    contracted = sp.zeros(4)
    coordinate_gamma = []
    for mu in range(4):
        curved_gamma = sum((inverse[mu, a]*gamma[a] for a in range(4)), sp.zeros(4))
        coordinate_gamma.append(curved_gamma)
        omega_mu = sum((omega[a][b][mu]*gamma[a]*gamma[b]/4
                        for a in range(4) for b in range(4)), sp.zeros(4))
        contracted += curved_gamma * omega_mu
    expected_contracted = (-gamma[0]*(beta_prime/2+beta*q) + gamma[1]*q
                           + gamma[2]*sp.cot(theta)/(2*radius))
    connection_residual = zero_matrix(contracted-expected_contracted)
    # r D (u/r): angular spin-connection term survives; both radial measure
    # terms cancel, including the term containing beta*r'/r.
    rescaled_connection = contracted-coordinate_gamma[1]*q
    rescaling_residual = zero_matrix(rescaled_connection
                                    + gamma[0]*beta_prime/2
                                    - gamma[2]*sp.cot(theta)/(2*radius))
    radial_principal_residual = zero_matrix(-sp.I*gamma[0]*coordinate_gamma[1]
                                           + sp.I*(gamma[0]*gamma[1]-beta*sp.eye(4)))
    shift_connection_residual = zero_matrix(-sp.I*gamma[0]*(rescaled_connection
                                             - gamma[2]*sp.cot(theta)/(2*radius))
                                            - sp.I*beta_prime*sp.eye(4)/2)

    # Generic complex amplitudes and derivatives, with real coefficient data.
    b, bp, w = sp.symbols("b bp w", real=True)
    x1, y1, x2, y2, dx1, dy1, dx2, dy2 = sp.symbols(
        "x1 y1 x2 y2 dx1 dy1 dx2 dy2", real=True)
    u = sp.Matrix([x1+sp.I*y1, x2+sp.I*y2])
    du = sp.Matrix([dx1+sp.I*dy1, dx2+sp.I*dy2])
    c = sigma2-b*sp.eye(2)
    potential = w*sigma1
    dt = -c*du+bp*u/2-sp.I*potential*u
    flux_derivative = (du.conjugate().T*c*u + u.conjugate().T*c*du
                       - bp*u.conjugate().T*u)[0]
    density_derivative = (dt.conjugate().T*u+u.conjugate().T*dt)[0]
    continuity = sp.simplify(density_derivative+flux_derivative)
    assert continuity == 0
    missing_connection_dt = -c*du-sp.I*potential*u
    false_source = sp.simplify((missing_connection_dt.conjugate().T*u
                               + u.conjugate().T*missing_connection_dt)[0]
                              + flux_derivative)
    density = sp.simplify((u.conjugate().T*u)[0])
    assert sp.simplify(false_source+bp*density) == 0
    assert false_source != 0
    explicit_current = 2*(x1*y2-y1*x2)-b*density
    assert sp.simplify((u.conjugate().T*c*u)[0]-explicit_current) == 0
    old_wall_current = sp.simplify((sp.Matrix([[0, 1]])*c*sp.Matrix([0, 1]))[0])
    assert old_wall_current == -b

    throat_beta = sp.sqrt(3*sp.pi/2)
    throat_c = sigma2-throat_beta*sp.eye(2)
    assert throat_c.det() == 3*sp.pi/2-1
    assert bool(throat_c.det() > 0) and bool(throat_c[0, 0] < 0)
    # Any exterior beta in (0,1) has the same inertia. beta=1/2 is an exact
    # representative, not a claim that the rho=8 probe has that exact value.
    representative_boundary_form = sp.diag(sigma2-sp.eye(2)/2, -throat_c)
    exact_boundary_eigenvalues = list(representative_boundary_form.eigenvals())
    assert all(value.is_real for value in exact_boundary_eigenvalues)
    positive = sum(int(bool(value > 0)) for value in exact_boundary_eigenvalues)
    negative = sum(int(bool(value < 0)) for value in exact_boundary_eigenvalues)
    assert (positive, negative) == (3, 1)

    # Independent numerical evaluation of the actual geometry and its horizon.
    def metric_a(value):
        return 1+3*value+3*(1+value*value)*(np.arctan(value)-np.pi/2)

    def derivative_a(value):
        return 6+6*value*(np.arctan(value)-np.pi/2)

    horizon = brentq(metric_a, 0., 8., xtol=5e-15)
    probes = []
    for point in (-8., 0., horizon, 8.):
        avalue = metric_a(point)
        beta_value = np.sqrt(1-avalue)
        probes.append({"rho": float(point), "A": float(avalue),
                       "beta": float(beta_value),
                       "beta_prime": float(-derivative_a(point)/(2*beta_value)),
                       "speed_plus": float(1-beta_value),
                       "speed_minus": float(-1-beta_value)})
    assert abs(metric_a(horizon)) < 1e-13
    assert all(row["speed_minus"] < 0 for row in probes)
    assert all(row["speed_plus"] < 0 for row in probes[:2])
    assert probes[-1]["speed_plus"] > 0

    # Analytic global-growth bound: x=pi/2-atan(rho) lies in (0,pi), and
    # 2|rho| <= 1+rho². The remaining nonnegative polynomial is a square.
    v = sp.symbols("v", nonnegative=True)
    square_residual = sp.expand(1+v*v-2*v-(v-1)**2)
    assert square_residual == 0
    return {
        "schema": "nsc-dirac-tetrad-v1", "artifact_id": "NSC-4-DIRAC-TETRAD",
        "classification": "exact_fixed_background_lorentzian_connection_current_and_boundary_domain_derivation",
        "source_hashes": {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                          for name in ("scripts/check_nsc_dirac_tetrad.py", "docs/nsc-dirac-tetrad.md")},
        "conventions": {"signature": "+---", "units": "c=hbar=L_throat=1",
                        "dimension": "four spacetime dimensions; no compact conformal warp",
                        "metric": "d_tau²-(d_rho+beta*d_tau)²-r²*d_Omega²",
                        "geometry": "r=sqrt(1+rho²); A=1+3rho+3(1+rho²)(atan(rho)-pi/2); beta=sqrt(1-A)",
                        "spin_connection": "nabla=d+omega_ab gamma^a gamma^b/4 with lowered antisymmetric Lorentz indices",
                        "angular_basis": "spinor-sphere kappa=plus_or_minus(ell+1); canonical radial Pauli basis",
                        "measure_before": "r² d_rho d_Omega", "rescaling": "u=r*psi",
                        "measure_after": "d_rho per normalized angular mode"},
        "cartan": {"metric_residual": metric_residual, "torsion_residuals": torsion,
                   "lowered_connection_forms": {"01": "-beta_prime e1", "02": "-beta r_prime/r e2",
                                                "03": "-beta r_prime/r e3", "12": "r_prime/r e2",
                                                "13": "r_prime/r e3", "23": "cot(theta)/r e3"},
                   "metric_compatibility_exact": True, "gamma_clifford_relations_exact": True},
        "operator": {"contracted_connection": "-gamma0(beta_prime/2+beta*r_prime/r)+gamma1*r_prime/r+gamma2*cot(theta)/(2r)",
                     "connection_matrix_residual": connection_residual,
                     "radial_rescaling_matrix_residual": rescaling_residual,
                     "hamiltonian_principal_matrix_residual": radial_principal_residual,
                     "hamiltonian_shift_matrix_residual": shift_connection_residual,
                     "radial_hamiltonian": "H=-i(sigma2-beta I)partial_rho+i beta_prime I/2+(kappa/r)sigma1",
                     "symmetrized_form": "H=-i(C partial_rho+C_prime/2)+V; C=sigma2-beta I; V=(kappa/r)sigma1"},
        "current": {"density": "u_dagger u", "flux": "u_dagger(sigma2-beta I)u=2 Im(conj(u1)u2)-beta u_dagger u",
                    "continuity_residual": str(continuity),
                    "green_form": "<u,Hv>-<Hu,v>=-i[u_dagger C v]_left^right",
                    "missing_spin_connection_source": "-beta_prime*u_dagger*u",
                    "negative_control_polynomial": str(false_source),
                    "negative_control_is_nonzero": True,
                    "old_u1_zero_wall_unit_trace_current": str(old_wall_current)},
        "boundary": {"throat_beta": str(throat_beta),
                     "throat_current_matrix_determinant": str(throat_c.det()),
                     "throat_current_inertia": {"positive": 0, "negative": 2, "zero": 0},
                     "joined_transmission": "u_child=u_parent; generalized U obeys U_dagger C(0)U=C(0)",
                     "interval_assumptions": "left endpoint trapped beta_left>1; right endpoint exterior 0<beta_right<1",
                     "endpoint_form_eigenvalues": ["1-beta_right", "-1-beta_right", "beta_left-1", "beta_left+1"],
                     "exact_inertia_representative": {"beta_left": str(throat_beta), "beta_right": "1/2",
                                                       "eigenvalues": [str(e) for e in exact_boundary_eigenvalues],
                                                       "positive": positive, "negative": negative},
                     "forward_incoming_characteristics": {"left": 0, "right": 1},
                     "forward_outgoing_characteristics": {"left": 2, "right": 1},
                     "nonzero_local_reflecting_line_at_trapped_cut": False,
                     "endpoint_only_self_adjoint_finite_interval_realization": False,
                     "horizon": float(horizon), "actual_geometry_probes": probes},
        "complete_line": {"domain": "initial core C_c^infinity(R;C²), each fixed finite kappa",
                          "growth_bound": "beta² <= (3pi+3/2)(1+rho²); hence |1-beta|,|-1-beta| <= C(1+|rho|)",
                          "bound_square_residual": str(square_residual),
                          "flow_argument": "smooth vector fields a_plus=1-beta and a_minus=-1-beta have complete flows by linear growth; zero at horizon is stationary and adds no finite-time boundary",
                          "transport_group": "U(t)f(rho)=sqrt(partial_rho Phi(-t,rho))*f(Phi(-t,rho))",
                          "generator": "-i(a partial_rho+a_prime/2)",
                          "bounded_potential": "||V|| <= |kappa| because r>=1",
                          "self_adjointness_argument": "complete smooth half-density transport generators on C_c^infinity, followed by bounded Hermitian perturbation",
                          "causal_completeness_argument": "|d_rho/d_tau+beta|<=1 prevents finite-tau escape; S² is compact; the declared boundaryless product manifold has Cauchy tau slices",
                          "claim_scope": "fixed unwarped background fermion evolution; no finite reflecting truncation, vacuum, metric variation or gravitational stability theorem"},
        "primary_sources": [
            {"url": "https://arxiv.org/html/gr-qc/0605031", "use": "Schwarzschild PG connection/shift control, section IV.2; not a proof for this geometry"},
            {"url": "https://www.sciencedirect.com/science/article/pii/0022123673900037", "use": "Chernoff hyperbolic-generator essential self-adjointness framework; explicit complete-flow argument supplied here"},
            {"url": "https://arxiv.org/abs/1512.00761", "use": "horizon nonellipticity can coexist with Dirac self-adjointness; its particular boundary theorem is not imported"}],
        "nonclaims": {"compact_warp_derived": False, "physical_vacuum_selected": False,
                      "full_covariant_action_derived": False, "stationary_geometry_derived": False,
                      "metric_hessian_health_verified": False, "particle_or_cosmological_scale_predicted": False,
                      "new_to_world_priority_established": False},
        "terminal": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = calculate()
    if args.check:
        compare(json.loads(OUTPUT.read_text()), result)
        print("PG tetrad, Lorentzian current and boundary-domain result reproduced; every field checked.")
    elif args.output:
        with args.output.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
        print(args.output)
    else:
        print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
