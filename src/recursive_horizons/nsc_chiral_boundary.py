"""Domain-aware paired-angular boundary response of the spatial throat Dirac.

The implicit-midpoint transfer regulator carries all four components at every
trace. It is deliberately distinct from NSC-3's staggered node/edge regulator.
The PG current/domain is derived here; this finite spatial boundary solve is
not a reflecting approximation to the Lorentzian trapped-region problem.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import splu

from .nsc_boundary import compare, jsonable
from .nsc_lorentzian import geometry, horizon

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = "NSC-8-CHIRAL-BOUNDARY-v1"
DEFAULT_RESULT = ROOT / "results/nsc-8-chiral-boundary.json"
SOURCE_PATHS = (
    "src/recursive_horizons/nsc_chiral_boundary.py",
    "scripts/check_nsc_chiral_boundary.py",
    "tests/test_nsc_chiral_boundary.py",
    "docs/nsc-chiral-boundary.md",
    "src/recursive_horizons/nsc_boundary.py",
    "src/recursive_horizons/nsc_lorentzian.py",
    "src/recursive_horizons/nsc_spinor_bridge.py",
)
ENERGIES = (.2+.25j, .6+.35j, 1.1+.3j)
RADIUS = 4.
PORTS = (1., 0., -1.)  # parent, eliminated throat, child
I2 = np.eye(2, dtype=complex)
I4 = np.eye(4, dtype=complex)
PAULI = (I2, np.array([[0, 1], [1, 0]], complex),
         np.array([[0, -1j], [1j, 0]], complex), np.diag([1., -1.]))
ALPHA = np.kron(I2, PAULI[2])
ANGULAR = np.kron(PAULI[3], PAULI[1])
BETA = np.kron(I2, PAULI[3])
GAMMA5 = -np.kron(PAULI[1], PAULI[2])
WALL_Q = np.kron(PAULI[3], PAULI[3])
REFLECTION = np.kron(PAULI[1], PAULI[3])


def norm(a):
    return float(np.linalg.norm(a))


def relative(a, reference):
    return norm(a-reference)/max(norm(reference), 1e-30)


def comm(a, b):
    return a@b-b@a


def wall_rows(kind="compatible"):
    """Rows annihilating the allowed Q=+1 wall data (two dimensions)."""
    if kind not in ("compatible", "legacy_radial"):
        raise ValueError("unknown spatial wall")
    q = WALL_Q if kind == "compatible" else -BETA
    return I4[np.real(np.diag(q)) < 0]


def potential(x, mass=0., kappa=1.):
    return kappa/np.sqrt(1+x*x)*ANGULAR + mass*BETA


def ode_matrix(x, energy, mass=0., kappa=1.):
    # (-i alpha d + V - z) psi=0; jump of G is +i alpha.
    return 1j*ALPHA@(energy*I4-potential(x, mass, kappa))


def midpoint_step(left, right, energy, mass=0., kappa=1.):
    half = .5*(right-left)*ode_matrix(.5*(left+right), energy, mass, kappa)
    return np.linalg.solve(I4-half, I4+half)


def transfer_data(intervals=160, energy=ENERGIES[0], mass=0., kappa=1.,
                  radius=RADIUS):
    if not isinstance(intervals, int) or intervals < 8 or intervals % 8:
        raise ValueError("intervals must be a positive multiple of eight, at least eight")
    grid = np.linspace(-radius, radius, intervals+1)
    values = [I4.copy()]
    for left, right in zip(grid[:-1], grid[1:]):
        values.append(midpoint_step(left, right, energy, mass, kappa)@values[-1])
    return grid, values


def continuum_data(energy, mass=0., kappa=1., radius=RADIUS):
    solution = solve_ivp(
        lambda x, y: (ode_matrix(x, energy, mass, kappa)@y.reshape(4, 4)).ravel(),
        (-radius, radius), I4.ravel(), method="DOP853", rtol=2e-12, atol=2e-14,
        dense_output=True)
    if not solution.success:
        raise RuntimeError(solution.message)
    return lambda x: solution.sol(x).reshape(4, 4)


def green_from_fundamental(fundamental, points, radius=RADIUS, wall="compatible"):
    """Symmetric coincident trace of (H-z)^-1, in one coordinate spin frame.

    Point sources have unit integrated weight. Consequently G's entries are
    dimensionless in c=hbar=L_star=1. Its inverse is a port kernel, not a local
    Hamiltonian mass. Noncoincident entries are ordinary continuum Green data.
    """
    w = wall_rows(wall)
    end = fundamental(radius)
    boundary = np.vstack((w, w@end))
    columns = []
    for y in points:
        jump = np.linalg.solve(fundamental(y), 1j*ALPHA)
        initial = np.linalg.solve(boundary, np.vstack((np.zeros((2, 4)), -w@end@jump)))
        columns.append([fundamental(x)@(initial+(1. if x > y else .5 if x == y else 0.)*jump)
                        for x in points])
    return np.block([[columns[j][i] for j in range(len(points))]
                     for i in range(len(points))])


def transfer_green(intervals=160, energy=ENERGIES[0], mass=0., wall="compatible",
                   points=PORTS, radius=RADIUS, kappa=1.):
    grid, values = transfer_data(intervals, energy, mass, kappa, radius)

    def fundamental(x):
        node = int(np.argmin(abs(grid-x)))
        if abs(grid[node]-x) > 1e-12:
            raise ValueError("every retained port must be a grid node")
        return values[node]

    return green_from_fundamental(fundamental, points, radius, wall)


def joined_green(intervals=160, energy=ENERGIES[0], mass=0., wall="compatible",
                 points=PORTS, radius=RADIUS, kappa=1.):
    """Independent sparse solve of every cell equation and every source jump.

    Each source node has distinct left/right traces. Interior cells impose
    implicit-midpoint equations; sources impose psi+ - psi- = i alpha f.
    There is no forward transfer product or dense inverse of a bulk matrix.
    """
    grid = np.linspace(-radius, radius, intervals+1)
    nodes = [int(np.argmin(abs(grid-x))) for x in points]
    if len(set(nodes)) != len(nodes) or any(abs(grid[n]-x) > 1e-12 or n in (0, intervals)
                                           for n, x in zip(nodes, points)):
        raise ValueError("distinct interior ports must lie on the grid")
    left_ids, right_ids = [], []
    count = 0
    for node in range(intervals+1):
        left_ids.append(count)
        count += 1
        if node in nodes:
            right_ids.append(count)
            count += 1
        else:
            right_ids.append(left_ids[-1])
    matrix = lil_matrix((4*count, 4*count), dtype=complex)
    rhs = np.zeros((4*count, 4*len(points)), complex)
    row = 0
    for node, (left, right) in enumerate(zip(grid[:-1], grid[1:])):
        half = .5*(right-left)*ode_matrix(.5*(left+right), energy, mass, kappa)
        matrix[row:row+4, 4*left_ids[node+1]:4*left_ids[node+1]+4] = I4-half
        matrix[row:row+4, 4*right_ids[node]:4*right_ids[node]+4] = -I4-half
        row += 4
    for p, node in enumerate(nodes):
        matrix[row:row+4, 4*right_ids[node]:4*right_ids[node]+4] = I4
        matrix[row:row+4, 4*left_ids[node]:4*left_ids[node]+4] = -I4
        rhs[row:row+4, 4*p:4*p+4] = 1j*ALPHA
        row += 4
    wall_data = wall_rows(wall)
    matrix[row:row+2, :4] = wall_data
    matrix[row+2:row+4, -4:] = wall_data
    matrix = matrix.tocsc()
    solution = splu(matrix).solve(rhs)
    green = np.vstack([.5*(solution[4*left_ids[n]:4*left_ids[n]+4]
                           +solution[4*right_ids[n]:4*right_ids[n]+4]) for n in nodes])
    residual = norm(matrix@solution-rhs)/norm(rhs)
    return green, residual


def eliminate_center(green):
    kernel = np.linalg.solve(green, np.eye(12))
    retained = [0, 1, 2, 3, 8, 9, 10, 11]
    removed = [4, 5, 6, 7]
    krr = kernel[np.ix_(retained, retained)]
    kre = kernel[np.ix_(retained, removed)]
    ker = kernel[np.ix_(removed, retained)]
    kee = kernel[np.ix_(removed, removed)]
    return krr-kre@np.linalg.solve(kee, ker)


def clifford_basis():
    """An algebraic +--- Clifford completion in the paired representation.

    alpha_1=I rho2, alpha_2=eta3 rho1, alpha_3=eta2 rho1; their
    orientation gives i gamma0 gamma1 gamma2 gamma3 = physical Gamma5.
    The remaining two axes are a matrix completion, not a separately derived
    angular tetrad projection. A spherical partial-wave port kernel remains
    nonlocal; these coefficients classify Clifford content, not a local action.
    """
    alphas = (ALPHA, ANGULAR, np.kron(PAULI[2], PAULI[1]))
    gamma = (BETA, *(BETA@a for a in alphas))
    basis = {"scalar": I4, "pseudoscalar": GAMMA5}
    basis.update({f"vector_{a}": g for a, g in enumerate(gamma)})
    basis.update({f"axial_{a}": g@GAMMA5 for a, g in enumerate(gamma)})
    basis.update({f"tensor_{a}{b}": .5j*comm(gamma[a], gamma[b])
                  for a in range(4) for b in range(a+1, 4)})
    return basis


def decompose_link(link):
    action = BETA@link
    basis = clifford_basis()
    coefficients = {key: np.trace(matrix.conj().T@action)/4 for key, matrix in basis.items()}
    rebuilt = sum(coefficients[key]*matrix for key, matrix in basis.items())
    forbidden = [key for key in basis if key in ("scalar", "pseudoscalar") or key.startswith("tensor")]
    even = .5*(link+GAMMA5@link@GAMMA5)
    odd = .5*(link-GAMMA5@link@GAMMA5)
    return {"hamiltonian_port_link": link, "action_port_kernel_beta_B": action,
            "action_clifford_coefficients": coefficients,
            "clifford_reconstruction_relative_error": relative(rebuilt, action),
            "link_norm": norm(link), "chiral_even_link_norm": norm(even),
            "chiral_odd_link_norm": norm(odd),
            "commutator_relative_norm": norm(comm(link, GAMMA5))/norm(link),
            "anticommutator_relative_norm": norm(link@GAMMA5+GAMMA5@link)/norm(link),
            "action_scalar_pseudoscalar_tensor_coefficient_norm": float(np.linalg.norm([coefficients[k] for k in forbidden]))}


def algebra_and_domain():
    spatial_allowed = I4[:, np.real(np.diag(WALL_Q)) > 0]
    legacy_allowed = I4[:, np.real(np.diag(BETA)) < 0]
    endpoint_rows = []
    for x, normal, label in ((-4., -1, "child_outer_cut"), (0., 1, "child_throat"),
                             (0., -1, "parent_throat"), (4., 1, "parent_outer_cut")):
        shift = float(geometry(x)[0])
        current = ALPHA-shift*I4
        speeds = [1-shift, -1-shift]
        endpoint_rows.append({"location": label, "rho": x, "normal": normal,
                              "shift": shift, "coordinate_speeds_per_angular_pair": speeds,
                              "incoming_full_spinor_channels": 2*sum(normal*v < 0 for v in speeds),
                              "outgoing_full_spinor_channels": 2*sum(normal*v > 0 for v in speeds),
                              "current_commutator_gamma5_norm": norm(comm(current, GAMMA5))})
    beta0 = float(geometry(0.)[0])
    gamma = (BETA, BETA@ALPHA, BETA@ANGULAR, BETA@np.kron(PAULI[2], PAULI[1]))
    return {
        "physical_gamma5": GAMMA5, "beta_mass_matrix": BETA,
        "radial_alpha": ALPHA, "angular_matrix": ANGULAR,
        "gamma5_from_Clifford_orientation_residual": norm(1j*gamma[0]@gamma[1]@gamma[2]@gamma[3]-GAMMA5),
        "gamma5_square_residual": norm(GAMMA5@GAMMA5-I4),
        "kinetic_commutator_residual": norm(comm(ALPHA, GAMMA5)),
        "angular_commutator_residual": norm(comm(ANGULAR, GAMMA5)),
        "mass_anticommutator_residual": norm(BETA@GAMMA5+GAMMA5@BETA),
        "compatible_wall_commutator_residual": norm(comm(WALL_Q, GAMMA5)),
        "compatible_wall_spatial_flux_residual": norm(spatial_allowed.conj().T@ALPHA@spatial_allowed),
        "legacy_wall_chirality_leakage_norm": norm(wall_rows("legacy_radial")@GAMMA5@legacy_allowed),
        "same_wall_PG_throat_flux_norm": norm(spatial_allowed.conj().T@(ALPHA-beta0*I4)@spatial_allowed),
        "PG_horizon_rho": horizon(), "oriented_PG_characteristics": endpoint_rows,
        "identity_transmission_oriented_flux_residual": norm(-(ALPHA-beta0*I4)+(ALPHA-beta0*I4)),
        "PG_generator": "-i I_eta tensor (rho2-beta_shift I) d_rho + (i/2) beta_shift_prime I4 + (kappa/r) eta3 tensor rho1",
        "PG_current": "psi_dagger [I_eta tensor (rho2-beta_shift I)] psi",
        "transmission_domain": "common coordinate spin frame: psi_parent(0)=psi_child(0); parent normal=-d_rho, child normal=+d_rho",
        "PG_outer_domain": "child cut all four channels outflow; parent cut two incoming channels fixed to zero; throat is transmission, not a reflecting wall",
        "schur_theorem": "If [K,Gamma5]=0 and retained/eliminated projections commute with Gamma5, invertible K_ee commutes with Gamma5 and S=K_rr-K_re K_ee^-1 K_er commutes with Gamma5_r. Nonzero common-frame offdiagonal B then cannot also anticommute with Gamma5.",
        "kernel_scope": "spatial finite interval on fixed r=sqrt(1+rho^2); no PG stationary map across horizon is evaluated",
    }


def frame_control(link):
    child_gamma = REFLECTION@GAMMA5@REFLECTION.conj().T
    local_link = link@REFLECTION.conj().T
    common_restored = local_link@REFLECTION
    return {"child_frame_map": REFLECTION, "child_frame_gamma5": child_gamma,
            "child_frame_wall_Q": REFLECTION@WALL_Q@REFLECTION.conj().T,
            "child_gamma_is_minus_common_residual": norm(child_gamma+GAMMA5),
            "misleading_raw_same_gamma_anticommutator_relative_norm": norm(local_link@GAMMA5+GAMMA5@local_link)/norm(link),
            "correct_intertwining_relative_norm": norm(GAMMA5@local_link-local_link@child_gamma)/norm(link),
            "restored_common_frame_relative_error": relative(common_restored, link),
            "common_frame_anticommutator_relative_norm": norm(link@GAMMA5+GAMMA5@link)/norm(link),
            "passive_frame_transport_only": True,
            "same_Q_plus_outer_walls_are_reflection_invariant": False,
            "reflection_alone_derives_scalar_mass": False}


def response_row(energy, intervals=160, mass=0., wall="compatible"):
    direct, residual = joined_green(intervals, energy, mass, wall)
    transfer = transfer_green(intervals, energy, mass, wall)
    continuum = green_from_fundamental(continuum_data(energy, mass), PORTS, wall=wall)
    response = eliminate_center(direct)
    port_green, _ = joined_green(intervals, energy, mass, wall, points=(1., -1.))
    reduced_direct = np.linalg.solve(port_green, np.eye(8))
    continuum_response = eliminate_center(continuum)
    chirality = np.kron(np.eye(2), GAMMA5)
    pi_sheet = .5*(np.eye(8)-np.kron(PAULI[3], GAMMA5))
    return {"energy": energy, "intervals": intervals, "mass_control": mass, "wall": wall,
            "three_port_green": direct, "two_port_response": response,
            "three_port_green_condition_number": float(np.linalg.cond(direct)),
            "sparse_linear_relative_residual": residual,
            "direct_versus_transfer_relative_error": relative(direct, transfer),
            "direct_versus_Schur_relative_error": relative(response, reduced_direct),
            "continuum_green_relative_error": relative(direct, continuum),
            "continuum_response_relative_error": relative(response, continuum_response),
            "continuum_response_chiral_commutator_relative_norm": norm(comm(continuum_response, chirality))/norm(continuum_response),
            "continuum_link_chiral_odd_norm": norm(.5*(continuum_response[:4, 4:]-GAMMA5@continuum_response[:4, 4:]@GAMMA5)),
            "response_chiral_commutator_relative_norm": norm(comm(response, chirality))/norm(response),
            "Pi_sheet_commutator_relative_norm": norm(comm(response, pi_sheet))/norm(response),
            "link": decompose_link(response[:4, 4:])}


def mass_control(phi=.35):
    kinetic = .7*ALPHA+.9*ANGULAR
    link = phi*BETA
    full = np.block([[kinetic, link], [link, kinetic]])
    pi_sheet = .5*(np.eye(8)-np.kron(PAULI[3], GAMMA5))
    return {"inserted_phi": phi, "inserted_link": link,
            "Pi_sheet_commutator_norm": norm(comm(full, pi_sheet)),
            "dispersion_square_residual": norm(full@full-(.7**2+.9**2+phi**2)*np.eye(8)),
            "action_scalar_coefficient": np.trace(BETA@link)/4,
            "derived_from_throat": False}


def convergence():
    energy = ENERGIES[0]
    continuum = green_from_fundamental(continuum_data(energy), PORTS)
    target = eliminate_center(continuum)
    rows = []
    for intervals in (40, 80, 160, 320):
        green, _ = joined_green(intervals, energy)
        response = eliminate_center(green)
        rows.append({"intervals": intervals, "step": 2*RADIUS/intervals,
                     "green_relative_error": relative(green, continuum),
                     "response_relative_error": relative(response, target),
                     "link_relative_error": relative(response[:4, 4:], target[:4, 4:])})
    return {"energy": energy, "rows": rows,
            "response_orders": [float(np.log2(a["response_relative_error"]/b["response_relative_error"]))
                                for a, b in zip(rows[:-1], rows[1:])],
            "continuum_two_port_response": target,
            "continuum_fullspinor_link": decompose_link(target[:4, 4:])}


def regulator_checks():
    energy = .6
    step = midpoint_step(-.1, .1, energy)
    complex_step = midpoint_step(-.1, .1, ENERGIES[0])
    # Discrete Green identity: T† alpha T-alpha = -2 Im(z) h M†M.
    average = .5*(I4+complex_step)
    current_change = complex_step.conj().T@ALPHA@complex_step-ALPHA
    expected = -2*ENERGIES[0].imag*.2*(average.conj().T@average)
    green_plus, _ = joined_green(80, ENERGIES[0])
    green_minus, _ = joined_green(80, ENERGIES[0].conjugate())
    imaginary_green = (green_plus-green_plus.conj().T)/(2j)
    return {"real_energy_cell_current_residual": norm(step.conj().T@ALPHA@step-ALPHA),
            "complex_energy_cell_Green_identity_residual": norm(current_change-expected),
            "cell_chiral_commutator_residual": norm(comm(step, GAMMA5)),
            "resolvent_Schwarz_reflection_relative_error": relative(green_plus.conj().T, green_minus),
            "upper_half_plane_green_imaginary_min_eigenvalue": float(np.min(np.linalg.eigvalsh(imaginary_green))),
            "all_components_share_trace_space": True,
            "legacy_staggered_gamma5_assumed": False,
            "absence_of_spectral_doubling_claimed": False}


def build_record():
    rows = [response_row(energy) for energy in ENERGIES]
    refinement = convergence()
    wall_control = response_row(ENERGIES[0], wall="legacy_radial")
    inserted_mass = response_row(ENERGIES[0], mass=.35)
    identities = algebra_and_domain()
    regulator = regulator_checks()
    frame = frame_control(rows[0]["link"]["hamiltonian_port_link"])
    gate = {
        "full_paired_angular_chirality_preserved": max(r["response_chiral_commutator_relative_norm"] for r in rows) < 1e-10,
        "direct_joined_equals_eliminated": max(r["direct_versus_Schur_relative_error"] for r in rows) < 1e-10,
        "independent_transfer_matches_joined": max(r["direct_versus_transfer_relative_error"] for r in rows) < 1e-10,
        "continuum_response_agrees_below_one_percent": max(r["continuum_response_relative_error"] for r in rows) < .01,
        "resolution_second_order": all(1.8 < p < 2.2 for p in refinement["response_orders"]),
        "nonzero_baseline_link": min(r["link"]["link_norm"] for r in rows) > .1,
        "baseline_scalar_pseudoscalar_tensor_absent": max(r["link"]["action_scalar_pseudoscalar_tensor_coefficient_norm"] for r in rows) < 1e-10,
        "mass_control_has_chiral_odd_response": min(inserted_mass["link"]["chiral_odd_link_norm"], inserted_mass["continuum_link_chiral_odd_norm"]) > .01,
        "legacy_wall_causes_chiral_breaking": min(wall_control["response_chiral_commutator_relative_norm"], wall_control["continuum_response_chiral_commutator_relative_norm"]) > .01,
        "frame_transport_restores_obstruction": frame["correct_intertwining_relative_norm"] < 1e-10 and frame["common_frame_anticommutator_relative_norm"] > 1.9,
        "cell_current_and_Green_identity": max(regulator["real_energy_cell_current_residual"], regulator["complex_energy_cell_Green_identity_residual"]) < 1e-12,
        "resolvent_has_correct_Schwarz_and_Herglotz_sign": regulator["resolvent_Schwarz_reflection_relative_error"] < 1e-12 and regulator["upper_half_plane_green_imaginary_min_eigenvalue"] > 0,
        "massless_response_does_not_preserve_candidate_Pi_sheet": min(r["Pi_sheet_commutator_relative_norm"] for r in rows) > .01,
    }
    if not all(gate.values()):
        raise RuntimeError(f"chiral boundary gate failed: {gate}")
    return {"schema": SCHEMA,
            "classification": "evaluated_finite_spatial_fullspinor_response_and_domain_aware_chiral_obstruction",
            "source_hashes": [{"path": p, "sha256": hashlib.sha256((ROOT/p).read_bytes()).hexdigest()} for p in SOURCE_PATHS],
            "conventions": {"tensor_order": "paired angular eta, radial rho; ports parent, throat, child",
                            "units": "hbar=c=L_star=1; point sources carry unit integrated radial weight",
                            "resolvent": "(H-z)^-1; trace diagonal is arithmetic mean of left/right limits",
                            "radius": RADIUS, "kappa_pair": [1., -1.], "ports": PORTS,
                            "spatial_domain": "separated Q=eta3 tensor rho3 eigenvalue +1 at both exterior endpoints; throat continuity",
                            "regulator": "implicit midpoint first-order cells with exact source jumps; four spinor traces at every node",
                            "action_decomposition": "beta times Hamiltonian-form port link, coefficients Tr(C_dagger beta B)/4 in 16 orthogonal Clifford matrices",
                            "Clifford_completion_scope": "alpha2=angular_matrix and alpha3=eta2 tensor rho1 complete the algebra; no independent local angular tetrad projection is inferred",
                            "pole_scope": "three nonreal energies; no pole search, local mass fit, or asymptotic spectral-gap inference"},
            "algebra_and_domain": identities, "regulator_checks": regulator,
            "energy_resolved_maps": rows, "resolution": refinement,
            "frame_transport_control": frame,
            "inserted_scalar_mass_control": {"local_link_algebra": mass_control(), "boundary_response": inserted_mass},
            "artificial_domain_breaking_wall_control": wall_control,
            "gate": gate,
            "nonclaims": {"smooth_massless_throat_generates_scalar_Dirac_mass": False,
                          "physical_PG_boundary_map_evaluated": False,
                          "reflection_selects_Pi_sheet_sector": False,
                          "local_4D_covariant_action_derived_from_port_kernel": False,
                          "stationary_common_action_solved": False,
                          "electron_mass_or_observational_target_fitted": False}}
