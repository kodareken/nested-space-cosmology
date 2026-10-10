"""One sourced de Sitter cap joined through a vacuum bridge to an AE end.

Known spherical constraints and the v3 wall are reused, with fixed source
parameters. This is initial-slice compatibility, not a formation calculation.
"""

import math

import numpy as np

from .nsc_atlas_moving_boundaries import VacuumWallParameters, vacuum_wall_state


SCHEMA = "nsc-atlas-connection-slice-v1"


def smooth_step(x):
    """C-infinity step and its first derivative, including flat endpoints."""
    x = np.asarray(x, dtype=float)
    if not np.all(np.isfinite(x)):
        raise ValueError("step arguments must be finite")
    value = np.where(x >= 1, 1.0, 0.0)
    derivative = np.zeros_like(x)
    mask = (x > 0) & (x < 1)
    z = x[mask]
    logit = -1/z + 1/(1-z)
    exponential = np.exp(-np.abs(logit))
    s = np.where(logit >= 0, 1/(1+exponential), exponential/(1+exponential))
    value[mask] = s
    derivative[mask] = exponential/(1+exponential)**2 * (1/z**2+1/(1-z)**2)
    return value, derivative


def wall_slice_data(length_scale=1.0):
    """The v3 source tuple at Rw/ell=1.5, boosted to a flat de Sitter slice."""
    ell = float(length_scale)
    if not math.isfinite(ell) or ell <= 0:
        raise ValueError("length_scale must be finite and positive")
    parameters = VacuumWallParameters(H=1/ell, k=2/ell, mass=.4*ell)
    R = 1.5*ell
    state = vacuum_wall_state(R, parameters)
    speed, bm, bp = (float(state[key]) for key in ("Rdot", "beta_minus", "beta_plus"))
    Vminus = parameters.H*R
    fminus = 1-Vminus**2
    gamma = (bm-Vminus*speed)/fminus
    gamma_v = (speed-Vminus*bm)/fminus
    v = gamma_v/gamma
    Rprime_plus = gamma*bp-gamma_v*speed
    Vplus = gamma*speed-gamma_v*bp
    if gamma < 1-1e-12 or abs(v) >= 1 or Rprime_plus >= 0:
        raise ValueError("wall does not provide the declared outward slice orientation")
    return {
        "ell": ell, "parameters": parameters, "R": R,
        "Rdot": speed, "beta_minus": bm, "beta_plus": bp,
        "gamma": gamma, "gamma_v": gamma_v, "v": v,
        "Rprime_minus": 1.0, "Rprime_plus": Rprime_plus,
        "V_minus": Vminus, "V_plus": Vplus,
        "energy_surface": parameters.sigma*gamma,
        "momentum_surface": parameters.sigma*gamma_v,
        "taper_stop": ell,
    }


def vacuum_transition(radius, wall):
    """Decreasing-radius spacelike Schwarzschild data for fixed slicing V(R)."""
    R = np.asarray(radius, dtype=float)
    ell, Rw, R0 = wall["ell"], wall["R"], wall["taper_stop"]
    if not np.all(np.isfinite(R)) or np.any(R < R0) or np.any(R > Rw):
        raise ValueError("transition is restricted to its recorded [R0,Rw] interval")
    M = wall["parameters"].mass
    if not 2*M < R0 < Rw:
        raise ValueError("taper must stop strictly outside the throat")
    step, step_prime = smooth_step((R-R0)/(Rw-R0))
    V = wall["V_plus"]*step
    V_R = wall["V_plus"]*step_prime/(Rw-R0)
    f = 1-2*M/R
    Rprime = -np.sqrt(f+V*V)
    Rsecond = M/R**2+V*V_R
    A, B = -V_R, -V/R
    Bprime = (-V_R/R+V/R**2)*Rprime
    scalar = -4*Rsecond/R+2*(1-Rprime**2)/R**2
    hamiltonian = scalar+4*A*B+2*B*B
    momentum = 2*Rprime/R*(A-B)-2*Bprime
    # Kuchar's embedding formula T_l=V/f; no coordinate age is inferred.
    T_l = V/f
    metric_norm = -f*T_l**2+Rprime**2/f
    return {
        "R": R, "Rprime": Rprime, "Rsecond": Rsecond, "V": V, "V_R": V_R,
        "Kll": A, "Ktheta": B, "scalar_curvature": scalar,
        "hamiltonian_residual": hamiltonian, "momentum_residual": momentum,
        "mass": R/2*(1+V*V-Rprime**2),
        "embedding_spatial_norm": metric_norm,
        "embedding_future_normal_norm": -f*(Rprime/f)**2+V*V/f,
    }


def transition_proper_distance(radii, wall, order=32):
    """Proper distance increasing outward from the cap centre."""
    radii = np.asarray(radii, dtype=float)
    if radii.ndim != 1 or len(radii) < 2 or np.any(np.diff(radii) >= 0):
        raise ValueError("transition radii must be a strictly decreasing vector")
    if not math.isclose(float(radii[0]), wall["R"], rel_tol=1e-13):
        raise ValueError("distance grid must start at the wall")
    if not isinstance(order, int) or not 4 <= order <= 64:
        raise ValueError("quadrature order must be an integer from 4 to 64")
    vacuum_transition(radii, wall)
    nodes, weights = np.polynomial.legendre.leggauss(order)
    distances = [wall["R"]]
    for upper, lower in zip(radii[:-1], radii[1:]):
        R = (upper+lower)/2+(upper-lower)/2*nodes
        radial = vacuum_transition(R, wall)["Rprime"]
        distances.append(distances[-1]+float((upper-lower)/2*np.dot(weights, -1/radial)))
    return np.asarray(distances)


def schwarzschild_bridge(u, mass, throat_distance=0.0):
    """Regular proper-distance parametrization, including the minimum sphere."""
    u = np.asarray(u, dtype=float)
    if not math.isfinite(mass) or mass <= 0 or not np.all(np.isfinite(u)):
        raise ValueError("bridge mass must be positive and parameters finite")
    if np.any(np.abs(u) > 10) or not math.isfinite(throat_distance):
        raise ValueError("bridge sample is outside the finite display adapter")
    R = 2*mass*np.cosh(u)**2
    l = throat_distance+2*mass*(u+np.sinh(u)*np.cosh(u))
    Rprime = np.tanh(u)
    Rsecond = mass/R**2
    scalar = -4*Rsecond/R+2*(1-Rprime**2)/R**2
    return {"u": u, "R": R, "l": l, "Rprime": Rprime,
            "Rsecond": Rsecond, "scalar_curvature": scalar,
            "mass": R/2*(1-Rprime**2)}


def _maxabs(array):
    return float(np.max(np.abs(array)))


def _finite_difference_checks(wall):
    ell = wall["ell"]
    R = ell*np.linspace(1.015, 1.485, 31)
    fields = vacuum_transition(R, wall)

    def derivative(key, h):
        def field(x):
            return vacuum_transition(x, wall)[key]
        return (field(R-2*h)-8*field(R-h)+8*field(R+h)-field(R+2*h))/(12*h)

    results = []
    for h in (1e-4*ell, 5e-5*ell):
        Rsecond = derivative("Rprime", h)*fields["Rprime"]
        Bprime = derivative("Ktheta", h)*fields["Rprime"]
        scalar = -4*Rsecond/R+2*(1-fields["Rprime"]**2)/R**2
        A, B = fields["Kll"], fields["Ktheta"]
        results.append({
            "step_over_ell": h/ell,
            "hamiltonian_max_ell2": ell**2*_maxabs(scalar+4*A*B+2*B*B),
            "momentum_max_ell2": ell**2*_maxabs(2*fields["Rprime"]/R*(A-B)-2*Bprime),
            "Rsecond_max_ell_error": ell*_maxabs(Rsecond-fields["Rsecond"]),
        })
    return results


def connection_slice_record(length_scale=1.0):
    """One finite atlas connection with an exact analytic AE continuation."""
    wall = wall_slice_data(length_scale)
    ell, parameters = wall["ell"], wall["parameters"]
    Rw, R0, M = wall["R"], wall["taper_stop"], parameters.mass
    radii = np.linspace(Rw, R0, 65)
    fields = vacuum_transition(radii, wall)
    distance32 = transition_proper_distance(radii, wall, 32)
    distance64 = transition_proper_distance(radii, wall, 64)
    u0 = -math.acosh(math.sqrt(R0/(2*M)))
    throat_distance = distance32[-1]-2*M*(u0+math.sinh(u0)*math.cosh(u0))
    uend = math.acosh(math.sqrt(6*ell/(2*M)))
    u = np.concatenate((np.linspace(u0, 0, 33), np.linspace(0, uend, 65)[1:]))
    bridge = schwarzschild_bridge(u, M, throat_distance)
    Rcap = np.linspace(0, Rw, 33)
    jump_Rprime = wall["Rprime_plus"]-1
    jump_V = wall["V_plus"]-wall["V_minus"]
    jump_B = -jump_V/Rw
    energy, momentum = wall["energy_surface"], wall["momentum_surface"]
    # Direct 1+1 tensor projection of the wall energy onto (T,s), including
    # delta(eta)|slice=delta(l-lw)/gamma. Angular tension has zero T,s projection.
    Ucontr = np.array([wall["gamma"], wall["gamma_v"]])
    Ulower = np.diag([-1., 1.])@Ucontr
    slice_tensor = parameters.sigma*np.outer(Ulower, Ulower)/wall["gamma"]
    transition = {"l_over_ell": (distance32/ell).tolist(), "R_over_ell": (radii/ell).tolist()}
    for key in ("Rprime", "V"):
        transition[key] = fields[key].tolist()
    for key in ("Kll", "Ktheta", "Rsecond"):
        transition["ell_"+key] = (ell*fields[key]).tolist()
    transition["mass_over_ell"] = (fields["mass"]/ell).tolist()
    transition["ell2_scalar_curvature"] = (ell**2*fields["scalar_curvature"]).tolist()
    throat_index = int(np.argmin(bridge["R"]))
    return {
        "schema": SCHEMA, "status": "CONNECTED_SOURCED_INITIAL_SLICE",
        "claim_type": "connection_of_imported_wall_and_spherical_initial_data",
        "admission": {"owner": "lab/docs/nsc-atlas-connection-slice.md",
                      "borrowed": ["v3 signed vacuum wall", "gr-qc/9403003v1 spherical constraints and Eq81-82",
                                   "gr-qc/9708008v2 shell frame and surface stress"],
                      "missing": "a sourced interior connected to an exact vacuum AE end",
                      "decision": "one cap/wall/vacuum-bridge/AE initial slice is explicit",
                      "unlocks": "gluing an unchanged sourced body through its vacuum outer end",
                      "cpu_seconds_max": 60, "quadrature_evaluations": 64*(32+64),
                      "stop": "failed spacelike, bulk-constraint or shell-jump check",
                      "priority_claim": False},
        "theory": {"H_ell": 1., "k_ell": 2., "M_over_ell": .4,
                   "sigma_ell": parameters.sigma*ell, "wall_law": "pressure=-sigma, fixed positive sigma",
                   "fixed_source_within_test": True, "K_convention": "K=-Lie_T(h)/2",
                   "momentum_convention": "J_i=-T_ab T^a e_i^b",
                   "taper_changes": "spacelike slicing only, not the vacuum geometry or source"},
        "wall": {"R_over_ell": Rw/ell, "l_over_ell": Rw/ell, "Rdot": wall["Rdot"],
                 "gamma": wall["gamma"], "v": wall["v"], "gamma_v": wall["gamma_v"],
                 "Rprime_minus": 1., "Rprime_plus": wall["Rprime_plus"],
                 "V_minus": wall["V_minus"], "V_plus": wall["V_plus"],
                 "beta_minus": wall["beta_minus"], "beta_plus": wall["beta_plus"],
                 "energy_surface_ell": energy*ell, "momentum_surface_ell": momentum*ell,
                 "jump_Rprime": jump_Rprime, "jump_V": jump_V,
                 "surface_delta_factor": "delta(eta)|slice=delta(l-lw)/gamma"},
        "cap": {"geometry": "flat expanding de Sitter slice, finite ball",
                "regular_centre": True, "l_over_ell": (Rcap/ell).tolist(),
                "R_over_ell": (Rcap/ell).tolist(), "V": (parameters.H*Rcap).tolist(),
                "Rprime": np.ones_like(Rcap).tolist(), "ell_Kll": (-np.ones_like(Rcap)).tolist(),
                "ell_Ktheta": (-np.ones_like(Rcap)).tolist(), "rho_ell2": 3/(8*math.pi),
                "hamiltonian_residual_ell2": 0., "momentum_residual_ell2": 0.,
                "contains_antitrapped_spheres": True,
                "antitrapped_condition": "HR>1 gives both future null expansions positive"},
        "transition": {"status": "VACUUM_SCHWARZSCHILD_SLICE", "R0_over_ell": R0/ell,
                       "profile": "C-infinity step with all endpoint derivatives zero",
                       "proper_length_over_ell": float((distance32[-1]-Rw)/ell),
                       "samples": transition},
        "bridge": {"geometry": "time-symmetric Schwarzschild bridge",
                   "throat_R_over_ell": 2*M/ell, "throat_l_over_ell": float(throat_distance/ell),
                   "throat_Rprime": float(bridge["Rprime"][throat_index]),
                   "throat_ell_Rsecond": float(ell*bridge["Rsecond"][throat_index]),
                   "regular_at_throat": True, "K_is_zero": True,
                   "samples": {"u": u.tolist(), "R_over_ell": (bridge["R"]/ell).tolist(),
                               "l_over_ell": (bridge["l"]/ell).tolist(),
                               "Rprime": bridge["Rprime"].tolist(),
                               "mass_over_ell": (bridge["mass"]/ell).tolist()}},
        "parent_end": {"geometry": "exact Schwarzschild AE end", "ADM_mass_over_ell": M/ell,
                       "sample_outer_R_over_ell": 6., "analytic_extension": "R to infinity, K=0, h_RR=(1-2M/R)^-1",
                       "finite_display_is_physical_cutoff": False,
                       "exact_vacuum_overlap_available": True},
        "matching": {"cap_to_wall": "induced sphere and common shell boost",
                     "wall_to_transition": "distributional Hamiltonian and momentum constraints",
                     "transition_to_bridge": "exact metric and K on a smooth vacuum overlap",
                     "bridge_to_parent_end": "same Schwarzschild initial data",
                     "shared_two_body_parent_established": False,
                     "causal_channel_established": False},
        "independent_checks": {
            "boost_unit_max_error": abs(wall["gamma"]**2-wall["gamma_v"]**2-1),
            "cap_reconstruction_max_error": max(abs(wall["gamma"]*wall["Rdot"]-wall["gamma_v"]*wall["beta_minus"]-wall["V_minus"]),
                                                abs(wall["gamma"]*wall["beta_minus"]-wall["gamma_v"]*wall["Rdot"]-1)),
            "surface_energy_projection_ell_error": abs(slice_tensor[0, 0]-energy)*ell,
            "surface_momentum_projection_ell_error": abs(-slice_tensor[0, 1]-momentum)*ell,
            "shell_hamiltonian_delta_ell_error": abs(-4*jump_Rprime/Rw-16*math.pi*energy)*ell,
            "shell_momentum_delta_ell_error": abs(-2*jump_B-8*math.pi*momentum)*ell,
            "vacuum_hamiltonian_max_ell2": ell**2*_maxabs(fields["hamiltonian_residual"]),
            "vacuum_momentum_max_ell2": ell**2*_maxabs(fields["momentum_residual"]),
            "schwarzschild_mass_max_error_over_ell": max(_maxabs(fields["mass"]-M),_maxabs(bridge["mass"]-M))/ell,
            "embedding_spatial_norm_max_error": _maxabs(fields["embedding_spatial_norm"]-1),
            "embedding_future_normal_norm_max_error": _maxabs(fields["embedding_future_normal_norm"]+1),
            "bridge_hamiltonian_max_ell2": ell**2*_maxabs(bridge["scalar_curvature"]),
            "transition_32_64_distance_max_error_over_ell": _maxabs(distance32-distance64)/ell,
            "quadrature_comparison_is_rigorous_bound": False,
            "finite_difference_constraints": _finite_difference_checks(wall),
        },
        "source_reading": [
            {"id": "gr-qc/9403003v1", "url": "https://arxiv.org/html/gr-qc/9403003v1",
             "depth": "Sections2.2,3.1-3.2,4.1 Eq14-18,23,27-33,81-82",
             "use": "spherical slice curvature, extrinsic curvature, constraints and vacuum mass"},
            {"id": "gr-qc/9708008v2", "url": "https://arxiv.org/html/gr-qc/9708008v2",
             "depth": "Section2 stress/constraint definitions and Section4.1 frame boosts",
             "use": "surface stress projection and frame conventions, not full vacuum-only reduction"},
            {"id": "0711.2198v1", "url": "https://arxiv.org/html/0711.2198v1",
             "depth": "v3 imported Section2 Eq2 and Section3", "use": "fixed signed vacuum-wall trajectory"},
        ],
        "limits": ["initial-slice connection, not formation or past completeness", "no perturbative stability claim",
                   "one spherical sourced body, not the BL two-hole data", "no causal sibling communication proved",
                   "no new-to-world theorem or priority claim", "no observer-age, inheritance or local dark-sector prediction"],
    }
