"""Named leading-Einstein EFT diagnostic on the existing carrier and source.

Four independent real geometry variables; no chi or p_chi degree of freedom.
The historical full auxiliary action is never evaluated by this branch.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import time

import numpy as np
from threadpoolctl import threadpool_limits

from . import nsc_discovery_backend as backend
from . import nsc_discovery_dynamic_preparation as preparation
from . import nsc_discovery_episode as episode
from . import nsc_discovery_extent as extent
from . import nsc_discovery_tidal as tidal
from . import nsc_nested_parent_child as nested
from . import nsc_regional_energy_exchange as regional
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_galerkin_coupling as galerkin

SCHEMA = "NSC-DISCOVERY-LEADING-EINSTEIN-v1"
FIELDS = ("Q", "r", "p_Q", "p_r", "phi0", "phi1")
REAL_FIELDS = FIELDS[:4]
LAB = episode.LAB
ROOT = episode.REPO
INITIAL = LAB / "results/development/nsc-discovery-dynamic-preparation-v2.json"
OUTPUT = LAB / "results/development/nsc-discovery-leading-einstein-v1"
STATIONS = (.3, 1., 3.)
MAX_CPU = 600.
MAX_CHUNK = 64*1024**2


@dataclass
class State:
    """Nested Q/r coefficients and pi=dx_g W.T p; columns stay in the AP band."""
    Q: np.ndarray
    r: np.ndarray
    p_Q: np.ndarray
    p_r: np.ndarray
    phi0: np.ndarray
    phi1: np.ndarray

    def copy(self):
        return State(*(np.array(getattr(self, name), copy=True) for name in FIELDS))


@dataclass
class Rate(State):
    force_L: np.ndarray
    force_Q: np.ndarray
    force_beta: np.ndarray
    fieldwork_power: float


def coefficients(system):
    a = -8*np.pi*system.A
    return a, 3*a, 2*np.pi*system.C_F*system.flux**2


def encode(pair, nodal):
    W, dx = pair.geometry_map, pair.grid.dx_g
    return State(W.T@nodal.Q, W.T@nodal.r, dx*(W.T@nodal.p_Q), dx*(W.T@nodal.p_r),
                 np.array(nodal.phi0, copy=True), np.array(nodal.phi1, copy=True))


def decode(pair, state):
    W, dx = pair.geometry_map, pair.grid.dx_g
    return State(W@state.Q, W@state.r, W@state.p_Q/dx, W@state.p_r/dx,
                 np.array(state.phi0, copy=True), np.array(state.phi1, copy=True))


def prolong(grid, nodal):
    return State(*(galerkin.prolong_geometry(grid, getattr(nodal, name)) for name in REAL_FIELDS),
                 grid.U_f@nodal.phi0, grid.U_f@nodal.phi1)


def fine_state(pair, state):
    return prolong(pair.grid, decode(pair, state))


def check_chart(pair, state):
    fine = fine_state(pair, state)
    if not all(np.isfinite(getattr(fine, name)).all() for name in FIELDS):
        raise coupling.PositiveChartExit("nonfinite_leading_state", np.nan, state.copy())
    if min(float(np.min(fine.Q)), float(np.min(fine.r))) <= 0:
        raise coupling.PositiveChartExit("nonpositive_leading_chart", np.nan, state.copy())
    return fine


def active_system(grid, fine):
    if grid.gauge != "conformal":
        raise ValueError("leading Einstein diagnostic requires conformal gauge L=Q,beta=0")
    return galerkin.active_fine_system(grid, fine)


def require_conformal(system, state):
    if not np.array_equal(np.asarray(system.length_density), state.Q) or np.any(np.asarray(system.shift) != 0):
        raise ValueError("leading fine evaluator requires L=Q,beta=0")


def geometric_density(system, state):
    """Actual finite SBP Hamiltonian; D(F) is never replaced by b D(r)."""
    a, Z, mag = coefficients(system)
    b, Q, p, v = a*state.r, state.Q, state.p_Q, state.p_r
    D = system.derivative
    F, V = a*state.r**2/2, -a*state.r**2-mag
    return Q*p*v/(2*b)-Z*Q**2*p**2/(4*b**2)+Z*(D@state.r)**2-Q**2*V+2*(D@F)*(D@Q)/Q


def fine_rates(system, state, source=None):
    """Canonical derivative of the finite Hamiltonian, including the lapse chain."""
    require_conformal(system, state)
    a, Z, mag = coefficients(system)
    Q, r, p, v, D = state.Q, state.r, state.p_Q, state.p_r, system.derivative
    b, F, V = a*r, a*r*r/2, -a*r*r-mag
    rx, qx, Fx = D@r, D@Q, D@F
    if source is None:
        source = coupling.source_from_columns(system, state)
    qdot = Q*v/(2*b)-Z*Q*Q*p/(2*b*b)
    rdot = Q*p/(2*b)
    pdot = (-p*v/(2*b)+Z*Q*p*p/(2*b*b)+2*Q*V
            +2*(D@(Fx/Q))+2*Fx*qx/(Q*Q)
            -(source["force_Q"]+source["force_L"])/system.dx)
    vdot = (a*Q*p*v/(2*b*b)-a*Z*Q*Q*p*p/(2*b**3)
            +2*Z*(D@rx)+Q*Q*(-2*a*r)+2*b*(D@(qx/Q)))
    return qdot, rdot, pdot, vdot, -1j*source["image0"], -1j*source["image1"]


def rates(pair, state, *, return_bundle=False, control_mode="coupled"):
    mode = episode.normalize_control_mode(control_mode)
    fine = check_chart(pair, state)
    system = active_system(pair.grid, fine)
    source = coupling.source_from_columns(system, fine)
    raw = fine_rates(system, fine, source)
    forcing = raw[:4]
    if mode=="frozen_geometry":
        raw = tuple(np.zeros_like(value) for value in raw[:4])+raw[4:]
    nodal = State(*(galerkin.pull_geometry(pair.grid, value) for value in raw[:4]),
                  pair.grid.U_f.conj().T@raw[4], pair.grid.U_f.conj().T@raw[5])
    encoded = encode(pair, nodal)
    lifted_Q = galerkin.prolong_geometry(pair.grid, nodal.Q)
    work = float(np.sum((source["force_Q"]+source["force_L"])*lifted_Q))
    rate = Rate(*(getattr(encoded, name) for name in FIELDS),
                source["force_L"], source["force_Q"], source["force_beta"], work)
    bundle = {"fine_state": fine, "fine_system": system, "source": source, "nodal_rate": nodal,
              "unprojected_rates": raw, "lifted_Q": lifted_Q, "lapse_dot": lifted_Q,
              "shift_dot": np.zeros(pair.grid.nq), "control_mode": mode,
              "coupled_geometry_forcing": forcing}
    return (rate, bundle) if return_bundle else rate


def fine_jvp(system, s, t):
    """Analytic derivative of the new four-variable SBP rate and unchanged field."""
    require_conformal(system, s)
    a, Z, mag = coefficients(system)
    Q, r, p, v = s.Q, s.r, s.p_Q, s.p_r
    q, u, dp, dv = t.Q, t.r, t.p_Q, t.p_r
    b, db, D = a*r, a*u, system.derivative
    F, dF, V, dV = a*r*r/2, b*u, -a*r*r-mag, -2*a*r*u
    Fx, dFx, Qx, dQx = D@F, D@dF, D@Q, D@q
    dqdot = ((q*v+Q*dv)/(2*b)-Q*v*db/(2*b*b)
             -Z*(2*Q*q*p+Q*Q*dp)/(2*b*b)+Z*Q*Q*p*db/b**3)
    drdot = (q*p+Q*dp)/(2*b)-Q*p*db/(2*b*b)
    ds = 2*np.sum(np.real(np.conj(t.phi0)*s.phi1+np.conj(s.phi0)*t.phi1)*system.occupations, axis=1)
    dpdot = (-(dp*v+p*dv)/(2*b)+p*v*db/(2*b*b)
             +Z*(q*p*p+2*Q*p*dp)/(2*b*b)-Z*Q*p*p*db/b**3+2*q*V+2*Q*dV
             +2*(D@(dFx/Q-Fx*q/Q**2))+2*(dFx*Qx+Fx*dQx)/Q**2-4*Fx*Qx*q/Q**3
             -system.multiplicity*system.kappa*ds/system.dx)
    dvdot = (a*(q*p*v+Q*dp*v+Q*p*dv)/(2*b*b)-a*Q*p*v*db/b**3
             -a*Z*(2*Q*q*p*p+2*Q*Q*p*dp)/(2*b**3)+3*a*Z*Q*Q*p*p*db/(2*b**4)
             +2*Z*(D@(D@u))-4*a*Q*r*q-2*a*Q*Q*u
             +2*db*(D@(Qx/Q))+2*b*(D@(dQx/Q-Qx*q/Q**2)))
    P, k = system.momentum, system.kappa
    dphi0 = -(P@t.phi1)-1j*k*(q[:, None]*s.phi1+Q[:, None]*t.phi1)
    dphi1 = P@t.phi0-1j*k*(q[:, None]*s.phi0+Q[:, None]*t.phi0)
    return State(dqdot, drdot, dpdot, dvdot, dphi0, dphi1)


def jvp(pair, state, tangent, *, control_mode="coupled"):
    fine = check_chart(pair, state)
    lifted = prolong(pair.grid, decode(pair, tangent))
    raw = fine_jvp(active_system(pair.grid, fine), fine, lifted)
    if episode.normalize_control_mode(control_mode)=="frozen_geometry":
        for name in REAL_FIELDS:
            setattr(raw, name, np.zeros_like(getattr(raw, name)))
    nodal = State(*(galerkin.pull_geometry(pair.grid, getattr(raw, name)) for name in REAL_FIELDS),
                  pair.grid.U_f.conj().T@raw.phi0, pair.grid.U_f.conj().T@raw.phi1)
    return encode(pair, nodal)


def energy(pair, state):
    fine = check_chart(pair, state)
    system = active_system(pair.grid, fine)
    gravity = float(pair.grid.dx_q*np.sum(geometric_density(system, fine)))
    field = coupling.field_energy(system, fine)
    return {"gravity": gravity, "field": float(field), "total": gravity+float(field)}


def constraint_arrays(grid, fine, system, source=None):
    if source is None:
        source = coupling.source_from_columns(system, fine)
    a, Z, mag = coefficients(system)
    b, Q, p, v, D = a*fine.r, fine.Q, fine.p_Q, fine.p_r, system.derivative
    F, V = a*fine.r**2/2, -a*fine.r**2-mag
    Cg = p*v/(2*b)-Z*Q*p*p/(4*b*b)+Z*(D@fine.r)**2/Q-Q*V-2*(D@((D@F)/Q))
    Dg = v*(D@fine.r)-Q*(D@p)
    return Cg+source["force_L"]/grid.dx_q, Dg+source["force_beta"]/grid.dx_q


def constraints(pair, state):
    fine = check_chart(pair, state)
    system = active_system(pair.grid, fine)
    source = coupling.source_from_columns(system, fine)
    C, shift = constraint_arrays(pair.grid, fine, system, source)
    h = fine.Q*C
    project = galerkin.pull_geometry(pair.grid, h)
    maximum = lambda value: float(np.max(abs(value)))
    return {"raw_C_max": maximum(C), "h_c_max": maximum(h), "D_max": maximum(shift),
            "projected_h_c_max": maximum(project),
            "held_out_h_c_max": maximum(h-galerkin.prolong_geometry(pair.grid, project)),
            "source_Q_rho_max": maximum(fine.Q*source["force_L"]/pair.grid.dx_q),
            "source_current_max": maximum(source["force_beta"]/pair.grid.dx_q),
            "mean_current": float(np.mean(source["force_beta"]/pair.grid.dx_q)),
            "observable_error_bound": None, "physical_instability_inferred": False}


def metric_jets(pair, state, rate=None, bundle=None, *, control_mode="coupled"):
    if rate is None or bundle is None:
        rate, bundle = rates(pair, state, return_bundle=True, control_mode=control_mode)
    mode = bundle.get("control_mode", episode.normalize_control_mode(control_mode))
    acceleration = prolong(pair.grid, decode(pair, jvp(pair, state, rate, control_mode=mode)))
    velocity = prolong(pair.grid, decode(pair, rate))
    fine = bundle["fine_state"]
    dx = lambda values: tidal.metric.spectral_dx(values, pair.grid.length)
    qx, rx = dx(fine.Q), dx(fine.r)
    curvature = tidal.curvature_from_local_jets(
        fine.Q, fine.r, velocity.Q, qx, acceleration.Q, dx(velocity.Q), dx(qx),
        velocity.r, rx, acceleration.r, dx(velocity.r), dx(rx))
    return curvature, velocity, acceleration


def combine(state, tangent, scale):
    return State(*(getattr(state, name)+scale*getattr(tangent, name) for name in FIELDS))


def rk4_step(pair, state, dt, *, control_mode="coupled"):
    mode = episode.normalize_control_mode(control_mode)
    k1 = rates(pair, state, control_mode=mode)
    k2 = rates(pair, combine(state, k1, dt/2), control_mode=mode)
    k3 = rates(pair, combine(state, k2, dt/2), control_mode=mode)
    k4 = rates(pair, combine(state, k3, dt), control_mode=mode)
    result = State(*(getattr(state, name)+dt*(getattr(k1, name)+2*getattr(k2, name)
                                             +2*getattr(k3, name)+getattr(k4, name))/6 for name in FIELDS))
    if mode=="frozen_geometry":
        for name in REAL_FIELDS:
            setattr(result, name, np.array(getattr(state, name), copy=True))
    check_chart(pair, result)
    return result


def frozen_geometry_step(pair, state, dt):
    return rk4_step(pair, state, dt, control_mode="frozen_geometry")


def step_restriction(pair, state, step_cap):
    """Unit-speed principal waves plus a conservative local coupled reaction bound."""
    fine = check_chart(pair, state)
    system = active_system(pair.grid, fine)
    a, Z, mag = coefficients(system)
    # Local reaction Jacobian: evaluate the analytic Jv with zero spatial
    # operators. Principal Fourier symbols are treated separately.
    class Zero:
        def __matmul__(self, values):
            return np.zeros_like(values)
    local_system = SimpleNamespace(A=system.A, C_F=system.C_F, flux=system.flux,
        derivative=Zero(), momentum=Zero(), kappa=system.kappa,
        occupations=system.occupations, multiplicity=system.multiplicity, dx=system.dx,
        length_density=fine.Q, shift=np.zeros(pair.grid.nq))
    row_sums = np.zeros((pair.grid.nq, 4))
    zeros = np.zeros(pair.grid.nq)
    zcol = np.zeros_like(fine.phi0)
    for name in REAL_FIELDS:
        tangent = State(*(np.ones_like(zeros) if field==name else zeros for field in REAL_FIELDS), zcol, zcol)
        image = fine_jvp(local_system, fine, tangent)
        row_sums += np.column_stack([abs(getattr(image, field)) for field in REAL_FIELDS])
    psi0, psi1 = fine.phi0/np.sqrt(pair.grid.dx_q), fine.phi1/np.sqrt(pair.grid.dx_q)
    row_sums[:, 2] += 2*system.multiplicity*abs(system.kappa)*np.sum(
        (abs(psi0.real)+abs(psi0.imag)+abs(psi1.real)+abs(psi1.imag))*system.occupations, axis=1)
    field_reaction = abs(system.kappa)*(np.max(abs(fine.Q))+max(np.max(abs(psi0)), np.max(abs(psi1))))
    reaction = max(float(np.max(row_sums)), float(field_reaction))
    k_geometry = episode.quadrature_principal_wavenumber(pair.grid)
    k_field = float(2*np.pi*np.max(abs(pair.grid.modes_f))/pair.grid.length)
    gradient = 2*max(float(np.max(abs(system.derivative@fine.Q/fine.Q))),
                     float(np.max(abs(system.derivative@fine.r/fine.r))))
    omega = max(k_geometry, k_field)+reaction+gradient
    dt = min(float(step_cap), episode.RK4_HALF_STABILITY/max(omega, 1.))
    return dt, {"dt": dt, "step_cap": float(step_cap), "combined_omega": omega,
                "principal_coordinate_speed": 1., "principal_AB_identity": True,
                "reaction_majorant": reaction, "gradient_majorant": gradient,
                "field_cfl_alone": False, "stability_certificate": False,
                "restriction_scope": "unit-speed frozen principal part plus local coupled/source and gradient majorants"}


def clock_rates(pair, state):
    fine = check_chart(pair, state)
    return extent.real_periodic_values(pair.grid, fine.r*fine.Q, pair.clock_locations)


def observe(pair, state, instant, control_mode="coupled"):
    mode = episode.normalize_control_mode(control_mode)
    rate, bundle = rates(pair, state, return_bundle=True, control_mode=mode)
    fine, system = bundle["fine_state"], bundle["fine_system"]
    curvature, velocity, _acceleration = metric_jets(pair, state, rate, bundle)
    ledger = regional.matter_ledger(system, fine)
    balance = regional.proper_balance_terms(system, fine, velocity, ledger)
    grid = pair.grid
    integral = lambda values, interval: extent.real_interval_integral(grid, values, interval)
    sample = lambda values, locations: extent.real_periodic_values(grid, values, locations)
    mass = np.sum((abs(fine.phi0)**2+abs(fine.phi1)**2)*pair.weights, axis=1)
    proper = fine.r*fine.Q
    density = mass/(grid.dx_q*proper)
    mean = float(np.sum(mass)/integral(proper, (0., grid.length)))
    contrast = float((np.max(density)-np.min(density))/mean) if mean else 0.
    gram = state.phi0.conj().T@state.phi0+state.phi1.conj().T@state.phi1
    eig = np.linalg.eigvalsh(np.sqrt(pair.weights[:, None]*pair.weights[None, :])*gram)
    child = pair.child_interval
    boundary = sample(balance["flux_nodal"]/grid.dx_q, child)
    parent_boundary = sample(balance["flux_nodal"]/grid.dx_q, pair.parent_interval)
    summary = lambda values: {"min": float(np.min(values)), "max": float(np.max(values)),
                              "max_abs": float(np.max(abs(values)))}
    return {"time": float(instant), "control_mode": mode, "model": "leading_Einstein_EFT_diagnostic",
        "child_proper_length": integral(proper, child), "child_r_proper_mean": integral(proper*fine.r, child)/integral(proper, child),
        "child_probability": integral(mass/grid.dx_q, child), "parent_probability": integral(mass/grid.dx_q, pair.parent_interval),
        "probability_contrast": contrast, "child_normal_energy": integral(ledger["normal_energy_nodal"]/grid.dx_q, child),
        "normal_clock_rates": clock_rates(pair, state).tolist(),
        "coordinate_fieldwork": rate.fieldwork_power,
        "pressure_work": float(np.sum(balance["proper_pressure_work"])),
        "lapse_work": float(np.sum(balance["momentum_lapse_work"])),
        "boundary_child": float(boundary[0]-boundary[1]),
        "boundary_parent": float(parent_boundary[0]-parent_boundary[1]),
        "energy": energy(pair, state), "constraints": constraints(pair, state),
        "Q": summary(fine.Q), "r": summary(fine.r), "metric_N": summary(proper),
        "tides": {name: summary(curvature["tides"][name]) for name in ("R_0101", "R_0202", "R4", "owned_W")},
        "curvature": {name: summary(curvature["base"][name]) for name in ("Ricci2", "K")},
        "Gram_max": float(np.max(abs(gram-np.eye(6)))), "CAR_min": float(eig.min()), "CAR_max": float(eig.max()),
        "source_current_mean_deleted": False, "independent_auxiliary_variables": False,
        "controlled_metric": mode=="frozen_geometry", "constraint_preservation_claim": False,
        "constraint_scope": "initial own leading constraint solution; later frozen residuals and inactive geometry forcing are measured"
                            if mode=="frozen_geometry" else "actual coupled leading constraint residuals",
        "inactive_geometry_forcing": None if mode=="coupled" else
            {name: summary(value) for name, value in zip(REAL_FIELDS, bundle["coupled_geometry_forcing"])},
        "coordinate_readout_is_equal_proper_time": False,
        "full_first_CW_correction_implemented": False, "holding_claim": False}


def source_hashes():
    paths = set(preparation.OWNERS[4:])
    paths.update("lab/"+directory+"/"+name for directory, name in (
        ("src/recursive_horizons", "nsc_discovery_leading_einstein.py"),
        ("scripts", "derive_nsc_discovery_leading_einstein.py"),
        ("tests", "test_nsc_discovery_leading_einstein.py"),
        ("docs", "nsc-discovery-leading-einstein.md")))
    modules = (backend, preparation, episode, extent, tidal, nested, regional, coupling, galerkin)
    paths.update(str(Path(module.__file__).resolve().relative_to(ROOT)) for module in modules)
    paths.update(str(Path(path).relative_to(ROOT)) for path in tidal.producer_hashes())
    return {name: episode.file_sha256(ROOT/name) for name in sorted(paths)}


def verify_sources(record):
    for name, digest in record["source_hashes"].items():
        if episode.file_sha256(ROOT/name) != digest:
            raise ValueError("leading diagnostic source changed: "+name)
    for name, digest in record.get("input_hashes", {}).items():
        if episode.file_sha256(ROOT/name) != digest:
            raise ValueError("frozen initial input changed: "+name)


def transfer_initial(case, nf, initial=INITIAL):
    """Frozen fields/Q only; verify removed variables vanish and reprepare C0."""
    pair, old, record = preparation.load_case(initial, case)
    old_nodal = nested.reconstruct_state(pair, old)
    for name in ("chi", "p_chi", "p_Q", "p_r"):
        if np.any(getattr(old_nodal, name) != 0):
            raise ValueError("leading turnover transfer requires zero "+name)
    if nf == pair.grid.nf:
        nodal = State(old_nodal.Q.copy(), old_nodal.r.copy(), old_nodal.p_Q.copy(),
                      old_nodal.p_r.copy(), old_nodal.phi0.copy(), old_nodal.phi1.copy())
        return pair, nodal, {"source_nf": nf, "AP_extension": False, "frozen_Q_and_columns_preserved": True}
    if nf < pair.grid.nf:
        raise ValueError("frozen source transfer only supports exact band extension")
    phi_map = backend.SpinorCarrier(pair.grid.nf, nf, pair.grid.length, canonical=True)
    phi0, phi1 = phi_map@old_nodal.phi0, phi_map@old_nodal.phi1
    new_pair = nested.build_pair(nf, coarse_modes=1, child_details=2,
        source_layout="override", columns_override=(phi0, phi1), occupations=pair.weights)
    modes = galerkin.geometry_modes(pair.grid.ng)
    analysis = np.exp(-2j*np.pi*pair.grid.xi_g[:, None]*modes/pair.grid.length)/pair.grid.ng
    synthesis = np.exp(2j*np.pi*new_pair.grid.xi_g[:, None]*modes/pair.grid.length)
    transform = lambda values: np.real(synthesis@(analysis.T@values))
    nodal = State(transform(old_nodal.Q), transform(old_nodal.r), np.zeros(new_pair.grid.ng),
                  np.zeros(new_pair.grid.ng), phi0, phi1)
    return new_pair, nodal, {"source_nf": pair.grid.nf, "AP_extension": True,
                           "no_new_eigenframe_or_occupation_selection": True}


def prepare_radius(pair, nodal, deadline):
    """Same source-consistent zero-momentum C0, evaluated by this branch only."""
    if np.any(nodal.p_Q != 0) or np.any(nodal.p_r != 0):
        raise ValueError("leading C0 radius preparation requires zero momenta")
    grid = pair.grid
    fine = prolong(grid, nodal)
    system = active_system(grid, fine)
    source = coupling.source_from_columns(system, fine)
    rho = source["force_L"]/grid.dx_q
    D, U = system.derivative, grid.A_g
    DU = D@U
    a, Z, _mag = coefficients(system)
    def oracle(radius):
        trial = nodal.copy()
        trial.r = radius
        lifted = prolong(grid, trial)
        full, _ = constraint_arrays(grid, lifted, system, source)
        b = a*lifted.r
        image = (2*Z*(D@lifted.r/lifted.Q)[:, None]*DU
                 +2*a*(lifted.Q*lifted.r)[:, None]*U
                 -2*D@((D@(b[:, None]*U))/lifted.Q[:, None]))
        return galerkin.pull_geometry(grid, full), galerkin.pull_geometry(grid, image)
    nodal.r, report = preparation._positive_newton(nodal.r, oracle,
        lambda value: U@value, deadline)
    state = encode(pair, nodal)
    measured = constraints(pair, state)
    report.update(leading_constraints=measured, old_exact_rates_called=False,
                  same_zero_momentum_lapse_equation=True)
    if not report["converged"]:
        raise ValueError("leading initial radius unresolved: "+str(report["blocker"]))
    return state, report


def plan():
    cases = [{"case_id": f"nf{nf}_{name}_dt{cap}", "nf": nf, "initial_case": name,
              "step_cap": float(cap), "control_mode": "coupled"}
             for nf, name in ((128, "uniform"), (128, "coherent"), (256, "coherent"))
             for cap in ("0.001", "0.0005")]
    return cases+[{"case_id": f"nf128_coherent_frozen_geometry_dt{cap}", "nf": 128,
                   "initial_case": "coherent", "step_cap": float(cap), "control_mode": "frozen_geometry"}
                  for cap in ("0.001", "0.0005")]


def state_from_arrays(arrays, representation=episode.CANONICAL_PI):
    if representation != episode.CANONICAL_PI:
        raise ValueError("leading checkpoints store canonical frame pi")
    if "chi" in arrays or "pi_chi" in arrays:
        raise ValueError("leading checkpoint has an independent auxiliary variable")
    return State(arrays["Q"], arrays["r"], arrays["pi_Q"], arrays["pi_r"], arrays["phi0"], arrays["phi1"])


def real_array_names(representation):
    if representation != episode.CANONICAL_PI:
        raise ValueError("leading checkpoints store canonical frame pi")
    return ("Q", "r", "pi_Q", "pi_r", "source_weights", "clock_rates", "normal_clocks", "W")


def normalize_stations(values=None):
    values = STATIONS if values is None else tuple(float(value) for value in values)
    if not values or list(values) != sorted(set(values)) or min(values) <= 0 or max(values) > 3:
        raise ValueError("leading stations must be ordered positive times through T3")
    return values


def scalar_row(pair, state, instant, mode="coupled"):
    started = time.process_time()
    row = observe(pair, state, instant, mode)
    fine = fine_state(pair, state)
    src = coupling.source_from_columns(active_system(pair.grid, fine), fine)
    row["source_rho"] = {"min": float(np.min(src["force_L"]/pair.grid.dx_q)),
                         "max": float(np.max(src["force_L"]/pair.grid.dx_q))}
    row.update(observation_cpu_seconds=time.process_time()-started,
        dense_propagator_stored=False, quadrature="trapezoid_coordinate_time",
        normal_clocks=None if preparation._observer_normal_clocks is None
                      else preparation._observer_normal_clocks.tolist(),
        clock_values_source="checkpoint-carried endpoint trapezoid; event wrapper updates endpoint rows")
    return row


def episode_observe(pair, state, instant, *, control_mode="coupled"):
    return {"stability": scalar_row(pair, state, instant, control_mode)}


def estimate_case_cpu(pair, state, *, step_cap, stations, current_time, probe_step_cpu, probe_diag_cpu):
    dt, _ = step_restriction(pair, state, step_cap)
    span = max(normalize_stations(stations)[-1]-current_time, 0.)
    steps = int(np.ceil(span/dt))
    diagnostics = int(np.ceil(span/episode.OBSERVATION_CADENCE))+len(stations)+1
    return {"dt": dt, "estimated_steps": steps, "estimated_cpu_seconds":
            steps*probe_step_cpu+diagnostics*probe_diag_cpu}


def execute_case(directory, case_id, *, cpu_allowance, forecast_factor=1.5,
                 memory_limit_bytes=episode.DEFAULT_MEMORY_BYTES, max_steps=None, fft=None,
                 backend=None, diagnostics=False, **kwargs):
    record, arrays = episode.load_checkpoint(directory, case_id)
    pair = episode.pair_from_arrays(arrays, record)
    state = state_from_arrays(arrays)
    pair, state, info = episode.resolve_pair(pair, state, fft=fft, backend=backend)
    mode = episode.normalize_control_mode(record.get("control_mode", "coupled"))
    stepper = frozen_geometry_step if mode=="frozen_geometry" else rk4_step
    started = time.process_time()
    result = episode.advance_case(pair, state, coordinate_time=record["coordinate_time"],
        steps=record["steps"], stations=record["stations"], step_cap=record["step_cap"],
        geometry_frozen=mode=="frozen_geometry", cpu_allowance=cpu_allowance, spent=0., forecast_factor=forecast_factor,
        memory_limit_bytes=memory_limit_bytes, max_steps=max_steps, stepper=stepper,
        normal_clocks=arrays["normal_clocks"], control_mode=mode, backend=backend,
        work_ledger=record.get("work_ledger"), channel_sample=record.get("channel_sample"),
        channel_time=record.get("channel_time"), **kwargs)
    elapsed = time.process_time()-started
    episode.append_observations(directory, case_id, result["observations"])
    chunk = None
    for snapshot in result["snapshots"]:
        at_end = abs(snapshot["time"]-result["time"]) <= 1e-12
        changed = snapshot["time"] != record["coordinate_time"] or any(
            not np.array_equal(getattr(snapshot["state"], name), getattr(state, name)) for name in FIELDS)
        values = dict(arrays)
        for name, key in (("Q", "Q"), ("r", "r"), ("p_Q", "pi_Q"), ("p_r", "pi_r"), ("phi0", "phi0"), ("phi1", "phi1")):
            values[key] = np.array(getattr(snapshot["state"], name), copy=True)
        values.update(normal_clocks=snapshot["clocks"], clock_rates=snapshot["clock_rates"])
        updated = dict(record, coordinate_time=snapshot["time"], steps=snapshot["steps"],
            stations_reached=snapshot["stations_reached"], work_ledger=snapshot["ledger"],
            channel_sample=snapshot["channel_sample"], channel_time=snapshot["channel_time"],
            snapshot_kind=snapshot["kind"], status=result["status"] if at_end else "station_retained",
            stop=result["stop"] if at_end else None, evolved=changed,
            model="leading_Einstein_EFT_diagnostic", independent_auxiliary_variables=False)
        updated["source_pins_after"] = episode._live_pins(values, snapshot["state"])
        episode._pins_unchanged(record["source_pins"], updated["source_pins_after"])
        committed = episode.commit_checkpoint(directory, updated, values, limit=MAX_CHUNK)
        chunk = committed["npz"]
    return {"case_id": case_id, "status": result["status"], "coordinate_time": result["time"],
            "steps": result["steps"], "stations_reached": result["stations_reached"],
            "stop": result["stop"], "child_cpu_seconds": elapsed, "chunk": chunk,
            "stability_certificate": False, "state_clamped": False, "physical_instability_claimed": False}


@contextmanager
def episode_adapter():
    """Scoped branch hooks; reuse the existing runner, event wrapper, pool and I/O."""
    if episode.execute_case is execute_case:
        yield episode
        return
    replacements = {
        "SCHEMA": SCHEMA, "model": SimpleNamespace(STATE_NAMES=FIELDS, NestedState=State,
          NestedPair=nested.NestedPair, reconstruct_state=decode),
        "_real_array_names": real_array_names, "state_from_arrays": state_from_arrays,
        "normalize_stations": normalize_stations, "evolving_step": rk4_step,
        "frozen_geometry_step": frozen_geometry_step,
        "step_restriction": step_restriction, "_clock_rates": clock_rates,
        "scalar_observation_row": scalar_row, "observe": episode_observe,
        "estimate_case_cpu": estimate_case_cpu, "execute_case": execute_case,
        "advance_case": preparation._event_aware_advance,
    }
    saved = {name: getattr(episode, name) for name in replacements}
    old_advance = getattr(preparation, "_advance_owner", None)
    old_clocks = preparation._observer_normal_clocks
    preparation._advance_owner = episode.advance_case
    for name, value in replacements.items():
        setattr(episode, name, value)
    try:
        yield episode
    finally:
        for name, value in saved.items():
            setattr(episode, name, value)
        preparation._advance_owner = old_advance
        preparation._observer_normal_clocks = old_clocks


def prepare(initial=INITIAL, output=OUTPUT, *, execute=False, producer_commit=None, cpu_budget=MAX_CPU):
    if not 0 < cpu_budget <= MAX_CPU:
        raise ValueError("leading preflight/preparation budget must be at most600 CPU seconds")
    report = {"schema": SCHEMA, "status": "PREFLIGHT", "cases": plan(), "stations": list(STATIONS),
        "scope": "leading Einstein diagnostic of curvature EFT; full first-CW correction and completed theory not claimed",
        "active_geometry_variables": list(REAL_FIELDS), "independent_chi_or_p_chi": False,
        "C_W_retained_only_in_input_history": True, "cpu_budget_seconds": float(cpu_budget),
        "workers": 4, "numerical_threads_per_worker": 1, "chunk_limit_bytes": MAX_CHUNK,
        "evolved": False, "source_hashes": source_hashes()}
    if not execute:
        return report
    if not 0 < cpu_budget <= MAX_CPU or not producer_commit:
        raise ValueError("explicit preparation requires a frozen producer commit and budget at most600 CPU seconds")
    destination = episode.assert_campaign_output(output)
    if destination.exists():
        raise FileExistsError("leading output already exists")
    commit = preparation._git_hashes(producer_commit, report["source_hashes"])
    binding = preparation.authenticate_preparation(initial)
    jp, payload = preparation.output_paths(initial)
    report.update(producing_commit=commit, initial_binding=binding,
        input_hashes={str(jp.resolve().relative_to(ROOT)): episode.file_sha256(jp),
                      str(payload.resolve().relative_to(ROOT)): episode.file_sha256(payload)})
    started = time.process_time()
    prepared = {}
    with threadpool_limits(limits=1), backend.fft_thread_limit(1):
        for nf, name in ((128, "uniform"), (128, "coherent"), (256, "coherent")):
            pair, nodal, transfer = transfer_initial(name, nf, initial)
            state, correction = prepare_radius(pair, nodal, started+cpu_budget)
            prepared[(nf, name)] = pair, state, transfer, correction
    verify_sources(report)
    cases = []
    with episode_adapter():
        for spec in plan():
            pair, state, transfer, correction = prepared[(spec["nf"], spec["initial_case"])]
            arrays = {name: np.array(getattr(state, field), copy=True) for field, name in
                      (("Q", "Q"), ("r", "r"), ("p_Q", "pi_Q"), ("p_r", "pi_r"), ("phi0", "phi0"), ("phi1", "phi1"))}
            arrays.update(W=pair.geometry_map, source_phi0=pair.source_phi0, source_phi1=pair.source_phi1,
                observer_columns=pair.reference_columns, source_weights=pair.weights,
                normal_clocks=np.zeros(3), clock_rates=clock_rates(pair, state))
            pins = episode._live_pins(arrays, state)
            record = dict(spec, coordinate_time=0., steps=0, ordinal=0, status="PREPARED", stations=list(STATIONS),
                control_mode=spec["control_mode"], momentum_representation=episode.CANONICAL_PI,
                snapshot_kind="handoff", source_pins=pins, verify_external_pins=False,
                coarse_indices=pair.geometry_coarse_indices.tolist(), child_indices=pair.geometry_child_indices.tolist(),
                parent_indices=pair.geometry_parent_indices.tolist(), source_metadata=pair.source_metadata,
                geometry_metadata=pair.geometry_metadata, clock_locations=list(pair.clock_locations),
                initial_state_called=False, leading_initial_constraints_prepared=True,
                active_geometry_variables=list(REAL_FIELDS), model="leading_Einstein_EFT_diagnostic",
                transfer=transfer, initial_constraint_correction=correction, evolved=False)
            chunk = episode.commit_checkpoint(destination, record, arrays, limit=MAX_CHUNK)
            cases.append({name: chunk[name] for name in ("case_id", "ordinal", "npz", "json", "coordinate_time", "steps")})
        used = time.process_time()-started
        if used > cpu_budget:
            raise RuntimeError("leading preparation exceeded shared CPU budget")
        report.update(status="PREPARED", stage=0, cases=cases, preparation_cpu_seconds=used,
                      child_cpu_seconds=0., forecast_factor=1.5)
        episode._write_json(episode._manifest_path(destination), report)
        episode._ledger_update(destination, pilot=used, budget=cpu_budget)
    return report


_worker_scope = None
_worker_threads = None
_worker_fft = None


def worker_initializer():
    global _worker_scope, _worker_threads, _worker_fft
    _worker_threads = threadpool_limits(limits=1)
    _worker_threads.__enter__()
    _worker_fft = backend.fft_thread_limit(1)
    _worker_fft.__enter__()
    _worker_scope = episode_adapter()
    _worker_scope.__enter__()


def pool(max_workers):
    return ProcessPoolExecutor(max_workers=max_workers, initializer=worker_initializer)


def comparisons(directory):
    rows = {}
    for spec in plan():
        rows[spec["case_id"]] = episode.read_observations(directory, spec["case_id"])
    channels = ("child_proper_length", "child_r_proper_mean", "child_probability", "probability_contrast",
                "radial_tide", "angular_tide", "clock_rate_x2", "energy_total", "Gram_max")
    def value(row, name):
        if name in ("radial_tide", "angular_tide"):
            return row["tides"]["R_0101" if name=="radial_tide" else "R_0202"]["max_abs"]
        if name=="clock_rate_x2":
            return row["normal_clock_rates"][1]
        if name=="energy_total":
            return row["energy"]["total"]
        return row[name]
    result = {}
    pairs = [(f"nf{nf}_{name}_dt0.001", f"nf{nf}_{name}_dt0.0005")
             for nf, name in ((128, "uniform"), (128, "coherent"), (256, "coherent"))]
    pairs += [(f"nf128_coherent_dt{cap}", f"nf256_coherent_dt{cap}") for cap in ("0.001", "0.0005")]
    pairs += [(f"nf128_coherent_dt{cap}", f"nf128_coherent_frozen_geometry_dt{cap}")
              for cap in ("0.001", "0.0005")]
    for left, right in pairs:
        a, b = {round(row["time"], 8): row for row in rows[left]}, {round(row["time"], 8): row for row in rows[right]}
        result[left+"__"+right] = {str(t): {name: {"difference": value(b[t], name)-value(a[t], name),
            "relative_to_second": None if name=="energy_total" or value(b[t], name)==0 else
                (value(b[t], name)-value(a[t], name))/abs(value(b[t], name))}
            for name in channels} for t in sorted(a.keys() & b.keys()) if t in STATIONS}
    return result


def run(output=OUTPUT, *, workers=4, cpu_budget=MAX_CPU, max_steps=None):
    if not 1 <= workers <= 4 or not 0 < cpu_budget <= MAX_CPU:
        raise ValueError("one pool of at most4 workers and at most600 CPU seconds required")
    manifest = episode.read_manifest(output)
    if manifest.get("schema") != SCHEMA:
        raise ValueError("not a leading-Einstein diagnostic checkpoint")
    verify_sources(manifest)
    preparation._git_hashes(manifest["producing_commit"], manifest["source_hashes"])
    with episode_adapter():
        result = episode.run(output, workers=workers, cpu_budget_seconds=cpu_budget, max_steps=max_steps,
                             backend="fft", executor=pool)
        result["comparisons"] = comparisons(output)
        result["comparison_time_domain"] = (
            "Matched coordinate times with separate checkpoint-carried normal clocks. "
            "Frozen/coupled readouts are coordinate-controlled; equal proper times are not asserted.")
        result["cumulative_cpu_seconds"] = json.loads(episode._ledger_path(output).read_text())["spent"]
        result["physical_scope"] = "generated leading-Einstein feedback comparison; exact auxiliary histories remain preserved"
        episode._write_json(episode._manifest_path(output), result)
    return result


def check(output=OUTPUT):
    manifest = episode.read_manifest(output)
    if manifest.get("schema") != SCHEMA:
        raise ValueError("not a leading-Einstein checkpoint")
    verify_sources(manifest)
    with episode_adapter():
        report = episode.check(output)
        for case in manifest["cases"]:
            record, arrays = episode.load_checkpoint(output, case["case_id"])
            pair = episode.pair_from_arrays(arrays, record)
            state = state_from_arrays(arrays)
            check_chart(pair, state)
            report.setdefault("latest_constraints", {})[case["case_id"]] = constraints(pair, state)
    return report
