"""Analytic Jacobian-vector action of the owned conformal Galerkin ODE.

All eight field blocks vary. The nonlinear evolution stays in its existing
owner; this module differentiates its L=Q, beta=0 Hamiltonian pullback.
It neither resets constraints nor assumes that later columns remain chiral.
"""
from __future__ import annotations

import numpy as np

from .nsc_spherical_coupling import CauchyState, PositiveChartExit, _combine, chart_failure
from .nsc_spherical_feedback_action import (
    alpha_of, feedback_F, feedback_V, feedback_Z, partial_F,
)
from .nsc_spherical_galerkin_coupling import (
    _check_state_dimensions, _constraint_directional, active_fine_system,
    compose_fine_hamiltonian, prolong_state, pull_geometry,
)

GEOMETRY = ("Q", "r", "chi", "p_Q", "p_r", "p_chi")
FIELDS = GEOMETRY + ("phi0", "phi1")


def conformal_jvp(grid, state, direction, *, bundle=None, checkpoint=None):
    """Return D X_H(state) direction, using the same maps as the full ODE.

    An optional bundle must come from compose_fine_hamiltonian at this very
    state. It avoids repeating the baseline prolongation during a paired
    evolution. Complex columns are differentiated as real coordinates.
    """
    if grid.gauge != "conformal":
        raise ValueError("this tangent requires the conformal gauge L=Q, beta=0")
    if int(grid.fine.multiplicity) != 4*int(grid.fine.kappa):
        raise ValueError("4 kappa multiplicity drifted")
    _check_state_dimensions(grid, state)
    _check_state_dimensions(grid, direction)
    for name in FIELDS:
        values = np.asarray(getattr(direction, name))
        if not np.isfinite(values).all():
            raise ValueError("tangent must be finite: " + name)
        if name in GEOMETRY and np.iscomplexobj(values):
            raise ValueError("geometry tangent must be real: " + name)

    def check():
        if checkpoint is not None:
            checkpoint()

    check()
    fine = prolong_state(grid, state) if bundle is None else bundle["fine_state"]
    system = active_fine_system(grid, fine) if bundle is None else bundle["fine_system"]
    failure = chart_failure(system, fine)
    if failure is not None:
        raise PositiveChartExit(failure, np.nan, state.copy())
    tangent = prolong_state(grid, direction)
    check()

    Q, r, chi = fine.Q, fine.r, fine.chi
    pQ, pr, pchi = fine.p_Q, fine.p_r, fine.p_chi
    dQ, dr, dchi = tangent.Q, tangent.r, tangent.chi
    dpQ, dpr, dpchi = tangent.p_Q, tangent.p_r, tangent.p_chi
    A = float(system.A)
    Z = float(feedback_Z(A))
    alpha = float(alpha_of(system.C_W))
    Fr, fchi = partial_F(r, A, system.C_W)
    fchi = float(fchi)
    dFr = -8*np.pi*A*dr
    coefficient = 8*np.pi*A/fchi
    pi = pr - Fr*pchi/fchi
    dpi = dpr - (dFr*pchi + Fr*dpchi)/fchi
    F = feedback_F(r, chi, A, system.C_W)
    dF = Fr*dr + fchi*dchi
    V = feedback_V(r, chi, A, system.C_W, system.C_F, system.flux)
    dV = 16*np.pi*A*r*dr - alpha*(4+2*chi)*dchi

    # Differentiate the actual D(F) node operation. Do not use D(F)=F_r D(r)+...
    derivatives = system.derivative @ np.column_stack((Q, F, dr, dQ, dF))
    Qx, Fx, drx, dQx, dFx = derivatives.T
    weight = 2*Qx/Q
    dweight = 2*(dQx/Q - Qx*dQ/Q**2)
    quotient_direction = dFx/Q - Fx*dQ/Q**2
    second = system.derivative @ np.column_stack((weight, dweight, quotient_direction, drx))
    weight_x, dweight_x, quotient_x, drxx = second.T
    check()

    # In the conformal pullback, F_Q+F_L=M*kappa*S1 exactly. The K terms
    # cancel algebraically, not by deleting their effect from the initial data.
    a0, a1 = fine.phi0, fine.phi1
    da0, da1 = tangent.phi0, tangent.phi1
    weights = system.occupations[None, :]
    dS1 = np.sum(weights*(
        da0.conj()*a1 + a0.conj()*da1 + da1.conj()*a0 + a1.conj()*da0
    ), axis=1).real

    q_rate = (dQ*pchi + Q*dpchi)/(2*fchi)
    r_rate = dpi/(2*Z)
    chi_rate = (dQ*pQ + Q*dpQ)/(2*fchi) - (dFr*pi + Fr*dpi)/(2*fchi*Z)
    pq_rate = (
        -(dpQ*pchi + pQ*dpchi)/(2*fchi)
        + 2*dQ*V + 2*Q*dV
        + 2*dFx*Qx/Q**2 + 2*Fx*dQx/Q**2 - 4*Fx*Qx*dQ/Q**3
        + 2*quotient_x - system.multiplicity*system.kappa*dS1/system.dx
    )
    pr_rate = (
        -coefficient*(dpi*pchi + pi*dpchi)/(2*Z)
        + 2*Z*drxx + 16*np.pi*A*(2*Q*dQ*r + Q**2*dr)
        + dFr*weight_x + Fr*dweight_x
    )
    pchi_rate = (
        -alpha*(2*Q*dQ*(4+2*chi) + 2*Q**2*dchi) + fchi*dweight_x
    )

    # H= sigma2 P + kappa Q sigma1; delta(L/Q)=0. There is no M on phi_dot.
    count = a0.shape[1]
    images = system.momentum @ np.column_stack((da0, da1))
    image0, image1 = images[:, :count], images[:, count:]
    phi0_rate = -image1 - 1j*system.kappa*(Q[:, None]*da1 + dQ[:, None]*a1)
    phi1_rate = image0 - 1j*system.kappa*(Q[:, None]*da0 + dQ[:, None]*a0)
    check()

    geometry = pull_geometry(grid, np.column_stack((
        q_rate, r_rate, chi_rate, pq_rate, pr_rate, pchi_rate,
    )))
    columns = grid.U_f.conj().T @ np.column_stack((phi0_rate, phi1_rate))
    result = CauchyState(
        **{name: geometry[:, index] for index, name in enumerate(GEOMETRY)},
        phi0=columns[:, :count], phi1=columns[:, count:],
    )
    if not all(np.isfinite(getattr(result, name)).all() for name in FIELDS):
        raise ValueError("coupled tangent left the finite chart")
    check()
    return result


def _constraint_mixed(system, state, first, second, *, checkpoint=None):
    """Mixed derivative of the actual quadrature (h=Q*C, D) constraints.

    The two directions are held fixed in this Hessian action. Dependence of
    a direction on the state belongs in the separate first derivative term.
    """
    Q, r, chi = state.Q, state.r, state.chi
    pQ, pchi = state.p_Q, state.p_chi
    u, v = first, second
    A = float(system.A)
    Z = float(feedback_Z(A))
    alpha = float(alpha_of(system.C_W))
    Fr, fchi = partial_F(r, A, system.C_W)
    fchi = float(fchi)
    coefficient = 8*np.pi*A/fchi
    pi = state.p_r-Fr*pchi/fchi
    upi = u.p_r+coefficient*(u.r*pchi+r*u.p_chi)
    vpi = v.p_r+coefficient*(v.r*pchi+r*v.p_chi)
    uvpi = coefficient*(u.r*v.p_chi+v.r*u.p_chi)
    V = feedback_V(r, chi, A, system.C_W, system.C_F, system.flux)
    uV = 16*np.pi*A*r*u.r-alpha*(4+2*chi)*u.chi
    vV = 16*np.pi*A*r*v.r-alpha*(4+2*chi)*v.chi
    uvV = 16*np.pi*A*u.r*v.r-2*alpha*u.chi*v.chi
    F = feedback_F(r, chi, A, system.C_W)
    uF, vF = Fr*u.r+fchi*u.chi, Fr*v.r+fchi*v.chi
    uvF = -8*np.pi*A*u.r*v.r
    rows = system.derivative @ np.column_stack((
        F, uF, vF, uvF, u.r, v.r, u.chi, v.chi, u.p_Q, v.p_Q,
    ))
    Fx, uFx, vFx, uvFx, urx, vrx, uchix, vchix, upQx, vpQx = rows.T
    uG, vG = uFx/Q-Fx*u.Q/Q**2, vFx/Q-Fx*v.Q/Q**2
    uvG = uvFx/Q-(uFx*v.Q+vFx*u.Q)/Q**2+2*Fx*u.Q*v.Q/Q**3
    uGx, vGx, uvGx = (system.derivative @ np.column_stack((uG, vG, uvG))).T
    if checkpoint is not None:
        checkpoint()

    h_geometry = (
        (u.Q*(v.p_Q*pchi+pQ*v.p_chi) + v.Q*(u.p_Q*pchi+pQ*u.p_chi)
         + Q*(u.p_Q*v.p_chi+v.p_Q*u.p_chi))/(2*fchi)
        + (upi*vpi+pi*uvpi)/(2*Z) + 2*Z*urx*vrx
        - 2*u.Q*v.Q*V - 2*Q*(u.Q*vV+v.Q*uV) - Q**2*uvV
        - 2*(u.Q*vGx+v.Q*uGx+Q*uvGx)
    )
    d_geometry = (
        u.p_r*vrx+v.p_r*urx+u.p_chi*vchix+v.p_chi*uchix
        - u.Q*vpQx-v.Q*upQx
    )
    count = state.phi0.shape[1]
    images = system.momentum @ np.column_stack((u.phi0, u.phi1, v.phi0, v.phi1))
    Pu0, Pu1, Pv0, Pv1 = (images[:, i*count:(i+1)*count] for i in range(4))
    weights = system.occupations[None, :]

    def real_sum(values):
        return np.sum(weights*values, axis=1).real

    uS = real_sum(u.phi0.conj()*state.phi1+state.phi0.conj()*u.phi1
                 +u.phi1.conj()*state.phi0+state.phi1.conj()*u.phi0)
    vS = real_sum(v.phi0.conj()*state.phi1+state.phi0.conj()*v.phi1
                 +v.phi1.conj()*state.phi0+state.phi1.conj()*v.phi0)
    uvS = real_sum(u.phi0.conj()*v.phi1+v.phi0.conj()*u.phi1
                  +u.phi1.conj()*v.phi0+v.phi1.conj()*u.phi0)
    uvK = real_sum(
        -1j*(u.phi0.conj()*Pv1+v.phi0.conj()*Pu1)
        +1j*(u.phi1.conj()*Pv0+v.phi1.conj()*Pu0))
    uvP = real_sum(u.phi0.conj()*Pv0+v.phi0.conj()*Pu0
                  +u.phi1.conj()*Pv1+v.phi1.conj()*Pu1)
    h = h_geometry+system.multiplicity/system.dx*(
        uvK+system.kappa*(u.Q*vS+v.Q*uS+Q*uvS))
    d = d_geometry-system.multiplicity/system.dx*uvP
    if checkpoint is not None:
        checkpoint()
    return h, d


def conformal_tangent_constraints(grid, state, direction, *, rate=None, bundle=None,
                                 tangent_rate=None, include_vectors=False, checkpoint=None):
    """Constraint variation and variation of the original finite-model forcing.

    delta_f = D2C[u,X] + DC[DX u] - spatial_wave_operator(DC[u]).
    No zero-constraint premise is used. Norms are on the full quadrature
    grid, including unresolved modes; a time/error enclosure is not claimed.
    """
    if grid.gauge != "conformal":
        raise ValueError("constraint tangent transport requires the conformal gauge")
    if rate is None or bundle is None:
        rate, bundle = compose_fine_hamiltonian(grid, state)
    if tangent_rate is None:
        tangent_rate = conformal_jvp(grid, state, direction, bundle=bundle, checkpoint=checkpoint)
    fine, system, source = (bundle[key] for key in ("fine_state", "fine_system", "source"))
    du, velocity, dvelocity = (prolong_state(grid, value) for value in (direction, rate, tangent_rate))

    def fields(value):
        return tuple(getattr(value, name) for name in FIELDS)

    _, _, dh, dd = _constraint_directional(system, fine, fields(du), source)
    _, _, ah, ad = _constraint_directional(system, fine, fields(dvelocity), source)
    hh, hd = _constraint_mixed(system, fine, du, velocity, checkpoint=checkpoint)
    dh_time, dd_time = hh+ah, hd+ad
    forcing_h = dh_time-system.derivative@dd
    forcing_d = dd_time-system.derivative@dh
    derivatives = system.derivative @ np.column_stack((
        fine.r, fine.chi, fine.p_Q, du.r, du.chi, du.p_Q,
    ))
    rx, chix, pQx, urx, uchix, upQx = derivatives.T
    geometric_momentum = (
        du.p_r*rx+fine.p_r*urx+du.p_chi*chix+fine.p_chi*uchix
        -du.Q*pQx-fine.Q*upQx
    )
    shift_source = dd-geometric_momentum
    if not all(np.isfinite(value).all() for value in (dh, dd, dh_time, dd_time)):
        raise ValueError("nonfinite linearized constraint transport")
    norm = lambda first, second: float(np.sqrt(grid.dx_q*np.sum(first**2+second**2)))
    energy_rate = float(2*grid.dx_q*np.sum(dh*dh_time+dd*dd_time))
    forcing_rate = float(2*grid.dx_q*np.sum(dh*forcing_h+dd*forcing_d))
    report = {
        "constraint_tangent_norm": norm(dh, dd),
        "forcing_tangent_norm": norm(forcing_h, forcing_d),
        "hamilton_tangent_max_abs": float(np.max(abs(dh))),
        "momentum_tangent_max_abs": float(np.max(abs(dd))),
        "shift_source_tangent_max_abs": float(np.max(abs(shift_source))),
        "constraint_tangent_energy_rate": energy_rate,
        "forcing_tangent_energy_rate": forcing_rate,
        "sbp_identity_error": abs(energy_rate-forcing_rate),
        "continuum_constraint_certified": False,
        "time_integral_certified": False,
    }
    if include_vectors:
        report.update(h=dh, d=dd, h_dot=dh_time, d_dot=dd_time,
                      forcing_h=forcing_h, forcing_d=forcing_d, shift_source=shift_source)
    if checkpoint is not None:
        checkpoint()
    return report


def variational_rk4_step(grid, state, direction, dt, *, checkpoint=None):
    """Advance the original RK4 map and its real derivative at identical stages."""
    step = float(dt)
    if not np.isfinite(step) or step <= 0:
        raise ValueError("variational RK4 requires a finite positive time step")

    def stage(z, u):
        if checkpoint is not None:
            checkpoint()
        rate, bundle = compose_fine_hamiltonian(grid, z)
        tangent_rate = conformal_jvp(grid, z, u, bundle=bundle, checkpoint=checkpoint)
        return rate, tangent_rate

    k1, j1 = stage(state, direction)
    k2, j2 = stage(_combine(state, k1, .5*step), _combine(direction, j1, .5*step))
    k3, j3 = stage(_combine(state, k2, .5*step), _combine(direction, j2, .5*step))
    k4, j4 = stage(_combine(state, k3, step), _combine(direction, j3, step))

    def average(first, second, third, fourth):
        return CauchyState(**{name: (
            getattr(first, name)+2*getattr(second, name)
            +2*getattr(third, name)+getattr(fourth, name)
        )/6 for name in FIELDS})

    return (_combine(state, average(k1, k2, k3, k4), step),
            _combine(direction, average(j1, j2, j3, j4), step))
