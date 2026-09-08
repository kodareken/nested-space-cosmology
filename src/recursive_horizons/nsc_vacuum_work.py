"""Prescribed smooth geometry pulse: vacuum response and external work.

r(t,x)=r0(x)/(1+epsilon*s(x)*g(t)), N=q=1 gives H(t)=H0+epsilon*g(t)*V.
The drive is an input experiment, not a solution of the metric equations.
"""
import numpy as np
from math import erfc, pi, sqrt
from scipy.integrate import solve_ivp
from .nsc_shape_response import staggered_operator
from .nsc_energy_transfer import dirac_hamiltonian, vacuum, positions


def pulse_profile(metric):
    return np.exp(2*(np.cos(2*np.pi*(metric.x+.6)/metric.length)-1))


def radius_vertex(metric, kappa=1, eta=.5, profile=None):
    """Exact variation for inverse-radius pulse; no independently fitted link."""
    _, coefficients = staggered_operator(metric, kappa, eta)
    profile = pulse_profile(metric) if profile is None else np.asarray(profile)
    if profile.shape != (metric.points,) or not np.isfinite(profile).all():
        raise ValueError("finite pulse profile matching the spatial grid required")
    dv = metric.lapse*kappa/metric.sphere_radius*profile
    edge = (dv+np.roll(dv,-1))/4
    i = np.arange(metric.points)
    a = np.diag(edge)
    a[i,(i+1)%metric.points] = edge*coefficients[-1]
    zero = np.zeros_like(a)
    return np.block([[zero,a.T],[a,zero]])


def leading_response(h, vertex, width=.5):
    """Coefficients of epsilon² for Gaussian pulse, integrated over all time.

    One radial representative channel; full Dirac multiplier is 4*kappa.
    Energy counts the positive particle AND its negative-sea hole.
    """
    if not np.isfinite(width) or width <= 0:
        raise ValueError("positive finite pulse duration required")
    e,v,_ = vacuum(h)
    positive,negative = e>0,e<0
    rotated = v[:,positive].T.conj()@vertex@v[:,negative]
    frequency = e[positive,None]-e[None,negative]
    weights = abs(rotated)**2*(2*np.pi*width**2)*np.exp(-(frequency*width)**2)
    return {'pair_number_coefficient': float(weights.sum()),
            'energy_coefficient': float(np.sum(frequency*weights)),
            'particle_energy_coefficient': float(np.sum(e[positive,None]*weights)),
            'hole_energy_coefficient': float(np.sum(-e[None,negative]*weights)),
            'minimum_pair_frequency': float(np.min(frequency)),
            'transition_vertex_norm': float(np.linalg.norm(rotated))}


def evolve_vacuum_pulse(metric, epsilon=.02, width=.5, rtol=2e-11):
    """Independent finite real-time occupied-mode evolution plus work ledger.

    Integrate from -7*width to +7*width (Gaussian endpoints nonzero but tiny).
    Subtract the fixed vacuum's vertex expectation in the work integrand;
    its net integrated contribution is exactly zero at equal pulse endpoints.
    """
    if not np.isfinite(epsilon) or epsilon <= -1 or not np.isfinite(width) or width <= 0:
        raise ValueError("pulse must keep radius positive and have positive duration")
    h = dirac_hamiltonian(metric)
    vertex = radius_vertex(metric)
    e,v,_ = vacuum(h)
    w = v.T.conj()@vertex@v
    projection = v.T.conj()@np.diag((positions(metric)>=0).astype(float))@v
    occupied = e<0
    initial = np.eye(len(h),dtype=complex)[:,occupied]
    count = int(occupied.sum())
    shape = initial.shape
    static_v = float(np.trace(w[np.ix_(occupied,occupied)]).real)
    static_va = float(np.trace((projection@w)[np.ix_(occupied,occupied)]).real)
    static_ea = float(np.sum(e[occupied]*np.diag(projection)[occupied]).real)
    def interaction_vertex(time):
        phase = np.exp(1j*e*time)
        return phase[:,None]*w*phase.conj()[None,:]
    def rhs(time, packed):
        modes = packed[:-3].reshape(shape)
        g = np.exp(-.5*(time/width)**2)
        gp = -time/width**2*g
        # Remove free sea phases exactly: the integrator resolves the physical
        # interaction, not an extensive vacuum subtraction polluted by phase error.
        wt = interaction_vertex(time)
        wm = wt@modes
        derivative = -1j*epsilon*g*wm
        power = epsilon*gp*(np.vdot(modes,wm).real-static_v)
        phase = np.exp(1j*e*time)
        pt = phase[:,None]*projection*phase.conj()[None,:]
        pm = pt@modes
        hm = e[:,None]*modes+epsilon*g*wm
        ham = .5*(pt@hm+e[:,None]*pm+epsilon*g*(wt@pm))
        inward_power = -2*np.vdot(hm,ham).imag
        regional_work = epsilon*gp*(np.vdot(pm,wm).real-static_va)
        return np.r_[derivative.ravel(),power,inward_power,regional_work]
    packed = np.r_[initial.ravel(),0j,0j,0j]
    sol = solve_ivp(rhs,(-7*width,7*width),packed,method='DOP853',rtol=rtol,atol=rtol*.05)
    if not sol.success:
        raise RuntimeError(sol.message)
    final = sol.y[:-3,-1].reshape(shape)
    pairs = float(np.sum(abs(final[e>0])**2))
    holes = float(count-np.sum(abs(final[e<0])**2))
    energy = float(np.sum(e[:,None]*abs(final)**2)-np.sum(e[occupied]))
    endpoint = np.exp(-49/2)
    end_correction = epsilon*endpoint*(np.vdot(final,interaction_vertex(7*width)@final).real-static_v)
    work = float(sol.y[-3,-1].real)
    phase = np.exp(1j*e*7*width)
    pt = phase[:,None]*projection*phase.conj()[None,:]
    pm = pt@final
    regional_energy = float(np.vdot(pm,e[:,None]*final).real-static_ea)
    regional_endpoint = epsilon*endpoint*(np.vdot(pm,interaction_vertex(7*width)@final).real-static_va)
    regional_inflow = float(sol.y[-2,-1].real)
    regional_work = float(sol.y[-1,-1].real)
    # Duhamel bound for unitary evolution: ||U_full-U_truncated|| is at
    # most the integral of the omitted perturbation's operator norm.
    tail_bound = abs(epsilon)*np.linalg.norm(vertex,2)*sqrt(2*pi)*width*erfc(7/sqrt(2))
    return {'epsilon': epsilon, 'width': width, 'pair_number': pairs,
            'hole_number': holes, 'excitation_energy': energy,
            'external_work': work, 'endpoint_interaction_correction': float(end_correction),
            'work_balance_residual': float(abs(work-energy-end_correction)),
            'regional_energy_change': regional_energy,
            'integrated_regional_inflow': regional_inflow,
            'regional_geometric_work': regional_work,
            'regional_endpoint_correction': float(regional_endpoint),
            'regional_balance_residual': float(abs(regional_energy+regional_endpoint-regional_inflow-regional_work)),
            'orthonormality_residual': float(np.max(abs(final.T.conj()@final-np.eye(count)))),
            'relative_integration_tolerance': rtol,
            'finite_lattice_omitted_pulse_unitary_bound': float(tail_bound),
            'initial_metric_pulse_magnitude': float(epsilon*endpoint)}
