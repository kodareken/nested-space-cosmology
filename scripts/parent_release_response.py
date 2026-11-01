"""Read-only local release diagnostics on the ORIGINAL measured GH history.

There is no history loader, new Cauchy state, physical ODE or PDE here. Retarded
integrals reconstruct the supplied history and separate its measured forces.
Conditional fixed-field/fixed-action equilibria are comparisons, not exits.
"""
from functools import lru_cache
import math

import numpy as np
from scipy.optimize import brentq

from recursive_horizons.nsc_atlas_parent_scale_selector import VolumeParent


K = 84.
SQRT_K = math.sqrt(K)
D = 12.


@lru_cache(maxsize=1)
def _reference():
    model = VolumeParent()
    minimum = next(row for row in model.vacuum_roots() if row['kind'] == 'minimum')
    maximum = next(row for row in model.vacuum_roots() if row['kind'] == 'maximum')
    fold_x = brentq(lambda x: 225*x**6-560*x+336, .99, 1., xtol=5e-15)
    fold_s = -.5*math.log(fold_x)
    loading = lambda s: -14*math.exp(-6*s)+17.5*math.exp(-8*s)-3.125*math.exp(-18*s)
    return dict(sstar=minimum['s'], bare_top=maximum['s'],
        bare_top_V=maximum['vacuum'], bare_minimum_V=minimum['vacuum'],
        fold_s=fold_s, loading_critical=loading(fold_s),
        bare_minimum_x=math.exp(-2*minimum['s']),
        bare_maximum_x=math.exp(-2*maximum['s']), fold_x=fold_x)


def _potential(s, order=0):
    s = np.asarray(s, dtype=float)
    return (7*(-12)**order*np.exp(-12*s)
        - 7.5*(-14)**order*np.exp(-14*s)
        + (25/32)*(-24)**order*np.exp(-24*s))


def potential_hessian(s, chi):
    """Exact canonical (S,chi) Hessian, S=sqrt(84)*(s-sstar).

    W=V(s)+mu²(s)chi²/2; returned shape is broadcast(s,chi).shape+(2,2).
    A local eigenvalue is neither a finite-time growth exponent nor a seed.
    """
    s, chi = np.broadcast_arrays(np.asarray(s, dtype=float), np.asarray(chi, dtype=float))
    if not (np.all(np.isfinite(s)) and np.all(np.isfinite(chi))):
        raise ValueError('finite log volume and scalar field required')
    mu2 = 400*np.exp(-12*s)
    matrix = np.empty(s.shape+(2, 2))
    matrix[..., 0, 0] = (_potential(s, 2)+72*mu2*chi*chi)/K
    matrix[..., 0, 1] = matrix[..., 1, 0] = -12*mu2*chi/SQRT_K
    matrix[..., 1, 1] = mu2
    return matrix


def instantaneous_release(s, chi, PiChi, gradient_chi_squared=None):
    """Scalar instantaneous comparisons at the actual measured event.

    PiChi is the Eulerian normal derivative. The oscillator-action comparison
    omits spatial transport, expansion and fast phase; gradient energy is
    exposed separately. No local loading threshold is an escape theorem.
    """
    s, chi, PiChi = map(float, (s, chi, PiChi))
    grad2 = None if gradient_chi_squared is None else float(gradient_chi_squared)
    if (not all(map(math.isfinite, (s, chi, PiChi)))
            or (grad2 is not None and (not math.isfinite(grad2) or grad2 < 0))):
        raise ValueError('finite fields/momentum and nonnegative gradient norm required')
    ref = _reference()
    mu2 = 400*math.exp(-12*s)
    rho = (PiChi*PiChi+mu2*chi*chi)/2
    loading = rho*math.exp(6*s)
    A = 7+200*chi*chi
    xturn = (14/15)**.2
    critical_A = (175/24)*xturn
    current_fixed_chi = A*math.exp(-12*s)-7.5*math.exp(-14*s)+(25/32)*math.exp(-24*s)
    fixed_chi = dict(A=A, fold_A=critical_A, extrema_exist=A < critical_A,
        status='fold_or_no_extrema', minimum_s=None, top_s=None, barrier_height=None,
        potential_at_event=current_fixed_chi, potential_margin_to_top=None)
    if A < critical_A:
        polynomial = lambda x: 6*A-52.5*x+9.375*x**6
        upper = 2*xturn
        while polynomial(upper) <= 0:
            upper *= 2
        xmax = brentq(polynomial, 0., xturn, xtol=5e-15)
        xmin = brentq(polynomial, xturn, upper, xtol=5e-15)
        top, bottom = -.5*math.log(xmax), -.5*math.log(xmin)
        value = lambda q: A*math.exp(-12*q)-7.5*math.exp(-14*q)+(25/32)*math.exp(-24*q)
        fixed_chi.update(status='two_extrema', minimum_s=bottom, top_s=top,
            barrier_height=value(top)-value(bottom),
            potential_margin_to_top=value(top)-current_fixed_chi)
    fixed_action = dict(loading=loading, critical_loading=ref['loading_critical'],
        loading_ratio=loading/ref['loading_critical'], extrema_exist=loading < ref['loading_critical'],
        status='fold_or_no_extrema', minimum_s=None, top_s=None, barrier_height=None,
        potential_at_event=float(_potential(s))+rho, potential_margin_to_top=None)
    if loading < ref['loading_critical']:
        polynomial = lambda x: x**3*(-14+17.5*x-3.125*x**6)-loading
        if loading == 0:
            xmax, xmin = ref['bare_maximum_x'], ref['bare_minimum_x']
        else:
            xmax = brentq(polynomial, 0., ref['fold_x'], xtol=5e-15)
            xmin = brentq(polynomial, ref['fold_x'], 2*ref['bare_minimum_x'], xtol=5e-15)
        top, bottom = -.5*math.log(xmax), -.5*math.log(xmin)
        value = lambda q: float(_potential(q))+loading*math.exp(-6*q)
        fixed_action.update(status='two_extrema', minimum_s=bottom, top_s=top,
            barrier_height=value(top)-value(bottom),
            potential_margin_to_top=value(top)-fixed_action['potential_at_event'])
    return dict(s=s, chi=chi, mu_squared=mu2, oscillator_density=rho,
        measured_gradient_density=grad2/2 if grad2 is not None else None,
        gradient_to_oscillator_density=grad2/(2*rho) if rho and grad2 is not None else None,
        instantaneous_to_averaged_force_ratio=mu2*chi*chi/rho if rho else None,
        fixed_chi=fixed_chi, fixed_action=fixed_action,
        bare=dict(top_s=ref['bare_top'], potential_margin_to_top=ref['bare_top_V']-float(_potential(s))),
        scope='Conditional local equilibria/loading, not a prepared state or a dynamical release prediction.')


def _retarded_integrals(t, force):
    """Exact first/double integrals of each piecewise-linear sampled force."""
    f = np.asarray(force, dtype=float)
    first, second = np.zeros_like(f), np.zeros_like(f)
    for i, h in enumerate(np.diff(t)):
        second[i+1] = second[i]+h*first[i]+h*h*(2*f[i]+f[i+1])/6
        first[i+1] = first[i]+h*(f[i]+f[i+1])/2
    return first, second


def _work_integrals(t, forces, velocity):
    """Integrate products of linear samples, splitting their actual signs."""
    f = np.asarray(forces, dtype=float)
    positive, negative = np.zeros_like(f), np.zeros_like(f)
    for i, h in enumerate(np.diff(t)):
        v0, dv = velocity[i], velocity[i+1]-velocity[i]
        for j in range(f.shape[1]):
            f0, df = f[i, j], f[i+1, j]-f[i, j]
            coefficients = (f0*v0, f0*dv+df*v0, df*dv)
            primitive = lambda x: h*(coefficients[0]*x+coefficients[1]*x*x/2+coefficients[2]*x**3/3)
            breaks = [0., 1.]
            for start, change in ((f0, df), (v0, dv)):
                if change and 0 < -start/change < 1:
                    breaks.append(-start/change)
            pos, neg = 0., 0.
            for left, right in zip(sorted(set(breaks))[:-1], sorted(set(breaks))[1:]):
                value = primitive(right)-primitive(left)
                pos += max(0., value)
                neg += min(0., value)
            positive[i+1, j] = positive[i, j]+pos
            negative[i+1, j] = negative[i, j]+neg
    return positive+negative, positive, negative


def _crossings(t, values, direction=0):
    events = []
    for i in range(len(t)-1):
        left, right = values[i], values[i+1]
        if left == right or not (left <= 0 < right or left >= 0 > right):
            continue
        sense = 1 if right > left else -1
        if direction and sense != direction:
            continue
        fraction = -left/(right-left)
        event = dict(tau=float(t[i]+fraction*(t[i+1]-t[i])), left_index=i,
                     right_index=i+1, direction=sense)
        if not events or event['tau'] != events[-1]['tau']:
            events.append(event)
    return events


def _positive_duration(t, values):
    duration = 0.
    for i, h in enumerate(np.diff(t)):
        a, b = values[i], values[i+1]
        if a > 0 and b > 0:
            duration += h
        elif (a > 0) != (b > 0):
            fraction = -a/(b-a)
            duration += h*(fraction if a > 0 else 1-fraction)
    return float(duration)


def _carrier_cycles(t, chi, velocity, mu2):
    """Measured positive maxima bound cycles; no fitted damping exponent."""
    maxima = [event for event in _crossings(t, velocity, -1)
              if np.interp(event['tau'], t, chi) > 0]
    energy = (velocity*velocity+mu2*chi*chi)/2
    amplitude = np.sqrt(chi*chi+velocity*velocity/mu2)
    cycles = []
    for first, last in zip(maxima[:-1], maxima[1:]):
        left, right = first['tau'], last['tau']
        grid = np.r_[left, t[(t > left) & (t < right)], right]
        mean = lambda values: float(np.trapezoid(np.interp(grid, t, values), grid)/(right-left))
        cycles.append(dict(tau_start=left, tau_end=right,
            mean_observer_oscillator_energy=mean(energy),
            mean_observer_quadrature_amplitude=mean(amplitude)))
    change = None
    if len(cycles) >= 2:
        change = {key: cycles[-1][key]-cycles[0][key] for key in (
            'mean_observer_oscillator_energy', 'mean_observer_quadrature_amplitude')}
    return dict(positive_maxima=maxima, complete_cycles=cycles,
        first_to_last_complete_cycle_change=change,
        scope='Actual local-field maxima and observer quadratures; no fitted decay law, no conserved mode action.')


def causal_release_response(samples, *, sstar=None):
    """Reconstruct the supplied actual history and its mechanical work ledger.

    Required equal-length arrays: tau,S,chi,uS,uChi,nS,nChi,HessS_uu,
    HessS_nn,spatial_lapS,lapse_accel_dot_DS,traceK. Here S is the GH canonical
    displacement sqrt84*(s-sstar), u is the tracked unit geodesic and n is
    the actual Eulerian normal. All Hessians are physical covariant Hessians.
    Optional gradient_chi_squared is measured using h^-1, not flat gradients.
    Chi KG diagnostics additionally accept HessChi_uu,HessChi_nn,
    spatial_lapChi,lapse_accel_dot_DChi. No initial seed is supplied or selected.
    """
    required = ('tau', 'S', 'chi', 'uS', 'uChi', 'nS', 'nChi', 'HessS_uu',
        'HessS_nn', 'spatial_lapS', 'lapse_accel_dot_DS', 'traceK')
    values = {key: np.asarray(samples[key], dtype=float) for key in required}
    t = values['tau']
    if (t.ndim != 1 or len(t) < 2 or np.any(np.diff(t) <= 0)
            or any(v.shape != t.shape or not np.all(np.isfinite(v)) for v in values.values())):
        raise ValueError('finite aligned one-dimensional samples and increasing proper times required')
    sstar = _reference()['sstar'] if sstar is None else float(sstar)
    if not math.isfinite(sstar):
        raise ValueError('finite reference minimum required')
    S, chi, uS = values['S'], values['chi'], values['uS']
    s = sstar+S/SQRT_K
    mu2 = 400*np.exp(-12*s)
    terms = dict(spatial_laplacian=values['spatial_lapS'],
        slice_expansion=values['traceK']*values['nS'],
        lapse_acceleration=values['lapse_accel_dot_DS'],
        observer_conversion=values['HessS_uu']-values['HessS_nn']-values['lapse_accel_dot_DS'],
        bare_potential=-_potential(s, 1)/SQRT_K,
        carrier_force=6*mu2*chi*chi/SQRT_K)
    labels = tuple(terms)
    forces = np.column_stack(list(terms.values()))
    total = forces.sum(axis=1)
    first, second = _retarded_integrals(t, forces)
    dt = t-t[0]
    reconstructed_S = S[0]+uS[0]*dt+second.sum(axis=1)
    reconstructed_uS = uS[0]+first.sum(axis=1)
    work, positive, negative = _work_integrals(t, forces, uS)
    transport_indices = [labels.index(key) for key in (
        'spatial_laplacian', 'slice_expansion', 'lapse_acceleration', 'observer_conversion')]
    source_index = labels.index('carrier_force')
    mechanical = uS*uS/2+_potential(s)
    transfer = work[:, transport_indices].sum(axis=1)+work[:, source_index]
    grad2 = samples.get('gradient_chi_squared')
    if grad2 is None:
        grad2 = [None]*len(t)
    else:
        grad2 = np.asarray(grad2, dtype=float)
        if grad2.shape != t.shape or not np.all(np.isfinite(grad2)) or np.any(grad2 < 0):
            raise ValueError('aligned nonnegative physical carrier gradient norm required')
    instantaneous = [instantaneous_release(*args) for args in zip(s, chi, values['nChi'], grad2)]
    loading_ratio = np.array([row['fixed_action']['loading_ratio'] for row in instantaneous])
    ref = _reference()
    crossings_top = _crossings(t, s-ref['bare_top'], 1)
    crossings_fold = _crossings(t, s-ref['fold_s'], 1)
    upward = lambda event: np.interp(event['tau'], t, uS) > 0
    crossings_top = [event for event in crossings_top if upward(event)]
    crossings_fold = [event for event in crossings_fold if upward(event)]
    loading_up = _crossings(t, loading_ratio-1, 1)
    first_loaded = float(t[0]) if loading_ratio[0] > 1 else (loading_up[0]['tau'] if loading_up else None)
    lag = lambda events: events[0]['tau']-first_loaded if events and first_loaded is not None else None
    series = dict(tau=t.tolist(), s=s.tolist(), S=S.tolist(), chi=chi.tolist(),
        uS=uS.tolist(), uChi=values['uChi'].tolist(),
        canonical_potential_hessian=potential_hessian(s, chi).tolist(),
        canonical_hessian_eigenvalues=np.linalg.eigvalsh(potential_hessian(s, chi)).tolist(),
        force_labels=list(labels), measured_forces=forces.tolist(),
        kg_normal_residual=(values['HessS_nn']-(total-terms['observer_conversion']-terms['lapse_acceleration'])).tolist(),
        kg_observer_residual=(values['HessS_uu']-total).tolist(),
        force_velocity_integrals=first.tolist(), force_displacement_integrals=second.tolist(),
        reconstructed_S=reconstructed_S.tolist(), reconstructed_uS=reconstructed_uS.tolist(),
        reconstruction_S_error=(reconstructed_S-S).tolist(),
        reconstruction_uS_error=(reconstructed_uS-uS).tolist(),
        cumulative_work=work.tolist(), cumulative_positive_work=positive.tolist(),
        cumulative_negative_work=negative.tolist(),
        mechanical_radion_energy=mechanical.tolist(),
        mechanical_work_residual=(mechanical-mechanical[0]-transfer).tolist(),
        bare_top_distance=(ref['bare_top']-s).tolist(),
        bare_mechanical_margin=(ref['bare_top_V']-mechanical).tolist(),
        instantaneous=instantaneous)
    if 't' in samples:
        coord = np.asarray(samples['t'], dtype=float)
        if coord.shape != t.shape or not np.all(np.isfinite(coord)):
            raise ValueError('aligned finite coordinate times required')
        series['coordinate_t'] = coord.tolist()
    chi_keys = ('HessChi_uu', 'HessChi_nn', 'spatial_lapChi', 'lapse_accel_dot_DChi')
    if any(key in samples for key in chi_keys):
        if not all(key in samples for key in chi_keys):
            raise ValueError('provide all carrier Hessian/KG terms together')
        ck = {key: np.asarray(samples[key], dtype=float) for key in chi_keys}
        if any(v.shape != t.shape or not np.all(np.isfinite(v)) for v in ck.values()):
            raise ValueError('aligned finite carrier Hessian/KG terms required')
        normal = ck['spatial_lapChi']+values['traceK']*values['nChi']-mu2*chi
        series['carrier_kg_normal_residual'] = (ck['HessChi_nn']-normal).tolist()
        series['carrier_kg_observer_residual'] = (
            ck['HessChi_uu']-normal-(ck['HessChi_uu']-ck['HessChi_nn'])).tolist()
    summary = dict(proper_interval=[float(t[0]), float(t[-1])], sample_count=len(t),
        maximum_proper_step=float(max(np.diff(t))),
        maximum_abs_kg_normal_residual=float(max(abs(np.asarray(series['kg_normal_residual'])))),
        maximum_abs_kg_observer_residual=float(max(abs(np.asarray(series['kg_observer_residual'])))),
        maximum_abs_S_reconstruction_error=float(max(abs(reconstructed_S-S))),
        maximum_abs_uS_reconstruction_error=float(max(abs(reconstructed_uS-uS))),
        maximum_abs_mechanical_work_residual=float(max(abs(mechanical-mechanical[0]-transfer))),
        final_work={key:dict(net=float(work[-1, i]), positive=float(positive[-1, i]),
                            negative=float(negative[-1, i])) for i, key in enumerate(labels)},
        duration_loading_above_adiabatic_fold=_positive_duration(t, loading_ratio-1),
        first_loading_above_adiabatic_fold=first_loaded,
        observed_upward_bare_top_crossings=crossings_top,
        observed_upward_adiabatic_fold_s_crossings=crossings_fold,
        duration_from_first_overloading_to_bare_top=lag(crossings_top),
        duration_from_first_overloading_to_fold_s=lag(crossings_fold),
        minimum_bare_top_distance=float(min(ref['bare_top']-s)),
        minimum_bare_mechanical_margin=float(min(ref['bare_top_V']-mechanical)),
        minimum_existing_fixed_chi_potential_margin=min(
            (r['fixed_chi']['potential_margin_to_top'] for r in instantaneous if r['fixed_chi']['extrema_exist']), default=None),
        minimum_existing_fixed_chi_top_distance=min(
            (r['fixed_chi']['top_s']-r['s'] for r in instantaneous if r['fixed_chi']['extrema_exist']), default=None),
        minimum_existing_fixed_action_top_distance=min(
            (r['fixed_action']['top_s']-r['s'] for r in instantaneous if r['fixed_action']['extrema_exist']), default=None),
        observed_uS_turns=[event for event in _crossings(t, uS) if event['tau'] > t[0]],
        carrier_cycles=_carrier_cycles(t, chi, values['uChi'], mu2),
        release_time_prediction=None, escape_proved=False, child_formation_proved=False)
    return dict(normalizations=dict(K=K, d=D, omega0=20., sstar=sstar,
        canonical_field='S=sqrt84*(s-sstar)', kappa=1., units='Existing model reference units'),
        summary=summary, series=series,
        scope='Actual-history local identity and finite-sampling work/reconstruction diagnostics only. '
            'Retarded integrals are exact for linearly interpolated measured forces; compare time refinements. '
            'Observer conversion is Hess_uu-Hess_nn-a.DS, converting n(nS) to the geodesic derivative. '
            'The two KG residuals share this measured conversion and are not independent equations. '
            'Mechanical radion energy omits spatial gradient energy and is not a cap energy or physical flux. '
            'Conditional equilibria omit transport/expansion and do not select a preparation, release time, escape or child.')
