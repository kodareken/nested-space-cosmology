"""Reconstruct the existing quiet walls from clocks; predict a later reception.

Only compact scalar evidence is retained. This is a conditional synthetic
observer experiment, with no P001 evolution or observational calibration.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from recursive_horizons.nsc_atlas_nested_slice import tension_wall
from recursive_horizons.nsc_atlas_nested_spacetime import (
    central_radar_clock, invert_central_radar_clock, radar_cap_input_bounds,
    radar_junction_branches,
    radar_clock_metric,
)


ROOT = Path(__file__).resolve().parents[1]
RECEIVER_INTERVAL = .04
LENGTH_METRES = 10000.
LIGHT_SPEED = 299792458.
INPUT_RELATIVE_BOUND = 1e-4
CASES = (('one_level', .5, .4, 0.), ('nested_inner', .8, .01, 1.))
CLOCK_KEYS = ('radar_distance', 'one_plus_redshift', 'echo_stretch', 'log_redshift_rate')


def inferred_clock(row):
    return invert_central_radar_clock(**{key:row[key] for key in CLOCK_KEYS})


def predict_reception(inferred, interval, order=32):
    """Use only first-event inference and the constant-material-law hypothesis."""
    R0, H = inferred['radius'], inferred['H']
    a, c = inferred['mass_over_tension'], inferred['junction_radius_coefficient']
    nodes, weights = np.polynomial.legendre.leggauss(order)
    count = 0

    def elapsed(R):
        nonlocal count
        radii = (R+R0)/2+(R-R0)/2*nodes
        beta = a/radii**2+c*radii
        f = 1-(H*radii)**2
        squared = beta**2-f
        if np.any(f <= 0) or np.any(squared <= 0) or np.any(beta <= 0):
            raise ValueError('prediction left its admitted future static-cap interval')
        speed = np.sqrt(squared)
        count += order
        return float((R-R0)/2*weights @ ((beta+speed)/(f*speed)))

    # Fixed physical interval, chosen before the application: these two
    # inputs move less than .02 ell during .04 ell of receiver proper time.
    R = float(brentq(lambda radius:elapsed(radius)-interval,R0,R0+.02,
                    xtol=8e-15,rtol=8e-15))
    beta = a/R**2+c*R
    f = 1-(H*R)**2
    speed = math.sqrt(beta*beta-f)
    return dict(radius=R,one_plus_redshift=(beta+speed)/f,
                scalar_evaluations=count)


def independent_wall_history(R0, mass, host_H, interval, *, change_law=0.):
    """Proper-time second-order dynamics and separate outgoing null clock.

    It uses the original signed bulk/junction data, not inferred a,c.
    A counterexample adds change_law*(R-R0)^2 to beta_minus; its first
    event and derivative agree, but its material stress subsequently differs.
    """
    initial = tension_wall(R0,0.,1.,mass,host_H,2.)
    optical0 = math.atanh(R0)
    calls = 0

    def geometry(R):
        nonlocal calls
        state = tension_wall(R,0.,1.,mass,host_H,2.)
        beta0 = state['beta_minus']
        beta_radius = (state['Rddot']-R)/beta0
        beta = beta0+change_law*(R-R0)**2
        derivative = beta_radius+2*change_law*(R-R0)
        calls += 1
        return beta, derivative

    def rhs(tau, y):
        R, velocity, _ = y
        beta, derivative = geometry(R)
        f = 1-R*R
        if f <= 0:
            raise ValueError('proper-time check reached the cap horizon')
        # Normalization supplies T_tau independently of the clock formula.
        return [velocity,beta*derivative+R,math.sqrt(f+velocity**2)/f]

    def arrival(tau, y):
        return y[2]+math.atanh(y[0])-interval
    arrival.terminal = True
    arrival.direction = 1
    v0 = initial['Rdot']
    solution = solve_ivp(rhs,(0.,.08),(R0,v0,-optical0),events=arrival,
                         method='DOP853',rtol=2e-13,atol=2e-14,dense_output=True,
                         max_step=.01)
    if not solution.success or len(solution.t_events[0]) != 1:
        raise RuntimeError('independent receiver-time event was not reached')
    tau = float(solution.t_events[0][0])
    R, velocity, T = map(float, solution.y_events[0][0])
    beta, derivative = geometry(R)
    normalized = math.sqrt(1-R*R+velocity*velocity)
    Z = (normalized+velocity)/(1-R*R)
    # For the counterexample the same vacuum bulks determine a compatible
    # variable wall stress through Israel; no bulk mass flux is inserted.
    original_plus = tension_wall(R,0.,1.,mass,host_H,2.)['beta_plus']
    plus_f = 1-2*mass/R-host_H*host_H*R*R
    plus = math.copysign(math.sqrt(velocity*velocity+plus_f),original_plus)
    surface_sigma = (normalized-plus)/(4*math.pi*R)
    # Independent temporal differences of the propagated clock, about a
    # strict interior proper-time point, check d(log Z)/dt_receive.
    mid = tau/2
    h = 1e-4*tau
    def logZ(t):
        r,v,_ = solution.sol(t)
        return math.log((math.sqrt(1-r*r+v*v)+v)/(1-r*r))
    derivative_tau = (logZ(mid-2*h)-8*logZ(mid-h)+8*logZ(mid+h)-logZ(mid+2*h))/(12*h)
    r,v,_ = map(float,solution.sol(mid))
    bm,bmR = geometry(r)
    zm = (math.sqrt(1-r*r+v*v)+v)/(1-r*r)
    drift_gap = derivative_tau/zm-(bmR/zm+r)
    return dict(radius=R,one_plus_redshift=Z,wall_proper_interval=tau,
        receiver_event_residual=T+math.atanh(R)-interval,
        normalization_residual=normalized-beta,
        temporal_drift_check=drift_gap,positive_wall_surface_density=surface_sigma>0,
        scalar_evaluations=calls)


def calculation():
    started = time.process_time()
    result = {}
    evaluations = 0
    for name,R,mass,host in CASES:
        state = tension_wall(R,0.,1.,mass,host,2.)
        clock = central_radar_clock(radius=R,H=1.,beta_minus=state['beta_minus'],
                                   Rdot=state['Rdot'],Rddot=state['Rddot'])
        inferred = inferred_clock(clock)
        branches = radar_junction_branches(inferred,host_H=host)
        if len(branches)!=1 or max(abs(branches[0]['tension']-2.),
                                  abs(branches[0]['mass']-mass),abs(inferred['H']-1.))>3e-12:
            raise RuntimeError('first-event geometry/junction reconstruction failed')
        query = {key:(clock[key]*(1-INPUT_RELATIVE_BOUND),clock[key]*(1+INPUT_RELATIVE_BOUND))
                 for key in CLOCK_KEYS[:3]}
        bounds = radar_cap_input_bounds(**query)
        predictions = [predict_reception(inferred,RECEIVER_INTERVAL,order) for order in (32,64)]
        direct = independent_wall_history(R,mass,host,RECEIVER_INTERVAL)
        changed = independent_wall_history(R,mass,host,RECEIVER_INTERVAL,change_law=.5)
        fine = predictions[1]
        radius_gap = fine['radius']-direct['radius']
        clock_gap = fine['one_plus_redshift']-direct['one_plus_redshift']
        if max(abs(radius_gap),abs(clock_gap))>2e-10:
            raise RuntimeError('null quadrature prediction and proper-time dynamics disagree')
        if abs(changed['one_plus_redshift']-fine['one_plus_redshift'])<1e-7:
            raise RuntimeError('finite clock prediction did not distinguish the changed wall law')
        evaluations += sum(row['scalar_evaluations'] for row in predictions)+direct['scalar_evaluations']+changed['scalar_evaluations']
        metric_readout=radar_clock_metric(one_plus_redshift=clock['one_plus_redshift'],
                                         echo_stretch=clock['echo_stretch'])
        result[name] = dict(inputs=dict(R_over_ell=R,H_ell=1.,k_ell=2.,M_over_ell=mass,
            host_H_ell=host,scalar_state='quiet'),first_event_clocks=clock,
            model_independent_radial_readout=metric_readout,
            inferred_geometry=inferred,junction_branches=branches,
            assumed_100ppm_clock_box=bounds,
            SI=dict(length_metres=LENGTH_METRES,radar_one_way_time_seconds=clock['radar_distance']*LENGTH_METRES/LIGHT_SPEED,
                log_redshift_rate_per_second=clock['log_redshift_rate']*LIGHT_SPEED/LENGTH_METRES,
                inferred_cap_Lambda_per_metre_squared=inferred['Lambda_cap']/LENGTH_METRES**2,
                assumed_cap_Lambda_range_per_metre_squared=[x/LENGTH_METRES**2 for x in bounds['Lambda_cap']]),
            held_out_reception=dict(receiver_interval_over_ell=RECEIVER_INTERVAL,
                receiver_interval_seconds=RECEIVER_INTERVAL*LENGTH_METRES/LIGHT_SPEED,
                predicted_radius_over_ell=fine['radius'],predicted_one_plus_redshift=fine['one_plus_redshift'],
                independent_proper_time_solution=direct,
                quadrature_64_minus_32={key:fine[key]-predictions[0][key] for key in ('radius','one_plus_redshift')},
                independent_radius_gap=radius_gap,independent_redshift_gap=clock_gap,
                counterexample_redshift_gap=changed['one_plus_redshift']-fine['one_plus_redshift'],
                counterexample_positive_surface_density=changed['positive_wall_surface_density']),
            limits='Synthetic quiet model and ideal clocks; supplied curvature is recovered, not generated. No observed Lambda or total physical error inferred.')
    elapsed = time.process_time()-started
    if elapsed>30 or evaluations>10000:
        raise RuntimeError('finite radar calculation exceeded its compute allowance')
    sources = ('src/recursive_horizons/nsc_atlas_nested_spacetime.py',
        'src/recursive_horizons/nsc_atlas_nested_slice.py',
        'scripts/check_nsc_atlas_radar_inverse.py','tests/test_nsc_atlas_nested_spacetime.py')
    return dict(schema='nsc-atlas-radar-inverse-v1',cases=result,
        counts=dict(scalar_evaluations=evaluations,new_field_grids=0,P001_evolutions=0),
        process_cpu_seconds=elapsed,
        source_sha256={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sources},
        reproduce='python scripts/run_python.py scripts/check_nsc_atlas_radar_inverse.py')


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    record=calculation()
    encoded=json.dumps(record,indent=2,allow_nan=False)+'\n'
    if len(encoded.encode())>24000:
        raise RuntimeError('summary exceeded its compact evidence size')
    if args.output:
        with args.output.open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(encoded)
    print(encoded)
