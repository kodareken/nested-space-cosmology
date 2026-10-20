"""Continue the existing two quiet walls in one host and locate incoming rays.

The finite application shares one scalar/CPU allowance across the 32/64 rules,
independent proper-time evolution and separate null-ray integrations. No
field grids or sample histories are written.
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

from recursive_horizons.nsc_atlas_nested_slice import tension_wall
from recursive_horizons.nsc_atlas_nested_spacetime import (
    TwoWallHistoryBudget, TwoWallHostHistory,
)


ROOT = Path(__file__).resolve().parents[1]
PROPER_TIMES = (0., .02, .04)
SCHEMA = 'nsc-atlas-two-wall-history-v1'


def independent_proper_history(name, initial, budget):
    """Second-order wall motion and the normalized signed host Killing clock."""
    inner = name == 'inner'
    args = (0., 1., .01, 1., 2.) if inner else (.01, 1., .4, 0., 2.)
    state = tension_wall(initial['radius'], *args)
    budget.add()
    sign = -1 if inner else 1

    def rhs(tau, y):
        R, speed, _ = y
        budget.add()
        wall = tension_wall(R, *args)
        f = 1-.02/R-R*R
        beta = sign*math.sqrt(f+speed*speed)
        return [speed, wall['Rddot'], beta/f]

    solution = solve_ivp(rhs, (0., .04),
        (initial['radius'], state['Rdot'], initial['T']), t_eval=PROPER_TIMES,
        method='DOP853', rtol=2e-13, atol=2e-14, max_step=.01)
    if not solution.success or solution.y.shape != (3, len(PROPER_TIMES)):
        raise RuntimeError('independent proper-time continuation failed')
    return solution.y.T


def independent_null_arrival(name, origin, event, budget):
    """Integrate the null equation in the independent static overlap."""
    direction = -1 if name == 'inner' else 1
    R0, R1 = origin['radius'], event['radius']

    def rhs(radius, _):
        budget.add()
        return [direction/(1-.02/radius-radius*radius)]

    solution = solve_ivp(rhs, (R0, R1), [origin['initial_event']['T']],
        method='DOP853', rtol=2e-13, atol=2e-14, max_step=.01)
    if not solution.success:
        raise RuntimeError('independent null propagation failed')
    return float(solution.y[0, -1])-event['T']


def _maximum(values):
    return max((abs(float(value)) for value in values), default=0.)


def calculation():
    budget = TwoWallHistoryBudget()
    histories = {}
    for order in (32, 64):
        model = TwoWallHostHistory(order=order, budget=budget)
        cases = {}
        for name in ('inner', 'outer'):
            events = [model.wall_event(name, tau) for tau in PROPER_TIMES]
            origins = [model.characteristic_origin(name, event) for event in events]
            cases[name] = (events, origins)
        histories[order] = (model, cases)
    fine, fine_cases = histories[64]
    coarse_cases = histories[32][1]
    result = {}
    all_differences = []
    for name, (events, origins) in fine_cases.items():
        coarse_events, coarse_origins = coarse_cases[name]
        direct = independent_proper_history(name, fine.initial_wall(name), budget)
        differences = dict(
            radius_64_minus_32=_maximum(a['radius']-b['radius']
                                        for a, b in zip(events, coarse_events)),
            shared_T_64_minus_32=_maximum(a['T']-b['T']
                                          for a, b in zip(events, coarse_events)),
            source_radius_64_minus_32=_maximum(a['radius']-b['radius']
                                              for a, b in zip(origins, coarse_origins)),
            source_label_64_minus_32=_maximum(a['proper_label']-b['proper_label']
                                             for a, b in zip(origins, coarse_origins)),
            arrival_derivative_64_minus_32=_maximum(
                a['arrival_d_tau_d_label']-b['arrival_d_tau_d_label']
                for a, b in zip(origins, coarse_origins)),
            radius_vs_second_order=_maximum(row['radius']-check[0]
                                             for row, check in zip(events, direct)),
            velocity_vs_second_order=_maximum(row['Rdot']-check[1]
                                               for row, check in zip(events, direct)),
            shared_T_vs_second_order=_maximum(row['T']-check[2]
                                              for row, check in zip(events, direct)),
            null_T_vs_independent_ODE=_maximum(
                independent_null_arrival(name, origin, event, budget)
                for origin, event in zip(origins, events)))
        all_differences.extend(differences.values())
        normalization = _maximum(row['lake']['normalization_residual'] for row in events)
        junction = _maximum(value for row in events for value in
                            (*row['junction_residuals'].values(),
                             row['angular_jump_residual'], row['temporal_jump_residual']))
        null_match = _maximum(row['null_matching_residual'] for row in origins)
        if max(normalization, junction, null_match) > 2e-10:
            raise RuntimeError('worldtube or null matching failed')
        if not all(row['lake']['u_tau'] < 0 and row['lake']['w_tau'] > 0 for row in events):
            raise RuntimeError('common host time orientation failed')
        sign = -1 if name == 'inner' else 1
        if not all(sign*row['lake']['u'] > 0 and sign*row['lake']['w'] > 0 for row in events):
            raise RuntimeError('worldtube left the admitted signed host patch')
        initial, endpoint = events[0], events[-1]
        labels = [row['proper_label'] for row in origins]
        radius_bound = (.8*math.exp(math.sqrt(2)*.04) if name == 'inner'
                        else 1.5*math.exp(1.25*.04))
        lower_source = .7 if name == 'inner' else 1.3
        lower_null = fine.initial_event(name, lower_source)['null_coordinate']
        target = endpoint['null_coordinate']
        bracket_product = (lower_null-target)*(initial['null_coordinate']-target)
        if not bracket_product < 0:
            raise RuntimeError('endpoint characteristic lacks a strict original-slice bracket')
        # These algebraic inequalities hold on the full radius intervals.
        # Positive motion preserves R>=R0; Gronwall's elementary exponential
        # bound then keeps the worldtubes in the declared signed patches.
        certificate = dict(
            analytic_radius_upper_bound=radius_bound,
            speed_squared_positive_lower_bound=2*.8**2-1 if name == 'inner'
                else (25/16)*1.5**2-1,
            subtracted_growth_remainder_upper_bound=.01/.8+.000025/.8**4
                if name == 'inner' else .3125/1.5+.038025/1.5**4,
            growth_inequality='Rdot^2=2R^2-1+.01/R+.000025/R^4<2R^2'
                if name == 'inner' else
                'Rdot^2=25R^2/16-1+.3125/R+.038025/R^4<(5R/4)^2',
            source_endpoint_bracket_product=bracket_product,
            source_monotonicity='p_R=1/[Rprime*(Rprime-V)]>0; p_tau<0'
                if name == 'inner' else 'q_R=-1/[Rprime*(Rprime+V)]<0; q_tau>0',
            patch_argument='R>=R0>rb; finite q preserves sign(u); w=(R-rb)/(kappa*u) has the same sign. Thus signs hold on the entire proper-time interval.',
            time_orientation='q_tau and u have the same sign, so u_tau<0. F>0 and proper normalization then imply w_tau>0.')
        if name == 'inner':
            certificate['host_f_positive_lower_bound'] = float(fine.lapse(radius_bound))
            certificate['host_beta_negative_reason'] = 'beta_host=-R+.005/R^2<0 for R>=.8'
            if not certificate['host_f_positive_lower_bound'] > 0:
                raise RuntimeError('inner analytic domain reaches the cosmological horizon')
        else:
            certificate['host_beta_positive_reason'] = 'beta_host=.195/R^2+.75R>0'
        if (certificate['speed_squared_positive_lower_bound'] <= 0
                or certificate['subtracted_growth_remainder_upper_bound'] >= 1):
            raise RuntimeError('analytic whole-interval growth certificate failed')
        result[name] = dict(
            host_null_family=endpoint['null_kind'],
            initial_event=initial,
            endpoint=endpoint,
            source_radius_interval=sorted([origins[0]['radius'], origins[-1]['radius']]),
            source_proper_label_interval=sorted([labels[0], labels[-1]]),
            arrival_derivatives_at_initial_midpoint_endpoint=[
                row['arrival_d_tau_d_label'] for row in origins],
            endpoint_source=origins[-1],
            comparisons=differences,
            residuals=dict(wall_normalization=normalization, junction=junction,
                           characteristic_matching=null_match),
            whole_interval_certificate=certificate,
            domain=dict(wall_radius_bracket=[.8,.85] if name == 'inner' else [1.5,1.6],
                original_source_bracket=[.7,.8] if name == 'inner' else [1.3,1.5],
                proper_time_interval=[0.,.04],
                continuation='Positive expanding root throughout the declared radius bracket.'))
    if max(all_differences) > 2e-9:
        raise RuntimeError('independent finite-history comparison exceeds 2e-9')
    # These initial checks also cover the neck and the actual cosmological
    # horizon, beyond the small source collars used by the finite rays.
    bridge = fine.bridge_event(0.)
    horizon = fine.initial_event('outer', fine.rc)
    initial_residual = _maximum((
        bridge['spatial_norm_residual'], bridge['normal_radius_residual'],
        horizon['lake']['spatial_norm_residual'], horizon['lake']['normal_norm_residual'],
        horizon['lake']['orthogonality_residual'], horizon['lake']['normal_radius_residual']))
    if initial_residual > 2e-10:
        raise RuntimeError('shared initial embedding failed at neck or horizon')
    budget.add(0)
    sources = ('src/recursive_horizons/nsc_atlas_nested_spacetime.py',
        'src/recursive_horizons/nsc_atlas_nested_slice.py',
        'scripts/check_nsc_atlas_two_wall_history.py',
        'tests/test_nsc_atlas_two_wall_history.py')
    return dict(schema=SCHEMA,
        inputs=dict(masses=[0.,.01,.4], H_values=[1.,1.,0.], common_tension=2.,
            initial_wall_radii=[.8,1.5], scalar='zero',
            gauge='Bridge T=0; rstar(.5 ell)=0; Lake L0=ell. Original slice fixes both wall origins.'),
        units='Lengths and proper times in ell, velocities in c=1 units.',
        cases=result,
        common_host=dict(black_hole_radius=fine.rb, cosmological_radius=fine.rc,
            surface_gravity=fine.kappa, outer_initial_proper_label=fine.outer_initial_label(),
            neck=bridge, initial_cosmological_horizon=horizon,
            initial_embedding_max_residual=initial_residual),
        causal_access=dict(inner_segment_to_outer_segment=False,
            outer_segment_to_inner_segment=False,
            reason='On the common retained host chart, future causal curves have du<=0 and dw>=0; F>=0 for R>0. The inner segment has u,w<0 and outer has u,w>0.',
            scope='Within the constructed retained host patch only; no unbuilt global identifications or formation history are assumed.'),
        counts=dict(scalar_geometry_or_RHS_evaluations=budget.count, new_field_grids=0,
                    counting_rule='Quadrature nodes, scalar wall queries, explicit initial/frame and Lake-chart queries, Kottler root-lapse evaluations and independent ODE RHS calls; constituent arithmetic and transcendental primitives are not counted separately.',
                    quadrature_orders=[32,64], wall_proper_time_cuts=list(PROPER_TIMES)),
        process_cpu_seconds=time.process_time()-budget.started,
        source_sha256={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sources},
        reproduce='python scripts/run_python.py scripts/check_nsc_atlas_two_wall_history.py',
        limits='Conditional quiet two-level state. Proper wall clocks are not a simultaneity slicing; no evolved volume average, new field state or generated vacuum curvature.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = calculation()
    encoded = json.dumps(result, indent=2, allow_nan=False)+'\n'
    if len(encoded.encode()) > 24000:
        raise RuntimeError('two-wall evidence exceeds compact 24KB allowance')
    if args.output:
        with args.output.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(encoded)
    print(encoded)


if __name__ == '__main__':
    main()
