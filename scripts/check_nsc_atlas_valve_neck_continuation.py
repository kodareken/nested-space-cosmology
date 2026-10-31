"""Bounded continuation diagnostic; no source or preparation changes."""
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
from scipy.integrate import DOP853
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_valve_spherical import SphericalValve


def diagnostic(model, t, state, band=False):
    row = model.diagnostic(t, state)
    neck = model.pieces(state, [0.])
    q, chi, ps, pc = neck['y'][[4, 5, 6, 7], 0]
    row.update(neck_S=float(q), neck_chi=float(chi), neck_Pi_S=float(ps),
               neck_Pi_chi=float(pc), neck_W=float(neck['W'][0]),
               neck_scalar_gradient_squared=float(neck['grad2'][0]),
               neck_normal_expansion=float(-neck['K'][0]),
               neck_radius_rate=float(-neck['kb'][0] * neck['b'][0]))
    kb0, K0 = float(neck['kb'][0]), float(neck['K'][0])
    kb_t = float(neck['Rt'][0])+K0*kb0-model.kappa*float(neck['W'][0])
    row.update(neck_Rll_over_R_H2=float(-neck['Rr'][0]/(2*model.H**2)),
               neck_radius_acceleration_over_radius_H2=(kb0*kb0-kb_t)/model.H**2,
               neck_theta_time_derivative_over_H2=-2*kb_t/model.H**2,
               neck_null_positive_flux=float((ps*ps+pc*pc)+neck['grad2'][0]))
    # Source-inferred curvature invariants, only meaningful while constraints hold.
    u = np.linspace(-.99, .99, 301)
    z = model.pieces(state, u)
    PS, PC = z['y'][[6, 7]]
    Su, Cu = z['yu'][[4, 5]]
    gS = np.sqrt(1-u*u)/z['q']*Su
    gC = np.sqrt(1-u*u)/z['q']*Cu
    pi2 = PS*PS+PC*PC
    kinetic_trace = -pi2 + z['grad2']
    ricci_scalar = model.kappa*(kinetic_trace+4*z['W'])
    Mss, Mcc, Msc = -PS*PS+gS*gS, -PC*PC+gC*gC, -PS*PC+gS*gC
    ricci_square = model.kappa**2*(Mss*Mss+Mcc*Mcc+2*Msc*Msc+
                                    2*z['W']*kinetic_trace+4*z['W']**2)
    R = z['radius']; gamma = z['rr']/z['q']; U = -z['kb']*R
    mass_over_R3 = (1+U*U-gamma*gamma)/(2*R*R)
    Weyl_E = mass_over_R3-model.kappa*(pi2-z['grad2']+2*z['W'])/12
    kretsch = 48*Weyl_E*Weyl_E+2*ricci_square-ricci_scalar*ricci_scalar/3
    row.update(sampled_max_abs_Ricci_scalar_over_H2=float(max(abs(ricci_scalar))/model.H**2),
               sampled_max_abs_source_inferred_Kretschmann_over_H4=float(max(abs(kretsch))/model.H**4),
               sampled_max_spatial_sectional_curvature_over_H2=float(max(max(abs(z['Rr']/2)),max(abs(z['Rt']-z['Rr']/2))))/model.H**2,
               curvature_scope='Einstein-source inferred spacetime invariants; sampled |u|<=.99, conditional on constraints')
    if band:
        row['trapped_band'] = model.trapped_band(state)
    return row


def run(modes):
    with threadpool_limits(limits=1):
        model = SphericalValve(modes, initial_degree=40, budget_seconds=30.)
        model.started = time.perf_counter()
        solver = DOP853(model.rhs, 0., model.initial.ravel(), 2.,
                        max_step=.002, rtol=3e-11, atol=1e-13)
        checkpoints = list(np.arange(.1, 2.001, .1)); rows = [diagnostic(model,0.,model.initial)]
        status, message, last_t, last_y = 'running', None, 0., model.initial.copy()
        accepted = 0
        try:
            while solver.status == 'running':
                result = solver.step()
                if solver.status == 'failed':
                    status, message = 'integrator_failed', result
                    break
                accepted += 1
                last_t, last_y = float(solver.t), solver.y.reshape(8,modes).copy()
                if checkpoints and solver.t >= checkpoints[0]-1e-13:
                    dense = solver.dense_output()
                    while checkpoints and solver.t >= checkpoints[0]-1e-13:
                        target=checkpoints.pop(0)
                        row=diagnostic(model,target,dense(target).reshape(8,modes),band=abs(round(target)-target)<1e-8)
                        rows.append(row)
                        print(modes, target, row['neck_radius'], row['theta_plus_over_H'],
                              row['Hamiltonian_RMS_over_H_squared'], flush=True)
                if time.perf_counter()-model.started > 30:
                    status, message='time_budget', '30 second wall-clock cap'
                    break
            else:
                status = 'finished'
        except (ValueError, RuntimeError, TimeoutError, FloatingPointError) as exc:
            status, message='exception', str(exc)
        if not rows or abs(rows[-1]['proper_time']-last_t)>1e-9:
            try: rows.append(diagnostic(model,last_t,last_y,band=True))
            except Exception as exc: message=f'{message}; final diagnostic: {exc}'
        out=dict(schema='nsc-one-mode-neck-to2-prototype-v1', modes=modes,
                 preparation='same production SphericalValve epsilon=1 initial_degree=40',
                 proper_endpoint_requested=2.,accepted_end=last_t,status=status,message=message,
                 elapsed_seconds=time.perf_counter()-model.started,accepted_steps=accepted,
                 nfev=model.calls,rows=rows,
                 source_sha256=hashlib.sha256(Path('src/recursive_horizons/nsc_atlas_valve_spherical.py').read_bytes()).hexdigest(),
                 scope='Bounded numerical continuation, not certified continuum, event horizon or regular child; reflection tube j=0')
        Path(f'.build/valve-neck-to2-{modes}.json').write_text(json.dumps(out,indent=2)+'\n')
        print('END',modes,status,last_t,out['elapsed_seconds'],message,flush=True)


if __name__ == '__main__': run(int(sys.argv[1]))
