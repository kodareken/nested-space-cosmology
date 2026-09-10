#!/usr/bin/env python3
"""Derive the spherical source/constraint map and bind its new ADM vertices."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_adm_source import (
    adm_hamiltonian, covariant_history_operator, canonical_force_gradients,
)
from recursive_horizons.nsc_covariant_operator import smooth_metric, euclidean_operator, SIGMA3
from recursive_horizons.nsc_influence import ground_covariance
from recursive_horizons.nsc_regulated import RegulatedOperator, OperatorConventions
from check_nsc_compact_boundary_action import authenticated_record
from check_nsc_compact_casimir import compare, native
from check_nsc_vacuum_charge_matching import hashes

OUTPUT = ROOT/'results/development/adm-source-constraints.json'
INPUTS = ('results/development/vacuum-matched-ctp.json',
          'results/nsc-4-dirac-tetrad.json', 'results/nsc-9-covariant-source.json')
SOURCES = ('scripts/derive_nsc_adm_source.py', 'src/recursive_horizons/nsc_adm_source.py',
           'src/recursive_horizons/nsc_covariant_operator.py',
           'src/recursive_horizons/nsc_influence.py',
           'src/recursive_horizons/nsc_regulated.py',
           'docs/nsc-adm-source-constraints.md',
           'scripts/check_nsc_compact_boundary_action.py',
           'scripts/check_nsc_compact_casimir.py',
           'scripts/check_nsc_vacuum_charge_matching.py')


def exact_identities():
    t, x = sp.symbols('t x', real=True)
    n, q, r, b = (sp.Function(s)(t, x) for s in ('N', 'q', 'r', 'beta'))
    zt, zx = (sp.Function(s)(t, x) for s in ('xi_t', 'xi_x'))
    fields = (n, q, r, b)
    advect = lambda a: zt*sp.diff(a, t)+zx*sp.diff(a, x)
    variations = (
        advect(n)+n*(sp.diff(zt, t)-b*sp.diff(zt, x)),
        advect(q)+q*(sp.diff(zx, x)+b*sp.diff(zt, x)),
        advect(r),
        advect(b)+sp.diff(zx, t)-b*sp.diff(zx, x)+b*sp.diff(zt, t)
                  -(b*b+n*n/(q*q))*sp.diff(zt, x),
    )
    exact = {}
    def zero(name, expression):
        assert sp.simplify(expression) == 0, name
        exact[name] = '0'
    g = sp.Matrix([[n*n-q*q*b*b, -q*q*b], [-q*q*b, -q*q]])
    coords, vector = (t, x), (zt, zx)
    for i in range(2):
        for j in range(2):
            induced = sum(sp.diff(g[i,j], a)*da for a, da in zip(fields, variations))
            lie = advect(g[i,j])+sum(g[k,j]*sp.diff(vector[k], coords[i])
                                    +g[i,k]*sp.diff(vector[k], coords[j]) for k in range(2))
            zero(f'ADM_Lie_metric_{i}{j}', induced-lie)
    zero('ADM_Lie_sphere_radius', 2*r*variations[2]-advect(r*r))

    fn, fq, fr, fb = (sp.Function(s)(t, x) for s in ('F_N', 'F_q', 'F_r', 'F_beta'))
    force = (fn, fq, fr, fb)
    variation = sum(f*da for f, da in zip(force, variations))
    euler = lambda z: (sp.diff(variation, z)-sp.diff(sp.diff(variation, sp.diff(z,t)),t)
                      -sp.diff(sp.diff(variation, sp.diff(z,x)),x))
    energy = n*fn+b*fb
    flux = -n*b*fn-(b*b+n*n/(q*q))*fb+q*b*fq
    rt = sum(f*sp.diff(a,t) for f,a in zip(force, fields))-sp.diff(energy,t)-sp.diff(flux,x)
    rx = sum(f*sp.diff(a,x) for f,a in zip(force, fields))-sp.diff(fb,t)-sp.diff(q*fq-b*fb,x)
    zero('time_Ward_from_Lie_variation', euler(zt)-rt)
    zero('radial_Ward_from_Lie_variation', euler(zx)-rx)
    transport_b = b*sp.diff(fb,x)+sp.diff(n,x)*fn+2*sp.diff(b,x)*fb
    transport_n = (b*sp.diff(fn,x)+sp.diff(b,x)*fn+n/q**2*sp.diff(fb,x)
                   +2*(sp.diff(n,x)/q**2-n*sp.diff(q,x)/q**3)*fb)
    for name, expr in [('energy',rt),('momentum',rx)]:
        reduced = expr.subs({fq:0,fr:0}).doit()
        zero(f'homogeneous_constraint_{name}',
             reduced.subs({sp.diff(fb,t):transport_b, sp.diff(fn,t):transport_n}))

    rho, j, pr, pt = sp.symbols('rho j p_radial p_sphere')
    coframe = sp.Matrix([[n,0],[q*b,q]])
    stress = coframe.T*sp.Matrix([[rho,j],[j,pr]])*coframe
    inverse = g.inv()
    volume = 4*sp.pi*n*q*r*r
    expected = {n:4*sp.pi*q*r*r*rho, b:4*sp.pi*q*q*r*r*j,
                q:-4*sp.pi*n*r*r*pr}
    for label, a in [('N',n),('beta',b),('q',q)]:
        value = -volume/2*sum(stress[i,k]*sp.diff(inverse[i,k],a)
                             for i in range(2) for k in range(2))
        zero(f'stress_frame_{label}',value-expected[a])
    zero('stress_frame_r', -volume/2*(4*pt/r)+8*sp.pi*n*q*r*pt)
    return exact


def calculate():
    for p in INPUTS:
        authenticated_record(p)
    exact = exact_identities()
    base = smooth_metric(12, general=True)
    nt, period = 4, 3.1
    tau = np.arange(nt)*period/nt
    modes = np.arange(nt)-nt//2+.5
    fourier = np.exp(2j*np.pi*tau[:,None]*modes[None,:]/period)/np.sqrt(nt)
    momentum = (fourier*(2*np.pi*modes/period))@fourier.conj().T
    angle = 2*np.pi*base.x/base.length
    metrics, shifts = [], []
    all_directions, all_shifts = [], []
    for j in range(nt):
        phase = 2*np.pi*j/nt
        metrics.append(replace(base, lapse=base.lapse*np.exp(.02*np.cos(phase)*np.cos(angle)),
                               radial_scale=base.radial_scale*np.exp(.015*np.sin(phase+angle)),
                               sphere_radius=base.sphere_radius*np.exp(.01*np.cos(phase-angle))))
        shifts.append(.08*np.cos(angle)+.01*np.sin(phase)*np.sin(angle))
        all_directions.append([.06*(1+.2*np.cos(angle))*np.cos(phase),
                               .04*np.sin(angle+phase), .05*np.cos(2*angle-phase)])
        all_shifts.append(.03*np.sin(angle+.7*phase))
    directions, shifts_d = np.array(all_directions), np.array(all_shifts)
    baseline = covariant_history_operator(metrics, shifts, momentum)['operator']
    assert np.max(abs(baseline-baseline.conj().T)) < 1e-12
    conventions = OperatorConventions(cutoff=3.)
    rows = []
    for label in ('N','q','r','beta','mixed'):
        ds = np.zeros_like(directions); db = np.zeros_like(shifts_d)
        if label in ('N','q','r'):
            k = ('N','q','r').index(label); ds[:,k]=directions[:,k]
        elif label == 'beta': db=shifts_d.copy()
        else: ds,db=directions,shifts_d
        owner = covariant_history_operator(metrics, shifts, momentum, ds, db)
        def build(s):
            changed = [replace(m, lapse=m.lapse*np.exp(s*d[0]),
                               radial_scale=m.radial_scale*np.exp(s*d[1]),
                               sphere_radius=m.sphere_radius*np.exp(s*d[2])) for m,d in zip(metrics,ds)]
            return covariant_history_operator(changed, np.array(shifts)+s*db, momentum)['operator']
        step=.02
        dm2,dm1,dp1,dp2 = (build(s*step) for s in (-2,-1,1,2))
        fd1=(dm2-8*dm1+8*dp1-dp2)/(12*step)
        fd2=(-dp2+16*dp1-30*baseline+16*dm1-dm2)/(12*step**2)
        first,second=RegulatedOperator(baseline,conventions).variation(owner['first'],owner['second'])
        actions=[RegulatedOperator(d,conventions).action() for d in (dm2,dm1,baseline,dp1,dp2)]
        action1=(actions[0]-8*actions[1]+8*actions[3]-actions[4])/(12*step)
        action2=(-actions[4]+16*actions[3]-30*actions[2]+16*actions[1]-actions[0])/(12*step**2)
        row={'direction':label,'first_operator_error':float(np.max(abs(fd1-owner['first']))),
             'second_operator_error':float(np.max(abs(fd2-owner['second']))),
             'first_action_variation':first,'second_action_variation':second,
             'first_action_difference':action1-first,'second_action_difference':action2-second}
        assert row['first_operator_error']<2e-9 and row['second_operator_error']<2e-7, row
        assert abs(row['first_action_difference'])<1e-8 and abs(row['second_action_difference'])<5e-7, row
        rows.append(row)
    beta=np.kron(SIGMA3,np.eye(base.points))
    fiber=covariant_history_operator([base],[np.zeros(base.points)],np.array([[.8]]))['operator']
    fiber_error=float(np.max(abs(fiber-beta@euclidean_operator(base,.8)@beta)))
    h=adm_hamiltonian(base,shifts[0]);_,_,c=ground_covariance(h)
    forces=canonical_force_gradients(base,shifts[0],c)
    energy_error=float(np.dot(base.lapse,forces['N'])+np.dot(shifts[0],forces['beta'])-forces['energy'])
    assert fiber_error<1e-12 and abs(energy_error)<1e-11
    return native({'schema':'NSC-ADM-SOURCE-CONSTRAINTS-v1',
        'status':'coupled spherical operator vertices, stress projection and continuum Ward/constraint equations derived',
        'source_hashes':hashes(SOURCES),'input_hashes':hashes(INPUTS),
        'exact_identities':exact,
        'domain':{'signature':'+---','metric':'N²dt²-q²(dx+beta dt)²-r²dOmega²',
                  'canonical_spinor':'u=r sqrt(q) psi','Euclidean_spinor':'chi=sqrt(N) u',
                  'spatial_scheme':'existing factorized Fourier scheme; no exact finite Leibniz rule claimed',
                  'numerical_control':'4 Euclidean time nodes by 12 spatial nodes, one radial block; no vacuum/thermal source interpretation',
                  'Euclidean_time_period':period,'control_cutoff':3.,'finite_difference_step':.02,
                  'physical_shift_continuation':'b_E=-i beta_L for t=-i tau',
                  'source_frame':'spherical angular average; full angular multiplicities must be restored'},
        'new_vertex_binding':rows,'old_static_fiber_unitary_error':fiber_error,
        'canonical_force_gradients':forces,'energy_constraint_projection_error':energy_error,
        'ward_equations':{
            'energy_density':'E=N F_N+beta F_beta',
            'energy_flux':'J=-N beta F_N-(beta²+N²/q²) F_beta+q beta F_q',
            'energy':'partial_t E+partial_x J=sum_A F_A partial_t A',
            'momentum':'partial_t F_beta+partial_x(q F_q-beta F_beta)=sum_A F_A partial_x A',
            'constraint_beta':'(partial_t-beta partial_x)F_beta=N_x F_N+2 beta_x F_beta when F_q=F_r=0',
            'constraint_N':'(partial_t-beta partial_x)F_N=beta_x F_N+(N/q²)partial_x F_beta+2(N_x/q²-N q_x/q³)F_beta when F_q=F_r=0'},
        'scope':{'full_absolute_state_source_evaluated':False,'joint_geometry_solution_found':False,
                 'continuum_Cauchy_wellposedness_proved':False,'quantum_constraint_algebra_proved':False,
                 'old_physics_generators_rerun':False},
        'comparison':{'fields':'all','float_atol':3e-9,'float_rtol':3e-8,
                      'exact':'structure, strings, non-float values and source/input hashes','exceptions':[]}})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true');parser.add_argument('--output',type=Path,default=OUTPUT)
    args=parser.parse_args();result=calculate()
    if args.check:
        compare(json.loads(args.output.read_text()),result)
        print('ADM source assembly reproduced: exact Ward identities and all record fields agree.')
    else:
        with args.output.open('x') as f:
            json.dump(result,f,sort_keys=True,indent=2,allow_nan=False);f.write('\n')
        print(f'Wrote {args.output}')


if __name__=='__main__': main()
