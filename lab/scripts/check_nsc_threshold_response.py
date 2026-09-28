#!/usr/bin/env python3
"""Positive zero-energy DtN jump does not imply a mass gap: explicit throat."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import sympy as sp
from check_nsc_scale_closure import compare

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'results/nsc-3-threshold-response.json'


def calculate():
    t=sp.symbols('t',positive=True)
    ip=sp.exp(3*t)/6+sp.exp(t)/2-sp.Rational(2,3)
    im=sp.Rational(2,3)-sp.exp(-t)/2-sp.exp(-3*t)/6
    assert sp.simplify((sp.diff(ip,t)-sp.exp(2*t)*sp.cosh(t)).rewrite(sp.exp))==0
    assert sp.simplify((sp.diff(im,t)-sp.exp(-2*t)*sp.cosh(t)).rewrite(sp.exp))==0
    limit=sp.limit(1/ip+1/im,t,sp.oo)
    assert limit==sp.Rational(3,2)
    rows=[]
    for radius in (4.,8.,16.,32.):
        depth=np.arcsinh(radius)
        plus=float(ip.subs(t,depth)); minus=float(im.subs(t,depth))
        analytic_parent=1+1/plus; analytic_child=-1+1/minus
        def equation(rho,state):
            # Upper component of D²: A†A=-d²+w²-w'.
            potential=1/(1+rho*rho)+rho/(1+rho*rho)**1.5
            return [state[1],potential*state[0]]
        numerical=[]
        for endpoint,normal in ((radius,-1),(-radius,1)):
            result=solve_ivp(equation,(endpoint,0),[0.,1.],method='DOP853',
                             rtol=2e-11,atol=2e-12,max_step=.1)
            if not result.success: raise RuntimeError(result.message)
            u,du=result.y[:,-1]
            numerical.append(float(normal*du/u))
        error=max(abs(numerical[0]-analytic_parent),abs(numerical[1]-analytic_child))
        assert error<1e-9
        rows.append({'R':radius,'parent_rho_positive_outward_minus':analytic_parent,
                     'child_rho_negative_outward_plus':analytic_child,
                     'jump':analytic_parent+analytic_child,'independent_ODE_maps':numerical,
                     'maximum_residual':error})
    return {'schema':'nsc-threshold-response-v1','artifact_id':'NSC-3-THRESHOLD-RESPONSE',
            'classification':'positive_zero_energy_throat_DtN_jump_coexists_with_gapless_bulk_spectrum',
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'operator':'upper D²=A†A=-d²+w²-w_prime, w=1/sqrt(1+rho²)',
            'boundary':'u(0)=1, u(plus_or_minus R)=0; outward normals at throat are minus on parent, plus on child',
            'exact':{'t':'asinh(R)','I_parent':str(ip),'I_child':str(im),
                     'N_parent':'1+1/I_parent','N_child':'-1+1/I_child',
                     'jump':'1/I_parent+1/I_child','infinite_box_threshold_jump':str(limit),
                     'boundary_unit':'1/L_star with L_star=1'},
            'finite_boxes':rows,
            'interpretation':'the threshold boundary stiffness is nonzero even though the separately proved essential spectrum of D is R; DtN(0) is not a constant Dirac mass Phi',
            'zero_energy_caution':'m(E)=v/u may have a pole at E=0; evaluate E*m(E) by the D² boundary problem or its limit, not substitution of m(0)',
            'nonclaims':{'physical_Phi_equals_3_over_2':False,'global_spectrum_gapped':False,
                         'full_warped_or_recursive_threshold_derived':False,'new_to_world_physics_established':False},
            'terminal':True}


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true');group.add_argument('--output',type=Path)
    args=parser.parse_args();record=calculate()
    if args.check:
        compare(json.loads(OUTPUT.read_text()),record)
        print('Threshold DtN counterexample reproduced; all fields checked.')
    elif args.output:
        with args.output.open('x') as stream:json.dump(record,stream,sort_keys=True,indent=2,allow_nan=False);stream.write('\n')
        print(args.output)
    else:print(json.dumps(record,indent=2))


if __name__=='__main__':main()
