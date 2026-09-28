"""Conditional nonzero-velocity analytic crossing; no physical tuple selection."""
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import sympy as sp

THEOREM={'author':'Feng Rong','title':'The Briot-Bouquet systems and the center families for holomorphic dynamical systems',
    'identifier':'arXiv:1307.3823v1','locator':'Proposition 2.1, page 2',
    'url':'https://arxiv.org/pdf/1307.3823v1#page=2',
    'published_url':'https://www.sciencedirect.com/science/article/pii/S0001870813002302',
    'hypotheses':['holomorphic f(x,y) near (0,0)','f(0,0)=0','no positive-integer eigenvalue of f_y(0,0)'],
    'conclusion':'unique holomorphic germ y(x) with y(0)=0',
    'zero_eigenvalue_permitted':True,'source_checked':'primary author PDF and published Proposition 2.1 on 2026-09-18'}


@lru_cache(maxsize=1)
def crossing_identities():
    x,y,s,zc,wD=sp.symbols('x y s z_c w_D',real=True)
    t,t0=sp.symbols('t t0',real=True,nonzero=True)
    delta1=sp.Symbol('Delta1',positive=True);R=sp.Function('R')
    tp=sp.Symbol('t_x',real=True)
    # d_z=1/t*d_x, w_z=1/t and w_zz=-t_x/t^3.
    flux=delta1/t**2-delta1*x*tp/t**3-R(wD+x,zc+y)
    transformed=x*tp-t+R(wD+x,zc+y)*t**3/delta1
    transformation=sp.factor(flux+delta1*transformed/t**3)
    f=sp.Matrix([x*(t0+s),t0+s-R(wD+x,zc+y)*(t0+s)**3/delta1])
    origin={x:0,y:0,s:0};compatible={R(wD,zc):delta1/t0**2}
    constant=f.subs(origin,simultaneous=True).subs(compatible)
    J=f.jacobian((y,s)).subs(origin,simultaneous=True).subs(compatible)
    Rz=sp.Symbol('R_z_at_crossing',real=True)
    J=J.xreplace({atom:Rz for atom in J.atoms(sp.Subs) if atom.expr.has(sp.Derivative)})
    expected=sp.Matrix([[0,0],[-Rz*t0**3/delta1,-2]])
    lam=sp.Symbol('lambda');characteristic=sp.factor(J.charpoly(lam).as_expr())
    # Original operators after solving the momentum primitive for U.
    A,B,C,d,e,v,acc,Q,Dval=sp.symbols('A B C d e v acceleration Q D_value',real=True)
    Delta=A*e-B*d
    dEN=A*(Q-e*acc)+B*d*acc+d*C*v*v+d*Dval
    fluxres=Delta*acc-d*C*v*v-(A*Q+d*Dval)
    lapse_reconstruction=sp.factor(dEN+fluxres)
    Sb,F,jerk,Qz=sp.symbols('S_beta F jerk Q_z',real=True)
    momentum=sp.factor((Qz-e*jerk)+e*jerk+F*v+Sb).subs(Qz,-Sb-F*v)
    Rw=sp.Symbol('R_w_at_crossing',real=True);vc=sp.Symbol('v_c',real=True,nonzero=True)
    curvature=(Rw+Rz/vc)/(3*delta1)
    differentiated_flux=sp.factor(3*delta1*vc*curvature-Rw*vc-Rz)
    t1=-(Rw*t0**3+Rz*t0**4)/(3*delta1)
    inverse_first_jet=sp.factor(-t1/t0**3-curvature.subs(vc,1/t0))
    residuals={'inverse_coordinate_flux':str(transformation),'centered_constant_first':str(constant[0]),
        'centered_constant_second':str(sp.factor(constant[1])),
        'centered_Jacobian':str(sum(sp.factor(q)**2 for q in J-expected)),
        'characteristic_polynomial':str(sp.factor(characteristic-lam*(lam+2))),
        'lapse_reconstruction':str(lapse_reconstruction),'momentum_reconstruction':str(momentum),
        'first_finite_curvature':str(differentiated_flux),'inverse_first_jet':str(inverse_first_jet)}
    if set(residuals.values())!={'0'}:raise ArithmeticError('crossing transformation identity failed: '+str(residuals))
    return {'necessary_compatibility':'R(w_D,z_c)=Delta1*v_c²>0; t0=1/v_c!=0',
        'centered_right_side':[str(q) for q in f],
        'Jacobian':[[str(J[i,j]) for j in range(2)] for i in range(2)],
        'eigenvalues':[0,-2],'positive_integer_eigenvalues':[],
        'first_curvature':str(curvature),'residuals':residuals,
        'root_symbol':'w_D is the exact affine principal root; no numerical midpoint substitution',
        'analyticity':'R is the committed real polynomial A(w)*(Pi0-S_beta*z-G(w))+d*D(w)',
        'reality':'conjugation invariance and uniqueness of the holomorphic germ',
        'invertibility':'Z_prime(0)=t0!=0; analytic inverse-function theorem',
        'uniqueness_scope':'analytic germs at fixed compatible real symbolic data',
        'finite_smooth_existence':'real analytic implies finite derivatives of every order locally'}


SOURCES=('src/recursive_horizons/nsc_incoming_surface_crossing.py','tests/test_nsc_incoming_surface_crossing.py','docs/nsc-incoming-surface-crossing.md')
INPUTS=('results/development/nsc-incoming-surface-coefficients.json','results/development/nsc-incoming-surface-regular-branch.json')
OUTPUT='results/development/nsc-incoming-surface-crossing.json'


def make_record(root):
    import mpmath as mp
    from .nsc_incoming_lapse_coefficient import _interval
    from .nsc_incoming_vacuum_tail_bound import _precision
    root=Path(root);digest=lambda p:sha256((root/p).read_bytes()).hexdigest()
    records=[json.loads((root/p).read_text()) for p in INPUTS]
    for record in records:
        for section in ('source_hashes','input_hashes'):
            for path,h in record[section].items():
                if digest(path)!=h:raise ValueError('crossing input changed: '+path)
        if 'payload' in record and digest(record['payload']['path'])!=record['payload']['sha256']:
            raise ValueError('crossing input payload changed')
    coefficient,branch=records;C=coefficient['coefficients']['C_interval'];d=branch['branch']['imported_unchanged']['d']
    if coefficient['exact_polynomials']['identities']['Delta_w_plus_dC']!='0':raise ValueError('exact affine slope identity required')
    if d[1]>=0 or C[0]<=0:raise ValueError('certified d<0 and C>0 required')
    with _precision(50):slope=_interval(-mp.iv.mpf(d)*mp.iv.mpf(C))
    if slope[0]<=0:raise ArithmeticError('positive affine principal slope not certified')
    proof=crossing_identities()
    return {'schema':'NSC-INCOMING-SURFACE-CROSSING-v1','accountable_author':'Douglas Ek',
        'status':'PASS: conditional real-analytic nonzero-velocity crossing germ',
        'theorem':THEOREM,'proof':proof,'imported_intervals':{'d':d,'C':C,'Delta1':slope,
            'w_D':branch['branch']['intervals']['w_principal']},
        'conditions':['all coefficients and z_c,Pi0,S_N,S_beta are finite real',
            'R(w_D,z_c)=Delta1*v_c² with v_c real and nonzero',
            'same included compatible source/operator and frozen intrinsic slice'],
        'conclusion':['a unique real-analytic crossing germ exists for each such fixed compatible tuple and slope sign',
            'w(z_c)=w_D and w_prime(z_c)=v_c; U=(Pi0-S_beta*z-G(w)-e*w_second)/d is analytic',
            'both included lapse and shift constraints hold at and near the crossing'],
        'residuals':{k:0. for k in proof['residuals']},'verification_tolerances':{'exact_symbolic_identities':0.},
        'scope':{'physical_tuple_selected':False,'source_constants_numerically_substituted':False,
            'particular_trajectory_reaches_root_proved':False,'arbitrary_smooth_uniqueness_claimed':False,
            'zero_velocity_or_zero_R_case_analyzed':False,'global_existence_claimed':False,
            'old_coefficient_or_source_generators_rerun':False,'metric_evolution':False,'pressure_equations_or_Gamma_rest_closed':False},
        'source_hashes':{p:digest(p) for p in SOURCES},'input_hashes':{p:digest(p) for p in INPUTS},
        'reproducer':'PYTHONPATH=src python3 -m recursive_horizons.nsc_incoming_surface_crossing --check'}


def main():
    import argparse
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true');args=parser.parse_args()
    root=Path(__file__).resolve().parents[2];record=make_record(root);path=root/OUTPUT
    if args.write:path.write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+'\n')
    elif json.loads(path.read_text())!=record:raise ValueError('crossing record differs from authenticated replay')
    print(json.dumps({'status':record['status'],'eigenvalues':record['proof']['eigenvalues'],
        'Delta1':record['imported_intervals']['Delta1'],'residuals':record['residuals']},indent=2))

if __name__=='__main__':main()
