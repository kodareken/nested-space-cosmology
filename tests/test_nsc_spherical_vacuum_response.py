"""Same-action constrained vacuum derivation; no data archive or trajectory."""
import time
import hashlib
import json

import pytest
import sympy as sp

from recursive_horizons.nsc_spherical_feedback_action import alpha_of, first_order_density
from recursive_horizons.nsc_spherical_timelike_boundary import timelike_gravitational_response
from recursive_horizons.nsc_spherical_vacuum_response import (
    vacuum_coefficients, vacuum_curvature_profile, vacuum_exterior_response,
)


ACTIVE_CPU=[0.0]


@pytest.fixture(autouse=True)
def bounded():
    started=time.process_time()
    yield
    ACTIVE_CPU[0]+=time.process_time()-started
    assert ACTIVE_CPU[0]<60,"vacuum response controls exhausted their symbolic allocation"


def zero(expr):
    assert sp.cancel(sp.expand(expr))==0


def el(density,field,t,R):
    return (sp.diff(density,field)-sp.diff(sp.diff(density,sp.diff(field,t)),t)
            -sp.diff(sp.diff(density,sp.diff(field,R)),R))


@pytest.fixture(scope='module')
def linear():
    t,R=sp.symbols('t R',real=True,positive=True)
    eps=sp.symbols('epsilon',real=True)
    A=sp.symbols('A',positive=True);CW=sp.symbols('C_W',nonzero=True,real=True)
    L,Q,r,chi,beta=[sp.Function(name)(t,R) for name in ('L','Q','r','chi','beta')]
    n,b,c=[sp.Function(name)(t,R) for name in ('n','b','c')]
    density=first_order_density(L,sp.diff(L,R),Q,sp.diff(Q,t),sp.diff(Q,R),
                               r,sp.diff(r,t),sp.diff(r,R),chi,sp.diff(chi,t),
                               sp.diff(chi,R),beta,sp.diff(beta,R),A,CW,0,0)
    sub={L:(1+eps*n)/R,Q:(1+eps*b)/R,r:R,chi:eps*c,beta:0}
    restricted=density.subs(sub,simultaneous=True).doit()
    quadratic=sp.expand(sp.diff(restricted,eps,2).subs(eps,0)/2)
    full={name:sp.expand(sp.diff(el(density,field,t,R).subs(sub,simultaneous=True).doit(),eps).subs(eps,0))
          for name,field in (('lapse',L),('radial',Q),('radius',r),('auxiliary',chi),('shift',beta))}
    return t,R,eps,A,CW,n,b,c,quadratic,full


def test_quadratic_action_and_all_pre_gauge_equations(linear):
    t,R,eps,A,CW,n,b,c,quadratic,full=linear
    alpha=alpha_of(CW);k=4*sp.pi*A/alpha
    expected=(4*alpha*sp.diff(b,t)*sp.diff(c,t)-16*sp.pi*A*R*b*sp.diff(n,R)+8*sp.pi*A*b*b
              -4*alpha*sp.diff(c,R)*sp.diff(n,R)+4*alpha*(n-b)*sp.diff(c,R)/R
              -alpha*c*c/R**2-4*alpha*c*(n+b)/R**2)
    zero(quadratic-expected)
    for name,field in (('lapse',n),('radial',b),('auxiliary',c)):
        factor=R if name in ('lapse','radial') else 1
        zero(full[name]/factor-el(quadratic,field,t,R))
    constraint=k*R*sp.diff(b,t)+sp.diff(c,t,R)+sp.diff(c,t)/R
    zero(full['shift']-4*alpha*constraint)
    # The full radius equation is retained, not dropped by areal gauge fixing.
    trace=(sp.diff(n,R,2)-sp.diff(b,t,2)+2*sp.diff(n,R)/R
           -2*sp.diff(b,R)/R-2*b/R**2)
    zero(full['radius']+16*sp.pi*A*R*trace)


def test_constraints_master_equation_and_schwarzschild_mass(linear):
    t,R,eps,A,CW,n,b,c,quadratic,full=linear
    alpha=alpha_of(CW);k=4*sp.pi*A/alpha;M=sp.symbols('M',real=True)
    bv=M/R-(sp.diff(c,R)+c/R)/(k*R)
    nr=M/R**2-(sp.diff(c,t,2)/R+2*sp.diff(c,R)/R**2+2*c/R**3)/k
    substitutions={sp.diff(n,R,2):sp.diff(nr,R),sp.diff(n,R):nr,
                   sp.diff(b,R):sp.diff(bv,R),sp.diff(b,t):sp.diff(bv,t),
                   sp.diff(b,t,2):sp.diff(bv,t,2),b:bv}
    for name in ('lapse','radial','shift','radius'):
        zero(full[name].subs(substitutions,simultaneous=True).doit())
    master=(sp.diff(c,t,2)-sp.diff(c,R,2)+2*sp.diff(c,R)/R+4*c/R**2-k*c/6-2*k*M/R)
    zero(full['auxiliary'].subs(substitutions,simultaneous=True).doit()-12*alpha*master/(k*R**2))
    static={n:-M/R,b:M/R,c:-12*M/R}
    for value in full.values():zero(value.subs(static,simultaneous=True).doit())


def test_reduced_master_action_including_boundary_flux(linear):
    t,R,eps,A,CW,n,b,c,quadratic,full=linear
    psi=sp.Function('psi')(t,R);alpha=alpha_of(CW);k=4*sp.pi*A/alpha
    cv=R*psi;bv=-(sp.diff(psi,R)+2*psi/R)/k
    reduced=quadratic.subs({c:cv,b:bv},simultaneous=True).doit()
    coefficients=vacuum_coefficients(A,CW)
    K=coefficients['quadratic_kinetic_coefficient'];mass2=coefficients['mass_squared']
    target=K*(sp.diff(psi,t)**2-sp.diff(psi,R)**2-(6/R**2+mass2)*psi**2)
    flux=4*alpha*n*psi-2*alpha*R*sp.diff(psi,t)**2/k+12*alpha*psi**2/(k*R)
    zero(reduced-target-sp.diff(flux,R))
    zero(mass2-A/(2*CW))
    zero(K+8*sp.pi*CW**2/(3*A))


def test_exact_outgoing_radial_solution_and_metric_reconstruction():
    R,q,k=sp.symbols('R q k',positive=True)
    wave=sp.exp(-q*R)*(1+3/(q*R)+3/(q*R)**2)
    zero(sp.diff(wave,R,2)-(q*q+6/R**2)*wave)
    ratio=-q-3*(q*R+2)/(R*((q*R)**2+3*q*R+3))
    zero(sp.diff(wave,R)-ratio*wave)
    integral=sp.exp(-q*R)*(1/q+3/(q*q*R))
    zero(sp.diff(integral,R)+wave)
    bv=-(sp.diff(wave,R)+2*wave/R)/k
    nv=bv+integral/6;squared_s=q*q+k/6
    zero(sp.diff(nv,R)+(squared_s*wave+2*sp.diff(wave,R)/R+4*wave/R**2)/k)


def test_face_response_linearization_respects_constraints():
    R,L,Q,A=sp.symbols('R L Q A',positive=True)
    CW=sp.symbols('C_W',nonzero=True,real=True)
    chi,lx,rx,cx,eps,n,b,nr,c,cr=sp.symbols('chi lx rx cx eps n b nr c cr',real=True)
    alpha=alpha_of(CW);k=4*sp.pi*A/alpha
    response=timelike_gravitational_response(L,Q,R,chi,lx,rx,cx,A,CW,orientation=-1)
    sub={L:(1+eps*n)/R,Q:(1+eps*b)/R,chi:eps*c,lx:-1/R**2+eps*(nr/R-n/R**2),rx:1,cx:eps*cr}
    linear={key:sp.diff(response[key].subs(sub,simultaneous=True),eps).subs(eps,0)
            for key in ('unreferenced_energy','unreferenced_area_pressure','auxiliary_response_per_proper_time')}
    zero(linear['unreferenced_energy']+16*sp.pi*A*R*b+4*alpha*cr)
    zero(linear['unreferenced_area_pressure']+2*A*(nr-b/R)+alpha*cr/(2*sp.pi*R**2))
    zero(linear['auxiliary_response_per_proper_time']-4*alpha*(nr+b/R))


@pytest.mark.parametrize('cw',[-0.1,0.1])
def test_laplace_branch_profile_and_response_coefficients(cw):
    R0=1.3;A=1.0;s=4+1j
    row=vacuum_exterior_response(s,R0,A,cw)
    opposite=vacuum_exterior_response(s,R0,A,cw,orientation=1)
    assert row['decay_q'].real>0
    assert vacuum_curvature_profile(s,R0,R0,A,cw)==pytest.approx(1)
    assert abs(vacuum_curvature_profile(s,10*R0,R0,A,cw))<1e-10
    conjugate=vacuum_exterior_response(s.conjugate(),R0,A,cw)
    assert conjugate['chi_R_per_chi']==pytest.approx(row['chi_R_per_chi'].conjugate())
    for key in ('unreferenced_energy_per_chi','unreferenced_area_pressure_per_chi','auxiliary_proper_response_per_chi'):
        assert row[key]==pytest.approx(-opposite[key])
    alpha=row['alpha']
    bv=row['radial_metric_b_per_chi'];nr=row['lapse_n_R_per_chi'];cr=row['chi_R_per_chi']
    assert row['unreferenced_energy_per_chi']==pytest.approx(-16*float(sp.pi)*A*R0*bv-4*alpha*cr)
    assert row['unreferenced_area_pressure_per_chi']==pytest.approx(-2*A*(nr-bv/R0)-alpha*cr/(2*float(sp.pi)*R0**2))
    assert row['auxiliary_proper_response_per_chi']==pytest.approx(4*alpha*(nr+bv/R0))


def test_negative_far_exterior_form_and_laplace_domain():
    gamma=sp.symbols('gamma',positive=True)
    width=4*sp.pi/gamma;far=4/gamma
    # A sine supported beyond max(R0,far) is in the Dirichlet form domain.
    # Its kinetic quotient is pi²/width² and 6/R²<=6/far² on the support.
    zero(sp.pi**2/width**2+6/far**2-gamma**2+9*gamma**2/16)
    row=vacuum_coefficients(1,-.1)
    assert row['mass_squared']<0 and row['quadratic_kinetic_coefficient']<0
    with pytest.raises(ValueError,match='growth bound'):
        vacuum_exterior_response(2,1,1,-.1)
    with pytest.raises(ValueError,match='C_W=0'):
        vacuum_exterior_response(4,1,1,0)
    with pytest.raises(ValueError,match='positive radius'):
        vacuum_exterior_response(4,0,1,-.1)


def test_scalar_producer_rejects_changed_couplings_or_target_selected_ruler(tmp_path):
    from measure_nsc_vacuum_gravitational_response import ACTION, ROOT, load_parameters
    action_hash=hashlib.sha256((ROOT/ACTION).read_bytes()).hexdigest()
    causal={'status':'completed_causal_photon_response','couplings':{'A':1.,'C_W':-.1},
            'input':{'payload_sha256':'fixture'},'source_sha256':{ACTION:action_hash}}
    resolution={'status':'measured_instantaneous_correlation','parameters':{'A':1.,'C_W':-.1},
                'domain':{'reference_length':1.3,'target_observation_used':False},
                'input':{'payload_sha256':'fixture'},'source_sha256':{ACTION:action_hash}}
    first=tmp_path/'first.json';second=tmp_path/'second.json'
    first.write_text(json.dumps(causal));second.write_text(json.dumps(resolution))
    assert load_parameters(first,second)==(1.,-.1,1.3)
    resolution['parameters']['C_W']=.1;second.write_text(json.dumps(resolution))
    with pytest.raises(ValueError,match='couplings differ'):load_parameters(first,second)
    resolution['parameters']['C_W']=-.1;resolution['domain']['target_observation_used']=True
    second.write_text(json.dumps(resolution))
    with pytest.raises(ValueError,match='ruler'):load_parameters(first,second)
