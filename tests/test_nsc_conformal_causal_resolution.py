"""Small manufactured checks of clocked canonical causal elimination."""
import time

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.linalg import expm
import sympy as sp
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_conformal_causal_resolution import (
    CANONICAL_FIELDS,CanonicalVariation,canonical_jvp,causal_split_rates,
    clocked_canonical_jvp,combine_variations,environment_projection,
    fixed_receiver_clock_redshift,initial_packet,observer_values,
    proper_time_reference,reconstructed_tangent,resolution_rk4_step,
    retained_projection,to_canonical_tangent,to_physical_tangent,
)
from recursive_horizons.nsc_conformal_resolution import from_log_state,log_state
from recursive_horizons.nsc_conformal_tangent_evolution import FIELDS,variational_rk4_step
from recursive_horizons.nsc_spherical_coupling import CauchyState
from recursive_horizons.nsc_spherical_galerkin_coupling import (
    blank_state,build_grid,compose_fine_hamiltonian,manufactured_columns,
)


CPU=[0.0]


@pytest.fixture(autouse=True)
def allocation():
    started=time.process_time()
    with threadpool_limits(limits=1):yield
    CPU[0]+=time.process_time()-started
    assert CPU[0]<60,'clocked resolution controls exhausted their allocation'


@pytest.fixture(scope='module')
def example():
    grid=build_grid(10,quadrature=40,gauge='conformal')
    state=blank_state(grid,*manufactured_columns(grid.nf))
    x=2*np.pi*grid.xi_g/grid.length
    state.Q=.8+.09*np.sin(x)+.01*np.cos(3*x)
    state.r=2.1+.13*np.cos(x)+.01*np.sin(4*x)
    state.chi=.04*np.sin(2*x)
    state.p_Q=.02*np.cos(2*x);state.p_r=.03*np.sin(x);state.p_chi=.01*np.cos(3*x)
    rng=np.random.default_rng(43)
    values={n:.03*rng.normal(size=grid.ng) for n in FIELDS[:6]}
    values.update({n:.03*(rng.normal(size=(grid.nf,6))+1j*rng.normal(size=(grid.nf,6)))
                   for n in FIELDS[6:]})
    return grid,state,CauchyState(**values),grid.A_g[5]


def close(first,second,atol=1e-11,rtol=1e-10):
    for n in CANONICAL_FIELDS:assert_allclose(getattr(first,n),getattr(second,n),atol=atol,rtol=rtol)


def test_log_tangent_roundtrip_and_canonical_pairing(example):
    _grid,state,u,_row=example
    mapped=to_canonical_tangent(state,u);back=to_physical_tangent(state,mapped)
    for n in FIELDS:assert_allclose(getattr(back,n),getattr(u,n),atol=1e-15,rtol=1e-14)
    other=CauchyState(**{n:np.roll(getattr(u,n),1,axis=0) for n in FIELDS})
    v=to_canonical_tangent(state,other)
    old=sum(np.sum(getattr(u,q)*getattr(other,p)-getattr(u,p)*getattr(other,q))
            for q,p in [('r','p_r'),('Q','p_Q'),('chi','p_chi')])
    new=sum(np.sum(getattr(mapped,q)*getattr(v,p)-getattr(mapped,p)*getattr(v,q))
            for q,p in [('a','pi_a'),('b','pi_b'),('chi','p_chi')])
    assert_allclose(old,new,atol=1e-15)


def test_log_and_clock_JVP_against_independent_coordinate_difference(example):
    grid,state,u,row=example;v=to_canonical_tangent(state,u);h=3e-6
    rate,bundle=compose_fine_hamiltonian(grid,state)
    exact=canonical_jvp(grid,state,v,rate=rate,bundle=bundle)
    clocked=clocked_canonical_jvp(grid,state,v,row,rate=rate,bundle=bundle)
    answers=[];proper=[]
    for sign in (1,-1):
        base=log_state(state)
        shifted=from_log_state({n:base[n]+sign*h*getattr(v,n) for n in base},state)
        for n in FIELDS[6:]:setattr(shifted,n,getattr(state,n)+sign*h*getattr(v,n))
        f,_=compose_fine_hamiltonian(grid,shifted)
        pushed=to_canonical_tangent(shifted,f)
        answers.append(pushed)
        proper.append(combine_variations((1/observer_values(shifted,row)[2],pushed)))
    close(exact,combine_variations((1/(2*h),answers[0]),(-1/(2*h),answers[1])),atol=2e-7,rtol=1e-6)
    close(combine_variations((1/observer_values(state,row)[2],clocked)),
          combine_variations((1/(2*h),proper[0]),(-1/(2*h),proper[1])),atol=2e-7,rtol=1e-6)
    assert max(np.max(abs(getattr(exact,n)-getattr(clocked,n))) for n in CANONICAL_FIELDS)>1e-5


def test_constant_projection_preserves_fermions_and_is_complementary(example):
    _grid,state,u,_row=example;v=to_canonical_tangent(state,u)
    p,q=retained_projection(v),environment_projection(v)
    close(combine_variations((1,p),(1,q)),v,atol=1e-16)
    close(retained_projection(p),p,atol=1e-16)
    for n in CANONICAL_FIELDS[:6]:
        assert np.max(abs(getattr(p,n)-np.mean(getattr(v,n))))<1e-16
        assert abs(np.mean(getattr(q,n)))<1e-16
    for n in ('phi0','phi1'):
        assert_allclose(getattr(p,n),getattr(v,n),atol=0)
        assert_allclose(getattr(q,n),0,atol=0)


def _matrix_split(matrix,initial,time_value=.2,dt=.005):
    basis=np.array([[1,1],[1,-1]])/np.sqrt(2)
    nodal=basis@matrix@basis.T
    def vector(v):
        values={n:np.zeros(2) for n in CANONICAL_FIELDS[:6]}
        values['chi']=np.asarray(v)
        values.update(phi0=np.zeros((2,1),complex),phi1=np.zeros((2,1),complex))
        return CanonicalVariation(**values)
    def apply(v):return vector(nodal@v.chi)
    v=vector(basis@initial)
    parts=(retained_projection(v),environment_projection(v),combine_variations((0,v)))
    for _ in range(round(time_value/dt)):
        k=[]
        for w,j in [(0,None),(.5,0),(.5,1),(1,2)]:
            stage=parts if j is None else tuple(combine_variations((1,v),(dt*w,r)) for v,r in zip(parts,k[j]))
            k.append(causal_split_rates(apply,*stage))
        parts=tuple(combine_variations((1,parts[j]),*(
            (dt*w,k[i][j]) for i,w in enumerate((1/6,1/3,1/3,1/6)))) for j in range(3))
    return parts,basis


def test_closed_retarded_split_against_independent_matrix_exponential():
    matrix=np.array([[.1,.4],[-.3,-.2]]);initial=np.array([.2,.4]);t=.2
    (c,h,f),basis=_matrix_split(matrix,initial,time_value=t)
    expected=basis@expm(t*matrix)@initial
    assert_allclose(c.chi+h.chi+f.chi,expected,atol=2e-13)
    assert_allclose(h.chi,basis[:,1]*initial[1]*np.exp(matrix[1,1]*t),atol=2e-13)
    assert np.max(abs(f.chi))>1e-4
    assert_allclose(c.chi,np.mean(expected),atol=2e-13)


def test_zero_environment_coupling_has_no_feedback():
    (c,h,f),basis=_matrix_split(np.diag([.1,-.2]),np.array([.2,.4]))
    assert_allclose(f.chi,0,atol=1e-16)
    assert_allclose(c.chi,basis[:,0]*.2*np.exp(.02),atol=2e-13)
    assert_allclose(h.chi,basis[:,1]*.4*np.exp(-.04),atol=2e-13)


def test_packet_recovery_matches_independent_original_tangent_and_clock(example):
    grid,state,u,row=example
    packet=initial_packet(state,u);original=state.copy();direction=u.copy()
    for _ in range(8):
        packet=resolution_rk4_step(grid,packet,2e-4,row)
        original,direction=variational_rk4_step(grid,original,direction,2e-4)
    for n in FIELDS:
        assert_allclose(getattr(packet.state,n),getattr(original,n),atol=1e-12,rtol=1e-12)
        assert_allclose(getattr(packet.coordinate_tangent,n),getattr(direction,n),atol=1e-12,rtol=1e-12)
    rate,_=compose_fine_hamiltonian(grid,packet.state)
    close(reconstructed_tangent(packet),proper_time_reference(packet,rate,row),atol=1e-10)
    assert packet.clock>0 and abs(packet.clock_tangent)>1e-8
    for block in (packet.initial_drive,packet.feedback):
        assert_allclose(block.phi0,0,atol=0);assert_allclose(block.phi1,0,atol=0)


def test_proper_clock_formula_at_fixed_time_is_an_exact_family_derivative():
    initial,T=sp.symbols('initial T',positive=True)
    z=initial/(1-initial*T);u=sp.diff(z,initial)
    tau=-sp.log(1-initial*T);delta_tau=sp.diff(tau,initial)
    proper=u-z*delta_tau
    assert sp.cancel(proper-1/(1-initial*T))==0
    assert sp.cancel(proper/z-1/initial)==0


def test_fixed_receiver_clock_photon_response_against_exact_event_family():
    e,t0,delay,a,b,c,d=sp.symbols('e t0 delay a b c d',real=True,nonzero=True)
    shifted=sp.log(1+sp.exp(-b*e)*(sp.exp(a*t0)-1))/a
    logz=a*shifted+b*e-c*(shifted-delay)-d*e
    exact=sp.simplify(sp.diff(logz,e).subs(e,0))
    formula=b-d-(a-c)*b*(sp.exp(a*t0)-1)/(a*sp.exp(a*t0))
    assert sp.simplify(exact-formula)==0
    values={a:.4,b:.3,c:-.2,d:.1,t0:.3,delay:.125}
    assert fixed_receiver_clock_redshift(.3,.1,.4,-.2,
        .3*np.expm1(.4*.3)/.4,np.exp(.4*.3))==pytest.approx(float(exact.subs(values)))


def test_invalid_clock_and_step_domains(example):
    grid,state,u,row=example
    for dt in (0,-1,float('nan')):
        with pytest.raises(ValueError):resolution_rk4_step(grid,initial_packet(state,u),dt,row)
    with pytest.raises(ValueError):observer_values(state,np.zeros_like(row))
    with pytest.raises(ValueError):observer_values(state,np.zeros(len(row)+1))
    bad=state.copy();bad.r[0]=0
    with pytest.raises(ValueError):to_canonical_tangent(bad,u)
    with pytest.raises(ValueError):fixed_receiver_clock_redshift(0,0,0,0,0,0)
