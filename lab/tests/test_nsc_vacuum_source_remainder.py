"""The covariance expansion has an equation remainder, not an order-difference error."""
import numpy as np
import pytest
import sympy as sp
from flint import arb,arb_series,ctx
from scipy.integrate import solve_ivp

from recursive_horizons.nsc_massive_jost_phase_bound import _series_context
from recursive_horizons.nsc_metric_horizon_frame import metric_horizon_frame
from recursive_horizons.nsc_paired_horizon_preparation import source_covariance
from recursive_horizons.nsc_vacuum_source_remainder import VacuumSourceExpansion,coefficient_jets

HINT=1.9006916054701435


def model(order=8):
    with ctx.workprec(192):
        return VacuumSourceExpansion(HINT,arb.pi()/2,arb(5).sqrt(),order=order)


def test_exact_nonconstant_manufactured_ode_and_normalization_residuals():
    t,e=sp.symbols('t e',positive=True)
    H=2+t*t/10;L=1+t*t/3+sp.I*(sp.Rational(2,5)+t*t/7)
    W=[sp.Integer(0)];Z=[sp.Integer(1)]
    derivative=lambda f:sp.expand(t*sp.diff(f,t)/2)
    for n in range(1,10):
        W.append(sp.expand(sp.I*H*derivative(W[-1])/2+t*L*Z[-1]))
        Z.append(sp.expand(-sum((sp.conjugate(W[j])*W[n-j]+Z[j]*Z[n-j]
                                for j in range(1,n)))/2))
    for M in (1,2,4,8):
        w=sum(W[j]/e**j for j in range(1,M+1))
        z=sum(Z[j]/e**j for j in range(M+1))
        rw=derivative(w)+2*sp.I*e/H*w-2*sp.I*t*L/H*z
        rz=derivative(z)+2*t/H*sp.im(sp.expand(sp.conjugate(L)*w))
        assert sp.cancel(rw+2*sp.I*W[M+1]/(H*e**M))==0
        assert sp.cancel(rz)==0
        norm=sp.expand(w*sp.conjugate(w)+z*z-1)
        for j in range(M+1):assert norm.coeff(e,-j)==0
    # Independent t-coordinate recurrence versus the implemented delta jets.
    with _series_context(192,9):
        x=arb_series([arb(1)/5,1],prec=10)
        ps,zs=coefficient_jets(x,2+x/10,(1+x/3,arb(2)/5+x/7),8)
        for j in range(1,10):
            exact=sp.expand(W[j]/t).subs(t,sp.sqrt(sp.Rational(1,5)))
            for k,part in enumerate((sp.re(exact),sp.im(exact))):
                value=arb(int(sp.numer(part)))/int(sp.denom(part))
                assert ps[j][k][0].contains(value)
            part=Z[j].subs(t,sp.sqrt(sp.Rational(1,5)))
            value=arb(int(sp.numer(part)))/int(sp.denom(part))
            assert zs[j][0].contains(value)


class ConstantFixture(VacuumSourceExpansion):
    def target_distance(self,rho):return arb(1)
    def jets(self,delta):
        with _series_context(self.bits,self.order+1):
            x=arb_series([delta,1],prec=self.order+2)
            H=arb_series([1],prec=self.order+2)
            ps,zs=coefficient_jets(x,H,(H,0*H),self.order)
            return H,ps,zs


def test_defect_integral_needs_its_factor_four():
    fixture=ConstantFixture(HINT,1,0,order=1)
    C=fixture.remainder_constant(1,cells=1)
    assert C==1
    # At E=.1 the first approximation has n=(10,0,1). Every exact Bloch
    # vector has norm one, so the error is at least sqrt(101)-1.
    lower=arb(101).sqrt()-1
    assert C/arb('.1')>lower
    assert C/(4*arb('.1'))<lower


def test_uniform_remainder_and_horizon_initial_limit():
    with ctx.workprec(192):
        m=model();_,ps,zs=m.jets(arb(0))
        for j in range(1,m.order+1):assert zs[j][0].contains(0)
        C=m.remainder_constant(1.03,cells=32)
        assert C>0
        bound=m.covariance_error(C,160,kappa=.23832579963401956,omega=3.9)
        assert bound['source_operator_error']<arb('3e-14')
        double=m.covariance_error(C,320,kappa=.23832579963401956,omega=3.9)
        assert double['vacuum_operator_error']<bound['vacuum_operator_error']/255
        zero=VacuumSourceExpansion(HINT,0,0,order=4)
        assert zero.remainder_constant(1.03,cells=4)==0
        delta=m.target_distance(1.03)
        _,short,_=m.jets(delta)
        _,long,_=m.jets(delta,remainder_jet_order=4)
        for k in (0,1):
            assert long[-1][k].prec==5
            assert long[-1][k][0].overlaps(short[-1][k][0])


def test_full_metric_vacuum_agrees_with_independent_ode():
    m=model(order=8);E=16.;rho=1.03;delta0=1e-8
    frame=metric_horizon_frame(E,np.pi/2,np.sqrt(5),HINT,order=32,bits=192)
    with ctx.workprec(192):
        F,tail=frame.evaluate(arb(delta0),interior=True)
        assert tail<arb('1e-20')
        column=np.array([complex(F[j][0].mid()) for j in range(2)])
        Q=np.outer(column,column.conj())
        initial=np.array([2*Q[0,1].real,-2*Q[0,1].imag,(Q[0,0]-Q[1,1]).real])
        target=float(m.target_distance(rho).mid())
        approx=np.array([float(v.mid()) for v in m.approximate_bloch(rho,E)])
        bound=float(m.remainder_constant(rho,cells=64)/arb(E)**m.order)
    q=float(m.q);cc=1.5*np.sin(2*q)+.5*np.cos(2*q);ss=1.5*np.cos(2*q)-.5*np.sin(2*q)
    def rhs(y,n):
        d=np.exp(y);u=-d
        # Stable closed metric, independently of the Taylor jet implementation.
        H=3+2*cc*np.sin(u)**2/u-2*ss*np.sinc(2*u/np.pi)
        scale=np.sqrt(d/H)
        h=np.array([-np.pi/2/np.sin(q-d)*scale,np.sqrt(5)*scale,-E/H])
        return 2*np.cross(h,n)
    run=solve_ivp(rhs,(np.log(delta0),np.log(target)),initial,method='DOP853',rtol=2e-11,atol=2e-13)
    assert run.success
    error=np.linalg.norm(run.y[:,-1]-approx)
    assert error<bound
    wrong=approx.copy();wrong[:2]*=-1
    assert np.linalg.norm(run.y[:,-1]-wrong)>bound


def test_source_thermal_and_coherent_bound_includes_open_infinity_port():
    m=model();E=.3;kappa=.23832579963401956;omega=3.9
    C=source_covariance(E,kappa,omega,0.)
    R=.6*np.exp(.73j);T=1-abs(R)**2
    sewing=np.array([[0,1,0],[R,0,np.sqrt(T)]],complex)
    Q=sewing@C@sewing.conj().T
    error=np.linalg.norm(Q-np.diag([1.,0.]),2)
    bound=m.covariance_error(0,E,kappa=kappa,omega=omega)['thermal_coherent_error']
    assert error<=float(bound)
    assert error>float(np.exp(-2*np.pi*E/kappa))
    assert error>float(np.exp(-2*np.pi*E/kappa)+np.exp(-np.pi*E/kappa))
    # With no transmitted infinity contribution, omitting coherence alone fails.
    C=source_covariance(E,kappa,.1,0.)
    sewing=np.array([[0,1,0],[1,0,0]],complex)
    error=np.linalg.norm(sewing@C@sewing.conj().T-np.diag([1.,0.]),2)
    no_coherence=np.exp(-2*np.pi*E/kappa)+np.exp(-2*np.pi*E/(.1*kappa))
    assert error>no_coherence
    assert error<=float(m.covariance_error(0,E,kappa=kappa,omega=.1)['thermal_coherent_error'])


def test_original_comparison_keeps_signed_gram_error():
    m=VacuumSourceExpansion(HINT,0,0)
    A=np.array([[1,0,0],[0,2,0]],complex);C=np.diag([1.,0.,0.])
    kwargs={'kappa':.23832579963401956,'omega':3.9}
    positive=m.compare_original(1.03,160,A,C,0,**kwargs)
    negative=m.compare_original(1.03,160,np.diag([1,-1])@A,np.eye(3)-C,0,
                                negative_partner=True,**kwargs)
    assert positive['polynomial_to_stored_operator_distance']==0
    assert negative['polynomial_to_stored_operator_distance']==3
    assert negative['original_source_operator_error']>=3


def test_original_high_energy_inventory_has_complete_selected_row_bounds():
    import runpy
    from pathlib import Path
    from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
    root=Path(__file__).resolve().parents[1]
    driver=runpy.run_path(str(root/'scripts/derive_nsc_vacuum_source_remainder.py'))
    record=driver['calculate']();coverage=record['archived_comparison']
    assert coverage['signed_rows']==len(coverage['rows'])==768
    assert coverage['positive_rows']==384
    assert coverage['remaining_signed_rows_in_family']==400
    assert len({(r['panel'],r['row'],r['energy_sign']) for r in coverage['rows']})==768
    assert restored_upper(coverage['maximum_original_source_operator_error_upper'])<arb('1.1e-12')
    assert coverage['source_quadrature_error_included'] is False
    assert record['physical_upstream_budget_component'] is None
    assert record['changed_history_C_M'] is None
    assert record['physical_local_gate']=='OPEN'


def test_invalid_domains_and_precision_are_rejected():
    before=ctx.prec,ctx.cap
    for kwargs in ({'order':0},{'mass':arb(0,'.1')},{'mass':None},
                   {'angular':None},{'metric_terms':3},{'bits':32}):
        args={'mass':1,'angular':2};args.update(kwargs)
        with pytest.raises(ValueError):VacuumSourceExpansion(HINT,**args)
    m=model()
    for rho in (0,2):
        with pytest.raises(ValueError):m.target_distance(rho)
    with pytest.raises(ValueError):m.remainder_constant(1.03,cells=0)
    with pytest.raises(ValueError):m.approximate_bloch(1.03,0)
    with pytest.raises(ValueError):m.covariance_error(None,160,kappa=.2,omega=3.)
    assert (ctx.prec,ctx.cap)==before
