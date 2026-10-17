"""Independent potential, wall-difference and manufactured envelope controls."""
from fractions import Fraction as Q

from flint import arb,ctx
import mpmath as mp
import pytest

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_coupled_tube_bound import (
    DECLARED_BOX,coupled_vacuum_tube_bound,potential_majorants,wall_difference_bounds)
from recursive_horizons.nsc_atlas_reference_gauge_bound import _ball,_record

CONTROLS=ControlBudget(maximum_evaluations=40000,cpu_seconds=20.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('coupled tube controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def _mpupper(value):
    n,e=value.upper().man_exp()
    return mp.mpf(str(n))*mp.mpf(2)**int(e)


def test_hyperbolic_majorants_and_regular_reaction_limit(budget):
    budget.check(2000)
    with mp.workdps(80),ctx.workprec(160):
        f,al=mp.mpf(1)/5,mp.mpf(3)/4;z=mp.sqrt(al/6)*f
        _,P,U=potential_majorants(_ball(Q(1,5)),_ball(Q(3,4)))
        exactP=f*f*(mp.sinh(z)/z)**2*(mp.cosh(z)**2+1)/2
        exactU=(2*f*mp.sinh(z)/z*mp.cosh(z)**3-2*f)/al
        assert exactP<=_mpupper(P) and exactU<=_mpupper(U)
        _,P0,U0=potential_majorants(_ball(Q(1,5)),_ball(0))
        assert f*f<=_mpupper(P0) and 5*f**3/9<=_mpupper(U0)


def _wall_case(*,alpha,field,outgoing,scalar_Q,mass,lapse,rx_change,budget):
    budget.check(8000)
    with mp.workdps(80),ctx.workprec(160):
        r=mp.mpf(1)/2;A0=1-r*r
        bm0=1/(5*r*r)+3*r/4;bp0=1/(5*r*r)-5*r/4
        v0=mp.sqrt(bp0*bp0+mp.mpf(4)/(5*r)-1)
        H0=bm0+v0;rx0=(H0-2*r)/2
        al,f,em,q,n,d,drx=map(mp.mpf,(alpha,field,outgoing,scalar_Q,mass,lapse,rx_change))
        a=mp.sqrt(al/6);root=mp.cosh(a*f);u=mp.sinh(a*f)/a if al else f
        Sw=2*root**3;gap=u-f*root;rx=rx0+drx
        A=A0-2*al*n/r;H=2*rx+Sw*r;Agap=A/H
        bbar=(2*rx+Agap)/2;bminus=(H+Agap)/2
        cx=root*q+gap*rx;cm=Agap*(root*em/A0-gap);nc=(2*cx-cm)/2
        k=-3*root/4-al*u*(nc-u*bbar)/(2*r*root**2)-2*bbar/r-Sw/2
        k0=-mp.mpf(2)/(5*r**3)-mp.mpf(5)/4
        pr=(Sw*r-Agap)/2;pr0=(2*r-A0/H0)/2
        Lx=mp.mpf('0.00001') if al else mp.mpf(0)
        pL_difference=abs(-2*(k-k0)+Lx)
        Tx=bminus/(A*mp.exp(d));Tx0=bm0/A0
        # All nonzero prescribed changes are strictly inside these envelopes.
        rootmax=mp.cosh(a*(abs(f)+mp.mpf('.001')));Swmax=2*rootmax**3
        droot=al*(abs(f)+mp.mpf('.001'))**2*rootmax/12
        nmax=abs(n)+mp.mpf('.00001');dmax=abs(d)+mp.mpf('.00001')
        uMax=(abs(f)+mp.mpf('.001'))*rootmax
        betaBarMax=(2*(rx0+abs(drx))+A0/H0)/2+1
        parameters=dict(eps_r=0,eps_L=0,Verror=abs(drx),Lxerror=abs(Lx),alpha=al,N=nmax,
            mu_error=dmax*mp.exp(dmax),Rlo=r,Rhi=r,R0lo=r,R0hi=r,
            Amin=A0-2*al*nmax/r,Hmin=H0,Upmin=1,Upmax=1,Vpmin=1,Vpmax=1,Wxmin=1,
            rx0max=rx0,Sw=Swmax,root_change=droot,Sw_change=6*rootmax**2*droot,
            u_wall=uMax,nc=abs(nc)+1,beta_bar_max=betaBarMax)
        bound=wall_difference_bounds(**{k:arb(str(v)) for k,v in parameters.items()})
        for key,value in dict(r_U=abs(pr-pr0),L_U=pL_difference,T_x=abs(Tx-Tx0),k_plus=abs(k-k0)).items():
            assert value<=_mpupper(bound[key]),key
        return abs(k-k0),_mpupper(bound['k_plus'])


def test_wall_curvature_keeps_the_factor_two(budget):
    exact,bound=_wall_case(alpha='0',field='0',outgoing='0',scalar_Q='0',mass='0',
                           lapse='0',rx_change='.0001',budget=budget)
    # A missing factor two in -2*beta_bar/R makes this envelope too small.
    assert exact>bound/2


def test_wall_mass_lapse_and_scalar_cross_terms(budget):
    _wall_case(alpha='.1',field='.002',outgoing='.001',scalar_Q='.003',mass='.00001',
                lapse='.00001',rx_change='0',budget=budget)


def test_manufactured_envelope_wiring_has_all_norms_and_guards(budget):
    # Broad synthetic coefficient inputs exercise the envelope algebra only.
    # They are not an authentication or calculation of the original source cell.
    with ctx.workprec(160):
        gauge=dict(inputs=dict(left='21/50',right='420001/1000000',radius_lower='47/100',radius_upper='12/25'),
                   wall=dict(W=_record(_ball(1)),W_x=_record(_ball(1).union(_ball(4))),
                             C_wall=_record(_ball(Q(1,2)).union(_ball(1)))))
        source=dict(duration='1/1000000',alpha_max='1/1000',forcing={k:'0' for k in
            ('r','r_x','log_C','log_C_x','c','c_x','X','Q','r_times_r_x_minus_x')})
        result=coupled_vacuum_tube_bound(source,gauge,z_star=Q(9,20),box=DECLARED_BOX,budget=budget)
    assert result['strict']
    assert set(result['image'])==set(DECLARED_BOX)
    assert all(result['margin_checks'].values())
    assert 'construction' in result['scope']
