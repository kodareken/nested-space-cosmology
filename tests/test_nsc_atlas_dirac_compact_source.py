"""Independent nonlinear mass-family controls; no physical wave/grid run."""
import numpy as np
import pytest
from scipy.integrate import quad,solve_ivp

from recursive_horizons.nsc_atlas_cap_wall_response import reference_coefficients
from recursive_horizons.nsc_atlas_dirac_radar import WallRadarChart
from recursive_horizons.nsc_atlas_dirac_compact_source import (
    CompactMassTangent,anchored_tortoise_mass_derivative,parent_mass_derivative)


@pytest.fixture(scope='module')
def model():
    return CompactMassTangent(WallRadarChart())


def full_nonlinear_family(chart,b,U,V,*,region):
    """Exact finite-b wall plus radial ray, without the optical tangent/K."""
    R0,k=.45,2.
    v0=float(reference_coefficients(R0)['v'])
    bm0=np.sqrt(v0*v0+1-2*b/R0-R0*R0)
    parent=R0*(1+v0*v0-(bm0-k*R0)**2)/2
    bm=lambda r:(parent-b)/(k*r*r)+3*r/4
    bp=lambda r:bm(r)-k*r
    bmr=lambda r:-2*(parent-b)/(k*r**3)+3/4
    bpr=lambda r:bmr(r)-k
    fm=lambda r:1-2*b/r-r*r
    fp=lambda r:1-2*parent/r
    fmr=lambda r:2*b/r**2-2*r
    fpr=lambda r:2*parent/r**2
    velocity=lambda r:np.sqrt(bm(r)**2-fm(r))
    acc=lambda r:bm(r)*bmr(r)-fmr(r)/2
    wall=solve_ivp(lambda t,y:velocity(y),(0.,chart.duration),
        np.array([R0],dtype=complex),method='DOP853',rtol=2e-12,atol=2e-14,
        max_step=chart.duration/30,dense_output=True)
    assert wall.success
    cap=region=='cap'
    beta,betar,f,fr=(bm,bmr,fm,fmr) if cap else (bp,bpr,fp,fpr)
    a=lambda u:1/(beta(wall.sol(u)[0])+velocity(wall.sol(u)[0]))
    initial=wall.sol(V)[0]
    ray=solve_ivp(lambda u,y:-a(u)*f(y)/2,(V,U),np.array([initial]),
        method='DOP853',rtol=2e-12,atol=2e-14,max_step=chart.duration/40)
    assert ray.success
    r=ray.y[0,-1]
    Ru,Rv=wall.sol(U)[0],wall.sol(V)[0]
    vu,vv=velocity(Ru),velocity(Rv)
    pu,qv=1/(beta(Ru)+vu),1/(beta(Rv)-vv)
    C=f(r)*pu*qv
    ru,rv=-f(r)*pu/2,f(r)*qv/2
    lu=-fr(r)*pu/2-(betar(Ru)*vu+acc(Ru))/(beta(Ru)+vu)
    lv=fr(r)*qv/2-(betar(Rv)*vv-acc(Rv))/(beta(Rv)-vv)
    W=np.sqrt(C)/r
    return dict(parent=parent,radius=r,C=C,radius_u=ru,radius_v=rv,
        log_C_u=lu,log_C_v=lv,potential=W,
        potential_u=W*(lu/2-ru/r),potential_v=W*(lv/2-rv/r))


def test_initial_mass_choice_preserves_physical_wall_speed(model):
    assert parent_mass_derivative()==pytest.approx(2755/8587,abs=1e-15)
    wall=model.mean.wall_state(0.)
    assert wall['displacement']==0.
    assert abs(wall['displacement_rate'])<2e-15
    assert wall['cap_mass']==1.
    assert wall['exterior_mass']==parent_mass_derivative()
    # Keeping the parent mass fixed would violate that initial speed.
    g=reference_coefficients(.45)
    wrong=-g['beta_plus']/(2*.45**2*g['v'])
    assert abs(wrong)>.5


def test_anchored_tortoise_mass_primitive_against_independent_quadrature():
    for side,r in (('cap',.43),('cap',.53),('exterior',.51)):
        F=(lambda x:1-x*x) if side=='cap' else (lambda x:1-.8/x)
        exact,error=quad(lambda x:2/(x*F(x)**2),.45,r,epsabs=2e-13)
        assert error<2e-12
        assert anchored_tortoise_mass_derivative(.45,region=side)==0.
        assert anchored_tortoise_mass_derivative(r,region=side)==pytest.approx(exact,abs=2e-14)


def test_optical_metric_and_Dirac_tangent_against_full_nonlinear_family(model):
    S=model.chart.duration
    for U,V,side in ((.8*S,.7*S,'cap'),(.2*S,.8*S,'exterior')):
        tangent=model.fields(U,V,region=side)
        step=1e-7
        exact=full_nonlinear_family(model.chart,step*1j,U,V,region=side)
        assert exact['parent'].imag/step==pytest.approx(parent_mass_derivative(),abs=3e-12)
        background=tangent['background']
        assert exact['radius'].real==pytest.approx(background['radius'],abs=3e-11)
        assert exact['C'].real==pytest.approx(background['conformal_factor'],abs=3e-10)
        pairs=(('radius_change','radius'),('radius_u_change','radius_u'),
            ('radius_v_change','radius_v'),('log_C_u_change','log_C_u'),
            ('log_C_v_change','log_C_v'),('potential_change','potential'),
            ('potential_u_change','potential_u'),('potential_v_change','potential_v'))
        for measured,key in pairs:
            assert tangent[measured]==pytest.approx(exact[key].imag/step,rel=2e-8,abs=2e-7)
        assert tangent['log_C_change']==pytest.approx(
            exact['C'].imag/(step*background['conformal_factor']),rel=2e-8,abs=2e-7)


def test_wall_matching_and_signed_beta_crossing(model):
    for r in (.45,.5,.16**(1/3),.55):
        t=float(model.chart.maps(r)[0])
        wall=model.mean.wall_state(t)
        sides=[model.fields(t,t,region=side) for side in ('cap','exterior')]
        for row in sides:
            assert row['radius_change']==pytest.approx(wall['displacement'],abs=2e-10)
            assert abs(row['log_C_change'])<2e-10
            assert row['radius_u_change']+row['radius_v_change']==pytest.approx(
                wall['displacement_rate'],abs=2e-9)
        for key in ('potential_change','potential_u_change','potential_v_change'):
            assert sides[0][key]==pytest.approx(sides[1][key],abs=2e-8)


def test_same_metric_from_existing_EF_reconstruction_and_bounded_queries(model):
    S=model.chart.duration
    for U,V in ((.8*S,.7*S),(.2*S,.8*S)):
        optical=model.fields(U,V)
        ef=model.mean.clock_response(U,V,ray_rtol=2e-11,ray_atol=2e-13)
        assert ef['radius_change']==pytest.approx(optical['radius_change'],abs=2e-9)
        assert ef['log_radar_C_change']==pytest.approx(optical['log_C_change'],abs=2e-9)
        assert ef['mass_path_gap']==0.
    u,v=np.meshgrid(np.array([.4,.8])*S,np.array([.3,.7])*S,indexing='ij')
    fields=model.fields(u,v)
    assert fields['radius_change'].shape==(2,2)
    with pytest.raises(ValueError,match='null clock labels'):
        model.fields(-.1*S,.5*S)
    with pytest.raises(ValueError,match='angular'):
        model.fields(.4*S,.5*S,angular_k=True)
