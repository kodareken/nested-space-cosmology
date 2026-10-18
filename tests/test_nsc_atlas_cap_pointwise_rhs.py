"""Independent polar-field and reduced-coordinate scalar-equation controls."""
import math

import pytest
import sympy as s
import numpy as np

from recursive_horizons.nsc_atlas_cap_pointwise_rhs import scalar_rhs,incoming_trace_defect


def test_reduced_zero_coupling_recovers_both_flat_characteristic_equations():
    # For arbitrary local chi jets, the conformal de Sitter reduction is
    # chi_T=k, k_T=chi_zz; this checks all curvature and potential terms.
    r=.53;chi=.17;k=-.23;cz=.8;czz=-1.7;kz=.9
    row=scalar_rhs(radius=r,chi=chi,k=k,chi_z=cz,chi_zz=czz,k_z=kz,
        mass=.3,lapse=0.,alpha=0.,potential_difference=(chi/r)**2,potential_derivative=2*chi/r)
    assert row['chi_T']==pytest.approx(k,abs=2e-15)
    assert row['k_T']==pytest.approx(czz,abs=2e-15)
    assert row['chi_T_z']==pytest.approx(kz,abs=2e-15)


def test_independent_polar_field_equation_at_nonzero_coupling():
    r=s.symbols('r',positive=True)
    f=s.Rational(1,7)+r*r/3-r**4/9
    Pi=s.Rational(1,5)-r*r/11
    A0=1-r*r
    chi=r*f;k=A0*r*Pi
    jets=[chi,A0*s.diff(chi,r),A0*s.diff(A0*s.diff(chi,r),r),k,A0*s.diff(k,r)]
    x=.51;chi,cz,czz,k,kz=[float(v.subs(r,x)) for v in jets]
    fv=float(f.subs(r,x));fr=float(s.diff(f,r).subs(r,x));frr=float(s.diff(f,r,2).subs(r,x))
    pv=float(Pi.subs(r,x));pr=float(s.diff(Pi,r).subs(r,x))
    alpha=.004;mass=.2;lapse=.013;q=math.sqrt(alpha/6)*fv
    P=fv*fv*(math.sinh(q)/q)**2*(math.cosh(q)**2+1)/2
    U=2*fv*(math.sinh(q)/q)*math.cosh(q)**3
    row=scalar_rhs(radius=x,chi=chi,k=k,chi_z=cz,chi_zz=czz,k_z=kz,
        mass=mass,lapse=lapse,alpha=alpha,potential_difference=P,potential_derivative=U)
    A=1-x*x-2*alpha*mass/x;mu=math.exp(lapse);B=A*mu
    rho=x*x*P/2+A*x*x*(pv*pv+fr*fr)/4
    Ar=-2*x-2*alpha*(rho/x-mass/x**2)
    dr=alpha*x*(pv*pv+fr*fr)/2
    Br=mu*(Ar+A*dr)
    assert row['chi_T']==pytest.approx(x*B*pv,abs=2e-15)
    # r^2 Pi_T = partial_r(r^2 B f_r)-r^2 exp(delta) U.
    assert row['k_T']==pytest.approx((1-x*x)*x*(B*frr+(2*B/x+Br)*fr-mu*U),abs=2e-15)
    assert row['chi_T_z']==pytest.approx((1-x*x)*((B+x*Br)*pv+x*B*pr),abs=2e-15)


def test_incoming_defect_is_the_direct_derivative_of_reconstructed_trace():
    z,t=s.symbols('z t',real=True)
    a=1+z+2*z**3;b=2-z*z
    r=s.tanh(z);A0=1-r*r
    # Difference family delta chi=t*a(z), delta k=t*b(z), at fixed z.
    delta_trace=(t*b/A0+s.diff(t*a,z)/A0-t*a/r)/r
    where=.55;R=math.tanh(where)
    v=float(a.subs(z,where));vz=float(s.diff(a,z).subs(z,where));d=float(b.subs(z,where))
    row=incoming_trace_defect(radius=R,continuum=dict(chi_T=0.,k_T=0.,chi_T_z=0.),
        represented_chi_T=v,represented_k_T=d,represented_chi_T_z=vz)
    assert row['incoming_cap_rate_defect']==pytest.approx(float(s.diff(delta_trace,t).subs(z,where)),abs=1e-13)
    channels=[row[k] for k in ('incoming_density_channel','incoming_velocity_slope_channel',
                               'incoming_velocity_value_channel')]
    assert abs(row['incoming_cap_rate_defect']-sum(channels))<=2*math.ulp(sum(map(abs,channels)))


@pytest.mark.parametrize('r',[0.,1.,-1.])
def test_invalid_cap_or_centre_needs_a_different_regular_operator(r):
    with pytest.raises(ValueError):
        scalar_rhs(radius=r,chi=0.,k=0.,chi_z=0.,chi_zz=0.,k_z=0.,mass=0.,lapse=0.,alpha=0.,
            potential_difference=0.,potential_derivative=0.)


def test_projected_equation_recovers_boundary_load_and_rejects_wrong_bulk_potential(monkeypatch):
    import diagnose_nsc_atlas_wall_rhs as diagnostic
    from recursive_horizons.nsc_atlas_cap_background_jets import RetainedCapBackground
    from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
    from recursive_horizons.nsc_atlas_cap_fixed_reference import fixed_cap_reference
    from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
    from recursive_horizons import nsc_atlas_cap_variational_reference as radial
    q=np.array([-.8,0.,.8]);positive=np.maximum(q,0.)
    bg=RetainedCapBackground(q,.1*positive**4,.4*positive**3,1.2*positive**2)
    budget=Budget(maximum_evaluations=10_000_000,cpu_seconds=30.)
    basis=FixedCapBasis([0,.05,.85],budget=budget);zero=np.zeros(basis.dimension)
    state=dict(cap=dict(chi=zero,canonical_density=zero))
    Z=math.atanh(.45);jet=bg.reduced_jets(.2,Z,budget=budget)
    incoming=(jet['p']+jet['h']+.17)/.45
    original=radial.solve_ivp;solutions=[]
    def capture(*a,**kw):
        result=original(*a,**kw);solutions.append(result);return result
    with monkeypatch.context() as routed:
        routed.setattr(radial,'solve_ivp',capture)
        cap=fixed_cap_reference(basis,zero,zero,Z,alpha=.001,T_x=1.,radius_rate=0.,
            incoming_w_plus=incoming,budget=budget,background=bg,background_time=.2,
            background_panel_sides=True,rtol=2.5e-14,atol=1e-15)
    result=diagnostic.projection_check(basis,state,cap,bg,solutions,budget)
    assert max(result['normalized_load_differences'])<1e-8
    assert max(map(abs,result['expected_density_boundary_load']))>1e-3
    true_rhs=diagnostic.scalar_rhs
    def wrong_potential(**kwargs):
        kwargs['potential_derivative']*=.5
        return true_rhs(**kwargs)
    with monkeypatch.context() as altered:
        altered.setattr(diagnostic,'scalar_rhs',wrong_potential)
        with pytest.raises(ValueError,match='projected scalar equation'):
            diagnostic.projection_check(basis,state,cap,bg,solutions,budget)
