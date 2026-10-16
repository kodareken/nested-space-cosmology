"""Exact moving-wall scalar and mass residuals of the fixed-space cap port."""
import sympy as s
import pytest

from control_budget import ControlBudget

CONTROLS=ControlBudget(maximum_evaluations=54000,cpu_seconds=20.)
R,A,H,mu,rootD=s.symbols('R A H mu rootD',positive=True)
P,Q,dp,dh=s.symbols('P Q dp dh',real=True)
beta=(H+A/H)/2
v=(H-A/H)/2
B=A*mu
p=P+dp
h=Q+dh
incoming=(P+Q)/R
outgoing=(p-h)/R
matched_p=R*(incoming+outgoing)/2
matched_h=R*(incoming-outgoing)/2


@pytest.fixture
def budget():
    with CONTROLS.active() as value: yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('fixed cap hybrid identities:',CONTROLS.budget.evaluations,
          'charged algebra values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def test_common_scalar_tangent_keeps_both_trace_gaps(budget):
    budget.check(8000)
    # Independently differentiate f=chi/r on the moving wall, with
    # chi_s=T_s*V_W+z_s*chi_z, T_s=beta/B and z_s=v/A0.
    chi,chi_z,A0=s.symbols('chi chi_z A0',real=True,nonzero=True)
    primitive_rate=(beta/B)*(B*P)+(v/A0)*chi_z
    cap_tangent=primitive_rate/R-chi*v/R**2
    cap_tangent=s.expand(cap_tangent.subs(chi_z,A0*(h+chi/R)))
    hybrid_tangent=(beta*matched_p+v*matched_h)/R
    residual=s.factor(cap_tangent-hybrid_tangent)
    assert s.simplify(residual-(-A*dp/H+H*dh)/(2*R))==0
    assert s.simplify(matched_p-P-(dp-dh)/2)==0
    assert s.simplify(matched_h-Q+(dp-dh)/2)==0
    assert s.simplify((p-matched_p)-(dp+dh)/2)==0
    assert s.simplify((h-matched_h)-(dp+dh)/2)==0


def test_frame_conversion_maps_the_canonical_tangent_residual(budget):
    budget.check(8000)
    f,qscale,fcap,ftarget=s.symbols('f qscale fcap ftarget',real=True,nonzero=True)
    # u=sinh(qscale*f)/qscale has derivative cosh(qscale*f)=sqrt(D).
    u=s.sinh(qscale*f)/qscale
    cap_c_rate=v*u+R*s.diff(u,f)*fcap
    target_c_rate=v*u+R*s.diff(u,f)*ftarget
    difference=s.expand(cap_c_rate-target_c_rate)
    canonical_gap=(-A*dp/H+H*dh)/(2*R)
    assert s.simplify(difference-R*s.cosh(qscale*f)*(fcap-ftarget))==0
    assert s.simplify((R*rootD*canonical_gap)-rootD*(-A*dp/H+H*dh)/2)==0
    assert s.limit(s.diff(u,f),qscale,0)==1


def test_temporal_mass_residual_distinguishes_actual_and_matched_pairs(budget):
    budget.check(8000)
    # Remove the same moving-domain rho*v term from the total mass rate.
    rho=s.symbols('rho',real=True)
    weak_T=A*A*mu*P*Q/2
    total_s=(beta/B)*weak_T+v*rho
    recovered_T=s.simplify((total_s-v*rho)/(beta/B))
    actual_T=A*A*mu*p*h/2
    matched_T=A*A*mu*matched_p*matched_h/2
    assert s.simplify(recovered_T-actual_T+A*A*mu*(P*dh+Q*dp+dp*dh)/2)==0
    a=(dp-dh)/2
    assert s.simplify(recovered_T-matched_T-A*A*mu*(a*(P-Q)+a*a)/2)==0


def test_weak_outgoing_pair_does_not_erase_moving_primitive_gap(budget):
    budget.check(6000)
    weak_outgoing=(P-Q)/R
    target_p=R*(incoming+weak_outgoing)/2
    target_h=R*(incoming-weak_outgoing)/2
    cap_tangent=(beta*P+v*h)/R
    weak_target=(beta*target_p+v*target_h)/R
    assert s.simplify(cap_tangent-weak_target-v*dh/R)==0
    assert s.simplify((cap_tangent-weak_target).subs(H,s.sqrt(A)))==0


def test_zero_weak_balance_alone_does_not_close_hybrid_scalar_or_actual_mass(budget):
    budget.check(6000)
    fixed={A:1,H:2,R:1,mu:1,P:1,Q:2,dp:1,dh:0}
    cap_tangent=(beta*P+v*h)/R
    hybrid_tangent=(beta*matched_p+v*matched_h)/R
    assert s.simplify((cap_tangent-hybrid_tangent).subs(fixed))==-s.Rational(1,4)
    weak_T=A*A*mu*P*Q/2
    actual_T=A*A*mu*p*h/2
    assert s.simplify((weak_T-actual_T).subs(fixed))==-1
    # The cap port's own mass identity remains exact in this counterexample.
    Vw=B*P;C=B*R*incoming
    assert s.simplify(weak_T-Vw*(C-Vw)/(2*mu))==0


def test_physical_mass_dictionary_keeps_the_moving_vacuum_density(budget):
    budget.check(4000)
    alpha,N,NT,rho=s.symbols('alpha N NT rho',real=True)
    # A0-2*alpha*N/R = 1-2*m/R determines the mass dictionary.
    mass=R**3/2+alpha*N
    assert s.simplify(1-2*mass/R-(1-R*R-2*alpha*N/R))==0
    rescaled_s=(beta/B)*NT+v*rho
    mass_s=s.diff(mass,R)*v+s.diff(mass,N)*rescaled_s
    physical_rho=3*R*R/2+alpha*rho
    recovered_T=s.simplify((mass_s-v*physical_rho)/(beta/B))
    assert s.simplify(recovered_T-alpha*NT)==0
    assert s.simplify(recovered_T-alpha*A*A*mu*p*h/2
                      -alpha*(NT-A*A*mu*p*h/2))==0
    assert s.simplify(mass_s.subs(alpha,0)-3*R*R*v/2)==0


def test_velocity_and_incoming_channels_from_independent_moving_trace(budget):
    budget.check(6000)
    rx,tx=s.symbols('rx tx',real=True)
    # Weak primitive velocity is B*P; the physical radial trace is h.
    cap_tangent=(tx*B*P+rx*h)/R
    hybrid_tangent=(tx*B*matched_p+rx*matched_h)/R
    velocity_gap=B*(P-p)
    incoming_gap=(p+h)/R-incoming
    velocity_channel=tx*velocity_gap/R
    incoming_channel=(rx+tx*B)*incoming_gap/2
    assert s.simplify(cap_tangent-hybrid_tangent-velocity_channel-incoming_channel)==0
    fixed={A:1,mu:1,R:1,rx:1,tx:2}
    assert s.simplify(velocity_channel.subs({**fixed,dp:1,dh:-1}))==-2
    assert s.simplify(incoming_channel.subs({**fixed,dp:1,dh:-1}))==0
    assert s.simplify(velocity_channel.subs({**fixed,dp:0,dh:1}))==0
    assert s.simplify(incoming_channel.subs({**fixed,dp:0,dh:1}))==s.Rational(3,2)


def test_mass_sensitivity_yields_only_metric_correction_projection(budget):
    budget.check(8000)
    r,A0,alpha,N,ps,hs,potential,phi=s.symbols('r A0 alpha N ps hs potential phi',positive=True)
    actual_A=A0-2*alpha*N/r
    rho=r*r*potential/2+actual_A*(ps*ps+hs*hs)/4
    damping=alpha*(ps*ps+hs*hs)/(2*r)
    assert s.simplify(s.diff(rho,N)+damping)==0
    assert s.simplify(s.diff(rho,ps)*phi/A0-actual_A*ps*phi/(2*A0))==0
    delta=s.Function('delta');radius=s.symbols('radius',positive=True)
    weight=s.exp(delta(r)-delta(radius))
    assert s.simplify(s.diff(weight,r)-s.diff(delta(r),r)*weight)==0
    assert weight.subs(r,radius)==1
    assert s.simplify(s.exp(delta(radius))*weight-s.exp(delta(r)))==0
    # Independent exact one-function Galerkin integral. The nonconstant
    # b below is a projection control, not an asserted metric solution.
    z,Z,beta,density=s.symbols('z Z beta density',positive=True)
    background_k=2*z**3;total_k=background_k+density*z;b=1+beta*z*z
    gram=s.integrate(z*z,(z,0,Z))
    gradient=s.integrate(z*b*total_k/2,(z,0,Z))
    projected_background=s.integrate(z*background_k,(z,0,Z))/gram
    weak_velocity=Z*(2*gradient/gram-projected_background)+background_k.subs(z,Z)
    physical_velocity=(b*total_k).subs(z,Z)
    defect=s.factor(weak_velocity-physical_velocity)
    independent=-beta*(8*Z**5/s.Integer(7)+2*density*Z**3/s.Integer(5))
    assert s.simplify(defect-independent)==0
    assert defect.subs(beta,0)==0
    assert defect.subs({beta:1,Z:1,density:1})==-s.Rational(54,35)
