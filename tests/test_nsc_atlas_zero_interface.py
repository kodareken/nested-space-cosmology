"""Independent identities for field zeros, conformal interfaces and winding.

These tests derive equations directly rather than inspect interpretation flags
in the adapter. A reciprocal scale-factor extension is a declared comparison,
and oscillator damping alone is not a microscopic entropy-production theorem.
"""
import sympy as s


def test_damped_offset_lobes_against_segmented_quadrature():
    import importlib.util
    import math
    from pathlib import Path
    from scipy.integrate import quad
    path=Path(__file__).resolve().parents[1]/'scripts/check_nsc_atlas_flow_cosmology.py'
    spec=importlib.util.spec_from_file_location('flow_damping_check',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    result=module.damped_offset_test()
    gamma,omega=result['gamma'],result['omega']
    half=math.pi/omega;areas=[0.,0.]
    # A finite quadrature with an independent absolute remainder bound.
    for j in range(512):
        value=quad(lambda t:abs(math.exp(-gamma*t)*math.sin(omega*t)),
                   j*half,(j+1)*half,epsabs=1e-15)[0]
        areas[j%2]+=value
    tail=math.exp(-gamma*512*half)/gamma
    assert abs(areas[0]-result['amplitude_lobes']['positive_area'])<=tail+1e-13
    assert abs(areas[1]-result['amplitude_lobes']['negative_magnitude_area'])<=tail+1e-13
    residue=(areas[0]-areas[1])/sum(areas)
    assert abs(residue-result['amplitude_lobes']['normalized_signed_residue'])<2*tail/sum(areas)


def test_field_equation_selects_amplitude_damping_not_friction_coefficient():
    t,H,mu,lam=s.symbols('t H mu lam',positive=True)
    y=s.Function('y')(t)
    a=s.cosh(H*t)/H;h=H*s.tanh(H*t)
    envelope=s.cosh(H*t)**(-s.Rational(3,2));q=envelope*y
    equation=s.diff(q,t,2)+3*h*s.diff(q,t)+(mu**2+lam/a**2)*q
    expected=s.diff(y,t,2)+(mu**2-s.Rational(9,4)*H**2+
        (lam+s.Rational(3,4))*H**2/s.cosh(H*t)**2)*y
    assert s.simplify(s.expand(equation/envelope-expected))==0
    assert s.limit(-s.diff(envelope,t)/envelope,t,s.oo)==3*H/2


def test_weight_one_klein_gordon_operator_on_the_conformal_cylinder():
    eta,lam,m2,xi=s.symbols('eta lam m2 xi',real=True)
    Q=s.Function('Q')(eta)
    # g=cos(eta)^(-2) gbar; a weight-one physical field is q=cos(eta)*Q.
    q=s.cos(eta)*Q
    original=s.diff(q,eta,2)+2*s.tan(eta)*s.diff(q,eta)+(lam+(m2+12*xi)/s.cos(eta)**2)*q
    cylinder=s.diff(Q,eta,2)+(lam+1+(m2+12*xi-2)/s.cos(eta)**2)*Q
    assert s.simplify(s.expand_trig(original/s.cos(eta)-cylinder))==0
    # The regular comparison is a parameter condition, not automatic mass fading.
    assert s.simplify((m2+12*xi-2).subs({m2:0,xi:s.Rational(1,6)}))==0
    assert s.simplify((m2+12*xi-2).subs({m2:2,xi:0}))==0


def test_null_observer_horizon_and_spacelike_conformal_zero_are_distinct():
    eta,chi=s.symbols('eta chi',real=True)
    H=s.symbols('H',positive=True)
    inverse=s.diag(-H*H,H*H)
    v=eta+chi
    Omega=s.cos(eta)
    dv=s.Matrix([s.diff(v,x) for x in (eta,chi)])
    dOmega=s.Matrix([s.diff(Omega,x) for x in (eta,chi)])
    assert (dv.T*inverse*dv)[0]==0
    norm=s.simplify((dOmega.T*inverse*dOmega)[0])
    assert norm==-H*H*s.sin(eta)**2
    assert norm.subs(eta,s.pi/2)==-H*H
    # Their common endpoint has chi=0; the whole future S3 is not that null sheet.
    assert s.solve((v-s.pi/2).subs(eta,s.pi/2),chi)==[0]


def test_chosen_reciprocal_de_sitter_extension_has_radiation_friedmann_source():
    tau,H,kappa=s.symbols('tau H kappa',positive=True)
    old=1/(H*s.sin(tau));new=s.sin(tau)/H
    assert s.simplify(H*H*old*new)==1
    proper=(1-s.cos(tau))/H
    assert s.diff(proper,tau)==new
    assert proper.subs(tau,0)==0
    hubble=s.diff(new,tau)/(new*new)
    rho=3*(hubble*hubble+1/new**2)/kappa
    assert s.trigsimp(rho*new**4)==3/(kappa*H*H)
    # Conformal scalar curvature R=6(a''/a^3+1/a^2) distinguishes the geometries.
    assert s.trigsimp(6*(s.diff(old,tau,2)/old**3+1/old**2))==12*H*H
    assert s.trigsimp(6*(s.diff(new,tau,2)/new**3+1/new**2))==0


def test_reciprocity_maps_density_exponent_and_equation_of_state():
    tau,H,K=s.symbols('tau H K',real=True,positive=True)
    a=s.Function('a',positive=True)(tau)
    b=1/(H*H*a)
    transformed=s.diff(b,tau)**2+K*b*b
    assert s.simplify(transformed-(s.diff(a,tau)**2+K*a*a)/(H**4*a**4))==0
    # Positive scales written logarithmically avoid assumptions about powers.
    x,h,n,C,kappa=s.symbols('x h n C kappa',real=True)
    new_log=-2*h-x
    transformed_source=kappa*C*s.exp(-n*x-4*h)/3
    new_source=kappa*C*s.exp((2*n-4)*h+n*new_log)/3
    assert s.simplify(transformed_source-new_source)==0
    w=s.symbols('w',real=True)
    new_n=4-3*(1+w)
    new_w=s.simplify(new_n/3-1)
    assert new_w==-s.Rational(2,3)-w
    assert new_w.subs(w,-1)==s.Rational(1,3)
    assert s.simplify(new_w.subs(w,new_w)-w)==0


def test_opposite_complex_windings_have_opposite_charge_and_positive_energy():
    t,A,w=s.symbols('t A w',real=True,positive=True)
    def charge(z):
        return s.simplify(-s.I*(s.conjugate(z)*s.diff(z,t)-z*s.conjugate(s.diff(z,t))))
    def energy(z):
        return s.simplify(s.diff(z,t)*s.conjugate(s.diff(z,t))+w*w*z*s.conjugate(z))
    plus=A*s.exp(s.I*w*t);minus=A*s.exp(-s.I*w*t)
    assert charge(plus)==2*A*A*w
    assert charge(minus)==-2*A*A*w
    assert energy(plus)==energy(minus)==2*A*A*w*w
    standing=plus+minus
    assert charge(standing)==0
    assert s.trigsimp(s.expand_complex(energy(standing)))==4*A*A*w*w
    # A zero signal at a crossing still carries kinetic energy.
    assert s.simplify(standing.subs(t,s.pi/(2*w)))==0
    assert s.simplify(energy(standing).subs(t,s.pi/(2*w)))==4*A*A*w*w


def test_real_phase_offset_does_not_damp_an_oscillator_amplitude():
    t,A,w,delta=s.symbols('t A w delta',real=True)
    z=A*s.exp(s.I*(w*t+delta))
    amplitude_squared=s.simplify(z*s.conjugate(z))
    assert amplitude_squared==A*A
    assert s.diff(amplitude_squared,t)==0
    assert s.diff(amplitude_squared,delta)==0
    q=A*s.cos(w*t+delta)
    E=(s.diff(q,t)**2+w*w*q*q)/2
    assert s.trigsimp(E)==A*A*w*w/2


def test_energy_loss_requires_the_specified_damping_term():
    t,gamma,w0=s.symbols('t gamma w0',real=True,positive=True)
    q=s.Function('q',real=True)(t)
    E=(s.diff(q,t)**2+w0*w0*q*q)/2
    damped_acceleration=-2*gamma*s.diff(q,t)-w0*w0*q
    derivative=s.diff(E,t).subs(s.diff(q,t,2),damped_acceleration)
    assert s.simplify(derivative+2*gamma*s.diff(q,t)**2)==0
    v=s.symbols('v',real=True)
    assert s.ask(s.Q.nonnegative(2*gamma*v*v)) is True


def test_stereographic_damped_spiral_has_unit_sphere_norm_and_constant_bearing():
    r,phi=s.symbols('r phi',real=True,positive=True)
    X=s.Matrix([2*r*s.cos(phi),2*r*s.sin(phi),1-r*r])/(1+r*r)
    assert s.trigsimp(X.dot(X))==1
    grr=s.trigsimp(s.diff(X,r).dot(s.diff(X,r)))
    gpp=s.trigsimp(s.diff(X,phi).dot(s.diff(X,phi)))
    assert grr==4/(1+r*r)**2
    assert gpp==4*r*r/(1+r*r)**2
    assert s.trigsimp(s.diff(X,r).dot(s.diff(X,phi)))==0
    gamma,w=s.symbols('gamma w',positive=True)
    # rdot=-gamma*r, phidot=w. The signed angle is relative to increasing r.
    bearing=s.simplify(s.sqrt(gpp)*w/(s.sqrt(grr)*(-gamma*r)))
    assert bearing==-w/gamma
    u,nu=s.symbols('u nu',positive=True)
    envelope=s.exp(-s.Rational(3,2)*u)
    assert s.simplify(envelope.subs(u,u+2*s.pi/nu)/envelope)==s.exp(-3*s.pi/nu)


def test_s3_mode_count_is_a_family_index_not_a_selected_uv_cutoff():
    ell,L=s.symbols('ell L',integer=True,nonnegative=True)
    count=s.summation((ell+1)**2,(ell,0,L))
    expected=(L+1)*(L+2)*(2*L+3)/6
    assert s.factor(count-expected)==0
    assert s.expand(expected.subs(L,L+1)-expected-(L+2)**2)==0
    assert s.limit(expected,L,s.oo)==s.oo
