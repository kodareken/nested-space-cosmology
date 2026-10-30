"""Independent mathematics for the smooth-parent characteristic continuation.

These checks deliberately do not call the numerical handoff implementation.
The second null trace is a restriction of a globally supplied parent state;
its same-time correlation with the horizon trace is not causal preparation.
"""
import sympy as s


def test_radial_klein_gordon_reduction_and_double_null_operator():
    eta,chi,H,m=s.symbols('eta chi H m',real=True,nonzero=True)
    psi=s.Function('Psi')(eta,chi)
    a=1/(H*s.cos(eta))
    field=psi/(a*s.sin(chi))
    kg=(s.diff(field,eta,2)+2*s.tan(eta)*s.diff(field,eta)
        -s.diff(field,chi,2)-2*s.cot(chi)*s.diff(field,chi)
        +(2/s.sin(chi)**2+m*m/s.cos(eta)**2)*field)
    potential=(m*m-2)/s.cos(eta)**2+2/s.sin(chi)**2
    target=s.diff(psi,eta,2)-s.diff(psi,chi,2)+potential*psi
    assert s.simplify(s.expand_trig(a*s.sin(chi)*kg-target))==0

    # Independent chain-rule check with an arbitrary function, not a solution.
    u,v=s.symbols('u v',real=True)
    f=s.Function('f')(u,v)
    transformed=f.subs({u:eta-chi,v:eta+chi})
    wave=s.diff(transformed,eta,2)-s.diff(transformed,chi,2)
    expected=4*s.diff(f,u,v).subs({u:eta-chi,v:eta+chi})
    assert s.simplify(wave-expected)==0


def test_horizon_potential_is_regular_and_has_no_surface_term():
    eta,m=s.symbols('eta m',real=True)
    chi_h=s.pi/2-eta
    potential=(m*m-2)/s.cos(eta)**2+2/s.sin(chi_h)**2
    assert s.trigsimp(potential-m*m/s.cos(eta)**2)==0
    # On the finite chart -pi/2<eta<pi/2 this has no horizon denominator.
    assert s.simplify(potential.subs(eta,0))==m*m


def test_reduced_wave_energy_includes_the_time_dependent_potential_work():
    eta,chi=s.symbols('eta chi',real=True)
    psi=s.Function('Psi')(eta,chi);V=s.Function('V')(eta,chi)
    energy=(s.diff(psi,eta)**2+s.diff(psi,chi)**2+V*psi**2)/2
    flux=-s.diff(psi,eta)*s.diff(psi,chi)
    balance=s.diff(energy,eta)+s.diff(flux,chi)
    reduced=balance.subs(s.diff(psi,eta,2),s.diff(psi,chi,2)-V*psi)
    assert s.simplify(reduced-s.diff(V,eta)*psi**2/2)==0


def test_characteristic_symplectic_flux_is_conserved_for_a_common_operator():
    u,v=s.symbols('u v',real=True)
    A=s.Function('A')(u,v);B=s.Function('B')(u,v);V=s.Function('V')(u,v)
    W_u=A*s.diff(B,u)-B*s.diff(A,u)
    W_v=A*s.diff(B,v)-B*s.diff(A,v)
    divergence=s.diff(W_u,v)+s.diff(W_v,u)
    divergence=divergence.subs({s.diff(A,u,v):-V*A/4,s.diff(B,u,v):-V*B/4})
    assert s.simplify(divergence)==0


def test_physical_killing_flux_divergence_vanishes_after_the_full_radial_pde():
    eta,chi,H,m=s.symbols('eta chi H m',real=True,nonzero=True)
    f,fu,fv=s.symbols('f fu fv',real=True)
    a=1/(H*s.cos(eta));radius_squared=a*a*s.sin(chi)**2
    B=2/s.sin(chi)**2+m*m/s.cos(eta)**2
    Ku=H*s.cos(eta-chi);Kv=H*s.cos(eta+chi)
    # The common angular norm 4*pi/3 is suppressed; it multiplies all terms.
    Tuu=fu*fu;Tvv=fv*fv;Tuv=B*f*f/4
    Ju=radius_squared*(Tuv*Ku+Tvv*Kv)
    Jv=radius_squared*(Tuu*Ku+Tuv*Kv)
    # Mixed derivative from the unrescaled radial KG equation.
    fuv=(-(s.tan(eta)+s.cot(chi))*fu
         -(s.tan(eta)-s.cot(chi))*fv)/2-B*f/4
    # Jet derivatives: Ju has no fu dependence and Jv has no fv dependence.
    Du=(s.diff(Ju,eta)-s.diff(Ju,chi))/2+fu*s.diff(Ju,f)+fuv*s.diff(Ju,fv)
    Dv=(s.diff(Jv,eta)+s.diff(Jv,chi))/2+fv*s.diff(Jv,f)+fuv*s.diff(Jv,fu)
    polynomial=s.Poly(s.expand(Du+Dv),f,fu,fv)
    assert all(s.simplify(s.trigsimp(s.expand_trig(c).rewrite(s.sin)))==0
               for c in polynomial.coeffs())


def test_diamond_algebra_and_constant_positive_potential_dispersion():
    c=s.symbols('c',nonnegative=True)
    p00,p10,p01,p11=s.symbols('p00 p10 p01 p11')
    equation=p11-p10-p01+p00+c*(p11+p10+p01+p00)
    update=((1-c)*(p10+p01)-(1+c)*p00)/(1+c)
    assert s.simplify(equation.subs(p11,update))==0

    omega,k,h=s.symbols('omega k h',real=True)
    tphase=s.exp(-s.I*omega*h/2);xphase=s.exp(s.I*k*h/2)
    plane=equation.subs({p11:tphase**2,p10:tphase/xphase,p01:tphase*xphase,p00:1})
    target=(1+c)*s.cos(omega*h/2)-(1-c)*s.cos(k*h/2)
    assert s.simplify(s.expand_complex(plane/(2*tphase))-target)==0
    ratio=(1-c)/(1+c)
    assert s.simplify(1-ratio**2-4*c/(1+c)**2)==0


def test_constant_kernel_volterra_iterates_have_two_factorials():
    u,v,x,y,M=s.symbols('u v x y M',positive=True)
    n=s.symbols('n',integer=True,nonnegative=True)
    # The positive majorant applies after |V|<=M. Both integration orders count.
    term=(M/4)**n*x**n*y**n/s.factorial(n)**2
    integrated=M*s.integrate(s.integrate(term,(x,0,u)),(y,0,v))/4
    expected=(M*u*v/4)**(n+1)/s.factorial(n+1)**2
    assert s.simplify(integrated/expected)==1


def test_second_characteristic_is_an_exact_same_time_parent_correlation():
    u0,v=s.symbols('u0 v',real=True)
    v0=s.pi/2;uh=u0+v-v0
    eta_second=(u0+v)/2;chi_second=(v-u0)/2
    eta_horizon=(uh+v0)/2;chi_horizon=(v0-uh)/2
    assert s.simplify(eta_second-eta_horizon)==0
    assert s.simplify(chi_second-chi_horizon-(v-v0))==0

    a=s.Function('a');q1=s.Function('q1');q2=s.Function('q2')
    C1,C2=s.symbols('C1 C2',nonzero=True)
    f1=a(eta_horizon)*C1*q1(eta_horizon)*s.sin(chi_horizon)**2
    f2=a(eta_horizon)*C2*q2(eta_horizon)*s.sin(chi_horizon)**2*s.cos(chi_horizon)
    g1=a(eta_second)*C1*q1(eta_second)*s.sin(chi_second)**2
    g2=a(eta_second)*C2*q2(eta_second)*s.sin(chi_second)**2*s.cos(chi_second)
    inherited1=f1*s.sin(chi_second)**2/s.sin(chi_horizon)**2
    inherited2=f2*s.sin(chi_second)**2*s.cos(chi_second)/(s.sin(chi_horizon)**2*s.cos(chi_horizon))
    assert s.simplify(g1-inherited1)==0
    assert s.simplify(g2-inherited2)==0
    # At the corner the two traces agree. Mode-2 quotient requires eta_horizon!=0.
    assert s.simplify((chi_second-chi_horizon).subs(v,v0))==0
    assert s.trigsimp(s.cos(chi_horizon)-s.sin(eta_horizon))==0
