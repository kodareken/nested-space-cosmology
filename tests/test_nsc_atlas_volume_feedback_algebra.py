"""Independent action, harmonic and local-ratio checks for volume feedback.

No numerical feedback module is imported. The static elimination is a
separate variational limit; it does not replace the causal volume equation.
"""
import sympy as s

from recursive_horizons.nsc_atlas_s3_polynomials import (
    X,laplacian,reduce_on_sphere,mean,inner_product)


def basis():
    x,y,_,w=X
    return s.Matrix([1,x*x-s.Rational(1,4),y*y+w*w-s.Rational(1,2),x*y*w,
                     y*y*w*w-(y*y+w*w)/8+s.Rational(1,48)])


def gram():
    B=basis()
    return s.Matrix([[inner_product(f,g) for g in B] for f in B])


def test_minimal_action_variation_fixes_the_volume_source_and_kinetic_factor():
    t,r=s.symbols('t r',real=True)
    d,K,m0=s.symbols('d K m0',positive=True)
    a=s.Function('a',positive=True)(t)
    volume=s.Function('s')(t,r);chi=s.Function('chi')(t,r);V=s.Function('V')
    mu2=m0*m0*s.exp(-d*volume)
    density=a**3*s.sin(r)**2*(K*(s.diff(volume,t)**2-s.diff(volume,r)**2/a**2)/2
                              -V(volume)-mu2*chi*chi/2)
    variation=(s.diff(s.diff(density,s.diff(volume,t)),t)
               +s.diff(s.diff(density,s.diff(volume,r)),r)-s.diff(density,volume))
    expected=(K*(s.diff(volume,t,2)+3*s.diff(a,t)*s.diff(volume,t)/a
                 -(s.diff(volume,r,2)+2*s.cot(r)*s.diff(volume,r))/a**2)
              +s.diff(V(volume),volume)-d*mu2*chi*chi/2)
    assert s.simplify(variation/(a**3*s.sin(r)**2)-expected)==0
    kinetic=d*(d+2)/2
    assert kinetic.subs(d,12)==84
    assert s.simplify((d/(2*kinetic)).subs(d,12))==s.Rational(1,14)


def test_stationary_background_expansion_has_a_nonzero_quadratic_source():
    eps,S0,S2,phi,d,m0,Vpp=s.symbols('eps S0 S2 phi d m0 Vpp',real=True)
    delta=eps**2*S2
    # V'(S0)=0; the next derivative is Vpp. The scalar itself is eps*phi.
    force=Vpp*delta-d*m0*m0*s.exp(-d*(S0+delta))*eps**2*phi*phi/2
    quadratic=s.series(force,eps,0,3).removeO().coeff(eps,2)
    assert s.expand(quadratic-(Vpp*S2-d*m0*m0*s.exp(-d*S0)*phi*phi/2))==0
    assert s.simplify(quadratic.subs(S2,0))==-d*m0*m0*phi*phi*s.exp(-S0*d)/2


def test_rescaled_volume_mode_has_the_correct_cylinder_operator_and_source_factor():
    u,lam,M,K,H=s.symbols('u lam M K H',real=True)
    b=s.Function('b')(u);J=s.Function('J')(u)
    mode=s.cosh(u)**(-s.Rational(3,2))*b
    original=s.diff(mode,u,2)+3*s.tanh(u)*s.diff(mode,u)+(M+lam/s.cosh(u)**2)*mode-J/(K*H*H)
    expected=(s.diff(b,u,2)+(M-s.Rational(9,4)+(lam+s.Rational(3,4))/s.cosh(u)**2)*b
              -s.cosh(u)**s.Rational(3,2)*J/(K*H*H))
    assert s.simplify(s.expand_trig(s.cosh(u)**s.Rational(3,2)*original-expected))==0


def test_source_square_decomposes_into_the_complete_five_scalar_basis():
    a,b=s.symbols('a b',real=True)
    phi=a*X[0]+b*X[1]*X[3]
    B=basis();eigenvalues=[0,8,8,15,24]
    coefficients=s.Matrix([a*a/4+b*b/24,a*a,b*b/8,2*a*b,b*b])
    assert s.expand(phi*phi-(B.T*coefficients)[0])==0
    for harmonic,eigenvalue in zip(B,eigenvalues):
        assert reduce_on_sphere(laplacian(harmonic)+eigenvalue*harmonic)==0


def test_harmonic_gram_matrix_is_positive_but_not_diagonal():
    G=gram()
    expected=s.Matrix([[1,0,0,0,0],
                       [0,s.Rational(1,16),-s.Rational(1,24),0,0],
                       [0,-s.Rational(1,24),s.Rational(1,12),0,0],
                       [0,0,0,s.Rational(1,192),0],
                       [0,0,0,0,s.Rational(19,11520)]])
    assert G==expected
    assert G[1,2]!=0
    assert all(G[:n,:n].det().is_positive for n in range(1,6))
    a,b=s.symbols('a b',real=True)
    coefficients=s.Matrix([a*a/4+b*b/24,a*a,b*b/8,2*a*b,b*b])
    assert s.simplify((coefficients.T*G*coefficients)[0]-mean((a*X[0]+b*X[1]*X[3])**4))==0


def test_static_elimination_uses_the_gradient_operator_and_full_gram_pairing():
    Vpp,K,R=s.symbols('Vpp K R',positive=True)
    eigenvalues=[0,8,8,15,24]
    L=s.diag(*[Vpp+K*l/R**2 for l in eigenvalues])
    G=gram();c=s.Matrix(s.symbols('c0:5',real=True));J=s.Matrix(s.symbols('j0:5',real=True))
    assert G*L==L*G
    energy=(c.T*G*L*c)[0]/2-(J.T*G*c)[0]
    solution=L.inv()*J
    stationary={ci:si for ci,si in zip(c,solution)}
    assert all(s.simplify(s.diff(energy,ci).subs(stationary))==0 for ci in c)
    eliminated=s.simplify(energy.subs(stationary))
    assert s.simplify(eliminated+(J.T*G*L.inv()*J)[0]/2)==0
    # A positive sum-of-squares form proves the static correction is nonpositive.
    positive=(J[0]**2/Vpp+(J[1]/4-J[2]/6)**2/L[1,1]
              +J[2]**2/(18*L[1,1])+J[3]**2/(192*L[3,3])
              +19*J[4]**2/(11520*L[4,4]))
    assert s.simplify((J.T*G*L.inv()*J)[0]-positive)==0


def test_local_quartic_coefficient_follows_from_stationary_volume_elimination():
    eps,c,phi,d,mu2,Vpp=s.symbols('eps c phi d mu2 Vpp',positive=True)
    potential=Vpp*(eps**2*c)**2/2+mu2*s.exp(-d*eps**2*c)*(eps*phi)**2/2
    quartic=s.series(potential,eps,0,5).removeO().coeff(eps,4)
    stationary=s.solve(s.diff(quartic,c),c)[0]
    assert stationary==d*mu2*phi*phi/(2*Vpp)
    assert s.simplify(quartic.subs(c,stationary)+d*d*mu2*mu2*phi**4/(8*Vpp))==0


def test_local_gauge_and_kk_mass_ratios_have_fixed_logarithmic_response():
    S,S0,c,eps,d,mu0,lam,Z0=s.symbols('S S0 c eps d mu0 lam Z0',real=True,positive=True)
    mu2=mu0*mu0*s.exp(-d*S)
    delta_mass=lam*s.exp(-(d+2)*S)
    # The inherited KK gauge normalization is Z_F=Z0 exp[(d+2)s].
    gauge=Z0**-s.Rational(1,2)*s.exp(-(d+2)*S/2)
    gauge_log_slope=s.simplify(s.diff(gauge,S)/gauge)
    ratio=delta_mass/mu2
    ratio_log_slope=s.simplify(s.diff(ratio,S)/ratio)
    assert gauge_log_slope.subs(d,12)==-7
    assert ratio_log_slope==-2
    gauge_shift=gauge_log_slope.subs(d,12)*eps**2*c
    ratio_shift=ratio_log_slope*eps**2*c
    assert s.simplify(gauge_shift-s.Rational(7,2)*ratio_shift)==0
    assert s.simplify(s.diff(gauge*gauge/delta_mass,S))==0


def test_constant_background_volume_makes_the_first_response_gauge_invariant():
    t,x,y,z=s.symbols('t x y z',real=True)
    Sstar=s.symbols('Sstar',real=True)
    xi=s.Matrix([s.Function('xi'+str(i))(t,x,y,z) for i in range(4)])
    constant_lie=sum(v*s.diff(Sstar,c) for v,c in zip(xi,(t,x,y,z)))
    assert constant_lie==0
    # For an evolving background this cancellation cannot be assumed.
    S0=s.Function('S0')(t)
    evolving_lie=sum(v*s.diff(S0,c) for v,c in zip(xi,(t,x,y,z)))
    assert evolving_lie==xi[0]*s.diff(S0,t)
