"""Exact ingredients of the finite, time-symmetric parent-data certificate.

The tests establish the conformal equation, physical integrands, symmetry
and strict variational-bound algebra. The concentration-compactness argument
is mathematical analysis, not replaced by a numerical-solver success flag.
"""
from itertools import product

import sympy as s


def test_full_hamiltonian_from_conformal_curvature_and_physical_density():
    H,kappa,eps,mu2,f2,P,psi,Lpsi=s.symbols('H kappa eps mu2 f2 P psi Lpsi',positive=True)
    vacuum=3*H**2/kappa
    scalar_curvature=H**2*psi**-5*(-8*Lpsi+6*psi)
    density=vacuum+eps**2*(H**2*psi**-4*P+mu2*f2)/2
    equation=-8*Lpsi+6*psi-6*psi**5-kappa*eps**2*(P*psi+mu2/H**2*f2*psi**5)
    assert s.expand((scalar_curvature-2*kappa*density)*psi**5/H**2-equation)==0


def test_quadratic_and_quartic_conformal_correctors():
    eps,u,v,Du,Dv,kappa,P,m2,f2=s.symbols('eps u v Du Dv kappa P m2 f2')
    psi=1+eps**2*u+eps**4*v
    equation=-8*(eps**2*Du+eps**4*Dv)+6*psi-6*psi**5-kappa*eps**2*(P*psi+m2*f2*psi**5)
    quadratic=-8*Du-24*u-kappa*(P+m2*f2)
    quartic=-8*Dv-24*v-60*u*u-kappa*(P*u+5*m2*f2*u)
    assert s.expand(s.expand(equation).coeff(eps,2)-quadratic)==0
    assert s.expand(s.expand(equation).coeff(eps,4)-quartic)==0


def test_exact_proper_volume_and_excitation_energy_integrands():
    a0,psi,eps,u,v,kappa,mu2,f2,P=s.symbols('a0 psi eps u v kappa mu2 f2 P',positive=True)
    volume_weight=a0**3*psi**6
    density_exc=eps**2*(P/(a0**2*psi**4)+mu2*f2)/2
    expected=eps**2*a0**3*(psi**2*P/a0**2+mu2*psi**6*f2)/2
    assert s.expand(volume_weight*density_exc-expected)==0
    expansion=1+eps**2*u+eps**4*v
    volume_series=s.expand(expansion**6)
    assert volume_series.coeff(eps,2)==6*u
    assert volume_series.coeff(eps,4)==6*v+15*u*u
    energy_series=s.expand(expected.subs(psi,expansion))
    assert s.expand(energy_series.coeff(eps,4)-a0**3*(u*P/a0**2+3*u*mu2*f2))==0


def flip_group():
    generators=[(-1,-1,1,1),(-1,1,1,-1),(1,1,-1,1)]
    return [s.diag(*(s.prod(generators[j][i]**bits[j] for j in range(3)) for i in range(4)))
            for bits in product((0,1),repeat=3)]


def test_existing_symmetry_preserves_sources_has_no_fixed_point_and_minimum_orbit_two():
    X=s.Matrix(s.symbols('x1:5',real=True));a,b=s.symbols('a b',real=True)
    f=a*X[0]+b*X[1]*X[3]
    D=a*a+b*b*(X[1]**2+X[3]**2)-(a*X[0]+2*b*X[1]*X[3])**2
    group=flip_group()
    assert len({tuple(g) for g in group})==8
    assert sum(group,s.zeros(4))==s.zeros(4)
    # The Reynolds projection on degree-one harmonics is zero; in particular
    # no nonzero ambient vector is fixed by the whole group.
    for g in group:
        assert g.T*g==s.eye(4)
        mapping=dict(zip(X,g*X))
        assert s.expand((f*f).subs(mapping,simultaneous=True)-f*f)==0
        assert s.expand(D.subs(mapping,simultaneous=True)-D)==0
    axis=s.Matrix([0,0,1,0])
    assert len({tuple(g*axis) for g in group})==2


def sphere_mean(polynomial,X):
    total=0
    for powers,coefficient in s.Poly(s.expand(polynomial),*X).terms():
        if any(p%2 for p in powers):continue
        half=[p//2 for p in powers];N=sum(half)
        moment=s.prod(s.factorial2(2*p-1) for p in half)/s.prod(4+2*j for j in range(N))
        total+=coefficient*moment
    return s.factor(total)


def test_source_normalization_and_mean_formulas_cancel_pi_exactly():
    z,H,kappa=s.symbols('z H kappa',positive=True)
    V3=2*s.pi**2;vacuum=3*H**2/kappa
    a2=(2/s.pi**2)*vacuum*V3/(H**2*(z+3))
    b2=(12/s.pi**2)*vacuum*V3/(H**2*(z+8))
    assert s.simplify(kappa*a2)==12/(z+3)
    assert s.simplify(kappa*b2)==72/(z+8)
    X=s.Matrix(s.symbols('x1:5',real=True));a,b=s.symbols('a b',real=True)
    f=a*X[0]+b*X[1]*X[3]
    D=a*a+b*b*(X[1]**2+X[3]**2)-(a*X[0]+2*b*X[1]*X[3])**2
    assert s.expand(sphere_mean(f*f,X)-a*a/4-b*b/24)==0
    assert s.expand(sphere_mean(D,X)-3*a*a/4-b*b/3)==0
    p=9/(z+3)+24/(z+8)
    m=3*z/(z+3)+3*z/(z+8)
    assert s.simplify(p+m)==6


def test_global_field_maximum_and_gradient_upper_bound_are_sums_of_squares():
    X=s.Matrix(s.symbols('x1:5',real=True));a,b=s.symbols('a b',positive=True)
    z,x,t,y=X;R2=X.dot(X);f=a*z+b*x*y
    maximum=(a*a+b*b)/(2*b)
    squares=b*((x-y)**2+t*t+(z-a/b)**2)/2
    assert s.expand(maximum-f-squares+b*(R2-1)/2)==0
    D=a*a+b*b*(x*x+y*y)-(a*z+2*b*x*y)**2
    remainder=b*b*(z*z+t*t)+(a*z+2*b*x*y)**2
    assert s.expand(a*a+b*b-D-remainder+b*b*(R2-1))==0
    # The maximum is attained when a/b<1, as ensured by the normalized family.
    mass=s.symbols('mass',positive=True)
    ratio=(mass+8)/(6*(mass+3))
    assert s.simplify(1-ratio-(5*mass+10)/(6*(mass+3)))==0


def test_rational_coercivity_and_bubbling_bounds_for_mass_squared_above_33():
    z=s.symbols('z',positive=True)
    p=9/(z+3)+24/(z+8)
    amplitude_bound=12/(z+3)+72/(z+8)
    assert s.factor(33/z-p-(219*z+792)/(z*(z+3)*(z+8)))==0
    assert s.factor(84/z-amplitude_bound-(612*z+2016)/(z*(z+3)*(z+8)))==0
    assert s.simplify(s.Rational(28,11)-84/z-28*(z-33)/(11*z))==0
    assert 6-s.Rational(28,11)==s.Rational(38,11)
    # Every displayed remainder has strictly positive numerator/denominator
    # for z>33; these are analytic inequalities, not a sampling test.
    ah,bh=12/(z+3),72/(z+8)
    M=s.factor(z*(ah+bh)**2/(4*bh))
    expected=z*(7*z+26)**2/(2*(z+3)**2*(z+8))
    assert s.simplify(M-expected)==0
    difference=s.expand(49*(z+3)**2*(z+8)-z*(7*z+26)**2)
    assert difference==322*z*z+2117*z+3528
    m=6-p
    assert s.factor(M-m-z*(37*z*z+262*z+478)/(2*(z+3)**2*(z+8)))==0
    tau=s.symbols('tau',nonnegative=True)
    ratio=(6+tau*M)/(6+tau*m)
    assert s.simplify(s.diff(ratio,tau)-6*(M-m)/(6+tau*m)**2)==0
    # At tau<=1: ratio<=ratio(1)<(6+49/2)/(6+5)=61/22.
    assert (6+s.Rational(49,2))/(6+5)==s.Rational(61,22)
    # <A>/6<=1 then bounds the CUBE of Q_constant / T by61/88<1.
    assert s.Rational(61,22)/4==s.Rational(61,88)
    assert 1-s.Rational(61,88)==s.Rational(27,88)


def test_the_constant_trial_has_the_stated_normalized_bubbling_ratio():
    V,Amean,Bmean,Bmax=s.symbols('V Amean Bmean Bmax',positive=True)
    constant_trial=V**s.Rational(2,3)*Amean/Bmean**s.Rational(1,3)
    threshold=2**s.Rational(2,3)*6*V**s.Rational(2,3)/Bmax**s.Rational(1,3)
    ratio_cubed=s.simplify((constant_trial/threshold)**3)
    assert s.simplify(ratio_cubed-(Amean/6)**3*Bmax/(4*Bmean))==0


def test_positive_principal_eigenvalue_is_necessary_but_not_the_nonlinear_jacobian():
    eps,kappa,D,m2,f2,psi=s.symbols('eps kappa D m2 f2 psi',positive=True)
    A=6-kappa*eps**2*D;B=6+kappa*eps**2*m2*f2
    # Integrating w*(Lpsi)=psi*(Lw) for the positive principal eigenfunction
    # gives lambda1 int(w psi)=int(w B psi5)>0.
    assert B>0
    jacobian_potential=s.diff(A*psi-B*psi**5,psi)
    assert s.expand(jacobian_potential-(A-5*B*psi**4))==0
    assert jacobian_potential.subs({eps:0,psi:1})==-24


def test_initial_geodesic_normal_raychaudhuri_acceleration_is_positive_potential():
    G,U,Vstar,mu2,eps,f2,kappa=s.symbols('G U Vstar mu2 eps f2 kappa',nonnegative=True)
    rho=G/2+U;spatial_trace=-G/2-3*U
    spacetime_trace=-rho+spatial_trace
    Ric_nn=s.expand(kappa*(rho+spacetime_trace/2))
    assert Ric_nn==-kappa*U
    theta_dot=-Ric_nn  # Kij=0: expansion, shear and normal acceleration vanish.
    potential=Vstar+mu2*eps**2*f2/2
    assert s.expand(theta_dot.subs(U,potential)-kappa*(Vstar+mu2*eps**2*f2/2))==0
