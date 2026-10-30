"""Independent Gaussian-normal and excitation-energy perturbation algebra.

No numerical backreaction producer is imported. Sphere tensors are checked
directly as tangent ambient polynomials; source conservation is kept explicit.
"""
import sympy as s


def symmetric_matrix(prefix):
    a,b,c,d,e,f=s.symbols(f'{prefix}11 {prefix}22 {prefix}33 {prefix}12 {prefix}13 {prefix}23')
    return s.Matrix([[a,d,e],[d,b,f],[e,f,c]])


def test_spatial_metric_equation_from_gaussian_normal_ricci():
    a,h,hd,kappa,mu2,phi=s.symbols('a h hd kappa mu2 phi',nonzero=True)
    k,v,w,R=symmetric_matrix('k'),symmetric_matrix('v'),symmetric_matrix('w'),symmetric_matrix('R')
    I=s.eye(3);g0=a*a*I;g2=a*a*k
    inv0=I/a**2;inv2=-k/a**2
    dot0=2*h*g0;dot2=a*a*(2*h*k+v)
    ddot2=a*a*((2*hd+4*h*h)*k+4*h*v+w)
    trace0=s.trace(inv0*dot0)
    trace2=s.expand(s.trace(inv0*dot2+inv2*dot0))
    assert trace0==6*h and trace2==s.trace(v)
    product2=dot2*inv0*dot0+dot0*inv2*dot0+dot0*inv0*dot2
    # R4_ij=R3_ij+gddot/2+tr(g^-1 gdot)gdot/4-gdot g^-1 gdot/2.
    ricci2=R+ddot2/2+(trace0*dot2+trace2*dot0)/4-product2/2
    gradient=s.Matrix(s.symbols('d1:4'))
    Lambda=hd+3*h*h+2/a**2
    field_equation=2*(ricci2-Lambda*g2-kappa*(gradient*gradient.T+mu2*phi*phi*g0/2))/a**2
    expected=w+3*h*v+h*I*s.trace(v)+2*(R-2*k)/a**2-2*kappa*gradient*gradient.T/a**2-kappa*mu2*phi*phi*I
    assert all(s.expand(z)==0 for z in field_equation-expected)


def test_scalar_curvature_momentum_and_pure_gauge_ricci_on_the_sphere():
    X=s.Matrix(s.symbols('x1:5',real=True));P=s.eye(4)-X*X.T
    sphere=s.Poly(X.dot(X)-1,*X)
    def mod(f):return s.rem(s.Poly(s.expand(f),*X),sphere).as_expr()
    def clean(M):return M.applyfunc(mod)
    def grad(f):return clean(P*s.Matrix([s.diff(f,x) for x in X]))
    def hess(f):
        ambient=s.hessian(f,X);euler=sum(x*s.diff(f,x) for x in X)
        return clean(P*ambient*P-euler*P)
    def lap(f):
        return mod(sum(P[i,j]*s.diff(f,X[i],X[j]) for i in range(4) for j in range(4))
                   -3*sum(x*s.diff(f,x) for x in X))
    def div(T):
        raw=s.Matrix([sum(P[i,a]*s.diff(T[i,j],X[a]) for i in range(4) for a in range(4)) for j in range(4)])
        return clean(P*raw)
    def covector_derivative(v):return clean(P*v.jacobian(X).T*P)
    def rough(T):return clean(P*T.applyfunc(lap)*P+2*T)
    def ricci_variation(T):
        D=covector_derivative(div(T))
        return clean((D+D.T-hess(mod(s.trace(T)))-rough(T))/2+3*T-P*mod(s.trace(T)))

    S=X[0]*X[1];lam=8;HS=hess(S)
    assert lap(S)==-lam*S
    assert clean(ricci_variation(S*P)+HS/2-lam*S*P/2)==s.zeros(4)
    assert clean(ricci_variation(HS)-2*HS)==s.zeros(4)
    A,B,Ad,Bd=s.symbols('A B Ad Bd')
    T=clean(A*S*P+B*HS);trace=mod(s.trace(T))
    curvature=mod(sum(s.diff(div(T)[i],X[j])*P[i,j] for i in range(4) for j in range(4))-lap(trace)-2*trace)
    assert mod(curvature-2*(lam-3)*A*S)==0
    velocity=clean(Ad*S*P+Bd*HS)
    momentum=clean(div(velocity)-grad(mod(s.trace(velocity))))
    assert clean(momentum-2*(Bd-Ad)*grad(S))==s.zeros(4,1)


def scalar_jet():
    a,h,lam,kappa,mu2=s.symbols('a h lam kappa mu2',nonzero=True)
    A,B,Ad,Bd,F,P,D,g,Ggamma,GH,N=s.symbols('A B Ad Bd F P D g Ggamma GH N')
    Add=-6*h*Ad+h*lam*Bd-(lam-4)*A/a**2+kappa*mu2*F+2*kappa*Ggamma/a**2
    Bdd=-3*h*Bd+A/a**2+2*kappa*GH/a**2
    Hres=2*(lam-3)*A/a**2+2*h*(3*Ad-lam*Bd)-kappa*N
    mres=Bd-Ad-kappa*g
    return (a,h,lam,kappa,mu2,A,B,Ad,Bd,F,P,D,g,Ggamma,GH,N,Add,Bdd,Hres,mres)


def test_constraint_propagation_from_metric_and_scalar_source_conservation():
    a,h,lam,kappa,mu2,A,B,Ad,Bd,F,P,D,g,Ggamma,GH,N,Add,Bdd,Hres,mres=scalar_jet()
    # These are the projected KG energy and momentum conservation identities.
    Ndot=-6*h*P-2*h*D/a**2-2*lam*g/a**2
    gdot=-3*h*g+(P-mu2*F)/2+(-Ggamma+(4-lam)*GH)/(2*a**2)
    Hd=(s.diff(Hres,a)*h*a+s.diff(Hres,h)/a**2+s.diff(Hres,A)*Ad
        +s.diff(Hres,Ad)*Add+s.diff(Hres,Bd)*Bdd+s.diff(Hres,N)*Ndot)
    md=Bdd-Add-kappa*gdot
    Hidentity=s.expand(Hd+6*h*Hres+2*lam*mres/a**2).subs(N,P+D/a**2+mu2*F).subs(D,3*Ggamma-lam*GH)
    midentity=s.expand(md-Hres/2+3*h*mres).subs(N,P+D/a**2+mu2*F).subs(D,3*Ggamma-lam*GH)
    assert s.simplify(Hidentity)==0
    assert s.simplify(midentity)==0


def test_trace_evolution_and_its_hamiltonian_residual_term():
    a,h,lam,kappa,mu2,A,B,Ad,Bd,F,P,D,g,Ggamma,GH,N,Add,Bdd,Hres,mres=scalar_jet()
    trace_dot=3*Ad-lam*Bd;trace_ddot=3*Add-lam*Bdd
    identity=trace_ddot+2*h*trace_dot-kappa*(mu2*F-2*P)+2*Hres
    identity=s.expand(identity).subs(N,P+D/a**2+mu2*F).subs(D,3*Ggamma-lam*GH)
    assert s.simplify(identity)==0


def test_excitation_density_pressure_measure_and_strain_work_to_fourth_order():
    eps,a,h,mu2,d,c,phi,phi3,p,p3,rhos,ps,Vstar,m4=s.symbols('eps a h mu2 d c phi phi3 p p3 rhos ps Vstar m4')
    k=symmetric_matrix('k');kd=symmetric_matrix('v')
    f=s.Matrix(s.symbols('f1:4'));g=s.Matrix(s.symbols('g1:4'))
    field=eps*phi+eps**3*phi3;velocity=eps*p+eps**3*p3;gradient=eps*f+eps**3*g
    inverse4=symmetric_matrix('L')
    inverse=(s.eye(3)-eps**2*k+eps**4*inverse4)/a**2
    potential=mu2*(1-d*eps**2*c)*field**2/2
    gradient_squared=(gradient.T*inverse*gradient)[0]
    rho=s.expand(velocity**2/2+gradient_squared/2+potential+eps**4*rhos)
    pressure=s.expand(velocity**2/2-gradient_squared/6-potential+eps**4*ps)
    rho2=(p*p+f.dot(f)/a**2+mu2*phi*phi)/2
    p2=(p*p-f.dot(f)/(3*a**2)-mu2*phi*phi)/2
    rho4=p*p3+f.dot(g)/a**2-(f.T*k*f)[0]/(2*a**2)+mu2*phi*phi3-d*mu2*c*phi*phi/2+rhos
    p4=p*p3-f.dot(g)/(3*a**2)+(f.T*k*f)[0]/(6*a**2)-mu2*phi*phi3+d*mu2*c*phi*phi/2+ps
    assert s.expand(rho.coeff(eps,4)-rho4)==0
    assert s.expand(pressure.coeff(eps,4)-p4)==0
    assert s.diff(rho.coeff(eps,4),inverse4[0,0])==0
    measure=1+eps**2*s.trace(k)/2+eps**4*m4
    # Subtraction is the constant vacuum stress on the actual metric, not
    # merely the background vacuum energy. The unknown measure m4 cancels.
    excitation=s.expand(measure*(Vstar+rho)-Vstar*measure).coeff(eps,4)
    assert s.expand(excitation-rho4-s.trace(k)*rho2/2)==0
    assert s.diff(excitation,m4)==0
    background_only=s.expand(measure*(Vstar+rho)-Vstar).coeff(eps,4)
    assert s.diff(background_only,m4)==Vstar
    S2=f*f.T/a**2+s.eye(3)*(p*p-f.dot(f)/a**2-mu2*phi*phi)/2
    strain=s.trace(kd*S2)/2
    expected=(f.T*kd*f)[0]/(2*a**2)+s.trace(kd)*(p*p-f.dot(f)/a**2-mu2*phi*phi)/4
    assert s.expand(strain-expected)==0
    W4=-3*h*(p4+s.trace(k)*p2/2)-strain
    assert s.expand(W4+3*h*(p4+s.trace(k)*p2/2)+expected)==0


def test_physical_inverse_and_volume_measure_have_the_expected_orders():
    t=s.symbols('t');k=symmetric_matrix('k');l=symmetric_matrix('l');I=s.eye(3)
    metric=I+t*k+t*t*l
    inverse=I-t*k+t*t*(k*k-l)
    product=(metric*inverse-I).applyfunc(s.expand)
    assert all(z.coeff(t,1)==0 and z.coeff(t,2)==0 for z in product)
    determinant=s.expand(metric.det())
    first=determinant.coeff(t,1);second=determinant.coeff(t,2)
    volume_second=s.expand(second/2-first*first/8)
    expected=s.trace(l)/2+s.trace(k)**2/8-s.trace(k*k)/4
    assert first==s.trace(k)
    assert s.expand(volume_second-expected)==0


def test_initial_conformal_excitation_energy_expansion():
    eps,a,u,w,mu2,phi,gradient_squared=s.symbols('eps a u w mu2 phi gradient_squared')
    Psi=1+eps**2*u+eps**4*w
    # Zero initial velocities, h=a^2 Psi^4 gamma, and exact vacuum subtraction.
    excitation=a**3*eps**2*(Psi**2*gradient_squared/a**2+Psi**6*mu2*phi*phi)/2
    coefficient=s.expand(excitation).coeff(eps,4)
    expected=a**3*(u*gradient_squared/a**2+3*u*mu2*phi*phi)
    assert s.expand(coefficient-expected)==0
    assert s.diff(coefficient,w)==0
