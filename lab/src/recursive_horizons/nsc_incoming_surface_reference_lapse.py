"""Closed B_ref coefficient of r_Tzz from the full order-four Weyl projector."""
from functools import lru_cache
from math import comb

import sympy as sp

I=sp.eye(2)
S1=sp.Matrix([[0,1],[1,0]])
S2=sp.Matrix([[0,-sp.I],[sp.I,0]])
S3=sp.diag(1,-1)


def _matrix(v):return v[0]*S1+v[1]*S2+v[2]*S3
def _zero(v):
    values=list(v) if isinstance(v,sp.MatrixBase) else [v]
    out=[str(sp.factor(sp.simplify(x))) for x in values]
    if set(out)!={'0'}:raise ArithmeticError('nonzero exact reference-lapse residual: '+str(out))
    return out


def _product(A,B):
    return [sum((A[j]*B[n-j] for j in range(n+1)),sp.zeros(2)) for n in range(3)]


def _odd_trace(order):
    A=[sp.Matrix(2,2,sp.symbols(f'A{j}_0:4')) for j in range(order+1)]
    B=[sp.Matrix(2,2,sp.symbols(f'B{j}_0:4')) for j in range(order+1)]
    expression=sum((-1)**j*comb(order,j)*sp.trace(A[j]*B[order-j]+B[j]*A[order-j]) for j in range(order+1))
    return _zero(sp.expand(expression))


def reference_lapse_group(m,L,a,r,Ha,Hr,D,pi):
    """Existing full-line D/(2pi), dk=a dp; only gapped retained groups."""
    M2=m*m+L*L
    return D*L*L*((3*Hr-2*Ha)*M2+8*Hr*m*m)/(120*pi*a*r*M2**3)


@lru_cache(maxsize=1)
def reference_lapse_principal():
    m,L,p,Ha,Hr=sp.symbols('m L p H_a H_r',real=True)
    a,r,D=sp.symbols('a r D',positive=True)
    Q=m*m+L*L+p*p;omega=sp.sqrt(Q)
    h=sp.Matrix([-m,L,p]);v=sp.Matrix([0,-L/r,0]);hd=sp.Matrix([0,-Hr*L,-Ha*p])
    H=_matrix(h);V=_matrix(v);P=(I-H/omega)/2
    # Exact Fourier-Weyl response, expanded only in the coefficient needed
    # here: one formal time frequency and two formal spatial frequencies.
    Pplus=[P,P.diff(p)/(2*a),P.diff(p,2)/(8*a*a)]
    Pminus=[P,-P.diff(p)/(2*a),P.diff(p,2)/(8*a*a)]
    Qplus=[I-Pplus[0],-Pplus[1],-Pplus[2]]
    Qminus=[I-Pminus[0],-Pminus[1],-Pminus[2]]
    A=_product([x*V for x in Pplus],Qminus)
    B=_product([x*V for x in Qplus],Pminus)
    epp=sp.diff(omega,p,2)
    X0=[-(x+y)/(2*omega) for x,y in zip(A,B)]
    X0[2]+=(A[0]+B[0])*epp/(16*a*a*omega**2)
    X1=[(x-y)/(4*omega**2) for x,y in zip(A,B)]
    X1[2]-=(A[0]-B[0])*epp/(16*a*a*omega**3)
    Hplus=[H,S3/(2*a),sp.zeros(2)];Hminus=[H,-S3/(2*a),sp.zeros(2)]
    transport=[];projector=[]
    for left,right,x0 in zip(_product(Hplus,X1),_product(X1,Hminus),X0):
        transport.extend(_zero(left-right-x0))
    for left,right,x1 in zip(_product(Pplus,X1),_product(X1,Pminus),X1):
        projector.extend(_zero(left+right-x1))
    # dT dz² exp(-i nu T+i q z)=i nu q² exp(...), so divide by i.
    deltaP3=-sp.I*X1[2]
    n=h/omega
    delta_p3=-n.diff(p,2).cross(v)/(32*a*a*omega**2)+n.cross(v)*epp/(16*a*a*omega**3)
    third_matrix=_zero(deltaP3-_matrix(delta_p3))

    # Homogeneous lower-order coefficients from the owned recursion.
    # d(r_T) changes h_T by v and h_TT by -4 Hr v while r_TT is fixed.
    eps,Y,Z=sp.symbols('eps h_TTy h_TTz',real=True)
    ht=hd+eps*v;htt=sp.Matrix([0,Y,Z])-4*Hr*eps*v
    p1=h.cross(hd)/(4*omega**3);p1w=h.cross(v)/(4*omega**3)
    p1t=h.cross(htt)/(4*omega**3)-3*h.cross(ht)*(h.dot(ht)/omega)/(4*omega**4)
    homogeneous_p2=h*h.cross(ht).dot(h.cross(ht))/(16*omega**7)-h.cross(p1t)/(2*omega**2)
    delta_p2=homogeneous_p2.diff(eps).subs(eps,0)
    if delta_p2.has(Y,Z):raise ArithmeticError('unowned higher normal background jet survived')
    expected_p2=(h*(2*h.cross(hd).dot(h.cross(v)))/(16*omega**7)
        -h.cross(h.cross(-4*Hr*v))/(8*omega**5)
        +3*((h.dot(v)/omega)*h.cross(h.cross(hd))+(h.dot(hd)/omega)*h.cross(h.cross(v)))/(8*omega**6))
    second_matrix=_zero(delta_p2-expected_p2)

    P1=_matrix(p1);P1zz=_matrix(p1w);P2zz=_matrix(delta_p2)
    # Full trace G4: ordinary P1/P3 pair and both surviving ell=2 pieces.
    ordinary=2*sp.trace(P1*deltaP3)
    moyal02=-sp.trace(P.diff(p,2)*P2zz)/(4*a*a)
    moyal11=-sp.trace(P1.diff(p,2)*P1zz)/(4*a*a)
    kernel=sp.factor(omega*(ordinary+moyal02+moyal11))
    expected=L*L*(-7*Ha*(m*m+L*L)*p*p+Hr*(5*m*m*(m*m+L*L)+(9*L*L+7*m*m)*p*p+2*p**4))/(32*a*a*r*Q**sp.Rational(9,2))
    kernel_residual=_zero(kernel-expected)
    # The third-order traced insertion is not pointwise zero; its full-line
    # integral vanishes by exact odd parity and an integrable tail.
    third_kernel=sp.factor(-omega*sp.trace(P.diff(p,2)*P1zz)/(4*a*a))
    third_expected=-L*m*p/(8*a*a*r*Q**sp.Rational(5,2))
    third_residual=_zero(third_kernel-third_expected)+_zero(third_kernel+third_kernel.subs(p,-p))
    # Reuse beta-function moments, not a new momentum quadrature.
    M=sp.Symbol('M',positive=True)
    moments={degree:sp.simplify((M**(degree-8)*sp.beta(sp.Rational(degree+1,2),sp.Rational(8-degree,2))).rewrite(sp.gamma)) for degree in (0,2,4)}
    wanted_moments={0:sp.Rational(32,35)/M**8,2:sp.Rational(16,105)/M**6,4:sp.Rational(4,35)/M**4}
    moment_residuals=[_zero(moments[j]-wanted_moments[j])[0] for j in (0,2,4)]
    integral=D*L*L/(64*sp.pi*a*r)*(-7*Ha*M*M*moments[2]+Hr*(5*m*m*M*M*moments[0]+(9*M*M-2*m*m)*moments[2]+2*moments[4]))
    closed=reference_lapse_group(m,L,a,r,Ha,Hr,D,sp.pi)
    moment_residuals+=_zero(integral.subs(M,sp.sqrt(m*m+L*L))-closed)
    # The projector diagonal blocks give Tr(H Pn)=omega Tr(Gn).
    g0,gx,gy,gz=sp.symbols('g0 gx gy gz');G=g0*I+gx*S1+gy*S2+gz*S3
    diagonal=-P*G*P+(I-P)*G*(I-P)
    trace_residual=_zero(sp.trace(H*diagonal)-omega*sp.trace(G))
    return {'symbols':(m,L,p,a,r,Ha,Hr,D),'kernel':kernel,'third_order_odd_kernel':third_kernel,
            'closed':closed,'trace_G4_pieces':{'ordinary_P1_deltaP3':sp.factor(ordinary),
                'Moyal2_P0_P2_pair':sp.factor(moyal02),'Moyal2_P1_P1':sp.factor(moyal11)},
            'residuals':{'linear_Weyl_transport':transport,'linear_Weyl_idempotence':projector,
                'third_order_matrix':third_matrix,'homogeneous_P2_variation':second_matrix,
                'full_fourth_order_kernel':kernel_residual,'third_order_and_oddness':third_residual,
                'standard_moments_and_integral':moment_residuals,'projector_energy_trace':trace_residual,
                'odd_Moyal1_pair':_odd_trace(1),'odd_Moyal3_pair':_odd_trace(3)},
            'grade':'r_Tzz grade3; only first-normal background H_a,H_r can enter grade4',
            'measure':'D/(2*pi) times full canonical dk=a*dp; D includes both angular signs',
            'scope':'B_ref only; formal response labels; no physical frequency, quadrature, state or boundary selected'}
