"""Physical fixed-transfer first-normal necessity: finite algebra/proof ledger.

This module does not evolve fields or scan energies. Its SU(2) frames belong
to mathematical Riccati approximants, not normalized archived physical modes.
Analytic uniform estimates are stated separately from exact symbolic checks.
"""
from functools import lru_cache
import sympy as sp


def zeros(matrix):return [str(sp.simplify(sp.factor(value))) for value in matrix]
def entries(matrix):return [[str(sp.factor(value)) for value in row] for row in matrix.tolist()]
def comm(a,b):return a*b-b*a


@lru_cache(maxsize=1)
def frame_projector_identity():
    x,y,xp,yp,g1,g2,g3=sp.symbols('x y xp yp g1 g2 g3',real=True)
    z=x+sp.I*y;zb=x-sp.I*y;den=1+x*x+y*y
    V=sp.Matrix([[1,-zb],[z,1]])/sp.sqrt(den);Vdag=sp.conjugate(V.T)
    P0=sp.diag(1,0);P=sp.Matrix([[1,zb],[z,z*zb]])/den
    S1=sp.Matrix([[0,1],[1,0]]);S2=sp.Matrix([[0,-sp.I],[sp.I,0]]);S3=sp.diag(1,-1)
    G=sp.I*(g1*S1+g2*S2+g3*S3)
    Vp=V.diff(x)*xp+V.diff(y)*yp;Pp=P.diff(x)*xp+P.diff(y)*yp
    L=Vdag*G*V-Vdag*Vp;R=Pp-comm(G,P)
    ar,ai=sp.symbols('alpha_re alpha_im',real=True);alpha=ar+sp.I*ai
    off=sp.Matrix([[0,alpha],[-sp.conjugate(alpha),0]])
    defect=-comm(off,P0);norm_square=(ar*ar+ai*ai)*sp.eye(2)
    return {'frame':entries(V),'unitarity_residuals':zeros(Vdag*V-sp.eye(2)),
        'determinant_residual':str(sp.factor(V.det()-1)),
        'projector_residuals':zeros(V*P0*Vdag-P),
        'transformed_defect_residuals':zeros(Vdag*R*V+comm(L,P0)),
        'transformed_generator_antiHermitian_residuals':zeros(L+sp.conjugate(L.T)),
        'offdiagonal_norm_square_residuals':zeros(sp.conjugate(off.T)*off-norm_square),
        'defect_norm_square_residuals':zeros(sp.conjugate(defect.T)*defect-norm_square),
        'identity':'V^dagger R16 V=-[L,P0]; norm(L_off)=norm(R16)',
        'frame_scope':'normalized mathematical P16 frame; archived columns and physical state are unchanged'}


@lru_cache(maxsize=1)
def inverse_gap_identity():
    rho=sp.Symbol('rho');f=sp.Function('f')(rho);lam=sp.Function('lambda')(rho)
    terms=[-f/lam]
    for _ in range(3):terms.append(sp.diff(terms[-1],rho)/lam)
    total=sum(terms);residual=sp.cancel(sp.diff(total,rho)-lam*total-f-sp.diff(terms[-1],rho))
    f0,fp,bp=sp.symbols('f0 fprime lambdaprime')
    b=sp.Symbol('lambda0',nonzero=True)
    q1=-fp/b**2+f0*bp/b**3
    return {'term_count':4,'terms':['q0=-f/lambda','q1=q0_prime/lambda','q2=q1_prime/lambda','q3=q2_prime/lambda'],
        'residual_equation':'(q0+q1+q2+q3)_prime-lambda*(q0+q1+q2+q3)-f=q3_prime',
        'telescoping_residual':str(residual),
        'endpoint_q1_when_f_zero':str(q1.subs(f0,0)),
        'endpoint_identity_residual':str(sp.factor(q1.subs(f0,0)+fp/b**2)),
        'uniform_orders':['q0=O(E^-1)','q1=O(E^-2)','q2=O(E^-3)','q3=O(E^-4)','residual=q3_prime=O(E^-4)'],
        'required_uniformity':'fixed transfer; smooth bounded coefficient derivatives through q3_prime on compact slab; |lambda|>=c*E for sufficiently large E'}


@lru_cache(maxsize=1)
def pair_transport_identity():
    P0=sp.diag(1,0)
    D=sp.Matrix(2,2,sp.symbols('D0:4'));Po=sp.Matrix(2,2,sp.symbols('Po0:4'));Pi=sp.Matrix(2,2,sp.symbols('Pi0:4'))
    source=D*Pi-Po*D
    compatibility=Po*source+source*Pi-source-(D*(Pi*Pi-Pi)-(Po*Po-Po)*D)
    z01,z10,lo01,lo10,li01,li10=sp.symbols('z01 z10 lo01 lo10 li01 li10')
    Z=sp.Matrix([[0,z01],[z10,0]]);Lo=sp.Matrix([[0,lo01],[lo10,0]]);Li=sp.Matrix([[0,li01],[li10,0]])
    omitted=-Lo*Z+Z*Li;F=D*P0-P0*D
    Co=sp.Matrix(2,2,sp.symbols('Co0:4'));Ci=sp.Matrix(2,2,sp.symbols('Ci0:4'))
    Ao=sp.Matrix(2,2,sp.symbols('Ao0:4'));Ai=sp.Matrix(2,2,sp.symbols('Ai0:4'));Zp=sp.Matrix(2,2,sp.symbols('Zp0:4'))
    transformed=(Co*Z+Zp-Z*Ci)-(Ao+Co)*Z+Z*(Ai+Ci)-F
    expected=Zp-Ao*Z+Z*Ai-F
    return {'source_projector_compatibility_residuals':zeros(compatibility),
        'offdiagonal_source_diagonal_residuals':[str(F[0,0]),str(F[1,1])],
        'approximate_projector_pair_constraint_residuals':zeros(P0*Z+Z*P0-Z),
        'ignored_generator_coupling_offdiagonal_residuals':[str(omitted[0,1]),str(omitted[1,0])],
        'two_frame_chain_rule_residuals':zeros(transformed-expected),
        'exact_equation':'Q_prime=G_o Q-Q G_i+D P_i-P_o D; Q(upstream)=0',
        'exact_projector_pair_constraint':'P_o Q+Q P_i=Q',
        'source_error_bound':'norm(D*(P_i-P16_i)-(P_o-P16_o)*D)<=norm(D)*(norm(P_i-P16_i)+norm(P_o-P16_o))',
        'left_right_bound':'norm(error)(rho)<=integral norm(residual) drho; left and right exact background transport are unitary'}


@lru_cache(maxsize=1)
def endpoint_pair_identity():
    a,r,m,ell=sp.symbols('a1 r1 m ell',positive=True)
    E,omega=sp.symbols('E omega',real=True)
    n,b,u,v=sp.symbols('deltaN deltabeta deltaa deltar')  # complex Fourier directions
    Ah,Rh=sp.symbols('Ahat Rhat')  # intentionally NOT real
    S1=sp.Matrix([[0,1],[1,0]]);S2=sp.Matrix([[0,-sp.I],[sp.I,0]]);S3=sp.diag(1,-1);I=sp.eye(2)
    A=-sp.I/a**2*S3;B=sp.I/a*(-m*S1+ell/r*S2);P0=sp.diag(1,0)
    P1=a*sp.Matrix([[0,m+sp.I*ell/r],[m-sp.I*ell/r,0]])/2
    D1=sp.I/a*((-n/a+u/a**2)*S3+b*I)
    D0=sp.I/a*(-m*n*S1+ell*(n/r-v/r**2)*S2)
    forcing0=comm(D0,P0)+comm(D1,P1)
    Q1=sp.Matrix([[0,((m+sp.I*ell/r)*u-sp.I*a*ell*v/r**2)/2],
                  [((m-sp.I*ell/r)*u+sp.I*a*ell*v/r**2)/2,0]])
    intrinsic={n:0,b:0,u:0,v:0}
    Q1p=Q1.diff(u)*(-Ah/a)+Q1.diff(v)*(-Rh/a)
    common=ell*a*Ah/r-ell*a*a*Rh/r**2
    L01=(common-sp.I*m*a*Ah)/4;L10=(common+sp.I*m*a*Ah)/4
    Q2=sp.Matrix([[0,L01],[L10,0]])
    P2=sp.Matrix(2,2,sp.symbols('P2_0:4'))
    forcing1=comm(D0,P1)+comm(D1,P2)+omega*(D1*P1+P1*D1)/2
    pair_solution=sp.solve((L01,L10),(Ah,Rh),dict=True)
    coefficient_matrix=sp.Matrix([[sp.diff(L01,Ah),sp.diff(L01,Rh)],[sp.diff(L10,Ah),sp.diff(L10,Rh)]])
    conjugate_opposite=sp.conjugate(L01.subs({Ah:sp.conjugate(Ah),Rh:sp.conjugate(Rh)},simultaneous=True))
    Eo,Ei=E+omega/2,E-omega/2
    return {'Weyl_midpoint_residual':str(sp.simplify(Ei+(Eo-Ei)/2-E)),
        'background_A':entries(A),'background_B':entries(B),'D1':entries(D1),'D0':entries(D0),
        'Q1':entries(Q1),'leading_transport_residuals':zeros(comm(A,Q1)+forcing0),
        'leading_projector_constraint_residuals':zeros(P0*Q1+Q1*P0-Q1),
        'Q1_intrinsic_zero_residuals':zeros(Q1.subs(intrinsic)),
        'omega_anticommutator_residuals':zeros(A*Q1+Q1*A),
        'next_source_intrinsic_zero_residuals':zeros(forcing1.subs(intrinsic)),
        'endpoint_transport_residuals':zeros(comm(A,Q2)-Q1p),
        'endpoint_diagonal_projector_residuals':zeros(P0*Q2+Q2*P0-Q2),
        'lapse_and_shift_Q1_residuals':zeros(Q1.diff(n))+zeros(Q1.diff(b)),
        'limit01':str(sp.factor(L01)),'limit10':str(sp.factor(L10)),
        'complex_Fourier_solution':[{str(k):str(val) for k,val in row.items()} for row in pair_solution],
        'coefficient_map_determinant':str(sp.factor(coefficient_matrix.det())),
        'opposite_transfer_Hermiticity_residual':str(sp.simplify(L10-conjugate_opposite)),
        'complex_amplitudes_retained':True,'intrinsic_requirement':'deltaN,deltabeta,deltaa,deltar vanish as functions of z on Sigma',
        'no_omega_lapse_or_shift_first_normal_term_at_E_minus2':True,
        'normal_derivative_conversion':'deltaa_rho(1)=-Ahat/a1; deltar_rho(1)=-Rhat/a1; all undifferentiated intrinsic directions zero'}


def uniform_proof():
    return {'fixed_transfer':'Eo=E+omega/2, Ei=E-omega/2; omega fixed finite real, E→+infinity',
        'background':'Paff(E±omega/2)-P16(E±omega/2)=O(E^-16) uniformly on the compact slab by the positive full-collar bound',
        'mathematical_frames':'V16 unitary, P16=V16 P0 V16^dagger; V16=I+O(E^-1), needed positive-order rho derivatives=O(E^-1), V_o-V_i and its needed rho derivatives=O(E^-2) for fixed omega',
        'transformed_generators':'L=V^dagger G V-V^dagger V_prime; offdiagonal norm equals projector-defect norm=O(E^-16)',
        'cross_band_gaps':'lambda01=(Lo)00-(Li)11=-2iE/a0²+O(E^-1); lambda10=+2iE/a0²+O(E^-1), uniformly separated from0 eventually',
        'transformed_forcing':'F=Dtilde P0-P0 Dtilde is offdiagonal; D=E D1+D0 with D1 diagonal and D0=O(1), hence F and required rho derivatives O(1)',
        'four_terms':'q0=-f/lambda; q(n+1)=q(n)_prime/lambda for n=0,1,2; residual=q3_prime=O(E^-4)',
        'frame_coupling_error':'ignored offdiagonal L times q=O(E^-16)*O(E^-1)=O(E^-17)',
        'physical_forcing_error':'D times physical-minus-P16 projectors=O(E)*O(E^-16)=O(E^-15)',
        'upstream':'all metric directions vanish in an open upstream neighborhood, so F and q0..q3 vanish there',
        'physical_pair_error':'unitary left/right propagation bounds Qexact-Vo*q*Vi^dagger=O(E^-4)',
        'endpoint':'intrinsic delta g is identically0 in z on Sigma, so F(Sigma)=0; Qapprox starts at E^-2 with the checked pair coefficients',
        'coherent_source_error':'exact full-minus-vacuum source in the pair equation is D(Ci-Pi)-(Co-Po)D; source occupations and coherence give O(E*exp(-cE)), c>0',
        'uniformity_hypotheses':['fixed compact slab away from horizon, positive background a,r',
            'smooth real energy-independent retarded metric directions, with sufficient bounded rho derivatives through q3_prime',
            'smooth compact axial first-normal tangents, and unchanged original upstream/past source law',
            'all intrinsic metric directions, including shift, vanish as functions of z on Sigma'],
        'Fourier_injectivity':'if matching holds at all pairs, both checked entries vanish for every fixed omega; m*ell!=0 implies Ahat=Rhat=0, hence deltaa_T=deltar_T=0 pointwise for smooth compact axial tangents',
        'analytic_not_numerical':'uniform orders and injectivity are analytic deductions; finite algebra checks do not supply constants or an energy threshold',
        'finite_energy_threshold':None,'uniform_large_transfer_estimate':None}


def symbolic_proof():
    proof={'frame_projector':frame_projector_identity(),'inverse_gap':inverse_gap_identity(),
        'pair_transport':pair_transport_identity(),'endpoint_pair':endpoint_pair_identity()}
    residuals=[]
    for row in proof.values():
        for name,value in row.items():
            if name.endswith('_residuals'):residuals.extend(value)
            elif name.endswith('_residual'):residuals.append(value)
    if set(residuals)!={'0'}:raise ArithmeticError('fixed-transfer finite identity failed')
    proof['exact_scalar_residual_count']=len(residuals)
    return proof
