"""Physical high-energy endpoint obstruction: algebra and analytic proof ledger.

No propagation or numerical spectral limit is evaluated. Exact symbolic checks
verify finite identities; the uniform asymptotic conclusion is a separate
analytic argument using the authenticated affine-vacuum defect bound.
"""
from fractions import Fraction as Q
from functools import lru_cache
import math
import sympy as sp


def _zeros(matrix):return [str(sp.factor(value)) for value in matrix]


@lru_cache(maxsize=1)
def scalar_identity():
    E,N,a0,a,ap,m,L=sp.symbols('E N a0 a ap m L',real=True,nonzero=True)
    S,Sprime=sp.symbols('S Sprime');U=L+sp.I*m;Ub=L-sp.I*m
    G=sp.I*N/a0*sp.Matrix([[-E/a,-m-sp.I*L],[-m+sp.I*L,E/a]])
    z=-sp.I*a*S;zprime=-sp.I*(ap*S+a*Sprime)
    riccati=G[1,0]+(G[1,1]-G[0,0])*z-G[0,1]*z*z
    p=a0*a/N;q=a0*ap/N;d=a*a
    F=2*E*S-U+sp.I*(p*Sprime+q*S)+d*Ub*S*S
    residual=sp.factor(zprime-riccati+N/a0*F)
    return {'p':str(p),'q':str(q),'d':str(d),'U':str(U),'generator':[[str(v) for v in row] for row in G.tolist()],
        'scalar_F':str(F),'z':'-i*a*S','ratio_residual_identity':'z_prime-Riccati_G(z)=-(N/a0)*F',
        'residual':str(residual),'antiHermitian_residuals':_zeros(G+sp.conjugate(G.T)),
        'energy_degree_of_G':max(int(sp.Poly(v,E).degree()) for v in G if v!=0),
        'clock_policy':'a0 is the unchanged reference dT/drho clock; N,a,r are varied raw fields'}


@lru_cache(maxsize=1)
def projector_identity():
    z,zb,e,eb,u,ub=sp.symbols('z zbar e ebar u ubar');omega=sp.Symbol('omega',real=True)
    G=sp.Matrix([[-sp.I*omega,u],[-ub,sp.I*omega]])
    P=sp.Matrix([[1,zb],[z,z*zb]])/(1+z*zb)
    rz=-ub+2*sp.I*omega*z-u*z*z
    rzb=-u-2*sp.I*omega*zb-ub*zb*zb
    residual=P.diff(z)*(rz+e)+P.diff(zb)*(rzb+eb)-(G*P-P*G)
    expected=sp.Matrix([[-(zb*e+z*eb),eb-zb*zb*e],
                        [e-z*z*eb,zb*e+z*eb]])/(1+z*zb)**2
    return {'matrix':[[str(v) for v in row] for row in expected.tolist()],
        'residuals':_zeros(residual-expected),
        'interpretation':'zbar/ebar/ubar denote conjugates on the real metric family; no scalar-only surrogate for projector defect',
        'denominator':'(1+abs(z)^2)^2 >= 1 for the real parameter family'}


@lru_cache(maxsize=1)
def recurrence_identity():
    x,p,q,d,U,Ub=sp.symbols('x p q d U Ubar');c=sp.symbols('c1:5');dc=sp.symbols('dc1:5')
    S=sum(c[j]*x**(j+1) for j in range(4));Sprime=sum(dc[j]*x**(j+1) for j in range(4))
    F=sp.expand(S/x-U+sp.I*(p*Sprime+q*S)+d*Ub*S*S)
    errors=[sp.factor(F.coeff(x,0).subs(c[0],U))];rules={'c1':str(U)}
    for n in range(1,4):
        rhs=-sp.I*(p*dc[n-1]+q*c[n-1])-d*Ub*sum(c[j-1]*c[n-j-1] for j in range(1,n))
        rules[str(c[n])]=str(rhs);errors.append(sp.factor(F.coeff(x,n).subs(c[n],rhs)))
    return {'small_parameter':'x=1/(2E)','rules':rules,'orders_0_through_3_residuals':[str(v) for v in errors],
        'remaining_coefficients':{str(n):str(sp.factor(F.coeff(x,n))) for n in range(4,9)},
        'uniformity_reason':'finite coefficient functions and their real-parameter derivatives are bounded on the fixed compact positive-metric slab; remaining powers begin at E^-4'}


@lru_cache(maxsize=1)
def endpoint_identity():
    a,r,m,ell=sp.symbols('a1 r1 m ell',positive=True)
    Amean,Rmean,W0=sp.symbols('Amean Rmean W0',real=True)
    x,eta=sp.symbols('x eta',real=True);c2=sp.Symbol('c2');U=ell/r+sp.I*m
    delta_ap=-Amean/a;delta_rp=-Rmean/a
    delta_Up=-ell*delta_rp/r**2
    delta_c2=sp.factor(-sp.I*(a*a*delta_Up+a*delta_ap*U))
    S=U*x+(c2+eta*delta_c2)*x*x
    P01=sp.I*a*sp.conjugate(S)/(1+a*a*S*sp.conjugate(S))
    coefficient=sp.simplify(sp.limit(sp.diff(P01,eta).subs(eta,0)/x**2,x,0)/4)
    expected=(ell*a*Amean/r-ell*a*a*Rmean/r**2-sp.I*m*a*Amean)/4
    original=sp.simplify(coefficient.subs({Amean:0,Rmean:W0}))
    radius_c2=sp.simplify(delta_c2.subs({Amean:0,Rmean:W0}))
    necessary=sp.solve((sp.re(expected),sp.im(expected)),(Amean,Rmean),dict=True)
    return {'delta_c2_general':str(delta_c2),'general_limit':str(coefficient),
        'general_limit_identity_residual':str(sp.simplify(coefficient-expected)),
        'original_delta_c2':str(radius_c2),'original_limit':str(original),
        'original_limit_identity_residual':str(sp.simplify(original+ell*a*a*W0/(4*r*r))),
        'flat_endpoint_limit_residual':str(sp.simplify(coefficient.subs({Amean:0,Rmean:0}))),
        'two_mean_necessary_conditions':[{str(k):str(v) for k,v in row.items()} for row in necessary],
        'means':'Amean=integral delta(a_T) dz; Rmean=integral delta(r_T) dz, real; intrinsic deltaN=deltaa=deltar=0',
        'original_case':'Amean=0, Rmean=W0>0; compensators vanish on an open neighborhood of Sigma',
        'no_pointwise_or_all_jet_claim':True}


@lru_cache(maxsize=1)
def source_identity():
    f,s,n,Rr,Ri,t=sp.symbols('f s n Rr Ri t',real=True)
    Ch=sp.Matrix([[f,-sp.I*s],[sp.I*s,1-f]])
    D=Ch-sp.diag(0,1)
    remainder=(D*D-f*sp.eye(2)).subs(s*s,f*(1-f))
    sewing=sp.Matrix([[0,1,0],[Rr+sp.I*Ri,0,t]])
    vacuum=sewing*sp.diag(0,1,0)*sp.conjugate(sewing.T)-sp.diag(1,0)
    coisometry=(sewing*sp.conjugate(sewing.T)-sp.eye(2)).subs(t*t,1-Rr*Rr-Ri*Ri)
    return {'horizon_remainder_square_residuals':_zeros(remainder),'sewing_vacuum_residuals':_zeros(vacuum),
        'sewing_coisometry_residuals':_zeros(coisometry),
        'positive_energy_full_source_remainder_bound':'norm(Csrc-diag(0,1,0))=max(sqrt(f_H),n_in) on an open fiber; <= same bound after canonical coisometry',
        'source_decay':'sqrt(f_H)<=exp(-pi*E/kappa); n_in<=exp(-2*pi*E/(Omega*kappa)) for sufficiently large positive E',
        'response_remainder':'norm([K,C_Sigma]-[K,Pv_Sigma])<=2*norm(K)*max(sqrt(f_H),n_in)=O(E*exp(-c*E)), c>0',
        'K_bound':'norm(K)<=integral norm(deltaG) drho=O(E) on the fixed slab',
        'coherence_and_Jost_phase':'retained in the bounded remainder; no phase is set to zero or fitted'}


@lru_cache(maxsize=1)
def linearized_bridge_identity():
    matrices=[sp.Matrix(2,2,sp.symbols(name+'0:4')) for name in ('G','dG','P','PN','X','Q','dR')]
    G,dG,P,PN,X,Q,dR=matrices
    comm=lambda a,b:a*b-b*a
    left=comm(G,X)+comm(dG,P)-comm(G,Q)-comm(dG,PN)-dR
    right=comm(G,X-Q)+comm(dG,P-PN)-dR
    return {'residuals':_zeros(left-right),
        'exact_error_equation':'(X-Q)_rho=[G,X-Q]+[deltaG,P-P4]-deltaR4',
        'initial_error':'X=Q=0 upstream because the physical family and all coefficient jets are unchanged there',
        'unitary_estimate':'norm(X-Q)(1)<=2 integral norm(deltaG)*norm(P-P4) drho+integral norm(deltaR4) drho',
        'background_orders':['Paff-P16=O(E^-16) uniformly from the positive full-collar bound',
            'P16-P4=O(E^-5) on the compact slab from the finite series and bounded coefficients',
            'hence Paff-P4=O(E^-5), stronger than the direct P4 defect bound'],
        'differentiated_defect_order':'deltaR4=O(E^-4) by differentiating its displayed finite rational expression, not an error inequality',
        'variation_size':'deltaG=O(E) for fixed energy-independent raw N,a,r directions',
        'physical_linearized_error_order':'X-deltaP4=O(E^-4) uniformly on the compact slab',
        'endpoint_expansion':'deltaP4_01=L/E²+O(E^-3); full coherent response has the same limit',
        'proof_kind':'analytic uniform Big-O proof; exact algebra checks its finite identities, not numerical constants or a finite-energy threshold'}


@lru_cache(maxsize=1)
def shift_identity():
    beta,k,E,a0=sp.symbols('beta k E a0',real=True,nonzero=True)
    P=sp.Matrix(2,2,sp.symbols('P0:4'))
    vertex=sp.diff(-beta*k*sp.eye(2),beta).subs(k,-E)
    scalar_G=sp.I*E/a0*sp.eye(2)
    return {'vertex':str(vertex),'vertex_residuals':_zeros(vertex-E*sp.eye(2)),
        'commutator_residuals':_zeros(scalar_G*P-P*scalar_G),
        'scope':'flat raw beta_K control at zero energy transfer only; the scalar phase has zero diagonal covariance commutator',
        'nonzero_transfer_claimed':False}


def rational_guard(channel,pulse,config,compensator):
    if (channel['index']!=14 or channel['compact_level']!=1 or channel['angular_level']!=1
            or not math.isfinite(channel['compact_mass']) or not math.isfinite(channel['angular_eigenvalue'])):
        raise ValueError('selected massive angular block14 required')
    if Q(channel['angular_eigenvalue'])<=2 or channel['compact_mass']<=0:
        raise ValueError('nonzero actual mass and angular labels required')
    if Q(pulse['axial_halfwidth'])<=Q(1,20) or Q(pulse['normal_inner'])<=0:
        raise ValueError('original positive compact axial bump and plateau required')
    if pulse['normal_outer']!=.03 or config['magnetic_flux']!=4 or min(config['surface_gravity'],config['omega'])<=0:
        raise ValueError('unchanged pulse and positive frozen source scales required')
    support=compensator['control']['support']
    if support['support_exact_rationals']!=[[203,200],[41,40]] or not compensator['selected_compensation']['same_coefficients_all_cases']:
        raise ValueError('common energy-independent flat compensator binding changed')
    if not all(math.isfinite(v) for v in compensator['selected_compensation']['coefficients_N_a_r']):
        raise ValueError('finite common compensator coefficients required')
    pi_lower=Q(157,50);a_lower=3*pi_lower/2-4
    W_lower=Q(2,3)*Q(1,20);magnitude=Q(2)*Q(16,25)*W_lower/(4*Q(2))
    if not a_lower>Q(16,25) or W_lower!=Q(1,30) or magnitude!=Q(2,375):
        raise ArithmeticError('rational nonzero endpoint margin failed')
    encode=lambda q:{'numerator':q.numerator,'denominator':q.denominator}
    return {'limit_upper_strict':encode(-magnitude),'radius_squared':2,
        'a_squared_lower_from_pi':encode(a_lower),'a_squared_used_strict_lower':encode(Q(16,25)),
        'angular_strict_lower':2,'W0_strict_lower':encode(W_lower),
        'bump_bound':'on |u|<=1/2, exp(1-1/(1-u²))>=exp(-1/3)>2/3; actual halfwidth>1/20',
        'pi_bound':'pi>157/50','intrinsic_geometry':'rho1, a1²=3*pi/2-4, r1²=2',
        'finite_energy_threshold':None,'Big_O_constant':None,
        'actual_common_coefficients_context_only':compensator['selected_compensation']['coefficients_N_a_r']}


def symbolic_proof():
    proof={'scalar':scalar_identity(),'projector':projector_identity(),'recurrence':recurrence_identity(),
           'endpoint':endpoint_identity(),'source':source_identity(),'linearized_bridge':linearized_bridge_identity(),
           'flat_shift_zero_transfer':shift_identity()}
    zeros=[proof['scalar']['residual'],*proof['scalar']['antiHermitian_residuals'],*proof['projector']['residuals'],
           *proof['recurrence']['orders_0_through_3_residuals'],proof['endpoint']['general_limit_identity_residual'],
           proof['endpoint']['original_limit_identity_residual'],proof['endpoint']['flat_endpoint_limit_residual'],
           *proof['source']['horizon_remainder_square_residuals'],*proof['source']['sewing_vacuum_residuals'],
           *proof['source']['sewing_coisometry_residuals'],*proof['linearized_bridge']['residuals'],
           *proof['flat_shift_zero_transfer']['vertex_residuals'],*proof['flat_shift_zero_transfer']['commutator_residuals']]
    if set(zeros)!={'0'}:raise ArithmeticError('endpoint obstruction finite algebra failed')
    proof['checked_scalar_residual_count']=len(zeros)
    return proof
