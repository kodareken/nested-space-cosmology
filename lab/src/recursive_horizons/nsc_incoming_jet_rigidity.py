"""Finite algebra for prepared fixed-C0 tangent rigidity of the twenty slots.

The physical P16 bridge is reused, not computed here. Six inverse-gap terms
provide an analytic remainder; exact identities below do not assign its
constants or a finite energy threshold. Fourier amplitudes remain complex.
"""
from functools import lru_cache
import sympy as sp

from .nsc_incoming_cauchy_jets import NORMAL_ENTRIES


def zero(value):return str(sp.factor(sp.cancel(value)))


@lru_cache(maxsize=1)
def inverse_gap_six():
    rho=sp.Symbol('rho');f=sp.Function('f')(rho);lam=sp.Function('lambda')(rho)
    # Telescoping uses independent q/dq symbols and the exact recurrence,
    # avoiding expansion of a sixth-order rational differential expression.
    q=sp.symbols('q0:6');dq=sp.symbols('dq0:6');F,L=sp.symbols('f lambda',nonzero=True)
    residual=sum(dq)-L*sum(q)-F-dq[5]
    rules={q[0]:-F/L,**{dq[j]:L*q[j+1] for j in range(5)}}
    telescoping=zero(residual.subs(rules,simultaneous=True))
    terms=[-f/lam]
    for _ in range(4):terms.append(sp.diff(terms[-1],rho)/lam)
    endpoints=[]
    for k in range(1,5):
        fk,l0=sp.symbols('f_k lambda_0',nonzero=True)
        replacements={sp.diff(f,rho,j):(fk if j==k else sp.Integer(0)) for j in range(k+1)}
        replacements.update({sp.diff(lam,rho,j):sp.Symbol('lambda_'+str(j),nonzero=True) for j in range(k+1)})
        endpoint=terms[k].xreplace(replacements)
        lower=[zero(term.xreplace(replacements)) for term in terms[:k]]
        endpoints.append({'normal_order':k,'endpoint_qk':str(sp.factor(endpoint)),
            'lower_q_endpoint_residuals':lower,'qk_endpoint_residual':zero(endpoint+fk/l0**(k+1))})
    return {'term_count':6,'recurrence':'q0=-f/lambda; q(j+1)=qj_prime/lambda for j=0,...,4',
        'residual_equation':'(sum(q0,...,q5))_prime-lambda*sum(q0,...,q5)-f=q5_prime',
        'telescoping_residual':telescoping,'endpoints':endpoints,
        'uniform_orders':['qj=O(E^(-j-1)), j=0,...,5','q5_prime=O(E^-6)'],
        'endpoint_assumption':'f^(j)(Sigma)=0 exactly for j<k; lambda(Sigma)!=0'}


@lru_cache(maxsize=1)
def normal_chain():
    T=sp.Symbol('T');g=sp.Function('g')(T);v=sp.Function('v')(T)
    a=sp.Symbol('a1',positive=True);gk=sp.Symbol('g_k')
    expression=g;rows=[]
    for k in range(1,5):
        expression=v*sp.diff(expression,T)
        rules={sp.diff(g,T,j):(gk if j==k else sp.Integer(0)) for j in range(k+1)}
        rules.update({sp.diff(v,T,j):(sp.Symbol('v_'+str(j)) if j else -1/a) for j in range(k)})
        endpoint=expression.xreplace(rules)
        rows.append({'normal_order':k,'rho_derivative_at_Sigma':str(sp.factor(endpoint)),
                     'normal_chain_residual':zero(endpoint-(-1/a)**k*gk)})
    return {'fixed_axial_Fourier_coordinate':'z; radial derivative is at fixed z after the owned clock-phase cancellation',
        'operator':'d_rho=(-1/a0)*d_T on the transformed metric Fourier functions',
        'lower_normal_functions_zero':'all g_T^j(z)=0 for j<k, so all their spatial derivatives vanish',
        'orders':rows}


@lru_cache(maxsize=1)
def endpoint_maps(q1_entries):
    """Differentiate only the imported leading Q1 coefficient, not P16 anew."""
    a,r=sp.symbols('a1 r1',positive=True);m,ell=sp.symbols('m ell',real=True,nonzero=True)
    da,dr,A,R=sp.symbols('deltaa deltar A_k R_k')
    context={'a1':a,'r1':r,'m':m,'ell':ell,'deltaa':da,'deltar':dr,'I':sp.I}
    Q1=sp.Matrix([[sp.sympify(entry,locals=context) for entry in row] for row in q1_entries])
    if Q1.shape!=(2,2) or Q1.free_symbols-set(context.values()):raise ValueError('owned two-field leading coefficient required')
    gaps=(-2*sp.I/a**2,2*sp.I/a**2)
    f0=(-gaps[0]*Q1[0,1],-gaps[1]*Q1[1,0])
    rows=[]
    for k in range(1,5):
        substitution={da:(-1/a)**k*A,dr:(-1/a)**k*R}
        calculated=[sp.factor(-forcing.subs(substitution)/gap**(k+1)) for forcing,gap in zip(f0,gaps)]
        expected=[-a**k/(2*sp.I)**(k+1)*((ell/r-sp.I*m)*A-ell*a*R/r**2),
                  (-1)**k*a**k/(2*sp.I)**(k+1)*((ell/r+sp.I*m)*A-ell*a*R/r**2)]
        matrix=sp.Matrix([[sp.diff(value,v) for v in (A,R)] for value in calculated])
        determinant=sp.factor(matrix.det());target=sp.I*m*ell*a**(2*k+1)/(2**(2*k+1)*r**2)
        adjoint_opposite=sp.conjugate(calculated[0].subs({A:sp.conjugate(A),R:sp.conjugate(R)},simultaneous=True))
        solved=sp.solve(calculated,(A,R),dict=True)
        rows.append({'normal_order':k,'scaled_physical_limit':'E^'+str(k+1)+' Q01,Q10',
            'limit01':str(calculated[0]),'limit10':str(calculated[1]),
            'coefficient_identity_residuals':[zero(x-y) for x,y in zip(calculated,expected)],
            'coefficient_map_determinant':str(determinant),'determinant_identity_residual':zero(determinant-target),
            'complex_kernel':[{str(key):str(value) for key,value in row.items()} for row in solved],
            'opposite_transfer_adjoint_residual':zero(calculated[1]-adjoint_opposite),
            'next_endpoint_order':'O(E^-'+str(k+2)+')','physical_bridge_remainder':'O(E^-6)'})
    return {'imported_Q1':[list(row) for row in q1_entries],'complex_Fourier_amplitudes':True,
        'leading_forcing_from_owned_Q1':'f0_entry=-(lambda_entry/E)*Q1_entry',
        'orders':rows,'nondegeneracy':'a1,r1>0 and m*ell!=0'}


def slot_inventory():
    expected=tuple((k,j) for k in range(1,5) for j in range(5-k))
    if NORMAL_ENTRIES!=expected:raise ValueError('owned normal-jet domain changed')
    slots=[{'field':field,'normal_order':k,'spatial_order':j,'necessary_tangent_change':0}
           for field in ('a','r') for k,j in NORMAL_ENTRIES]
    return {'slots':slots,'count':len(slots),'per_field_count':len(NORMAL_ENTRIES),
        'mixed_slot_count_residual':str(len(slots)-20),'per_normal_order_per_field':[4,3,2,1],
        'justification':'each normal function vanishes for all z; differentiate in z only afterward'}


def uniform_induction():
    return {'hypotheses':{
        'family':'first derivative at the fixed reference of a C-infinity real metric history, energy independent, with compact axial support and a fixed compact radial slab away from horizon',
        'upstream':'all variations vanish on an open upstream neighborhood; unchanged original past, source law and physical affine component',
        'intrinsic':'deltaN=deltabeta=deltaa=deltar=0 AS FUNCTIONS of z on Sigma',
        'lapse_shift':'delta partial_T^j N(z)=delta partial_T^j beta(z)=0 for j=1,...,4 AS FUNCTIONS on Sigma',
        'regularity':'bounded forcing derivatives through six, with the corresponding smooth background-frame and gap derivatives, uniformly at each fixed finite omega; a bare twenty-number germ is insufficient',
        'spectrum':'existing unbounded signed-real-energy continuum; Eo=E+omega/2,Ei=E-omega/2 go to positive infinity at each fixed real omega',
        'nondegenerate_block':'reuse authenticated group14 with nonzero m and ell'},
        'reused_physical_bridge':'Paff-P16=O(E^-16) uniformly on compact slab, from positive full-collar bound; no derivative of a numerical error inequality',
        'six_term_residual':'f and required derivatives O(1), cross-band gap lambda=+/-2iE/a0^2+O(E^-1); q5_prime=O(E^-6)',
        'remaining_errors':{'offdiagonal_frame_generator':'O(E^-16)*O(E^-1)=O(E^-17)',
            'physical_source_projector_replacement':'O(E)*O(E^-16)=O(E^-15)',
            'full_coherent_source_minus_vacuum':'O(E exp(-cE)), c>0, including occupations and coherence'},
        'physical_pair_error':'exact unitary left/right transport and zero upstream approximate data give Qphysical-Vo Z Vi^dagger=O(E^-6)',
        'induction':'if lower normal a/r functions vanish, all lower radial metric derivatives vanish; fixed lapse/shift jets give f^(j)(Sigma)=0 for j<k',
        'endpoint':'q0,...,q(k-1) vanish exactly; qk=-f^(k)/lambda^(k+1); higher q terms contribute O(E^-k-2)',
        'frames':'Vo,Vi=I+O(E^-1), so endpoint frame corrections are O(E^-k-2) after lower endpoint terms vanish',
        'injectivity':'full matching makes BOTH complex entries vanish for every fixed omega; nonzero determinant gives A_k=R_k=0; Fourier injectivity yields normal functions zero for k=1,...,4',
        'twenty_slots':'then all partial_z^j delta partial_T^k a,r vanish for k>=1,k+j<=4',
        'higher_derivatives':'needed for the finite remainder estimate but neither required to vanish nor classified by this result',
        'finite_energy_threshold':None,'uniform_large_transfer_bound':None,'sufficiency_for_C0_matching':False,
        'scope':'first-order tangent necessity in the stated functional class; no finite-amplitude exclusion, noncompact-class result, new boundary prescription or global stationarity claim'}


def symbolic_proof(q1_entries):
    proof={'six_inverse_gap':inverse_gap_six(),'normal_chain':normal_chain(),
           'endpoint_maps':endpoint_maps(tuple(tuple(row) for row in q1_entries)),'normal_slot_inventory':slot_inventory()}
    residuals=[]
    def collect(value):
        if isinstance(value,dict):
            for key,item in value.items():
                if key.endswith('_residuals'):residuals.extend(item)
                elif key.endswith('_residual'):residuals.append(item)
                else:collect(item)
        elif isinstance(value,list):
            for item in value:collect(item)
    collect(proof)
    if set(residuals)!={'0'}:raise ArithmeticError('jet-rigidity finite identity failed')
    if any(row['complex_kernel']!=[{'A_k':'0','R_k':'0'}] for row in proof['endpoint_maps']['orders']):
        raise ArithmeticError('both complex Fourier functions must vanish at each order')
    proof['exact_scalar_residual_count']=len(residuals)
    return proof
