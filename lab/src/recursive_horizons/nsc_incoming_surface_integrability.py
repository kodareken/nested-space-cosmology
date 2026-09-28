"""Exact constraint contractions on the compatible incoming r-normal family."""
from functools import lru_cache

import sympy as sp

I=sp.eye(2)
S1=sp.Matrix([[0,1],[1,0]])
S2=sp.Matrix([[0,-sp.I],[sp.I,0]])
S3=sp.diag(1,-1)


def zeros(value):
    entries=list(value) if isinstance(value,sp.MatrixBase) else [value]
    result=[str(sp.simplify(v)) for v in entries]
    if set(result)!={'0'}:raise ArithmeticError('nonzero exact incoming contraction: '+str(result))
    return result


def star1(Az,Ak,Bz,Bk):return sp.I*(Az*Bk-Ak*Bz)/2


@lru_cache(maxsize=1)
def first_order_certificate():
    a,r=sp.symbols('a r',positive=True)
    m,ell,k,w=sp.symbols('m ell k w',real=True)
    az,rz,Nz,bz=sp.symbols('a_z r_z N_z beta_z',real=True)
    H=-m*S1+ell/r*S2+k/a*S3
    gap2=m*m+ell*ell/(r*r)+k*k/(a*a)
    P0=(I-H/sp.sqrt(gap2))/2
    Hz=Nz*H+H.diff(a)*az+H.diff(r)*rz-bz*k*I
    Pz=P0.diff(a)*az+P0.diff(r)*rz
    intrinsic={az:0,rz:0,Nz:0,bz:0}
    G1=star1(Pz,P0.diff(k),Pz,P0.diff(k)).subs(intrinsic)
    transport=(star1(Hz,H.diff(k),Pz,P0.diff(k))-star1(Pz,P0.diff(k),Hz,H.diff(k))).subs(intrinsic)
    deltaF1=sp.I*w*P0.diff(r)
    deltaP1=sp.simplify((H*deltaF1-deltaF1*H)/(4*gap2))
    wanted=ell*w/(4*r*r*gap2**sp.Rational(3,2))*(k/a*S1+m*S3)
    M=sp.Matrix(2,2,sp.symbols('M00 M01 M10 M11'))
    vertices={'N':H,'beta':-k*I}
    residuals={'intrinsic_H_z':zeros(Hz.subs(intrinsic)),
               'intrinsic_P0_z':zeros(Pz.subs(intrinsic)),
               'first_Weyl_idempotence':zeros(G1),
               'first_Weyl_transport_commutator':zeros(transport),
               'first_projector_formula':zeros(deltaP1-wanted),
               'first_projector_trace':zeros(sp.trace(deltaP1))}
    for name,vertex in vertices.items():
        residuals[name+'_first_order_insertion']=zeros(sp.trace(deltaP1*vertex))
        # M is an arbitrary Pj_z, including w', w'', U' dependence.
        sym_moyal=sp.I*sp.trace(M*vertex.diff(k)-vertex.diff(k)*M)/4
        residuals[name+'_symmetric_first_Moyal_trace']=zeros(sym_moyal)
        residuals[name+'_vertex_second_k_derivative']=zeros(vertex.diff(k,2))
    if deltaP1==sp.zeros(2):raise ArithmeticError('changed P1 must not be silently set to zero')
    return {'P0':'(I-H/omega)/2; omega^2=m^2+ell^2/r^2+k^2/a^2',
            'delta_P1':'ell*w*(k*sigma1/a+m*sigma3)/(4*r^2*omega^3)',
            'first_projector_generically_nonzero':True,
            'raw_vertices':{'N':'H','beta':'-k*I'},'residuals':residuals,
            'intrinsic_assumption':'all pure spatial metric jets unchanged; N=1,beta=0,a,r constant on T=0',
            'mixed_jets_retained':'w,w_z,w_zz,w_zzz,U,U_z; Pj_z need not vanish',
            'Moyal_scope':'bare vertex trace; nonconstant test-envelope derivatives use the existing compact-test bulk identity'}


def certificate(uv_certificate):
    rows=uv_certificate['result']['recursion_power_certificate']['rows']
    powers=[]
    for order in (2,3,4):
        row=next(r for r in rows if r['formal_order']==order)
        if row['P_power']!=-order-1:raise ValueError('owned generic UV power changed')
        power=row['P_power']+1
        if power>=-1:raise ArithmeticError('constraint tail is not absolutely integrable')
        powers.append({'formal_order':order,'projector_power':row['P_power'],
                       'raw_vertex_power_upper':1,'constraint_insertion_power_upper':power})
    return {'exact':first_order_certificate(),'imported_tail_powers':powers,
            'constraints':['N','beta'],'reference_change_absolutely_integrable':True,
            'zeroth_reference_change':0,'first_order_constraint_insertions':{'N':0,'beta':0},
            'finite_momentum_domain':'32 retained non-LLL groups, each m^2+ell^2/r^2>0',
            'LLL':'owned constant nonzero-k chiral charts; geometric allocation remains separate once',
            'compatible_family':'delta r=T*w(z)+T^3*U(z)/6; delta a=0',
            'scope':{'local_finite_compatible_jets':True,'source_C0_changed':False,
                     'occupation_symmetry_assumed':False,'spatial_domain_or_boundary_selected':False,
                     'pressure_direction_integrability_claimed':False,'full4D_Hadamard_gate_added':False,
                     'physical_surface_solution':False,'reference_quadrature_run':False,
                     'endpoint_terms_removed':False,'rho0_seam_changed':False,
                     'Gamma_rest_assigned':False,'metric_evolution':False}}
