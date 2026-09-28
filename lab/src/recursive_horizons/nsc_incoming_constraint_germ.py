"""Structural two-direction relation for the retained incoming constraints.

This owns a conditional local-germ relation, not selected initial data,
a source-accuracy certificate, a Cauchy-surface solution or metric evolution.
"""
import sympy as sp


def structural_constraint_plane():
    """Grade and reflection of the EXISTING four-derivative bulk owners."""
    u,v=sp.symbols('u v',real=True)
    allowed=[(i,j) for i in range(2) for j in range(3) if 3*i+2*j<=4]
    even=[pair for pair in allowed if pair[1]%2==0]
    odd=[pair for pair in allowed if pair[1]%2==1]
    if set(even)!={(0,0),(1,0),(0,2)} or odd!=[(0,1)]:
        raise ArithmeticError('retained derivative grade produced unexpected monomials')
    N,beta,a,r,m,ell,k=sp.symbols('N beta a r m ell k',real=True,nonzero=True)
    s1=sp.Matrix([[0,1],[1,0]]);s2=sp.Matrix([[0,-sp.I],[sp.I,0]]);s3=sp.diag(1,-1);I=sp.eye(2)
    H=-N*m*s1+N*ell/r*s2+N*k/a*s3-beta*k*I
    reflected=H.subs({ell:-ell,k:-k,beta:-beta},simultaneous=True)
    residual=sp.simplify(reflected-s1*H*s1)
    vertexN=H.diff(N);vertexBeta=H.diff(beta)
    reflectedN=vertexN.subs({ell:-ell,k:-k,beta:-beta},simultaneous=True)
    reflectedBeta=vertexBeta.subs({ell:-ell,k:-k,beta:-beta},simultaneous=True)
    vertex_residuals=[sp.simplify(reflectedN-s1*vertexN*s1),
                      sp.simplify(reflectedBeta+s1*vertexBeta*s1)]
    if residual!=sp.zeros(2) or any(x!=sp.zeros(2) for x in vertex_residuals):
        raise ArithmeticError('Hamiltonian/vertex reflection identity failed')
    SN,Sb,cvv=sp.symbols('S_N S_beta c_v2',real=True)
    cu,cv=sp.symbols('c_u c_v',real=True,nonzero=True)
    EN=SN+cu*u+cvv*v*v;Eb=Sb+cv*v
    formal_v=-Sb/cv;formal_u=-(SN+cvv*formal_v**2)/cu
    substitution={u:formal_u,v:formal_v}
    identities=[sp.simplify(EN.subs(substitution)),sp.simplify(Eb.subs(substitution))]
    if identities!=[0,0]:raise ArithmeticError('conditional algebraic relation failed')
    return {'variables':{'u':'delta r_TTT','v':'delta r_Tz'},'derivative_grades':{'u':3,'v':2},
            'retained_maximum_grade':4,'allowed_lapse_monomials':[list(x) for x in even],
            'allowed_shift_response_monomials':[list(x) for x in odd],
            'constraint_form':{'N':'S_N+c_u*u+c_v2*v^2','beta':'S_beta+c_v*v'},
            'reflection':{'coordinates':'z -> -z; k -> -k; ell -> -ell; beta -> -beta',
                          'matrix':'sigma1','Hamiltonian_residual':'0','vertex_residuals':['0','0'],
                          'state_reflection_required':False,'state_enters_fixed_constants':True},
            'conditional_relation':{'v':'-S_beta/c_v','u':'-(S_N+c_v2*v^2)/c_u',
                                    'requires':['c_u != 0','c_v != 0','finite retained source and reference insertion'],
                                    'substitution_residuals':[str(x) for x in identities]},
            'coefficients_are_action_responses_not_adjustable_couplings':True,
            'physical_initial_data_or_surface_solution':False}
