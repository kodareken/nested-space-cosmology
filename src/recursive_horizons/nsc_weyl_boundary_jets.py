"""Boundary jets of the already owned NSC homogeneous Weyl action.

Imports the higher-derivative first-variation identity.  This module applies
it to the existing density and decomposes its stored discrete node gradient;
it does not add an action term or solve a metric field equation.
"""
import numpy as np

FIELDS=('lapse','beta','q_ADM','r')


def weyl_density(jets,coefficient):
    """Existing local-owner L_W, with its eight local jets independent.

    Jet order: N,Ndot,a,adot,addot,r,rdot,rddot.  q_ADM=a.
    """
    N,Nd,a,ad,add,r,rd,rdd=np.asarray(jets).T
    A=(add-ad*Nd/N)/(N*N*a)
    B=ad*rd/(N*N*a*r)
    C=(rdd-rd*Nd/N)/(N*N*r)
    D=(1+(rd/N)**2)/(r*r)
    return -(16*np.pi/3)*coefficient*N*a*r*r*(A-B-C+D)**2


def density_partials(jets,coefficient):
    x=np.asarray(jets,float)
    if x.ndim!=2 or x.shape[1]!=8 or not np.isfinite(x).all() or np.any(x[:,(0,2,5)]<=0):
        raise ValueError('finite independent Weyl jets with positive N,a,r required')
    out=np.empty_like(x);h=1e-28
    for j in range(8):
        z=x.astype(complex);z[:,j]+=1j*h
        out[:,j]=weyl_density(z,coefficient).imag/h
    return out


def nodal_decomposition(history,coefficient):
    """Separate bulk, canonical end jets and the owned stencil's defect.

    G=W E+S P+D^T S Q+(B-S)P+D^T(B-S)Q, with
    B=WD+D^TW, P=L_gdot-D L_gddot. This is the exact chain rule for the
    existing trapezoid/np.gradient discretization. B is not silently set to S.
    """
    t=np.asarray(history.time,float);n=len(t)
    if n<7 or np.any(np.diff(t)<=0):raise ValueError('owned increasing KS node grid required')
    D=np.gradient(np.eye(n),t,axis=0,edge_order=2);D2=D@D
    w=np.empty(n);w[0]=(t[1]-t[0])/2;w[-1]=(t[-1]-t[-2])/2;w[1:-1]=(t[2:]-t[:-2])/2
    W=np.diag(w);S=np.diag(np.r_[-1.,np.zeros(n-2),1.]);B=W@D+D.T@W
    N,a,r=map(np.asarray,(history.lapse,history.a_parallel,history.radius))
    Nd,ad,rd=(D@x for x in (N,a,r));add,rdd=D2@a,D2@r
    jets=np.column_stack((N,Nd,a,ad,add,r,rd,rdd));partial=density_partials(jets,coefficient)
    first=np.zeros((4,n));second=first.copy();value=first.copy()
    for field,j0,j1,j2 in ((0,0,1,None),(2,2,3,4),(3,5,6,7)):
        value[field]=partial[:,j0];first[field]=partial[:,j1]
        if j2 is not None:second[field]=partial[:,j2]
    P=first-second@D.T;Q=second
    bulk=(value-P@D.T)@W
    canonical=P@S+Q@S@D
    defect=P@(B-S).T+Q@(B-S).T@D
    chain=value@W+first@W@D+second@W@D2
    reconstructed=bulk+canonical+defect
    # In proper normal variables delta adot=N*a*delta H_a + ... .
    proper_N=P[0]+ad/N*Q[2]+rd/N*Q[3]
    proper_a=P[2]+ad/a*Q[2];proper_r=P[3]+rd/r*Q[3]
    shear= N*a*Q[2]
    k=-(16*np.pi/3)*coefficient
    L=weyl_density(jets,coefficient)
    # Independent explicit highest-jet coefficients from the same density.
    F=((add-ad*Nd/N)/(N*N*a)-ad*rd/(N*N*a*r)
       -(rdd-rd*Nd/N)/(N*N*r)+(1+(rd/N)**2)/(r*r))
    expected_Qa=2*k*r*r*F/N;expected_Qr=-2*k*a*r*F/N
    return {'time':t,'jets':jets,'density':L,'density_partials':partial,
      'gradient':chain,'bulk_component':bulk,'canonical_boundary_component':canonical,
      'stencil_component':defect,'reconstructed_gradient':reconstructed,
      'P':P,'Q':Q,'proper_boundary_N':proper_N,'proper_boundary_a':proper_a,
      'proper_boundary_r':proper_r,'proper_boundary_shear':shear,
      'quadrature_weights':w,'derivative_matrix':D,'boundary_matrix':B,'ideal_boundary_matrix':S,
      'residuals':{'gradient_decomposition':float(np.max(abs(chain-reconstructed))),
          'highest_jet_Qa':float(np.max(abs(Q[2]-expected_Qa))),
          'highest_jet_Qr':float(np.max(abs(Q[3]-expected_Qr))),
          'Weyl_normal_trace':float(np.max(abs(a*Q[2]+r*Q[3]))),
          'proper_boundary_lapse':float(np.max(abs(proper_N))),
          'stencil_SBP_defect_norm':float(np.linalg.norm(B-S,2))}}
