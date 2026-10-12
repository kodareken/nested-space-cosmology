"""Closed reference coefficient of (r_Tz)^2, retaining the full Weyl product."""
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import sympy as sp

from .nsc_incoming_surface_reference_lapse import I,S1,S2,S3,_matrix,_zero


def reference_quadratic_group(m,L,a,r,D,pi):
    """Coefficient of v², not its second derivative; full-line group measure."""
    M2=m*m+L*L
    return D*L*L*(L*L+3*m*m)/(60*pi*a*r*r*M2**3)


@lru_cache(maxsize=1)
def reference_quadratic_lapse():
    m,L,p,tau,zeta,v=sp.symbols('m L p r_T r_z v',real=True)
    a,r,D=sp.symbols('a r D',positive=True)
    W=sp.sqrt(m*m+L*L+p*p);h=sp.Matrix([-m,L,p]);V=sp.Matrix([0,-L/r,0])
    H=_matrix(h);P0=(I-H/W)/2;Q0=I-P0;Hk=S3/a
    # L=ell/r at fixed ell; k=a*p and a is fixed in this grade-four extraction.
    Dr=lambda expression:expression.diff(r)-L/r*expression.diff(L)
    P0r=Dr(P0);P0k=P0.diff(p)/a
    P0z=P0r*zeta;Hz=_matrix(V)*zeta
    G1=sp.I*(P0z*P0k-P0k*P0z)/2
    F1=sp.I*P0r*tau-sp.I*(Hz*P0k+P0k*Hz-Hk*P0z-P0z*Hk)/2
    rawP1=-P0*G1*P0+Q0*G1*Q0+(H*F1-F1*H)/(4*W*W)
    pv=h.cross(V)/(4*W**3);M=_matrix(pv)
    g=-m*L/(4*a*r*W**3)
    P1=tau*M+zeta*g*I
    residuals={'full_first_Weyl_transport_term':_zero(F1-sp.I*P0r*tau),
               'full_first_projector':_zero(rawP1-P1),
               'first_idempotence':_zero(P0*P1+P1*P0-P1+G1),
               'first_transport':_zero(H*P1-P1*H-F1)}

    # On T=0 the radius is spatially constant, but r_T=tau(z), r_Tz=v.
    # P1's spatial-Weyl scalar term has a nonzero time derivative g*v.
    P1surface=tau*M;P1T=tau*tau*Dr(M)+g*v*I;P1z=v*M
    G2=P1surface*P1surface+sp.I*(P1z*P0k-P0k*P1z)/2
    F2=sp.I*P1T+sp.I*(Hk*P1z+P1z*Hk)/2
    P2surface=-P0*G2*P0+Q0*G2*Q0+(H*F2-F2*H)/(4*W*W)
    P2surface=P2surface.applyfunc(sp.factor)
    scalar2=-L*L*p/(8*a*r*W**5)
    hom2=tau*tau*(pv.dot(pv)*h/W-h.cross(Dr(pv))/(2*W*W))
    residuals.update(second_idempotence=_zero(P0*P2surface+P2surface*P0-P2surface+G2),
                     second_transport=_zero(H*P2surface-P2surface*H-F2),
                     second_surface_decomposition=_zero(P2surface-_matrix(hom2)-v*scalar2*I),
                     second_point_transport=_zero(F2.subs(tau,0)),
                     second_point_scalar=_zero(P2surface.subs(tau,0)-v*scalar2*I))
    # tau=v*z: its second spatial derivative is zero, but the quadratic
    # homogeneous P2 dependence supplies P2_zz=v²*d_tau²(P2).
    P2zz=v*v*P2surface.diff(tau,2)
    P1zk=v*M.diff(p)/a
    P2point=v*scalar2*I
    ordinary=sp.trace(P2point*P2point)
    moyal02=-sp.trace(P0.diff(p,2)*P2zz)/(4*a*a)
    moyal11=sp.trace(P1zk*P1zk)/4
    energy=W*(ordinary+moyal02+moyal11)
    kernel=sp.factor(energy.diff(v,2).subs(v,0)/2)
    expected=L*L*(L**4+L*L*m*m+(12*m*m-8*L*L)*p*p+12*p**4)/(32*a*a*r*r*W**9)
    residuals['full_fourth_order_kernel']=_zero(kernel-expected)
    residuals['coefficient_not_second_derivative']=_zero(energy-v*v*kernel)
    # Reuse Tr(H P4)=omega Tr(G4) and odd-Moyal trace cancellation from
    # the independently proved B_ref owner; no B or cu/cv evaluation occurs.
    massgap=sp.Symbol('M_gap',positive=True)
    moments={j:sp.simplify((massgap**(j-8)*sp.beta(sp.Rational(j+1,2),sp.Rational(8-j,2))).rewrite(sp.gamma)) for j in (0,2,4)}
    exact={0:sp.Rational(32,35)/massgap**8,2:sp.Rational(16,105)/massgap**6,4:sp.Rational(4,35)/massgap**4}
    residuals['standard_moments']=[_zero(moments[j]-exact[j])[0] for j in (0,2,4)]
    integral=D*L*L/(64*sp.pi*a*r*r)*(L*L*massgap**2*moments[0]+(12*m*m-8*L*L)*moments[2]+12*moments[4])
    closed=reference_quadratic_group(m,L,a,r,D,sp.pi)
    residuals['full_line_integral']=_zero(integral.subs(massgap,sp.sqrt(m*m+L*L))-closed)
    return {'symbols':(m,L,p,a,r,D),'kernel':kernel,'closed':closed,
            'P1_scalar_spatial_coefficient':g,'P2_point_scalar_coefficient':scalar2,
            'trace_G4_pieces':{'ordinary_P2_squared':sp.factor(ordinary/v**2),
                              'Moyal2_P0_P2_pair':sp.factor(moyal02/v**2),
                              'Moyal2_P1_P1':sp.factor(moyal11/v**2)},
            'residuals':residuals,'coefficient_definition':'[v²] E_N; v=r_Tz',
            'grade':'v² has grade4; no other background derivative can multiply it at this order',
            'measure':'D/(2pi), dk=a dp; D includes angular signs once',
            'scope':'C_ref only; full spatial Weyl recursion; no new state or reference quadrature'}


def archived_v_control(root):
    """Compare paired existing n24 nodes; never evaluate the old projector."""
    root=Path(root);record=json.loads((root/'results/development/nsc-incoming-joint-constraints.json').read_text())
    for field in ('source_hashes','input_hashes'):
        for path,digest in record[field].items():
            if sha256((root/path).read_bytes()).hexdigest()!=digest:raise ValueError('stored control owner changed: '+path)
    spec=record['payload']
    if sha256((root/spec['path']).read_bytes()).hexdigest()!=spec['sha256']:raise ValueError('stored n24 payload changed')
    change=record['reference_controls']['n24']['changed_normal_entries']
    if len(change)!=1 or (change[0]['field'],change[0]['normal_order'],change[0]['spatial_order'],change[0]['delta'])!=('r',1,1,.01):
        raise ValueError('original v=.01 control required')
    channels=json.loads((root/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels']
    a=np.sqrt(3*np.pi/2-4);r=np.sqrt(2);v=.01
    proof=reference_quadratic_lapse();m,L,p,aa,rr,D=proof['symbols']
    kernel=sp.lambdify((m,L,p,aa,rr),proof['kernel'],'numpy')
    maximum=0.;groups={};momentum={}
    with np.load(root/spec['path'],allow_pickle=False) as archive:
        meta=json.loads(archive['metadata_json'].tobytes())
        for i,row in enumerate(meta['controls']['n24']['signed_inventory']):
            k=archive[f'n24/{i}/momenta'];values=archive[f'n24/{i}/insertions'][:,:,0].sum(axis=0);n=len(k)//2
            if not np.array_equal(k[:n],-k[n:]):raise ValueError('actual both-sign momentum grid required')
            group=row['group'];folded=(values[:n]+values[n:])/2
            if group in momentum and not np.array_equal(momentum[group],k[n:]):raise ValueError('angular grids differ')
            momentum[group]=k[n:];groups.setdefault(group,[]).append(folded)
        for group,rows in groups.items():
            c=channels[group];expected_signs=1 if c['angular_eigenvalue']==0 else 2
            if len(rows)!=expected_signs:raise ValueError('complete signed family required before comparison')
            paired=sum(rows)/expected_signs
            expected=v*v*kernel(c['compact_mass'],c['angular_eigenvalue']/r,momentum[group]/a,a,r)
            maximum=max(maximum,float(np.max(abs(paired-expected))))
    closed=sum(reference_quadratic_group(c['compact_mass'],c['angular_eigenvalue']/r,a,r,c['copy_count']*c['degeneracy'],np.pi) for c in channels[1:])
    old=record['reference_controls']['n24']['reference_action_gradient_change'][0]
    return {'paired_node_absolute_error':maximum,'aggregate_action_error':abs(v*v*closed-old),
            'closed_reference_coefficient':closed,'old_paired_control_coefficient':old/(v*v),
            'control_tolerance':record['numerical_tolerance'],
            'angular_and_momentum_signs_combined_before_comparison':True,
            'new_reference_evaluations':0,'new_momentum_quadratures':0,'payload':spec}
