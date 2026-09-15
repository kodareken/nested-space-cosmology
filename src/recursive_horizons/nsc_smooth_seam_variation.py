"""Geometric boundary variation on the actual shared smooth NSC seam.

The normal reversal and boundary first-variation identities are imported.
This owner applies the locked local action to the same metric jet on the
two sides of rho=0. It neither assigns Gamma_rest=0 nor solves bulk sourcing.
"""
import numpy as np


def reverse_normal(jet):
    result=dict(jet)
    for key in ('Nd','ad','addd','rd','rddd'):
        result[key]=-result[key]
    return result


def boundary_data(jet,coefficients):
    """Outward-normal coefficients per common axial coordinate.

    Derivatives in `jet` refer to an outward-oriented coordinate with lapse N.
    The returned proper Weyl covector acts on (delta a,delta r,delta sigma),
    sigma=H_a-H_r. No second orientation sign is applied to this covector.
    """
    N,Nd,Ndd,a,ad,add,addd,r,rd,rdd,rddd=(jet[k] for k in
         ('N','Nd','Ndd','a','ad','add','addd','r','rd','rdd','rddd'))
    if min(N,a,r)<=0 or not np.isfinite(list(jet.values())).all():
        raise ValueError('finite normal metric jets with positive N,a,r required')
    A=(add-ad*Nd/N)/(N*N*a)
    Ad=(addd-add*Nd/N-ad*Ndd/N+ad*Nd*Nd/(N*N))/(N*N*a)-A*(2*Nd/N+ad/a)
    B=ad*rd/(N*N*a*r)
    Bd=(add*rd+ad*rdd)/(N*N*a*r)-B*(2*Nd/N+ad/a+rd/r)
    C=(rdd-rd*Nd/N)/(N*N*r)
    Cd=(rddd-rdd*Nd/N-rd*Ndd/N+rd*Nd*Nd/(N*N))/(N*N*r)-C*(2*Nd/N+rd/r)
    D=(1+(rd/N)**2)/(r*r)
    Dd=2*rd*rdd/(N*N*r*r)-2*rd*rd*Nd/(N**3*r*r)-2*D*rd/r
    F=A-B-C+D;Fd=Ad-Bd-Cd+Dd
    scalar=-2*(A+2*B+2*C+D);scalar_d=-2*(Ad+2*Bd+2*Cd+Dd)
    k=-16*np.pi*coefficients['C_Weyl']/3
    Qa=2*k*r*r*F/N;Qr=-2*k*a*r*F/N
    Qad=2*k*((2*r*rd*F+r*r*Fd)/N-r*r*F*Nd/N**2)
    Qrd=-2*k*((ad*r*F+a*rd*F+a*r*Fd)/N-a*r*F*Nd/N**2)
    PN=-(ad*Qa+rd*Qr)/N
    Pa=-Qa*(Nd/N+rd/r)-Qad
    Pr=Qr*(ad/a-Nd/N-2*rd/r)-Qrd
    weyl=np.array([Pa+ad/a*Qa,Pr+rd/r*Qr,N*a*Qa])
    einstein=np.array([-16*np.pi*coefficients['A']*r*rd/N,
                       -16*np.pi*coefficients['A']*(r*ad+a*rd)/N])
    euler=-32*np.pi*coefficients['C_Euler']*(ad/N)*(1+(rd/N)**2)
    box=-4*np.pi*coefficients['C_boxR']*a*r*r*scalar_d/N
    return {'weyl':weyl,'einstein_ghy':einstein,'sigma':ad/(N*a)-rd/(N*r),
            'Q_a':Qa,'Q_r':Qr,'proper_lapse':PN+ad/N*Qa+rd/N*Qr,
            'Weyl_scalar_F':F,'normal_F_derivative':Fd/N,
            'scalar_curvature':scalar,'normal_scalar_derivative':scalar_d/N,
            'euler_boundary_density':euler,'boxR_boundary_density':box}


def smooth_seam_match(common_jet,coefficients):
    """Parent future-normal / child past-normal pullback on ONE smooth cut."""
    parent=boundary_data(common_jet,coefficients)
    child=boundary_data(reverse_normal(common_jet),coefficients)
    # Intrinsic metric variations agree; outward-normal shear variations oppose.
    child_pullback=np.diag([1.,1.,-1.])
    residual=parent['weyl']+child_pullback@child['weyl']
    return {'parent':parent,'child':child,'child_variation_pullback':child_pullback,
      'residuals':{'weyl_common_covector':float(np.max(abs(residual))),
        'einstein_ghy_common_covector':float(np.max(abs(parent['einstein_ghy']+child['einstein_ghy']))),
        'normal_shear_matching':float(abs(parent['sigma']+child['sigma'])),
        'euler_boundary_sum':float(abs(parent['euler_boundary_density']+child['euler_boundary_density'])),
        'boxR_boundary_sum':float(abs(parent['boxR_boundary_density']+child['boxR_boundary_density']))},
      'quantum_boundary_remainder':None,'bulk_metric_residual':None}
