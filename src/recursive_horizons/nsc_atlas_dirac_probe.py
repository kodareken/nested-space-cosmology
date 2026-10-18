"""Geometry inputs for a massless Dirac probe on the original cap and wall.

This extends the classical scalar candidate by a declared probe. It does
not choose a quantum state, compute a loop polarization or predict Psi/Phi.
The smooth control continues the same de Sitter cap across the cut. A flat
control would change the cap curvature as well as remove the interface.
"""
from __future__ import annotations

from math import atanh,isfinite,sqrt

from .nsc_atlas_cap_wall_response import (
    H_ELL,K_ELL,M_OVER_ELL,reference_coefficients,
)


def cap_radial_dirac(radius, *, hubble=H_ELL, angular_k=1):
    """H_k=-i sigma_2 d/dz+sigma_1 W, z=atanh(H*r)/H.

    The physical spinor is A^(-1/4)*r^(-1) times its radial column, and
    W=k*sqrt(A)/r. The two eigenchannels of H_k^2 have W^2 +/- dW/dz.
    All lengths here are in ell; k=+/-1,+/-2,... is a spinor angular label.
    The flat values are at the same areal radius and central-clock rate.
    """
    r,h=float(radius),float(hubble)
    if (not isfinite(r) or not isfinite(h) or r<=0 or h<0 or h*r>=1
            or isinstance(angular_k,bool) or not isinstance(angular_k,int) or angular_k==0):
        raise ValueError('require 0<r, 0<=H*r<1 and a nonzero integer spinor angular label')
    k=angular_k
    A=1-(h*r)**2
    W=k*sqrt(A)/r
    derivative=-k*sqrt(A)/r**2
    return dict(radius=r,hubble=h,angular_k=k,A=A,
        tortoise_z=atanh(h*r)/h if h else r,
        W=W,dW_dz=derivative,V_plus=W*W+derivative,V_minus=W*W-derivative,
        flat_W=k/r,flat_V_plus=k*(k-1)/r**2,flat_V_minus=k*(k+1)/r**2,
        W_over_flat_at_same_radius=sqrt(A),
        spinor_rescaling=A**(-.25)/r,
        bulk_operator_same_as_smooth_cap_control=True)


def reference_dirac_interface(radius, *, angular_k=1):
    """Transparent, zero-surface-Dirac-interaction domain in a common frame.

    A common normal points from cap to exterior and [X]=X_plus-X_minus.
    The first-order spinor is continuous. In the squared operator's domain,
    [d_n psi]=-[K]*psi/2 in a frame transported along that normal. The
    scalar coefficient is not a reflection amplitude or gravitational slip.
    """
    if isinstance(angular_k,bool) or not isinstance(angular_k,int) or angular_k==0:
        raise ValueError('a nonzero integer spinor angular label is required')
    g=reference_coefficients(radius)
    r=float(g['x'])
    bp,bm=float(g['beta_plus']),float(g['beta_minus'])
    kp=float(g['beta_plus_x'])+2*bp/r
    km=float(g['beta_minus_x'])+2*bm/r
    jump=kp-km
    return dict(radius=r,beta_plus=bp,beta_minus=bm,
        trace_K_plus=kp,trace_K_minus=km,trace_K_jump=jump,
        squared_Dirac_normal_derivative_jump_coefficient=-jump/2,
        angular_inverse_radius_normal_jump=-angular_k*(bp-bm)/r**2,
        smooth_control_trace_K_jump=0.,smooth_control_normal_derivative_jump=0.,
        exterior_f=1-2*M_OVER_ELL/r,
        exterior_static_time_is_timelike=1-2*M_OVER_ELL/r>0,
        domain='spinor and D_spinor continuous in a common nonsingular wall frame',
        added_premise='massless probe with no direct surface spinor interaction')


def quiet_reference_curvatures(radius):
    """One-sided smooth bulk invariants only; no products of surface deltas."""
    reference_coefficients(radius)
    r=float(radius)
    return dict(cap=dict(R=12*H_ELL**2,Ricci_squared=36*H_ELL**4,
                        Riemann_squared=24*H_ELL**4,Weyl_squared=0.),
                exterior=dict(R=0.,Ricci_squared=0.,
                    Riemann_squared=48*M_OVER_ELL**2/r**6,
                    Weyl_squared=48*M_OVER_ELL**2/r**6))


def reference_probe_report(radius=.5):
    return dict(radius=radius,units='ell=1',cap=cap_radial_dirac(radius),
        interface=reference_dirac_interface(radius),curvatures=quiet_reference_curvatures(radius),
        material_parameter_k_ell=K_ELL,nested_Psi_over_Phi=None,
        reason='No Dirac two-point state, renormalized interface response or conserved metric-probe/observer protocol is specified for this cap-wall model.',
        calculation_scope='Original zero-scalar reference and a declared transparent massless probe; not the unconverged finite-amplitude P001 metric.',
        nesting_changes='The common-frame interface domain and exterior causal response. The same cap bulk operator is shared with the smooth de Sitter control.',
        quantum_inputs_still_needed=[
            'Initial/boundary Dirac two-point function on the joined finite causal domain, including incoming exterior occupations.',
            'Renormalization conditions for bulk and interface terms in the probe extension.',
            'A conserved perturbing source, wall displacement response and gauge-invariant clock/motion/lensing readout on this time-dependent background.',
        ])
