"""Weyl-disk control of the selected smooth-motif spectral tail.

This is a fixed finite spatial discretization with H_n=Omega^n H and
B_n=Omega^n B in an unweighted direct-sum Hilbert space. It does not identify
the projected resolvent with a dark-energy density or a cosmological fraction.
"""
from __future__ import annotations

import mpmath as mp
import numpy as np
from .nsc_smooth_geometry import smooth_motif


def motif_hoppings(intervals=16):
    motif=smooth_motif(1.,2.,intervals)
    forward=motif['H'][intervals:,:intervals]
    positive=np.empty(2*intervals)
    positive[::2]=np.abs(np.diag(forward))
    positive[1:-1:2]=np.diag(forward,1)
    positive[-1]=motif['B'][-1,0]
    if np.any(positive<=0):raise ValueError('strictly positive Jacobi couplings required')
    return motif,positive


def coefficients(base,omega,depth):
    if omega<=0 or depth<1:raise ValueError('positive scale ratio and depth required')
    size=len(base)*depth
    # mp.mpf(float) retains the exact binary64 input at sufficient precision.
    # High-precision bounds concern this declared lattice, not exact continuum
    # geometry or the physical precision of any cosmological measurement.
    return [mp.mpf(float(base[j%len(base)]))*mp.mpf(omega)**(j//len(base))
            for j in range(size-1)]


def finite_response(hoppings,z,terminal_shift=0):
    response=1/(z-terminal_shift)
    for a in reversed(hoppings):response=1/(z-a*a*response)
    return response


def weyl_disk(hoppings,z):
    """Exact radius formula evaluated with the caller's mpmath precision.

    p_-1=0,p_0=1; s_0=0,s_1=1/a_0. The second solution is the
    inhomogeneous resolvent fundamental solution, not an extra left boundary.
    For Im(z)>0 the convention (z-J)^-1 has negative imaginary part.
    """
    if mp.im(z)<=0:raise ValueError('upper-half-plane energy required')
    pprev=mp.mpc(0);p=mp.mpc(1)
    sprev=mp.mpc(0);s=mp.mpc(0)
    psum=mp.mpf(1);ssum=mp.mpf(0);cross=mp.mpc(0)
    for j,a in enumerate(hoppings):
        prev_a=hoppings[j-1] if j else 0
        pn=(z*p-prev_a*pprev)/a
        sn=1/a if j==0 else (z*s-prev_a*sprev)/a
        pprev,p=p,pn;sprev,s=s,sn
        psum+=abs(p)**2;ssum+=abs(s)**2;cross+=p*mp.conj(s)
    radius=1/(2*mp.im(z)*psum)
    center=mp.conj(cross)/psum-1j*radius
    return center,radius,(psum,ssum,cross)


def record_family(omega,depths=(1,2,4,8,16,32,64),intervals=16,precision=220):
    motif,base=motif_hoppings(intervals)
    rows=[]
    with mp.workdps(precision):
        z=mp.mpc('.3','.4')
        for depth in depths:
            a=coefficients(base,omega,depth)
            center,radius,_=weyl_disk(a,z)
            response=finite_response(a,z)
            circumference_residual=abs(abs(response-center)-radius)/radius
            # A finite real end potential and a passive terminal self-energy
            # exercise the same disk independently of the Dirichlet choice.
            shifted=finite_response(a,z,mp.mpf(2))
            passive=finite_response(a,z,-mp.j)
            shifted_residual=abs(abs(shifted-center)-radius)/radius
            passive_radius_ratio=abs(passive-center)/radius
            if max(circumference_residual,shifted_residual)>mp.mpf('1e-35'):
                raise RuntimeError('Weyl circumference identity failed')
            if passive_radius_ratio>1+mp.mpf('1e-35'):
                raise RuntimeError('passive self-energy escaped the Weyl disk')
            rows.append({'depth':depth,'sites':len(base)*depth,
                         'dirichlet_response':[float(mp.re(response)),float(mp.im(response))],
                         'center':[float(mp.re(center)),float(mp.im(center))],
                         'log10_diameter':float(mp.log10(2*radius)),
                         'diameter_numeric':float(2*radius),
                         'dirichlet_circle_relative_residual':float(circumference_residual),
                         'real_terminal_circle_relative_residual':float(shifted_residual),
                         'passive_terminal_radius_ratio':float(passive_radius_ratio)})
    q=motif['zero_energy_transfer']
    return {'Omega':omega,'q':q,'critical_Omega':1/q,
            'zero_energy_even_multiplier_abs':q,
            'zero_energy_odd_multiplier_abs':1/(omega*q),
            'endpoint_class':'limit_point' if omega*q<=1 else 'limit_circle',
            'self_adjoint_condition_at_infinite_depth_needed':omega*q>1,
            'precision_decimal_digits':precision,'rows':rows}
