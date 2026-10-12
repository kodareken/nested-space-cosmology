"""Reference-column error from a causally unaffected point of the same solve."""
import numpy as np
from flint import arb,ctx

from .nsc_ks_ball_trajectory import complex_ball,exact_upper,restored_upper


def spectator_reference_error(initial_rho, period_origin, period_length,
                              support_lower, support_upper, endpoint_difference,
                              column_weights, energies, row_field_error, *, reference_pair=None,bits=120):
    """Bound A_ref error without another ODE or an empirical drift correction.

    The exact comparison solution at the period origin equals the homogeneous
    reference when its entire backward cone misses every periodic copy of
    the compact perturbation. The available field error is uniform over the
    numerical period. Physical source preparation error is separate.
    """
    difference=np.asarray(endpoint_difference,complex)
    weights=np.asarray(column_weights,float);labels=np.asarray(energies,float)
    if difference.shape!=(2,len(weights)) or labels.shape!=weights.shape or not len(weights):
        raise ValueError('one weighted source fiber per reference endpoint column required')
    if not all(np.isfinite(v).all() for v in (difference,weights,labels)) or np.any(weights<=0):
        raise ValueError('finite columns/labels and positive original weights required')
    pair=None
    if reference_pair is not None:
        pair=tuple(np.asarray(v,complex) for v in reference_pair)
        if len(pair)!=2 or any(v.shape!=difference.shape or not np.isfinite(v).all() for v in pair):
            raise ValueError('matching joint and cached reference arrays required')
    with ctx.workprec(bits):
        up,origin,length,lo,hi=map(arb,(initial_rho,period_origin,period_length,support_lower,support_upper))
        if not up>=1 or not up<=arb(33)/32 or not length>0 or not origin<lo or not lo<hi or not hi<origin+length:
            raise ValueError('owned trapped slab and interior compact support required')
        distance=(up-1)/(arb(4)/5)**2
        margin=arb.min(lo-origin,origin+length-hi)-distance
        if not margin>0:raise ValueError('spectator backward cone reaches a metric perturbation')
        field_error=arb(row_field_error)
        if not field_error.is_finite() or not field_error>=0:raise ValueError('uniform numerical row error bound required')
        leakage=[]
        for spin,row in enumerate(difference):
            squares=[]
            for source,(value,weight) in enumerate(zip(row,weights)):
                value=complex_ball(value)
                if pair is not None:
                    value+=complex_ball(pair[0][spin,source])-complex_ball(pair[1][spin,source])
                squares.append((value*arb(float(weight))).abs_upper()**2)
            leakage.append(sum(squares,arb(0)).upper().sqrt().upper())
        error=(field_error+arb.max(*leakage)).upper()
        epsilon=arb(2).sqrt()*error
        return {'row_reference_error_upper':exact_upper(error),
                'weighted_reference_F_error_upper':exact_upper(epsilon),
                'weighted_reference_F_z_error_upper':exact_upper(float(np.max(abs(labels)))*epsilon),
                'numerical_difference_at_spectator_upper':exact_upper(arb.max(*leakage)),
                'causal_margin_lower':{'mantissa':str(margin.lower().man_exp()[0]),'exponent':int(margin.lower().man_exp()[1])},
                'reference_pair_difference_retained':pair is not None,
                'physical_source_error_included':False,'stress_drift_subtracted':False}


def homogeneous_reference_norms(amplitudes,weights,energies,*,bits=120):
    A=np.asarray(amplitudes,complex);w=np.asarray(weights,float);E=np.asarray(energies,float)
    if A.shape!=(2,len(w)) or E.shape!=w.shape or not len(w):raise ValueError('matching reference columns required')
    if not all(np.isfinite(v).all() for v in (A,w,E)) or np.any(w<=0):raise ValueError('finite weighted reference required')
    with ctx.workprec(bits):
        norm,axial=arb(0),arb(0)
        for row in A:
            for value,weight,energy in zip(row,w,E):
                upper=(complex_ball(value)*arb(float(weight))).abs_upper()
                norm+=upper**2;axial+=(arb(float(energy))*upper)**2
        return norm.upper().sqrt().upper(),axial.upper().sqrt().upper()
