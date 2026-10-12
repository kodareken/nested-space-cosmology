"""Numerical insertion errors for one unchanged finite coherent source table."""
from functools import reduce
from flint import acb,arb,ctx
import numpy as np

from .nsc_ks_ball_trajectory import complex_ball,exact_upper,restored_upper


def source_norm_upper(covariance, *, bits=120):
    """sqrt(||C||_1 ||C||_infinity), retaining all source coherences."""
    matrix=np.asarray(covariance,complex)
    if matrix.ndim!=2 or matrix.shape[0]!=matrix.shape[1] or not matrix.shape[0] or not np.isfinite(matrix).all():
        raise ValueError('finite square coherent source matrix required')
    with ctx.workprec(bits):
        magnitude=[[complex_ball(v).abs_upper() for v in row] for row in matrix]
        row=reduce(arb.max,(sum(values,arb(0)).upper() for values in magnitude))
        column=reduce(arb.max,(sum((magnitude[i][j] for i in range(len(matrix))),arb(0)).upper() for j in range(len(matrix))))
        return (row*column).sqrt().upper()


def endpoint_norms(field,energies):
    """Uniform weighted Frobenius norms of the reconstructed F and actual F_z."""
    if len(energies)!=field.source_count:raise ValueError('matching source energies required')
    with ctx.workprec(field.bits):
        rows=[arb(0),arb(0)];derivative=[arb(0),arb(0)]
        for spin in range(2):
            for k,wave in enumerate(field.waves):
                values=[sum((field.coefficients[p][spin][s][k] for p in range(8)),acb(0)) for s in range(field.source_count)]
                rows[spin]+=sum((v.abs_upper()**2 for v in values),arb(0)).upper().sqrt()
                derivative[spin]+=sum(((acb(0,wave-arb(float(E)))*v).abs_upper()**2
                                      for E,v in zip(energies,values)),arb(0)).upper().sqrt()
        return (sum((v*v for v in rows),arb(0)).sqrt().upper(),
                sum((v*v for v in derivative),arb(0)).sqrt().upper())


def finite_matter_error(field_norm,axial_norm,field_error,axial_error,source_norm,*,
                        mass,absolute_angular,axial_lower,radius_lower,multiplicity,bits=120):
    """Same archived C for exact/reconstructed columns; source error is separate."""
    with ctx.workprec(bits):
        values=list(map(arb,(field_norm,axial_norm,field_error,axial_error,source_norm,
                             mass,absolute_angular,axial_lower,radius_lower,multiplicity)))
        if any(not v.is_finite() or not v>=0 for v in values) or any(not v>0 for v in values[7:]):
            raise ValueError('finite nonnegative bounds and positive geometry/multiplicity required')
        f,z,e,ez,C,m,ell,a,r,mu=values
        density=C*(2*f*e+e*e)
        current=C*(z*e+f*ez+e*ez)
        potential=(m*m+(ell/r)**2).sqrt()
        return {'N':exact_upper(mu*(potential*density+current/a)),
                'beta':exact_upper(mu*current),
                'source_accuracy_included':False,'physical_local_gate':'OPEN'}
