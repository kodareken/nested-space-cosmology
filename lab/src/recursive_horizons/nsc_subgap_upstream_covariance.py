"""Covariance accuracy of an unchanged, archived homogeneous preparation.

The signed partner construction applies only before the inhomogeneous history.
It is not a pointwise unitary assumption for the later evolved field.
"""
import numpy as np
from flint import acb, arb, ctx

from .nsc_subgap_source_covariance import original_covariance_error


def negative_subgap_covariance_error(columns, covariance, positive_validation):
    """Compare original negative columns with the exact complementary state.

    In the homogeneous subgap preparation the exact two-row mode map has
    Gram I, including the matched unitary interior frame. Hence
    Qminus = I - S3 conjugate(Qplus) S3 at opposite angular sign.
    Its Bloch vector is (nx,-ny,-nz), an error-norm isometry.
    The numerical negative covariance is evaluated independently; its Gram
    or source-complement arithmetic is never assumed exact.
    Caller binds the original opposite-energy/opposite-angular pair.
    """
    with ctx.workprec(positive_validation['bits']):
        negative = dict(positive_validation)
        endpoint = positive_validation['endpoint']
        if len(endpoint) != 3:
            raise ValueError('complete positive Bloch endpoint required')
        negative['endpoint'] = (endpoint[0], -endpoint[1], -endpoint[2])
        return original_covariance_error(columns, covariance, negative)


def numerical_gram_error(columns, *, bits=192):
    """Directed Frobenius upper bound on A A^dagger-I, for diagnostics."""
    values = np.asarray(columns, complex)
    if values.shape != (2, 3) or not np.isfinite(values).all():
        raise ValueError('finite two-row three-column preparation required')
    with ctx.workprec(bits):
        matrix = [[acb(v) for v in row] for row in values]
        squared = arb(0)
        for i in range(2):
            for j in range(2):
                entry = sum((matrix[i][k]*matrix[j][k].conjugate()
                             for k in range(3)), acb(0)) - int(i == j)
                squared += entry.abs_upper()**2
        return squared.sqrt().upper()
