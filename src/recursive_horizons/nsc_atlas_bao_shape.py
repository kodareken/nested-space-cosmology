"""Sound-ruler-free Alcock--Paczynski null contrasts of Gaussian BAO data."""
import math

import numpy as np
from scipy.linalg import cho_factor,cho_solve,solve

from .nsc_atlas_parent_supernova import flrw_distance


def ap_prediction(redshifts,*,matter,curvature=0.):
    z=np.asarray(redshifts,float)
    DM=flrw_distance(z,z,matter=matter,curvature=curvature)/(1+z)
    E=np.sqrt(matter*(1+z)**3+curvature*(1+z)**2+1-matter-curvature)
    return DM*E


def ap_null_contrast(mean,covariance,pairs,predicted):
    """Test y_M-F*y_H=0 without Gaussianising the observed ratio y_M/y_H.

    All within/between-bin covariance is retained. F is fixed independently
    of these BAO observations. This is conditional on that model calibration.
    """
    y=np.asarray(mean,float);C=np.asarray(covariance,float);F=np.asarray(predicted,float)
    pairs=np.asarray(pairs,int)
    if (y.ndim!=1 or C.shape!=(len(y),len(y)) or pairs.shape!=(len(F),2)
            or np.any(pairs<0) or np.any(pairs>=len(y)) or len(set(pairs.ravel()))!=pairs.size
            or np.any(F<=0) or not np.isfinite(F).all() or not np.isfinite(C).all()
            or not np.isfinite(y).all() or np.max(abs(C-C.T))>1e-12):
        raise ValueError('aligned finite mean, symmetric covariance, unique pairs and positive AP predictions required')
    W=np.zeros((len(F),len(y)))
    W[np.arange(len(F)),pairs[:,0]]=1.
    W[np.arange(len(F)),pairs[:,1]]=-F
    residual=W@y;V=W@C@W.T
    factor=cho_factor(V,lower=True)
    value=float(residual@cho_solve(factor,residual))
    independent=float(residual@solve(V,residual,assume_a='sym'))
    if abs(value-independent)>1e-9*(1+value):
        raise ValueError('independent AP covariance solve disagrees')
    return dict(chi_squared=value,degrees_of_freedom=len(F),
                linear_residuals=residual.tolist(),standardized_marginal_residuals=(residual/np.sqrt(np.diag(V))).tolist(),
                marginal_standard_uncertainties=np.sqrt(np.diag(V)).tolist(),
                independent_chi_squared_difference=abs(value-independent))
