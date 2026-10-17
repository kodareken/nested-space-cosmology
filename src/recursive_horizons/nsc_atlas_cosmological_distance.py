"""A scale-free distance contrast for the existing quiet de Sitter cap.

This is a fixed linear statistic on standardized apparent magnitudes, not a
cosmological parameter fit or the nonzero-source P001 forward model.
"""
import math

import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve

LOW_BAND = (.025, .075)
HIGH_BAND = (.15, .25)


def quiet_cap_distance_shape(z_hd, z_hel):
    """D_L/(c/H) with the released heliocentric prefactor convention."""
    z, hel = np.broadcast_arrays(np.asarray(z_hd, dtype=float),
                                 np.asarray(z_hel, dtype=float))
    if (not np.all(np.isfinite(z)) or not np.all(np.isfinite(hel))
            or np.any(z <= 0) or np.any(hel <= -1)):
        raise ValueError('finite positive zHD and zHEL>-1 required')
    return z*(1+hel)


def quiet_cap_domain():
    """Finite causal-domain argument at t=0, Rwall/ell=1.5; no wall solve.

    On 1<=r<=1.5 the separate terms of Rdot^2 give 4569/1600<4.
    Since dt/dtau>=1, the wall takes >1/4 cap-clock units to reach r=1
    backwards. The tested null cone has lookback <1/4 and radius <=1/5.
    """
    maximum_z = HIGH_BAND[1]
    return dict(observer_cap_time_over_ell=0., reference_wall_radius_over_ell=1.5,
                reference_parameters=dict(H_ell=1., k_ell=2., M_over_ell=.4),
                physical_frame='Jordan and Einstein metrics coincide at zero scalar',
                wall_speed_squared_upper=4569/1600,
                wall_horizon_crossing_lookback_lower=.25,
                emission_lookback_upper=math.log1p(maximum_z),
                emission_radius_upper_over_ell=maximum_z/(1+maximum_z),
                source='zero scalar mean; comoving test emitters; no probe backreaction')


def contrast_weights(covariance, groups):
    """Joint-covariance GLS contrast, independently checked by a KKT solve.

    Each row belongs to band0 or band1. Covariance is symmetrized only for
    antisymmetric floating-point discrepancies <=1e-12 times max(1,max|C|).
    No diagonal jitter or decorrelation is applied.
    """
    C = np.asarray(covariance, dtype=float)
    group = np.asarray(groups)
    count = group.size
    if (group.ndim != 1 or count < 2 or not np.all(np.isin(group, (0, 1)))
            or not np.any(group == 0) or not np.any(group == 1)
            or C.shape != (count, count) or not np.all(np.isfinite(C))):
        raise ValueError('two populated bands and their finite square covariance required')
    asymmetry = float(np.max(np.abs(C-C.T)))
    if asymmetry > 1e-12*max(1., float(np.max(np.abs(C)))):
        raise ValueError('covariance is not symmetric to numerical precision')
    C = (C+C.T)/2
    design = np.column_stack((group == 0, group == 1)).astype(float)
    target = np.array([-1., 1.])
    try:
        precision_design = cho_solve(cho_factor(C, lower=True), design)
        normal = design.T @ precision_design
        dual = solve(normal, target, assume_a='pos')
        weights = precision_design @ dual
        variance = float(target @ dual)
        kkt = np.block([[C, design], [design.T, np.zeros((2, 2))]])
        independent = solve(kkt, np.r_[np.zeros(count), target], assume_a='sym')[:count]
    except np.linalg.LinAlgError as error:
        raise ValueError('unresolved positive covariance or contrast solve') from error
    weight_difference = float(np.max(np.abs(weights-independent)))
    group_residual = float(np.max(np.abs(design.T @ weights-target)))
    variance_difference = abs(float(weights @ C @ weights)-variance)
    if (not np.isfinite(variance) or variance <= 0
            or weight_difference > 1e-10*(1+float(np.max(np.abs(weights))))
            or group_residual > 1e-10 or variance_difference > 1e-10*(1+variance)):
        raise ValueError('independent contrast or zero-point check failed')
    checks = dict(covariance_asymmetry=asymmetry,
                  KKT_weight_maximum_difference=weight_difference,
                  band_normalization_maximum_residual=group_residual,
                  variance_identity_difference=variance_difference,
                  common_zero_point_weight_sum=float(weights.sum()))
    return weights, variance, checks


def project_distance_contrast(weights, z_hd, z_hel, distances_over_common_scale):
    """Project a forward model at the SAME selected rows, not band midpoints."""
    w = np.asarray(weights, dtype=float)
    distance = np.asarray(distances_over_common_scale, dtype=float)
    quiet = quiet_cap_distance_shape(z_hd, z_hel)
    if (w.ndim != 1 or distance.shape != w.shape or quiet.shape != w.shape
            or not np.all(np.isfinite(w)) or not np.all(np.isfinite(distance))
            or np.any(distance <= 0) or abs(float(w.sum())) > 1e-10):
        raise ValueError('positive model distances and zero-sum weights in selected-row order required')
    return float(w @ (5*np.log10(distance/quiet)))


def distance_contrast(z_hd, z_hel, magnitude, covariance, calibrator):
    """Select the declared bands and measure their contrast in magnitudes.

    Returned indices/weights are for in-memory reuse. The producer saves
    aggregate statistics and input hashes, rather than those raw vectors.
    """
    z, hel, mag = (np.asarray(value, dtype=float) for value in (z_hd, z_hel, magnitude))
    flags = np.asarray(calibrator)
    C = np.asarray(covariance, dtype=float)
    if (z.ndim != 1 or z.shape != hel.shape or z.shape != mag.shape or flags.shape != z.shape
            or not np.all(np.isfinite(mag)) or not np.all(np.isin(flags, (0, 1)))
            or C.shape != (z.size, z.size)):
        raise ValueError('aligned redshifts, magnitudes, calibrator flags and covariance required')
    quiet_cap_distance_shape(z, hel)
    low = (z >= LOW_BAND[0]) & (z < LOW_BAND[1])
    high = (z >= HIGH_BAND[0]) & (z <= HIGH_BAND[1])
    indices = np.flatnonzero((low | high) & (flags == 0))
    groups = high[indices].astype(int)
    weights, variance, checks = contrast_weights(C[np.ix_(indices, indices)], groups)
    y = mag[indices]-5*np.log10(quiet_cap_distance_shape(z[indices], hel[indices]))
    delta = float(weights @ y)
    sigma = math.sqrt(variance)
    ratio = 10**(delta/5)
    interval = [10**((delta-sigma)/5), 10**((delta+sigma)/5)]
    statistics = dict(selected_light_curves=int(indices.size),
                      band_light_curve_counts=[int(np.sum(groups == band)) for band in (0, 1)],
                      contrast_magnitudes=delta, standard_uncertainty_magnitudes=sigma,
                      exponential_contrast=ratio, exponential_contrast_68_percent_interval=interval,
                      quiet_cap_predicted_contrast_magnitudes=0.,
                      quiet_cap_predicted_exponential_contrast=1.,
                      quiet_cap_standardized_contrast=delta/sigma,
                      quiet_cap_two_sided_Gaussian_p=math.erfc(abs(delta)/(math.sqrt(2)*sigma)),
                      checks=checks)
    return indices, weights, statistics
