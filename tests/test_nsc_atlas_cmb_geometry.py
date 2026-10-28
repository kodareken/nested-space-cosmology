import math
import numpy as np
import pytest
from scipy.integrate import quad
from recursive_horizons.nsc_atlas_cmb_geometry import (
    early_parent_bound,dust_transverse_distance,profiled_distance_prior,PRIOR_MEAN,PRIOR_COVARIANCE,
)


def test_high_redshift_integral_matches_direct_independent_form():
    m=.304;z=1089.92
    distance,_=dust_transverse_distance('flat_LCDM',m,z)
    direct=quad(lambda x:1/math.sqrt(m*(1+x)**3+1-m),0,z,epsabs=1e-11)[0]
    assert distance==pytest.approx(direct,abs=2e-11)


def test_profile_sound_scale_retains_correct_marginal_R_variance():
    for shift in (-.15,0,.02):
        row=profiled_distance_prior(PRIOR_MEAN[0]+shift)
        assert row['conditional_chi_squared_lower_bound']==pytest.approx(shift**2/PRIOR_COVARIANCE[0,0],abs=1e-10)
        inverse=np.linalg.inv(PRIOR_COVARIANCE)
        residual=np.array([shift,row['profiled_lA']-PRIOR_MEAN[1]])
        assert (inverse@residual)[1]==pytest.approx(0.,abs=1e-10)


def test_early_scalar_control_closes_without_an_early_thermal_claim():
    bound=early_parent_bound(1.4656003508,1089.92)
    assert bound['minimum_scale']<.002
    assert bound['source_density_relative_error']<.0002
    assert bound['causal_margin_lower']>.1
