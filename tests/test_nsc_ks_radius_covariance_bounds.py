"""The spatial profile binding alone cannot authorize a wider normal window."""
import importlib.util
from pathlib import Path
import pytest
pytest.importorskip('flint')
from recursive_horizons.nsc_fourier_covariance_moments import radius_control as module


def test_actual_normal_window_must_fit_the_bound_used_for_spatial_norms():
    from fractions import Fraction
    geometry={'normal_support':{'exact_rational':str(Fraction(.03))}}
    direction={'normal_inner':float(.007).hex(),'normal_outer':float(.03).hex()}
    history={'history_description':{'directions':[direction]}}
    module.require_normal_support(geometry,history)
    direction['normal_outer']=float(.035).hex()
    with pytest.raises(ValueError,match='normal support'):
        module.require_normal_support(geometry,history)
