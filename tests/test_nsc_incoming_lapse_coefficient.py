import json
from pathlib import Path

import numpy as np

from recursive_horizons.nsc_incoming_lapse_coefficient import reference_identity,closed_lapse_coefficient,local_density_stencil
from recursive_horizons.nsc_incoming_cauchy_jets import incoming_cauchy_jets
from recursive_horizons.nsc_incoming_local_constraints import scalar_taylor_coefficients,taylor_two_jets,local_action_densities
from recursive_horizons.nsc_light_restoration_action import LightRestorationAction
from recursive_horizons.nsc_magnetic_light_reference import magnetic_light_spectrum
from recursive_horizons.nsc_spherical_local_history import LockedSphericalLocalAction

ROOT=Path(__file__).resolve().parents[1]


def inputs():
    c=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels']
    l=json.loads((ROOT/'results/development/nsc-incoming-local-constraints.json').read_text())['locked_inputs']
    return c,l


def test_closed_reference_moments_and_positive_interval():
    identity=reference_identity();assert identity['cross_product_residual']=='0'
    assert identity['integral_residuals']==['0','0']
    c,l=inputs();r=closed_lapse_coefficient(c,l)
    assert 0.05739927758627<r['total_coefficient_interval'][0]<r['total_coefficient_interval'][1]<0.05739927758629
    assert r['strictly_positive']
    for row in r['per_group_reference']:
        if c[row['group']]['angular_eigenvalue']==0:
            assert row['coefficient_interval'][0]<=0<=row['coefficient_interval'][1]


def test_closed_local_channels_match_original_density_mixed_derivative():
    c,l=inputs();closed=closed_lapse_coefficient(c,l)
    for step in(.5,1.):
        direct=local_density_stencil(l,step)
        for name,bounds in closed['local_coefficient_intervals'].items():
            assert bounds[0]-3e-13<=direct[name]<=bounds[1]+3e-13
        assert max(abs(direct[name]) for name in closed['zero_channels'])<3e-13


def test_homogeneous_density_has_no_lapse_second_normal_derivative():
    _,l=inputs();domain=incoming_cauchy_jets()
    raw=taylor_two_jets(scalar_taylor_coefficients(domain.fields),0.,0.)[0]
    light=LightRestorationAction(magnetic_light_spectrum(4));compact=LockedSphericalLocalAction(l)
    base=local_action_densities(raw,light,compact)
    for step in(-1.,1.):
        altered=raw.copy();altered[0,3]+=step
        new=local_action_densities(altered,light,compact)
        assert max(float(abs(new[k]-base[k])) for k in base)<3e-13
