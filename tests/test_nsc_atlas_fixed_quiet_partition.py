"""A certified zero extension has no integration break at its old endpoint."""
import math

import numpy as np

from recursive_horizons.nsc_atlas_cap_background_jets import RetainedCapBackground
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis, exact_number
from recursive_horizons.nsc_atlas_cap_fixed_reference import _background_partition
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from run_nsc_atlas_fixed_endpoint import cap_knots


def test_sub_ulp_quiet_endpoint_is_not_a_physical_formula_break():
    budget=Budget(maximum_evaluations=100000, cpu_seconds=10.)
    Z=math.atanh(.45)
    qa=np.nextafter(Z,-np.inf)
    history=RetainedCapBackground([qa,qa+.1,qa+.2],[0.,0.,1.],
        [0.,0.,0.],[0.,0.,0.],onset=dict(index=1,quartic_coefficient=1.))
    assert history.profile.quiet_extension
    for level in (0,1,2):
        basis=FixedCapBasis(cap_knots(level),budget=budget)
        breaks,owners=_background_partition(basis,exact_number(Z,'Z'),0.,history)
        assert np.all(np.diff(breaks)>0)
        assert len(owners)==len(breaks)-1
        assert breaks[-1] == basis.cells-1+float(
            (exact_number(Z,'Z')-basis.exact_knots[-2])/
            (basis.exact_knots[-1]-basis.exact_knots[-2]))
        for z in (0.,qa,Z):
            jet=history.reduced_jets(0.,z,budget=budget)
            assert all(jet[k]==0. for k in ('chi','k','f','p','h'))
