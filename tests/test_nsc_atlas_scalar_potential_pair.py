"""Scalar hot-path arithmetic agrees with the full independent array route."""
import numpy as np
import pytest

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_cap_radial_reference import reference_potential
from recursive_horizons.nsc_atlas_cap_potential_pair import reference_potential_pair

CONTROLS = ControlBudget(maximum_evaluations=1_000_000, cpu_seconds=15.)


def test_scalar_pair_agrees_with_array_route_and_parity():
    with CONTROLS.active() as budget:
        for alpha in (0., 0.00012983127142551443, .002, .5):
            fields=np.array([0., 1e-150, 1e-8, .01, .1, 1., 4., 50., 100.])
            array=reference_potential(fields,alpha,budget=budget)
            for i,f in enumerate(fields):
                p,u=reference_potential_pair(float(f),alpha,budget=budget)
                assert p == pytest.approx(array['Pdiff'][i],rel=3e-15,abs=0.)
                assert u == pytest.approx(array['Pdiff_f'][i],rel=3e-15,abs=0.)
                pm,um=reference_potential_pair(-float(f),alpha,budget=budget)
                assert p==pm and u==-um


def test_zero_coupling_known_limit_and_charge():
    with CONTROLS.active() as budget:
        before=budget.evaluations
        assert reference_potential_pair(3.,0.,budget=budget)==(9.,6.)
        assert budget.evaluations-before==12


@pytest.mark.parametrize('field,alpha',[(float('inf'),.1),(1.,-1.),(1.,float('nan')),
    (1e200,1.),(1.,np.nextafter(0.,1.))])
def test_invalid_or_unrepresentable_arithmetic_fails(field,alpha):
    with CONTROLS.active() as budget:
        with pytest.raises(ValueError):
            reference_potential_pair(field,alpha,budget=budget)
