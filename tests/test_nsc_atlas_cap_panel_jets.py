"""Panel endpoints use their own polynomial, including rounded knot jumps."""
import numpy as np
import pytest

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile
from recursive_horizons.nsc_atlas_cap_panel_jets import history_piece_indices,panel_profile_jet

CONTROLS=ControlBudget(maximum_evaluations=100000,cpu_seconds=15.)


def profile():
    result=CharacteristicProfile([-1.,0.,1.],[-1.,0.,1.],[3.,0.,3.],[-6.,0.,6.])
    # A deliberate representation mismatch at zero, with a known left limit.
    # This tests endpoint ownership, not a physical source construction.
    result._spline.c[-3,0]+=1e-5
    result._derivatives=[result._spline.derivative(k) for k in range(1,len(result._spline.c))]
    return result


def test_left_polynomial_limit_at_unchanged_knot_coordinates():
    with CONTROLS.active() as budget:
        p=profile();labels=np.array([0.,.4]);before=labels.copy()
        indices=history_piece_indices(p,[-.1,.4],budget=budget)
        ordinary=p.jet(labels,budget=budget)
        sided=panel_profile_jet(p,labels,indices,budget=budget)
        np.testing.assert_array_equal(labels,before)
        assert ordinary['F_second'][0]==0.
        assert sided['F_second'][0]==pytest.approx(2e-5,rel=0,abs=2e-15)
        for name,polynomial in zip(('F','F_prime','F_second','F_third'),
                                  (p._spline,*p._derivatives[:3])):
            assert sided[name][0]==float(np.polyval(polynomial.c[:,0],1.))
            assert sided[name][1]==ordinary[name][1]


def test_interiors_unchanged_and_resolved_crossings_rejected():
    with CONTROLS.active() as budget:
        p=profile();labels=np.array([-.3,.4]);indices=np.array([0,1])
        ordinary=p.jet(labels,budget=budget)
        sided=panel_profile_jet(p,labels,indices,budget=budget)
        for name in ordinary:np.testing.assert_array_equal(ordinary[name],sided[name])
        with pytest.raises(ValueError,match='resolved knot'):
            panel_profile_jet(p,[.1,.4],indices,budget=budget)
        with pytest.raises(ValueError,match='integer'):
            panel_profile_jet(p,labels,[0.,1.],budget=budget)
        with pytest.raises(ValueError,match='outside'):
            panel_profile_jet(p,labels,[9,1],budget=budget)


def test_quiet_left_limit_and_original_finite_domain():
    with CONTROLS.active() as budget:
        p=CharacteristicProfile([0.,1.],[0.,1.],[0.,3.],[0.,6.])
        result=panel_profile_jet(p,[0.,.3],[-1,0],budget=budget)
        assert all(values[0]==0. for values in result.values())
        with pytest.raises(ValueError,match='finite labels'):
            panel_profile_jet(p,[0.,1.01],[-1,0],budget=budget)
        with pytest.raises(ValueError,match='quiet extension'):
            panel_profile_jet(profile(),[-1.,.3],[-1,1],budget=budget)


def test_certified_zero_prefix_can_span_omitted_quiet_knots():
    with CONTROLS.active() as budget:
        p=CharacteristicProfile([0.,1.,2.,3.],[0.,0.,0.,1.],[0.,0.,0.,4.],
            [0.,0.,0.,12.],onset=dict(index=2,quartic_coefficient=1.))
        indices=history_piece_indices(p,[-.5,.5],budget=budget)
        np.testing.assert_array_equal(indices,[-1,-1])
        zero=panel_profile_jet(p,[.25,1.75],indices,budget=budget)
        assert all(np.all(row==0.) for row in zero.values())
        endpoint=panel_profile_jet(p,[2.,2.3],[-1,2],budget=budget)
        assert all(row[0]==0. for row in endpoint.values())
        with pytest.raises(ValueError,match='resolved knot'):
            panel_profile_jet(p,[2.1,2.3],[-1,2],budget=budget)
