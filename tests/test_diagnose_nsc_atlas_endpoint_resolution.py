"""Controls for locating the incoming support independently of cell counts."""
from fractions import Fraction
import math

import pytest

from diagnose_nsc_atlas_endpoint_resolution import support_at_cut


def test_front_and_clock_origin_match_characteristic_labels():
    # F(q)=max(q-q0,0)^4. The incoming field is nonzero iff T+z>q0;
    # the outgoing field is still zero because T-z<q0.
    central,offset,q0=.2,-1.,-.5
    row=support_at_cut([0,.1,.2,1.],radius=math.tanh(.6),
        central_time=central,clock_offset=offset,onset_label=q0)
    for z in (.05,.25,.35,.55):
        f=lambda q:max(q-q0,0.)**4
        assert bool(f(central+offset+z)-f(central+offset-z))==(z>row['incoming_front_z'])
    assert row['incoming_supported_cell_indices']==[2]
    assert row['last_cell_clipped_width']==pytest.approx(.4)
    assert row['front_to_last_inner_knot']==pytest.approx(.1)


def test_actual_subdivision_of_active_layer_is_detected():
    kwargs=dict(radius=math.tanh(.6),central_time=.2,clock_offset=-1.,onset_label=-.5)
    original=support_at_cut([0,.1,.2,1.],**kwargs)
    split=support_at_cut([0,.1,.2,.4,.5,1.],**kwargs)
    assert original['incoming_supported_cell_indices']==[2]
    assert split['incoming_supported_cell_indices']==[2,3,4]
    assert not split['all_inner_cells_before_onset']
    assert original['incoming_supported_width']==split['incoming_supported_width']


def test_quiet_cap_has_no_active_cells_and_rational_knots_work():
    row=support_at_cut(['0','1/10','1'],radius=.45,central_time=0,
        clock_offset=-1.,onset_label=-.2)
    assert row['incoming_supported_cell_indices']==[]
    assert row['incoming_supported_width']==0


@pytest.mark.parametrize('change',[
    {'central_time':.6}, {'radius':1.}, {'knots':[0,.9,1.]}, {'knots':[0,.1,.1,1.]},
])
def test_invalid_or_outgoing_regime_rejected(change):
    kwargs=dict(knots=[Fraction(0),Fraction(1,5),Fraction(1)],radius=.5,
        central_time=.2,clock_offset=-1.,onset_label=-.5)
    kwargs.update(change)
    with pytest.raises(ValueError): support_at_cut(**kwargs)
