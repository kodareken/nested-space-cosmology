"""Changing an analytic profile cannot reuse an amplitude-only cache key."""
from dataclasses import replace
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_compatible_history_geometry import CompatibleIncomingMetric
from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction
from test_nsc_ks_difference_envelope import inputs


def test_same_amplitude_different_profile_changes_identity_even_at_zero():
    family=inputs(0.)[2]
    direction=family.directions[0]
    changed=replace(direction,w=LocalAxialFunction((0.,1.,.2),direction.w.center))
    other=CompatibleIncomingMetric(family.amplitudes,(changed,))
    assert profile_identity(family)!=profile_identity(other)


def test_axial_table_can_be_reused_across_normal_windows_but_history_cannot():
    family=inputs()[2]
    changed=replace(family.directions[0],inner_radius=.006)
    other=CompatibleIncomingMetric(family.amplitudes,(changed,))
    assert profile_identity(family)!=profile_identity(other)
    assert profile_identity(family,include_normal_window=False)==profile_identity(other,include_normal_window=False)
