from copy import deepcopy
import math

import pytest

from recursive_horizons.nsc_atlas_connection_slice import connection_slice_record
from recursive_horizons.nsc_atlas_connection_observer import (
    causal_window_record, protected_window, require_compatible_slice,
    require_protected_clock, tortoise,
)


def test_null_arrival_and_local_lapse():
    w = protected_window(.4, 1.6, 4.)
    assert w["coordinate_duration"] == pytest.approx(.4*(6+2*math.log(4)))
    assert w["proper_duration"] == pytest.approx(w["coordinate_duration"]/math.sqrt(2))
    assert w["coordinate_duration"] == pytest.approx(tortoise(4., .4)-tortoise(1.6, .4))


@pytest.mark.parametrize("factor", [.01, 7., 1e4])
def test_change_of_length_unit_does_not_select_physics(factor):
    a = protected_window(.4, 1.6, 4.)
    b = protected_window(.4*factor, 1.6*factor, 4.*factor)
    for key in ("coordinate_duration", "proper_duration"):
        assert b[key] == pytest.approx(a[key]*factor)
    assert a["observer_lapse"] == b["observer_lapse"]


@pytest.mark.parametrize("args", [(.4, .8, 4.), (.4, .7, 4.), (.4, 4., 4.),
                                   (.4, 5., 4.), (0, 1., 2.), (.4, 1., float("inf"))])
def test_incompatible_radial_domains_rejected(args):
    with pytest.raises(ValueError):
        protected_window(*args)


def test_open_clock_exit_is_not_a_sibling_signal():
    args = (.4, 1.6, 4.)
    stop = protected_window(*args)["proper_duration"]
    assert require_protected_clock(stop/2, *args) > 0
    for t in (-1., stop, stop+1., float("nan")):
        with pytest.raises(ValueError, match="outside"):
            require_protected_clock(t, *args)


def test_larger_retained_region_increases_window():
    assert protected_window(.4, 1.6, 5.)["proper_duration"] > protected_window(.4, 1.6, 4.)["proper_duration"]


@pytest.mark.parametrize("section,key,value", [("parent_end", "ADM_mass_over_ell", .5),
    ("parent_end", "exact_vacuum_overlap_available", False), ("bridge", "K_is_zero", False),
    ("bridge", "throat_R_over_ell", .9)])
def test_incompatible_parent_slice_rejected(section, key, value):
    record = deepcopy(connection_slice_record())
    record[section][key] = value
    with pytest.raises(ValueError, match="supplied slice"):
        require_compatible_slice(record)


def test_record_quadrature_and_coordinate_checks():
    record = causal_window_record(connection_slice_record())
    for key, value in record["independent_checks"].items():
        if key.endswith("error") or key.endswith("over_ell"):
            assert value < 1e-12
    assert record["window"]["upper_endpoint_included"] is False
    assert record["consequence"]["new_physical_mechanism"] is False


@pytest.mark.parametrize("mass,cut", [(float("nan"), 4.), (.4, float("inf")), (.4, float("nan")), (-1., 4.)])
def test_nonfinite_or_nonpositive_compatibility_inputs(mass, cut):
    with pytest.raises(ValueError, match="selected mass"):
        require_compatible_slice(connection_slice_record(), mass, cut)
