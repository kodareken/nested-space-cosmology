"""The final endpoint report must preserve all accumulated numeric diagnostics."""
import copy

import pytest

from run_nsc_atlas_fixed_endpoint import TransportComparison, transport_summary
from repair_nsc_atlas_live_summary import wrap_transport


def test_summary_keeps_the_flag_once_and_all_numeric_values():
    transport=TransportComparison()
    transport.start(.1, .2)
    transport.advance(.03, -.04, .14, .18)
    summary=transport_summary(transport)
    assert summary['state_corrected'] is False
    assert all(summary[key] == value for key,value in transport.last.items())
    assert summary['max_abs_running_mass_gap'] == transport.max_abs_mass_gap


def test_live_serialization_wrapper_does_not_change_numeric_diagnostics():
    expected=TransportComparison(); expected.start(.1, .2)
    expected.advance(.03, -.04, .14, .18)
    patched=copy.deepcopy(expected)
    report={}
    wrap_transport(patched, report, {'controller_commit':'test'})
    with pytest.raises(ValueError, match='already'):
        wrap_transport(patched, report, {})
    for increments in [(0., 0., .13, .18), (.04, -.03, .16, .12)]:
        expected.advance(*increments); patched.advance(*increments)
        assert dict(state_corrected=False, **patched.last) == expected.last
        assert patched.max_abs_mass_gap == expected.max_abs_mass_gap
        assert patched.max_abs_scalar_gap == expected.max_abs_scalar_gap
