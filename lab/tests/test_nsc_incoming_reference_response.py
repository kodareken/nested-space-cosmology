import json
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons.nsc_incoming_cauchy_jets import incoming_cauchy_jets, IncomingNormalJetChange
from recursive_horizons.nsc_incoming_reference_response import (
    signed_inventory, reference_vertex_change, incoming_reference_response,
)
from recursive_horizons.nsc_incoming_joint_constraints import baseline_matter_source, source_action_gradient

ROOT = Path(__file__).resolve().parents[1]


def record(name):
    return json.loads((ROOT/'results/development'/f'{name}.json').read_text())


def test_exact_zero_change_preserves_every_formal_order():
    result = reference_vertex_change(incoming_cauchy_jets(), [-4., -.2, .2, 4.], np.pi/2, 2.)
    np.testing.assert_array_equal(result['orders'], 0.)
    assert result['leading_projector_change'] == 0.


def test_multiplicity_and_lll_allocation_are_not_doubled():
    channels = record('nsc-mode-resolved-cauchy-state')['channels']
    rows = signed_inventory(channels)
    assert len(rows) == 63
    for g, c in enumerate(channels):
        assert sum(r['multiplicity'] for r in rows if r['group'] == g) == c['copy_count']*c['degeneracy']
    domain = incoming_cauchy_jets((IncomingNormalJetChange('r', 1, 1, .01),))
    lll = incoming_reference_response(domain, channels[:1], nodes=8)
    np.testing.assert_array_equal(lll['reference_action_gradient_change'], 0.)
    assert not lll['scope']['LLL_geometric_allocation_included']


def test_mixed_normal_current_changes_sign_under_spatial_reflection():
    channel = record('nsc-mode-resolved-cauchy-state')['channels'][13:14]
    values = []
    for sign in (1, -1):
        domain = incoming_cauchy_jets((IncomingNormalJetChange('r', 1, 1, sign*.01),))
        result = incoming_reference_response(domain, channel, nodes=16)
        values.append(result['reference_action_gradient_change'])
        assert result['maxima']['formal_projector_residual'] < 3e-12
    # The zero-angular channel has no r dependence in the canonical operator.
    np.testing.assert_allclose(values, 0., rtol=0, atol=3e-13)
    channel = record('nsc-mode-resolved-cauchy-state')['channels'][14:15]
    values = []
    for sign in (1, -1):
        domain = incoming_cauchy_jets((IncomingNormalJetChange('r', 1, 1, sign*.01),))
        values.append(incoming_reference_response(domain, channel, nodes=16)['reference_action_gradient_change'])
    assert abs(values[0][1]) > 1e-5
    np.testing.assert_allclose(values[0]*[1, -1, 1, 1], values[1], rtol=0, atol=3e-12)


def test_joint_allocation_equals_full_source_plus_compact_action_at_baseline():
    finite = record('nsc-incoming-spectral-source')
    subgap = record('nsc-incoming-subgap-source')
    tail = record('nsc-incoming-source-tail')
    local = record('nsc-incoming-local-constraints')['baseline']
    matter = baseline_matter_source(finite, subgap, tail)
    domain = incoming_cauchy_jets()
    joined = source_action_gradient(domain, matter['stress_approximant'])[:2]+local['local_action_gradient']
    full = np.array(finite['finite_quantum_source_approximant'])+matter['parts']['group13_replacement']+matter['parts']['paired_approximate_tail']
    other_route = source_action_gradient(domain, full)[:2]+local['compact_action_gradient']
    np.testing.assert_allclose(joined, other_route, rtol=0, atol=3e-12)
    assert matter['full_source_error_bound'] is None
    assert matter['parts']['group13_replacement'][2] == 0.
    assert matter['group13_raw_current_zero_residual'] != 0.


def test_reject_changed_intrinsic_or_ungapped_label():
    domain = incoming_cauchy_jets()
    with pytest.raises(ValueError, match='non-LLL'):
        reference_vertex_change(domain, [1.], 0., 0.)
    domain.fields[2].data[0] += .1*np.eye(2)
    with pytest.raises(ValueError, match='intrinsic'):
        reference_vertex_change(domain, [1.], np.pi/2, 0.)
