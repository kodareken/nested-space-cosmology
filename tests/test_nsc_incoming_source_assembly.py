"""Panel coverage and independent current normalization of the new observable."""
from pathlib import Path
import json

import numpy as np
import pytest

from recursive_horizons.nsc_incoming_source_assembly import assemble, read_panels, selected_pieces

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def assembled():
    return assemble(ROOT)


def test_each_signed_frequency_interval_is_counted_once(assembled):
    assert len(assembled['groups']) == 32
    assert assembled['residuals']['spectral_interval_measure'] < 3e-11
    for group in assembled['groups'].values():
        for row in group['signed_inventory'].values():
            names = [p['panel'] for p in row['pieces']]
            assert len(names) == len(set(names))
            assert abs(row['integrated_dE_measure']-row['upper_energy_endpoint']) < 3e-11
            for family in ('low', 'mid'):
                assert not (any(n.startswith(family+'/') for n in names)
                            and any(n.startswith(family+'_ref/') for n in names))
    data, meta = read_panels(ROOT)
    pieces = selected_pieces(data, meta['panels'], 14, 1)
    old_mask = dict(pieces)['group14/low16_1']
    E = data['group14/low16_1/energies']
    assert not np.any(old_mask & (((E > .25) & (E < 1)) | ((E > np.pi/2) & (E < 8))))
    assert any(n == 'group14/low32_1' for n, _ in pieces)


def test_current_matches_imported_power_normalization(assembled):
    certificate = json.loads((ROOT/'results/development/nsc-incoming-constraint-gate.json').read_text())
    power = certificate['source_bound']['LLL_Killing_power']
    geometry = certificate['source_bound']['incoming_geometry']
    expected = -power/(4*np.pi*geometry['r_squared']*geometry['a_squared'])
    assert abs(assembled['light_allocation']['LLL_state'][2]-expected) < 3e-18
    # In this already-owned preparation n_in=0 below the compact threshold;
    # above it both thermal occupations are tiny. No numerical current from
    # subtracting unit traces may swamp that physical compact current.
    compact_current = sum(assembled['groups'][str(j)]['finite_moments'][2] for j in range(13, 33))
    assert 0 <= compact_current < 1e-16
    angular_current = sum(assembled['groups'][str(j)]['finite_moments'][2] for j in range(1, 13))
    assert angular_current > 0
    assert abs(assembled['net_quadrature_change'][2]) < 3e-17
    assert assembled['light_allocation']['LLL_count'] == 1
    assert not assembled['light_allocation']['compact_complement_included']


def test_finite_source_does_not_claim_unevaluated_tail(assembled):
    assert all(not row['full_source_converged'] for row in assembled['groups'].values())
    assert all(not p['whole_spectrum'] for p in assembled['payload_metadata']['panels'].values())
    assert not assembled['payload_metadata']['history_selected']
    assert not assembled['payload_metadata']['state_changed']
    assert np.isfinite(assembled['finite_quantum_source_approximant']).all()
    # Total CAR applies to the source C, not the formal subtraction P_ad4.
    assert assembled['low_modal_diagnostics']['CAR_lower'] < 3e-11
    assert assembled['low_modal_diagnostics']['CAR_upper'] < 3e-11
