"""Read archived middle-panel observations; no mode/packet ODE reruns."""
import copy
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons.nsc_pg_archived_high_energy_modes import archived_middle_boundary_modes
from recursive_horizons.nsc_pg_fast_packets import _normalized_boundary
from recursive_horizons.nsc_pg_high_energy import branch_riccati_residual


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def artifacts():
    records = [json.loads((ROOT/('results/development/'+name+'.json')).read_text())
               for name in ('nsc-pg-retained-covariance', 'nsc-pg-spectral-mode-recovery')]
    files = []
    for record in records:
        path = ROOT/record['payload']['path']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record['payload']['sha256']
        files.append(np.load(path, allow_pickle=False))
    yield files
    for file in files:
        file.close()


@pytest.mark.parametrize('panel', ['mid/12_1', 'mid/32_-1'])
def test_archived_fields_cover_rejected_rows_and_accepted_overlap(artifacts, panel):
    retained, recovery = artifacts
    metadata = json.loads(retained[panel+'/metadata_json'].tobytes())
    energies = recovery[panel+'/energies']
    mask = recovery[panel+'/accepted']
    rejected = np.flatnonzero(~mask)
    accepted = np.flatnonzero(mask)
    indices = np.r_[accepted[-1], rejected[[0, len(rejected)//2, -1]]]
    result = archived_middle_boundary_modes(energies[indices], metadata)
    comparison = result.compare_observations(recovery[panel+'/observation'][indices],
                                             recovery[panel+'/projection'][indices])
    overlap = np.linalg.norm(result.mode_at_one[0]-recovery[panel+'/mode_at_one'][indices[0]])
    assert max(comparison['projection_absolute_residual']) < 3e-11
    assert overlap < 3e-8
    assert max(result.current_coisometry_residual) < 3e-11
    assert np.max(result.riccati_equation_residual) < 3e-11
    assert result.provenance['order'] == 16
    assert result.provenance['boundary_phase'].endswith('no added clock factor')
    mass = metadata['channel']['compact_mass']
    angular = metadata['angular_sign']*metadata['channel']['angular_eigenvalue']
    for branch, column in ((1, 1), (-1, 2)):
        expected = _normalized_boundary(1., energies[indices[-1]], mass, angular, branch, 16)
        assert np.linalg.norm(result.mode_at_one[-1, :, column]-expected) < 3e-14
        existing = branch_riccati_residual([1.], energies[indices[-1]], mass, angular, branch, 16)[0]
        assert abs(existing-result.riccati_equation_residual[-1, 0 if branch == 1 else 1]) < 3e-13
    print(json.dumps({'panel': panel, 'energies': result.energies.tolist(),
                      'rejected_rows_checked': 3,
                      'projection_absolute_residual': float(max(comparison['projection_absolute_residual'])),
                      'accepted_inverse_overlap': float(overlap),
                      'current_coisometry_residual': float(max(result.current_coisometry_residual)),
                      'riccati_equation_residual': float(np.max(result.riccati_equation_residual)),
                      'order12_16_difference': float(max(result.order12_16_field_difference))}))


def test_changed_producer_hash_or_out_of_panel_frequency_is_rejected(artifacts):
    retained, _ = artifacts
    metadata = json.loads(retained['mid/12_1/metadata_json'].tobytes())
    wrong = copy.deepcopy(metadata)
    wrong['sources']['src/recursive_horizons/nsc_pg_fast_packets.py'] = '0'*64
    with pytest.raises(ValueError, match='dependency changed'):
        archived_middle_boundary_modes([250.], wrong)
    with pytest.raises(ValueError, match='declared middle panel'):
        archived_middle_boundary_modes([5.], metadata)
    with pytest.raises(ValueError, match='positive real middle-panel frequencies'):
        archived_middle_boundary_modes([250.+.1j], metadata)


def test_covariance_is_not_a_field_observation(artifacts):
    retained, _ = artifacts
    metadata = json.loads(retained['mid/12_1/metadata_json'].tobytes())
    result = archived_middle_boundary_modes([250.], metadata)
    with pytest.raises(ValueError, match='source columns required'):
        result.compare_observations(np.zeros((1, 6, 2)), np.eye(8)[None])
