"""Recover actual stationary source fields without a horizon/scattering solve."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons.nsc_pg_batch_projection import project_mode_family
from recursive_horizons.nsc_pg_mode_recovery import stationary_interior_observation


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def stored_fields():
    record = json.loads((ROOT/'results/development/nsc-pg-massive-mode-resolution.json').read_text())
    path = ROOT/record['payload']['path']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == record['payload']['sha256']
    with np.load(path, allow_pickle=False) as payload:
        yield payload


@pytest.mark.parametrize('group', [13, 14])
def test_recover_archived_massive_source_columns(stored_fields, group):
    """Independent forward projection starts at rho=0, inverse at rho=1."""
    data = stored_fields
    label = data['label']
    candidates = np.flatnonzero((label[:, 0] == group) & (label[:, 1] == 1))
    indices = candidates[[np.argmin(abs(label[candidates, 5]-value))
                          for value in (.2, 2., 9.)]]
    energies = label[indices, 5]
    mass, angular = label[indices[0], 3:5]
    projection = project_mode_family(energies, mass, angular,
                                    data['interior'][indices, 1],
                                    data['exterior'][indices, 1])
    owner = stationary_interior_observation(energies, mass, angular)
    recovered = owner.recover(projection, data['source_projector'][indices])
    difference = np.linalg.norm(recovered.mode_at_one-data['interior'][indices, 2], axis=(1, 2))
    assert max(difference) < 3e-9
    assert max(row['projection_absolute_residual'] for row in recovered.diagnostics) < 3e-11
    assert max(row['current_coisometry_residual'] for row in recovered.diagnostics) < 3e-9
    assert max(owner.current_residual) < 3e-11
    for index, error, row in zip(indices, difference, recovered.diagnostics):
        print(json.dumps({'group': group, 'archived_index': int(index),
                          'field_error_vs_archived': float(error), **row}))


def test_unresolved_projection_and_covariance_are_rejected():
    owner = stationary_interior_observation([.5], np.pi/2, 0.)
    with pytest.raises(ValueError, match='covariance is not input'):
        owner.recover(np.eye(8)[None], np.eye(3)[None])
    # The declared upstream error alone prevents a resolved inverse; no
    # arbitrary singular-value floor, state completion, or finite stress.
    result = owner.recover(np.zeros((1, 8, 3)), np.eye(3)[None],
                           input_absolute_error=1., field_absolute_tolerance=1e-12,
                           require_all=False)
    assert not result.accepted.any()
    assert np.isnan(result.mode_at_one).all()
    assert result.diagnostics[0]['projection_absolute_residual'] is None
    assert result.diagnostics[0]['status'].startswith('OPEN')
    with pytest.raises(ValueError, match='unresolved stationary mode recovery'):
        result.require_all()


def test_complex_contour_nodes_are_not_real_physical_state_inputs():
    with pytest.raises(ValueError, match='real positive frequencies'):
        stationary_interior_observation([.5+.01j], np.pi/2, 0.)
