"""Clipped n=16 Newton proposals stay inside the owned damper."""
import importlib
from pathlib import Path

import numpy as np

from recursive_horizons.nsc_local_history_newton import LocalHistoryNewtonSettings
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily


def _module():
    import sys
    scripts = str(Path(__file__).resolve().parents[1]/'scripts')
    if scripts not in sys.path:
        sys.path.insert(0, scripts)
    return importlib.import_module('derive_nsc_ks_coupled_newton_n16_damped')


def test_clip_and_line_search_select_an_improving_positive_radius_scale():
    module = _module()
    family = LocalIncomingFamily(np.zeros((2, 16)))
    delta = np.zeros(32)
    delta[0] = 20.0
    gradient = np.ones((5, 2))
    jacobian = np.zeros((32, 5, 2))
    jacobian[0, :, 0] = -0.6
    jacobian[0, :, 1] = -0.6
    measured = np.array([0.6, 0.6])
    settings = LocalHistoryNewtonSettings()
    clipped, clip_norm, rows, selected = module.damped_scale_rows(
        family, delta, gradient, jacobian, measured, settings)
    assert clip_norm == settings.max_step_norm
    assert float(np.linalg.norm(clipped)) == settings.max_step_norm
    assert selected is not None
    assert selected['step_norm'] <= settings.max_step_norm
    assert selected['radius_lower_bound'] > 0
    assert selected['candidate_profile_identity'] not in module.FORBIDDEN_IDENTITIES
    assert module.beats_measured(selected['predicted_all_node_residual_maxima'], measured)
    assert rows[0]['step_scale'] == 1.0
    assert rows[-1]['step_scale'] == settings.min_line_search_scale


def test_forbidden_identity_is_never_selected(monkeypatch):
    module = _module()
    family = LocalIncomingFamily(np.zeros((2, 16)))
    delta = np.array([0.1]+[0.0]*31)
    gradient = np.full((3, 2), 0.2)
    jacobian = np.zeros((32, 3, 2))
    measured = np.array([1.0, 1.0])
    settings = LocalHistoryNewtonSettings()

    def fake_identity(_candidate):
        return module.FORBIDDEN_IDENTITIES[0]

    monkeypatch.setattr(module, 'profile_identity', fake_identity)
    _, _, rows, selected = module.damped_scale_rows(
        family, delta, gradient, jacobian, measured, settings)
    assert selected is None
    assert all(row['forbidden_identity'] for row in rows)
