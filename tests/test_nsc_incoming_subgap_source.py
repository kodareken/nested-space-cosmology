"""Replay the bounded group13 input artifact; never restart its mode solves."""
from pathlib import Path
from hashlib import sha256
import json

import numpy as np
import pytest

from recursive_horizons.nsc_incoming_subgap_source import (
    read_prepared_incoming_subgap, recompose_incoming_subgap, real_reference_remainder,
    integrate_local_subgap, MASS,
)
from recursive_horizons.nsc_subgap_history_response import AnalyticResponsePanel


@pytest.fixture(scope='module')
def prepared():
    root = Path(__file__).resolve().parents[1]
    payload = json.loads((root/'results/development/nsc-incoming-subgap-source.json').read_text())['payload']
    path = root/payload['path']
    assert sha256(path.read_bytes()).hexdigest() == payload['sha256']
    arrays = read_prepared_incoming_subgap(root, path)
    return arrays, recompose_incoming_subgap(arrays)


def test_same_surface_bases_and_local_normalization(prepared):
    _, report = prepared
    assert report['short_stationary_solves'] == 19
    assert report['reflection_samples_used'] == 360
    for key in ('cross_current', 'canonical_frame_vs_real_owner', 'real_basis_unitarity',
                'complex_cross_Gram', 'complex_bilinear_adjoint', 'normalization_control_at_zero',
                'normalization_control_reflection_modulus', 'real8_formula_vs_archived', 'real8_formula_vs_direct'):
        assert report['diagnostics'][key] < 3e-10, (key, report['diagnostics'][key])
    assert report['diagnostics']['T01_coherence_overlap'] < 3e-11
    assert report['diagnostics']['T01_diagonal_difference'] < 3e-11
    assert abs(report['physical_group13_contribution'][2]) < 3e-11
    assert np.isfinite(report['refined_minus_coarse8']).all()


def test_independent_complex_node_and_real_contour_refinements(prepared):
    arrays, report = prepared
    assert report['diagnostics']['complex_interpolation'] < 3e-10
    assert report['diagnostics']['real8_bilinear_interpolation'] < 3e-10
    assert max(np.max(v) for v in report['error_indicators'].values()) < 3e-10
    assert report['diagnostics']['smooth_imaginary'] < 3e-11
    with pytest.raises(ValueError, match='real axis'):
        real_reference_remainder(np.array([1.2+.03j]))
    panel = AnalyticResponsePanel(arrays['energies'], arrays['bilinears'], 1., MASS)
    with pytest.raises(ValueError, match='archived'):
        integrate_local_subgap(panel, .238, lambda z: 1., points=32)
