"""Integration of the real difference-state type with the KS assembly."""
import numpy as np

from recursive_horizons.nsc_ks_local_constraints import assemble_local_incoming, _REFERENCE_CACHE
from recursive_horizons.nsc_ks_matter_difference import source_fixed_matter_difference
from recursive_horizons.nsc_local_incoming_constraints import local_error_budget
from test_nsc_ks_difference_envelope import evolve, inputs
from test_nsc_ks_local_constraints import coefficients_for, channel


def assemble(amplitude):
    prepared = evolve(amplitude)
    return assemble_local_incoming({'14_1': prepared}, inputs(amplitude)[2], np.array([-.7, -.05]),
        coefficients_for(), {'14_1': channel(14, prepared.mass, prepared.angular, 1)},
        (float(prepared.z[0]), float(prepared.z[-1])), local_error_budget())


def test_stable_difference_is_used_with_the_full_constraint_jacobian():
    h = 1e-4
    base, plus, minus = [assemble(a) for a in (.002, .002 + h, .002 - h)]
    assert base['family_records']['14_1']['stable_difference']
    np.testing.assert_allclose((plus['action_gradient'] - minus['action_gradient']) / (2 * h),
                               base['history_jacobian'][0], atol=3e-8, rtol=0)
    direct = base['family_current_matter']['14_1'] - base['family_reference_matter']['14_1']
    np.testing.assert_allclose(base['family_corrections']['14_1'], direct, rtol=0, atol=1e-14)
    assert base['physical_constraint_status'] == 'OPEN'


def test_cached_reference_difference_is_retained_instead_of_removed():
    state = evolve()
    coeff = coefficients_for()
    key = (state.fixed_preparation_digest, state.binding.rtol, state.binding.atol, state.binding.max_step)
    altered = state.reference_amplitudes + 1e-7
    previous = _REFERENCE_CACHE.get(key)
    _REFERENCE_CACHE[key] = altered
    try:
        parameters = dict(axial_scale=coeff['a'], radius=coeff['r'], multiplicity=12.)
        joint = source_fixed_matter_difference(state, **parameters)
        cached = source_fixed_matter_difference(state, reference='cached', **parameters)
        assert cached['reference_integration_difference_retained'] > 9e-8
        assert np.max(abs(cached['action_gradient_change'] - joint['action_gradient_change'])) > 1e-9
    finally:
        if previous is None:
            _REFERENCE_CACHE.pop(key)
        else:
            _REFERENCE_CACHE[key] = previous


def test_authenticated_nonzero_history_arrays_reach_both_constraints():
    """Replay saved actual source fields; only the short homogeneous reference ODE runs."""
    import json
    from pathlib import Path
    import sys
    from recursive_horizons import nsc_ks_source_envelope as E
    from recursive_horizons.nsc_ks_difference_envelope import KSDifferenceIncoming

    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root / 'scripts'))
    import derive_nsc_ks_difference_control as control
    saved = control.check()
    with np.load(root / saved['payload']['path'], allow_pickle=False) as f:
        arrays = {k: f[k] for k in f.files}
    meta = json.loads(arrays['metadata_json'].tobytes())
    name, points, options = control.CASES[-1]
    family = control.R.C.P.family(control.R.C.P.ALPHA)
    mesh = E.computational_z_grid(points)
    w, U = E._sample_axial_profiles(family.directions, mesh)
    binding = E.KSEnvelopeBinding(mesh, w, U, family.amplitudes,
        tuple(d.inner_radius for d in family.directions), tuple(d.outer_radius for d in family.directions),
        E.usual_axial_support(), meta['rho_up'], 1., **options)
    assert binding.fingerprint == meta['cases'][name]['binding']
    state = KSDifferenceIncoming(arrays['target_z'], arrays[name + '/F'], arrays[name + '/dF'],
        arrays[name + '/Fz'], arrays[name + '/dFz'], arrays['source/covariance'],
        arrays['source/column_weights'], arrays['source/energies'], meta['parameters']['mass'],
        meta['parameters']['angular'], arrays['initial_canonical_columns'], meta['rho_up'], binding,
        meta['cases'][name]['preparation'], None, {}, arrays[name + '/A'], arrays[name + '/D'], arrays[name + '/Dz'])
    read = lambda path: json.loads((root / path).read_text())
    baseline = read('results/development/nsc-incoming-source-update-v5.json')
    coeff = control.R.C.P.coefficients_from_records(
        read('results/development/nsc-incoming-surface-coefficients.json'),
        read('results/development/nsc-incoming-surface-regular-branch.json'))
    ch = read('results/development/nsc-mode-resolved-cauchy-state.json')['channels'][14]
    budget = local_error_budget(baseline_covered_regions=
        baseline['partial_error_budget']['covered_spectral_regions_action_error_upper'])
    result = assemble_local_incoming({'14_1': state}, family,
        baseline['baseline']['action_gradient_approximant'], coeff,
        {'14_1': {**ch, 'angular_sign': 1}}, E.physical_incoming_interval(), budget)
    difference = result['family_current_matter']['14_1'] - result['family_reference_matter']['14_1']
    parity = float(np.max(abs(result['family_corrections']['14_1'] - difference)))
    assert parity < 3e-13
    assert result['history_jacobian'].shape == (1, 21, 2)
    assert result['scope']['physical_constraint_status'] == 'OPEN'
    assert result['family_records']['14_1']['stable_difference']
    assert result['error_budget']['changed_history_tail'] is None
    print({'assembly_parity': parity, 'partial_constraint_max': np.max(abs(result['action_gradient']), axis=0).tolist(),
           'full_retarded_matter_tangent_max': np.max(abs(result['family_matter_tangents']['14_1']), axis=(0, 1)).tolist()})
