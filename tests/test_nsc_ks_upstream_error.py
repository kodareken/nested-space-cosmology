"""A-posteriori upstream continuation and low/subgap decision controls."""
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import derive_nsc_ks_source_control_v2 as C
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_ks_source_inventory import RetainedSourceInventory
from recursive_horizons.nsc_ks_upstream_error import (
    classify_low_subgap_windows,
    combine_upstream_row_error,
    continue_columns,
    dop853_power_coefficients,
    eval_power_polynomial,
    independent_constant_matrix_control,
    nine_component_upstream_feed,
)
from scipy.integrate import DOP853


def family14_inputs():
    archive = RetainedUpstreamArchive(ROOT)
    source, _negative, archived, first, channel, selections = C.select_source(
        archive, (14, 1))
    unique = source.energies.reshape(-1, 3)[:, 0]
    at_one = []
    with RetainedSourceInventory(ROOT) as inventory:
        panels = {panel.name: panel for panel in inventory.positive_panels(groups=(14,))
                  if panel.angular_sign == 1}
        for energy in unique:
            panel = next(value for value in panels.values()
                         if np.any(value.energies == energy))
            index = np.flatnonzero(panel.energies == energy)
            assert len(index) == 1
            at_one.append(panel.amplitudes_at_one[index[0]])
    return source, archived, first, channel, selections, unique, np.asarray(at_one)


def test_constant_matrix_majorant_contains_independent_solution():
    result = independent_constant_matrix_control()
    assert result['majorant_contains_independent_error']
    assert result['refinement_movement_is_not_the_bound']
    assert not result['physical_source_error_included']


def test_dop853_power_coefficients_reproduce_dense_output():
    matrix = np.array([[0, 1j], [1j, 0]], complex)
    solver = DOP853(lambda _t, value: matrix @ value, 0., np.array([1.+0j, 0j]),
                    .1, rtol=1e-10, atol=1e-12, max_step=.05)
    while solver.status == 'running':
        solver.step()
        dense = solver.dense_output()
        coefficients = dop853_power_coefficients(dense)
        for fraction in (0., .2, .7, 1.):
            np.testing.assert_allclose(
                eval_power_polynomial(coefficients, fraction),
                dense(dense.t_old + fraction * dense.h), atol=2e-17, rtol=0)


def test_family14_directed_continuation_contains_archived_endpoint_indicator():
    _source, archived, first, _channel, _selections, energies, at_one = family14_inputs()
    result = continue_columns(at_one, energies, first.mass, first.angular,
                              max_step=1/4096, bits=90)
    bound = result['unweighted_frobenius_error_upper']
    numeric = int(bound['mantissa']) * 2.0 ** int(bound['exponent'])
    difference = np.linalg.norm(result['amplitudes_at_rho_up'].reshape(2, -1) - archived)
    assert numeric >= difference
    assert numeric < 1e-11
    assert result['logarithmic_norm_integral']['mantissa'] == '0'
    assert not result['physical_source_error_included']


def test_missing_physical_source_error_keeps_upstream_open():
    combined = combine_upstream_row_error(1e-13, None)
    assert combined['combined_row_error_upper'] is None
    feed = nine_component_upstream_feed([1e-13, 2e-13], None)
    assert feed['combined_upstream_bound_N_beta'] is None
    assert feed['combined_status'].startswith('OPEN')


def test_low_subgap_classification_keeps_reconstruction_unbounded():
    record = json.loads((ROOT / 'results/development/nsc-ks-source-inventory.json').read_text())
    v5 = json.loads((ROOT / 'results/development/nsc-incoming-source-update-v5.json').read_text())
    with RetainedSourceInventory(ROOT) as inventory:
        result = classify_low_subgap_windows(record, inventory.positive_panels(), v5['coverage'])
    assert result['ell0_zero_response_groups'] == [0, 13, 23]
    assert result['unbounded_reconstruction_windows'] > 0
    assert result['unbounded_reconstruction_rows'] > 0
    assert result['remaining_low_subgap_source_error_bound'] is None
    assert not result['pauli_or_nested_node_relabelled_as_bound']
