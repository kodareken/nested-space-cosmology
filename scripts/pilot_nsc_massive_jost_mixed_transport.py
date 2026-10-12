#!/usr/bin/env python3
"""One original low-row signed phase transport, then the archived covariance.

The reference selection is /tmp/nsc-full-low-row-pilot.py: group14/low16_1
row 15, both energy signs. This pilot does not write a record, edit an index,
or replace original source columns. One row is not a campaign.
"""
import json
from pathlib import Path
import time

import numpy as np
from flint import arb, ctx

from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_massive_jost_mixed_transport import (
    subgap_mixed_phase_transport_bound)
from recursive_horizons.nsc_massive_jost_modes import solve_jost
from recursive_horizons.nsc_massive_jost_transport_bound import (
    original_dense_solution, phase_transport_segments)
from recursive_horizons.nsc_metric_horizon_frame import (
    metric_horizon_frame, reflection_from_phase)
from recursive_horizons.nsc_paired_horizon_preparation import PairedHorizonSeedMap
from recursive_horizons.nsc_subgap_source_covariance import (
    BlochSource, capture_bloch, initial_bloch, original_covariance_error,
    validate_bloch)
from recursive_horizons.nsc_subgap_upstream_covariance import (
    negative_subgap_covariance_error)


ROOT = Path(__file__).resolve().parents[1]
ENERGY = 0.24867511687395624
ROW = 15
# A covariance upper bound at or above this scale does not separate the
# archived matrix. Crossing it marks the method too loose. It is not a
# certificate that would close the local gate.
LOOSE_COVARIANCE = 1e-3


def _float_upper(value):
    return float(value.upper() if hasattr(value, 'upper') else value)


def _emit(payload):
    print(json.dumps(payload, sort_keys=True), flush=True)


def run():
    started, wall = time.process_time(), time.perf_counter()
    archive = RetainedUpstreamArchive(ROOT)
    config = archive.meta['config']
    preparation = PairedHorizonSeedMap(
        config['horizon_rho'], config['surface_gravity'], config['omega'],
        config['horizon_offset'], config['scattering_tolerance'], config['outer_floor'])
    batches = [
        batch for batch, _ in archive.family_entries((14, 1))
        if batch.original_panel == 'group14/low16_1' and batch.rows == (0, 16)]
    positive = next(batch for batch in batches if batch.energy_sign > 0)
    negative = next(batch for batch in batches if batch.energy_sign < 0)
    energy = float(positive.source.energies[3 * ROW])
    if energy != ENERGY or positive.mass != np.pi / 2 or positive.angular != np.sqrt(5.):
        raise ValueError('original group14/low16_1 row 15 channel changed')
    solver = {'order': 8, 'radial_collar': 1e-10, 'rtol': 2e-13, 'atol': 2e-15}
    mode = solve_jost(
        preparation.background, energy, positive.mass, positive.angular, **solver)
    near = preparation.horizon_rho + 1.01e-4
    segments = phase_transport_segments(
        original_dense_solution(mode.run), preparation.horizon_rho, near)
    _emit({
        'phase': 'captured', 'energy_hex': energy.hex(), 'cells': len(segments),
        'outer_radius': mode.outer_radius, 'cpu_seconds': time.process_time() - started,
    })
    with ctx.workprec(192):
        mass = arb(positive.mass).union(arb.pi() / 2)
        angular = arb(positive.angular).union(arb(5).sqrt())
        bound = subgap_mixed_phase_transport_bound(
            energy, mass, angular, preparation.background, mode=mode,
            inner_radius=near, degree=16, metric_terms=48, defect_subdivisions=4)
    if float(segments[-1].y_end).hex() != bound['inner_y_node_hex']:
        raise ArithmeticError('transported node does not match the captured endpoint')
    phase_error = restored_upper(bound['phase_error_inner_upper'])
    _emit({
        'phase': 'transported', 'cells': bound['cells'], 'rate_signs': bound['rate_signs'],
        'phase_error_upper': _float_upper(phase_error),
        'cpu_seconds': time.process_time() - started,
    })
    with ctx.workprec(192):
        rho = arb(preparation.horizon_rho) + arb(float(segments[-1].y_end)).exp()
        theta = arb(float(segments[-1].theta_end)) + arb(0, phase_error)
        frame = metric_horizon_frame(energy, mass, angular, preparation.horizon_rho, bits=192)
        reflection, distance, tail = reflection_from_phase(frame, rho, theta)
        y_start = -18.
        initial = initial_bloch(frame, reflection, config['surface_gravity'], y_start)
        model = BlochSource(frame.q, energy, mass, angular, bits=192)
        target = (frame.q - arb.pi() / 2 - arb(positive.rho_up).atan()).log()
    trace, evaluations = capture_bloch(
        model, initial, y_start, float(target.mid()), max_step=.025)
    validated = validate_bloch(model, initial, trace, target)
    columns = slice(3 * ROW, 3 * ROW + 3)
    positive_error = original_covariance_error(
        positive.initial_columns[:, columns], positive.source.covariance[columns, columns],
        validated)
    negative_error = negative_subgap_covariance_error(
        negative.initial_columns[:, columns], negative.source.covariance[columns, columns],
        validated)
    phase_float = _float_upper(phase_error)
    positive_float = _float_upper(positive_error)
    negative_float = _float_upper(negative_error)
    loose = not (
        phase_float < 1e-4 and positive_float < LOOSE_COVARIANCE
        and negative_float < LOOSE_COVARIANCE)
    return {
        'status': 'TOO_LOOSE' if loose else 'BOUNDED_ROW_OPEN_GATE',
        'panel': 'group14/low16_1', 'row': ROW, 'energy_hex': energy.hex(),
        'mass_hex': float(positive.mass).hex(), 'angular_hex': float(positive.angular).hex(),
        'cells': bound['cells'], 'rate_signs': bound['rate_signs'],
        'min_transport_rate_lower': bound['min_transport_rate_lower'],
        'phase_error_upper': phase_float,
        'phase_error_dyadic': bound['phase_error_inner_upper'],
        'max_cell_defect_upper': bound['max_cell_defect_upper'],
        'initializer_phase_error_upper': bound['initializer_phase_error_upper'],
        'inner_offset_certified_lower': bound['inner_offset_certified_lower'],
        'reflection_tail_upper': _float_upper(tail),
        'compact_distance_upper': _float_upper(distance),
        'bloch_cells': validated['cells'], 'bloch_nfev': evaluations,
        'bloch_error_upper': _float_upper(validated['bloch_error']),
        'positive_covariance_error_upper': positive_float,
        'negative_covariance_error_upper': negative_float,
        'method_too_loose': loose,
        'loose_covariance_scale': LOOSE_COVARIANCE,
        'source_columns_replaced': False,
        'physical_rho1_source_error': None,
        'physical_local_gate': 'OPEN',
        'cpu_seconds': time.process_time() - started,
        'wall_seconds': time.perf_counter() - wall,
    }


if __name__ == '__main__':
    started, wall = time.process_time(), time.perf_counter()
    try:
        _emit(run())
    except (ArithmeticError, ValueError) as exc:
        _emit({
            'status': 'FAILED', 'error': str(exc), 'physical_local_gate': 'OPEN',
            'source_columns_replaced': False, 'row': ROW,
            'cpu_seconds': time.process_time() - started,
            'wall_seconds': time.perf_counter() - wall,
        })
        raise
