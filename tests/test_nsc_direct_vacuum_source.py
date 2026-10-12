"""Per-row homogeneous vacuum Bloch transport against one archived source row."""
import time
from pathlib import Path

import numpy as np
import pytest
from flint import arb, arb_series, ctx

from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_massive_jost_phase_bound import _series_context
from recursive_horizons.nsc_metric_horizon_frame import horizon_q
from recursive_horizons.nsc_source_occupation_enclosure import occupation_vacuum_distance
from recursive_horizons.nsc_subgap_source_covariance import (
    BlochSource, original_covariance_distance_bounds,
)
from recursive_horizons.nsc_subgap_upstream_covariance import negative_subgap_covariance_error
from recursive_horizons.nsc_vacuum_source_correction import VacuumCorrection
from recursive_horizons.nsc_vacuum_source_remainder import VacuumSourceExpansion
from recursive_horizons import nsc_direct_vacuum_source as direct
from recursive_horizons.nsc_direct_vacuum_source import (
    DirectVacuumBloch, archived_vacuum_distance, capture_direct_vacuum,
    compare_signed_archives, validate_direct_vacuum,
)

ROOT = Path(__file__).resolve().parents[1]
START = -18.
STEP = .00625


def _moderate_selection(archive):
    entries = archive.family_entries((14, 1))
    choice = min(
        ((float(energy), batch, index)
         for batch, _channel in entries if batch.energy_sign > 0
         for index, energy in enumerate(batch.source.energies[::3]) if energy >= 2),
        key=lambda row: row[0])
    energy, positive, index = choice
    negative = next(
        batch for batch, _channel in entries
        if batch.energy_sign < 0 and batch.original_panel == positive.original_panel
        and batch.rows == positive.rows)
    if not np.all(negative.source.energies[3*index:3*index+3] == -energy):
        raise AssertionError('the signed group14 partner does not carry -E')
    return energy, positive, negative, index


@pytest.fixture(scope='module')
def moderate():
    archive = RetainedUpstreamArchive(ROOT)
    energy, positive, negative, index = _moderate_selection(archive)
    config = archive.meta['config']
    if positive.mass != np.pi/2 or positive.angular != np.sqrt(5):
        raise AssertionError('original group14 mass and angular label required')
    with ctx.workprec(192):
        mass = arb(positive.mass).union(arb.pi()/2)
        angular = arb(positive.angular).union(arb(5).sqrt())
        model = DirectVacuumBloch(energy, mass, angular, config['horizon_rho'])
    sl = slice(3*index, 3*index+3)
    started = time.process_time()
    trace, nfev = capture_direct_vacuum(
        model, START, positive.rho_up, max_step=STEP, cpu_limit=30)
    cpu = time.process_time()-started
    proof = validate_direct_vacuum(model, trace, START, positive.rho_up)
    compared = compare_signed_archives(
        model, proof, positive.initial_columns[:, sl], positive.source.covariance[sl, sl],
        negative.initial_columns[:, sl], negative.source.covariance[sl, sl],
        kappa=config['surface_gravity'], omega=config['omega'])
    return {
        'model': model, 'trace': trace, 'proof': proof, 'compared': compared,
        'energy': energy, 'positive': positive, 'negative': negative, 'index': index,
        'slice': sl, 'hint': config['horizon_rho'], 'cpu': cpu, 'nfev': nfev,
        'kappa': config['surface_gravity'], 'omega': config['omega'],
    }


def test_group14_moderate_signed_row_stays_inside_demonstrated_scope(moderate):
    model, proof, compared = moderate['model'], moderate['proof'], moderate['compared']
    initial = model.vacuum_initial(START)
    assert moderate['cpu'] < 30
    assert proof['cells'] == 2680
    assert proof['physical_local_gate'] == 'OPEN'
    assert proof['initial_error'] > 0
    assert proof['initial_error'] <= proof['bloch_error']
    assert max(component.rad() for component in initial) > 0
    assert proof['initial_error'] >= max(component.rad() for component in initial)
    assert proof['bloch_error'] > arb('1.58e-12')
    assert proof['bloch_error'] < arb('1.60e-12')
    direct_positive = original_covariance_distance_bounds(
        moderate['positive'].initial_columns[:, moderate['slice']],
        moderate['positive'].source.covariance[moderate['slice'], moderate['slice']], proof)
    assert direct_positive['upper'] == compared['archived_positive']['upper']
    assert direct_positive['lower'] == compared['archived_positive']['lower']
    occupation = occupation_vacuum_distance(
        model.energy, model.mass, moderate['kappa'], moderate['omega'], bits=model.bits)
    assert occupation['operator_distance_upper'] == compared['occupation']['operator_distance_upper']
    assert occupation['source_occupation_law_changed'] is False
    assert occupation['incoming_gap_used'] is True
    once = (direct_positive['upper']+occupation['operator_distance_upper']).upper()
    assert once > direct_positive['upper']
    assert once > arb('2.38e-12')
    assert once < arb('2.40e-12')
    manual = dict(proof)
    manual['endpoint'] = (proof['endpoint'][0], -proof['endpoint'][1], -proof['endpoint'][2])
    manual_negative = original_covariance_distance_bounds(
        moderate['negative'].initial_columns[:, moderate['slice']],
        moderate['negative'].source.covariance[moderate['slice'], moderate['slice']], manual)
    assert manual_negative['upper'] == compared['archived_negative']['upper']
    assert manual_negative['lower'] == compared['archived_negative']['lower']
    assert compared['archived_negative']['upper'] == negative_subgap_covariance_error(
        moderate['negative'].initial_columns[:, moderate['slice']],
        moderate['negative'].source.covariance[moderate['slice'], moderate['slice']], proof)
    assert proof['endpoint'] == (
        manual['endpoint'][0], -manual['endpoint'][1], -manual['endpoint'][2])
    changed = moderate['negative'].initial_columns[:, moderate['slice']].copy()
    changed[0, 1] += 1e-6
    moved = archived_vacuum_distance(
        changed, moderate['negative'].source.covariance[moderate['slice'], moderate['slice']],
        proof, complement=True)
    assert moved['upper'] != compared['archived_negative']['upper']
    assert set(compared) == {'archived_positive', 'archived_negative', 'occupation'}
    print({
        'E': moderate['energy'], 'cells': proof['cells'], 'nfev': moderate['nfev'],
        'cpu': moderate['cpu'], 'vacuum_error': float(proof['bloch_error']),
        'initial_error': float(proof['initial_error']),
        'positive_archived_upper': float(direct_positive['upper']),
        'occupation_upper': float(occupation['operator_distance_upper']),
        'positive_source_upper': float(once),
        'negative_archived_upper': float(compared['archived_negative']['upper']),
        'negative_archived_lower': float(compared['archived_negative']['lower']),
    })


def test_saved_trace_replay_does_not_solve(moderate, monkeypatch):
    def boom(*_args, **_kwargs):
        raise AssertionError('replay attempted another solve')

    monkeypatch.setattr(direct, 'solve_ivp', boom)
    monkeypatch.setattr(direct, 'capture_direct_vacuum', boom)
    monkeypatch.setattr(moderate['model'], 'rhs_numeric', boom)
    replayed = validate_direct_vacuum(
        moderate['model'], moderate['trace'], START, moderate['positive'].rho_up)
    assert replayed['bloch_error'] == moderate['proof']['bloch_error']
    assert replayed['initial_error'] == moderate['proof']['initial_error']
    assert replayed['cells'] == moderate['proof']['cells']
    assert replayed['endpoint'] == moderate['proof']['endpoint']


def test_wrong_start_broken_join_and_capture_guards(moderate, monkeypatch):
    model, trace, rho = moderate['model'], moderate['trace'], moderate['positive'].rho_up
    with pytest.raises(ValueError, match='start'):
        validate_direct_vacuum(model, trace, -18.1, rho)
    shifted = trace.copy()
    shifted[0, 0] = -17.5
    with pytest.raises(ValueError, match='start'):
        validate_direct_vacuum(model, shifted, START, rho)
    broken = trace.copy()
    broken[1, 2] += 1e-8
    with pytest.raises(ValueError, match='join'):
        validate_direct_vacuum(model, broken, START, rho)
    gap = trace.copy()
    gap[1, 0] = gap[0, 1]+1e-12
    with pytest.raises(ValueError, match='join'):
        validate_direct_vacuum(model, gap, START, rho)

    def boom(*_args, **_kwargs):
        raise AssertionError('solver ran without an explicit budget')

    monkeypatch.setattr(direct, 'solve_ivp', boom)
    for bad_budget in (None, False, True, 0, -1):
        with pytest.raises(ValueError, match='CPU budget'):
            capture_direct_vacuum(model, START, rho, max_step=STEP, cpu_limit=bad_budget)
    for bad_step in (None, False, 0):
        with pytest.raises(ValueError, match='max step'):
            capture_direct_vacuum(model, START, rho, max_step=bad_step, cpu_limit=30)
    with pytest.raises(ValueError, match='target rho'):
        capture_direct_vacuum(model, START, None, max_step=STEP, cpu_limit=30)
    with pytest.raises(TypeError):
        capture_direct_vacuum(object(), START, rho, max_step=STEP, cpu_limit=30)


def test_numeric_rhs_excludes_correction_forcing(moderate):
    model = moderate['model']
    y = -1.5
    state = np.array([.2, -.4, .3])
    delta = np.exp(y)
    metric = np.polynomial.polynomial.polyval(-delta, model.numeric_coefficients)
    radius = 1/np.sin(float(model.q)-delta)
    scale = np.sqrt(delta/metric)
    generator = np.array([
        -float(model.mass)*radius*scale,
        float(model.angular)*scale,
        -float(model.energy)/metric])
    np.testing.assert_array_equal(model.rhs_numeric(y, state), 2*np.cross(generator, state))
    with _series_context(model.bits, 0):
        series = model.hamiltonian_series(arb_series([arb(y)], prec=1), 0)
    series_generator = np.array([float(component[0].mid()) for component in series])
    np.testing.assert_allclose(
        model.rhs_numeric(y, state), 2*np.cross(series_generator, state),
        rtol=1e-12, atol=1e-14)
    expansion = VacuumSourceExpansion(
        moderate['hint'], model.mass, model.angular, order=2, bits=model.bits)
    correction = VacuumCorrection(expansion, model.energy)
    with _series_context(correction.bits, 0):
        defect = correction.defect_series(arb_series([arb(y)], prec=1), 0)
    forcing = np.array([float(component[0].mid()) for component in defect])
    assert np.linalg.norm(forcing) > 1e-8
    forced = correction.rhs_numeric(y, state)
    np.testing.assert_allclose(model.rhs_numeric(y, state)-forced, forcing, rtol=1e-12, atol=1e-12)
    assert not np.allclose(model.rhs_numeric(y, state), forced)
    assert not issubclass(DirectVacuumBloch, VacuumCorrection)
    assert not issubclass(DirectVacuumBloch, BlochSource)
    with pytest.raises(ValueError, match='subgap'):
        BlochSource(model.q, model.energy, model.mass, model.angular)


def test_initial_projector_norm_domain_and_signed_parameters(moderate):
    model = moderate['model']
    before = ctx.prec, ctx.cap
    initial = model.vacuum_initial(START)
    norm = sum((component*component for component in initial), arb(0))
    assert norm.contains(1)
    assert initial[2] > .5
    assert all(component.rad() > 0 for component in initial)
    with ctx.workprec(model.bits):
        expected_q = horizon_q(moderate['hint'], bits=model.bits)
        assert model.q.lower().man_exp() == expected_q.lower().man_exp()
        assert model.q.upper().man_exp() == expected_q.upper().man_exp()
    assert model.energy == model._frame.energy
    assert not model.mass > model.energy
    with pytest.raises(ValueError):
        model.vacuum_initial(-1)
    with pytest.raises(ValueError):
        model.vacuum_initial(None)
    with pytest.raises(ValueError):
        model.vacuum_initial(False)
    with pytest.raises(ValueError):
        model.target_log_distance(None)
    with pytest.raises(ValueError):
        model.target_log_distance(True)
    with pytest.raises(ValueError):
        model.target_log_distance(.5)
    with pytest.raises(ValueError):
        model.target_log_distance(100)
    hint = moderate['hint']
    for energy, mass, angular, horizon in (
            (None, 0, 0, hint), (True, 0, 0, hint), (0, 0, 0, hint), (-1, 0, 0, hint),
            (1, None, 0, hint), (1, False, 0, hint), (1, -.1, 0, hint),
            (1, 0, None, hint), (1, 0, True, hint), (1, 0, 0, None), (1, 0, 0, False),
            (1, 0, 0, .5)):
        with pytest.raises((ValueError, ArithmeticError)):
            DirectVacuumBloch(energy, mass, angular, horizon, bits=80, frame_order=2)
    for kwargs in ({'bits': True}, {'bits': 40}, {'frame_order': True}, {'frame_order': 1},
                   {'metric_terms': True}, {'metric_terms': 8}, {'analytic_radius': None}):
        with pytest.raises(ValueError):
            DirectVacuumBloch(.4, 0, -1, hint, **kwargs)
    outward = DirectVacuumBloch(.5, 0, -np.sqrt(5), hint, bits=80, frame_order=2)
    assert outward.energy > outward.mass
    assert outward.angular < 0
    state = np.array([0., 0., 1.])
    assert outward.rhs_numeric(-2., state)[0] < 0
    inward = DirectVacuumBloch(.5, 0, np.sqrt(5), hint, bits=80, frame_order=2)
    assert inward.rhs_numeric(-2., state)[0] > 0
    enclosed = outward.vacuum_initial(START)
    assert sum((component*component for component in enclosed), arb(0)).contains(1)
    assert (ctx.prec, ctx.cap) == before


def test_timeout_does_not_return_a_trajectory(moderate, monkeypatch):
    calls = {'n': 0}

    def clock():
        calls['n'] += 1
        return 0. if calls['n'] == 1 else 1e9

    monkeypatch.setattr(direct.time, 'process_time', clock)
    caught = {}
    try:
        trace, _nfev = capture_direct_vacuum(
            moderate['model'], START, moderate['positive'].rho_up, max_step=STEP, cpu_limit=1)
        caught['trace'] = trace
    except TimeoutError as exc:
        caught['error'] = str(exc)
    assert 'trace' not in caught
    assert 'CPU budget' in caught['error']
    assert calls['n'] >= 2


def test_missing_or_boolean_archives_are_not_zero():
    proof = {
        'endpoint': (arb(0), arb(0), arb(1)),
        'bloch_error': arb(0),
        'bits': 64,
    }
    columns = np.array([[1, 0, 0], [0, 1, 0]], complex)
    covariance = np.diag([1., 0., 0.]).astype(complex)
    positive = archived_vacuum_distance(columns, covariance, proof, complement=False)
    complemented = archived_vacuum_distance(columns, covariance, proof, complement=True)
    assert complemented['upper'] > positive['upper']+arb('.1')
    assert proof['endpoint'] == (arb(0), arb(0), arb(1))
    for columns_bad, covariance_bad in (
            (None, covariance), (columns, None), (False, covariance),
            (columns, True), (np.ones((2, 3), dtype=bool), covariance),
            (np.array([[None, 0, 0], [0, 1, 0]], dtype=object), covariance)):
        with pytest.raises(ValueError):
            archived_vacuum_distance(columns_bad, covariance_bad, proof)
    with pytest.raises(ValueError):
        archived_vacuum_distance(columns, covariance, None)
    with pytest.raises(ValueError):
        archived_vacuum_distance(columns, covariance, {'bits': 64, 'bloch_error': arb(0)})
    with pytest.raises(ValueError):
        archived_vacuum_distance(columns, covariance, proof, complement=None)
    model = DirectVacuumBloch(.4, 0, 1, 1.9006916054701435, bits=80, frame_order=2)
    with pytest.raises(ValueError):
        compare_signed_archives(
            model, {'bits': model.bits, 'endpoint': proof['endpoint'], 'bloch_error': arb(0)},
            columns, covariance, None, covariance, kappa=.2, omega=4)
    for kappa, omega in ((None, 4), (.2, False), (True, 4)):
        with pytest.raises(ValueError):
            compare_signed_archives(
                model, {'bits': model.bits, 'endpoint': proof['endpoint'], 'bloch_error': arb(0)},
                columns, covariance, columns, covariance, kappa=kappa, omega=omega)
    with pytest.raises(TypeError):
        compare_signed_archives(
            model, {'bits': model.bits}, columns, covariance, energy=.4, mass=0, kappa=.2, omega=4)


def test_complex_trace_and_missing_error_are_rejected(moderate):
    bad=moderate['trace'].astype(complex);bad[0,2]+=1j
    with pytest.raises(ValueError,match='real'):
        validate_direct_vacuum(moderate['model'],bad,START,moderate['positive'].rho_up)
    for missing in (None,False,True):
        proof=dict(moderate['proof']);proof['bloch_error']=missing
        with pytest.raises(ValueError,match='explicit'):
            archived_vacuum_distance(np.zeros((2,3),complex),np.eye(3),proof)
