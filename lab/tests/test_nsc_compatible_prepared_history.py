"""Integrated control receipt; scientific calculation is not rerun here."""
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    'compatible_prepared_control', ROOT/'scripts/derive_nsc_compatible_prepared_history.py')
CONTROL = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CONTROL)


def test_integrated_receipt_and_provenance():
    record = json.loads(CONTROL.OUTPUT.read_text())
    assert record['source_sha256'] == CONTROL.provenance()
    residuals, settings = record['residuals'], record['control']
    assert residuals['complete_field_derivative'] < settings['field_derivative_tolerance']
    assert residuals['complete_coherent_covariance_derivative'] < settings['covariance_derivative_tolerance']
    assert min(record['omitted_term_effects'].values()) > 1e-6
    assert record['scope']['physical_parent_preparation'] == 'OPEN'
    assert record['scope']['incoming_C0_match'] == 'OPEN'
    assert record['scope']['extended_stationarity'] == 'OPEN'
    assert record['scope']['continuum_error_bound'] is None
    assert record['scope']['metric_timestep'] is False
    assert record['scope']['source_stress_computed'] is False
    assert settings['synthetic_columns_only'] is True
    assert settings['selected_physical_duration'] is None
    assert 'actual rho1' in settings['covariance_projection']
    assert settings['stationary_phase_reapplied'] is False
    assert settings['incoming_trace_interpolation'] is None
    assert residuals['incoming_trace_derivative'] < settings['field_derivative_tolerance']
