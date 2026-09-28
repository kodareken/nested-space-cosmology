"""Exact new slope and distinct lapse-pivot/principal normal-data roots."""
import importlib.util
import json
from pathlib import Path

from recursive_horizons.nsc_incoming_surface_regular_branch import affine_identities,enclose_branch

ROOT=Path(__file__).resolve().parents[1]


def script():
    spec=importlib.util.spec_from_file_location('surface_branch',ROOT/'scripts/derive_nsc_incoming_surface_regular_branch.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def test_new_affine_derivative_and_independent_C_identities_are_exact():
    result=affine_identities()
    assert set(result['residuals'].values())=={'0'}
    assert not result['Cu_Cv_producers_called']
    assert result['residuals']['determinant_slope_from_C']=='0'


def test_distinct_roots_and_safe_component_are_directed():
    module=script();lapse,principal,channels,cauchy=module.inputs()
    result=enclose_branch(lapse,principal,channels,cauchy);bounds=result['intervals']
    assert bounds['A1_total'][1]<0 and bounds['R1'][1]<0
    assert 0<bounds['w_principal'][0]<bounds['w_principal'][1]<bounds['w_lapse_pivot'][0]
    assert bounds['A_at_principal_root'][0]>0
    assert bounds['determinant_at_lapse_pivot'][0]>0 and bounds['B_at_lapse_pivot'][0]>0
    assert result['regular_component_containing_zero']['guaranteed_numeric_subinterval']['upper_exclusive']==bounds['w_principal'][0]
    assert result['imported_unchanged']['A0']==lapse['coefficient']['total_coefficient_interval']
    assert result['imported_unchanged']['d']==principal['matrix']['matrix_intervals'][1][0]
    assert result['imported_unchanged']['e']==principal['matrix']['matrix_intervals'][1][1]
    assert not result['spacetime_singularity_claimed']


def test_record_uses_saved_baseline_and_never_runs_old_producers(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('old coefficient producer called')
    monkeypatch.setattr('recursive_horizons.nsc_incoming_lapse_coefficient.closed_lapse_coefficient',forbidden)
    monkeypatch.setattr('recursive_horizons.nsc_incoming_surface_principal.directed_principal_matrix',forbidden)
    module=script();result=module.make_record()
    cross=result['independent_C_crosscheck']
    assert cross['intersection'][0]<=cross['intersection'][1]
    assert not cross['C_used_to_compute_roots'] and not cross['C_D_F_producer_called']
    assert not result['scope']['C_D_F_producers_rerun']
    assert not result['scope']['physical_initial_data_selected']
    if module.OUTPUT.exists():assert result==json.loads(module.OUTPUT.read_text())
