"""No-evolution diagnosis replay and the owned interior symbol discriminator."""
import importlib.util
import json
from pathlib import Path
import numpy as np
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('spatial_diagnosis',ROOT/'scripts/derive_nsc_retarded_spatial_diagnosis.py')
D=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(D)


def test_owned_fourth_order_group_factor_is_symbol_derivative():
    theta=sp.Symbol('theta',real=True)
    symbol=(8*sp.sin(theta)-sp.sin(2*theta))/6
    assert sp.trigsimp(sp.diff(symbol,theta)-(4*sp.cos(theta)-sp.cos(2*theta))/3)==0


def test_saved_and_local_diagnosis_replays_without_any_evolution(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('propagation/producer entered a diagnosis')
    monkeypatch.setattr(D.pilot,'run',forbidden)
    monkeypatch.setattr(D.pilot,'one_case',forbidden)
    monkeypatch.setattr(D.pilot,'evolve_prepared_field_jets',forbidden)
    monkeypatch.setattr(D.pilot.FourthOrderModePropagator,'evolve',forbidden)
    record=D.make_record();saved=json.loads((ROOT/D.OUTPUT).read_text())
    assert record==saved
    assert record['scope']['field_propagations']==record['scope']['radial_solves']==0
    assert record['decision']['existing_three_solve_response_status'].startswith('OPEN')
    assert record['residuals']['common_node_forcing_identity']==0
    local=record['local_discriminator']['times']
    assert np.allclose([r['LF_relative_difference'] for r in local],[.5239659969886,.45211274657919764,.4935501286168124],rtol=0,atol=2e-12)
    assert np.allclose([r['grids']['401']['geometry_power_fraction_group_error_above_0_25'] for r in local],
                       [.3371253835272251,.3650634154225146,.3386758253347509],rtol=0,atol=2e-12)
    localization=record['localization']
    assert localization['spatial_field']['unweighted_error_squared_norm_fraction_by_spin'][0]>.999
    assert localization['spatial_field']['spin0_error_fraction_in_rho_0_95_to_1_05']>.90
    assert localization['spatial_trace']['location']['PG_time']==.253125
    assert localization['spatial_field']['unweighted_error_fraction_by_source'][2]==0
    assert not record['decision']['pulse_widening_or_source_alteration_selected']
    assert record['scope']['continuum_error_bound'] is None


def test_conditional_fast_slow_support_identities_and_saved_leakage():
    beta=sp.Symbol('beta',positive=True)
    assert sp.factor(1/(beta+1)-beta/(beta**2-1)+1/(beta**2-1))==0
    assert sp.factor(1/(beta-1)-beta/(beta**2-1)-1/(beta**2-1))==0
    saved=json.loads((ROOT/D.OUTPUT).read_text());support=saved['conditional_causal_support']
    assert np.allclose(support['PG_time_window'],[.05404795482456391,.24595204517543606],rtol=0,atol=2e-14)
    assert support['per_case']['401_64']['location']['PG_time']>support['PG_time_window'][1]
    assert not support['directed_edge_certificate']
    assert not support['discrete_finite_propagation_claimed']
