"""Strict directed predicate, exact sign composition and preserved scope."""
import importlib.util
import json
from pathlib import Path
import sys
from fractions import Fraction

import pytest

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
SPEC=importlib.util.spec_from_file_location('fixed_c0_certificate',ROOT/'scripts/derive_nsc_retarded_fixed_c0_certificate.py')
M=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(M)


def inputs():
    v=json.loads((ROOT/M.INPUTS[0]).read_text());p=json.loads((ROOT/M.INPUTS[4]).read_text())
    return v,p['control']['support']


def test_binary_lower_endpoint_and_exact_coherent_sign():
    v,p=inputs()
    r=M.combine(v['result']['final_Re_P01']['binary_interval'],v['energy'],v['channel'],v['source_config'],p)
    assert r['strict_projector_predicate']
    q=r['full_coherent_response_00_strict_upper']
    assert Fraction(q['numerator'],q['denominator'])==-Fraction(69,875000000)
    assert r['exact_constants']['full_coherent_remainder_allowance_over_I']=={'numerator':1,'denominator':500}
    assert Fraction(r['response_upper_outward_float'])>=-Fraction(69,875000000)
    # A lower bound equal to1/4 is insufficient, even with a high upper bound.
    inconclusive=[[0,1,-2,1],[0,1,-1,1]]
    r=M.combine(inconclusive,v['energy'],v['channel'],v['source_config'],p)
    assert not r['strict_projector_predicate'] and r['full_coherent_response_00_strict_upper'] is None


def test_changed_source_or_pulse_cannot_inherit_the_certificate():
    v,p=inputs();interval=v['result']['final_Re_P01']['binary_interval']
    with pytest.raises(ValueError,match='guard'):
        M.combine(interval,v['energy'],v['channel'],{**v['source_config'],'surface_gravity':.3},p)
    with pytest.raises(ValueError,match='guard'):
        M.combine(interval,v['energy'],v['channel'],v['source_config'],{**p,'normal_outer':.04})


def test_receipt_requires_physical_bound_and_retains_full_scope():
    r=json.loads(M.OUTPUT.read_text())
    assert r['source_and_input_sha256']==M.signature()
    assert r['certificate']['strict_projector_predicate']
    assert r['physical_scope']['projector_is_a_component_not_a_replacement_state']
    assert r['physical_scope']['linearized_fixed_C0_matching']=='EXCLUDED for this specified tangent'
    assert r['physical_scope']['extended_stationarity']=='OPEN'
    assert r['physical_scope']['other_parent_extensions']=='OPEN'
    assert r['physical_scope']['new_field_or_radial_solves']==0
    assert not r['physical_scope']['metric_timestep']
    assert not r['physical_scope']['source_or_coupling_change']
