import importlib.util
import json
from pathlib import Path

import numpy as np

from recursive_horizons.nsc_incoming_constraint_germ import structural_constraint_plane

ROOT=Path(__file__).resolve().parents[1]


def test_exact_grade_reflection_and_conditional_relation():
    certificate=structural_constraint_plane()
    assert certificate['allowed_lapse_monomials']==[[0,0],[0,2],[1,0]]
    assert certificate['allowed_shift_response_monomials']==[[0,1]]
    assert certificate['reflection']['Hamiltonian_residual']=='0'
    assert certificate['reflection']['vertex_residuals']==['0','0']
    assert certificate['conditional_relation']['substitution_residuals']==['0','0']
    assert not certificate['reflection']['state_reflection_required']


def test_independent_mixed_response_is_resolved_without_root_claim():
    r=json.loads((ROOT/'results/development/nsc-incoming-constraint-germ.json').read_text())
    np.testing.assert_allclose(r['mixed_control_response'],r['mixed_control_prediction'],rtol=0,atol=3e-11)
    assert abs(r['numerical_Jacobian_determinant'])>.9
    assert r['residuals']['first_reference_order_unchanged']==0.
    assert not r['scope']['source_accuracy_certified']
    assert not r['scope']['numerical_root_selected']
    assert not r['scope']['Cauchy_surface_constraints_solved']
    assert r['scope']['extended_stationarity']=='OPEN'


def test_replay_contracts_saved_nodes_only(monkeypatch):
    r=json.loads((ROOT/'results/development/nsc-incoming-constraint-germ.json').read_text())
    spec=importlib.util.spec_from_file_location('germ_replay',ROOT/'scripts/derive_nsc_incoming_constraint_germ.py')
    script=importlib.util.module_from_spec(spec);spec.loader.exec_module(script)
    def forbidden(*args,**kwargs):raise AssertionError('replay cannot generate reference or local response')
    monkeypatch.setattr('recursive_horizons.nsc_incoming_reference_response.incoming_reference_response',forbidden)
    monkeypatch.setattr('recursive_horizons.nsc_incoming_local_constraints.incoming_local_constraints',forbidden)
    assert script.replay(ROOT/r['payload']['path'])==r
