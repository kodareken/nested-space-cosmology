"""Inverse-coordinate identities and the exact hypotheses of the imported theorem."""
import json
from pathlib import Path
import pytest
import recursive_horizons.nsc_incoming_surface_crossing as owner
ROOT=Path(__file__).resolve().parents[1]


def test_centered_system_has_no_forbidden_positive_integer_eigenvalue():
    proof=owner.crossing_identities()
    assert set(proof['residuals'].values())=={'0'}
    assert proof['Jacobian']==[['0','0'],['-R_z_at_crossing*t0**3/Delta1','-2']]
    assert proof['eigenvalues']==[0,-2]
    assert not proof['positive_integer_eigenvalues']
    assert owner.THEOREM['zero_eigenvalue_permitted']
    assert owner.THEOREM['locator']=='Proposition 2.1, page 2'


def test_crossing_receipt_reuses_frozen_coefficients_without_producers(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('scientific coefficient producer entered crossing proof')
    import recursive_horizons.nsc_incoming_surface_coefficients as coefficient
    import recursive_horizons.nsc_incoming_surface_regular_branch as branch
    for name in ('polynomial_proof','directed_coefficients','new_algebra_controls'):
        monkeypatch.setattr(coefficient,name,forbidden)
    # Block every public callable from the branch owner: the crossing reads its receipt only.
    for name in dir(branch):
        if not name.startswith('_') and callable(getattr(branch,name)):
            monkeypatch.setattr(branch,name,forbidden)
    record=owner.make_record(ROOT)
    assert record==json.loads((ROOT/owner.OUTPUT).read_text())
    assert record['imported_intervals']['Delta1'][0]>0
    assert record['imported_intervals']['d'][1]<0
    assert record['scope']['physical_tuple_selected'] is False
    assert record['scope']['arbitrary_smooth_uniqueness_claimed'] is False
    assert record['scope']['zero_velocity_or_zero_R_case_analyzed'] is False


def test_crossing_receipt_rejects_modified_input(tmp_path):
    # An altered parent receipt is rejected before its coefficient values are used.
    for path in owner.INPUTS:
        target=tmp_path/path;target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes((ROOT/path).read_bytes())
    record=json.loads((tmp_path/owner.INPUTS[0]).read_text())
    path=next(iter(record['source_hashes']));target=tmp_path/path
    target.parent.mkdir(parents=True,exist_ok=True);target.write_text('changed owner')
    with pytest.raises(ValueError,match='crossing input changed'):
        owner.make_record(tmp_path)
