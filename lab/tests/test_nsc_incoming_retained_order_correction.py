import importlib.util
import json
from pathlib import Path

import mpmath as mp
import pytest

from recursive_horizons.nsc_incoming_retained_order_correction import delta_kernel,owned_middle_source
from recursive_horizons.nsc_incoming_source_tail import _riccati_at_one

ROOT=Path(__file__).resolve().parents[1]


def test_stable_increment_matches_independent_rank_one_matrices():
    c=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels'][12]
    with mp.workdps(70):
        a,r=mp.sqrt(3*mp.pi/2-4),mp.sqrt(2);m=mp.mpf(c['compact_mass']);ell=mp.mpf(c['angular_eigenvalue'])
        for sign in(1,-1):
            cs=_riccati_at_one(m,sign*ell,order=24)
            for E in(40,160,320):
                matrices=[]
                for n in(16,24):
                    S=sum(cs[j]/(2*E)**(j+1) for j in range(n))
                    v=mp.matrix([1,-1j*a*S]);matrices.append(v*v.H/(1+a*a*abs(S)**2))
                delta=matrices[1]-matrices[0];p=-E/a
                H=mp.matrix([[p,-m-1j*sign*ell/r],[-m+1j*sign*ell/r,-p]])
                got=delta_kernel(E,m,sign*ell,cs)[0]
                expected=(H*delta)[0,0]+(H*delta)[1,1]
                assert abs(got-expected)<mp.mpf('1e-55')


def test_scope_rejects_unowned_group_or_coefficients():
    with pytest.raises(ValueError,match='group12'):
        owned_middle_source(ROOT,group=14)
    with pytest.raises(ValueError,match='explicit24'):
        delta_kernel(40,0,1,[1]*16)


def test_artifact_replay_cannot_regenerate_coefficients(monkeypatch):
    record=json.loads((ROOT/'results/development/nsc-incoming-retained-order-correction.json').read_text())
    spec=importlib.util.spec_from_file_location('retained_correction_replay',ROOT/'scripts/derive_nsc_incoming_retained_order_correction.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    def forbidden(*args,**kwargs):raise AssertionError('replay cannot generate a source')
    monkeypatch.setattr(module,'prepare_retained_correction',forbidden)
    assert module.replay(ROOT/record['payload']['path'])==record
    assert record['physical_order24_error_bound'] is None
    assert not record['scope']['original_source_replaced']
