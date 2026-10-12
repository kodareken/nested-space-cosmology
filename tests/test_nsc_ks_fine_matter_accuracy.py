"""Completed-coverage guard and independent normal/shift contraction signs."""
import json
from pathlib import Path
import sys

import pytest
pytest.importorskip('flint')
from flint import arb,acb,acb_mat,ctx

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import derive_nsc_ks_fine_matter_accuracy as M


def test_partial_history_cannot_emit_a_numerical_matter_certificate(monkeypatch,tmp_path):
    path=tmp_path/'partial.json'
    path.write_text(json.dumps({'all_segments_covered':False}))
    monkeypatch.setattr(M.V,'OUTPUT',path)
    with pytest.raises(ValueError,match='partial history'):M.compute()


def test_coherent_pauli_and_momentum_contraction_signs():
    with ctx.workprec(120):
        E=arb(7)/10
        C=acb_mat([[arb(3)/5,acb(0,arb(1)/10)],[-acb(0,arb(1)/10),arb(2)/5]])
        F=acb_mat([[1,0],[0,1]])
        Fz=-acb(0,E)*F
        N,beta=M.insertion(F,Fz,C,arb(1),arb(2),arb(3))
        a=(3*arb.pi()/2-4).sqrt();r=arb(2).sqrt()
        expected_N=3*(arb(2)/(5*r)+E/(5*a))
        assert N.real.overlaps(expected_N) and N.imag.contains(0)
        assert beta.real.overlaps(-3*E) and beta.imag.contains(0)


def test_support_hull_accepts_the_actual_single_direction_control():
    single=[(arb(1),arb(2))]
    assert M.support_hull(single)==single[0]
    assert M.support_hull(single+[(arb(0),arb(3))])==(arb(0),arb(3))
    with pytest.raises(ValueError,match='support'):
        M.support_hull([])
