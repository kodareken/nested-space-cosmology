"""Physical source-mode derivative controls, never a packet-state replacement."""
import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from recursive_horizons.nsc_pg_ctp_mode_jets import PGSourceModeHistoryDerivative,continue_mode_columns
from derive_nsc_pg_ctp_mode_jets import evaluate,load


def payload():
    r=json.loads((ROOT/'results/development/nsc-pg-ctp-mode-jets.json').read_text())
    return load(ROOT/r['payload']['path'])


def family(a,name):return {k[len(name)+1:]:v for k,v in a.items() if k.startswith(name+'/')}


def test_actual_mode_vertices_agree_with_direct_hamiltonian_variation():
    rows,maxima=evaluate(payload())
    assert len(rows)==63 and {r['group'] for r in rows.values()}==set(range(33))
    assert maxima['weak_vs_direct']<3e-8 and maxima['spatial_refinement']<3e-8
    assert maxima['energy_exchange_hermiticity']<3e-11
    assert maxima['CTP_trace']<3e-11 and maxima['CTP_tangent']<3e-11
    assert min(r['positive_frequency_transfer_norm'] for r in rows.values())>1e-3
    assert maxima['seam']==0


def test_conditional_time_kernel_is_tangent_not_a_selected_history():
    d=family(payload(),'32_-1')
    owner=PGSourceModeHistoryDerivative(d['energies'],d['vertex_72'],d['source'],d['projector'])
    # Coordinates used only to verify the exact phase identity. No metric or
    # state is stepped and these numbers are never recorded as a physical time.
    ti,tf,tau=-.2,.7,.1
    U=np.diag(np.exp(-1j*np.repeat(d['energies'],3)*(tf-ti)))
    for derivative in owner.derivative_kernel(ti,tf,tau):
        assert np.linalg.norm(U.conj().T@derivative+derivative.conj().T@U)<3e-11
    with pytest.raises(ValueError,match='raw KS'):
        owner.as_endpoint_branch_jets()
    with pytest.raises(ValueError,match='ordered'):
        owner.derivative_kernel(1.,0.,.5)


def test_source_fibers_cannot_be_replaced_by_seed_or_packet_covariance():
    with pytest.raises(ValueError,match='mode columns'):
        continue_mode_columns(np.array([.2]),1.,1.,np.eye(2))
    d=family(payload(),'13_1')
    with pytest.raises(ValueError,match='matching physical source'):
        PGSourceModeHistoryDerivative(d['energies'],d['vertex_72'],np.eye(8),d['projector'])
    M=d['vertex_72'].copy();M[0,2,2]=1 # closed infinity source at |E|<m
    with pytest.raises(ValueError,match='closed infinity'):
        PGSourceModeHistoryDerivative(d['energies'],M,d['source'],d['projector'])
