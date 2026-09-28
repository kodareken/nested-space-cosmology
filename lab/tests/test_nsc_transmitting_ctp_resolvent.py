"""Focused checks of the new memory convention and archived-field bridge."""
import json
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons.nsc_transmitting_ctp_resolvent import (
    ArchivedResolventField, CTPRetardedMemory, TransmittingCrossResolvent,
    ctp_completion_inventory, frequency_transfer_driven, frequency_transfer_weak,
)

ROOT=Path(__file__).resolve().parents[1]


def inputs():
    r=json.loads((ROOT/'results/development/nsc-transmitting-cross-resolvent.json').read_text())
    a=np.load(ROOT/r['payload']['path'],allow_pickle=False)
    return r,a


def test_memory_preserves_direction_and_rejects_full_jet_claim():
    _,a=inputs()
    m=CTPRetardedMemory(a['upper_response'][0],a['lower_response'][0])
    assert np.linalg.norm(m.G_retarded[:2,2:])==0
    assert np.linalg.norm(m.G_retarded[2:,:2])>0.1
    assert np.linalg.norm(m.G_retarded@m.D_retarded-np.eye(4))<3e-11
    with pytest.raises(ValueError,match='preparation'): m.as_endpoint_branch_jets()
    assert ctp_completion_inventory()['lesser_greater_and_Keldysh'] is None


def test_temporal_vertex_uses_distinct_frequencies_on_locked_domain():
    r,a=inputs();index=12;c=r['channels'][index]['channel']
    o=TransmittingCrossResolvent(c['compact_mass'],c['angular_eigenvalue'])
    zi=complex(*r['parameters']['z']);zo=complex(.6,zi.imag)
    inc=ArchivedResolventField(o,zi,a['rho_quadrature'],a['rho_weights'],a['upper_fields'][index])
    out_adj=o.solve(zo.conjugate(),with_jets=False)
    d,seam=frequency_transfer_driven(o,zo,inc)
    w=frequency_transfer_weak(o,out_adj,inc,a['rho_quadrature'],a['rho_weights'])
    assert np.linalg.norm(d-w)<3e-11
    assert np.linalg.norm(d-a['upper_metric_jets'][index])>1e-4
    assert np.linalg.norm(d[:,:2,2:])==0 and seam==0
    with pytest.raises(ValueError,match='same nonzero imaginary'):
        frequency_transfer_driven(o,complex(.6,.35),inc)
