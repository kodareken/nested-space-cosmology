"""Exact retained-source bound, without stress or metric integration."""
import json
from fractions import Fraction
from pathlib import Path

import pytest
import sympy as sp

from recursive_horizons.nsc_incoming_constraint_gate import retained_incoming_flux_bound

ROOT=Path(__file__).resolve().parents[1]


def input_data():
    load=lambda p:json.loads((ROOT/p).read_text())
    return (load('results/development/nsc-pg-retained-covariance.json')['locked_inputs'],
            load('results/development/nsc-mode-resolved-cauchy-state.json')['channels'],
            load('results/development/nsc-compact-matched-restart.json')['scattering_provenance']['config'])


def sewing_current_identity():
    f,n,s,rr,ri,T=sp.symbols('f n s rr ri T',real=True)
    S=sp.Matrix([[0,1,0],[rr+sp.I*ri,0,sp.sqrt(T)]])
    C=sp.Matrix([[f,-sp.I*s,0],[sp.I*s,1-f,0],[0,0,n]])
    # Use a real transmission square root to avoid irrelevant symbolic
    # branch choices; the inherited physical T is in [0,1].
    t=sp.symbols('t',nonnegative=True)
    S=S.subs(sp.sqrt(T),t)
    current=sp.expand(sp.trace(S*C*S.conjugate().T))
    expected=1+t*t*(n-f)
    return {'horizon_coherence_overlap':str(sp.simplify((S.conjugate().T*S)[0,1])),
            'trace_mod_current_unitarity':str(sp.simplify(current-expected-f*(rr*rr+ri*ri+t*t-1)))}


def test_source_coherence_and_signed_normalization():
    assert sewing_current_identity()=={'horizon_coherence_overlap':'0','trace_mod_current_unitarity':'0'}
    r=retained_incoming_flux_bound(*input_data())
    assert r['Killing_power_upper_bound'] < -.022
    assert sorted(v['multiplicity'] for v in r['massive_emission_bounds'].values())==[512,512]
    witness=r['rational_witness']
    assert Fraction(witness['strict_negative_power_margin'])>Fraction(22,1000)
    assert Fraction(witness['T01_strict_lower'])>Fraction(12,10000)
    assert witness['normalized_shift_miss_lower_decimal']>1e8*3e-11


def test_does_not_generalize_past_declared_inventory_or_parameters():
    ledger,channels,prep=input_data()
    with pytest.raises(ValueError,match='inventory'):retained_incoming_flux_bound(ledger,channels[:-1],prep)
    with pytest.raises(ValueError,match='ledger'):retained_incoming_flux_bound({**ledger,'Omega':1.},channels,prep)
