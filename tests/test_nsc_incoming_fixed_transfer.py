"""Complex fixed-transfer algebra and physical-proof receipt; no propagation."""
import importlib.util
import json
from pathlib import Path
import sympy as sp
import recursive_horizons.nsc_incoming_fixed_transfer as M
ROOT=Path(__file__).resolve().parents[1]


def test_frames_four_inverse_gap_terms_and_endpoint_algebra():
    proof=M.symbolic_proof()
    assert proof['exact_scalar_residual_count']==81
    assert proof['inverse_gap']['term_count']==4
    assert proof['inverse_gap']['telescoping_residual']=='0'
    assert proof['frame_projector']['determinant_residual']=='0'
    assert proof['endpoint_pair']['complex_amplitudes_retained']
    assert proof['endpoint_pair']['complex_Fourier_solution']==[{'Ahat':'0','Rhat':'0'}]
    assert proof['endpoint_pair']['opposite_transfer_Hermiticity_residual']=='0'


def test_complex_Fourier_amplitudes_require_both_entries():
    proof=M.endpoint_pair_identity();symbols={k:sp.Symbol(k) for k in ('a1','r1','m','ell','Ahat','Rhat')}
    L01=sp.sympify(proof['limit01'],locals=symbols);L10=sp.sympify(proof['limit10'],locals=symbols)
    a,r,m,ell,A,R=[symbols[k] for k in ('a1','r1','m','ell','Ahat','Rhat')]
    values={a:1,r:1,m:1,ell:1,A:1+2*sp.I,R:(1-sp.I)*(1+2*sp.I)}
    # A nonzero complex pair can cancel one entry, but cannot cancel both.
    assert sp.simplify(L01.subs(values))==0
    assert sp.simplify(L10.subs(values))!=0


def test_receipt_authentication_and_no_scientific_producer(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('scientific producer entered fixed-transfer proof')
    import scipy.integrate
    import recursive_horizons.nsc_incoming_vacuum_tail_bound as tail
    import recursive_horizons.nsc_retarded_radial_response as radial
    import recursive_horizons.nsc_pg_high_energy as high
    monkeypatch.setattr(scipy.integrate,'solve_ivp',forbidden)
    monkeypatch.setattr(tail.IncomingVacuumTailBoundGeometry,'__init__',forbidden)
    monkeypatch.setattr(radial,'short_response',forbidden)
    monkeypatch.setattr(high,'riccati_coefficients',forbidden)
    spec=importlib.util.spec_from_file_location('fixed_transfer_proof',ROOT/'scripts/derive_nsc_incoming_fixed_transfer.py')
    script=importlib.util.module_from_spec(spec);spec.loader.exec_module(script)
    record=script.make_record()
    assert record==json.loads((ROOT/script.OUTPUT).read_text())
    assert record['assumptions']['intrinsic_fields_identically_zero_on_Sigma']==['deltaN(z)','deltabeta(z)','deltaa(z)','deltar(z)']
    assert record['result']['requires_matching_at_every_fixed_transfer']
    assert not record['result']['finite_set_of_pairs_suffices']
    assert not record['result']['sufficient_for_full_C0_matching']
    assert record['result']['finite_energy_threshold'] is None
    assert record['result']['uniform_large_omega_error_bound'] is None
    assert not record['scope']['noncompact_axial_tangents_covered']
    assert not record['scope']['higher_normal_jets_classified']
    assert not record['scope']['finite_amplitude_exclusion_claimed']
    assert record['scope']['field_radial_source_old_producer_or_energy_scan_runs']==0
