"""Six-term endpoint identities, complex maps and the owned twenty slots."""
import importlib.util
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
SPEC=importlib.util.spec_from_file_location('jet_rigidity',ROOT/'scripts/derive_nsc_incoming_jet_rigidity.py')
R=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(R)
from recursive_horizons.nsc_incoming_jet_rigidity import inverse_gap_six,normal_chain,endpoint_maps,slot_inventory
from recursive_horizons.nsc_incoming_cauchy_jets import NORMAL_ENTRIES


def test_six_term_telescoping_and_vanishing_lower_endpoint_terms():
    proof=inverse_gap_six()
    assert proof['term_count']==6 and proof['telescoping_residual']=='0'
    for row in proof['endpoints']:
        assert row['qk_endpoint_residual']=='0'
        assert row['lower_q_endpoint_residuals']==['0']*row['normal_order']
    assert all(row['normal_chain_residual']=='0' for row in normal_chain()['orders'])


def test_all_four_complex_maps_are_injective_and_retain_adjoint_parity():
    parent,q1=R.inputs();proof=endpoint_maps(tuple(tuple(row) for row in q1))
    assert proof['complex_Fourier_amplitudes']
    for row in proof['orders']:
        assert row['coefficient_identity_residuals']==['0','0']
        assert row['determinant_identity_residual']=='0'
        assert row['complex_kernel']==[{'A_k':'0','R_k':'0'}]
        assert row['opposite_transfer_adjoint_residual']=='0'
        assert row['normal_order'] in range(1,5)


def test_twenty_slots_exactly_match_the_existing_domain():
    proof=slot_inventory();expected={(field,k,j) for field in ('a','r') for k,j in NORMAL_ENTRIES}
    actual={(row['field'],row['normal_order'],row['spatial_order']) for row in proof['slots']}
    assert actual==expected and len(actual)==proof['count']==20
    assert proof['per_normal_order_per_field']==[4,3,2,1]
    assert all(row['necessary_tangent_change']==0 for row in proof['slots'])


def test_authenticated_receipt_keeps_functional_and_first_order_scope():
    record=R.check()
    assert record['result']['slot_count']==20
    assert record['result']['physical_pair_remainder']=='O(E^-6) for six inverse-gap terms'
    assert record['scope']['N_beta_normal_jets_fixed_as_functions_through_four']
    assert not record['scope']['bare_point_germ_sufficient']
    assert not record['scope']['finite_amplitude_exclusion']
    assert not record['result']['sufficient_for_full_matching']
    assert record['scope']['physical_field_source_radial_producer_or_energy_scan_runs']==0
    assert record['residuals']['exact_symbolic_identities']==0
