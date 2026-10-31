"""Independent pressure, spectral, normalization and Clifford controls."""
import math

import numpy as np
import pytest

from recursive_horizons.dirac.charged_throat import CommonCoefficients
from recursive_horizons.nsc_atlas_inherited_charge import harmonic_sector
from recursive_horizons.nsc_atlas_valve_source_gate import (
    radial_source_gate,geometric_spinor_candidate,mmp_massless_field_gate,
    dirac_clifford_control,neutral_thermal_ec_benchmark,
    axial_pair_correlation_witness,
    parent_mode_axial_kernel,
)


def test_closed_de_sitter_belly_bounces_without_negative_nec():
    H=.3;kappa=.7;R=1/H;rho=3*H*H/kappa
    row=radial_source_gate(radius=R,radius_rate=0.,density=rho,radial_pressure=-rho,
                          radial_spacetime_Hessian_over_radius=-H*H,kappa=kappa,
                          actual_acceleration=H*H*R)
    assert row['positive_acceleration_supported']
    assert row['pressure_capacity']==pytest.approx(3.)
    assert row['radial_equation_residual']<1e-15
    assert row['classical_radial_NEC']==0.
    assert row['extra_total_NEC_at_fixed_geometry_must_be_less_than']>0


def test_contracting_receiver_requires_more_than_turning_point_pressure():
    row=radial_source_gate(radius=.5,radius_rate=-2.,density=1.,radial_pressure=.25,
                          radial_spacetime_Hessian_over_radius=3.)
    assert not row['positive_acceleration_supported']
    assert row['pressure_for_positive_acceleration_must_be_less_than']==-20.
    assert row['pressure_at_turning_point_must_be_less_than']==-4.
    assert row['extra_total_NEC_at_fixed_geometry_must_be_less_than']==-7.25
    with pytest.raises(ValueError):radial_source_gate(radius=0.,radius_rate=0.,density=1.,radial_pressure=-1.,radial_spacetime_Hessian_over_radius=0.)


def test_warped_spinor_gap_matches_independent_scalar_and_curvature_owners():
    row=geometric_spinor_candidate();scalar=harmonic_sector()
    np.testing.assert_allclose(row['lowest_mass_squared'],3*scalar['KK_mass_squared_increment'],rtol=2e-15)
    assert row['Friedrich_internal_eigenvalue_squared_lower_bound']==pytest.approx(45/11)
    assert row['modes'][0]['internal_absolute_eigenvalue']**2==pytest.approx(45/11)
    assert [m['multiplicity_per_internal_Dirac_sign'] for m in row['modes']]==[64,768,4992]
    assert row['lowest_four_dimensional_Dirac_fields_both_internal_signs']==128
    assert row['ground_charge_multiplicities_per_Dirac_sign']=={'minus_half':32,'plus_half':32}
    assert row['mass_radion_derivative']==pytest.approx(row['logarithmic_mass_radion_derivative']*row['lowest_mass'])
    assert row['source_population_selected'] is False
    assert row['SI_mass_selected'] is False


def test_charge_exchange_is_a_proper_spin_rotation():
    X=np.array([[0.,1.],[1.,0.]],complex);Y=np.array([[0.,-1j],[1j,0.]],complex);Z=np.diag([1.,-1.]);I=np.eye(2)
    def tensor(factors):
        out=np.array([[1.]],complex)
        for f in factors:out=np.kron(out,f)
        return out
    gamma=[]
    for j in range(6):
        gamma.extend(tensor([Z]*j+[p]+[I]*(5-j)) for p in (X,Y))
    gamma.append(tensor([Z]*6))
    J=-.5j*gamma[0]@gamma[1];U=gamma[1]@gamma[12]
    np.testing.assert_allclose(U.conj().T@U,np.eye(64),atol=0)
    np.testing.assert_allclose(U@J@U.conj().T,-J,atol=0)
    assert np.count_nonzero(np.diag(J).real==.5)==32
    assert np.count_nonzero(np.diag(J).real==-.5)==32


def test_mmp_mass_gate_reuses_existing_common_action_normalization():
    row=geometric_spinor_candidate();m2=row['lowest_mass_squared']
    # A formal common-unit check, not selection of the physical free g.
    e2=row['minimal_charge_coupling_squared_over_g_grav']
    common=CommonCoefficients(einstein=.5,gauge=1/(4*e2),vacuum=0.)
    for Q in (1,2,7):
        mapped=common.mmp_map(Q);gate=mmp_massless_field_gate(row,Q)
        assert m2*mapped['r_e_squared']==pytest.approx(gate['lowest_mass_times_extremal_radius_squared'])
        assert gate['lowest_mass_times_extremal_radius_squared']==pytest.approx(18*Q*Q/13)
        assert gate['original_MMP_massless_field_assumption_satisfied'] is False


def test_literal_thermal_ec_wick_trace_and_full_nec_normalization():
    clifford=dirac_clifford_control();row=neutral_thermal_ec_benchmark()
    assert clifford['Clifford_defect']==0.
    assert clifford['axial_pair_trace']==16.
    assert clifford['literal_Wick_A_squared_over_m_squared_sigma_squared']==-1.
    ratio=row['contact_NEC_hat_over_g_squared_T_hat_fourth']/row['free_NEC_hat_over_g_T_hat_fourth']
    assert ratio==pytest.approx(row['contact_to_free_NEC_ratio_over_g'])
    assert ratio==pytest.approx(2*row['reference_half_Wick_ratio_over_g'])
    assert row['contact_NEC_sign']=='positive'
    assert row['total_NEC_sign_in_controlled_regime']=='positive'
    assert row['g_selected'] is False and row['temperature_selected'] is False
    assert row['supports_minimum_neck_turnaround'] is False


def test_equal_free_covariance_does_not_select_connected_torsion_source():
    row=axial_pair_correlation_witness()
    assert row['CAR_exact'] and row['full_normal_order_identity_exact']
    a=row['states']['singlet'];b=row['states']['isotropic_triplet']
    for state in (a,b):
        assert state['N']=='2' and state['Q']=='0'
        assert state['free_vacuum_subtracted_energy_over_m']=='2'
        assert state['rotation_invariant'] and state['C_invariant']
        assert state['one_body_covariance']=='diag(1/2,1/2,1/2,1/2)'
        assert state['anomalous_covariance']=='zero'
    assert a['full_vacuum_Wick_AA_times_V_squared']=='-10'
    assert b['full_vacuum_Wick_AA_times_V_squared']=='2'
    assert row['triplet_minus_singlet_AA_times_V_squared']=='12'
    assert row['triplet_minus_singlet_contact_H_times_V_over_kappa']=='9/4'


def test_parent_kinetic_mode_kernel_differs_from_rest_projection():
    # Independent exact Fock/Dirac audit at m/E=3/5, lambda/E=4/5 gives
    # raw AA singlet=-36/25, triplet=8; their common-subtraction gap=236/25.
    row=parent_mode_axial_kernel(mass=3.,spatial_eigenvalue=4.,density_overlap=.2)
    assert row['energy']==5.
    assert row['local_quartic_gap_over_density_squared']==pytest.approx(236/25)
    assert row['triplet_minus_singlet_energy_hat_over_g']==pytest.approx((3/16)*(236/25)*.2)
    rest=parent_mode_axial_kernel(mass=3.,spatial_eigenvalue=0.,density_overlap=.2)
    assert rest['local_quartic_gap_over_density_squared']==12.
    assert rest['triplet_minus_singlet_energy_hat_over_g']==pytest.approx(9*.2/4)
    assert row['source_state_selected'] is False and row['physical_pressure_or_NEC_derived'] is False
