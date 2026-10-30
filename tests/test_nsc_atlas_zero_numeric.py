"""Fixed-model tails, retained state, pressure-work and reference-only quantities."""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest
import sympy as s

from recursive_horizons.nsc_atlas_horizon_response import solve_source
from recursive_horizons.nsc_atlas_zero_interface import (
    tail_state,thermodynamic_record,winding_and_charge,ccc_record,spectrum_record,pet_comparison_only)


@pytest.fixture(scope='module')
def source():return solve_source()


def test_nonzero_state_survives_zero_excitation_density(source):
    a=tail_state(source,6.);b=tail_state(source,8.)
    for name in ('amplitude_quadrature','phase_quadrature'):
        np.testing.assert_array_less(abs(np.array(a[name])-b[name]),np.array(a['each_quadrature_future_error_upper_bound'])+1e-10)
    assert np.linalg.norm(b['amplitude_quadrature'])>.1
    assert b['late_phase_average_density_estimate_over_epsilon_squared_vacuum']<4e-10
    assert min(b['rescaled_oscillator_energy'])>29
    assert b['Omega']<a['Omega']


def test_tail_spectral_and_density_exponents_are_distinct(source):
    a=tail_state(source,6.);b=tail_state(source,8.)
    log_a=math.log(math.cosh(8.)/math.cosh(6.))
    assert math.log(b['H_times_spatial_wavelengths'][0]/a['H_times_spatial_wavelengths'][0])/log_a==pytest.approx(1.)
    slope=math.log(b['late_phase_average_density_estimate_over_epsilon_squared']/a['late_phase_average_density_estimate_over_epsilon_squared'])/log_a
    assert slope==pytest.approx(-3.,abs=1e-7)
    assert source.parent.mass_squared-2.25*source.parent.H**2>0


def test_field_thermodynamic_work_identity(source):
    r=thermodynamic_record(source)
    assert max(x['work_identity_absolute_residual'] for x in r['samples'])<1e-10
    assert r['phase_shift_alone_is_amplitude_decay'] is False
    assert 0<r['late_spiral_radius_ratio_per_turn']<1
    assert r['late_spiral_inverse_amplitude_ratio_per_turn']*r['late_spiral_radius_ratio_per_turn']==pytest.approx(1.)
    response=r['late_retarded_response']
    assert response['mean_absorbed_power_over_F0_squared']>0
    assert response['power_balance_absolute_residual']<1e-14
    assert response['phase_lag_radians']==pytest.approx(math.pi/2,abs=1e-13)


def test_conformal_matching_benchmark_is_not_an_inferred_child(source):
    r=ccc_record(source.parent)
    assert r['horizon_zero']['normal_norm']==0
    assert r['conformal_zero']['normal_norm']<0
    assert r['conformal_zero']['same_surface_as_observer_horizon'] is False
    assert r['reciprocal_comparison']['dimensionless_kappa_H_squared_rho_a_fourth']==3
    assert r['reciprocal_comparison']['same_action_matter_continuation'] is False
    assert r['actual_scalar']['conformal_singular_coefficient']>8300
    assert [(x['crossover_value'],x['crossover_eta_derivative']) for x in r['regular_operator_comparisons']['regular_real_unit_modes']]==[(-1.,0.),(0.,3.)]


def test_nonzero_circulation_and_energy_do_not_choose_charge():
    r=winding_and_charge();q=r['state_charge'];w=r['reversed_winding']
    gauge=r['inherited_KK_connection']
    assert q['same_real_circulation_helicity']=='0'
    assert q['same_real_circulation_mean_norm_squared']=='8/(3*pi**4)'
    assert q['neutral_state_expected_internal_charge']==0
    assert q['physical_baryon_photon_eta'] is None
    assert w['standing_wave_E_per_A_squared']>0
    assert w['opposite_winding_Q_per_A_squared'][0]==-w['opposite_winding_Q_per_A_squared'][1]
    assert w['each_positive_E_per_A_squared']/abs(w['opposite_winding_Q_per_A_squared'][0])==pytest.approx(w['model_mass'])
    assert gauge['gauge_coupling_squared_over_kappa_mass_increment']==pytest.approx(13/12,abs=5e-15)
    assert gauge['gauge_coupling_squared_over_g_grav']==pytest.approx(2.433310295207691)
    assert gauge['neutral_preparation_gauge_invariant_horizon_potential']==0


def test_detector_constants_remain_comparison_only_and_do_not_select_spectrum():
    r=pet_comparison_only();spectrum=spectrum_record()
    assert r['reference_only'] is True and r['NSC_prediction'] is False
    assert r['photon_reference_energy_keV']==pytest.approx(510.99895069)
    assert r['photon_reference_frequency_Hz']*r['photon_reference_wavelength_m']==pytest.approx(299792458.)
    assert r['example_200ps_contains_reference_carrier_cycles']>1e10
    assert spectrum['highest_l_selected'] is False
    assert spectrum['Planck_length_proved_as_hard_wavelength_cutoff'] is False


def test_internal_gauge_gap_relation_preserves_the_unselected_scale():
    d,r,xi,kappa,g,L=s.symbols('d r xi kappa g L',positive=True)
    g2=kappa*(d+1)*s.exp(-(d+2)*xi)/r**2
    gap=d*s.exp(-(d+2)*xi)/r**2
    assert s.simplify(g2/(kappa*gap))==(d+1)/d
    restored=g2.subs({kappa:g*L**2,r:r*L},simultaneous=True)
    assert s.simplify(restored/g-(d+1)*s.exp(-(d+2)*xi)/r**2)==0


def test_charge_blind_regional_selection_preserves_zero_mean_charge():
    C=s.kronecker_product(s.eye(2),s.Matrix([[0,1],[1,0]]))
    Q=s.kronecker_product(s.eye(2),s.diag(1,-1))
    P=s.kronecker_product(s.diag(1,0),s.eye(2))
    a,b,c=s.symbols('a b c',real=True)
    rho=s.kronecker_product(s.Matrix([[a,c],[c,b]]),s.eye(2))
    assert C*P-P*C==s.zeros(4)
    assert C*rho*C-rho==s.zeros(4)
    assert C*Q*C+Q==s.zeros(4)
    assert s.trace(P*rho*Q)==0


def test_horizon_can_be_a_signal_flow_zero_without_zero_vacuum_energy():
    H,R,velocity=s.symbols('H R velocity',positive=True)
    null_equation=-1+(velocity-H*R)**2
    assert s.expand(null_equation.subs(velocity,H*R+1))==0
    assert s.expand(null_equation.subs(velocity,H*R-1))==0
    assert (1-H*R).subs(R,1/H)==0
    kappa=s.symbols('kappa',positive=True)
    assert s.ask(s.Q.positive(3*H*H/kappa)) is True


def test_circle_minimizes_stationary_mode_energy_at_fixed_charge():
    x,y,px,py,w=s.symbols('x y px py w',real=True)
    E=(px*px+py*py+w*w*(x*x+y*y))/2
    Q=x*py-y*px
    assert s.expand(E-w*Q-((px+w*y)**2+(py-w*x)**2)/2)==0
    assert s.expand(E+w*Q-((px-w*y)**2+(py+w*x)**2)/2)==0


def test_retarded_phase_lag_has_the_positive_work_channel():
    omega,mu,gamma=s.symbols('omega mu gamma',positive=True,real=True)
    response=1/(mu*mu-omega*omega-s.I*gamma*omega)
    absorptive=s.simplify(s.im(s.expand_complex(response)))
    expected=gamma*omega/((mu*mu-omega*omega)**2+gamma*gamma*omega*omega)
    assert s.simplify(absorptive-expected)==0
    power=omega*absorptive/2
    damping=gamma*omega*omega*s.simplify(response*s.conjugate(response))/2
    assert s.simplify(power-damping)==0


def test_saved_source_bindings():
    root=Path(__file__).resolve().parents[1];path=root/'results/summaries/nsc-atlas-zero-interface-v1.json'
    if not path.exists():pytest.skip('source record follows final freeze')
    r=json.loads(path.read_text(encoding='utf-8'))
    for name,digest in r['source_sha256'].items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest
