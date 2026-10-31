"""Constraint, geometry and energy checks for the one-mode initial neck."""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_parent_circulation import VOLUME
from recursive_horizons.nsc_atlas_valve_neck import (
    solve_necks,neck_diagnostics,critical_amplitude,zonal_rule,
)


def test_exact_vacuum_limit_and_round_s3_normalization():
    u,w,labels,Y=zonal_rule(16,90)
    np.testing.assert_allclose(w.sum(),VOLUME,atol=1e-14)
    np.testing.assert_allclose(Y.T@(w[:,None]*Y),np.eye(len(labels)),atol=5e-14)
    data=solve_necks(16,(0.,))[0];row=neck_diagnostics(data);p=data.parent
    np.testing.assert_allclose(data.values(np.linspace(-1,1,13)),1.,atol=1e-15)
    assert row['full_conformal_PDE_RMS']<1e-13
    assert row['dot_theta_pm_over_H_squared_direct']==2.
    assert row['dot_theta_pm_over_H_squared_constraint']==2.
    assert row['physical_excitation_energy']==0.
    np.testing.assert_allclose(row['neck_areal_radius'],1/p.H)
    np.testing.assert_allclose(row['proper_total_volume'],VOLUME/p.H**3)
    assert row['neck_radial_pressure']==row['neck_tangential_pressure']==-p.vacuum


def test_finite_neck_constraint_slope_refinement_and_energy_accounting():
    medium,fine=(neck_diagnostics(solve_necks(L,(1.,))[0]) for L in (24,32))
    assert fine['full_conformal_PDE_RMS']<1e-8
    assert fine['full_conformal_PDE_RMS']<medium['full_conformal_PDE_RMS']/100
    assert fine['slope_identity_residual']<3e-8
    assert fine['dot_theta_pm_over_H_squared_direct']<-1.3
    assert abs(fine['dot_theta_pm_over_H_squared_direct']-medium['dot_theta_pm_over_H_squared_direct'])<7e-6
    assert fine['analytic_positivity_lower_bound']>0
    assert fine['psi_neck']<fine['psi_neck_strict_trapping_threshold']-.06
    for total,half in (('proper_total_volume','proper_half_volume'),
                       ('physical_excitation_energy','physical_half_excitation_energy'),
                       ('physical_vacuum_energy','physical_half_vacuum_energy'),
                       ('physical_total_energy','physical_half_total_energy')):
        assert fine[total]==2*fine[half]
    np.testing.assert_allclose(fine['physical_total_energy'],fine['physical_excitation_energy']+fine['physical_vacuum_energy'])
    assert fine['physical_excitation_over_reference_round_vacuum_energy']<.43
    assert fine['supplied_round_reference_excitation_over_reference_vacuum']==1.
    np.testing.assert_allclose(fine['neck_radial_pressure']-fine['neck_tangential_pressure'],2*fine['neck_gradient_energy_density'])
    np.testing.assert_allclose(fine['neck_isotropic_pressure'],(fine['neck_radial_pressure']+2*fine['neck_tangential_pressure'])/3)


def test_numerical_critical_amplitude_and_independent_slopes():
    low,high=(neck_diagnostics(row) for row in solve_necks(24,(.75,1.)))
    assert low['dot_theta_pm_over_H_squared_direct']>0
    assert high['dot_theta_pm_over_H_squared_direct']<0
    row=critical_amplitude(32)
    assert .75<row['epsilon_crit']<1.
    assert abs(row['dot_theta_pm_over_H_squared_constraint'])<1e-9
    assert abs(row['dot_theta_pm_over_H_squared_direct'])<1e-8
    assert row['full_conformal_PDE_RMS']<1e-9


def test_continuation_order_derivatives_and_input_rejection():
    data=solve_necks(16,(.5,.25,.5),max_epsilon_step=.125)
    np.testing.assert_allclose(data[0].coefficients,data[-1].coefficients,atol=3e-13)
    assert not data[0].coefficients.flags.writeable
    assert data[0].values(0.,1)==0.
    # Independent centered derivative check away from the reflection neck.
    x=.3;h=1e-4
    second=(data[0].values(x+h)-2*data[0].values(x)+data[0].values(x-h))/h**2
    np.testing.assert_allclose(second,data[0].values(x,2),atol=3e-8)
    for degree in (3,5,66,True):
        with pytest.raises(ValueError,match='degree'):solve_necks(degree)
    for eps in (-.1,1.1,math.nan,True):
        with pytest.raises(ValueError,match='epsilon'):solve_necks(8,(eps,))
    with pytest.raises(ValueError,match='validation order'):neck_diagnostics(data[0],validation_order=30)
