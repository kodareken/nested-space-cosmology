"""Full constraint, parity and represented-GH controls for the same kick."""
import numpy as np
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_finite_parent_data import solve_initial_data
from recursive_horizons.nsc_atlas_gravitating_kick import solve_kicked_initial_data,validate_kicked_data
from recursive_horizons.nsc_atlas_parent_gh import ParentGH


def test_zero_kick_recovers_old_constraint_and_sign_is_even():
    with threadpool_limits(limits=1):
        old=solve_initial_data(4,epsilon=.25)
        zero=solve_kicked_initial_data(4,epsilon=.25,kick_amplitude=0)
        np.testing.assert_array_equal(old.coefficients,zero.coefficients)
        plus=solve_kicked_initial_data(4,epsilon=.25,kick_amplitude=1/32)
        minus=solve_kicked_initial_data(4,epsilon=.25,kick_amplitude=-1/32)
        np.testing.assert_array_equal(plus.coefficients,minus.coefficients)
        assert np.linalg.norm(plus.coefficients-old.coefficients)>0
        assert plus.metadata['sampled_psi_min']>0
        assert plus.metadata['stages'][-1]['projected_residual_norm']<2e-11


def test_independent_rational_constraint_refinement_and_positive_supply():
    with threadpool_limits(limits=1):
        cases=[validate_kicked_data(solve_kicked_initial_data(L,epsilon=.25,kick_amplitude=1/32)) for L in (4,8)]
    assert cases[1]['full_conformal_RMS_estimate']<.1*cases[0]['full_conformal_RMS_estimate']
    assert all(v['initial_normal_kinetic_over_round_vacuum_energy']>0 for v in cases)
    assert all(v['integrated_shear_constraint_over_round_vacuum']>0 for v in cases)


def test_represented_gh_constraints_and_source_signs():
    with threadpool_limits(limits=1):
        diagnostics=[]
        for L in (8,12):
            system=ParentGH(L,.25,kick_amplitude=1/32)
            y=system.initial();J=system.jets(y)
            point=system.diagnostic(0,y);diagnostics.append(point)
            baseline=ParentGH(L,.25)
            reference=baseline.diagnostic(0,baseline.initial())
            # The coarse lapse-gauge projection has the same defect before
            # the kick; check the change and convergence, not an invented floor.
            assert point['GH_C_RMS_over_H']<1.001*reference['GH_C_RMS_over_H']+1e-12
            opposite=ParentGH(L,.25,kick_amplitude=-1/32)
            other=opposite.initial()
            C=system.unpack(y);D=opposite.unpack(other)
            np.testing.assert_array_equal(C[:10],D[:10])
            np.testing.assert_array_equal(C[20:22],D[20:22])
            np.testing.assert_allclose(C[23],-D[23],atol=1e-18,rtol=0)
            # The spatial time-derivative trace is zero at represented level.
            assert np.max(abs(J.v[:,1,1]+J.v[:,2,2]+J.v[:,3,3]))<1e-14
        assert diagnostics[1]['ADM_H_RMS_over_H_squared']<.2*diagnostics[0]['ADM_H_RMS_over_H_squared']
        assert diagnostics[1]['ADM_M_RMS_over_H_squared']<.2*diagnostics[0]['ADM_M_RMS_over_H_squared']
        assert diagnostics[1]['GH_C_RMS_over_H']<.1*diagnostics[0]['GH_C_RMS_over_H']
