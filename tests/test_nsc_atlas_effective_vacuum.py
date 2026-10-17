"""Check the Step-7 trace independently against the implemented bulk equations."""
import numpy as np
import pytest

from recursive_horizons.nsc_atlas_effective_vacuum import classical_vacuum_calculation
from recursive_horizons.nsc_atlas_nonlinear_bulk import bulk_operator


def test_exact_jordan_and_einstein_action_calculation():
    result = classical_vacuum_calculation()
    assert all(value == 0 for value in result['checks'].values())
    assert result['induced_lambda_Lref_squared'] == 0
    assert result['omitted_improvement_control'] != 0


@pytest.mark.parametrize('region,expected', [('cap', 12.), ('exterior', 0.)])
@pytest.mark.parametrize('alpha', [0., 0.00012983127142551443, .02])
def test_geometric_ricci_from_actual_mixed_operator(region, expected, alpha):
    # Independent point jets, not an evolved state or a PDE existence claim.
    r, C, c = .49, 1.3, .21
    rU, rV, cU, cV = -.37, .42, -.53, .19
    row = bulk_operator(dict(r=r, C=C, c=c),
                        dict(r_U=rU, r_V=rV, c_U=cU, c_V=cV),
                        alpha=alpha, region=region)
    rUV, cUV = row['r_UV'], row['c_UV']
    # Ricci of -C dU dV+r^2 dOmega^2, computed from its metric jets.
    R_E = 4*row['log_C_UV']/C+2/r**2+8*rU*rV/(C*r**2)+16*rUV/(C*r)
    u = c/r
    uU, uV = (cU-u*rU)/r, (cV-u*rV)/r
    uUV = cUV/r-(cU*rV+cV*rU+c*rUV)/r**2+2*c*rU*rV/r**3
    D = 1+alpha*u*u/6
    lU, lV = -alpha*u*uU/(3*D), -alpha*u*uV/(3*D)
    lUV = -alpha*(uU*uV+u*uUV)/(3*D)+alpha**2*u*u*uU*uV/(9*D**2)
    box_l = -4*(lUV+(rU*lV+rV*lU)/r)/C
    grad_l_squared = -4*lU*lV/C
    R_J = (R_E+3*box_l-1.5*grad_l_squared)/D
    assert np.isfinite(R_J)
    assert R_J == pytest.approx(expected, rel=0, abs=2e-14)
