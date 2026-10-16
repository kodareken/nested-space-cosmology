"""Independent power-polynomial check of the local DOP853 derivative audit."""
import numpy as np
from numpy.polynomial import Polynomial
import pytest

from control_budget import ControlBudget
from check_nsc_atlas_prepared_endpoint_connection import dense_polynomial_jet


def test_nested_dense_jet_matches_independent_power_derivatives():
    control=ControlBudget(maximum_evaluations=4000,cpu_seconds=10.)
    with control.active() as budget:
        coefficients=np.array([1.,-2.,3.,-4.,5.,-6.,7.])/8
        x=Polynomial([0,1]);power=Polynomial([0])
        # Expand all products in power coordinates; runtime uses the
        # differentiated recurrence and no polynomial derivative library.
        for i,coefficient in enumerate(coefficients[::-1]):
            power=(power+coefficient)*(x if i%2==0 else 1-x)
        power+=.25
        for coordinate in (0.,.375,1.):
            got=dense_polynomial_jet(coefficients,.25,2.,2.5,2.+.5*coordinate,budget=budget)
            for key,order in (('value',0),('first',1),('second',2)):
                assert got[key]==pytest.approx(power.deriv(order)(coordinate)/.5**order,abs=2e-13)
        with pytest.raises(ValueError,match='own interval'):
            dense_polynomial_jet(coefficients,.25,2.,2.5,2.6,budget=budget)
        budget.check(1000)
    print('dense derivative control:',control.budget.evaluations,'values;',control.active_cpu_seconds,'CPU seconds')
