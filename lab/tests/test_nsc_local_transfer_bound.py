from fractions import Fraction as F

import pytest
import sympy as sp

from recursive_horizons.nsc_local_transfer_bound import (
    two_energy_tail_coefficients,normal_insertion_bounds,linear_large_transfer_bound,
)


def test_double_integral_powers_and_constants_independently():
    E,E0,omega=sp.symbols('E E0 omega',positive=True)
    for p in (4,6):
        density=2*sp.integrate(2*sp.integrate(omega**(-p),(omega,E/2,sp.oo)),(E,E0,sp.oo))
        current=2*sp.integrate(2*sp.integrate((2*E+omega)*omega**(-p),(omega,E/2,sp.oo)),(E,E0,sp.oo))
        result=two_energy_tail_coefficients(p)
        assert sp.simplify(density-sp.Rational(result['density'])*E0**result['density_power'])==0
        assert sp.simplify(current-sp.Rational(result['current'])*E0**result['current_power'])==0
    assert two_energy_tail_coefficients(4)['current']==F(88,3)
    assert two_energy_tail_coefficients(6)['current']==F(112,5)


def test_collar_and_retained_zero_angular_response():
    assert normal_insertion_bounds()=={'w':F(9,40000),'U':F(27,1600000000)}
    result=linear_large_transfer_bound(order=6,energy_cut=32,w_derivative_l1=10,U_derivative_l1=20,
        mass_upper=2,angular_abs_upper=0,multiplicity=12,axial_lower=F(4,5),radius_lower=F(7,5))
    assert result['N_upper']==result['beta_upper']==0
    assert result['finite_history_remainder_bound'] is None
    assert result['physical_local_gate']=='OPEN'


def test_required_norms_and_weighted_integrability():
    for p in (0,1,2,3,True,4.5):
        with pytest.raises(ValueError):two_energy_tail_coefficients(p)
    with pytest.raises(ValueError):
        linear_large_transfer_bound(order=6,energy_cut=32,w_derivative_l1=None,U_derivative_l1=0,
            mass_upper=2,angular_abs_upper=3,multiplicity=12,axial_lower=F(4,5),radius_lower=F(7,5))
