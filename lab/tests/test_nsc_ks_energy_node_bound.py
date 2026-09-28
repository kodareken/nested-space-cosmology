"""Exact polynomial controls for actual-node, whole-interval remainders."""
import numpy as np
import pytest
pytest.importorskip('flint')
from flint import arb,ctx
from recursive_horizons.nsc_ks_energy_node_bound import actual_node_remainder_factor,multiply_derivative_bound
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_ks_energy_propagator import chebyshev_energy_nodes


def test_nodal_quadratic_has_exact_sharp_interval_bound():
    r=actual_node_remainder_factor(np.array([-.5,.5]),(-1.,1.))
    with ctx.workprec(160):
        assert restored_upper(r['node_polynomial_sup_upper'])==arb(3)/4
        assert restored_upper(r['remainder_factor_upper'])==arb(3)/8
        assert restored_upper(multiply_derivative_bound(r,2))==arb(3)/4
    assert r['node_displacement_included']
    assert r['node_solve_error_bound'] is None
    assert r['interpolation_arithmetic_error_bound'] is None


def test_rounded_degree48_roots_have_near_ideal_bound_over_whole_interval():
    nodes=chebyshev_energy_nodes((-320.,320.),48)
    r=actual_node_remainder_factor(nodes,(-320.,320.))
    with ctx.workprec(160):
        ideal=arb(2)**-48
        upper=restored_upper(r['normalized_node_polynomial_sup_upper'])
        # Sorted factor multiplication has cancellation, enclosed at160 bits.
        assert upper>=ideal
        assert upper<ideal*(1+arb('1e-10'))
    with pytest.raises(ValueError,match='distinct'):
        actual_node_remainder_factor([.2,.2],(0.,1.))


def test_scaled_and_perturbed_nodes_are_bounded_without_ideal_root_assumption():
    interval=(2.,6.)
    nodes=np.array([2.5,3.7,5.1,5.9])
    r=actual_node_remainder_factor(nodes,interval)
    with ctx.workprec(160):
        upper=restored_upper(r['node_polynomial_sup_upper'])
        for x in np.linspace(*interval,501):
            value=arb(1)
            for node in nodes:value*=arb(float(x))-arb(float(node))
            assert value.abs_upper()<=upper
    assert not r['ideal_Chebyshev_nodes_assumed']
    assert r['physical_local_gate']=='OPEN'
