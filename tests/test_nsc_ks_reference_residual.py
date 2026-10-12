import numpy as np
import pytest
from scipy.integrate import DOP853, quad

from recursive_horizons.nsc_ks_reference_residual import reference_segment_defect
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_ks_trajectory import TrajectorySegment
from recursive_horizons.nsc_retarded_radial_response import ks_generator


def test_constant_interpolant_encloses_independent_integrated_defect():
    state = np.array([1., .5j])
    segment = TrajectorySegment(1.01,1.009,state,state,np.zeros((6,2),complex))
    value = reference_segment_defect(segment,[.7],[.4],.2,.8)
    actual = quad(lambda rho: .4*np.linalg.norm(
        ks_generator(rho,np.array([.7]),.2,.8)[0]@state),1.009,1.01,
        epsabs=1e-14)[0]
    bound = float(restored_upper(value['weighted_row_residual_integral']))
    assert actual <= bound < 1.02*actual
    assert value['physical_source_error_included'] is False


def test_actual_decreasing_dop853_reference_defect_and_weight_scaling():
    energies = np.array([.7,1.3])
    initial = np.array([[1.,.2j],[.1j,.7]])
    def rhs(t,y):
        g = ks_generator(t,energies,.2,.8)
        return np.einsum('sab,bs->as',g,y.reshape(2,2)).ravel()
    solver = DOP853(rhs,1.01,initial.ravel(),1.0099,rtol=1e-13,atol=1e-15)
    solver.step()
    dense = solver.dense_output()
    segment = TrajectorySegment(solver.t_old,solver.t,solver.y_old,solver.y,dense.F[1:])
    v = reference_segment_defect(segment,energies,[.4,.5],.2,.8)
    v2 = reference_segment_defect(segment,energies,[.8,1.],.2,.8)
    upper = float(restored_upper(v['weighted_frobenius_defect_sup']))
    for x in np.linspace(0,1,11):
        rho = segment.rho_start+x*(segment.rho_end-segment.rho_start)
        r = (segment.evaluate(x,derivative=True)-rhs(rho,segment.evaluate(x))).reshape(2,2)
        assert np.linalg.norm(r*np.array([.4,.5])) < upper+2e-11
    assert float(restored_upper(v2['weighted_row_residual_integral'])) == pytest.approx(
        2*float(restored_upper(v['weighted_row_residual_integral'])),rel=1e-12)
    assert upper < 1e-8


def test_bad_weights_are_not_a_zero_error():
    s = TrajectorySegment(1.01,1.,np.ones(2),np.ones(2),np.zeros((6,2)))
    with pytest.raises(ValueError):reference_segment_defect(s,[1.],[0.],0.,0.)
