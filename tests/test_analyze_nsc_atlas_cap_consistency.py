"""Independent manufactured cap/curve controls, without retained replay."""
import copy
import time

import numpy as np
import pytest

import analyze_nsc_atlas_cap_consistency as tool
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget

from control_budget import ControlBudget

_CONTROL_BUDGET = ControlBudget(maximum_evaluations=1000000, cpu_seconds=30.)


@pytest.fixture
def budget():
    # Pytest may register the imported fixture separately in the recorder
    # module. Both registrations intentionally share this one accumulator.
    with _CONTROL_BUDGET.active() as shared:
        yield shared


def manufactured_curve():
    """Exact jets of a moving curve in an explicitly time-dependent metric."""
    alpha = .07
    R, T = .46, .12
    Rx, Tx, Rxx, Txx = .17, 1.3, -.09, .11
    a, aT = .2+.3*T+.1*T*T, .3+.2*T
    n, nr, nT = a*R**3, 3*a*R**2, aT*R**3
    A = 1-R*R-2*alpha*n/R
    Ar, AT = -2*R-2*alpha*(nr/R-n/R**2), -2*alpha*nT/R
    delta, dr, dT = .05*T+.08*R*R, .16*R, .05
    Ax, dx = AT*Tx+Ar*Rx, dT*Tx+dr*Rx
    sigma2 = A*np.exp(2*delta)*Tx**2-Rx**2/A
    sigma = np.sqrt(sigma2)
    sigmax = (np.exp(2*delta)*((Ax+2*A*dx)*Tx**2+2*A*Tx*Txx)
              -2*Rx*Rxx/A+Rx**2*Ax/A**2)/(2*sigma)
    u = np.array([Tx, Rx])/sigma
    udot = np.array([Txx, Rxx])/sigma**2-np.array([Tx, Rx])*sigmax/sigma**3
    # Christoffels from g_TT=-A exp(2 delta), g_rr=1/A, computed independently
    # of the compact extrinsic-curvature formula used by the diagnostic.
    gamma_T = np.array([[AT/(2*A)+dT, Ar/(2*A)+dr],
                        [Ar/(2*A)+dr, -AT/(2*A**3*np.exp(2*delta))]])
    gamma_r = np.array([[A*np.exp(2*delta)*(Ar/2+A*dr), -AT/(2*A)],
                        [-AT/(2*A), -Ar/(2*A)]])
    accel = udot+np.array([u@gamma_T@u, u@gamma_r@u])
    normal = np.array([-np.exp(delta)*Rx/sigma, np.exp(delta)*Tx/sigma])
    expected_K = normal@accel
    wall = dict(R=R, T=T, A=A, delta=delta, n=n, f=.2, Pi=.4, Phi=-.1)
    first = dict(R=Rx, T=Tx, s=sigma+.003, n=nT*Tx+nr*Rx,
                 A=Ax, delta=dx, f=.23, chi=.23*R+.2*Rx)
    second = dict(R=Rxx, T=Txx)
    radial = dict(n=nr, delta=dr)
    return wall, first, second, radial, alpha, nT, expected_K


def test_actual_curve_acceleration_against_christoffels(budget):
    wall, first, second, radial, alpha, nT, expected_K = manufactured_curve()
    before = copy.deepcopy((wall, first, second, radial))
    row = tool.curve_quantities(wall, first, second, radial, alpha=alpha, budget=budget)
    np.testing.assert_allclose(row["n_T"], nT, rtol=2e-14)
    np.testing.assert_allclose(row["K_geom"], expected_K, rtol=2e-14, atol=2e-15)
    assert abs(row["proper_clock_derivative_defect"]) == pytest.approx(.003)
    assert before == (wall, first, second, radial)


def test_curve_uses_measured_metric_derivative_and_rejects_spacelike(budget):
    wall, first, second, radial, alpha, _, _ = manufactured_curve()
    reference = tool.curve_quantities(wall, first, second, radial, alpha=alpha, budget=budget)
    changed = dict(first, A=first["A"]+.04)
    got = tool.curve_quantities(wall, changed, second, radial, alpha=alpha, budget=budget)
    assert got["sigma_x"] != reference["sigma_x"]
    assert got["n_T"] == reference["n_T"]
    with pytest.raises(ValueError):
        tool.curve_quantities(wall, dict(first, T=-1.), second, radial, alpha=alpha, budget=budget)
    with pytest.raises(ValueError):
        tool.curve_quantities(wall, dict(first, R=10.), second, radial, alpha=alpha, budget=budget)


def test_temporal_support_and_polynomial_stencils(budget):
    expected_neighbours = set(range(61, 71)) | {122, 124, 125, 126, 127, 128}
    used = set()
    dx = .001
    for cut in (63, 68, 128):
        for stride in (1, 2):
            for derivative in (1, 2):
                indices, weights = tool.temporal_stencil(cut, stride, derivative)
                indices = np.asarray(indices)
                used.update(indices.tolist())
                y = (indices-cut)*dx
                # Second derivatives are exact on cubic polynomials for both
                # centered-three and terminal-backward-four formulas.
                polynomial = .3+.2*y+.7*y*y
                if derivative == 2:
                    polynomial += .11*y**3
                got = tool.differentiate_samples(polynomial, indices, cut, stride, derivative, dx, budget)
                assert got == pytest.approx(.2 if derivative == 1 else 1.4, abs=2e-9)
                assert len(indices) == (4 if cut == 128 and derivative == 2 else 3)
    assert used == expected_neighbours
    for cut, stride, derivative in ((62, 1, 1), (63, 3, 1), (63, 1, 3)):
        with pytest.raises(ValueError):
            tool.temporal_stencil(cut, stride, derivative)


def test_radial_backward_quadratic_and_analytic_vacuum(budget):
    r = np.linspace(0, .52, 65)
    for stride in (1, 2):
        assert tool.radial_first(.3+.1*r+.4*r*r, r[1], stride, budget) == pytest.approx(.1+.8*r[-1], abs=2e-14)
    wall = dict(R=.52, A=1-.52**2, delta=0., n=.13, f=.2, Pi=.4, Phi=-.1)
    for alpha in (0., .03):
        row = tool.residual_quantities(wall, n_T=.17, n_r=.2, delta_r=.01, alpha=alpha, budget=budget)
        assert row["m_r"] == pytest.approx(1.5*.52**2+alpha*.2)
        assert row["E_m"] == pytest.approx(alpha*row["e_m"])
        assert row["E_r"] == pytest.approx(alpha*row["e_r"])
        expected_em = .17-.52**2*wall["A"]**2*.4*(-.1)/2
        assert row["e_m"] == pytest.approx(expected_em)
        if alpha == 0:
            assert row["E_m"] == row["E_r"] == 0.
            assert row["e_m"] != 0.


def test_stable_potential_quotient_and_lapse_residual(budget):
    R, A, delta, f, Pi, Phi = .52, .71, .03, .7, .4, -.1
    wall = dict(R=R, A=A, delta=delta, n=.13, f=f, Pi=Pi, Phi=Phi)
    for alpha in (0., 1e-20, .3):
        pdiff = f*f if alpha < 1e-15 else 3*np.expm1(4*np.log(np.cosh(np.sqrt(alpha/6)*f)))/alpha
        expected_er = .2-R*R*pdiff/2-R*R*A*(Pi*Pi+Phi*Phi)/4
        row = tool.residual_quantities(wall, n_T=.17, n_r=.2, delta_r=.01, alpha=alpha, budget=budget)
        assert row["e_r"] == pytest.approx(expected_er, abs=2e-15)
        assert row["E_delta"] == pytest.approx(.01-alpha*R*(Pi*Pi+Phi*Phi)/2)


def test_static_de_sitter_curve_has_nonzero_geometric_acceleration(budget):
    R = .52
    A = 1-R*R
    wall = dict(R=R, A=A, delta=0., n=0., f=0., Pi=0., Phi=0.)
    first = dict(R=0., T=1., s=np.sqrt(A), n=0., A=0., delta=0., f=0., chi=0.)
    row = tool.curve_quantities(wall, first, dict(R=0., T=0.), dict(n=0., delta=0.), alpha=0., budget=budget)
    assert row["K_geom"] == pytest.approx(-R/np.sqrt(A), rel=2e-15)
    assert row["proper_clock_derivative_defect"] == 0.


def test_helper_budget_exhaustion_debits_common_control_allowance(budget):
    class Child(Budget):
        def check(self, count=0):
            budget.check(count)
            super().check(count)
    child = Child(maximum_evaluations=1, cpu_seconds=30., started=budget.started)
    before = budget.evaluations
    with pytest.raises(RuntimeError, match="exhausted"):
        tool.radial_first(np.linspace(0, 1, 65), 1/64, 1, child)
    assert budget.evaluations-before == child.evaluations > 1


def test_distinct_incoming_traces_change_flux_without_replacing_left_residual(budget):
    wall, first, second, radial, alpha, _, _ = manufactured_curve()
    wall["alpha"] = alpha
    wm, wp = wall["Pi"]-wall["Phi"], wall["Pi"]+wall["Phi"]
    # Arrange an independently chosen kinematic mismatch, and a different
    # hybrid mismatch. All three traces must remain visible in the result.
    wave_speed = wall["A"]*np.exp(wall["delta"])
    kinematic_shift, imposed_shift = .06, -.04
    first["f"] = ((wave_speed*first["T"]+first["R"])*(wp+kinematic_shift)
                   +(wave_speed*first["T"]-first["R"])*wm)/2
    first["chi"] = wall["R"]*first["f"]+first["R"]*wall["f"]
    curve = tool.curve_quantities(wall, first, second, radial, alpha=alpha, budget=budget)
    hybrid = dict(w_plus=wp+imposed_shift, w_minus=wm, k_minus=.2,
                  r_s=curve["v"]+.02, T_s=curve["T_s"]-.01, Vdot_plus=.8)
    before = tool.residual_quantities(wall, n_T=curve["n_T"], n_r=radial["n"],
                                     delta_r=radial["delta"], alpha=alpha, budget=budget)
    trace = tool.trace_quantities(wall, curve, hybrid, budget)
    after = tool.residual_quantities(wall, n_T=curve["n_T"], n_r=radial["n"],
                                    delta_r=radial["delta"], alpha=alpha, budget=budget)
    assert before == after
    assert trace["kinematic_minus_left"] == pytest.approx(kinematic_shift)
    assert trace["hybrid_minus_left"] == pytest.approx(imposed_shift)
    assert trace["hybrid_minus_kinematic"] == pytest.approx(imposed_shift-kinematic_shift)
    expected = alpha*wall["R"]**2*wall["A"]**2*np.exp(wall["delta"])*(2*wp*imposed_shift+imposed_shift**2)/8
    assert trace["actual_flux_hybrid_minus_left"] == pytest.approx(expected)
    assert trace["normalized_flux_left"] == before["normalized_mass_flux"]
    assert trace["outgoing_identity_defect"] == 0.
    assert trace["primitive_product_derivative_defect"] == pytest.approx(0., abs=1e-16)
    assert trace["v_curve_minus_hybrid"] == pytest.approx(-.02)
    assert trace["T_s_curve_minus_hybrid"] == pytest.approx(.01)
