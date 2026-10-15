"""Independent manufactured moving-wall equations; no physical preparation."""
from types import SimpleNamespace
import time

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_bulk_response import BulkBudget
from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget
from recursive_horizons.nsc_atlas_scalar_evolution import solve_scalar_correction


@pytest.fixture(autouse=True)
def analytic_cpu_budget():
    started = time.process_time()
    yield
    assert time.process_time()-started < 30


class ManufacturedInputs:
    def __init__(self, mode="forced"):
        self.domain = dict(U0=0., Ua=.4, Ub=1., Va=2., Vb=2.7)
        external = PreparedFieldBudget()
        self.budget = BulkBudget(started=external.started, external=external)
        self.training_times = np.array([0., .25, .5, .75, 1.])
        self.geometry = {"wall": self.wall}
        self.mode, self.queries = mode, []
        self.lift = SimpleNamespace(domain=self.domain, budget=self.budget,
            knot_times=np.array([0., .25, .75, 1.]), B=self.B, b_at_U=self.b)

    def wall(self, V, derivative=0):
        t = (np.asarray(V)-2.)/.7
        return .6*(.7+.6*t)/.7 if derivative == 1 else .4+.6*(.7*t+.3*t*t)

    def B(self, U):
        return np.maximum(np.asarray(U)-.4, 0.)**4 if self.mode != "zero" else np.zeros_like(np.asarray(U))

    def b(self, U):
        return 4*np.maximum(np.asarray(U)-.4, 0.)**3 if self.mode != "zero" else np.zeros_like(np.asarray(U))

    def source(self, t):
        return (t**3, 3*t*t/.7) if self.mode != "zero" else (0., 0.)

    def kernel(self, U, V):
        return np.zeros_like(U) if self.mode in ("zero", "lift") else .2+np.asarray(U)+.1*V

    def exact(self, U, t):
        V = 2.+.7*t
        W, WV = self.wall(V), self.wall(V, 1)
        q, qV = t*t, 2*t/.7
        # z=q*(2U/W-U^2/W^2) has z(0,V)=0 and z_U(W,V)=0.
        z = q*(2*U/W-U*U/W**2)
        zUV = qV*(2/W-2*U/W**2)+q*(-2*WV/W**2+4*U*WV/W**3)
        return z, zUV

    def eval(self, U, t):
        U = np.asarray(U, dtype=float)
        V = 2.+.7*t
        W = float(self.wall(V))
        assert np.all(U <= W+4e-15), "a coefficient was queried in the excised cap"
        self.queries.append((t, U.copy()))
        g1, g1V = self.source(t)
        if self.mode in ("zero", "lift"):
            F = G = np.zeros_like(U)
            F0 = 0.
        else:
            a = .3+t*t
            F, F0 = a*(W-U), a*W
            z, zUV = self.exact(U, t)
            # Independent strong-form identity: G=F_U-z_UV-K*y.
            G = -a-zUV-self.kernel(U, V)*(z+g1+self.B(U))
        return dict(U=U, F=F, G=G, Fsource=F0, g1=g1, g1_V=g1V)


@pytest.mark.parametrize("mode", ("zero", "lift"))
def test_zero_kernel_quiet_and_nonconstant_source_lift_controls(mode):
    data = ManufacturedInputs(mode)
    grid = np.linspace(0., 1., 17)
    solution = solve_scalar_correction(data, data.lift, grid, kernel=data.kernel)
    for t in np.linspace(0., 1., 13):
        row = solution.sample(t)
        W = data.wall(2.+.7*t)
        g, gp = data.source(t)
        assert row["z_wall"] == 0
        assert row["y_wall"] == pytest.approx(g+data.B(W), abs=2e-15)
        assert row["y_wall_V"] == pytest.approx(gp+data.wall(2.+.7*t, 1)*data.b(W), abs=3e-15)
        assert row["source_defect"] == row["source_rate_defect"] == 0
        assert row["inactive_state_spread"] == 0
    initial = solution.bulk(np.linspace(0., .4, 11), 0.)
    assert np.array_equal(initial["y"], np.zeros(11))
    assert solution.kernel_values == data.budget.external.kernel_values


def test_nonzero_flux_density_kernel_and_nonlinear_moving_wall_exact_solution():
    data = ManufacturedInputs()
    grid = np.linspace(0., 1., 19)
    solution = solve_scalar_correction(data, data.lift, grid, kernel=data.kernel,
        tolerance=3e-11, atol=2e-13)
    for t in np.linspace(0., 1., 15):
        row = solution.sample(t)
        W = data.wall(2.+.7*t)
        g, gp = data.source(t)
        assert row["z_wall"] == pytest.approx(t*t, abs=3e-11)
        assert row["z_wall_V"] == pytest.approx(2*t/.7, abs=1e-10)
        assert row["y_wall"] == pytest.approx(t*t+g+data.B(W), abs=3e-11)
        assert row["y_wall_V"] == pytest.approx(2*t/.7+gp+data.wall(2.+.7*t, 1)*data.b(W), abs=1e-10)
        assert row["source_defect"] == row["source_rate_defect"] == 0
        assert abs(row["wall_flux"]) < 2e-15
        active = grid[grid <= W]
        result = solution.bulk(active, t)
        expected, _ = data.exact(active, t)
        assert result["z"] == pytest.approx(expected, abs=3e-11)
    assert all(any(abs(cut-edge) < 2e-15 for cut in solution.cuts) for edge in (0., .25, .5, .75, 1.))
    for birth in solution.wall_birth_times:
        assert min(abs(solution.cuts-birth)) < 2e-15
    assert solution.rhs_calls < 6000
    assert solution.summary()["maximum_absolute_source_defect"] == 0


def test_constant_kernel_prefix_matches_independent_semidiscrete_sum():
    # No physical reference preparation is invoked. Compare an independent
    # cell-by-cell constant-kernel integral with the vectorized clipped RHS.
    data = ManufacturedInputs("lift")
    solution = solve_scalar_correction(data, data.lift, np.linspace(0., 1., 11), kernel=data.kernel)
    t = .61
    V, W = 2.+.7*t, data.wall(2.+.7*t)
    U = np.minimum(solution.grid, W)
    state = np.linspace(.1, .9, len(U))
    solution.kernel = lambda U, V: np.full_like(U, .7)
    rate = solution._rate(t, state, counted=False)
    field = np.where(solution.grid < W, state, state[-1])+t**3+data.B(U)
    expected = np.zeros(len(U))
    for i in range(1, len(U)):
        expected[i] = -.7*.7*sum((field[j]+field[j-1])*(U[j]-U[j-1])/2 for j in range(1, i+1))
    assert rate == pytest.approx(expected, abs=2e-16)


def test_budget_rhs_cap_source_subtraction_and_incompatible_inputs_fail():
    data = ManufacturedInputs()
    with pytest.raises(RuntimeError, match="RHS budget"):
        solve_scalar_correction(data, data.lift, np.linspace(0., 1., 9), rhs_calls_max=1, kernel=data.kernel)
    bad = ManufacturedInputs()
    original = bad.eval
    def wrong_source(U, t):
        row = original(U, t)
        row["Fsource"] += 1.
        return row
    bad.eval = wrong_source
    with pytest.raises(ValueError, match="same cached"):
        solve_scalar_correction(bad, bad.lift, np.linspace(0., 1., 9), kernel=bad.kernel)
    bad = ManufacturedInputs()
    bad.lift.budget = SimpleNamespace()
    with pytest.raises(ValueError, match="share"):
        solve_scalar_correction(bad, bad.lift, kernel=bad.kernel)
    expired = ManufacturedInputs()
    expired.budget.started -= 121
    expired.budget.external.started = expired.budget.started
    with pytest.raises(RuntimeError, match="budget"):
        solve_scalar_correction(expired, expired.lift, kernel=expired.kernel)


@pytest.mark.parametrize("grid", ([0., 1., .5], [0., .9], np.linspace(0., 1., 258), [[0., 1.]], [0., np.nan, 1.]))
def test_invalid_grids_fail_before_coefficient_evaluations(grid):
    data = ManufacturedInputs()
    with pytest.raises(ValueError):
        solve_scalar_correction(data, data.lift, grid, kernel=data.kernel)
    assert not data.queries


def test_finite_readout_domain_and_nonzero_bottom_source_are_rejected():
    data = ManufacturedInputs("lift")
    solution = solve_scalar_correction(data, data.lift, np.linspace(0., 1., 9), kernel=data.kernel)
    for bad in (True, -.1, 1.1, np.nan):
        with pytest.raises(ValueError):
            solution.sample(bad)
    with pytest.raises(ValueError, match="retained exterior"):
        solution.bulk([.9], .1)
    bad = ManufacturedInputs()
    bad.source = lambda t: (1.+t, 1/.7)
    with pytest.raises(ValueError, match="quiet bottom"):
        solve_scalar_correction(bad, bad.lift, kernel=bad.kernel)
