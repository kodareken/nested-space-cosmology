"""Fixed-label coefficients: exact manufactured families, no physical solve."""
from types import SimpleNamespace

import numpy as np
import pytest
import sympy as sp

from control_budget import ControlBudget
from recursive_horizons import nsc_atlas_cap_response as model


CONTROLS = ControlBudget(maximum_evaluations=20000, cpu_seconds=20.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope="module", autouse=True)
def report():
    yield
    print("wall coefficient pullback controls:", CONTROLS.budget.evaluations,
          "charged values;", CONTROLS.active_cpu_seconds, "active CPU seconds")


def test_source_ode_variation_recovers_the_original_affine_constraint(budget):
    budget.check(6000)
    x = sp.symbols("x", positive=True)
    p = sp.Function("p")(x)
    g = x*p
    # Differentiate the exact radius ODE at alpha=0, rho=x, rho'=1.
    varied = (g*sp.diff(g, x, 2)-2*sp.diff(g, x)**2)/(6*x)
    varied += g*sp.diff(g, x)/(3*x*x)
    original = x*p*sp.diff(p, x, 2)/6-x*sp.diff(p, x)**2/3
    assert sp.simplify(varied-original) == 0
    # The apparent reduced-field terms cancel only with the second source term.
    omitted = varied-g*sp.diff(g, x)/(3*x*x)
    assert sp.simplify(omitted-original) != 0


def test_inverse_event_and_conformal_frame_match_an_exact_family(budget):
    budget.check(7000)
    alpha, s = sp.symbols("alpha s", real=True)
    s0, x0, v = sp.Rational(1, 10), sp.Rational(1, 2), sp.Rational(1, 5)
    beta = sp.Rational(6, 5)
    radius = x0+v*(s-s0)
    rho, scalar = sp.Rational(3, 10)+s/10, sp.Rational(2, 5)+s/5
    clock0 = beta/(1-x0*x0)*s+sp.Rational(7, 10)*(s-s0)**2
    clock1 = sp.Rational(2, 5)+3*s/10
    # V(alpha,s)=s+alpha*(1+2s) has an exact inverse at V=s0.
    event = (s0-alpha)/(1+2*alpha)
    assert sp.simplify((s+alpha*(1+2*s)).subs(s, event)-s0) == 0
    exact_radius = ((radius+alpha*rho)*sp.sqrt(1-alpha*scalar**2/6)).subs(s, event)
    exact_clock = (clock0+alpha*clock1).subs(s, event)
    e = (rho-radius*scalar**2/12).subs(s, s0)
    out = model.fixed_label_wall_coefficients(
        x=float(x0), v=float(v), beta_minus=float(beta), e=float(e),
        T1=float(clock1.subs(s, s0)), V1=float(1+2*s0), V0_s=1.)
    assert out["proper_time_shift"] == pytest.approx(float(sp.diff(event, alpha).subs(alpha, 0)))
    assert out["e_fixed_V"] == pytest.approx(
        float(sp.diff(exact_radius, alpha).subs(alpha, 0)), abs=2e-16)
    assert out["T1_fixed_V"] == pytest.approx(
        float(sp.diff(exact_clock, alpha).subs(alpha, 0)), abs=2e-15)
    assert out["e_fixed_V"] != pytest.approx(float(e))
    assert out["T1_fixed_V"] != pytest.approx(float(clock1.subs(s, s0)))


def test_label_rescaling_and_quiet_limits(budget):
    budget.check(100)
    data = dict(x=.5, v=.2, beta_minus=1.2, e=.3, T1=.7, V1=.4, V0_s=1.2)
    expected = model.fixed_label_wall_coefficients(**data)
    scaled = dict(data, V1=-3*data["V1"], V0_s=-3*data["V0_s"])
    assert model.fixed_label_wall_coefficients(**scaled) == pytest.approx(expected)
    quiet = model.fixed_label_wall_coefficients(**dict(data, e=0., T1=0., V1=0.))
    assert set(quiet) == {"proper_time_shift", "e_fixed_V", "T1_fixed_V"}
    assert all(value == 0 for value in quiet.values())
    assert data == dict(x=.5, v=.2, beta_minus=1.2, e=.3, T1=.7, V1=.4, V0_s=1.2)


@pytest.mark.parametrize("changed", [
    {"x": .44}, {"x": .56}, {"V0_s": 0.}, {"T1": np.nan},
    {"V1": np.inf}, {"T1": True}, {"v": [1.]},
    {"V1": 1e308, "V0_s": 1e-308},
])
def test_undefined_pullbacks_are_rejected(changed, budget):
    budget.check(10)
    data = dict(x=.5, v=.2, beta_minus=1.2, e=.3, T1=.7, V1=.4, V0_s=1.2)
    data.update(changed)
    with pytest.raises(ValueError):
        model.fixed_label_wall_coefficients(**data)


def manufactured_solution(monkeypatch):
    """Only polynomial clock data and supplied rows; no field queries."""
    calls = {"checks": 0, "cap": 0, "wall": 0, "ray": 0}

    def check():
        calls["checks"] += 1

    def cap_row(t):
        calls["cap"] += 1
        return dict(t=t, x=.5, e=.3, V0_s=1.2, q0=.8, q0_s=2.)

    def wall_row(t):
        calls["wall"] += 1
        return dict(V1=.4, v=.2, beta_minus=1.2)

    def forbidden_ray(*args, **kwargs):
        calls["ray"] += 1
        raise AssertionError("wall-only query evaluated a central-arrival ray")

    solution = object.__new__(model.CapSolution)
    solution.budget = SimpleNamespace(check=check)
    solution.settings = dict(time_panels=4, radial_order=4, ray_order=4)
    solution._wall = cap_row
    solution.wall_solution = SimpleNamespace(sample=wall_row)
    solution.cap_field = object()
    solution.ray_cache, solution.ray_nodes = {}, 0
    solution.maxima = dict(fixed_V_conversion_defect=0.)
    nodes, solution.time_basis = model._gauss_antiderivatives()
    solution.time_cuts = np.linspace(0, 1, 5)
    times = solution.time_cuts[:-1, None]+(1+nodes)/8
    solution.time_values = 2+3*times
    solution.time_prefix = 2*solution.time_cuts+1.5*solution.time_cuts**2
    monkeypatch.setattr(model, "ray_integral", forbidden_ray)
    return solution, calls


def test_wall_query_uses_known_clock_integral_without_a_ray(monkeypatch, budget):
    budget.check(100)
    solution, calls = manufactured_solution(monkeypatch)
    t = .37  # A partial panel, not a precomputed cut.
    out = solution.wall_sample(t)
    expected_clock = 2*t+1.5*t*t
    assert out["T1"] == pytest.approx(expected_clock, abs=2e-15)
    assert out["e_fixed_V"] == pytest.approx(.3-.2/3, abs=2e-16)
    assert out["T1_fixed_V"] == pytest.approx(expected_clock-1.6/3, abs=2e-15)
    assert calls == dict(checks=1, cap=1, wall=1, ray=0)
    assert solution.ray_cache == {} and solution.ray_nodes == 0
    assert "Delta_C_fixed_V" not in out


def test_historical_sample_still_uses_the_arrival_shift(monkeypatch, budget):
    budget.check(100)
    solution, calls = manufactured_solution(monkeypatch)

    def ray(*args, **kwargs):
        calls["ray"] += 1
        return dict(integral=.125, positive_integral=.2, negative_integral=-.075,
                    integrand_sample_min=-.1, integrand_sample_max=.3, ray_nodes=4)

    monkeypatch.setattr(model, "ray_integral", ray)
    out = solution.sample(.37)
    expected_clock = 2*.37+1.5*.37**2
    expected_arrival = expected_clock+.3/.75+.125-2*.4/1.2
    assert out["T1"] == pytest.approx(expected_clock, abs=2e-15)
    assert out["Delta_C_fixed_V"] == pytest.approx(expected_arrival, abs=2e-15)
    assert out["exterior_label_shift"] == pytest.approx(2*.4/1.2)
    assert calls == dict(checks=1, cap=1, wall=1, ray=1)
    assert solution.ray_nodes == 4
    assert not {"proper_time_shift", "e_fixed_V", "T1_fixed_V"}.intersection(out)
