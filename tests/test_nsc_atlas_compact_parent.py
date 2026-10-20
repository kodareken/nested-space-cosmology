"""Independent gluing, volume and causal controls for compact completion."""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_compact_parent import (
    CompactParent, reflected_initial_wall, quiet_throat_continuation,
)
from recursive_horizons.nsc_atlas_geometric_average import (
    KAPPA, _bridge_primitive, geometric_slice_average, averaged_gravity_balance,
)
from recursive_horizons.nsc_atlas_parent_evolution import normal_bulk
from recursive_horizons.nsc_atlas_nested_slice import bridge_chart


@pytest.fixture(scope="module")
def compact():
    return CompactParent(64)


def test_right_orientation_preserves_positive_membrane_and_future(compact):
    w = compact.regular.initial_walls["outer"]
    r = reflected_initial_wall(w)
    for key, value in r.items():
        if key.endswith("_residual"):
            assert abs(value) < 3e-14
    for side in ("minus", "plus"):
        radial_rate = r["gamma"]*r["V_"+side]+r["b"]*r["Rprime_"+side]
        beta = r["gamma"]*r["Rprime_"+side]+r["b"]*r["V_"+side]
        assert radial_rate == pytest.approx(w["Rdot"], abs=3e-14)
        assert beta == pytest.approx(r["beta_"+side], abs=3e-14)
    assert r["Rprime_minus"] > 0 > r["Rprime_plus"]
    assert r["gamma"] > 1 and r["b"] < 0


def test_bridge_open_overlap_reflects_smoothly_at_finite_time():
    y = np.array([-.1, 0., .1])
    state = bridge_chart(y, .8, 0., .4)
    rays = np.column_stack((state["radius"], state["velocity"], np.ones(3),
                           np.zeros(3), state["Rprime"], state["Rsecond"]))
    evolved = normal_bulk(rays, .01, .4, 0.)
    for column in (0, 1, 2, 3, 5):
        assert evolved[0, column] == pytest.approx(evolved[2, column], abs=2e-14)
    assert evolved[0, 4] == pytest.approx(-evolved[2, 4], abs=2e-14)
    assert evolved[1, 4] == 0 and evolved[1, 2] > 0
    # R_x/A tends to zero from both sides: no curvature delta at the join.
    assert state["Rprime"][1] == 0


def test_regular_side_matches_independent_original_body_and_bridge_primitive(compact):
    old = geometric_slice_average()
    q = old["integrals_in_ell_units"]
    positive_end = 64*math.pi*.4**3*_bridge_primitive(math.acosh(math.sqrt(4/.8)))
    reference = [old["volume_over_ell3"]-positive_end, q["theta"], q["theta_squared"],
                 q["shear_squared"], q["curvature_total"],
                 q["energy_cap"]+q["energy_wall"],
                 q["stress_trace_cap"]+q["stress_trace_wall"]]
    assert compact.side("regular", 0.)["moments"] == pytest.approx(reference, abs=2e-10)


def test_nested_side_is_the_original_body_without_its_ae_end(compact):
    old = compact.nested.snapshot(.005)
    reference = np.asarray(old["moments"])-old["pieces"]["outer_bridge_right"]
    assert compact.side("nested", .005)["moments"] == pytest.approx(reference, abs=2e-10)


def test_symmetric_double_has_same_means_and_twice_the_actual_sources(compact):
    half = compact.side("nested", .005)
    full = compact.snapshot(.005, right="nested")
    assert full["moments"] == pytest.approx(2*half["moments"], abs=2e-12)
    assert full["between_side_Q"] == pytest.approx(0., abs=1e-14)
    assert full["wall_count"] == 4
    assert full["side_volume_fractions"] == {"left": .5, "right": .5}


@pytest.mark.parametrize("t", [0., .005, .01])
@pytest.mark.parametrize("right", ["regular", "nested"])
def test_common_kruskal_clock_and_causal_quadrants(compact, t, right):
    left = compact.connector_event(t, side="left", right=right)
    other = compact.connector_event(t, side="right", right=right)
    assert left["U"] > 0 > left["V"]
    assert other["U"] < 0 < other["V"]
    for event in (left, other):
        assert event["U_tau"] > 0 and event["V_tau"] > 0
        assert abs(event["radius_identity"]) < 2e-13
        assert event["tangent_norm"] == pytest.approx(-1., abs=3e-14)
    # Even their inward null fronts have no common future event in the
    # physical Schwarzschild connector UV<1 for this emission pair.
    assert left["U"]*other["V"] > 1.


def test_nonzero_global_nesting_increment_with_identical_local_signal(compact):
    a = compact.snapshot(.01, right="regular")
    b = compact.snapshot(.01, right="nested")
    assert a["left_flash"] == b["left_flash"]
    assert a["left_walls"] == b["left_walls"]
    assert b["sources"]["total_trace_projection"]-a["sources"]["total_trace_projection"] > .04
    for row in (a, b):
        assert .48 < row["curved_region_volume_fraction"] < .53
        assert row["sources"]["geometric_vacuum_condition"] > 6.
        decomposition = (row["supplied_bulk_trace"]+row["material_wall_trace"]
                         +row["sources"]["geometric_trace_projection"])
        assert decomposition == pytest.approx(row["sources"]["total_trace_projection"], abs=3e-14)
        assert abs(row["sources"]["hamiltonian_residual"]) < 4e-14
        assert row["controls"]["max_bulk_constraint"] < 2e-8


def test_actual_compact_volume_evolution_and_conservation(compact):
    h = 5e-5
    for right in ("regular", "nested"):
        mid = compact.snapshot(.005, right=right)
        lo, hi = (compact.snapshot(.005+sign*h, right=right) for sign in (-1, 1))
        a = mid["sources"]
        rate = lambda k: (hi["sources"][k]-lo["sources"][k])/(2*h)
        balance = averaged_gravity_balance(H=a["H_normal"], Q=a["Q"],
            curvature=a["mean_intrinsic_curvature"], Q_rate=rate("Q"),
            curvature_rate=rate("mean_intrinsic_curvature"),
            mean_energy=a["kappa_mean_energy"]/KAPPA,
            mean_stress_trace=3*a["kappa_mean_pressure"]/KAPPA,
            mean_energy_rate=rate("kappa_mean_energy")/KAPPA)
        assert abs(balance["integrability_residual"]) < 2e-7
        assert abs(rate("H_normal")-(a["three_volume_acceleration"]/3-a["H_normal"]**2)) < 3e-8
        assert abs((hi["total_volume"]-lo["total_volume"])/(2*h)-mid["moments"][1]) < 2e-6


def test_named_completions_and_finite_time_only(compact):
    with pytest.raises(ValueError):
        compact.snapshot(.005, right="fitted_fraction")
    with pytest.raises(ValueError):
        compact.snapshot(.04)


def test_quiet_throat_lifetime_and_pressure_threshold():
    flat = quiet_throat_continuation(.4, 0.)
    assert flat["normal_lifetime"] == pytest.approx(math.pi*.4, abs=2e-15)
    host = quiet_throat_continuation(.01, 1.)
    assert host["analytic_lower_bound"] < host["normal_lifetime"] < host["analytic_upper_bound"]
    assert .0314 < host["normal_lifetime"] < .0315
    assert host["supplied_kappa_radial_pressure"] > host["kappa_radial_pressure_bounce_threshold"]


def test_general_spherical_radial_identity_from_warped_einstein_tensor():
    import sympy as sp
    R, m, Hnn, Hss, kappa, rho, pressure = sp.symbols("R m Hnn Hss kappa rho pressure")
    # Warped-product G_ab = -2 H_ab/R + [2 box(R)/R-2m/R^3] g_ab.
    common = 2*(-Hnn+Hss)/R-2*m/R**3
    Gnn, Gss = -2*Hnn/R-common, -2*Hss/R+common
    solved = sp.solve([Gnn-kappa*rho, Gss-kappa*pressure], (Hnn, Hss))
    assert sp.simplify(solved[Hnn]+m/R**2+kappa*R*pressure/2) == 0
    assert sp.simplify(solved[Hnn]+solved[Hss]+kappa*R*(rho+pressure)/2) == 0
