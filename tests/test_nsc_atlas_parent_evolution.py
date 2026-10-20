"""Independent controls for the transported two-wall averaging domain."""
import math

import numpy as np
import pytest
from scipy.optimize import brentq

from recursive_horizons.nsc_atlas_parent_evolution import (
    NormalParent, ParentBudget, normal_bulk, normal_moments,
)
from recursive_horizons.nsc_atlas_nested_slice import nested_slice_record
from recursive_horizons.nsc_atlas_geometric_average import KAPPA, averaged_gravity_balance


@pytest.fixture(scope="module")
def parent():
    return NormalParent(64)


@pytest.mark.parametrize("sign", [-1., 1.])
def test_schwarzschild_normal_against_parametric_freefall(sign):
    M, R0, t = .4, 1.2, .03
    E = sign*math.sqrt(1-2*M/R0)
    ray = normal_bulk([[R0, 0., 1., 0., E, M/R0**2]], t, M, 0.)[0]

    def radius(r):
        eta = brentq(lambda x: math.sqrt(r**3/(8*M))*(x+math.sin(x))-t, 0., .5)
        return r*(1+math.cos(eta))/2

    d = 1e-5
    assert ray[0] == pytest.approx(radius(R0), abs=2e-13)
    assert ray[2] == pytest.approx((radius(R0+d)-radius(R0-d))/(2*d), abs=3e-11)
    assert ray[4] == E


def test_exact_de_sitter_volume_and_raychaudhuri():
    t = .01
    rays = normal_bulk([[.8, .8, 1., 1., 1., 0.]], t, 0., 1.)
    assert rays[0, :4] == pytest.approx([.8*math.exp(t)]*2+[math.exp(t)]*2, abs=3e-14)
    m, c = normal_moments(rays, 0., 1.)
    assert c == pytest.approx([0.], abs=3e-14)
    assert m[0, 1]/m[0, 0] == pytest.approx(3.)


def test_bridge_zero_signed_energy_is_regular(parent):
    initial = [[parent.rb, 0., 1., 0., 0., .01/parent.rb**2-parent.rb]]
    ray = normal_bulk(initial, .01, .01, 1.)[0]
    assert parent.rb/2 < ray[0] < parent.rb
    assert ray[2] > 1
    _, c = normal_moments([ray], .01, 1.)
    assert abs(c[0]) < 2e-8


@pytest.mark.parametrize("name", ["inner", "outer"])
def test_crossed_spatial_metric_from_neighbouring_geodesics(parent, name):
    # This differentiates the physical map at fixed final time. It does not
    # use the Jacobi equation or the Hamiltonian to define its reference.
    s, t, h = .004, .009, 2e-6
    c = parent.crossing(name, s)
    final = normal_bulk([c["outgoing"]], t-s, .01, 1.)[0]
    radii, energies = [], []
    for si in (s-h, s+h):
        ci = parent.crossing(name, si)
        r = normal_bulk([ci["outgoing"]], t-si, .01, 1.)[0]
        radii.append(r[0])
        energies.append(r[4])
    reference = (radii[1]-radii[0])/(2*h*c["label_rate"])
    assert reference == pytest.approx(final[4]*final[2], abs=2e-8)
    El = (energies[1]-energies[0])/(2*h*c["label_rate"])
    assert El == pytest.approx(final[5], abs=2e-7)
    _, constraint = normal_moments([final], .01, 1.)
    assert abs(constraint[0]) < 2e-10


def test_initial_slice_recovery(parent):
    original = nested_slice_record(include_partition=True)
    # Recover the independently integrated historical source quantities.
    # Locate the existing projection without relying on a rewritten record.
    rows = original["partition"]["bulk_regions_in_ell_units"]
    moments = np.sum(list(rows.values()), axis=0)
    for s in original["partition"]["surfaces_in_ell_units"].values():
        moments[4:] += [s["curvature"], s["energy"], s["stress_trace"]]
    current = parent.snapshot(0.)
    assert current["moments"] == pytest.approx(moments, rel=2e-12, abs=2e-10)


def test_actual_evolved_balance_and_volume_derivatives(parent):
    t, h = .005, 5e-5
    mid = parent.snapshot(t)
    lo, hi = parent.snapshot(t-h), parent.snapshot(t+h)
    a = mid["sources"]
    d = lambda k: (hi["sources"][k]-lo["sources"][k])/(2*h)
    balance = averaged_gravity_balance(H=a["H_normal"], Q=a["Q"],
        curvature=a["mean_intrinsic_curvature"], Q_rate=d("Q"),
        curvature_rate=d("mean_intrinsic_curvature"),
        mean_energy=a["kappa_mean_energy"]/KAPPA,
        mean_stress_trace=3*a["kappa_mean_pressure"]/KAPPA,
        mean_energy_rate=d("kappa_mean_energy")/KAPPA)
    assert abs(balance["integrability_residual"]) < 3e-9
    assert abs(d("H_normal")-(a["three_volume_acceleration"]/3-a["H_normal"]**2)) < 2e-9
    assert abs((hi["moments"][0]-lo["moments"][0])/(2*h)-mid["moments"][1]) < 6e-7
    assert balance["Lambda_trace_rate"] > .18
    assert balance["explicit_source_balance"] > .04


def test_endpoint_is_evolved_nonvacuum_geometry(parent):
    r = parent.snapshot(.01)
    assert r["sources"]["geometric_trace_projection"] == pytest.approx(-.01346818189, abs=2e-9)
    assert r["sources"]["geometric_vacuum_condition"] > 2.
    assert r["controls"]["min_Jacobi_scale"] > 0
    assert r["controls"]["max_bulk_constraint"] < 2e-8
    assert r["walls"]["inner"]["proper_time"] != r["walls"]["outer"]["proper_time"]
    assert r["walls"]["inner"]["v"] < 0 < r["walls"]["outer"]["v"]
    assert r["pieces"]["inner_crossed"][0] > 0
    assert r["pieces"]["outer_crossed"][0] > 0


def test_explicit_domain_and_budget(parent):
    with pytest.raises(ValueError):
        parent.snapshot(.04)
    with pytest.raises(ValueError):
        NormalParent(128)
    with pytest.raises(RuntimeError):
        normal_bulk([[1., 0., 1., 0., 1., 0.]], .01, 0., 0., budget=ParentBudget(limit=1))


def test_received_frequency_is_derivative_of_same_wall_clock(parent):
    t, h = .005, 1e-6
    lo, hi = parent.central_flash(t-h), parent.central_flash(t+h)
    clock_ratio = ((hi["central_reception_time"]-lo["central_reception_time"])
                   /(hi["emission_wall_proper_time"]-lo["emission_wall_proper_time"]))
    assert clock_ratio == pytest.approx(parent.central_flash(t)["one_plus_redshift"], abs=3e-9)
    assert parent.central_flash(.01)["central_reception_time"] > parent.endpoint
