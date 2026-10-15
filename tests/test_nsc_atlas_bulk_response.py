"""Field-independent checks of the exterior rectangle's metric propagation."""
import time

import numpy as np
import pytest
import sympy as sp

from recursive_horizons import nsc_atlas_bulk_response as model
from recursive_horizons.nsc_atlas_characteristic_initialization import source_geometry, constraint_samples


@pytest.fixture(scope="module")
def geometry():
    return source_geometry()


class ZeroField:
    def __init__(self, geometry):
        self.domain = dict(geometry)

    def bulk(self, U, t, order=32):
        return {key: np.zeros_like(U) for key in ("p", "p_U", "p_V")}


def test_mixed_equations_follow_from_independent_linearization():
    alpha, r, C, ru, rv, H, c, Hu, Hv, Huv, cuv, pu, pv = sp.symbols(
        "alpha r C ru rv H c Hu Hv Huv cuv pu pv")
    ruv = -ru*rv/r-C/(4*r)
    radial_rhs = -(ru+alpha*Hu)*(rv+alpha*Hv)/(r+alpha*H)-C*(1+alpha*c)/(4*(r+alpha*H))
    Huv_rhs = sp.diff(radial_rhs, alpha).subs(alpha, 0)
    lapse_rhs = 2*(ru+alpha*Hu)*(rv+alpha*Hv)/(r+alpha*H)**2+C*(1+alpha*c)/(2*(r+alpha*H)**2)-alpha*pu*pv
    cuv_rhs = sp.diff(lapse_rhs, alpha).subs(alpha, 0)
    yuv = ruv*H+ru*Hv+rv*Hu+r*Huv_rhs
    zuv = cuv_rhs+Huv_rhs/r-(Hu*rv+Hv*ru+H*ruv)/r**2+2*H*ru*rv/r**3
    y, z = r*H, c+H/r
    first, second = model.mixed_sources(r, C, y, z, pu, pv)
    assert sp.factor(yuv-first) == 0
    assert sp.factor(zuv-second) == 0
    assert sp.factor(zuv-(second+pu*pv)) != 0


def test_source_trace_uses_original_gauge_and_existing_constraint(geometry):
    saved = constraint_samples(geometry, samples=9, order=64)
    for index, V in enumerate(saved["V"]):
        src = model.source_traces(V, geometry)
        x, p = src["x"], src["p"]
        assert src["h_J"] == pytest.approx(saved["h"][index], abs=2e-12)
        assert src["h_J_x"] == pytest.approx(saved["h_prime"][index], abs=2e-12)
        assert src["H_E"] == pytest.approx(src["h_J"]-x*p*p/12)
        assert src["c_E"] == pytest.approx(-p*p/6)
        assert src["z"] == pytest.approx(src["h_J"]/x-p*p/4)
        # The affine Einstein coordinate was NOT substituted for original V.
        assert src["H_E"] == pytest.approx(saved["h_einstein"][index]-src["integral_p_squared"]/6)


def test_source_mass_and_areal_null_lapse_are_independent(geometry):
    src = model.source_traces((geometry["Va"]+geometry["Vb"])/2, geometry)
    bg = model.background(geometry["U0"], (geometry["Va"]+geometry["Vb"])/2)
    reconstructed = (model.M*src["H_E"]/src["x"]
                     +2*src["x"]*(bg["r_V"]*src["k_E"]+bg["r_U"]*src["H_E_V"]/bg["C"]))
    assert reconstructed == pytest.approx(src["mu_E"], rel=2e-12, abs=1e-12)


def test_zero_source_rectangle_has_exact_zero_metric_and_constraints(geometry):
    started = time.process_time()
    solution = model.solve_bulk_rectangle(ZeroField(geometry), geometry, source_scale=0., nodes=33)
    for t in (0., .5, 1.):
        sample = solution.sample(t, independent_order=32)
        for key in ("H_E", "c_E", "h_J", "c_J", "mu_E", "mu_J", "mass_defect_U",
                    "mass_defect_V", "null_U_integral_defect", "null_V_integral_defect"):
            assert np.count_nonzero(sample[key]) == 0
    assert time.process_time()-started < model.CPU_SECONDS


def test_actual_source_cannot_be_combined_with_zero_scalar(geometry):
    with pytest.raises(ValueError, match="incompatible"):
        model.solve_bulk_rectangle(ZeroField(geometry), geometry)


def test_domain_and_budget_fail_before_evolution(geometry):
    changed = ZeroField(geometry)
    changed.domain["Ua"] += .1
    with pytest.raises(ValueError, match="domains"):
        model.solve_bulk_rectangle(changed, geometry, source_scale=0.)
    with pytest.raises(ValueError, match="grid"):
        model.solve_bulk_rectangle(ZeroField(geometry), geometry, nodes=513)
    exhausted = model.BulkBudget(started=time.process_time()-121)
    with pytest.raises(RuntimeError, match="budget"):
        model.solve_bulk_rectangle(ZeroField(geometry), geometry, source_scale=0., budget=exhausted)


def test_source_data_fail_closed(geometry):
    with pytest.raises(ValueError, match="fixed characteristic"):
        model.source_traces(geometry["Va"]-1, geometry)
    with pytest.raises(ValueError, match="source scale"):
        model.source_traces(geometry["Va"], geometry, source_scale=2.)
    with pytest.raises(ValueError, match="order"):
        model.source_traces(geometry["Va"], geometry, order=16)


def test_source_compatibility_is_checked_during_evolution(geometry):
    class AlteredSource(ZeroField):
        def bulk(self, U, t, order=32):
            values = super().bulk(U, t, order)
            # Passes all three initial compatibility probes, then changes.
            if 0 < t < .5:
                values["p"][0] = 1.
            return values
    with pytest.raises(ValueError, match="incompatible"):
        model.solve_bulk_rectangle(AlteredSource(geometry), geometry, source_scale=0.)


@pytest.mark.parametrize("invalid", [-1, 0, 121, np.inf, np.nan, True])
def test_budget_cannot_relax_admitted_cpu_limit(invalid):
    with pytest.raises(ValueError):
        model.BulkBudget(cpu_seconds=invalid)


@pytest.mark.parametrize("invalid", [-1, 1.5, True])
def test_budget_counts_cannot_be_subtracted(invalid):
    with pytest.raises(ValueError, match="nonnegative integers"):
        model.BulkBudget(source_nodes=invalid)
    budget = model.BulkBudget()
    with pytest.raises(ValueError, match="nonnegative integers"):
        budget.check(metric=invalid)


def test_noninteger_grid_and_order_fail_cleanly(geometry):
    with pytest.raises(ValueError, match="grid"):
        model.solve_bulk_rectangle(ZeroField(geometry), geometry, nodes=33.)
    with pytest.raises(ValueError, match="orders"):
        model.solve_bulk_rectangle(ZeroField(geometry), geometry, source_order=64.)
