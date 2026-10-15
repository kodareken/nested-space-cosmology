"""Manufactured forcing cache checks; no physical field preparation."""
import copy
import time
from types import SimpleNamespace

import numpy as np
import pytest

from recursive_horizons import nsc_atlas_forcing_cache as model
from recursive_horizons.nsc_atlas_bulk_response import BulkBudget
from recursive_horizons.nsc_atlas_field_preparation import radius_from_null
from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget


@pytest.fixture(autouse=True)
def bounded_manufactured_check():
    started = time.process_time()
    yield
    assert time.process_time()-started < 30


@pytest.fixture
def geometry():
    U0, Ua, Ub = -.32130105902052203, -.2992702973503126, -.2747468030926288
    Va, Vb = -2.565695357068247, -2.262011613052105
    support = [Va+(Vb-Va)/4, Va+3*(Vb-Va)/4]
    return dict(order=64, U0=U0, Ua=Ua, Ub=Ub, Va=Va, Vb=Vb,
                x_corner=.4, x_last=float(radius_from_null(U0, Vb)),
                support_V=support, support_x=radius_from_null(U0, support).tolist())


class ManufacturedWedge:
    def __init__(self, geometry):
        self.domain = dict(geometry)
        self.budget = PreparedFieldBudget()
        self.geometry = {"wall": self.wall}

    def wall(self, V):
        d = self.domain
        t = (np.asarray(V)-d["Va"])/(d["Vb"]-d["Va"])
        return np.where(V == d["Va"], d["Ua"], np.where(V == d["Vb"], d["Ub"],
                        d["Ua"]+(d["Ub"]-d["Ua"])*(t+t*(1-t)/5)))


@pytest.fixture
def setup(geometry, monkeypatch):
    wedge = ManufacturedWedge(geometry)
    budget = BulkBudget(started=wedge.budget.started, external=wedge.budget)
    calls = []
    parameters = dict(time_power=3, space_power=1, bad_wall=False, bad_quiet=False)

    class ManufacturedForcing:
        def __init__(self, supplied_wedge, supplied_geometry, U, *, budget, **settings):
            assert supplied_wedge is wedge and supplied_geometry == geometry
            self.U = np.asarray(U)
            assert self.U[0] == geometry["U0"] and self.U[-1] == geometry["Ub"]
            assert np.all(np.diff(self.U) > 0)
            self.budget = budget

        def slice(self, t):
            dV = geometry["Vb"]-geometry["Va"]
            V = geometry["Va"]+t*dV
            wall = float(wedge.geometry["wall"](V))
            U = np.minimum(self.U, wall)
            U[-1] = wall
            eta = (U-geometry["U0"])/(wall-geometry["U0"])
            q = max(0., t-.25)
            h = q**parameters["time_power"]
            hV = parameters["time_power"]*q**(parameters["time_power"]-1)/dV
            F, G = h*(1-eta), h*(2+eta**parameters["space_power"])
            if parameters["bad_wall"]:
                F[-1] = 1e-5
            if parameters["bad_quiet"] and t <= .25:
                G[0] = 1e-5
            calls.append((t, eta.copy()))
            wedge.budget.check(U.size)
            self.budget.check(metric=U.size, source=64)
            return dict(t=t, V=V, U_wall=wall, U=U, F=F, G=G, g1=2*h, g1_V=2*hV,
                        source_p_defect=0., source_p_V_defect=0., radial_nodes=2*(U.size-1),
                        field_points=U.size, field_batches=1)

    monkeypatch.setattr(model, "ScalarForcing", ManufacturedForcing)
    return SimpleNamespace(wedge=wedge, budget=budget, geometry=geometry,
                           calls=calls, parameters=parameters)


def cache(setup, **kwargs):
    return model.ForcingTimeCache(setup.wedge, setup.geometry,
                                 kwargs.pop("eta_grid", np.linspace(0, 1, 9)),
                                 kwargs.pop("training_times", np.linspace(0, 1, 17)),
                                 budget=setup.budget, **kwargs)


def targets(setup, t, eta):
    d = setup.geometry
    wall = setup.wedge.geometry["wall"](d["Va"]+t*(d["Vb"]-d["Va"]))
    return d["U0"]+np.asarray(eta)*(wall-d["U0"])


def test_cubic_time_linear_eta_and_hermite_source_are_exact(setup):
    result = cache(setup, validation_times=(.03125, .40625, .84375))
    for t in (.1, .25, .37, .75, .89, 1.):
        eta = np.array([0., .113, .63, 1.])
        row = result.eval(targets(setup, t, eta), t)
        h = max(0., t-.25)**3
        dV = setup.geometry["Vb"]-setup.geometry["Va"]
        assert row["F"] == pytest.approx(h*(1-eta), abs=2e-15)
        assert row["G"] == pytest.approx(h*(2+eta), abs=2e-15)
        assert row["g1"] == pytest.approx(2*h, abs=2e-15)
        assert row["g1_V"] == pytest.approx(6*max(0., t-.25)**2/dV, abs=2e-14)
        assert row["F"][0] == row["Fsource"]
        assert row["F"][-1] == 0
    assert result.eval(targets(setup, .9, [0, 1]), .9)["g1"] > 0
    for row in result.validation_rows:
        assert max(abs(row[key]) for key in model.DEFECT_KEYS) < 2e-13


def test_fixed_U_composition_needs_no_advection_and_clips_inactive_nodes(setup):
    result = cache(setup)
    U = np.linspace(setup.geometry["U0"], setup.geometry["Ub"], 31)
    for t in (.37, .61, .91):
        row = result.eval(U, t)
        wall = targets(setup, t, 1.)
        eta = (np.minimum(U, wall)-setup.geometry["U0"])/(wall-setup.geometry["U0"])
        h = (t-.25)**3
        assert row["F"] == pytest.approx(h*(1-eta), abs=2e-15)
        assert row["G"] == pytest.approx(h*(2+eta), abs=2e-15)
        assert np.all(row["U"][U >= wall] == wall)
        assert np.all(row["F"][U >= wall] == 0)


def test_eval_cannot_call_the_field_metric_or_source(setup, monkeypatch):
    result = cache(setup)
    before = (len(setup.calls), setup.wedge.budget.kernel_values,
              setup.budget.metric_values, setup.budget.source_nodes)
    def forbidden(*args, **kwargs):
        raise AssertionError("physical slice called during cached evaluation")
    monkeypatch.setattr(model, "ScalarForcing", forbidden)
    for t in np.linspace(0, 1, 101):
        result.eval(np.linspace(setup.geometry["U0"], setup.geometry["Ub"], 129), t)
    assert before == (len(setup.calls), setup.wedge.budget.kernel_values,
                      setup.budget.metric_values, setup.budget.source_nodes)


def test_one_ulp_terminal_wall_keeps_exact_zero_flux(setup):
    original = setup.wedge.wall
    terminal = np.nextafter(setup.geometry["Ub"], np.inf)
    setup.wedge.geometry["wall"] = lambda V: np.where(
        V == setup.geometry["Vb"], terminal, original(V))
    result = cache(setup)
    row = result.eval([setup.geometry["U0"], setup.geometry["Ub"]], 1.)
    assert row["U"][-1] == terminal
    assert row["F"][-1] == 0


def test_withheld_midpoints_separate_temporal_from_combined_spatial_error(setup):
    setup.parameters["space_power"] = 2
    result = cache(setup, validation_times=(.40625, .84375))
    for row in result.validation_rows:
        h = (row["t"]-.25)**3
        assert row["temporal_G_defect"] < 2e-15
        assert row["combined_G_defect"] == pytest.approx(h*(1/8)**2/4, abs=2e-15)
        assert row["combined_G_defect"] > 1e-6
    summary = result.summary()
    assert summary["uniform_error_enclosure"] is False
    assert summary["withheld_used_for_training"] is False
    assert "combined" in summary["midpoint_defect_scope"]


def test_nested_training_and_shared_withheld_slices_have_65_unique_calls(setup):
    eta = np.linspace(0, 1, 129)
    withheld = (np.arange(32)+.5)/32
    provider = model.ForcingSliceProvider(setup.wedge, setup.geometry, budget=setup.budget)
    coarse = cache(setup, eta_grid=eta, validation_times=withheld, slice_provider=provider)
    fine = cache(setup, eta_grid=eta, training_times=np.linspace(0, 1, 33),
                 validation_times=withheld, slice_provider=provider)
    assert len(setup.calls) == len(provider.cache) == 65
    assert sorted(len(row["eta"]) for row in provider.cache.values()) == [129]*33+[257]*32
    assert {row["t"] for row in fine.training_rows}.isdisjoint(withheld)
    assert all(row["eta"].size == 129 for row in coarse.training_rows+fine.training_rows)
    assert provider.cache_hits == 49
    assert provider.summary()["kernel_values"] == setup.wedge.budget.kernel_values
    assert provider.summary()["metric_values"] == setup.budget.metric_values
    assert provider.summary()["source_nodes"] == setup.budget.source_nodes


def test_validation_never_refits_and_scales_the_source_slope_by_its_cell_width(setup):
    setup.parameters["time_power"] = 4
    result = cache(setup)
    query = targets(setup, .43, [0., .4, 1.])
    before = result.eval(query, .43)
    result.validate((.40625, .84375))
    after = result.eval(query, .43)
    for key in before:
        assert np.array_equal(before[key], after[key])
    for row in result.validation_rows:
        assert row["cell_width_V"] == pytest.approx(result.dV/16)
        assert row["scaled_g1_V_defect"] == row["cell_width_V"]*row["g1_V_defect"]
    # A non-midpoint exposes the derivative error of a quartic Hermite fit.
    result.validate((.42,))
    assert abs(result.validation_rows[-1]["g1_V_defect"]) > 1e-5
    assert result.summary()["blocks"][1]["maximum_absolute_defects"]["scaled_g1_V_defect"] > 0
    count = len(setup.calls)
    result.validate((.42,))
    assert len(setup.calls) == count


def test_saved_slices_replay_without_physical_fallback(setup, monkeypatch):
    result = cache(setup, validation_times=(.40625, .84375))
    records = list(result.slice_provider.cache.values())
    monkeypatch.setattr(model, "ScalarForcing", lambda *a, **k: pytest.fail("physical replay"))
    provider = model.ForcingSliceProvider.from_records(setup.wedge, setup.geometry, records,
                                                       budget=setup.budget)
    replay = cache(setup, validation_times=(.40625, .84375), slice_provider=provider)
    for old, new in zip(result.validation_rows, replay.validation_rows):
        for key in model.DEFECT_KEYS:
            assert old[key] == new[key]
    assert provider.summary()["physical_slices"] == len(records)
    for key in provider.counts:
        assert provider.counts[key] == result.slice_provider.counts[key]
    with pytest.raises(ValueError, match="not retained"):
        provider.slice(.4, np.linspace(0, 1, 9))
    assert not provider.cache[next(iter(provider.cache))]["F"].flags.writeable


@pytest.mark.parametrize("mutation", ["U", "eta", "nan", "cost", "duplicate", "quiet", "wall"])
def test_saved_slice_tampering_is_rejected(setup, mutation):
    result = cache(setup)
    records = copy.deepcopy(list(result.slice_provider.cache.values()))
    if mutation == "U":
        records[0]["U"][1] += .001
    elif mutation == "eta":
        records[0]["eta"][0] = .01
    elif mutation == "nan":
        records[0]["G"][0] = np.nan
    elif mutation == "cost":
        records[0]["slice_costs"]["kernel_values"] = -1
    elif mutation == "duplicate":
        records.append(records[0])
    elif mutation == "quiet":
        records[0]["G"][0] = 1e-4
    else:
        records[-1]["F"][-1] = 1e-4
    with pytest.raises(ValueError):
        model.ForcingSliceProvider.from_records(setup.wedge, setup.geometry, records,
                                                budget=setup.budget)


@pytest.mark.parametrize("kwargs", [
    {"eta_grid": [0., 0., 1.]}, {"eta_grid": [0., 1e-16, 1.]},
    {"eta_grid": [0., np.nan, 1.]}, {"eta_grid": np.linspace(0, 1, 258)},
    {"eta_grid": [False, True]}, {"eta_grid": [.1, .5, 1.]},
    {"training_times": np.linspace(0, 1, 16)}, {"training_times": [0., .25, .75, 1.]},
    {"validation_times": [.25]}, {"validation_times": [.4, .3]},
    {"validation_times": [np.nan]}, {"radial_order": True}, {"source_order": 16},
])
def test_invalid_grids_source_blocks_orders_and_withheld_times_fail_before_solving(setup, kwargs):
    with pytest.raises(ValueError):
        cache(setup, **kwargs)
    assert not setup.calls


@pytest.mark.parametrize("which", ["bad_wall", "bad_quiet"])
def test_exact_boundary_and_quiet_input_contract_are_checked(setup, which):
    setup.parameters[which] = True
    with pytest.raises(ValueError, match="wall flux|quiet interval"):
        cache(setup)


def test_same_preparation_budget_and_provider_settings_are_required(setup):
    bad_budget = BulkBudget()
    with pytest.raises(ValueError, match="shared budget"):
        model.ForcingSliceProvider(setup.wedge, setup.geometry, budget=bad_budget)
    provider = model.ForcingSliceProvider(setup.wedge, setup.geometry, budget=setup.budget)
    with pytest.raises(ValueError, match="same forcing slice"):
        cache(setup, slice_provider=provider, radial_order=4)
    changed = dict(setup.geometry, Ub=setup.geometry["Ub"]+.0001)
    with pytest.raises(ValueError):
        model.ForcingSliceProvider(setup.wedge, changed, budget=setup.budget)
    assert not setup.calls


def test_eval_domain_shapes_and_cached_budget_remain_guarded(setup):
    result = cache(setup)
    for U, t in (([setup.geometry["U0"]-.001], .5), ([setup.geometry["Ub"]+.001], .5),
                 ([[setup.geometry["U0"]]], .5), ([setup.geometry["U0"]], True),
                 ([setup.geometry["U0"]], 1.001)):
        with pytest.raises(ValueError):
            result.eval(U, t)
    scalar = result.eval(setup.geometry["U0"], .8)
    assert np.asarray(scalar["F"]).shape == ()
    setup.budget.started -= 121
    with pytest.raises(RuntimeError, match="compute budget"):
        result.eval(setup.geometry["U0"], .8)


def test_shared_slice_limit_is_checked_before_any_additional_queries(setup):
    result = cache(setup)
    initial_calls = len(setup.calls)
    with pytest.raises(ValueError, match="finite slice count"):
        result.validate((np.arange(129)+.37)/129)
    assert len(setup.calls) == initial_calls
