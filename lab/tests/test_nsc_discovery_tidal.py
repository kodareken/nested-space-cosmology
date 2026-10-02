"""Observer tides of the realized conformal metric. Saved states are read only."""
import json
from pathlib import Path

import pytest

import derive_nsc_discovery_tidal as cli
import derive_nsc_spherical_conformal_curvature as helper
from recursive_horizons import nsc_discovery_episode as episode
from recursive_horizons import nsc_discovery_response as response
from recursive_horizons import nsc_discovery_tidal as tidal


@pytest.fixture(scope="module")
def report():
    measured = tidal.measurement_report()
    assert measured["cpu_seconds"] < 30.0
    return measured


def _station(report, case_id, time_value):
    cases = {case["case_id"]: case for case in report["episode_v1"]["cases"]}
    for row in cases[case_id]["stations"]:
        if abs(row["coordinate_time"] - time_value) <= 1e-8:
            return row
    raise AssertionError((case_id, time_value))


def test_flat_charts_have_vanishing_four_curvature():
    exact = tidal.exact_suite()
    assert exact["minkowski"]["flat_curvature_max_abs"] == 0.0
    assert exact["minkowski"]["R_h"] == 2.0
    assert exact["minkowski"]["grad2"] == -1.0
    assert exact["milne"]["flat_curvature_max_abs"] < 1e-12
    assert exact["milne"]["R_h"] == pytest.approx(2.0)
    assert exact["milne"]["grad2"] == pytest.approx(-1.0)
    assert abs(exact["minkowski"]["algebra_gap"]) == 0.0
    assert abs(exact["milne"]["algebra_gap"]) < 1e-15


def test_bertotti_robinson_sign_and_weyl_identity():
    exact = tidal.exact_suite()["bertotti_robinson"]
    assert exact["R_h"] == pytest.approx(2.0)
    assert exact["R4"] == pytest.approx(0.0, abs=1e-12)
    assert exact["owned_W"] == pytest.approx(0.0, abs=1e-15)
    assert exact["R_0202"] == pytest.approx(0.0, abs=1e-15)
    assert exact["R_0101"] == pytest.approx(-1.0 / 9.0)
    assert exact["Ricci2"] == pytest.approx(4.0 / 81.0)
    assert exact["K"] == pytest.approx(8.0 / 81.0)
    assert exact["Ricci2"] > 0.0 and exact["K"] > 0.0
    assert abs(exact["algebra_gap"]) < 1e-12
    product = tidal.exact_suite()["product_time_sphere"]
    assert product["R_0101"] == 0.0 and product["R_0202"] == 0.0
    assert product["R4"] == pytest.approx(-0.5)
    assert product["K"] == pytest.approx(0.25)
    assert product["owned_W"] == pytest.approx(product["K"] - 2.0 * product["Ricci2"] + product["R4"] ** 2 / 3.0)


def test_saved_hashes_and_no_evolution(report):
    assert report["saved_state_hashes_unchanged"] is True
    assert report["input_hashes"]["unchanged"] is True
    assert report["chi_substituted_for_curvature"] is False
    assert report["new_ode"] is False
    assert report["geometry_evolved"] is False
    assert report["production_record_written"] is False
    assert report["singularity_declared"] is False
    assert report["theory_declared_dead"] is False
    assert report["finite_window"]["chart_positive"] is True
    assert report["finite_window"]["values_finite"] is True
    assert report["confirmation_nf512"]["status"] == "measured"
    for campaign, prefix in tidal.FROZEN_PRODUCING_COMMITS.items():
        binding = report["frozen_producer_bindings"][campaign]
        assert binding["producing_commit"].startswith(prefix)
        assert binding["physics_matches_saved_producer"] is True
        assert tidal.sha256_file(binding["envelope"]) == binding["envelope_sha256"]
    sample = next(iter(report["input_hashes"]["before"]))
    assert tidal.sha256_file(sample) == report["input_hashes"]["before"][sample]


def test_frozen_jets_are_zero_while_the_field_flows(report):
    frozen = report["comparisons"]["frozen"]
    assert frozen["geometry_jets_zero"] is True
    assert frozen["tide_constant"] is True
    assert frozen["field_still_flows"] is True
    coupled = _station(report, "nf256_coupled_dt0.0005", 0.3)
    frozen_row = _station(report, "nf256_frozen_geometry_dt0.0005", 0.3)
    assert frozen_row["Q_t_max"] == 0.0 and frozen_row["r_tt_max"] == 0.0
    assert frozen_row["field_rate_max"] > 1.0
    assert frozen_row["summaries"]["owned_W"]["max_abs"] < 0.01
    assert coupled["summaries"]["owned_W"]["max_abs"] == pytest.approx(0.112, rel=0.02)
    assert coupled["summaries"]["owned_W"]["max_abs"] > 10.0 * frozen_row["summaries"]["owned_W"]["max_abs"]


def test_q_acceleration_matches_the_owned_helper():
    gap = tidal.static_conformal_section_gap()
    assert max(gap.values()) < 1e-12
    assert tidal.directional_radius_chain_gap()["gap"] < 1e-10
    directory = tidal.EPISODE_DIR
    record, arrays = episode.load_checkpoint(directory, "nf128_coupled_dt0.0005", 2)
    pair = episode.pair_from_arrays(arrays, record)
    state = episode.state_from_arrays(arrays, record["momentum_representation"])
    jets = tidal.analytic_accelerations(pair, state, "coupled")
    _q_dot, q_ddot, before = helper.actual_q_second_rate(
        pair.grid, jets["nodal"], jets["rate"], jets["bundle"],
    )
    assert float(abs(q_ddot - jets["Q_tt"]).max()) == 0.0
    assert float(abs(before - jets["Q_tt_before_projection"]).max()) == 0.0
    # Independent full analytic Jv along the full actual state rate.
    full_rate = tidal.nested.rates(pair, state)
    tangent = response.StateTangent(
        *(getattr(full_rate, name) for name in tidal.nested.STATE_NAMES),
        tidal.np.zeros(pair.weights.size),
    )
    jv = response.nested_rate_jacobian_vector(pair, state, tangent)
    for name in ("Q", "r"):
        lifted = tidal.galerkin.prolong_geometry(pair.grid, pair.geometry_map @ getattr(jv, name))
        assert lifted == pytest.approx(jets[name + "_tt"], rel=1e-10, abs=1e-10)
    gap = abs(jets["unprojected_r_formula"] - jets["unprojected_r_owner"]).max()
    assert float(gap) < 1e-12
    # Mixed Hessian contractions use the Lorentzian sign, including -8 B01^2.
    package = tidal.curvature_from_local_jets(
        0.8, 2.3, 0.04, -0.03, -0.06, 0.02, 0.05,
        0.1, -0.2, 0.07, 0.08, -0.09,
    )
    base, tides = package["base"], package["tides"]
    A, b, c, d, S = (base[name] for name in ("A", "B00", "B11", "B01", "S"))
    assert float(base["Ricci2"]) == pytest.approx((A-2*b)**2+(-A-2*c)**2-8*d**2+2*(S+b-c)**2)
    assert float(base["K"]) == pytest.approx(4*A**2+8*(b*b+c*c-2*d*d)+4*S*S)
    assert float(base["factorized_W"]) == pytest.approx(tides["owned_W"])
    assert float(tides["R_0101"]) == pytest.approx(-A)
    assert float(tides["R_0202"]) == pytest.approx(b)


def test_late_tide_is_stable_while_weyl_cancels(report):
    late = _station(report, "nf256_coupled_dt0.0005", 3.0)
    early = _station(report, "nf256_coupled_dt0.0005", 0.3)
    middle = _station(report, "nf256_coupled_dt0.0005", 1.0)
    assert 2.0e-5 < late["Q_min"] < late["Q_max"] < 5.0e-5
    assert 20.0 < late["r_min"] < late["r_max"] < 22.0
    assert late["conditioning"]["temporal_raw"]["numerator_max_abs"] > 1.0e10
    assert late["conditioning"]["full_Rh"]["condition"] > 1.0e6
    assert 1.0 < late["conditioning"]["angular_tide"]["condition"] < 2.0
    assert late["summaries"]["R_h"]["max"] < 1.0e5
    assert middle["conditioning"]["algebra_gap_max"] < 1e-8
    assert early["conditioning"]["algebra_gap_max"] < 1e-9
    assert late["conditioning"]["algebra_gap_max"] < 1.0
    assert late["conditioning"]["algebra_gap_max"] / late["summaries"]["owned_W"]["max_abs"] < 1e-3
    assert late["worldline_x2"]["R_0101"] == pytest.approx(-524203.28, rel=1e-6)
    assert late["worldline_x2"]["R_0202"] == pytest.approx(5632447.06, rel=1e-6)
    assert late["worldline_x2"]["Ricci2"] == pytest.approx(2.545895e14, rel=1e-6)
    assert late["worldline_x2"]["K"] == pytest.approx(5.091790e14, rel=1e-6)
    space = report["comparisons"]["space_nf256_versus_nf128_dt0.0005"]["3.0"]
    tide_gap = abs(space["R_0202"]["max_abs"]["relative_difference"])
    weyl_gap = abs(space["owned_W"]["max_abs"]["relative_difference"])
    scalar_gap = abs(space["R4"]["max_abs"]["relative_difference"])
    assert tide_gap < 0.01
    assert abs(space["K"]["max_abs"]["relative_difference"]) < 0.01
    assert abs(space["R_0101"]["max_abs"]["relative_difference"]) < 0.01
    assert weyl_gap > 0.2
    assert scalar_gap > 0.5
    assert tide_gap < weyl_gap
    for channel in ("R_0101", "R_0202", "K", "Ricci2"):
        assert space[channel]["shape"]["l2_relative"] < 0.01
    step = report["comparisons"]["timestep_nf256_dt0.001_versus_dt0.0005"]["3.0"]
    assert abs(step["R_0202"]["relative_difference"]) < 1e-8
    assert abs(step["K"]["relative_difference"]) < 1e-8
    world = report["comparisons"]["growth"]["worldline_x2"]
    assert world["x"]["3.0"] == pytest.approx(2.0)
    ratio = world["channels"]["R_0202"]["T0.3_to_3"]["ratio"]
    assert ratio > 1.0e6
    assert world["channels"]["K"]["values"]["3.0"] > 1.0e12
    clock = report["comparisons"]["growth"]["clock_comparison"]["worldline_x2"]["T1_to_3"]
    assert 0.0 < clock["proper_dtau_x2"] < clock["coordinate_dt"]
    assert abs(clock["per_proper_time_x2"]) > abs(clock["per_coordinate_time"])
    peak = report["comparisons"]["growth"]["peak_R0202"]["channels"]["R_0202"]
    assert abs(peak["values"]["3.0"]) / abs(peak["values"]["0.3"]) > 1.0e6
    assert peak["values"]["0.3"] < 0.0 < peak["values"]["3.0"]
    assert all(abs(value) < 1.0e20 for value in peak["values"].values())


def test_nf512_tide_moves_less_than_the_weyl_scalar(report):
    finer = report["comparisons"]["space_nf512_versus_nf256_dt0.0005"]
    assert finer is not None
    tide = abs(finer["3.0"]["R_0202"]["max_abs"]["relative_difference"])
    weyl = abs(finer["3.0"]["owned_W"]["max_abs"]["relative_difference"])
    assert tide < 1e-4
    assert weyl < 2e-3
    assert tide < weyl
    world = abs(finer["3.0"]["R_0202"]["worldline_x2"]["relative_difference"])
    assert world < 1e-4


def test_cli_records_once_and_checks_hashes(report, tmp_path, monkeypatch):
    monkeypatch.setattr(cli, "measurement_report", lambda **_kwargs: report)
    destination = tmp_path / "tidal.json"
    assert cli.main(["--write", str(destination)]) == 0
    assert destination.is_file()
    saved = json.loads(destination.read_text(encoding="utf-8"))
    assert saved["singularity_declared"] is False
    assert saved["saved_state_hashes_unchanged"] is True
    with pytest.raises(SystemExit):
        cli.main(["--write", str(destination)])
    assert cli.main(["--check", str(destination)]) == 0
    before = tidal.sha256_file(episode.LAB / "results/development/nsc-discovery-episode-v1/manifest.json")
    assert before == saved["input_hashes"]["before"][str((episode.LAB / "results/development/nsc-discovery-episode-v1/manifest.json").resolve())]
    saved["producer_hashes"][next(iter(saved["producer_hashes"]))] = "0" * 64
    destination.write_text(json.dumps(saved), encoding="utf-8")
    with pytest.raises(SystemExit, match="hash changed"):
        cli.main(["--check", str(destination)])
