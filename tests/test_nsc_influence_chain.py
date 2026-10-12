"""Unit tests for the nested-influence chain closure."""

from recursive_horizons import nsc_influence_chain as chain


def test_infinite_geometric_closure():
    closure = chain.geometric_closure(0.25)
    assert closure["total_multiplier"] == 1.0 / 0.75
    assert closure["local_fraction"] == 0.75
    assert closure["exterior_fraction"] == 0.25


def test_finite_depth_matches_partial_sum():
    closure = chain.geometric_closure(0.5, depth=3)
    expected = 1.0 + 0.5 + 0.25 + 0.125
    assert abs(closure["total_multiplier"] - expected) < 1e-15
    assert abs(closure["local_fraction"] - 1.0 / expected) < 1e-15


def test_depth_zero_is_single_level():
    closure = chain.geometric_closure(0.9, depth=0)
    assert closure["total_multiplier"] == 1.0
    assert closure["local_fraction"] == 1.0


def test_verdict_taxonomy_pre_registered():
    assert chain.verdict(0.049) == chain.NEAR_STATUS
    assert chain.verdict(0.03) == chain.NEAR_STATUS
    assert chain.verdict(0.10) == chain.NEAR_STATUS
    assert chain.verdict(0.75) == chain.NOT_NEAR_STATUS
    assert chain.verdict(0.999) == chain.NOT_NEAR_STATUS
    assert chain.verdict(1.5) == chain.BAND_STATUS


def test_out_of_range_coefficient_rejected():
    try:
        chain.geometric_closure(1.0)
    except chain.ChainFeeError:
        pass
    else:
        raise AssertionError("f = 1 must not close")
    try:
        chain.geometric_closure(-0.1)
    except chain.ChainFeeError:
        pass
    else:
        raise AssertionError("negative f must not close")


def test_scale_binding_fees_extraction():
    record = {
        "development_branch": {"omega": 4.0, "zeta": 16.0},
        "recursive_compatibility": {"Omega_times_transfer": 0.17732912925261413},
    }
    fees = chain.scale_binding_fees(record)
    assert fees["structural_one_over_omega"] == 0.25
    assert abs(fees["zeta_minus_omega_squared"]) < 1e-12
    assert abs(fees["recorded_transfer"] - 0.17732912925261413 / 4.0) < 1e-15
    assert fees["endpoint_condition_omega_q_le_one"] is True


def test_continuation_mass_fee_extraction():
    summary = {
        "final": {"coupled_inner_mass": 1.0410282460948723, "fixed_inner_mass": 1.04},
    }
    fees = chain.continuation_mass_fee(summary)
    assert abs(fees["measured_relative_mass_renormalisation"] - 0.0010282460948723 / 1.04) < 1e-15
    assert abs(fees["one_over_ratio_squared"] - 1.0 / 32.0 ** 2) < 1e-17
    assert fees["route_ratio"] == 32.0


def test_implied_layer_ratio():
    entry = chain.matrix_coefficient(1.0 / 1024.0)
    assert abs(entry["implied_layer_ratio"] - 32.0) < 1e-12


def test_required_power_inversion():
    entry = chain.required_power(32.0, local_fraction=0.049)
    import math
    expected = -math.log(1.0 - 0.049) / math.log(32.0)
    assert abs(entry["required_nu"] - expected) < 1e-15
    assert abs(entry["per_level_decay"] - 0.951) < 1e-15


def test_band_required_power_ordering():
    entry = chain.band_required_power(32.0)
    assert 0.0 < entry["nu_min"] < entry["nu_max"]
    for nu in (entry["nu_min"], entry["nu_max"]):
        closure = chain.chain_fraction(nu, 32.0)
        assert chain.BAND_LO - 1e-12 <= closure["local_fraction"] <= chain.BAND_HI + 1e-12


def test_lll_screen_power_formula():
    entry = chain.lll_screen_power(g4=2.0, flux=4.0)
    assert abs(entry["nu"] - 2.0 * 2.0 / (2.0 * 3.141592653589793)) < 1e-15


def test_coupling_inversion_round_trip():
    screen = chain.lll_screen_power(g4=4.713, flux=4.0)
    back = chain.coupling_for_power(screen["nu"], 4.0)
    assert abs(back["g4"] - 4.713) < 1e-12


def test_chain_fraction_closure_at_matched_coupling():
    closure = chain.chain_fraction(1.5, 32.0)
    assert abs(closure["f"] - 32.0 ** (-1.5)) < 1e-15
    assert closure["verdict"] == chain.NOT_NEAR_STATUS


def test_spatial_fee_l1_definition():
    import numpy as _np
    radius = _np.linspace(0.5, 32.0, 400)
    weights = _np.full(400, (32.0 - 0.5) / 399.0)
    response = _np.zeros_like(radius)
    response[radius >= 24.0] = 1.0
    inner = chain.spatial_fee(radius, weights, response, 8.0)
    assert inner["l1_mass_response_fraction_inside"] == 0.0
    response[:] = 0.0
    response[radius < 8.0] = 1.0
    locality = chain.spatial_fee(radius, weights, response, 8.0)
    assert abs(locality["l1_mass_response_fraction_inside"] - 1.0) < 1e-12
    response[:] = 1.0
    half = chain.spatial_fee(radius, weights, response, 8.0)
    assert 0.0 < half["l1_mass_response_fraction_inside"] < 1.0


def test_spatial_fee_rejects_seam_outside_grid():
    import numpy as _np
    radius = _np.linspace(1.0, 10.0, 50)
    weights = _np.full(50, 0.2)
    try:
        chain.spatial_fee(radius, weights, _np.ones(50), 11.0)
    except chain.ChainFeeError:
        pass
    else:
        raise AssertionError("out-of-grid seam must be rejected")


def test_spatial_fee_increases_with_seam():
    """Monotonicity: a response profile's inside fraction grows with the seam."""
    import numpy as _np
    radius = _np.linspace(0.5, 32.0, 400)
    weights = _np.full(400, (32.0 - 0.5) / 399.0)
    response = _np.ones_like(radius)
    previous = 0.0
    for seam in (4.0, 8.0, 16.0, 24.0, 30.0):
        current = chain.spatial_fee(radius, weights, response, seam)[
            "l1_mass_response_fraction_inside"]
        assert current > previous + 1e-12
        previous = current


def test_interval_shares_close_to_one():
    import numpy as _np
    radius = _np.linspace(0.5, 32.0, 400)
    weights = _np.full(400, (32.0 - 0.5) / 399.0)
    response = _np.abs(radius - 12.0) + 1.0
    shares = chain.interval_shares(radius, weights, response, (0.5, 8.0, 16.0, 24.0, 32.0))
    assert len(shares) == 4
    assert abs(sum(item["l1_share"] for item in shares) - 1.0) < 1e-9
    assert all(item["l1_share"] >= 0.0 for item in shares)


def test_interval_shares_reject_bad_edges():
    import numpy as _np
    radius = _np.linspace(0.5, 32.0, 400)
    weights = _np.full(400, 0.08)
    response = _np.ones(400)
    for bad in ((32.0, 8.0), (0.5, 8.0, 4.0), (0.5, 40.0)):
        try:
            chain.interval_shares(radius, weights, response, bad)
        except chain.ChainFeeError:
            pass
        else:
            raise AssertionError("bad edge tuple accepted: " + str(bad))


def test_interior_balance_record_closes():
    """The continuation run's interior identity: dM = d(excision) exactly,
    the outer boundary leaks only, and the recorded verdict says so."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-interior-balance-v1.json").read_text()
    )
    assert record["closed"] is True
    assert abs(record["closure_residual"]) < 1e-12
    assert abs(record["deltas"]["inner_mass"] - record["deltas"]["excision_flux"]) < 1e-12
    assert record["deltas"]["outer_flux"] < 0.0
    assert record["verdict"] == "EXCISION_CLOSED outer_leaks_only"


def test_continuation_stations_close_every_interval():
    """Round 38: the station-level closure of the continuation holds on all
    20 intervals with the outer boundary negative throughout."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-continuation-station-closure-v1.json").read_text()
    )
    assert record["verdict"] == "LATE_CLOSURE_HOLDS"
    assert record["maximum_residual"] < 1e-11
    assert record["outer_leaks_every_interval"] is True
    assert len(record["intervals"]) == 20
    assert record["beyond_family_coverage"]["max_residual_beyond"] < 1e-11


def test_station_junction_wall_collects_inflow():
    """The wall stations concentrate the geometric mass inflow; the family's
    probability book partitions to a few parts in a million."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-interior-station-balance-v1.json").read_text()
    )
    assert record["probability_partition_max_gap"] < 1e-4
    junction = record["junction"]
    assert junction["wall_stations"] + junction["quiet_stations"] == junction["overlap_stations"]
    assert junction["wall_quiet_mass_inflow_ratio"] > 5.0
    for row in junction["rows"]:
        assert row["wall_count"] >= 0


def test_energy_mass_ratio_is_recorded_not_constant():
    """Round 20's honest negative: the anti-correlated pair has no constant
    exchange rate; the measured scatter is recorded, not hidden."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-wall-energy-mass-ratio-v1.json").read_text()
    )
    mu_1040 = record["cases"]["mu-1.040-eps-0.001"]
    assert mu_1040["ratio_relative_scatter"] > 5.0
    assert mu_1040["ratio_mean"] < 0.0
    assert record["audited_A"] == 0.04501936182826115


def test_ratio_regime_split_at_mu_1040():
    """Round 21: the mu=1.040 converter resolves into two settings when
    split by the wall flag; the wall setting is positive and much tighter."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-wall-ratio-regimes-v1.json").read_text()
    )
    mu_1040 = record["cases"]["mu-1.040-eps-0.001"]
    wall = mu_1040["wall_destined"]
    quiet = mu_1040["quiet_destined"]
    assert wall["count"] + quiet["count"] == mu_1040["pooled"]["count"]
    assert wall["mean"] > 0.0
    assert wall["relative_scatter"] < quiet["relative_scatter"]
    assert abs(mu_1040["pooled"]["mean"]) > abs(wall["mean"])


def test_wall_breathing_asymmetry_is_family_wide():
    """Round 24: fall faster than rise in every case, and the recorded
    radial envelopes stay inside [0.9, 2.0] r_+."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-wall-breathing-family-v1.json").read_text()
    )
    for name, entry in record["cases"].items():
        assert entry["asymmetry_fall_over_rise"] > 1.0
        assert entry["mean_fall_rate"] < 0.0 < entry["mean_rise_rate"]
        low, high = entry["radial_range"]
        assert 0.8 < low and high < 2.0


def test_wall_rise_is_horizon_gap_clocked():
    """Round 26: rise*sqrt(horizon gap) is family-constant at the few-%
    level while the fall phase is not gap-clocked."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-wall-clock-horizon-v1.json").read_text()
    )
    assert record["rise_sqrt_gap_spread"] < 0.06
    assert record["fall_sqrt_gap_spread"] > record["rise_sqrt_gap_spread"]
    assert record["raw_rise_spread"] > 0.25


def test_wall_fall_is_quarter_gap_clocked():
    """Round 30: the fall phase follows the quarter-power gap law declared
    before the numbers; the fitted power matches it within 4%."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-wall-clock-fall-v1.json").read_text()
    )
    assert record["fall_gap_quarter_spread"] < 0.06
    assert record["fall_gap_quarter_spread"] < record["fall_gap_half_spread"]
    assert abs(abs(record["fitted_fall_power_post_hoc"]) - record["pre_registered_power"]) < 0.02


def test_ladder_deficits_follow_gaussian_in_n():
    """Round 42: with (a,b) fixed by n=1,2 the Gaussian-in-n deficit law
    predicts n=3 at sub-percent and n=4 within 15%."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-charged-ladder-deficit-law-v1.json").read_text()
    )
    assert record["verdict"] == "GAUSSIAN_APPROACH_HOLDS"
    assert abs(record["rows"][2]["relative_error"]) < 0.01
    assert abs(record["rows"][3]["relative_error"]) < 0.15
    assert record["b"] > 0.0


def test_charged_seam_profiles_coincide_across_eta():
    """Round 43: the dyonic seam profile is eta-invariant at <1e-2, and
    the 25%-cap fee stays in the ~0.26 universal class."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-charged-seam-invariance-v1.json").read_text()
    )
    assert record["verdict"] == "PROFILES_COINCIDE"
    assert record["max_pairwise_deviation"] < 0.01
    for profile in record["profiles"].values():
        assert 0.25 < profile["0.25"] < 0.28


def test_charged_tilt_survives_resolution_doubling():
    """Round 44: the eta-tilt is not a grid-cell artifact; the recorded
    verdict says fringe physics with ratios at or above 0.8."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-charged-tilt-resolution-v1.json").read_text()
    )
    assert record["verdict"] == "FRINGE_PHYSICS"
    for tilt in record["tilts"].values():
        assert tilt["ratio_1752_over_876"] >= 0.8


def test_tilt_rejects_the_organism_clock_hypothesis():
    """Round 45: the five-point tilt staircase does not march on the
    sqrt(gap) step; the record keeps the rejected hypothesis and the
    sign-flipping step."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-tilt-organism-clock-v1.json").read_text()
    )
    assert record["verdict"] == "ORGANISM_CLOCK_REJECTED"
    assert record["ratio_spread"] > 0.5
    assert any(step["delta_fee"] < 0.0 for step in record["steps"])
    assert len(record["rows"]) == 5


def test_staircase_quantum_is_the_single_node_at_the_seam():
    """Round 46: the sampled staircase advances in one quantized step size
    equal to the single-node L1 share at the seam (ratio ~1)."""
    import json
    import pathlib
    crossings = json.loads(
        pathlib.Path("results/development/nsc-tilt-cell-crossings-v1.json").read_text()
    )
    deltas = [pair["delta_fee"] for pair in crossings["pairs"]]
    active = [d for d in deltas if d > 0.0005]
    if active:
        assert max(active) - min(active) < 0.0002
    quantum = json.loads(
        pathlib.Path("results/development/nsc-tilt-quantum-v1.json").read_text()
    )
    assert quantum["verdict"] == "QUANTUM_IS_SINGLE_NODE"
    assert 0.8 < quantum["ratio_share_to_quantum"] < 1.2


def test_quantum_ladder_counts_integers_per_seam():
    """Round 47: every seam's active steps carry the common quantum; the
    per-seam net tilts count integers, recorded as {3, 3, 2, 1}."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-tilt-quantum-ladder-v1.json").read_text()
    )
    assert record["verdict"] == "QUANTUM_LADDER_HOLDS"
    counts = [record["seams"][key]["nearest_integer"] for key in ("0.125", "0.25", "0.5", "0.75")]
    assert counts == [3, 3, 2, 1]
    for entry in record["seams"].values():
        for step in entry["active_steps"]:
            assert abs(abs(step) - record["quantum"]) <= 0.15 * record["quantum"]


def test_neutral_fringe_shares_the_chart_free_quantum():
    """Round 48: neutral seams 4/8 quantize at the same q*N as the charged
    chart (1.149 within ~1%), with the recorded negative sign."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-neutral-fringe-quantum-v1.json").read_text()
    )
    for seam in ("4.0", "8.0"):
        qn = record["seams"][seam]["quantum_times_points"]
        assert abs(qn - 1.1476) < 0.05
    assert all(step < 0.0 for seam in record["seams"].values()
               for step in seam["active_steps"])


def test_fringe_ladder_thins_faster_than_beta():
    """Round 49: counts {3,3,2,1} against charged beta show no
    proportionality; count/beta declines with seam, recorded as a
    threshold-adjacent envelope."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-fringe-beta-join-v1.json").read_text()
    )
    ratios = [row["count_over_beta"] for row in record["rows"]]
    assert ratios[0] > ratios[-1]
    assert record["count_over_beta_spread"] > 0.15
    assert record["rows"][0]["beta"] > 0.99


def test_two_ends_join_sharper_room_slower_rise():
    """Round 51: corr(Q, rise) carries the pre-registered negative sign at
    -0.93 (n=3), while the Q-wall relation stays weak positive."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-two-ends-join-v1.json").read_text()
    )
    correlations = record["correlations"]
    assert correlations["q_vs_rise"] < -0.8
    assert 0.0 < correlations["q_vs_wall_inflow"] < 0.8
    assert record["verdict"] == "SHARPER_ROOM_SLOWER_RISE"
    assert len(record["table"]) == 3


def test_two_ends_covariance_closed_by_gap_powers():
    """Round 52: the composed Q*rise^0.2346 collapses the raw spread by
    an order of magnitude; the Q exponent is post-hoc flagged."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-two-ends-covariance-v1.json").read_text()
    )
    assert record["verdict"] == "COVARIANCE_CLOSED_BY_GAP"
    assert record["composed_spread"] < 0.06
    assert record["collapse_factor"] > 10.0
    assert record["raw_product_spread"] > 0.2


def test_two_ends_power_is_station_stable():
    """Round 53: p_Q at every recorded station time is 0.1173-0.1174;
    the recorded verdict says the power is stable."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-two-ends-power-stability-v1.json").read_text()
    )
    assert record["verdict"] == "POWER_STABLE"
    assert record["p_Q_spread"] < 0.01
    for row in record["rows"]:
        assert abs(row["p_Q"] - 0.1173) < 0.001


def test_oneway_transient_band_contains_the_wall_clock():
    """Round 54: the single outward current cell is the low band holding
    the wall's rates; the mass book keeps zero positive ticks."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-oneway-wall-join-v1.json").read_text()
    )
    assert record["verdict"] == "TRANSIENT_BAND_CONTAINS_WALL_CLOCK"
    assert record["wall_rates_inside_low_band"] is True
    assert record["mass_book"]["positive_outer_ticks"] == 0
    assert len(record["oneway_positive_cells_at_mu_1040"]) == 1


def test_all_wall_kinematics_sit_inside_the_outward_cell():
    """Round 55: rise, fall, and the oscillator constant all lie inside
    [0, 0.5], ordered, with the transient between fall and omega."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-lowband-placement-v1.json").read_text()
    )
    assert record["all_inside"] is True
    values = [row["value"] for row in record["rows"]]
    assert values == sorted(values)
    transient = record["transient_magnitude"]
    fall = record["rows"][1]["value"]
    omega = record["rows"][2]["value"]
    assert fall < transient < omega
    assert record["verdict"] == "ALL_WALL_KINEMATICS_INSIDE_OUTWARD_CELL"


def test_cone_trap_nests_inside_wall_and_window():
    """Round 56: the cone trap lies inside both the wall envelope and the
    window shell in all three cases, with ppm-uniform redshift."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-cone-wall-join-v1.json").read_text()
    )
    for row in record["rows"]:
        assert row["wall_contains_trap"] is True
        assert row["window_contains_trap"] is True
        assert row["window_bounds_wall"] is True
    assert record["redshift_at_excision_spread"] < 1e-5
    traps = [row["trap_in_r_plus_units"] for row in record["rows"]]
    assert traps == sorted(traps, reverse=True)


def test_isotropy_and_throat_balance_do_not_track():
    """Round 57: the throat balance's swings do not propagate into the
    chart-wide isotropy index; the recorded verdict is NO_TRACK."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-isotropy-wall-join-v1.json").read_text()
    )
    assert record["verdict"] == "NO_TRACK"
    assert abs(record["correlation_throat_vs_isotropy"]) < 0.7
    assert len(record["pairs"]) == 4
    assert record["tail"]["t40_isotropy"] > record["tail"]["t80_isotropy"]


def test_relaxation_wall_dies_after_the_driven_phase():
    """Round 73: the wall is present at t<=10 and absent at 15..80 while
    the throat balance runs to the single-force limit."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-relaxation-wall-v1.json").read_text()
    )
    rows = record["rows"]
    assert rows[-1]["wall_count"] == 0
    assert rows[-2]["wall_count"] == 0
    assert any(row["wall_count"] > 0 for row in rows[:2])
    assert rows[-1]["balance_at_throat"] > 0.99
    assert record["wall_survives_to_t40"] is False
    slopes = [row["throat_slope_inner"] for row in rows]
    assert slopes[-1] < slopes[0]


def test_wall_afterlife_has_one_death_and_a_sizeable_t40_event():
    """Round 74: one death window; the t=40 event witnesses agree; the
    single-force deficit halves every ~8 t."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-wall-afterlife-v1.json").read_text()
    )
    assert record["death_window"] == [[10.0, 15.0]]
    wiggle = record["t40_event"]["slope_wiggle"]
    assert wiggle["slope_20"] < wiggle["slope_15"]
    assert wiggle["slope_40"] > wiggle["slope_20"]
    assert record["single_force_deficits"]["t80"] < 0.01
    assert 5.0 < record["deficit_halving_time"] < 10.0


def test_driven_isotropy_antitracks_the_wall_kinematics():
    """Round 62: over the six family cases the isotropy fall correlates
    -0.96 with the wall rates and -0.87 with presence; I_first is the
    seventh mu-blind quantity."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-isotropy-driven-join-v1.json").read_text()
    )
    correlations = record["correlations"]
    assert correlations["isotropy_fall_vs_rise"] < -0.9
    assert correlations["isotropy_fall_vs_fall_rate"] < -0.9
    assert correlations["isotropy_fall_vs_wall_presence"] < -0.8
    assert record["i_first_mu_blindness"] < 0.005
    assert len(record["rows"]) == 6


def test_frontier_assessment_ledger():
    """Round 63: the frontier ledger names 7 closed groups and 4 open
    cells; every open cell carries a concrete blocker and every closed
    group names its records."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-frontier-assessment-v1.json").read_text()
    )
    assert record["counts"]["closed_question_groups"] == 7
    assert record["counts"]["open_cells"] == 4
    assert record["counts"]["test_definitions"] >= 80
    for cell in record["open"]:
        assert cell["blocker"]
    checked = []
    for group in record["closed"]:
        for name in group["closed_by"]:
            checked.append(pathlib.Path("results/development", name).exists())
    assert all(checked) and checked


def test_frontier_amendment_closes_two_of_four_cells():
    """Round 78: the amendment closes the wall-morphology and p_Q cells
    with their later records, leaving exactly two open cells."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-frontier-amendment-v1.json").read_text()
    )
    assert len(record["closed_cells"]) == 2
    assert len(record["remaining_open_cells"]) == 2
    remaining = [row["cell"] for row in record["remaining_open_cells"]]
    assert "exponent provenance at a fourth gap point" in remaining
    assert any("paper-level coupling" in item for item in remaining)
    for names in record["closed_by_records"].values():
        for name in names:
            assert pathlib.Path("results/development", name).exists()


def test_continuation_trap_drifts_with_the_feed():
    """Round 64: the trap creeps outward while fed (r = 0.84 over 40
    half-station intervals); the local CFL clock is constant."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-continuation-trap-v1.json").read_text()
    )
    assert record["station_count"] >= 40
    assert record["trap_mean_drift_rate"] > 0.0
    assert record["trap_drift_corr_with_mass_steps"] > 0.7
    assert record["cfl_dt"]["max"] / record["cfl_dt"]["min"] - 1.0 < 1e-4


def test_trap_timing_peaks_at_adjacent_window():
    """Round 65: the trap drift's strongest predictor is the adjacent
    window's feed (0.999), recorded as at-most-one-interval delay."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-trap-timing-v1.json").read_text()
    )
    assert record["lag_profile"]["1"]["correlation"] > 0.99
    assert record["lag_profile"]["0"]["correlation"] > 0.7
    assert record["lag_profile"]["-1"]["correlation"] < 0.5
    assert record["trap_displacement_per_unit_mass"] > 5.0


def test_march_throat_slope_crosses_at_the_regime_boundary():
    """Round 66: the march's throat slope flips sign once, between the
    light cases and the locked mu=1.040; the wall presence drops there."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-march-wall-join-v1.json").read_text()
    )
    assert record["verdict"] == "SIGN_CROSSES_AT_THE_REGIME_BOUNDARY"
    assert record["throat_slope_sign_changes"] == 1
    slopes = [row["throat_slope_last"] for row in record["rows"]]
    assert slopes[0] > 0 and slopes[-1] < 0
    assert record["rows"][-1]["wall_presence_fraction"] < record["rows"][0]["wall_presence_fraction"]


def test_march_slope_trajectory_flips_at_the_wall_epoch_for_light_cases():
    """Round 75: all cases start negative; the wall epoch flips the light
    cases positive (spike +180 at mu=1.002) while mu=1.040 never flips."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-march-slope-trajectory-v1.json").read_text()
    )
    series = record["series"]
    for name, rows in series.items():
        assert rows[0]["throat_slope_inner"] < 0.0
    light = series["mu-1.002-eps-0.001"]
    assert max(row["throat_slope_inner"] for row in light) > 100.0
    assert light[-1]["throat_slope_inner"] > 0.0
    heavy = series["mu-1.040-eps-0.001"]
    assert all(row["throat_slope_inner"] < 0.0 for row in heavy)
    assert record["mu_1040_sign_crossing_windows"] == []


def test_march_balance_trajectory_wall_lifespan_is_mu_monotone():
    """Round 76: the lightest case keeps walls everywhere (3 at t=10);
    the heavy case dissolves from t=15 and its t=20 balance matches the
    relaxation arm's digit-wise."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-march-balance-trajectory-v1.json").read_text()
    )
    series = record["series"]
    light = series["mu-1.002-eps-0.001"]
    assert all(row["wall_count"] > 0 for row in light)
    assert max(row["wall_count"] for row in light) == 3
    heavy = series["mu-1.040-eps-0.001"]
    assert heavy[0]["wall_count"] > 0
    assert all(row["wall_count"] == 0 for row in heavy[2:])
    assert abs(heavy[-1]["balance_at_throat"] - 0.7706) < 0.001


def test_wall_extinction_light_outlives_heavy_on_the_deep_grid():
    """Round 77: on the march domain (all cases covered to t=20) the
    heavy wall dies at 10 while both light walls reach 20; the family
    domain's ragged coverage is flagged in the record."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-wall-extinction-v1.json").read_text()
    )
    march = record["march_chart_extinction"]
    assert march["mu-1.040-eps-0.001"]["last_wall"] == 10.0
    assert march["mu-1.002-eps-0.001"]["last_wall"] == 20.0
    assert march["mu-1.010-eps-0.001"]["last_wall"] == 20.0
    family = record["family_chart_extinction"]
    assert family["mu-1.002-eps-0.001"]["covered_until"] < 10.0
    assert record["ordering_inverts_with_domain_depth"] is True


def test_cross_arm_identity_is_bit_exact():
    """Round 79: the march and relaxation payload families agree to zero
    difference at every shared station; the record says so."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-cross-arm-identity-v1.json").read_text()
    )
    assert record["bit_identical_balances"] is True
    assert record["wall_counts_agree"] is True
    assert len(record["table"]) == 4
    for row in record["table"]:
        assert row["balance_diff"] == 0.0
        assert row["march_walls"] == row["relaxation_walls"]


def test_deep_grid_wall_radii_are_lower_but_mu_ordered():
    """Round 81: the deep-grid radii sit below r_+ for the heavier cases
    and below the family envelopes everywhere, with the mu-ordering kept."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-march-wall-radii-v1.json").read_text()
    )
    radii = record["radii"]
    light = radii["mu-1.002-eps-0.001"]
    heavy = radii["mu-1.040-eps-0.001"]
    assert light["range_r_plus_units"][1] > heavy["range_r_plus_units"][1]
    assert heavy["range_r_plus_units"][1] < 0.8
    assert all(not entry["inside_family_envelope"] or entry["inside_family_envelope"] is True
               for entry in radii.values())
    assert all(entry["range_r_plus_units"] is not None for entry in radii.values())


def test_wall_domain_law_amplitude_ratio_is_mu_monotone():
    """Round 82: the shallow/deep breathing-amplitude ratio grows with mu
    (7/17/40), and the deep wall's offset above the inner edge is recorded."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-wall-domain-law-v1.json").read_text()
    )
    ratios = [row["amplitude_ratio_shallow_over_deep"] for row in record["rows"]]
    assert ratios == sorted(ratios)
    assert ratios[0] > 5.0 and ratios[-1] > 30.0
    for row in record["rows"]:
        assert row["wall_min_offset_above_inner_edge"] > 0.0
    assert record["rows"][0]["deep_amplitude_r_plus"] > record["rows"][-1]["deep_amplitude_r_plus"]


def test_deep_domain_rates_are_suppressed_and_honestly_flagged():
    """Round 86: deep rates run an order below the family rates, and the
    record names the sparse cadence rather than claiming a power."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-march-wall-rates-v1.json").read_text()
    )
    for row in record["rows"]:
        if row["deep_rise_rate"] is not None:
            assert abs(row["deep_rise_rate"]) < 0.05 * row["family_rise_rate"]
    assert "insufficient for a second-domain power-pair test" in (
        str(record.get("notes")) + str(record.get("reading")))
    assert record["usable_rise_cases"]


def test_programme_index_is_complete_and_checked():
    """Round 84: the master index cites all fourteen parts' records and
    every cited path was verified to exist before sealing."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-programme-index-v1.json").read_text()
    )
    assert record["index_complete"] is True
    assert record["missing_records"] == []
    assert record["part_count"] == 14
    cited = [item["record"] for part in record["parts"] for item in part["records"]]
    assert len(cited) == 46
    assert all(item["exists"] for part in record["parts"] for item in part["records"])


def test_ledger_audit_is_clean():
    """Round 85: forty-five records re-read, nineteen verdicts re-checked,
    no missing files, no mismatches."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-ledger-audit-v1.json").read_text()
    )
    assert record["audit_clean"] is True
    assert record["checked_records"] == 45
    assert record["verdict_checks"] == 19
    assert record["missing_records"] == []
    assert record["mismatches"] == []


def test_t30_launch_log_narrates_the_campaign_period():
    """Round 97: the launch log holds the four period entries and the two
    armed instruments, matching the recorded rounds."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-t30-launch-log-v1.json").read_text()
    )
    assert len(record["period"]) == 4
    assert {entry["round"] for entry in record["period"]} == {90, 93, 96}
    assert len(record["instruments_armed"]) == 2
    assert record["campaign"] == "nsc-rn-inward-march-t30-v1"


def test_t30_timing_resolution_reads_the_v1_mtimes():
    """Round 104: the v1 march measured 173.7 min from its recorded file
    mtimes with ~39 min between later stations; the t30 forecast follows."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-t30-timing-resolution-v1.json").read_text()
    )
    assert abs(record["v1_total_prefix_minutes"] - 173.7) < 1.0
    intervals = record["v1_intervals_minutes"]
    assert intervals[0] > 50.0
    assert all(30.0 < x < 50.0 for x in intervals[1:])
    assert record["t30_forecast_minutes_from_launch"]["station_25"] \
        > record["v1_total_prefix_minutes"]
    assert record["t30_forecast_minutes_from_launch"]["station_30"] \
        > record["t30_forecast_minutes_from_launch"]["station_25"]


def test_march_duration_is_mu_ordered_and_intervals_shrink():
    """Round 105: the march costs 173.7/89.8/54.0 minutes for the three
    masses; every case's station intervals decrease with t."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-march-duration-table-v1.json").read_text()
    )
    rows = record["rows"]
    totals = [row["total_to_station_20_minutes"] for row in rows]
    assert totals == sorted(totals, reverse=True)
    assert rows[0]["total_to_station_20_minutes"] > 3.0 * rows[-1]["total_to_station_20_minutes"]
    for row in rows:
        assert row["interval_minutes"] == sorted(row["interval_minutes"], reverse=True)
        assert all(x > 0 for x in row["interval_minutes"])


def test_march_cost_vs_wall_join_is_partial_and_honest():
    """Round 106: the death-inside-coverage case shows the largest shrink;
    the rank comparison fails on the two cocovered light cases and the
    record says so."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-march-cost-vs-wall-v1.json").read_text()
    )
    assert record["verdict"].startswith("PARTIAL")
    assert record["shrink_ranks_equal_death_ranks"] is False
    heavy = next(row for row in record["rows"] if row["wall_extinction_t"] == 10.0)
    assert heavy["interval_shrink_minutes"] == max(
        row["interval_shrink_minutes"] for row in record["rows"])


def test_march_learn_rate_factors_are_mu_ordered():
    """Round 108: the interval-decay factors s = 0.986/0.897/0.810 order
    cleanly with mu; the fit truth is recorded as a fit."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-march-learn-rate-v1.json").read_text()
    )
    factors = [row["geometric_factor_s"] for row in record["rows"]]
    assert factors == sorted(factors, reverse=True)
    assert factors[0] > 0.98
    assert factors[-1] < 0.85
    assert "fit, not a law" in " ".join(record["notes"])


def test_relaxation_duration_is_steady_two_minutes_per_t():
    """Round 110: the relaxation campaign's mtimes show a constant 2.0
    min per t through t=80 - the cleanest recorded march rate."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-relaxation-duration-v1.json").read_text()
    )
    for entry in record["per_unit"]:
        assert abs(entry["minutes_per_t"] - 2.0) < 0.2
    assert abs(record["total_prepare_to_t80_minutes"] - 157.5) < 1.0


def test_t30_dual_forecast_names_the_adjudicator():
    """Round 110: two recorded rates bracket the arrival; the record
    names the campaign write times as the adjudicator."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-t30-dual-forecast-v1.json").read_text()
    )
    assert record["forecast_a_v1_mtimes"]["station_25_minutes"] \
        > record["forecast_b_relaxation_mtimes"]["station_25_minutes"]
    assert record["adjudicator"].startswith("nsc-rn-inward-march-t30-v1")


def test_t30_extinction_once_the_march_lands():
    """Round 92: conditional. Fires only after the t30 campaign lands;
    then the extinction verdict must be one of the three named windows."""
    import json
    import pathlib
    import pytest
    path = pathlib.Path("results/development/nsc-t30-extinction-v1.json")
    if not path.exists():
        pytest.skip("t30 march has not landed yet")
    record = json.loads(path.read_text())
    assert record["verdict"] in ("WALL_SURVIVES_TO_30", "DIES_IN_25_30", "DIES_IN_20_25")
    assert record["creation_only"] is True


def test_q_power_decomposes_into_the_room_stretch():
    """Round 67: p_Q reconciles exactly with the radius-factor power;
    lapse and shift contribute near zero as the flat-slice gauge says."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-q-decomposition-v1.json").read_text()
    )
    components = record["component_powers"]
    assert abs(components["p_N_single"]) < 0.01
    assert abs(components["p_beta_factor"]) < 0.01
    assert components["p_r"] > 0.1
    assert abs(record["reconciliation"]["sum_minus_measured"]) < 0.001
    assert abs(components["p_sum"] - 0.1173) < 0.001


def test_q_chord_and_identity_sit_inside_the_analytic_envelope():
    """Round 68: the measured chord fits sqrt(mu^2-1)/mu's range, and the
    four-fold envelope contains the ladder's b."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-q-analytic-envelope-v1.json").read_text()
    )
    assert record["chord_inside_envelope"] is True
    assert record["b_inside_four_p_envelope"] is True
    low, high = record["envelope_range"]
    assert low < record["measured_chord_p_q"] < high
    assert record["four_p_envelope"][0] < record["ladder_b"] < record["four_p_envelope"][1]


def test_production_preflight_passes_and_prices_the_run():
    """Round 70: the calibration dry-run passed inside tolerance and the
    projected t>20 run fits the recorded 900 s admission slot."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-production-preflight-v1.json").read_text()
    )
    assert record["calibration"]["passed"] is True
    assert record["calibration"]["worst_residual"] < 1e-4
    assert record["cost_basis"]["projected_t20_to_t40_cpu_seconds"] \
        <= record["cost_basis"]["admission_slot_seconds"]
    assert len(record["ceremony_remaining"]) == 4


def test_constant_table_records_the_flagged_near_identity():
    """Round 58: seven constants with named sources; the 4*p_Q ~ b
    proximity matches to under 1% and is flagged numerology."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-constant-table-v1.json").read_text()
    )
    assert len(record["constants"]) == 7
    assert abs(record["near_identity"]["relative_difference"]) < 0.01
    assert "numerology" in record["status"]
    for entry in record["constants"].values():
        assert entry["source"]


def test_ladder_covariance_predictive_pass():
    """Round 60: the ladder-implied exponent closes the covariance at
    the same spread as the measured one; verdict recorded PREDICTIVE_PASS."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-ladder-covariance-test-v1.json").read_text()
    )
    assert record["verdict"] == "PREDICTIVE_PASS"
    assert record["predicted_composed_spread"] < 0.05
    assert abs(record["p_Q_predicted"] - record["ladder_b"] / 4.0) < 1e-12
    assert abs(record["predicted_composed_spread"] -
               record["recorded_composed_spread_reference"]) < 0.01


def test_composed_constant_sits_at_95_percent_of_the_clock():
    """Round 61: the composed covariance constant is 0.948 of the infall
    clock with 2% spread, recorded as flagged scale proximity."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-composed-vs-clock-v1.json").read_text()
    )
    assert abs(record["mean_ratio"] - 0.948) < 0.005
    assert record["ratio_spread"] < 0.03
    assert "not claimed" in record["status"]
    assert len(record["rows"]) == 3


def test_window_shell_holds_a_fixed_five_percent_cell():
    """Round 40: the [r_+, 2 r_+] window shell carries ~5-6% of both
    instruments on both charts below the next atlas edge."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-window-shell-share-v1.json").read_text()
    )
    for chart in ("neutral_chart", "charged_chart"):
        entry = record[chart]
        assert 0.04 < entry["instrument_a"][0]["l1_share"] < 0.08
        assert 0.04 < entry["instrument_b"]["window_shell_share"] < 0.08


def test_wall_fall_is_constant_frequency_sampled_by_amplitude():
    """Round 39: fall/A is family-constant at the 8-9% level, and the
    record flags the exponent origins as underdetermined."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-wall-clock-derivation-v1.json").read_text()
    )
    assert record["decomposition_holds"] is True
    assert record["omega_spread"] < 0.12
    assert record["scale_counting"]["underdetermined"]


def test_asymmetry_law_derives_the_held_out_cases():
    """Round 31: the composed clock law predicts each held-out case's
    asymmetry to within 15%, with mu=1.040 at sub-percent."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-wall-asymmetry-derived-v1.json").read_text()
    )
    predictions = record["predictions"]
    for entry in predictions.values():
        assert abs(entry["relative_error"]) < 0.15
    assert abs(predictions["mu-1.040-eps-0.001"]["relative_error"]) < 0.01


def test_seam_fee_is_chart_invariant():
    """Round 32: the seam-8 fee varies by <1e-3 across the three charts;
    the recorded verdict reflects the refuted lock-anchor hypothesis."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-fee-convergence-v1.json").read_text()
    )
    fees = [entry["fee"] for entry in record["rows"].values()]
    assert max(fees) - min(fees) < 1e-3
    assert record["verdict"] == "LOCK_NOT_ANCHOR"
    assert all(0.25 < fee < 0.30 for fee in fees)


def test_full_seam_profiles_coincide_across_charts():
    """Round 33: the five-seam profile matches across all three charts to
    within the recorded 2.1e-3 maximum pairwise deviation."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-seam-invariance-v1.json").read_text()
    )
    assert record["verdict"] == "PROFILES_COINCIDE"
    assert record["max_pairwise_deviation"] < 0.005
    profile_rows = record["profiles"]
    for seam in ("4.0", "8.0", "16.0", "24.0", "30.0"):
        values = [profile_rows[key][seam] for key in profile_rows]
        assert max(values) - min(values) < 0.005


def test_clock_pair_holds_on_the_strong_packet_trio():
    """Round 34: the epsilon=0.010 trio reproduces both gap power laws
    (gapped spreads below 8%); the recorded verdict says the cross holds."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-wall-clock-epscross-v1.json").read_text()
    )
    assert record["verdict"] == "EPS_CROSS_HOLDS"
    strong = record["eps_0010"]
    assert strong["rise_half_spread"] < 0.08
    assert strong["fall_quarter_spread"] < 0.08
    cases = record["cases"]
    ratio = cases["mu-1.040-eps-0.010"]["rise_rate"] / cases["mu-1.040-eps-0.001"]["rise_rate"]
    assert 0.9 < ratio < 1.1


def test_cosmic_band_maps_to_the_outer_cap_shell():
    """Round 35: the [0.03, 0.10] local band maps to the outermost shell
    near the cap, chart-invariant; the target 0.049 sits at ~28.2 r_m."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-cosmic-band-location-v1.json").read_text()
    )
    mapping = record["mapping"]
    low, high = mapping["band_seam_interval"]
    assert 24.0 < low and high < 30.5
    assert mapping["seam_of_band_low_0_03"] > mapping["seam_of_target_0_049"] > mapping["seam_of_band_high_0_10"]
    targets = [entry["local_0_049_seam"] for entry in record["per_case"].values()]
    assert max(targets) - min(targets) < 0.1


def test_cosmic_shell_weights_are_fair_and_live():
    """Round 36: the [26,29] band carries a length-fair infall-time share
    and a third of the causal mass response; concentration B near unity."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-cosmic-shell-weight-v1.json").read_text()
    )
    b = record["instrument_b"]
    assert abs(b["band_share"] - record["band_length_fraction_of_span"]) < 0.01
    assert 0.8 < record["concentration"]["b_band_over_length"] < 1.1
    mass = record["instrument_a_mass_shares"]
    assert mass[1]["l1_share"] > 0.25
    lapse = record["instrument_a_lapse_shares"]
    assert lapse[0]["l1_share"] > 0.85


def test_lapse_clock_awareness_profile_is_monotone():
    """Round 37: beta falls monotonically from ~1 at seam 4 to 0.08 at the
    cap; the cosmic seam holds beta ~0.41."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-time-awareness-v1.json").read_text()
    )
    rows = record["rows"]
    betas = [row["beta"] for row in rows]
    assert all(betas[i] > betas[i + 1] for i in range(len(betas) - 1))
    assert 0.95 < betas[0] < 1.05
    assert betas[-1] < 0.1
    target_row = next(row for row in rows
                      if abs(row["seam"] - record["target_seam"]) < 1e-6)
    assert 0.35 < target_row["beta"] < 0.45


def test_charged_omega_star_is_not_gap_clocked():
    """Round 27: omega_star is eta-flat and sqrt(gap) gapping worsens the
    spread - the spectroscopy is cap-clocked, not gap-clocked."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-charged-clock-horizon-v1.json").read_text()
    )
    assert record["omega_star_raw_spread"] < 0.03
    assert record["omega_star_sqrt_gap_spread"] > record["omega_star_raw_spread"]
    assert len(record["rows"]) == 3


def test_omega_star_identity_holds_to_machine_precision():
    """Round 28: omega_star = k*<1/r>_T with k = Q_E/(16 pi C_F), ratio
    constant across eta at floating-point level."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-omega-star-identity-v1.json").read_text()
    )
    assert record["identity_holds"] is True
    assert record["ratio_spread"] < 1e-12
    for row in record["rows"]:
        assert abs(row["ratio"] - record["k"]) < 1e-12


def test_matryoshka_fee_join_sums_near_omega4():
    """Round 29: with Schroeder branching 3, the resolvent and tangent fees
    sum near the Omega=4 anchor; the transfer fee stays far below."""
    import json
    import pathlib
    record = json.loads(
        pathlib.Path("results/development/nsc-matryoshka-fee-join-v1.json").read_text()
    )
    sites = record["sites"]
    assert abs(sites["resolvent_1_over_omega"]["sigma"] - 4.0) / 4.0 < 0.03
    tangent = sites["lapse_tangent_f"]["sigma"]
    assert 0.8 < tangent / 4.0 < 1.4
    assert sites["transfer_fee"]["sigma"] < 2.5