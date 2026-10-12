"""Same-source inventory, constrained preparation, and coarse source effects."""
import json
import numpy as np

import derive_nsc_spherical_conformal_episode as episode
import derive_nsc_spherical_conformal_source_controls as source


def test_parameter_changes_only_regional_pair_weights_and_preserves_total_inventory():
    for alpha in (0., -1., .95, 1., 1.05):
        weights = source.source_weights(alpha)
        assert abs(np.sum(weights) - 3) < 1e-14
        assert np.array_equal(weights[::2], weights[1::2])
        assert np.all((weights >= 0) & (weights <= 1))
    assert np.array_equal(source.source_weights(1), episode.galerkin.OCCUPATIONS)
    assert np.array_equal(source.source_weights(-1), source.source_weights(1)[::-1])


def test_saved_controls_reuse_columns_and_keep_the_original_witness_bytes():
    record = json.loads(source.OUT.read_text())
    assert record["sealed_sources_unchanged"] is True
    assert record["source_bindings_before"]["v1_source_npz"] == episode.sha256(episode.NPZ)
    assert record["mean_current_deleted"] is False
    assert record["radius_imposed"] is False
    assert record["midtrajectory_reset"] is False
    assert record["all_target_reached"] is True
    assert record["cpu_seconds"] < source.CPU_BUDGET
    with np.load(source.NPZ, allow_pickle=False) as saved, np.load(episode.NPZ, allow_pickle=False) as baseline:
        assert np.array_equal(saved["original_T0_phi0"], baseline["nf256_initial_phi0"])
        assert np.array_equal(saved["original_T0_phi1"], baseline["nf256_initial_phi1"])
        for label, alpha in source.RUN_PLAN:
            for phi in ("phi0", "phi1"):
                assert np.array_equal(saved[label + "_initial_" + phi], baseline["nf256_initial_" + phi])
            assert np.array_equal(saved[label + "_weights"], source.source_weights(alpha))
            report = record["results"][label]
            constructor = report["constructor"]
            assert constructor["source_derived"] is True
            assert constructor["columns_unchanged"] is True
            assert constructor["current_unchanged"] is True
            assert constructor["mean_retained_in_shift_residual"] is True
            assert constructor["source_mean_subtracted"] is False
            assert report["positive_chart"] is True
            assert report["gram_gap_max"] < 1e-8
            assert report["constraint_end_margin"] > 0
            assert report["observable_error_certified"] is False


def test_signed_surface_and_forcing_integrals_recompute_from_saved_times():
    record = json.loads(source.OUT.read_text())
    keys = list(record["series_schema"])
    window_keys = list(record["window_schema"])
    with np.load(source.NPZ, allow_pickle=False) as saved:
        for label, alpha in source.RUN_PLAN:
            series = saved[label + "_series"]
            times = series[:, keys.index("time")]
            forcing = series[:, keys.index("forcing_norm")]
            integral = float(np.sum(.5 * np.diff(times) * (forcing[:-1] + forcing[1:])))
            reported = record["results"][label]["integrals"]["forcing_norm"]["trapezoid"]
            assert abs(integral - reported) < 1e-14
            physical_times = saved[label + "_physical_times"]
            windows = saved[label + "_window_series"]
            for index in range(4):
                left = windows[:, index, window_keys.index("normal_left_outward_flux")]
                right = windows[:, index, window_keys.index("normal_right_outward_flux")]
                net = windows[:, index, window_keys.index("normal_boundary_flux")]
                assert np.max(abs(net + left + right)) < 1e-12
                for name, values in (("normal_left_outward_flux", left), ("normal_right_outward_flux", right)):
                    expected = np.sum(.5 * np.diff(physical_times) * (values[:-1] + values[1:]))
                    assert abs(expected - record["results"][label]["surfaces"][index][name]["signed_integral"]) < 1e-14
    reverse = record["results"]["reversed_alpha_minus1"]
    assert reverse["physical"]["final"]["shell_leader"] == 1
    uniform = record["results"]["uniform_alpha0"]
    assert uniform["surfaces"][0]["normal_left_outward_flux"]["absolute_integral"] > .01
    assert uniform["surfaces"][1]["normal_right_outward_flux"]["absolute_integral"] > .1
