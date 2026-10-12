"""Independent metric controls and saved-sample checks. No evolution."""
import hashlib
import json
import os
import subprocess
import zlib
from pathlib import Path

for _thread_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_var, "1")

import numpy as np
import pytest

from recursive_horizons.nsc_spherical_coupling import CALIBRATION, PERIOD, periodic_derivative
from recursive_horizons.provenance import resolve_pinned_source_bytes
from recursive_horizons.nsc_spherical_null_expansion import (
    CPU_BUDGET_S,
    EPISODE_JSON,
    EPISODE_NPZ,
    JSON_BYTE_LIMIT,
    RECORD_PATH,
    boosted_expansions,
    build_record,
    circular_intervals,
    expansions_from_stored_normal,
    null_expansions,
    spectral_derivative,
    static_clock,
)

_MODULE = Path(__file__).resolve().parents[1] / "src" / "recursive_horizons" / "nsc_spherical_null_expansion.py"
_REVIEW = Path(__file__).resolve().parents[1] / "results" / "development" / "nsc-local-boundary-review-v1.json"
_REVIEW_V2 = Path(__file__).resolve().parents[1] / "results" / "development" / "nsc-local-boundary-review-v2.json"
_SOURCE_HISTORY = Path(__file__).resolve().parents[1] / ".source-history"
_V2_SOURCE_COMMIT = "5f10ecd365843d1616e50eb16a20d7acd8377e2c"
_V1_REVIEW_SHA256 = "0d43accbc5d46407c094890a9722c1e78842224313580044df0d4b641dff8d9b"
_OLD_MODULE_SHA256 = "1aa47f0b3d680621d23ac9329e0d9e1e38070b7df42fc2ec48848ea367fd82d5"
_SEALED_INDEPENDENT_SHA256 = "68c2fa1559d705c3aae2d3b64d693af72311802e51cb5914432e1422063b4ea9"
_STORED_OBSERVER = "initial region-0 plus and minus source columns, unphased"
_LATER_INDEPENDENT_TESTS = (
    "test_control_policy_keeps_g_initial_and_proxies_off_the_gate",
    "test_v2_reads_stored_endpoint_without_touching_v1",
)
_ORIGINAL_SOURCE_OBJECTS = (
    ("09f3faa364ca88b8ff482c413022ce38f40bef0d", "commit", 339, "f17ca2a8a430a11161f005be9dac2178f391565c884d3a4475abcee7db6f1f98"),
    ("6be0057f68cc1ca3fc021b0a2865ad87ea79f5c1", "tree", 847, "62fa898de8369496dd1ef1df07113bee809954000a74e71dbe7534176d16ec08"),
    ("22993e882ceb7032db8b233fb71231ee3b226316", "tree", 45, "08ee451400f015721b2ef91e55dcfc0f94eadab268d1618bc4fd6ea4bc93a6a7"),
    ("61286d365dbda4cf1e19b8f154d6f2a298bf92de", "tree", 9458, "2ef8afced2e0123b7a1245674fd78af1f58456751fcd544c3ea509e0539ee0a4"),
    ("a5a85acef3d98a8d6ff66f928a30725856a1ac6c", "blob", 7350, "d6f31c307705b49371242f04c4cd44ae326593b068904cd3d71506751b0a42d9"),
    ("750f08abd50ef461600baf9dd319d1ac38b50ca3", "commit", 317, "6b9db11f0a897997fd1532d39e5062584dafc601376b7f89eb3b591c2192ebb8"),
    ("157a55cdf6ceaf2f69350914046c54fe059855e1", "tree", 882, "8532e5f38ada1e2c1312d753781d953850512c5ba7661034c7fdde06b5c64f87"),
    ("f7b6d43cae388b9e87a586e1bc87ee62327ab96d", "blob", 1262, "935a7d3b693692d5e91d3710934ec68af4e3a532dfba6e905fc5985b89b356c1"),
    ("c7dc0f98cacb28afd030b089bdb9e1369288120c", "commit", 334, "b08c61f8d20dde6afb8eea725ba7ed4ad3d012ee509004e4b1593298316fbadb"),
    ("201a906a41d70ca033fb19f4bf7f1483f1fa2552", "tree", 902, "529fb8a47bc619d56cc636a3c3a1392020bc22c4973edf80bff6b8167636e9cc"),
    ("d778878f61150cfa0dd2cd3b4024aea276001d46", "tree", 45, "0d60fb8d81e087c2e4c67afdd49c86fcf9660dcad1d329ae65d4365e465d7abb"),
    ("4dc21108424b29558c9092c703184543721f5248", "tree", 11006, "db08f2b7c71150db61d7f4f08b3c9760e15fc05e484e5007355430525c0341ff"),
    ("31818fdd6b8cc9c4d8b9d33f41e0f241b2dd1227", "blob", 10229, "42cc1f09dd49d9c3ee0eeab22fe3b778f0562729e05414e82e0ef19658c82971"),
    ("e62d28f3230ae6ddb2712df71eb7c61a8ed65ad6", "commit", 341, "6b51fb912f469919df458c2a14b97c8a5e2ce69ff7f8929126cfae19a98d38ef"),
    ("d482613c0af362b6b7440788b335b136e24b8b4f", "tree", 847, "a123fa73f023a45a743a9a7d0712d28038e090a7118cccb0d746b948e0867e5a"),
    ("f2556bf248ef23192a4c864cb816090cb8641931", "tree", 45, "675dfff3c855ac4228ad398fbf28962e7b169d2fb770f0fb5d38b6ff1f38d23b"),
    ("826473ae9a2ff7a3da07718738ae48be7badcd68", "tree", 9798, "c360866cd2d48b03a62f39b9f376a619c24596cd2169862764988ceb2e813e99"),
    ("9efa09bd82b7003391890f1a19828de04148173f", "blob", 9693, "dfb09f21ba206c4b7272d98dc3270ae902f6ade26c2b203596723c0e1b218cfd"),
)
_ORIGINAL_PINS = (
    ("09f3faa364ca88b8ff482c413022ce38f40bef0d", "src/recursive_horizons/nsc_local_incoming_family.py", "a5a85acef3d98a8d6ff66f928a30725856a1ac6c", "9d2b9fc57e2057d6c34d8943384428341b314370021bcd4156af88320c3de04b"),
    ("750f08abd50ef461600baf9dd319d1ac38b50ca3", "pyproject.toml", "f7b6d43cae388b9e87a586e1bc87ee62327ab96d", "b5bfe97f7225d19d6036e05d3199376885d4d626db05d49013495c94ec3e400e"),
    ("c7dc0f98cacb28afd030b089bdb9e1369288120c", "src/recursive_horizons/nsc_ks_energy_propagator.py", "31818fdd6b8cc9c4d8b9d33f41e0f241b2dd1227", "5f75023efc8c2b637fc4983b414cab50f15f90f8abc1a24dbc1d23ecef9ecf94"),
    ("e62d28f3230ae6ddb2712df71eb7c61a8ed65ad6", "src/recursive_horizons/nsc_ks_difference_envelope.py", "9efa09bd82b7003391890f1a19828de04148173f", "77c89c8e720ca97336ccdee441a645b5c73472cb6f15b09e8120f68ab53db196"),
)
_HISTORICAL_PIN_PATHS = {
    "src/recursive_horizons/nsc_spherical_null_expansion.py": _OLD_MODULE_SHA256,
    "tests/test_nsc_local_boundary_independent.py": _SEALED_INDEPENDENT_SHA256,
}
_RESPONSE_JSON = Path(__file__).resolve().parents[1] / "results" / "development" / "nsc-coupled-local-response-v1.json"
_RESPONSE_NPZ = Path(__file__).resolve().parents[1] / "results" / "development" / "nsc-coupled-local-response-v1.npz"
_RESPONSE_MODULE = Path(__file__).resolve().parents[1] / "src" / "recursive_horizons" / "nsc_coupled_local_response.py"
_NULL_DOC = Path(__file__).resolve().parents[1] / "docs" / "nsc-spherical-null-expansion.md"
_RESPONSE_DOC = Path(__file__).resolve().parents[1] / "docs" / "nsc-coupled-local-response.md"


@pytest.fixture(scope="module")
def measurement():
    return build_record(source_ref=_V2_SOURCE_COMMIT)


@pytest.fixture(scope="module")
def saved_record():
    return json.loads(RECORD_PATH.read_text())


def _independent_inverse(r, length_density, radial_factor, shift):
    """2x2 inverse from det(g), written apart from the module's inverse helper."""
    area = r * r
    conformal_square = radial_factor * radial_factor
    g_tt = area * (length_density ** 2 - conformal_square * shift ** 2)
    g_tx = -area * conformal_square * shift
    g_xx = -area * conformal_square
    determinant = g_tt * g_xx - g_tx ** 2
    inverse_tt = g_xx / determinant
    inverse_xx = g_tt / determinant
    inverse_tx = -g_tx / determinant
    return g_tt, g_tx, g_xx, inverse_tt, inverse_tx, inverse_xx


def _independent_norms(r, length_density, radial_factor, shift):
    g_tt, g_tx, g_xx, _inverse_tt, _inverse_tx, _inverse_xx = _independent_inverse(
        r, length_density, radial_factor, shift
    )
    lapse = r * length_density
    radial_metric = r * radial_factor
    normal_t = 1.0 / lapse
    normal_x = -shift / lapse
    radial_x = 1.0 / radial_metric
    plus_t, plus_x = normal_t, normal_x + radial_x
    minus_t, minus_x = normal_t, normal_x - radial_x

    def form(left_t, left_x, right_t, right_x):
        return (
            g_tt * left_t * right_t
            + g_tx * (left_t * right_x + left_x * right_t)
            + g_xx * left_x * right_x
        )

    return {
        "nn": form(normal_t, normal_x, normal_t, normal_x),
        "ee": form(0.0, radial_x, 0.0, radial_x),
        "plus": form(plus_t, plus_x, plus_t, plus_x),
        "minus": form(minus_t, minus_x, minus_t, minus_x),
        "cross": form(plus_t, plus_x, minus_t, minus_x),
    }


def test_postprocess_does_not_own_evolution():
    text = _MODULE.read_text()
    for forbidden in ("rk4_step", "evolve_episode", "solve_initial_radius", "full_source_rho"):
        assert forbidden not in text
    assert "frame_chi" not in text
    assert "weyl" not in text


def test_inverse_metric_product_null_norms_and_boost():
    cases = (
        (2.0, 1.0, 1.0, 0.0, -0.5, 0.0, -0.25, -0.25),
        (2.0, 1.0, 1.0, 0.0, 0.5, 0.0, 0.25, 0.25),
        (2.0, 1.0, 1.0, 0.0, 0.0, 0.5, 0.25, -0.25),
        (2.0, 3.0, 4.0, 0.5, -1.25, 0.25, None, None),
    )
    generator = np.random.default_rng(19)
    extra = [
        tuple(float(value) for value in row) + (None, None)
        for row in generator.uniform(
            [0.3, 0.2, 0.2, -1.0, -2.0, -2.0],
            [3.0, 2.0, 2.0, 2.0, 2.0, 2.0],
            size=(6, 6),
        )
    ]
    for radius, lapse_factor, conformal, shift, r_t, r_x, expected_plus, expected_minus in list(cases) + extra:
        theta_plus, theta_minus, product = null_expansions(
            radius, r_t, r_x, lapse_factor, conformal, shift
        )
        _g_tt, _g_tx, _g_xx, inverse_tt, inverse_tx, inverse_xx = _independent_inverse(
            radius, lapse_factor, conformal, shift
        )
        gradient = inverse_tt * r_t ** 2 + 2.0 * inverse_tx * r_t * r_x + inverse_xx * r_x ** 2
        assert abs(product - 4.0 / radius ** 2 * gradient) < 1e-12
        norms = _independent_norms(radius, lapse_factor, conformal, shift)
        assert abs(norms["nn"] - 1.0) < 1e-12
        assert abs(norms["ee"] + 1.0) < 1e-12
        assert abs(norms["plus"]) < 1e-12
        assert abs(norms["minus"]) < 1e-12
        assert abs(norms["cross"] - 2.0) < 1e-12
        for factor in (0.5, 3.0, 8.0):
            scaled_plus, scaled_minus, scaled_product = boosted_expansions(
                theta_plus, theta_minus, factor
            )
            assert abs(scaled_product - product) < 1e-12
            if theta_plus != 0.0:
                assert np.sign(scaled_plus) == np.sign(theta_plus)
            if theta_minus != 0.0:
                assert np.sign(scaled_minus) == np.sign(theta_minus)
        if expected_plus is not None:
            assert abs(float(theta_plus) - expected_plus) < 1e-12
            assert abs(float(theta_minus) - expected_minus) < 1e-12
        flipped_product = (-theta_plus) * (-theta_minus)
        assert abs(flipped_product - product) < 1e-15
        if theta_plus != 0.0:
            assert (-theta_plus) * theta_plus < 0.0


def test_derivative_matches_owner_and_kills_nyquist():
    points = 64
    dense = periodic_derivative(points, PERIOD)
    coordinate = np.arange(points) * (PERIOD / points)
    field = np.sin(4.0 * np.pi * coordinate / PERIOD) + 0.3 * np.cos(np.pi * np.arange(points))
    assert np.max(np.abs(spectral_derivative(field, PERIOD) - dense @ field)) < 1e-11
    nyquist = np.cos(np.pi * np.arange(points))
    assert np.max(np.abs(spectral_derivative(nyquist, PERIOD))) < 1e-12


def test_clock_arc_uses_calibration_literals():
    _coordinate, length_density, shift = static_clock(8)
    assert length_density[0] == pytest.approx(CALIBRATION["b0"])
    assert shift[0] == pytest.approx(CALIBRATION["beta0"])
    assert length_density[2] == pytest.approx(CALIBRATION["b0"] * CALIBRATION["Omega"] ** 2)
    assert shift[2] == pytest.approx(CALIBRATION["beta0"] * CALIBRATION["Omega"] ** 2)


def test_circular_intervals_wrap_and_split():
    assert circular_intervals(np.array([False, False, False])) == []
    split = circular_intervals(np.array([False, True, True, False, True]))
    assert [(item["start_index"], item["end_index"], item["count"]) for item in split] == [
        (1, 2, 2),
        (4, 4, 1),
    ]
    wrapped = circular_intervals(np.array([True, True, False, True]))
    assert wrapped[0]["wraps"] is True
    assert wrapped[0]["count"] == 3


def test_saved_cases_change_character_under_positive_chart(measurement, saved_record):
    assert measurement["evolution_performed"] is False
    assert measurement["inheritance_recomputed"] is False
    assert measurement["chi_used"] is False
    assert measurement["historical_1e-8_used_as_veto"] is False
    assert measurement["global_event_horizon_claimed"] is False
    assert measurement["black_hole_birth_claimed"] is False
    assert measurement["child_region_claimed"] is False
    assert measurement["regeneration_claimed"] is False
    assert measurement["sampling_is_continuous_horizon_proof"] is False
    assert measurement["cpu_budget_exceeded"] is False
    assert measurement["cpu_seconds"] <= CPU_BUDGET_S
    assert measurement["status"] == "MEASURED_SAMPLE_TRAPPING_CHARACTER_CHANGES"
    assert measurement["pattern"]["positive_geometry"] is True
    assert measurement["pattern"]["sample_character_changes"] is True
    assert measurement["pattern"]["arcs_track_areal_extrema"] is True
    assert measurement["pattern"]["strict_margin_agrees"] is True
    agreement = measurement["source_bindings"]["hash_agreement"]
    assert all(agreement.values())
    assert measurement["source_bindings"]["episode_files_unchanged"] is True
    assert measurement["metadata_normal_gap"] < 1e-12
    assert measurement["controls"]["inverse_metric_product_gap"] < 1e-12
    assert measurement["controls"]["null_norm_gap"] < 1e-12
    assert measurement["controls"]["boost_product_gap"] < 1e-12
    assert measurement["controls"]["analytic_sign_gap"] < 1e-12
    assert measurement["controls"]["positive_boost_preserves_each_sign"] is True
    assert measurement["controls"]["derivative_owner"]["field_gap"] < 1e-11
    comparison = measurement["comparisons"]
    assert comparison["dt_raw_counts_identical"] is True
    assert comparison["shared_node_sign_disagreements_max"] == 0
    assert comparison["dt_nf512_theta_plus_max_abs_gap"] < 1e-8
    assert measurement["clearance"]["clearance_ratio"] > 50.0
    by_name = {run["name"]: run for run in measurement["runs"]}
    fine = by_name["nf512_dt_0_0005"]
    coarse = by_name["nf256_dt_0_0005"]
    assert fine["frames"][0]["raw_counts"] == {
        "trapped": 0,
        "anti_trapped": 0,
        "untrapped": 2048,
        "marginal_indicator": 0,
    }
    assert fine["frames"][1]["raw_counts"]["trapped"] == 4
    assert fine["frames"][1]["raw_counts"]["anti_trapped"] == 8
    assert fine["frames"][-1]["raw_counts"]["trapped"] == 45
    assert fine["frames"][-1]["raw_counts"]["anti_trapped"] == 182
    assert fine["frames"][-1]["raw_counts"]["untrapped"] == 1821
    assert coarse["frames"][-1]["raw_counts"]["trapped"] == 23
    assert coarse["frames"][-1]["raw_counts"]["anti_trapped"] == 91
    assert fine["frames"][0]["product_max"] < 0.0
    assert fine["frames"][-1]["product_max"] > 0.0
    assert fine["frames"][0]["critical_point_estimates"]
    assert all(
        point["label"] == "marginal_within_margin"
        for point in fine["frames"][0]["critical_point_estimates"]
    )
    assert fine["frames"][-1]["product_crossing_count"] == 4
    assert fine["positive_r"] is True and fine["positive_Q"] is True
    assert fine["observer_proper_gap"] < 1e-14
    assert fine["high_mode_derivative_bound"] < 1e-8
    final_norm = fine["frames"][-1]["norm_gaps"]["g_nn_minus_one"]
    assert final_norm < 1e-12
    initial = fine["frames"][0]
    final = fine["frames"][-1]
    margin = measurement["margin"]["absolute_expansion_margin"]
    assert margin == pytest.approx(1.8325977887805045e-06)
    saved_fine = next(run for run in saved_record["runs"] if run["name"] == "nf512_dt_0_0005")
    peak = initial["areal_maximum"]
    saved_peak = saved_fine["frames"][0]["areal_maximum"]
    assert peak["index"] == saved_peak["index"] == 395
    assert abs(peak["r_x"]) > 1e-4
    recomputed_plus, recomputed_minus, recomputed_product, _spatial = expansions_from_stored_normal(
        peak["r"],
        peak["normal_velocity"],
        peak["r_x"],
        initial["Q_min"],
    )
    assert float(recomputed_plus) == pytest.approx(saved_peak["theta_plus"], abs=1e-14)
    assert float(recomputed_minus) == pytest.approx(saved_peak["theta_minus"], abs=1e-14)
    assert float(recomputed_plus) == pytest.approx(3.113501917378642e-4, rel=1e-9)
    assert float(recomputed_minus) == pytest.approx(-3.113501506382507e-4, rel=1e-9)
    assert float(recomputed_plus) * float(recomputed_minus) < 0.0
    assert float(recomputed_product) < 0.0
    critical = initial["critical_point_estimates"]
    theta_abs = max(abs(point["theta_both"]) for point in critical)
    assert theta_abs == pytest.approx(1.9636126344543498e-11, rel=1e-9)
    assert abs(theta_abs - 1.96e-11) < 5e-14
    assert theta_abs < margin
    assert abs(initial["normal_max"]) < 1e-9
    assert abs(final["normal_max"]) > 1e-2
    assert abs(final["areal_maximum"]["theta_plus"]) > margin
    assert abs(final["areal_maximum"]["theta_minus"]) > margin
    assert final["areal_maximum"]["theta_plus"] < 0.0 and final["areal_maximum"]["theta_minus"] < 0.0
    assert final["areal_minimum"]["theta_plus"] > margin
    assert final["areal_minimum"]["theta_minus"] > margin
    assert final["same_sign_min_abs_expansion"] > margin
    assert measurement["black_hole_birth_claimed"] is False
    assert measurement["global_event_horizon_claimed"] is False


def test_record_distinguishes_marginal_indicator_from_global_horizon(measurement, saved_record):
    assert saved_record["schema"] == measurement["schema"]
    assert saved_record["status"] == measurement["status"]
    assert saved_record["finding"] == measurement["finding"]
    assert saved_record["global_event_horizon_claimed"] is False
    assert saved_record["sampling_is_continuous_horizon_proof"] is False
    assert "global event horizon" in " ".join(saved_record["open"]).lower() or any(
        "event horizon" in item.lower() for item in saved_record["open"]
    )
    assert RECORD_PATH.stat().st_size < JSON_BYTE_LIMIT
    episode_hash = hashlib.sha256(EPISODE_NPZ.read_bytes()).hexdigest()
    assert episode_hash == saved_record["source_bindings"]["hashes"]["episode_npz"]
    assert episode_hash == saved_record["source_bindings"]["episode_npz_sha256_recorded_in_metadata"]
    json_hash = hashlib.sha256(EPISODE_JSON.read_bytes()).hexdigest()
    assert json_hash == saved_record["source_bindings"]["hashes"]["episode_json"]
    frozen_module = saved_record["source_bindings"]["module_sha256"]
    review = json.loads(_REVIEW.read_text())
    assert review["immutable_v1"]["null_expansion_recorded_module_sha256"] == frozen_module
    assert _source_matches(
        _MODULE, review["current_owner_sha256"]["nsc_spherical_null_expansion.py"],
        _V2_SOURCE_COMMIT)
    assert hashlib.sha256(_MODULE.read_bytes()).hexdigest() != frozen_module


def test_flipped_boost_return_fails_the_sign_flag(monkeypatch):
    import recursive_horizons.nsc_spherical_null_expansion as owner

    real = owner.boosted_expansions

    def flipped(theta_plus, theta_minus, scale):
        plus, minus, product = real(theta_plus, theta_minus, scale)
        return -plus, -minus, product

    monkeypatch.setattr(owner, "boosted_expansions", flipped)
    controls = owner.independent_controls()
    assert controls["positive_boost_preserves_each_sign"] is False
    assert controls["boost_product_gap"] < 1e-12


def test_observer_label_names_initial_mode_columns():
    from recursive_horizons.nsc_coupled_local_response import finite_domain

    label = finite_domain(
        {
            "case": "nf512_dt_0_0005",
            "frame_time": np.zeros(11),
            "interpolation": {"midpoint_linear_scale": 0.0},
            "time_resolution_geometry_gap": 0.0,
        },
        {"t_start": 0.0, "t_end": 0.05, "steps": 10},
        {"nf": 512, "nq": 2048, "length": 8.0},
        0.0,
        True,
    )["observer"]
    assert label == "initial mode columns 0 and 1, weights 0.75 and 0.75, no phase QR"
    assert "plus and minus" not in label


def _file_sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _owner_paths():
    return {
        "nsc_spherical_null_expansion.py": _MODULE,
        "test_nsc_spherical_null_expansion.py": Path(__file__),
        "nsc-spherical-null-expansion.md": _NULL_DOC,
        "nsc_coupled_local_response.py": _RESPONSE_MODULE,
        "nsc-coupled-local-response.md": _RESPONSE_DOC,
    }


def _payload_paths():
    return {
        "episode_json_sha256": EPISODE_JSON,
        "episode_npz_sha256": EPISODE_NPZ,
        "null_expansion_json_sha256": RECORD_PATH,
        "coupled_response_json_sha256": _RESPONSE_JSON,
        "coupled_response_npz_sha256": _RESPONSE_NPZ,
    }


def _show_pinned_source(commit, path):
    env = os.environ.copy()
    env["GIT_ALTERNATE_OBJECT_DIRECTORIES"] = str(_SOURCE_HISTORY / "objects")
    return subprocess.check_output(
        ["git", "-C", str(_SOURCE_HISTORY.parent), "show", f"{commit}:{path}"],
        env=env,
        stderr=subprocess.PIPE,
    )


def _source_matches(path, expected, historical_commit):
    if historical_commit is None:
        return path.is_file() and _file_sha256(path) == expected
    root = _SOURCE_HISTORY.parent.parent
    try:
        resolve_pinned_source_bytes(
            root, path.resolve().relative_to(root).as_posix(), expected,
            commit=historical_commit)
    except (ValueError, RuntimeError):
        return False
    return True


def _successor_binding_errors(review, *, historical_commit=None):
    """Strict current binding, or explicit replay of a sealed review's sources."""
    errors = []
    immutable = review.get("immutable_v1")
    if not isinstance(immutable, dict):
        return ["immutable_v1"]
    for key, path in _payload_paths().items():
        if key not in immutable or not path.is_file() or _file_sha256(path) != immutable[key]:
            errors.append(key)
    if immutable.get("null_expansion_recorded_module_sha256") != _OLD_MODULE_SHA256:
        errors.append("null_expansion_recorded_module_sha256")
    if immutable.get("stored_observer_label") != _STORED_OBSERVER:
        errors.append("stored_observer_label")
    owners = review.get("current_owner_sha256")
    if not isinstance(owners, dict):
        errors.append("current_owner_sha256")
        owners = {}
    for name, path in _owner_paths().items():
        if name not in owners or not _source_matches(path, owners[name], historical_commit):
            errors.append(name)
    independent = _SOURCE_HISTORY.parent / "tests" / "test_nsc_local_boundary_independent.py"
    current = review.get("reference_current")
    if not isinstance(current, dict) or not _source_matches(
            independent, current.get("independent_test_sha256"), historical_commit):
        errors.append("reference_current.independent_test_sha256")
    sealed = review.get("reference_at_v1_seal")
    if not isinstance(sealed, dict) or sealed.get("independent_test_sha256") != _SEALED_INDEPENDENT_SHA256:
        errors.append("reference_at_v1_seal")
    historical = review.get("historical_sources")
    if not isinstance(historical, dict):
        errors.append("historical_sources")
        return errors
    if historical.get("recovered") is not True or historical.get("archived") is not True:
        errors.append("historical_recovery")
    if historical.get("carrier_is_historical_commit") is not False:
        errors.append("historical_carrier")
    if historical.get("null_expansion_module_sha256") != _OLD_MODULE_SHA256:
        errors.append("historical_module")
    if historical.get("independent_test_sha256") != _SEALED_INDEPENDENT_SHA256:
        errors.append("historical_test")
    pins = historical.get("pins")
    if not isinstance(pins, list):
        errors.append("historical_pins")
        return errors
    seen = {}
    for pin in pins:
        if not isinstance(pin, dict):
            errors.append("historical_pin")
            continue
        path = pin.get("path")
        digest = pin.get("sha256")
        commit = pin.get("commit")
        blob = pin.get("blob")
        if path not in _HISTORICAL_PIN_PATHS or digest != _HISTORICAL_PIN_PATHS[path] or not commit or not blob:
            errors.append("historical_pin:" + str(path))
            continue
        try:
            raw = _show_pinned_source(commit, path)
            framed = zlib.decompress((_SOURCE_HISTORY / "objects" / blob[:2] / blob[2:]).read_bytes())
            header, blob_raw = framed.split(b"\0", 1)
        except (OSError, subprocess.CalledProcessError, zlib.error, ValueError):
            errors.append("historical_replay:" + path)
            continue
        if hashlib.sha1(framed).hexdigest() != blob or not header.startswith(b"blob "):
            errors.append("historical_replay:" + path)
        elif hashlib.sha256(raw).hexdigest() != digest or hashlib.sha256(blob_raw).hexdigest() != digest:
            errors.append("historical_replay:" + path)
        else:
            seen[path] = digest
    if set(seen) != set(_HISTORICAL_PIN_PATHS):
        errors.append("historical_pin_paths")
    return errors


def test_successor_review_binds_immutable_v1(measurement, saved_record):
    v1_bytes = _REVIEW.read_bytes()
    assert _file_sha256(_REVIEW) == _V1_REVIEW_SHA256
    v1 = json.loads(v1_bytes)
    v2 = json.loads(_REVIEW_V2.read_text())
    assert v2["schema"] == "NSC-LOCAL-BOUNDARY-REVIEW-SUCCESSOR-v2"
    assert v2["predecessor"]["path"] == "lab/results/development/nsc-local-boundary-review-v1.json"
    assert v2["predecessor"]["sha256"] == _V1_REVIEW_SHA256
    assert v2["predecessor"]["immutable"] is True
    assert v2["predecessor"]["rebound"] is False
    assert v2["immutable_v1"] == v1["immutable_v1"]
    assert v2["controls"] == v1["controls"]
    assert v2["domain"]["later_independent_tests"] == list(_LATER_INDEPENDENT_TESTS)
    independent = _SOURCE_HISTORY.parent / "tests" / "test_nsc_local_boundary_independent.py"
    independent_text = independent.read_text()
    for name in _LATER_INDEPENDENT_TESTS:
        assert f"def {name}(" in independent_text
    assert _successor_binding_errors(v2, historical_commit=_V2_SOURCE_COMMIT) == []
    sealed_pin = next(
        pin for pin in v2["historical_sources"]["pins"]
        if pin["path"] == "tests/test_nsc_local_boundary_independent.py"
    )
    sealed_text = _show_pinned_source(sealed_pin["commit"], sealed_pin["path"]).decode()
    for name in _LATER_INDEPENDENT_TESTS:
        assert f"def {name}(" not in sealed_text
    owners = v2["current_owner_sha256"]
    assert owners["nsc_spherical_null_expansion.py"] == v1["current_owner_sha256"]["nsc_spherical_null_expansion.py"]
    assert owners["nsc_coupled_local_response.py"] == v1["current_owner_sha256"]["nsc_coupled_local_response.py"]
    assert owners["test_nsc_spherical_null_expansion.py"] != v1["current_owner_sha256"]["test_nsc_spherical_null_expansion.py"]
    assert owners["nsc-spherical-null-expansion.md"] != v1["current_owner_sha256"]["nsc-spherical-null-expansion.md"]
    assert owners["nsc-coupled-local-response.md"] != v1["current_owner_sha256"]["nsc-coupled-local-response.md"]
    assert v1["reference_unchanged"]["independent_test_sha256"] == _SEALED_INDEPENDENT_SHA256
    assert v2["reference_current"]["independent_test_sha256"] != _SEALED_INDEPENDENT_SHA256
    immutable = v2["immutable_v1"]
    assert immutable["null_expansion_recorded_module_sha256"] == saved_record["source_bindings"]["module_sha256"]
    stored = json.loads(_RESPONSE_JSON.read_text())
    assert stored["finite_approximation_domain"]["observer"] == immutable["stored_observer_label"]
    saved_fine = next(run for run in saved_record["runs"] if run["name"] == "nf512_dt_0_0005")
    peak = saved_fine["frames"][0]["areal_maximum"]
    assert v2["controls"]["positive_boost_preserves_each_sign"] is True
    assert v2["controls"]["wrong_return_preserves_each_sign"] is False
    assert v2["controls"]["boost_product_gap"] == measurement["controls"]["boost_product_gap"]
    assert v2["controls"]["t0_areal_maximum_theta_plus"] == peak["theta_plus"]
    assert v2["controls"]["t0_areal_maximum_theta_minus"] == peak["theta_minus"]
    assert v2["controls"]["absolute_expansion_margin"] == saved_record["margin"]["absolute_expansion_margin"]
    manifest = json.loads((_SOURCE_HISTORY / "manifest.json").read_text())
    assert manifest["schema"] == "NSC-PINNED-SOURCE-OBJECTS-v1"
    rows = manifest["objects"]
    assert len(_ORIGINAL_SOURCE_OBJECTS) == 18
    assert [(row["oid"], row["type"], row["bytes"], row["sha256"]) for row in rows[:18]] == list(_ORIGINAL_SOURCE_OBJECTS)
    assert [(row["commit"], row["path"], row["blob"], row["sha256"]) for row in manifest["pins"][:4]] == list(_ORIGINAL_PINS)
    for row in rows[:18]:
        oid = row["oid"]
        framed = zlib.decompress((_SOURCE_HISTORY / "objects" / oid[:2] / oid[2:]).read_bytes())
        assert hashlib.sha1(framed).hexdigest() == oid
        assert hashlib.sha256(framed).hexdigest() == row["sha256"]
    object_ids = {row["oid"] for row in rows}
    for pin in v2["historical_sources"]["pins"]:
        assert pin["commit"] in object_ids
        assert pin["blob"] in object_ids
        match = next(row for row in manifest["pins"] if row["path"] == pin["path"] and row["commit"] == pin["commit"])
        assert match["blob"] == pin["blob"]
        assert match["sha256"] == pin["sha256"]


def test_sealed_successor_review_rejects_changed_or_missing_dependency():
    review = json.loads(_REVIEW_V2.read_text())
    assert _successor_binding_errors(review, historical_commit=_V2_SOURCE_COMMIT) == []
    changed = json.loads(_REVIEW_V2.read_text())
    changed["current_owner_sha256"]["nsc_coupled_local_response.py"] = "0" * 64
    assert "nsc_coupled_local_response.py" in _successor_binding_errors(changed, historical_commit=_V2_SOURCE_COMMIT)
    missing = json.loads(_REVIEW_V2.read_text())
    del missing["current_owner_sha256"]["test_nsc_spherical_null_expansion.py"]
    assert "test_nsc_spherical_null_expansion.py" in _successor_binding_errors(missing, historical_commit=_V2_SOURCE_COMMIT)
    dropped_payload = json.loads(_REVIEW_V2.read_text())
    del dropped_payload["immutable_v1"]["episode_npz_sha256"]
    assert "episode_npz_sha256" in _successor_binding_errors(dropped_payload, historical_commit=_V2_SOURCE_COMMIT)
    bad_pin = json.loads(_REVIEW_V2.read_text())
    bad_pin["historical_sources"]["pins"][0]["sha256"] = "1" * 64
    assert any(item.startswith("historical_pin") for item in _successor_binding_errors(bad_pin, historical_commit=_V2_SOURCE_COMMIT))
    missing_object = json.loads(_REVIEW_V2.read_text())
    missing_object["historical_sources"]["pins"][0]["commit"] = "0" * 40
    assert any(item.startswith("historical_replay") for item in _successor_binding_errors(missing_object, historical_commit=_V2_SOURCE_COMMIT))


def test_current_dependency_binding_does_not_fall_back_to_sealed_sources():
    review = json.loads(_REVIEW_V2.read_text())
    assert "test_nsc_spherical_null_expansion.py" in _successor_binding_errors(review)
    review["current_owner_sha256"] = {
        name: _file_sha256(path) for name, path in _owner_paths().items()
    }
    independent = _SOURCE_HISTORY.parent / "tests/test_nsc_local_boundary_independent.py"
    review["reference_current"]["independent_test_sha256"] = _file_sha256(independent)
    assert _successor_binding_errors(review) == []
    review["current_owner_sha256"]["nsc_coupled_local_response.py"] = "0" * 64
    assert "nsc_coupled_local_response.py" in _successor_binding_errors(review)


def test_saved_sample_formula_binding_names_its_historical_domain(measurement):
    assert measurement["source_bindings"]["source_ref"] == _V2_SOURCE_COMMIT
    assert measurement["source_bindings"]["hash_agreement"]["galerkin_matches_episode_hashes_after"] is True
    current = build_record()
    assert current["status"] == "SOURCE_BINDING_MISMATCH"
    assert current["source_bindings"]["source_ref"] is None
    assert current["source_bindings"]["hash_agreement"]["galerkin_matches_episode_hashes_after"] is False
    missing = build_record(source_ref="0" * 40)
    assert missing["status"] == "SOURCE_BINDING_MISMATCH"
    assert missing["source_bindings"]["historical_source_hashes"] == {}
