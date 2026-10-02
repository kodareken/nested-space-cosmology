"""Read-only historical authentication of the sealed feedback episode.

The current producer pin is not reapplied to the saved JSON or NPZ.
Historical source is parsed and is not imported. The episode generator is
not run.
"""
import hashlib
import os
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import pytest

import check_nsc_spherical_feedback_episode_historical as historical

ROOT = Path(__file__).resolve().parents[2]
EPISODE_JSON = ROOT / "lab/results/development/nsc-spherical-feedback-episode-v1.json"
EPISODE_NPZ = ROOT / "lab/results/development/nsc-spherical-feedback-episode-v1.npz"
CAUCHY_JSON = ROOT / "lab/results/development/nsc-spherical-cauchy-data-v1.json"
GALERKIN = ROOT / "lab/src/recursive_horizons/nsc_spherical_galerkin_coupling.py"
EPISODE_JSON_SHA256 = "9fc9f460aa076db5eb67a3136d3e53797160f2885e124b62e6fb61ee247d9954"
EPISODE_NPZ_SHA256 = "2864d3a8d3ba413e840c395961f14a3d71eabee702c27fd9f47c7d333eb24b77"
CAUCHY_JSON_SHA256 = "2c6a9ede084e3a56ae7798e44ac017d1086e591f4a2165d26e75a71ff98150f7"
UNPINNED = {
    "nsc_adm_source",
    "nsc_covariant_operator",
    "nsc_influence",
    "nsc_nested_qualities",
}


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_historical_consumer_authenticates_context_and_saved_arrays():
    assert _sha256(EPISODE_JSON) == EPISODE_JSON_SHA256
    assert _sha256(EPISODE_NPZ) == EPISODE_NPZ_SHA256
    assert _sha256(CAUCHY_JSON) == CAUCHY_JSON_SHA256
    assert _sha256(GALERKIN) == historical.CURRENT_GALERKIN_SHA256
    report = historical.authenticate()
    assert report["historical_commit"] == historical.HISTORICAL_COMMIT
    assert report["historical_galerkin_sha256"] == historical.HISTORICAL_GALERKIN_SHA256
    assert report["historical_galerkin_bytes"] == 52028
    assert report["source_origin"] == "resolve_pinned_source_bytes"
    assert report["current_galerkin_sha256"] == historical.CURRENT_GALERKIN_SHA256
    assert report["current_code_compatible"] is False
    assert report["current_code_replay_claimed"] is False
    assert report["current_direct_modules"] == [
        "nsc_conformal_adm_source",
        "nsc_spherical_coupling",
        "nsc_spherical_feedback_action",
    ]
    inspection = report["dependency_inspection"]
    assert inspection["inspected_before_import"] is True
    assert inspection["historical_direct_modules"] == [
        "nsc_conformal_adm_source",
        "nsc_spherical_coupling",
    ]
    assert inspection["missing_names"] == []
    assert inspection["module_level_calls"] == []
    assert inspection["import_allowed"] is False
    assert set(inspection["unpinned_transitive_modules"]) == UNPINNED
    assert "nsc_spherical_feedback_action" in inspection["authenticated_import_modules"]
    assert report["historical_source_imported"] is False
    assert report["producer_invoked"] is False
    assert report["generation_invoked"] is False
    assert inspection["producer_defines_run_episode"] is True
    assert report["saved_arrays_authenticated"] is True
    assert report["driver_before_matches_commit"] is False
    assert report["driver_after_sha256"] != report["driver_before_sha256"]
    by_name = {item["name"]: item for item in report["dependencies"]}
    assert set(by_name) == {name for name, _path, _role in historical.DEPENDENCIES}
    galerkin = by_name["galerkin"]
    assert galerkin["origin"] == "resolve_pinned_source_bytes"
    assert galerkin["commit"] == historical.HISTORICAL_COMMIT
    assert galerkin["authenticated_sha256"] == historical.HISTORICAL_GALERKIN_SHA256
    assert galerkin["working_tree_matches_historical"] is False
    for name, item in by_name.items():
        if name == "galerkin":
            continue
        if name == "driver":
            assert item["working_tree_sha256"] == item["declared_after"]
            assert item["declared_before"] != item["declared_after"]
            assert item["origin"] == "resolve_pinned_source_bytes"
            continue
        assert item["declared_before"] == item["declared_after"]
        assert item["working_tree_matches_historical"] is True
        if item["role"] == "scientific_source":
            assert item["origin"] == "resolve_pinned_source_bytes"
    assert report["sealed_bytes_unchanged"] is True
    assert report["episode_json_sha256"] == EPISODE_JSON_SHA256
    assert report["episode_npz_sha256"] == EPISODE_NPZ_SHA256
    assert _sha256(EPISODE_JSON) == EPISODE_JSON_SHA256
    assert _sha256(EPISODE_NPZ) == EPISODE_NPZ_SHA256
    assert _sha256(CAUCHY_JSON) == CAUCHY_JSON_SHA256
    assert _sha256(GALERKIN) == historical.CURRENT_GALERKIN_SHA256


def test_pinned_source_failure_does_not_heal_or_generate(monkeypatch):
    before = {
        EPISODE_JSON: EPISODE_JSON.read_bytes(),
        EPISODE_NPZ: EPISODE_NPZ.read_bytes(),
        CAUCHY_JSON: CAUCHY_JSON.read_bytes(),
        GALERKIN: GALERKIN.read_bytes(),
    }

    def unavailable(*_args, **_kwargs):
        raise RuntimeError("pinned source drift: nsc_spherical_galerkin_coupling.py")

    monkeypatch.setattr(historical, "resolve_pinned_source_bytes", unavailable)
    with pytest.raises(RuntimeError, match="pinned source drift"):
        historical.authenticate()
    for path, raw in before.items():
        assert path.read_bytes() == raw
