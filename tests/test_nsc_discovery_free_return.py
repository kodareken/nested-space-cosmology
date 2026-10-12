"""Exact AP identities, saved-state decoding and bounded saved-field comparison."""
import copy
import json

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

from recursive_horizons import nsc_discovery_free_return as free


def test_exact_AP_period_is_minus_identity_and_every_window_probability_returns():
    with threadpool_limits(limits=1):
        grid = free.galerkin.build_grid(32, gauge="conformal")
        backend = free.free_backend(grid)
        assert backend["AP_spectrum_gap"] < 1e-10
        eye = np.eye(2*grid.nf, dtype=complex)
        assert np.max(abs(free.propagate(backend, eye, 8)+eye)) < 1e-10
        rng = np.random.default_rng(37)
        phi = rng.normal(size=(2*grid.nf, 6))+1j*rng.normal(size=(2*grid.nf, 6))
        weights = np.array([.75, .75, .5, .5, .25, .25])
        out = free.propagate(backend, phi, 8)
        for interval in ((1, 3), (0, 4), (0, 8)):
            assert free.fixed_window_probability(grid, out, weights, interval) == pytest.approx(
                free.fixed_window_probability(grid, phi, weights, interval), abs=1e-10)
        assert np.max(abs(out.conj().T@out-phi.conj().T@phi)) < 1e-10
        assert np.max(abs(free.propagate(backend, np.zeros_like(phi), 5))) == 0


@pytest.fixture(scope="module")
def report():
    return free.measure()


def test_saved_join_decoding_operator_and_CAR_controls(report):
    assert report["cpu_seconds"] < 10
    assert report["input_sources_unchanged"]
    assert all(not free.Path(name).is_absolute() for name in report["source_hashes"] | report["input_hashes"])
    for row in report["cases"]:
        assert row["join"]["exact_T3_handoff"]
        assert min(row["Q_sup_samples"].values()) > 0
        assert "reconstruct_state before" in row["momentum_geometry_decoding"]
        assert row["source_cross_terms_omitted"] is False
        assert row["controls"]["actual_H_minus_free_minus_mass_max_gap"] < 1e-9
        assert row["controls"]["zero_mass_operator_max_gap"] < 1e-9
        assert row["controls"]["free_Gram_change_max"] < 1e-10
        assert row["fixed_child_window"]["saved_ledger_gap"] < 1e-10
        assert row["Duhamel"]["certified_error_upper_bound"] is None
        assert row["Duhamel"]["endpoint_samples_are_integral_bound"] is False
        eigen = row["gram_CAR"]["free"]["occupied_covariance_eigenvalues"]
        assert min(eigen) >= 0 and max(eigen) < 1
    old, _ = free.episode.load_checkpoint(free.PREDECESSOR, free.PRIMARY, 2)
    handoff, _ = free.episode.load_checkpoint(free.CROSSING, free.PRIMARY, 0)
    endpoint, _ = free.episode.load_checkpoint(free.CROSSING, free.PRIMARY, 1)
    bad = copy.deepcopy(handoff)
    bad["array_sha256"].pop("Q")
    with pytest.raises(ValueError, match="continuation array dictionaries"):
        free.validate_join(old, bad, endpoint)


def test_nearly_free_primary_relation_is_distinct_from_frozen_return(report):
    primary = report["cases"][0]
    assert primary["weighted_field_relative_error"] == pytest.approx(2.982e-6, rel=0.005)
    assert primary["fixed_child_window"]["absolute_free_saved_gap"] < 1e-10
    assert primary["source_covariance_relative_error"] < 1e-5
    frozen = next(row for row in report["cases"] if row["case_id"] == "nf256_frozen_geometry_dt0.0005")
    assert frozen["weighted_field_relative_error"] > 1000*primary["weighted_field_relative_error"]
    assert report["interpretation"]["renewal_demonstrated"] is False
    assert report["interpretation"]["frozen_return_alone_proves_free_limit"] is False
    assert free.duhamel_bound(1, np.sqrt(3), .001) == pytest.approx(np.sqrt(3)*.001)
    with pytest.raises(ValueError):
        free.duhamel_bound(1, 1, -1)


def test_creation_only_and_readonly_exact_replay(report, tmp_path, monkeypatch):
    destination = tmp_path / "free-return.json"
    monkeypatch.setattr(free, "OUTPUT", destination)
    free.write_record(report)
    before = destination.read_bytes()
    assert free.check_record()["recomputed_exact_free_operator"]
    assert destination.read_bytes() == before
    with pytest.raises(FileExistsError):
        free.write_record(report)
    changed = json.loads(before)
    changed["input_hashes"][next(iter(changed["input_hashes"]))] = "0"*64
    destination.write_text(json.dumps(changed))
    with pytest.raises(ValueError, match="bound hash changed"):
        free.check_record()
