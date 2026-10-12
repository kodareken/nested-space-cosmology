"""Independent finite-Fock/counting controls and saved-source geometric vertices."""
from copy import deepcopy
import json
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.linalg import expm

from recursive_horizons import nsc_discovery_geometric_vertices as vertices
from recursive_horizons import nsc_influence as influence


def test_gaussian_bilinear_third_cumulant_matches_independent_small_fock_oracle():
    covariance = np.array([[.7, .08 + .04j], [.08 - .04j, .3]])
    g = np.array([[.8, .3 - .2j], [.3 + .2j, -.4]])
    annihilators = influence.fock_annihilators(2)
    density = influence.gaussian_fock_state(covariance, annihilators)
    operator = influence.second_quantize(g, annihilators)
    mean = np.trace(density @ operator)
    centered = operator - mean * np.eye(4)
    variance = np.trace(density @ centered @ centered)
    third = np.trace(density @ centered @ centered @ centered)
    calculated = vertices.bilinear_cumulants(covariance, g)
    assert calculated["kappa1"] == pytest.approx(mean.real, abs=1e-14)
    assert calculated["kappa2"] == pytest.approx(variance.real, abs=1e-14)
    assert calculated["kappa3"] == pytest.approx(third.real, abs=1e-14)
    assert abs(third.real) > .01
    assert calculated["elementary_fermion_state_is_Gaussian"]
    assert calculated["bilinear_force_is_not_assumed_Gaussian"] is True


def test_thin_counting_determinant_keeps_full_vertex_images_and_matches_fock():
    c = np.array([.75, .35])
    columns = np.eye(3, dtype=complex)[:, :2]
    covariance = (columns * c) @ columns.conj().T
    g = np.array([[.4, .2j, .3], [-.2j, -.6, .5j], [.3, -.5j, .9]])
    lam = .13
    thin = vertices.counting_logdet(columns, c, g, lam)
    full = np.linalg.det(np.eye(3) - covariance + covariance @ expm(-1j * lam * g))
    assert np.exp(thin) == pytest.approx(full, abs=2e-14)
    # Mutation: merely compressing G into occupied columns drops real ambient paths.
    truncated = np.linalg.det(np.eye(2) - np.diag(c) + np.diag(c) @ expm(-1j * lam * g[:2, :2]))
    assert abs(full - truncated) > 1e-4


def test_occupied_empty_cut_matches_independent_fock_commutator():
    covariance = np.array([[.65, .08j], [-.08j, .25]])
    h = np.array([[.4, .3 - .1j], [.3 + .1j, -.7]])
    child = np.array([[.1, .4j], [-.4j, .7]])
    parent = np.array([[.8, .2], [.2, -.3]])
    result = vertices.occupied_empty_cut(covariance, h, child, parent)
    a = 1j * (h @ child - child @ h)
    operators = influence.fock_annihilators(2)
    rho = influence.gaussian_fock_state(covariance, operators)
    af = influence.second_quantize(a, operators)
    bf = influence.second_quantize(parent, operators)
    exact = np.trace(rho @ (af @ bf - bf @ af))
    expected = complex(result["commutator"]["real"], result["commutator"]["imag"])
    assert expected == pytest.approx(exact, abs=2e-14)
    assert result["cut_identity_gap"] < 1e-14
    assert result["instantaneous_retarded_cross_slope"] == pytest.approx((-1j * exact).real, abs=1e-14)
    assert not result["stationarity_assumed"]
    assert not result["physical_pole_claimed"]


@pytest.fixture(scope="module")
def actual():
    return vertices.calculate()


def test_actual_saved_state_has_nonzero_third_vertex_and_common_action_mean(actual):
    assert actual["nf"] == 128 and actual["coordinate_time"] == .3
    assert actual["case"] == "nf128_baseline_dt0.0005"
    assert actual["CAR_admissible"] and actual["full_ambient_band_retained"]
    assert not actual["evolution_ran"] and not actual["initial_state_called"]
    assert actual["covariance_stationarity_gap"] > 1e-4
    for index in ("0", "5", "6"):
        row = actual["vertices"][index]
        raw = row["representative_channel"]
        assert abs(raw["kappa3"]) > 1e-7
        assert raw["kappa3_imaginary_residual"] < 1e-15
        assert row["mean_force_mapping_gap"] < 1e-12
        assert row["physical_operator"]["operator_dimension"] == 256
        assert not row["physical_operator"]["representative_vertex_has_M"]
        counting = row["counting_field_check"]
        assert counting["held_out_not_a_derivative_step"]
        assert counting["held_out_gap_over_cubic_effect"] < .01
        for estimate in counting["derivative_steps"]:
            assert estimate["gap_from_trace_kappa3"] < .01 * abs(raw["kappa3"])
    assert actual["vertices"]["6"]["representative_channel"]["kappa3"] == pytest.approx(1.1498940792281297e-5, rel=1e-8)
    assert actual["occupied_empty_cut"]["cut_identity_gap"] < 1e-13
    assert actual["occupied_empty_cut"]["instantaneous_retarded_cross_slope"] == pytest.approx(-.015091406833530857, rel=1e-8)
    assert actual["cpu_seconds"] < actual["cpu_limit_seconds"] <= 10


def test_total_angular_noise_is_not_invented_from_the_mean_mapping(actual):
    scope = actual["multiplicity_scope"]
    assert scope["M"] == 4
    assert scope["angular_product_state_derived_here"] is False
    assert scope["total_stress_noise_mapping"] is None
    assert "M^3" in scope["scaled_observable_assumption"]
    assert "M kappa3" in scope["independent_angular_copy_assumption"]
    assert not actual["new_dynamics"] and not actual["cosmohedron_equivalence_claimed"]


def test_frozen_git_authentication_rejects_missing_or_different_bytes(monkeypatch):
    commit = "a" * 40
    blob = b"finite physical source\n"
    expected = __import__("hashlib").sha256(blob).hexdigest()
    oid = "b" * 40
    answer = oid.encode() + b" blob " + str(len(blob)).encode() + b"\n" + blob + b"\n"
    monkeypatch.setattr(vertices.subprocess, "run", lambda *args, **kwargs: SimpleNamespace(stdout=answer))
    authenticated = vertices.authenticate_git(commit, {"lab/owner.py": expected})
    assert authenticated["all_declared_bytes_authenticated"]
    assert authenticated["blob_ids"]["lab/owner.py"] == oid
    with pytest.raises(ValueError, match="differ"):
        vertices.authenticate_git(commit, {"lab/owner.py": "0" * 64})
    monkeypatch.setattr(vertices.subprocess, "run", lambda *args, **kwargs: SimpleNamespace(stdout=b"request missing\n"))
    with pytest.raises(ValueError, match="absent"):
        vertices.authenticate_git(commit, {"lab/owner.py": expected})


def test_exclusive_temp_write_and_readonly_replay_reject_numeric_mutation(actual, tmp_path, monkeypatch):
    bindings = actual["source_bindings_before"]
    frozen = {"commit": "a" * 40, "all_declared_bytes_authenticated": True, "blob_ids": {}}
    monkeypatch.setattr(vertices, "calculate", lambda: deepcopy(actual))
    monkeypatch.setattr(vertices, "source_bindings", lambda: bindings.copy())
    monkeypatch.setattr(vertices, "authenticate_git", lambda *args: frozen.copy())
    path = tmp_path / "vertices.json"
    vertices.write_record(path, "a" * 40)
    before = path.read_bytes()
    checked = vertices.check_record(path)
    assert checked["vertices"]["6"]["representative_channel"]["kappa3"] != 0
    assert path.read_bytes() == before
    with pytest.raises(FileExistsError):
        vertices.write_record(path, "a" * 40)
    changed = json.loads(path.read_text())
    changed["vertices"]["6"]["representative_channel"]["kappa3"] += 1e-5
    path.write_text(json.dumps(changed))
    with pytest.raises(ValueError, match="number changed"):
        vertices.check_record(path)
