"""Independent Fock oracles and the finite Dirac-contact benchmark."""
import hashlib
import io
import json

import numpy as np
import pytest
import sympy as sp
from scipy.linalg import block_diag, expm

from recursive_horizons import nsc_discovery_dirac_contact as contact
from recursive_horizons import nsc_influence as influence
from recursive_horizons.nsc_conformal_adm_source import direct_hamiltonian
from recursive_horizons.nsc_covariant_operator import smooth_metric
from recursive_horizons.nsc_curvature_eft import spherical_contact, stress_contact


def test_one_mode_fock_splits_the_product_of_means_from_the_connected_piece():
    occupation, left, right = contact.ONE_MODE_OCCUPATION, contact.ONE_MODE_LEFT, contact.ONE_MODE_RIGHT
    density = np.diag([1 - occupation, occupation]).astype(complex)
    number = np.diag([0.0, 1.0]).astype(complex)
    full = np.trace(density @ (left * number) @ (right * number))
    product = np.trace(density @ (left * number)) * np.trace(density @ (right * number))
    connected = occupation * (1 - occupation) * left * right
    result = contact.one_mode_oracle()
    assert full == pytest.approx(product + connected, abs=1e-14)
    assert result["connected"]["real"] == pytest.approx(connected, abs=1e-14)
    assert result["product_of_means"]["real"] == pytest.approx(product, abs=1e-14)
    assert result["trace_gap"] < 1e-12
    assert result["means_entered_the_connected_piece"] is False


def test_two_mode_fock_matches_the_wick_trace_and_not_the_mean_product():
    oracle = contact.two_mode_wick_oracle()
    operators = influence.fock_annihilators(2)
    density = influence.gaussian_fock_state(oracle["covariance"], operators)
    left = influence.second_quantize(oracle["left"], operators)
    right = influence.second_quantize(oracle["right"], operators)
    full = np.trace(density @ left @ right)
    product = (np.trace(density @ left) * np.trace(density @ right))
    wick = influence.connected(oracle["covariance"], oracle["left"], oracle["right"])
    assert full == pytest.approx(product + wick, abs=1e-12)
    assert abs(wick) > 1e-4
    assert abs(full - product) > 1e-4
    assert oracle["gap"] < 1e-12


def test_four_copies_scale_noise_by_four_and_a_scaled_observable_by_sixteen():
    covariance = np.array([[contact.ONE_MODE_OCCUPATION]], dtype=complex)
    vertex = np.array([[contact.ONE_MODE_LEFT]], dtype=complex)
    one = influence.connected(covariance, vertex, vertex)
    summed = influence.connected(
        block_diag(*([covariance] * 4)), block_diag(*([vertex] * 4)), block_diag(*([vertex] * 4)))
    scaled = influence.connected(covariance, 4 * vertex, 4 * vertex)
    result = contact.multiplicity_noise()
    assert summed == pytest.approx(4 * one, abs=1e-12)
    assert scaled == pytest.approx(16 * one, abs=1e-12)
    assert result["product_sum_over_one"] == pytest.approx(4)
    assert result["scaled_over_one"] == pytest.approx(16)
    assert result["noise_set_to_M_squared_without_a_declared_state"] is False


def test_three_mode_quartic_is_not_a_gaussian_preserving_mean_force():
    covariance, vertex = contact._quartic_inputs()
    operators = influence.fock_annihilators(3)
    density = influence.gaussian_fock_state(covariance, operators)
    bilinear = influence.second_quantize(vertex, operators)
    quartic = bilinear @ bilinear
    evolved = expm(-1j * contact.QUARTIC_EPSILON * quartic) @ density @ expm(1j * contact.QUARTIC_EPSILON * quartic)
    matrix = np.zeros((3, 3), dtype=complex)
    for i in range(3):
        for j in range(3):
            matrix[i, j] = np.trace(evolved @ operators[j].conj().T @ operators[i])
    matrix = 0.5 * (matrix + matrix.conj().T)
    rebuilt = influence.gaussian_fock_state(matrix, operators)
    result = contact.quartic_witness()
    assert np.linalg.norm(evolved - rebuilt) == pytest.approx(result["exact_state_minus_gaussian_rebuild"], abs=1e-12)
    assert result["exact_state_minus_gaussian_rebuild"] > 1e-4
    assert result["scalar_expectation_moves_state"] == pytest.approx(0, abs=1e-12)
    assert result["hartree_covariance_minus_exact"] > 1e-4
    assert result["four_fermion_installed_as_gaussian_mean_force"] is False
    assert result["P_minus_used"] is False


def test_pure_radial_current_is_outside_the_collective_specialization():
    eta = sp.diag(1, -1, -1, -1)
    current = sp.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    pure = stress_contact(current, eta, 1, 1, -2)
    collective = spherical_contact(2, 0, 0, 1, 0, 1, 0, 1, -2)
    witness = contact.specialization_witness()
    assert sp.simplify(collective) == 0
    assert sp.simplify(pure) != 0
    assert witness["scalar_collective_specialization_used_as_dirac_stress"] is False
    assert witness["difference"] != "0"


def test_projection_solves_the_owned_inverse_metric_variation():
    identity = contact.projection_identity()
    assert identity["gaps"] == {"N": "0", "beta": "0", "q": "0", "r": "0"}
    projected = contact.orthonormal_from_force_densities(
        1.0, 2.0, 3.0, 4.0, 1.5, 0.8, 1.2)
    rho, current, radial, angular = (
        projected["rho"], projected["current"], projected["p_radial"], projected["p_angular"])
    assert rho == pytest.approx(1 / (4 * np.pi * 0.8 * 1.2**2))
    assert current == pytest.approx(4 / (4 * np.pi * 0.8**2 * 1.2**2))
    assert radial == pytest.approx(-2 / (4 * np.pi * 1.5 * 1.2**2))
    assert angular == pytest.approx(-3 / (8 * np.pi * 1.5 * 0.8 * 1.2))


@pytest.fixture(scope="module")
def benchmark():
    record, arrays, payload = contact.assemble()
    return record, arrays, payload


def test_benchmark_variations_wick_completion_and_open_scope(benchmark):
    record, arrays, payload = benchmark
    assert record["schema"] == contact.SCHEMA
    assert record["scope_open"] is True
    assert record["blocks_parent_construction"] is False
    assert record["missing"]["four_dimensional_angular_vertex"] is None
    assert record["missing"]["angular_current_map"] is None
    assert record["missing"]["renormalized_contact_expectation"] is None
    assert record["off_shell_identity"]["gap"] == "0"
    assert record["off_shell_identity"]["boundary_variation_discarded"] is False
    assert record["off_shell_identity"]["euler_box_and_boundary_counted_again"] is False
    assert record["off_shell_identity"]["jacobian_set_to_one"] is False
    assert record["off_shell_identity"]["ghost_initial_data_added"] is False
    assert record["contact_polynomial"]["scalar_collective_relations_imposed"] is False
    assert record["raw_projection"]["identified_with_minus_delta_Gamma"] is False
    assert record["raw_projection"]["coefficients_chosen"] is False
    assert record["raw_projection"]["renormalized"] is False
    assert record["raw_projection"]["raw_trace_is_the_massless_ward_identity"] is True
    assert record["raw_projection"]["raw_trace_is_not_the_renormalized_dirac_trace"] is True
    assert record["raw_projection"]["ward_trace_gap"] < 1e-12
    assert record["raw_projection"]["weyl_structure"]["max_abs"] > 1e-4
    assert record["new_action_coefficient_selected"] is False
    assert record["evolution_ran"] is False
    assert record["ultraviolet_campaign_ran"] is False
    assert record["sea_comparison"]["P_minus_installed_as_physical_vacuum"] is False
    assert record["sea_comparison"]["P_minus_frobenius_distance"] > 1
    for kind in ("lapse", "Q", "r", "shift"):
        assert record["geometric_variation"]["directions"][kind]["max_abs_gap"] < 1e-8
        assert abs(record["geometric_wick"][kind]["connected"]["real"]) > 0
        assert record["geometric_wick"][kind]["means_added_into_connected"] is False
    assert record["geometric_variation"]["multiplied_energy_over_one_block"] == pytest.approx(4)
    assert record["geometric_variation"]["multiplicity_applied_to_noise"] is False
    assert record["angular_completion"]["unique"] is False
    assert record["angular_completion"]["multiplied_again_by_M"] is False
    assert record["angular_completion"]["degeneracy_2_abs_kappa_state_constructed"] is False
    assert record["angular_completion"]["four_dimensional_angular_vertex"] is None
    assert record["angular_completion"]["conjugation_gap"] < 1e-10
    assert record["angular_completion"]["completed_energy_over_one_block"] == pytest.approx(2)
    assert record["angular_completion"]["completed_states_differ_by"] > 1e-6
    metric = smooth_metric(contact.GRID_POINTS, general=True)
    assert np.allclose(metric.lapse, arrays["lapse"])
    negative = direct_hamiltonian(metric, arrays["shift"], -contact.ANGULAR_LABEL)
    mapped = arrays["sign_map"] @ arrays["hamiltonian"] @ arrays["sign_map"].conj().T
    assert np.max(np.abs(mapped - negative)) < 1e-10
    assert hashlib.sha256(payload).hexdigest() == record["npz_sha256"]
    with np.load(io.BytesIO(payload)) as loaded:
        for name, values in arrays.items():
            assert contact.array_digest(loaded[name]) == record["array_sha256"][name]
            assert np.allclose(loaded[name], values)
    encoded = json.dumps(record, allow_nan=False).encode()
    assert len(payload) + len(encoded) < contact.OUTPUT_LIMIT_BYTES


def test_sign_pair_noise_is_two_copies_not_sixteen(benchmark):
    record, arrays, _payload = benchmark
    covariance = arrays["covariance"]
    kernel = arrays["kernel_lapse"]
    unitary = arrays["sign_map"]
    one = influence.connected(covariance, kernel, kernel)
    completed_covariance = block_diag(covariance, unitary @ covariance @ unitary.conj().T)
    completed_kernel = block_diag(kernel, unitary @ kernel @ unitary.conj().T)
    paired = influence.connected(completed_covariance, completed_kernel, completed_kernel)
    assert record["angular_completion"]["sign_pair_factor"] == 2
    assert paired == pytest.approx(2 * one, abs=1e-10)
    assert paired != pytest.approx(16 * one, abs=1e-6)


def test_publish_is_exclusive_and_a_bad_commit_writes_nothing(tmp_path, monkeypatch):
    destination = tmp_path / "benchmark.json"
    contact.publish_exclusive(destination, b'{"schema": "NSC-DISCOVERY-DIRAC-CONTACT-v1"}\n', b"npz")
    assert destination.is_file() and destination.with_suffix(".npz").is_file()
    with pytest.raises(FileExistsError):
        contact.publish_exclusive(destination, b"{}\n", b"npz")
    with pytest.raises(ValueError):
        contact.write_record(tmp_path / "refused.json", "not-a-commit")
    assert not (tmp_path / "refused.json").exists()
    monkeypatch.setattr(contact, "OUTPUT_LIMIT_BYTES", 8)
    with pytest.raises(ValueError):
        contact.publish_exclusive(tmp_path / "large.json", b"01234567", b"01234567")
    assert not (tmp_path / "large.json").exists()


def test_replay_agrees_except_for_cpu_time():
    first = contact.calculate()
    second = contact.calculate()
    contact.same_calculation(first, second)
