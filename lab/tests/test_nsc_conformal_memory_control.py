"""Independent triangular memory-control algebra and recorded reference errors."""
import json

import numpy as np
from scipy.integrate import solve_ivp

from recursive_horizons import nsc_conformal_memory_control as memory
from recursive_horizons.nsc_coupled_local_response import sha256_file


def test_triangular_generator_removes_only_retained_to_exterior_response():
    rng = np.random.default_rng(729)
    n, retained = 8, 2
    raw = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    h0 = .5 * (raw + raw.conj().T)
    raw = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    h1 = .5 * (raw + raw.conj().T)
    observer = np.eye(n, dtype=complex)[:, :retained]
    initial = rng.normal(size=(n, 3)) + 1j * rng.normal(size=(n, 3))

    def hamiltonian(mark):
        return h0 + mark * h1

    g = memory.triangular_generator(hamiltonian(.02), observer)
    expected = hamiltonian(.02).copy()
    expected[retained:, :retained] = 0
    np.testing.assert_allclose(g, expected, rtol=0, atol=1e-14)
    assert np.linalg.norm(g - g.conj().T) > 1

    def block_ode(mark, flat):
        state = flat.reshape(initial.shape)
        h = hamiltonian(mark)
        x, y = state[:retained], state[retained:]
        return np.vstack((-1j * (h[:retained, :retained] @ x + h[:retained, retained:] @ y),
                          -1j * (h[retained:, retained:] @ y))).ravel()

    def triangular_ode(mark, flat):
        g = memory.triangular_generator(hamiltonian(mark), observer)
        return (-1j * g @ flat.reshape(initial.shape)).ravel()

    options = {"method": "DOP853", "rtol": 1e-12, "atol": 1e-13, "t_eval": [.02]}
    independent = solve_ivp(block_ode, (0, .02), initial.ravel(), **options)
    actual = solve_ivp(triangular_ode, (0, .02), initial.ravel(), **options)
    assert independent.success and actual.success
    np.testing.assert_allclose(actual.y, independent.y, rtol=0, atol=1e-12)


def test_actual_memory_reference_has_own_error_and_substep_control():
    record = json.loads(memory.OUTPUT_JSON.read_text())
    production = json.loads(memory.PRODUCTION_JSON.read_text())
    assert record["negative_control"] is True
    assert record["physical_Hermitian_model_claimed"] is False
    assert record["geometry_rerun"] is False
    assert record["dense_propagator_history_bytes"] == 0
    assert record["cpu_seconds"] < record["cpu_budget_seconds"] <= 30
    assert record["error_within_one_percent"] is True
    assert record["reference_indicator_within_one_percent"] is True
    assert record["active_initial_C_AE"] == production["initial_blocks"]["C_AE_frobenius"]
    assert record["lower_exterior_from_retained_block_max_abs"] < 1e-12
    assert record["retained_rows_unchanged_max_abs"] < 1e-12
    assert record["source_code_sha256"] == sha256_file(memory.__file__)
    assert record["payload_sha256"] == sha256_file(memory.OUTPUT_NPZ)
    with np.load(memory.PRODUCTION_NPZ, allow_pickle=False) as source, np.load(memory.OUTPUT_NPZ, allow_pickle=False) as reference:
        np.testing.assert_array_equal(source["times"], reference["times"])
        error = float(np.max(abs(source["occupation_memory_off"] - reference["occupation_reference_substeps_4"])))
        indicator = float(np.max(abs(reference["occupation_reference_substeps_2"] - reference["occupation_reference_substeps_4"])))
    assert error == record["streamed_memory_off_error_against_independent_reference"]
    assert indicator == record["reference_substeps_2_to_4_indicator"]
    assert record["error_fraction_of_memory_effect"] == error / record["memory_effect"]
    assert record["reference_indicator_fraction_of_memory_effect"] == indicator / record["memory_effect"]
    for path in (memory.PRODUCTION_JSON, memory.PRODUCTION_NPZ):
        assert sha256_file(path) == record["production_bindings"][path.name]
