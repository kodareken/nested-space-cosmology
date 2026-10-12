"""Direct controls for the incoherent imbalance witness.

The shared witness is computed once. Numerical diagnostics use an upper
tolerance. The physical gate stays open.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest
import sympy as sp

from recursive_horizons.nsc_finite_turnover import _onsite_projectors
from recursive_horizons.nsc_imbalance_turnover import (
    FAIL,
    PASS,
    SCHEMA,
    _matrix_is_zero,
    commutant_projection,
    compute_result,
)

ROOT = Path(__file__).resolve().parents[1]
MODULE_SHA256 = "565221169e8d188eb464928f51f980ae4d96974a635e0662b00110c24039ab32"
MEAN_POPULATIONS = (
    1.3837895948404104,
    1.0058995585808204,
    0.6103108465787691,
)
CHANNEL_CURRENT = 0.0012511140479781876
REPROJECTION_SHIFT = 0.004920298239642384
FROZEN_SHORTCUT = 0.00030421076749148856


def load_driver():
    path = ROOT / "scripts" / "derive_nsc_imbalance_turnover.py"
    spec = importlib.util.spec_from_file_location("derive_nsc_imbalance_turnover", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def algebraic_float(text):
    value = sp.N(sp.sympify(text), 40)
    assert abs(float(sp.im(value))) <= 1e-28
    return float(sp.re(value))


def pinch_error(operator, state, mean):
    def convert(matrix):
        data = [
            [complex(sp.expand(matrix[row, col]).evalf(40)) for col in range(matrix.cols)]
            for row in range(matrix.rows)
        ]
        return np.array(data, dtype=np.complex128)

    generated = convert(operator)
    generated = (generated + generated.conj().T) / 2
    _values, vectors = np.linalg.eigh(generated)
    rotated = vectors.conj().T @ convert(state) @ vectors
    pinched = vectors @ np.diag(np.diag(rotated)) @ vectors.conj().T
    return float(np.linalg.norm(convert(mean) - pinched, ord="fro"))


def test_core_bytes_stay_at_the_accepted_sha():
    digest = hashlib.sha256(
        (ROOT / "src/recursive_horizons/nsc_imbalance_turnover.py").read_bytes()
    ).hexdigest()
    assert digest == MODULE_SHA256


def test_qsqrt29_projection_regression():
    local = sp.Matrix([[1, sp.I / 5], [-sp.I / 5, 2]])
    _larger, _smaller, plus, minus = _onsite_projectors(local)
    operator = sp.Matrix([
        [sp.Rational(1, 2), sp.Rational(1, 3)],
        [sp.Rational(1, 3), sp.Rational(5, 2)],
    ])
    state = sp.simplify(sp.Rational(3, 5) * plus + sp.Rational(1, 7) * minus)
    projected = commutant_projection(operator, state)
    expected = [
        sp.Rational(13, 35) - sp.Rational(54, 1015) * sp.sqrt(29),
        sp.Rational(36, 1015) * sp.sqrt(29),
    ]
    assert projected["gram_determinant"] != 0
    assert projected["moments_matched"] is True
    assert projected["commutes"] is True
    assert state.has(sp.sqrt(29))
    for index, want in enumerate(expected):
        assert sp.expand(projected["alpha"][index, 0] - want) == 0
    assert _matrix_is_zero(operator * projected["mean"] - projected["mean"] * operator)
    assert pinch_error(operator, state, projected["mean"]) < 1e-8


@pytest.fixture(scope="module")
def result():
    value = compute_result()
    json.dumps(value, sort_keys=True, allow_nan=False)
    return value


def test_witness_invariants_reversal_and_generated_mean(result):
    metrics = result["metrics"]
    checks = result["checks"]
    assert result["schema"] == SCHEMA
    assert result["verdict"] == PASS
    assert result["answer"] == "yes"
    assert result["physical_local_gate"] == "OPEN"
    assert len(checks) == 34
    assert all(value is True for value in checks.values())
    assert metrics["total_content_exact"] == "3"
    assert Fraction(metrics["total_energy_exact"]) == Fraction(99, 16)
    assert metrics["initial_region_populations_exact"] == ["3/2", "1", "1/2"]
    assert checks["initial_interregion_currents_zero"] is True
    assert checks["initial_current_derivative_nonzero"] is True
    assert checks["uniform_preparation_is_stationary_and_current_free"] is True
    assert checks["sign_reversal_negates_adjacent_currents"] is True
    assert checks["currents_scale_exactly_with_delta"] is True
    assert metrics["uniform_adjacent_currents_max_abs_exact"] == "0"
    assert metrics["delta_scale_factor_exact"] == "1/4"
    populations = [Fraction(item) for item in metrics["mean_region_populations_exact"]]
    assert sum(populations) == 3
    assert populations != [Fraction(3, 2), Fraction(1), Fraction(1, 2)]
    for got, want in zip(populations, MEAN_POPULATIONS):
        assert float(got) == pytest.approx(want, rel=1e-12, abs=1e-15)
    current = algebraic_float(metrics["mean_adjacent_currents_exact"]["0:plus->1:plus"])
    assert current == pytest.approx(CHANNEL_CURRENT, rel=1e-9, abs=1e-15)
    assert current > 0
    dephasing = metrics["spectral_dephasing_frobenius_error"]
    assert isinstance(dephasing, float)
    assert 0 <= dephasing < 1e-30
    assert metrics["spectral_dephasing_digits"] == 157
    assumptions = " ".join(result["assumptions"])
    assert "Cesaro average" in assumptions
    assert "not a claim that C(t) itself relaxes pointwise" in assumptions
    assert "not a stress tensor" in assumptions
    assert "geometry is an input" in assumptions


def test_finite_time_bound_domain(result):
    bound = result["finite_time_bound"]
    assert result["finite_time_bound_error"] is None
    assert bound["obtained"] is True
    assert bound["reading"] == "finite-window averaging error bound, not a physical cooling law"
    gap = Fraction(bound["gamma_lower_exact"])
    deviation = Fraction(bound["deviation_frobenius_upper_exact"])
    operator = Fraction(bound["current_operator_frobenius_upper_exact"])
    lower = Fraction(bound["current_abs_lower_exact"])
    threshold = Fraction(bound["sign_separation_if_T_greater_than"])
    square = Fraction(bound["deviation_frobenius_square_exact"])
    assert gap > 0 and lower > 0 and threshold > 0
    assert deviation * deviation >= square
    assert threshold == 2 * deviation * operator / (gap * lower)
    assert 2 * deviation * operator / (gap * (2 * threshold)) < lower
    intervals = [(Fraction(left), Fraction(right)) for left, right in bound["isolating_intervals"]]
    assert len(intervals) == 6
    assert all(left <= right for left, right in intervals)
    separations = [intervals[index + 1][0] - intervals[index][1] for index in range(5)]
    assert min(separations) == gap
    assert all(item > 0 for item in separations)
    coefficients = [Fraction(item) for item in bound["characteristic_polynomial_high_to_low"]]
    assert coefficients[0] == 1 and len(coefficients) == 7
    assert bound["state_error_at_most"] == "2*deviation_upper/(gamma*T)"
    assert bound["mean_current_error_at_most"] == "2*deviation_upper*current_operator_upper/(gamma*T)"
    channel = bound["current_channel"]
    assert result["metrics"]["mean_adjacent_currents_exact"][channel] == bound["current_exact"]


def test_frozen_state_is_not_the_reprojected_response(result):
    metrics = result["metrics"]
    probe = result["inputs"]["exterior_probe"]
    assert probe == {"index": [2, 2], "delta": "1/8"}
    reprojected = metrics["perturbed_versus_baseline_response_frobenius"]
    frozen = metrics["frozen_cbar_shortcut_versus_baseline_frobenius"]
    frozen_gap = metrics["frozen_cbar_shortcut_versus_recomputed_frobenius"]
    assert reprojected == pytest.approx(REPROJECTION_SHIFT, rel=1e-9, abs=1e-12)
    assert frozen == pytest.approx(FROZEN_SHORTCUT, rel=1e-9, abs=1e-12)
    assert frozen_gap > 1e-3
    assert reprojected > frozen


def test_filling_failure_retains_primary_results(result, monkeypatch):
    import recursive_horizons.nsc_imbalance_turnover as imbalance

    def fail_filling(*_args, **_kwargs):
        raise RuntimeError("forced filling failure")

    monkeypatch.setattr(imbalance, "_fill_currents", fail_filling)
    failed = imbalance.compute_result()
    json.dumps(failed, sort_keys=True, allow_nan=False)
    assert failed["schema"] == SCHEMA
    assert failed["verdict"] == FAIL
    assert failed["answer"] == "unresolved"
    assert failed["verdict"] != PASS
    assert failed["physical_local_gate"] == "OPEN"
    assert failed["filling_control_error"].startswith("RuntimeError:")
    assert failed["checks"]["filling_map_reproduces_imbalance_currents"] is False
    assert failed["metrics"]["mode_filling_to_mean_current_map"] == {}
    for key in (
        "total_content_exact",
        "total_energy_exact",
        "initial_region_populations_exact",
        "mean_region_populations_exact",
        "mean_adjacent_currents_exact",
        "projection_coefficients_exact",
    ):
        assert failed["metrics"][key] == result["metrics"][key]
    assert failed["finite_time_bound"]["gamma_lower_exact"] == result["finite_time_bound"]["gamma_lower_exact"]
    assert failed["finite_time_bound"]["current_exact"] == result["finite_time_bound"]["current_exact"]
    assert "optional_pass" not in failed


def test_record_rejects_modified_or_missing_source(result, tmp_path, monkeypatch):
    driver = load_driver()
    import recursive_horizons.nsc_imbalance_turnover as imbalance

    monkeypatch.setattr(imbalance, "compute_result", lambda: result)
    target = tmp_path / "nsc-imbalance-turnover-v1.json"
    monkeypatch.setattr(driver, "OUTPUT", target)
    assert driver.SCHEMA == SCHEMA
    assert driver.main(["--record"]) == 0
    assert driver.main(["--check"]) == 0
    with pytest.raises(FileExistsError, match="use --check"):
        driver.main(["--record"])

    saved = json.loads(target.read_text())
    assert saved["source_bindings"]["src/recursive_horizons/nsc_imbalance_turnover.py"] == MODULE_SHA256
    saved["source_bindings"]["src/recursive_horizons/nsc_imbalance_turnover.py"] = "0" * 64
    target.write_bytes(driver.canonical(saved))
    with pytest.raises(ValueError, match="source hash mismatch"):
        driver.main(["--check"])

    fresh = driver.payload()
    missing = json.loads(driver.canonical(fresh))
    del missing["source_bindings"]["src/recursive_horizons/nsc_nested_qualities.py"]
    target.write_bytes(driver.canonical(missing))
    with pytest.raises(FileNotFoundError, match="missing bound source"):
        driver.main(["--check"])

    shifted = json.loads(driver.canonical(fresh))
    shifted["metrics"]["total_content_exact"] = "4"
    target.write_bytes(driver.canonical(shifted))
    with pytest.raises(ValueError, match="differs from recomputation"):
        driver.main(["--check"])
