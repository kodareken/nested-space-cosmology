"""FFT carriers against the untouched dense Galerkin and nested-pair oracle."""
import numpy as np
import pytest

from recursive_horizons import nsc_discovery_backend as backend
from recursive_horizons import nsc_nested_parent_child as nested
from recursive_horizons import nsc_spherical_coupling as coupling
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin


def test_periodic_derivative_matches_dense_and_kills_nyquist():
    points, length = 32, 8.0
    dense = coupling.periodic_derivative(points, length)
    carrier = backend.PeriodicDerivative(points, length)
    rng = np.random.default_rng(4)
    real = rng.normal(size=points)
    imag = rng.normal(size=(points, 3)) + 1j * rng.normal(size=(points, 3))
    assert np.max(np.abs(carrier @ real - dense @ real)) < 1e-12
    assert np.max(np.abs(carrier @ imag - dense @ imag)) < 1e-12
    nyquist = np.cos(np.pi * np.arange(points))
    assert np.max(np.abs(carrier @ nyquist)) < 1e-12
    assert np.max(np.abs(dense @ nyquist)) < 1e-12
    mode = 3
    samples = np.sin(2 * np.pi * mode * np.arange(points) / points)
    expected = (2 * np.pi * mode / length) * np.cos(2 * np.pi * mode * np.arange(points) / points)
    assert np.max(np.abs(carrier @ samples - expected)) < 1e-12
    with pytest.raises(TypeError, match="refusing to densify"):
        np.array(backend.PeriodicDerivative(1024, length))
    realized = carrier.to_dense()
    assert realized.shape == dense.shape
    assert np.max(np.abs(realized - dense)) < 1e-11


def test_antiperiodic_momentum_matches_dense_half_integer_symbol():
    points, length = 32, 8.0
    dense, _metric = coupling.antiperiodic_momentum(points, length)
    carrier = backend.AntiperiodicMomentum(points, length)
    rng = np.random.default_rng(5)
    values = rng.normal(size=(points, 4)) + 1j * rng.normal(size=(points, 4))
    assert np.max(np.abs(carrier @ values - dense @ values)) < 1e-11
    assert np.max(np.abs(carrier.adjoint() @ values - dense.conj().T @ values)) < 1e-11
    mode = 2.5
    wave = np.exp(2j * np.pi * mode * np.arange(points) * (length / points) / length)
    assert np.max(np.abs(carrier @ wave - (2 * np.pi * mode / length) * wave)) < 1e-11
    with pytest.raises(TypeError, match="refusing to densify"):
        np.asarray(backend.AntiperiodicMomentum(1024, length))


def test_geometry_and_spinor_adjoints_match_dense_maps():
    nf, length = 16, 8.0
    ng, nq = nf - 1, 40
    geometry = galerkin.periodic_interpolation(ng, nq, length)
    fermions = galerkin.antiperiodic_interpolation(nf, nq, length)
    columns = galerkin.canonical_column_map(fermions, nf, nq)
    weight = ng / nq
    fft_geometry = backend.GeometryProlongation(ng, nq, length)
    fft_columns = backend.SpinorCarrier(nf, nq, length, canonical=True)
    fft_interpolation = backend.SpinorCarrier(nf, nq, length, canonical=False)
    rng = np.random.default_rng(6)
    coarse = rng.normal(size=(ng, 2))
    fine = rng.normal(size=nq)
    phi = rng.normal(size=(nf, 6)) + 1j * rng.normal(size=(nf, 6))
    lifted = rng.normal(size=(nq, 6)) + 1j * rng.normal(size=(nq, 6))
    assert np.max(np.abs(fft_geometry @ coarse - geometry @ coarse)) < 1e-11
    assert np.max(np.abs(weight * fft_geometry.T @ fine - weight * geometry.T @ fine)) < 1e-11
    assert np.max(np.abs(fft_interpolation @ phi - fermions @ phi)) < 1e-11
    assert np.max(np.abs(fft_columns @ phi - columns @ phi)) < 1e-11
    assert np.max(np.abs(fft_columns.conj().T @ lifted - columns.conj().T @ lifted)) < 1e-11
    restored = weight * fft_geometry.T @ (fft_geometry @ coarse)
    assert np.max(np.abs(restored - coarse)) < 1e-11
    assert np.max(np.abs(fft_columns.conj().T @ (fft_columns @ phi) - phi)) < 1e-11


def test_large_operator_product_does_not_densify_through_array():
    carrier = backend.PeriodicDerivative(1024, 8.0)
    with pytest.raises(TypeError, match="refusing to densify"):
        np.array(carrier, copy=False)
    with pytest.raises(TypeError, match="refusing to realize"):
        carrier.to_dense()
    with pytest.raises(TypeError, match="do not broadcast"):
        _ = carrier * np.ones(1024)
    image = carrier @ np.zeros(1024)
    assert image.shape == (1024,)
    assert np.max(np.abs(image)) == 0.0


def test_fft_grid_and_pair_keep_metadata_state_w_and_source():
    grid = galerkin.build_grid(32, quadrature=128, gauge="conformal")
    before = backend._shared_grid_identity(grid)
    fft_grid = backend.make_fft_grid(grid)
    assert backend.operator_backend(grid) == "dense"
    assert backend.operator_backend(fft_grid) == "fft"
    assert fft_grid.fine.derivative is not grid.fine.derivative
    assert fft_grid.fine.length_density is grid.fine.length_density
    assert fft_grid.fine.shift is grid.fine.shift
    assert fft_grid.fine.occupations is grid.fine.occupations
    assert fft_grid.fine.coefficients is grid.fine.coefficients
    assert fft_grid.fine.calibration is grid.fine.calibration
    assert fft_grid.fine.preparation is grid.fine.preparation
    assert fft_grid.xi_q is grid.xi_q
    assert fft_grid.modes_f is grid.modes_f
    assert fft_grid.dx_q == grid.dx_q
    assert fft_grid.weight == grid.weight
    assert before == backend._shared_grid_identity(grid)
    assert backend.normalization_record(fft_grid)["multiplicity_is_4kappa"] is True
    assert backend.normalization_record(fft_grid)["occupation_rank"] == 6
    pair = nested.build_pair(48, source_layout="separated")
    fft_pair = backend.make_fft_pair(pair)
    assert fft_pair.geometry_map is pair.geometry_map
    assert fft_pair.source_phi0 is pair.source_phi0
    assert fft_pair.source_phi1 is pair.source_phi1
    assert fft_pair.source_columns is pair.source_columns
    assert fft_pair.weights is pair.weights
    assert fft_pair.source_metadata is pair.source_metadata
    assert fft_pair.geometry_metadata is pair.geometry_metadata
    assert backend.resolve_pair(pair, "dense") is pair
    assert backend.resolve_grid(grid, "dense") is grid
    assert backend.operator_backend(backend.resolve_pair(pair, "fft").grid) == "fft"


def test_manufactured_stage_matches_dense_oracle_including_one_rk4():
    grid = galerkin.build_grid(32, quadrature=128, gauge="conformal")
    fft_grid = backend.make_fft_grid(grid)
    state = backend.manufactured_state(grid, seed=3)
    original = state.copy()
    dt, _omega, _quadrature_omega = galerkin.stable_timestep(grid, state)
    report = backend.compare_stages(
        backend.evaluate_stage(grid, state, dt),
        backend.evaluate_stage(fft_grid, state, dt),
    )
    assert report["equivalent"] is True
    assert report["oracle_normalization"]["multiplicity"] == 4 * report["oracle_normalization"]["kappa"]
    assert report["force_Q_over_dx"] < 1e-8
    assert report["jet_R_h"] < 1e-8
    assert report["step_phi0"] < 1e-8
    for name in ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1"):
        assert np.array_equal(getattr(state, name), getattr(original, name))
    pair = nested.build_pair(48, source_layout="separated")
    fft_pair = backend.resolve_pair(pair, "fft")
    encoded, _report = nested.initial_state(pair, solve_constraints=False)
    pair_report = backend.compare_pair_stages(
        backend.evaluate_pair_stage(pair, encoded, dt),
        backend.evaluate_pair_stage(fft_pair, encoded, dt),
    )
    assert pair_report["equivalent"] is True


def test_saved_nf256_frame_matches_dense_oracle():
    grid = galerkin.build_grid(256, quadrature=1024, gauge="conformal")
    fft_grid = backend.make_fft_grid(grid)
    saved, time_value = backend.load_saved_conformal_state()
    manufactured = backend.manufactured_state(grid, seed=3)
    assert time_value == 0.3
    for state in (manufactured, saved):
        dt, _omega, _quadrature_omega = galerkin.stable_timestep(grid, state)
        oracle = backend.evaluate_stage(grid, state, dt)
        candidate = backend.evaluate_stage(fft_grid, state, dt)
        report = backend.compare_stages(oracle, candidate)
        assert report["failures"] == []
        assert report["equivalent"] is True
        assert report["oracle_normalization"]["multiplicity_is_4kappa"] is True
        assert report["oracle_normalization"]["dx_q"] == grid.length / grid.nq
        spectral_gap = float(np.max(np.abs(candidate["jets"]["R_h"] - candidate["jets"]["R_h_spectral_dx"])))
        assert spectral_gap < 1e-10
        assert report["force_Q_over_dx"] < 1e-8
        assert report["rate_force_L"] < 1e-9


def test_auto_backend_follows_the_recorded_production_gate():
    assert backend.select_backend("auto", speedup=2.0, equivalent=True) == "fft"
    assert backend.select_backend("auto", speedup=1.49, equivalent=True) == "dense"
    assert backend.select_backend("auto", speedup=3.0, equivalent=False) == "dense"
    assert backend.select_backend("dense") == "dense"
    assert backend.select_backend("fft") == "fft"
    selected = "fft" if backend.PRODUCTION_FFT_SELECTED else "dense"
    assert backend.select_backend("auto") == selected
    if backend.PRODUCTION_FFT_SELECTED:
        assert backend.PRODUCTION_FFT_EQUIVALENT is True
        assert backend.PRODUCTION_FFT_SPEEDUP >= backend.PRODUCTION_MIN_SPEEDUP
    else:
        speed = backend.PRODUCTION_FFT_SPEEDUP or 0.0
        assert not (backend.PRODUCTION_FFT_EQUIVALENT and speed >= backend.PRODUCTION_MIN_SPEEDUP)


def test_dense_derivative_builder_remains_the_callable_oracle():
    assert coupling.periodic_derivative is coupling.periodic_derivative
    matrix = coupling.periodic_derivative(16, 8.0)
    assert matrix.shape == (16, 16)
    assert np.allclose(matrix, -matrix.T)
