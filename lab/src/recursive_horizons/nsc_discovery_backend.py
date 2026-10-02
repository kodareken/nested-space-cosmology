"""Matrix-free FFT carriers for the existing spherical Galerkin action.

The sealed dense builders stay the oracle. This module does not modify
Galerkin, coupling, or nested-pair source, and it does not assemble a dense
rank-6 covariance. Quadrature maps are applied with ``@``. Real geometry
pullback is ``weight * A.T``, and the spinor column adjoint is ``U.conj().T``.

``__array__`` raises for every carrier. A high ``__array_priority__`` makes
NumPy use those reflected operations instead of materializing the map.
Explicit ``to_dense`` is only for small reference checks.

``resolve_pair(pair, backend)`` accepts ``dense``, ``fft``, and ``auto``.
``auto`` keeps the dense oracle unless the recorded production-resolution
gate says the FFT step is at least ``PRODUCTION_MIN_SPEEDUP`` times faster
and the observable comparison passed. That record is filled from
``lab/scripts/benchmark_nsc_discovery_backend.py``; this module does not
write a scientific campaign.
"""
from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import replace
from pathlib import Path

import numpy as np
import scipy.fft

from . import nsc_nested_parent_child as nested
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_episode_assessment as metric
from . import nsc_spherical_galerkin_coupling as galerkin
from .nsc_spherical_feedback_action import partial_F


PRODUCTION_RESOLUTION_NF = 256
PRODUCTION_QUADRATURE_FACTOR = 4
PRODUCTION_MIN_SPEEDUP = 1.5
# One-thread wall clock at nf=256, nq=1024, three RK4 steps plus observation.
# Confirmation run: dense 0.357 s, FFT 0.037 s, speedup 9.55. Threads 2 and 4
# were 8.58 and 9.07. Peak RSS was about 278 MB. Episode drivers pin one thread,
# so auto follows that sample. No campaign file is written.
PRODUCTION_FFT_SPEEDUP = 9.55
PRODUCTION_FFT_EQUIVALENT = True
PRODUCTION_FFT_SELECTED = True
PRODUCTION_THREAD_COUNT = 1

EQUIVALENCE_ATOL = 1e-7
EQUIVALENCE_RTOL = 1e-8
# R_h divides second derivatives by powers of Q. At nf=256 the explicit dense
# Fourier matrix and the Nyquist-zero spectral derivative differ by a few
# 1e-7 here; the FFT carrier tracks spectral_dx much more tightly. The sealed
# curvature note records that dense-versus-FFT gap near 2.6e-6.
CURVATURE_ATOL = 1e-6
_FFT_WORKERS: ContextVar[int] = ContextVar("nsc_discovery_fft_workers", default=1)
_SAVED_EPISODE = (
    Path(__file__).resolve().parents[2]
    / "results"
    / "development"
    / "nsc-spherical-conformal-episode-v2.npz"
)


def _workers():
    return _FFT_WORKERS.get()


@contextmanager
def fft_thread_limit(workers):
    """Limit carrier FFTs in this context. Does not change BLAS."""
    count = int(workers)
    if count < 1:
        raise ValueError("FFT worker count must be positive")
    token = _FFT_WORKERS.set(count)
    try:
        yield
    finally:
        _FFT_WORKERS.reset(token)


def _scalar(value):
    if isinstance(value, (bool, np.bool_)):
        return None
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, (bool, np.bool_)):
        return None
    if isinstance(value, (int, np.integer)):
        number = float(value)
    elif isinstance(value, (float, np.floating)):
        number = float(value)
    elif isinstance(value, (complex, np.complexfloating)):
        number = float(value.real) if value.imag == 0.0 else complex(value)
    elif isinstance(value, np.ndarray) and value.shape == ():
        return _scalar(value.item())
    else:
        return None
    if not np.isfinite(number):
        raise ValueError("operator scale must be finite")
    return number


def _prepare(values, length):
    array = np.asarray(values)
    if array.ndim == 1:
        if array.shape[0] != length:
            raise ValueError(f"spectral operator expected length {length}, got {array.shape[0]}")
        return array[:, None], True
    if array.ndim != 2 or array.shape[0] != length:
        raise ValueError(f"spectral operator expected leading dimension {length}, got {array.shape}")
    return array, False


def _finish(result, squeeze):
    return result[:, 0] if squeeze else result


def _real_carrier(values):
    imag = float(np.max(np.abs(np.imag(values)))) if np.iscomplexobj(values) else 0.0
    scale = max(1.0, float(np.max(np.abs(np.real(values)))))
    if imag > 1e-8 * scale:
        raise ValueError("FFT carrier left the real trigonometric subspace")
    return np.real(values)


class SpectralOperator:
    """A square or rectangular Fourier map with matrix-style ``@`` and adjoint."""

    __array_priority__ = 1000.0

    def __init__(self, shape, dtype):
        rows, cols = (int(shape[0]), int(shape[1]))
        if rows < 1 or cols < 1:
            raise ValueError("spectral operator shape must be positive")
        self.shape = (rows, cols)
        self.dtype = np.dtype(dtype)
        self.ndim = 2

    def apply(self, values):
        raise NotImplementedError

    def apply_adjoint(self, values):
        raise NotImplementedError

    def apply_transpose(self, values):
        return np.conjugate(self.apply_adjoint(np.conjugate(values)))

    def __matmul__(self, other):
        if isinstance(other, SpectralOperator):
            return Composed(self, other)
        array, squeeze = _prepare(other, self.shape[1])
        return _finish(self.apply(array), squeeze)

    def __rmatmul__(self, other):
        if isinstance(other, SpectralOperator):
            return Composed(other, self)
        array = np.asarray(other)
        squeezed = array.ndim == 1
        if squeezed:
            array = array[None, :]
        if array.shape[-1] != self.shape[0]:
            raise ValueError(
                f"right-multiply expected trailing dimension {self.shape[0]}, got {array.shape}"
            )
        moved = np.moveaxis(array, -1, 0)
        image = self.apply_adjoint(np.conjugate(moved))
        result = np.conjugate(np.moveaxis(image, 0, -1))
        return result[0] if squeezed else result

    def __mul__(self, other):
        scale = _scalar(other)
        if scale is None:
            raise TypeError("spectral operators do not broadcast; apply them with @")
        return Scaled(self, scale)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __array__(self, dtype=None, copy=None):
        raise TypeError(
            "refusing to densify a spectral operator of shape "
            f"{self.shape}; apply it with @ or call to_dense explicitly"
        )

    def to_dense(self, *, allow_large=False):
        entries = self.shape[0] * self.shape[1]
        if entries > 65_536 and not allow_large:
            raise TypeError(
                f"refusing to realize a spectral operator of shape {self.shape}"
            )
        identity = np.eye(self.shape[1], dtype=self.dtype)
        return self.apply(identity)

    @property
    def T(self):
        return Transposed(self)

    def conj(self):
        return Conjugated(self)

    def adjoint(self):
        return self.conj().T

    @property
    def H(self):
        return self.adjoint()

    def __repr__(self):
        return f"{type(self).__name__}(shape={self.shape}, dtype={self.dtype.name})"


class Scaled(SpectralOperator):
    def __init__(self, base, scale):
        if isinstance(base, Scaled):
            scale = scale * base.scale
            base = base.base
        super().__init__(base.shape, np.result_type(base.dtype, np.asarray(scale).dtype))
        self.base = base
        self.scale = scale

    def apply(self, values):
        return self.scale * self.base.apply(values)

    def apply_adjoint(self, values):
        return np.conjugate(self.scale) * self.base.apply_adjoint(values)

    def apply_transpose(self, values):
        return self.scale * self.base.apply_transpose(values)


class Transposed(SpectralOperator):
    def __init__(self, base):
        super().__init__((base.shape[1], base.shape[0]), base.dtype)
        self.base = base

    def apply(self, values):
        return self.base.apply_transpose(values)

    def apply_adjoint(self, values):
        return np.conjugate(self.base.apply(np.conjugate(values)))

    @property
    def T(self):
        return self.base


class Conjugated(SpectralOperator):
    def __init__(self, base):
        super().__init__(base.shape, base.dtype)
        self.base = base

    def apply(self, values):
        return np.conjugate(self.base.apply(np.conjugate(values)))

    def apply_adjoint(self, values):
        return self.base.apply_transpose(values)

    def conj(self):
        return self.base


class Composed(SpectralOperator):
    def __init__(self, left, right):
        if left.shape[1] != right.shape[0]:
            raise ValueError(f"cannot compose shapes {left.shape} and {right.shape}")
        super().__init__((left.shape[0], right.shape[1]), np.result_type(left.dtype, right.dtype))
        self.left = left
        self.right = right

    def apply(self, values):
        return self.left.apply(self.right.apply(values))

    def apply_adjoint(self, values):
        return self.right.apply_adjoint(self.left.apply_adjoint(values))


class PeriodicDerivative(SpectralOperator):
    """Real periodic derivative. The even Nyquist symbol is zero, as in the dense oracle."""

    def __init__(self, points, length):
        points = int(points)
        length = float(length)
        if points < 4 or points % 2:
            raise ValueError("periodic derivative expects an even point count")
        if not np.isfinite(length) or length <= 0.0:
            raise ValueError("positive period required")
        super().__init__((points, points), np.float64)
        self.points = points
        self.length = length
        symbol = 2.0 * np.pi * np.arange(points // 2 + 1) / length
        symbol[-1] = 0.0
        self._symbol = symbol

    def apply(self, values):
        array, squeeze = _prepare(values, self.points)
        if np.isrealobj(array):
            spectrum = scipy.fft.rfft(np.real(array), axis=0, workers=_workers())
            symbol = self._symbol.reshape((-1, 1))
            derived = scipy.fft.irfft(spectrum * (1j * symbol), n=self.points, axis=0, workers=_workers())
            return _finish(derived, squeeze)
        spectrum = scipy.fft.fft(array, axis=0, workers=_workers())
        modes = np.fft.fftfreq(self.points) * self.points
        wavenumber = 2.0 * np.pi * modes / self.length
        wavenumber = np.array(wavenumber, dtype=float)
        wavenumber[self.points // 2] = 0.0
        derived = scipy.fft.ifft(spectrum * wavenumber[:, None] * 1j, axis=0, workers=_workers())
        return _finish(derived, squeeze)

    def apply_adjoint(self, values):
        return -self.apply(values)


class AntiperiodicMomentum(SpectralOperator):
    """Half-integer AP momentum, the eta=1/2 twist of the dense covariant symbol."""

    def __init__(self, points, length):
        points = int(points)
        length = float(length)
        if points < 10 or points % 2:
            raise ValueError("AP fermion grid must be even and at least 10 points")
        if not np.isfinite(length) or length <= 0.0:
            raise ValueError("positive period required")
        super().__init__((points, points), np.complex128)
        self.points = points
        self.length = length
        index = np.arange(points)
        self._twist = np.exp(-1j * np.pi * index / points)
        self._untwist = np.exp(1j * np.pi * index / points)
        modes = np.fft.fftfreq(points) * points
        self._symbol = 2.0 * np.pi * (modes + 0.5) / length

    def apply(self, values):
        array, squeeze = _prepare(values, self.points)
        spectrum = scipy.fft.fft(array * self._twist[:, None], axis=0, workers=_workers())
        twisted = scipy.fft.ifft(spectrum * self._symbol[:, None], axis=0, workers=_workers())
        return _finish(twisted * self._untwist[:, None], squeeze)

    def apply_adjoint(self, values):
        return self.apply(values)


class GeometryProlongation(SpectralOperator):
    """Odd-band periodic interpolation A_g, shape (nq, ng). No Nyquist mode."""

    def __init__(self, ng, nq, length):
        ng, nq = int(ng), int(nq)
        length = float(length)
        if ng < 3 or ng % 2 == 0:
            raise ValueError("geometry degree count must be odd so the Nyquist mode is absent")
        if nq < ng or nq % 2:
            raise ValueError("quadrature grid must be even and at least as fine as the geometry")
        if not np.isfinite(length) or length <= 0.0:
            raise ValueError("positive period required")
        super().__init__((nq, ng), np.float64)
        self.ng = ng
        self.nq = nq
        self.length = length
        self.half = ng // 2

    def apply(self, values):
        array, squeeze = _prepare(values, self.ng)
        coefficients = scipy.fft.fft(array, axis=0, workers=_workers()) / self.ng
        spectrum = np.zeros((self.nq, array.shape[1]), dtype=np.complex128)
        spectrum[: self.half + 1] = coefficients[: self.half + 1]
        spectrum[-self.half :] = coefficients[-self.half :]
        synthesized = scipy.fft.ifft(spectrum * self.nq, axis=0, workers=_workers())
        return _finish(_real_carrier(synthesized), squeeze)

    def apply_adjoint(self, values):
        """Unweighted A_g.T. The Galerkin weight is a separate scalar factor."""
        array, squeeze = _prepare(values, self.nq)
        spectrum = scipy.fft.fft(array, axis=0, workers=_workers())
        band = np.zeros((self.ng, array.shape[1]), dtype=np.complex128)
        band[: self.half + 1] = spectrum[: self.half + 1]
        band[-self.half :] = spectrum[-self.half :]
        pulled = scipy.fft.ifft(band, axis=0, workers=_workers())
        return _finish(_real_carrier(pulled), squeeze)


class SpinorCarrier(SpectralOperator):
    """Half-integer AP interpolation. ``canonical`` includes sqrt(nf/nq) so U is an isometry."""

    def __init__(self, nf, nq, length, *, canonical):
        nf, nq = int(nf), int(nq)
        length = float(length)
        if nf < 10 or nf % 2:
            raise ValueError("AP fermion count must be even and at least 10")
        if nq < nf or nq % 2:
            raise ValueError("AP quadrature must be even and finer than the fermion band")
        if not np.isfinite(length) or length <= 0.0:
            raise ValueError("positive period required")
        super().__init__((nq, nf), np.complex128)
        self.nf = nf
        self.nq = nq
        self.length = length
        self.canonical = bool(canonical)
        self.scale = float(np.sqrt(nf / nq) if canonical else 1.0)
        coarse = np.arange(nf)
        fine = np.arange(nq)
        self._twist_coarse = np.exp(-1j * np.pi * coarse / nf)
        self._untwist_coarse = np.exp(1j * np.pi * coarse / nf)
        self._twist_fine = np.exp(-1j * np.pi * fine / nq)
        self._untwist_fine = np.exp(1j * np.pi * fine / nq)
        self._modes = np.arange(-nf // 2, nf // 2)

    def apply(self, values):
        array, squeeze = _prepare(values, self.nf)
        coefficients = scipy.fft.fft(array * self._twist_coarse[:, None], axis=0, workers=_workers()) / self.nf
        spectrum = np.zeros((self.nq, array.shape[1]), dtype=np.complex128)
        spectrum[self._modes % self.nq] = coefficients[self._modes % self.nf]
        synthesized = scipy.fft.ifft(spectrum * self.nq, axis=0, workers=_workers())
        return _finish(self.scale * synthesized * self._untwist_fine[:, None], squeeze)

    def apply_adjoint(self, values):
        array, squeeze = _prepare(values, self.nq)
        spectrum = scipy.fft.fft(array * self._twist_fine[:, None], axis=0, workers=_workers())
        gathered = np.zeros((self.nf, array.shape[1]), dtype=np.complex128)
        gathered[self._modes % self.nf] = spectrum[self._modes % self.nq]
        returned = scipy.fft.ifft(gathered, axis=0, workers=_workers())
        return _finish(self.scale * returned * self._untwist_coarse[:, None], squeeze)


def operator_backend(grid):
    momentum = grid.fine.momentum
    if isinstance(momentum, SpectralOperator):
        return "fft"
    return "dense"


def select_backend(backend="auto", *, speedup=None, equivalent=None):
    """Resolve ``dense``, ``fft``, or ``auto`` without touching the oracle builders."""
    name = str(backend).lower()
    if name == "dense":
        return "dense"
    if name == "fft":
        return "fft"
    if name != "auto":
        raise ValueError("backend must be dense, fft, or auto")
    measured = PRODUCTION_FFT_SPEEDUP if speedup is None else speedup
    passed = PRODUCTION_FFT_EQUIVALENT if equivalent is None else equivalent
    if passed and measured is not None and float(measured) >= PRODUCTION_MIN_SPEEDUP:
        return "fft"
    return "dense"


def _shared_grid_identity(reference):
    return {
        "length_density": reference.fine.length_density,
        "shift": reference.fine.shift,
        "occupations": reference.fine.occupations,
        "coefficients": reference.fine.coefficients,
        "calibration": reference.fine.calibration,
        "preparation": reference.fine.preparation,
        "xi_q": reference.xi_q,
        "xi_g": reference.xi_g,
        "xi_f": reference.xi_f,
        "modes_g": reference.modes_g,
        "modes_f": reference.modes_f,
    }


def make_fft_grid(reference_grid):
    """Quadrature operators as FFT carriers. Samples, weights, and calibration are shared."""
    if operator_backend(reference_grid) == "fft":
        return reference_grid
    if not isinstance(reference_grid, galerkin.GalerkinGrid):
        raise TypeError("make_fft_grid expects a GalerkinGrid")
    grid = reference_grid
    if grid.fine.occupations.shape != (6,):
        raise ValueError("source_from_columns stays rank 6; refusing a different occupation width")
    if int(grid.fine.multiplicity) != 4 * int(grid.fine.kappa):
        raise ValueError("multiplicity must remain M=4 kappa")
    if int(grid.fine.points) != int(grid.nq):
        raise ValueError("fine system points and quadrature count differ")
    derivative = PeriodicDerivative(grid.nq, grid.length)
    momentum = AntiperiodicMomentum(grid.nq, grid.length)
    geometry = GeometryProlongation(grid.ng, grid.nq, grid.length)
    columns = SpinorCarrier(grid.nf, grid.nq, grid.length, canonical=True)
    interpolation = SpinorCarrier(grid.nf, grid.nq, grid.length, canonical=False)
    fine = replace(grid.fine, derivative=derivative, momentum=momentum)
    fft_grid = replace(
        grid,
        A_g=geometry,
        A_f=interpolation,
        U_f=columns,
        derivative_on_geometry=Composed(derivative, geometry),
        fine=fine,
    )
    _assert_carrier_isometry(fft_grid)
    return fft_grid


def _assert_carrier_isometry(grid):
    geometry = np.eye(grid.ng)
    geometry_pull = grid.weight * (grid.A_g.T @ (grid.A_g @ geometry))
    geometry_error = float(np.max(np.abs(geometry_pull - geometry)))
    columns = np.eye(grid.nf, dtype=np.complex128)
    column_pull = grid.U_f.conj().T @ (grid.U_f @ columns)
    column_error = float(np.max(np.abs(column_pull - columns)))
    if geometry_error > 1e-8 or column_error > 1e-8:
        raise ValueError(
            f"FFT interpolation isometry failed: geometry {geometry_error}, columns {column_error}"
        )


def make_fft_pair(reference_pair):
    """Same W, source columns, weights, and metadata. Only the grid carriers change."""
    if not isinstance(reference_pair, nested.NestedPair):
        raise TypeError("make_fft_pair expects a NestedPair")
    if operator_backend(reference_pair.grid) == "fft":
        return reference_pair
    return replace(reference_pair, grid=make_fft_grid(reference_pair.grid))


def resolve_grid(grid, backend="auto"):
    selected = select_backend(backend)
    if selected == "dense":
        return grid
    return make_fft_grid(grid)


def resolve_pair(pair, backend="auto"):
    """Return ``pair`` for the dense oracle, or a carrier-replaced pair for FFT.

    ``auto`` follows the recorded production gate. Worker drivers can pass
    ``dense`` or ``fft`` explicitly. State arrays are not stored on the pair
    and are not copied here.
    """
    selected = select_backend(backend)
    if selected == "dense":
        return pair
    return make_fft_pair(pair)


def projected_curvature_jets(grid, state, rate=None, bundle=None):
    """Projected Q jet and R_h used by the sealed conformal curvature postprocessor.

    The formula applies this grid's derivative. The spectral copy uses the
    owned Nyquist-zero helper and does not replace the projected jet.
    """
    if grid.gauge != "conformal":
        raise ValueError("projected curvature jets use the conformal chart L=Q, beta=0")
    if rate is None or bundle is None:
        rate, bundle = galerkin.compose_fine_hamiltonian(grid, state)
    fine = bundle["fine_state"]
    _force_r, f_chi = partial_F(fine.r, grid.fine.A, grid.fine.C_W)
    f_chi = float(f_chi)
    q_dot = galerkin.prolong_geometry(grid, rate.Q)
    p_chi_dot = galerkin.prolong_geometry(grid, rate.p_chi)
    directional = (q_dot * fine.p_chi + fine.Q * p_chi_dot) / (2.0 * f_chi)
    q_ddot = galerkin.prolong_geometry(grid, galerkin.pull_geometry(grid, directional))
    q_x = grid.derivative @ fine.Q
    q_xx = grid.derivative @ q_x
    formula = 2.0 / fine.Q ** 2 * (
        q_xx / fine.Q - q_x ** 2 / fine.Q ** 2 - q_ddot / fine.Q + q_dot ** 2 / fine.Q ** 2
    )
    zero = np.zeros(grid.nq)
    spectral = metric.direct_rh_grid(
        fine.Q, fine.Q, zero, q_dot, q_ddot, grid.length, L_dot=q_dot, beta_dot=zero,
    )["R_h"]
    return {
        "Q_dot": q_dot,
        "Q_ddot": q_ddot,
        "Q_x": q_x,
        "Q_xx": q_xx,
        "R_h": formula,
        "R_h_spectral_dx": spectral,
        "projection_gap": q_ddot - directional,
    }


def normalization_record(grid):
    system = grid.fine
    dx = float(grid.dx_q)
    return {
        "kappa": int(system.kappa),
        "multiplicity": int(system.multiplicity),
        "multiplicity_is_4kappa": bool(int(system.multiplicity) == 4 * int(system.kappa)),
        "dx_q": dx,
        "dx_system": float(system.dx),
        "dx_matches_length_over_nq": bool(np.isclose(system.dx, grid.length / grid.nq)),
        "occupation_rank": int(system.occupations.size),
        "source_layout": "rank-6 columns; dense C is not built",
    }


def evaluate_stage(grid, state, dt):
    """Rates, forces, energy, constraints, curvature jets, and one RK4 step."""
    if not np.isfinite(dt) or dt <= 0.0:
        raise ValueError("a positive finite timestep is required")
    rate, bundle = galerkin.compose_fine_hamiltonian(grid, state)
    source = bundle["source"]
    if source["image0"].shape != (grid.nq, 6) or source["image1"].shape != (grid.nq, 6):
        raise ValueError("column source left the rank-6 carrier")
    diagnostics = galerkin.constraint_diagnostics(grid, state, source, bundle["fine_state"])
    jets = projected_curvature_jets(grid, state, rate, bundle) if grid.gauge == "conformal" else None
    observed = galerkin.observation(grid, state)
    stepped = galerkin.rk4_step(grid, state, float(dt))
    return {
        "rate": rate,
        "source": source,
        "diagnostics": diagnostics,
        "jets": jets,
        "energy": float(galerkin.energy(grid, state)),
        "observation": observed,
        "stepped": stepped,
        "normalization": normalization_record(grid),
        "dt": float(dt),
    }


def evaluate_pair_stage(pair, state, dt):
    """Nested-pair rates, energy, source forces, and one RK4 step. No dense C."""
    if not np.isfinite(dt) or dt <= 0.0:
        raise ValueError("a positive finite timestep is required")
    rate, bundle = nested.rates(pair, state, return_bundle=True)
    source = bundle["source"]
    if source["image0"].shape[-1] != 6:
        raise ValueError("nested source left the rank-6 carrier")
    stepped = nested.rk4_step(pair, state, float(dt))
    return {
        "rate": rate,
        "source": source,
        "energy": float(nested.energy(pair, state)),
        "stepped": stepped,
        "normalization": normalization_record(pair.grid),
        "dt": float(dt),
    }


def _scale(reference):
    array = np.asarray(reference)
    if array.size == 0:
        return 1.0
    return max(1.0, float(np.max(np.abs(array))))


def _gap(actual, reference):
    difference = float(np.max(np.abs(np.asarray(actual) - np.asarray(reference))))
    return difference, _scale(reference)


def within_equivalence(difference, scale, *, atol=EQUIVALENCE_ATOL, rtol=EQUIVALENCE_RTOL):
    return bool(difference <= atol + rtol * scale)


def compare_stages(oracle, candidate):
    """Max-abs movement from the dense oracle. Forces are also reported over dx."""
    report = {}
    failures = []
    equivalent = True

    def keep(name, difference, scale, atol=EQUIVALENCE_ATOL):
        nonlocal equivalent
        accepted = within_equivalence(difference, scale, atol=atol)
        if not accepted:
            failures.append(name)
            equivalent = False
        return accepted

    rate_names = ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1", "force_L", "force_Q", "force_beta")
    for name in rate_names:
        difference, scale = _gap(getattr(candidate["rate"], name), getattr(oracle["rate"], name))
        report[f"rate_{name}"] = difference
        keep(f"rate_{name}", difference, scale)
    dx = float(oracle["normalization"]["dx_q"])
    for name in ("force_L", "force_Q", "force_beta"):
        difference, scale = _gap(
            getattr(candidate["rate"], name) / dx,
            getattr(oracle["rate"], name) / dx,
        )
        report[f"{name}_over_dx"] = difference
        keep(f"{name}_over_dx", difference, scale)
    for name in ("K", "S1", "Pmom", "image0", "image1", "force_L", "force_Q", "force_beta"):
        difference, scale = _gap(candidate["source"][name], oracle["source"][name])
        report[f"source_{name}"] = difference
        keep(f"source_{name}", difference, scale)
    energy_gap = abs(candidate["energy"] - oracle["energy"])
    report["energy"] = float(energy_gap)
    keep("energy", energy_gap, max(1.0, abs(oracle["energy"])))
    for name, value in oracle["diagnostics"].items():
        if isinstance(value, (bool, np.bool_)):
            continue
        if isinstance(value, (float, int, np.floating, np.integer)):
            difference = abs(float(candidate["diagnostics"][name]) - float(value))
            report[f"diagnostic_{name}"] = float(difference)
            keep(f"diagnostic_{name}", difference, max(1.0, abs(float(value))))
    if oracle["jets"] is not None or candidate["jets"] is not None:
        for name in ("Q_dot", "Q_ddot", "Q_x", "Q_xx", "R_h", "R_h_spectral_dx", "projection_gap"):
            difference, scale = _gap(candidate["jets"][name], oracle["jets"][name])
            report[f"jet_{name}"] = difference
            atol = CURVATURE_ATOL if name in ("R_h", "R_h_spectral_dx") else EQUIVALENCE_ATOL
            keep(f"jet_{name}", difference, scale, atol=atol)
    for name in ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1"):
        difference, scale = _gap(getattr(candidate["stepped"], name), getattr(oracle["stepped"], name))
        report[f"step_{name}"] = difference
        keep(f"step_{name}", difference, scale)
    for name, value in oracle["observation"].items():
        if isinstance(value, (bool, np.bool_)):
            matched = bool(value) == bool(candidate["observation"][name])
            report[f"observe_{name}"] = matched
            if not matched:
                failures.append(f"observe_{name}")
                equivalent = False
        elif isinstance(value, (float, int, np.floating, np.integer)):
            difference = abs(float(candidate["observation"][name]) - float(value))
            report[f"observe_{name}"] = float(difference)
            keep(f"observe_{name}", difference, max(1.0, abs(float(value))))
    report["normalization_match"] = candidate["normalization"] == oracle["normalization"]
    if not report["normalization_match"]:
        failures.append("normalization")
    report["failures"] = failures
    report["equivalent"] = bool(equivalent and report["normalization_match"])
    report["oracle_normalization"] = oracle["normalization"]
    report["candidate_normalization"] = candidate["normalization"]
    return report


def compare_pair_stages(oracle, candidate):
    report = {}
    equivalent = True
    for name in nested.STATE_NAMES:
        difference, scale = _gap(getattr(candidate["rate"], name), getattr(oracle["rate"], name))
        report[f"rate_{name}"] = difference
        equivalent = equivalent and within_equivalence(difference, scale)
        difference, scale = _gap(getattr(candidate["stepped"], name), getattr(oracle["stepped"], name))
        report[f"step_{name}"] = difference
        equivalent = equivalent and within_equivalence(difference, scale)
    for name in ("force_L", "force_Q", "force_beta"):
        difference, scale = _gap(getattr(candidate["rate"], name), getattr(oracle["rate"], name))
        report[name] = difference
        equivalent = equivalent and within_equivalence(difference, scale)
    energy_gap = abs(candidate["energy"] - oracle["energy"])
    report["energy"] = float(energy_gap)
    equivalent = equivalent and within_equivalence(energy_gap, max(1.0, abs(oracle["energy"])))
    report["normalization_match"] = candidate["normalization"] == oracle["normalization"]
    report["equivalent"] = bool(equivalent and report["normalization_match"])
    return report


def step_movement(before, after):
    """Infinity-norm change of one RK4 step, per state field."""
    movement = {}
    for name in ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1"):
        movement[name] = float(np.max(np.abs(getattr(after, name) - getattr(before, name))))
    return movement


def manufactured_state(grid, seed=3):
    """Band-limited geometry plus orthonormal AP columns. Not a saved physical packet."""
    phi0, phi1 = galerkin.manufactured_columns(grid.nf, seed=seed)
    state = galerkin.blank_state(grid, phi0, phi1)
    angle = 2.0 * np.pi * grid.xi_g / grid.length
    state.Q = state.Q * (1.0 + 0.01 * np.cos(angle))
    state.r = 1.0 + 0.02 * np.sin(angle)
    state.chi = 0.01 * np.cos(angle)
    state.p_Q = 0.004 * np.sin(2.0 * angle)
    state.p_r = 0.01 * np.cos(2.0 * angle)
    state.p_chi = 0.005 * np.sin(angle)
    return state


def load_saved_conformal_state(nf=PRODUCTION_RESOLUTION_NF, dt=0.0005, time=0.3, path=None):
    """Read one sealed conformal frame. Does not write or evolve the payload."""
    payload = Path(path) if path is not None else _SAVED_EPISODE
    prefix = f"nf{int(nf)}_dt_{float(dt):.4f}"
    with np.load(payload, allow_pickle=False) as stored:
        times = np.asarray(stored[prefix + "_frame_times"], dtype=float)
        matches = np.flatnonzero(np.isclose(times, float(time)))
        if matches.size == 0:
            raise ValueError(f"{prefix} has no frame at T={time}")
        index = int(matches[-1])
        values = {
            name: np.array(stored[f"{prefix}_frames_{name}"][index], copy=True)
            for name in ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1")
        }
    return coupling.CauchyState(**values), float(times[index])
