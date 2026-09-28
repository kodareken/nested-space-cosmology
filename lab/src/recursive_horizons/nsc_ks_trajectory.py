"""Capture the owned KS solve as anchored local reconstruction polynomials.

SciPy supplies DOP853 dense coefficients. The same evolution code object runs
in a private namespace with a recording solver factory; no global is patched.
Dense polynomials are approximations, not continuous residual enclosures.
"""
from dataclasses import dataclass
from types import FunctionType

import numpy as np
from scipy.integrate import DOP853

from .nsc_evolved_incoming_state import _finite_array
from .nsc_ks_difference_envelope import evolve_ks_difference_envelope
from . import nsc_ks_source_envelope as E


@dataclass(frozen=True)
class TrajectorySegment:
    rho_start: float
    rho_end: float
    start: object
    end: object
    coefficients: object

    def __post_init__(self):
        if not np.isfinite([self.rho_start, self.rho_end]).all() or self.rho_start <= self.rho_end:
            raise ValueError('decreasing finite rho segment required')
        for name in ('start', 'end', 'coefficients'):
            object.__setattr__(self, name, _finite_array(getattr(self, name), complex, name))
        if self.start.ndim != 1 or self.end.shape != self.start.shape or self.coefficients.shape != (6, len(self.start)):
            raise ValueError('six DOP853 correction coefficients and matching endpoint vectors required')

    def evaluate(self, fraction, *, derivative=False):
        """Degree-seven polynomial; mathematical endpoints are stored exactly."""
        x = float(fraction)
        if not np.isfinite(x) or not 0 <= x <= 1:
            raise ValueError('segment fraction must lie in [0,1]')
        if derivative:
            value = (self.end - self.start).copy()
        else:
            value = (1 - x) * self.start + x * self.end
        for j, coefficient in enumerate(self.coefficients, 1):
            p, q = (j + 2) // 2, (j + 1) // 2
            factor = (p * x**(p-1) * (1-x)**q - q * x**p * (1-x)**(q-1)
                      if derivative else x**p * (1-x)**q)
            value += factor * coefficient
        return value / (self.rho_end - self.rho_start) if derivative else value


def capture_ks_trajectory(*args, **kwargs):
    """Same upstream/history solve with local dense-output observations."""
    segments = []

    class Recorder(DOP853):
        def step(self):
            message = super().step()
            if self.status != 'failed':
                dense = self.dense_output()
                if dense.F.shape != (7, self.n):
                    raise ValueError('unsupported SciPy DOP853 reconstruction layout')
                segments.append(TrajectorySegment(self.t_old, self.t,
                    self.y_old, self.y, dense.F[1:]))
            return message

    original = evolve_ks_difference_envelope
    namespace = {**original.__globals__, 'DOP853': Recorder}
    runner = FunctionType(original.__code__, namespace, original.__name__,
                          original.__defaults__, original.__closure__)
    runner.__kwdefaults__ = dict(original.__kwdefaults__ or {})
    prepared = runner(*args, **kwargs)
    if not segments or segments[0].rho_start != prepared.rho_up or segments[-1].rho_end != 1.:
        raise ValueError('captured trajectory does not cover the prepared slab')
    for previous, following in zip(segments, segments[1:]):
        if previous.rho_end != following.rho_start or not np.array_equal(previous.end, following.start):
            raise ValueError('captured reconstruction endpoints are not continuous')
    return prepared, tuple(segments)


class TrajectoryResidualSampler:
    """Actual PDE residual and its first two axial derivatives at sample nodes.

    Results are diagnostics. The sampled periodic quadrature is not promoted
    to a continuous norm enclosure.
    """
    def __init__(self, prepared, family, oversampling=2):
        prepared.require_history(family)
        if isinstance(oversampling, bool) or not isinstance(oversampling, int) or oversampling < 1:
            raise ValueError('positive integer oversampling required')
        self.prepared = prepared
        _, self.metric = E._as_metric(family)
        self.native = np.asarray(prepared.binding.computational_z)
        self.length = float((self.native[1] - self.native[0]) * len(self.native))
        self.z = E.computational_z_grid(len(self.native) * oversampling, self.length)
        self.profiles = np.array([[[[callback(float(z), j) for z in self.z]
                                   for callback in (d.w, d.U)] for d in self.metric.directions] for j in range(3)])
        self.T1 = E.chart_coordinates(1.)[0]

    def _envelope(self, vector):
        nsrc, nz = len(self.prepared.source_energies), len(self.native)
        size = 2 * nsrc
        reference = vector[:size].reshape(2, nsrc)
        difference = vector[size:size + size*nz].reshape(2, nsrc, nz)
        return reference, difference

    def sample(self, segment, fraction):
        A, D = self._envelope(segment.evaluate(fraction))
        Arho, Drho = self._envelope(segment.evaluate(fraction, derivative=True))
        rho = segment.rho_start + float(fraction) * (segment.rho_end - segment.rho_start)
        dz = float(self.native[1] - self.native[0])
        wave = 2*np.pi*np.fft.fftfreq(len(self.native), d=dz)
        phase = np.exp(1j*np.multiply.outer(self.z - self.native[0], wave))

        def evaluate(field, order):
            modes = np.fft.fft(field, axis=-1) * (1j*wave)**order
            return np.tensordot(modes, phase, axes=(-1, -1)) / len(wave)

        X = [evaluate(D, j) for j in range(4)]
        X[0] += A[:, :, None]
        Xrho = [evaluate(Drho, j) for j in range(3)]
        Xrho[0] += Arho[:, :, None]
        s = E.chart_coordinates(float(rho))[0] - self.T1
        chi = np.array([float(E._plateau(np.array([s]), d.inner_radius, d.outer_radius)[0])
                        for d in self.metric.directions])
        radius = np.array([np.einsum('d,dz->z', np.asarray(self.metric.amplitudes)*chi,
                           s*row[:, 0] + s**3*row[:, 1]/6) for row in self.profiles])
        r_ref = np.sqrt(1 + rho*rho)
        delta_r = radius[0].copy()
        radius[0] += r_ref
        r, rz, rzz = radius
        beta = float(E.geometry(rho)[0])
        axial = np.sqrt(beta*beta - 1)
        G = E.ks_generator(rho, self.prepared.source_energies, self.prepared.mass, self.prepared.angular)
        k = 1j*self.prepared.angular/axial
        potential = (-k*delta_r/(r_ref*r), -k*rz/r**2, k*(2*rz**2/r**3 - rzz/r**2))
        spin = [np.einsum('ab,bsz->asz', E.S2, x) for x in X[:3]]
        residuals = []
        for j in range(3):
            result = Xrho[j] - np.einsum('ab,bsz->asz', E.S3, X[j+1])/axial**2 - np.einsum('sab,bsz->asz', G, X[j])
            result -= potential[0]*spin[j]
            if j >= 1:
                result -= j*potential[1]*spin[j-1]
            if j == 2:
                result -= potential[2]*spin[0]
            residuals.append(result * self.prepared.column_weights[None, :, None])
        norms = [float(np.sqrt(self.length/len(self.z)*np.sum(abs(v)**2))) for v in residuals]
        return {'rho': float(rho), 'weighted_L2_indicators': norms,
                'continuous_integral_bounds': None, 'is_error_certificate': False}
