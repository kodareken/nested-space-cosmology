"""DOP853 step control on the homogeneous-plus-difference block only.

The joint solver concatenates retarded tangents onto (A, D). Its RMS error
norm divides by the square root of that whole length, and the initial step
uses the same RMS. Extra tangent directions then loosen the declared rtol on
the primal block and can drop an accepted step. Matter N amplifies that
primal change. evolve_primal_controlled requests the envelope's primal-block
controller. Optional tangent_rtol/tangent_atol are a separate accuracy
control and cannot relax the primal prefix.
"""
import numpy as np
from scipy.integrate import DOP853

from .nsc_ks_difference_envelope import evolve_ks_difference_envelope


def _rms(value):
    value = np.asarray(value)
    return float(np.linalg.norm(value) / value.size ** 0.5)


def _free(func, name):
    found = {} if func.__closure__ is None else dict(zip(func.__code__.co_freevars, func.__closure__))
    if name not in found:
        raise ValueError('difference-envelope primal block is not visible on this right-hand side')
    return found[name].cell_contents


def primal_block_length(fun, y0):
    """Prefix length of (A, D) in the difference-envelope state."""
    unpack = _free(fun, 'unpack')
    size_a = int(_free(unpack, 'size_A'))
    size_d = int(_free(unpack, 'size_D'))
    ndir = int(_free(unpack, 'ndir'))
    if size_a <= 0 or size_d <= 0 or ndir < 0:
        raise ValueError('difference-envelope block sizes must be nonnegative')
    n_primal = size_a + size_d
    if size_a + (1 + ndir) * size_d != int(np.asarray(y0).size) or n_primal > int(np.asarray(y0).size):
        raise ValueError('difference-envelope state does not split into A, D and tangents')
    return n_primal


def primal_initial_step(fun, t0, y0, t_bound, max_step, f0, direction, order, rtol, atol, n_primal):
    """Hairer initial step with the RMS taken on the primal prefix."""
    y0 = np.asarray(y0)
    f0 = np.asarray(f0)
    n_primal = int(n_primal)
    y0p, f0p = y0[:n_primal], f0[:n_primal]
    interval_length = abs(t_bound - t0)
    if interval_length == 0.0:
        return 0.0
    scale = atol + np.abs(y0p) * rtol
    d0 = _rms(y0p / scale)
    d1 = _rms(f0p / scale)
    if d0 < 1e-5 or d1 < 1e-5:
        h0 = 1e-6
    else:
        h0 = 0.01 * d0 / d1
    h0 = min(h0, interval_length)
    y1 = y0 + h0 * direction * f0
    f1 = np.asarray(fun(t0 + h0 * direction, y1))
    d2 = _rms((f1[:n_primal] - f0p) / scale) / h0
    if d1 <= 1e-15 and d2 <= 1e-15:
        h1 = max(1e-6, h0 * 1e-3)
    else:
        h1 = (0.01 / max(d1, d2)) ** (1 / (order + 1))
    return min(100 * h0, h1, interval_length, max_step)


class PrimalBlockDOP853(DOP853):
    """Same tableau; the primal error norm and initial step exclude tangents.

    Optional tangent_rtol/tangent_atol never relax the primal prefix. They can
    only shrink a step when the tangent block itself exceeds that separate scale.
    """

    def __init__(self, fun, t0, y0, t_bound, n_primal, max_step=np.inf,
                 rtol=1e-3, atol=1e-6, vectorized=False, first_step=None,
                 tangent_rtol=None, tangent_atol=None, tangent_block_size=None,
                 **extraneous):
        y0 = np.asarray(y0)
        n_primal = int(n_primal)
        if n_primal <= 0 or n_primal > y0.shape[0]:
            raise ValueError('primal block must be a nonempty prefix of the state')
        if (tangent_rtol is None) != (tangent_atol is None):
            raise ValueError('tangent_rtol and tangent_atol must be supplied together')
        self._n_primal = n_primal
        rest = y0.shape[0] - n_primal
        if tangent_block_size is None:
            tangent_block_size = rest or 1
        if (isinstance(tangent_block_size, bool) or int(tangent_block_size) != tangent_block_size
                or tangent_block_size <= 0 or rest % int(tangent_block_size)):
            raise ValueError('positive tangent block size must divide the tangent state')
        self._tangent_block_size = int(tangent_block_size)
        self._tangent_rtol = None if tangent_rtol is None else float(tangent_rtol)
        self._tangent_atol = None if tangent_atol is None else float(tangent_atol)
        if self._tangent_rtol is not None and (self._tangent_rtol <= 0 or self._tangent_atol < 0
                                               or not np.isfinite([self._tangent_rtol, self._tangent_atol]).all()):
            raise ValueError('finite tangent_rtol>0 and tangent_atol>=0 required')
        self.primal_error_norm = None
        self.tangent_error_norm = None
        if first_step is None:
            f0 = np.asarray(fun(t0, y0))
            direction = np.sign(t_bound - t0) if t_bound != t0 else 1
            first_step = primal_initial_step(
                fun, t0, y0, t_bound, max_step, f0, direction,
                self.error_estimator_order, rtol, atol, n_primal)
        super().__init__(
            fun, t0, y0, t_bound, max_step=max_step, rtol=rtol, atol=atol,
            vectorized=vectorized, first_step=first_step, **extraneous)

    def _estimate_error_norm(self, K, h, scale):
        n = self._n_primal
        primal = float(DOP853._estimate_error_norm(self, K[:, :n], h, scale[:n]))
        self.primal_error_norm = primal
        rest = K.shape[1] - n
        if rest <= 0:
            self.tangent_error_norm = 0.0
            return primal
        if self._tangent_rtol is None:
            self.tangent_error_norm = None
            return primal
        atol = np.asarray(self.atol)
        rtol = np.asarray(self.rtol)
        atol_tail = float(atol) if atol.ndim == 0 else atol[n:]
        rtol_tail = float(rtol) if rtol.ndim == 0 else rtol[n:]
        mag = (np.asarray(scale[n:], float) - atol_tail) / rtol_tail
        mag = np.maximum(mag, 0.0)
        tangent_scale = self._tangent_atol + mag * self._tangent_rtol
        # Each retarded direction gets its own error norm. Appending zero
        # directions must not dilute the tolerance of a nonzero direction.
        width = self._tangent_block_size
        tangent = max(float(DOP853._estimate_error_norm(
            self, K[:, start:start + width], h,
            tangent_scale[start-n:start-n + width]))
            for start in range(n, K.shape[1], width))
        self.tangent_error_norm = tangent
        return max(primal, tangent)


def _factory(fun, t0, y0, t_bound, **kwargs):
    y0 = np.asarray(y0)
    return PrimalBlockDOP853(fun, t0, y0, t_bound, n_primal=primal_block_length(fun, y0), **kwargs)


def evolve_primal_controlled(*args, **kwargs):
    """evolve_ks_difference_envelope with primal-block step control."""
    if kwargs.get('step_control', 'primal') == 'joint':
        raise ValueError('primal-controlled evolution cannot use joint step control')
    kwargs['step_control'] = 'primal'
    return evolve_ks_difference_envelope(*args, **kwargs)
