"""Compiled point coefficients for one fixed-cap radial stage.

The DOP853 driver, panel breaks, and adjoint sum stay with the caller.
This module replaces the repeated Python body of that callback for points
outside the near-centre background formula. A batch evaluates adjoint nodes
independently and does not add them.

Small-z background points, where the retained history uses split integrals,
return no row. The caller keeps those points on the existing reference path.
"""
import math

import numpy as np
from numba import njit, prange

from .nsc_atlas_cap_background_jets import SMALL_Z

_MAX_LOCAL = 6
_OK = 0
_NONFINITE = 2
_METRIC = 3
_POTENTIAL = 4
_DOMAIN = 5
_PIECE = 6
_AMBIGUOUS = 7
_CROSSED = 8
_QUIET = 9
_MESSAGES = {
    _NONFINITE: 'nonfinite one-sided history jet',
    _METRIC: 'numerical fixed cap left the positive metric domain',
    _POTENTIAL: 'nonfinite reference potential',
    _DOMAIN: 'profile query exits the supplied finite labels',
    _PIECE: 'history-piece index outside the finite profile',
    _AMBIGUOUS: 'ambiguous history-piece ownership at radial endpoint',
    _CROSSED: 'history-piece selection crosses a resolved knot',
    _QUIET: 'history piece requires a certified quiet extension',
}


def accounted_evaluations(local_count, derivative_count, degree):
    """Python checks replaced by one compiled point, counted before the call."""
    return int(32+18*local_count+160+24+12+derivative_count+48+4*degree)


@njit(cache=True, nogil=True, fastmath=False)
def _positive_gap(z):
    x = 2.0*z
    if x <= 1.0:
        term = 1.0
        tail = 1.0
        for k in range(1, 10):
            term *= x*x/((2*k+2)*(2*k+3))
            tail += term
        return (2.0*z/3.0)*tail/(math.sinh(x)/x)
    return 1.0/z-2.0/math.sinh(x)


@njit(cache=True, nogil=True, fastmath=False)
def _horner_low(coeff, x):
    acc = 0.0
    for i in range(coeff.shape[0]-1, -1, -1):
        acc = acc*x+coeff[i]
    return acc


@njit(cache=True, nogil=True, fastmath=False)
def _horner_high(coeff, length, x):
    acc = 0.0
    for i in range(length):
        acc = acc*x+coeff[i]
    return acc


@njit(cache=True, nogil=True, fastmath=False)
def _actual_piece(knots, onset_index, label):
    found = int(np.searchsorted(knots, label, side='right')-1)
    last = knots.shape[0]-2
    if found > last:
        found = last
    if onset_index >= 1 and 0 <= found < onset_index:
        found = -1
    return found


@njit(cache=True, nogil=True, fastmath=False)
def _piece_value(coeff, length, knots, piece, label):
    if piece < 0:
        return 0.0
    return _horner_high(coeff[:, piece], length, label-knots[piece])


@njit(cache=True, nogil=True, fastmath=False)
def _owned_value(coeff, length, knots, onset_index, quiet, selected, label):
    if label > knots[-1] or (quiet == 0 and label < knots[0]):
        return _DOMAIN, 0.0
    found = _actual_piece(knots, onset_index, label)
    if selected < -1 or selected > knots.shape[0]-2:
        return _PIECE, 0.0
    if selected == found:
        return _OK, _piece_value(coeff, length, knots, selected, label)
    # {-1, onset} is the certified quiet edge. Any other pair must be adjacent.
    if selected == -1 and found == onset_index and onset_index >= 1:
        quiet_onset = True
    elif found == -1 and selected == onset_index and onset_index >= 1:
        quiet_onset = True
    else:
        quiet_onset = False
    if not quiet_onset and abs(selected-found) != 1:
        return _AMBIGUOUS, 0.0
    knot = knots[max(selected, found)]
    tolerance = 64.0*np.spacing(max(1.0, abs(label), abs(knot)))
    if abs(label-knot) > tolerance:
        return _CROSSED, 0.0
    if selected == -1 and quiet == 0:
        return _QUIET, 0.0
    return _OK, _piece_value(coeff, length, knots, selected, label)


@njit(cache=True, nogil=True, fastmath=False)
def _background_limits(coeff, lengths, knots, onset_index, quiet, selected, time, z):
    """F and its first two derivatives at T-z and T+z, with panel ownership."""
    out = np.empty(6)
    labels = (time-z, time+z)
    for j in range(2):
        for derivative in range(3):
            status, value = _owned_value(
                coeff[derivative], int(lengths[derivative]), knots, onset_index, quiet,
                int(selected[j]), labels[j])
            if status != _OK:
                return status, out
            out[derivative*2+j] = value
    if not np.isfinite(out).all():
        return _NONFINITE, out
    return _OK, out


@njit(cache=True, nogil=True, fastmath=False)
def _potential(f, alpha):
    scale = math.sqrt(alpha/6.0)
    if alpha > 0.0 and scale == 0.0:
        return _POTENTIAL, 0.0, 0.0
    q = scale*f
    ch = math.cosh(q)
    sh = math.sinh(q)
    sinc = sh/q if q != 0.0 else 1.0
    P = f*f*sinc*sinc*((ch*ch+1.0)/2.0)
    U = 2.0*f*sinc*ch*ch*ch
    if not math.isfinite(P) or not math.isfinite(U):
        return _POTENTIAL, 0.0, 0.0
    return _OK, P, U


@njit(cache=True, nogil=True, fastmath=False)
def _point(panel, point, mass, alpha, count, cut, has_background, time, onset_index, quiet,
           use_selected, breaks, owners, left, width, n_local, indices, powers, deltas,
           shapes, shape_deltas, knots, selected, coeff, lengths, out):
    """One coefficient row. ``out`` is r, jacobian, A, a, rho, then partials, then saw."""
    last = breaks.shape[0]-1
    if point == breaks[last]:
        xi = cut
    else:
        xi = point-owners[panel]
    cell = owners[panel]
    z = left[cell]+width[cell]*xi
    if not (0.0 <= xi <= 1.0) or not math.isfinite(z):
        return _NONFINITE
    local = int(n_local[cell])
    value = np.empty(_MAX_LOCAL)
    fb = np.empty(_MAX_LOCAL)
    hb = np.empty(_MAX_LOCAL)
    for j in range(local):
        value[j] = _horner_low(shapes[cell, j], xi)
    chi = _horner_low(powers[cell, 0], xi)
    k = _horner_low(powers[cell, 1], xi)
    r = math.tanh(z)
    A0 = 1.0/(math.cosh(z)*math.cosh(z))
    if cell == 0:
        qfield = powers[cell, 0, 1]+xi*xi*(powers[cell, 0, 3]+xi*xi*powers[cell, 0, 5])
        ratio = z/r if z != 0.0 else 1.0
        f = qfield*ratio/width[cell]
        gap = _positive_gap(z) if z != 0.0 else 0.0
        h = ((2.0*powers[cell, 0, 3]*xi*xi+4.0*powers[cell, 0, 5]*xi**4)/width[cell]+gap*chi)/A0
        for j in range(local):
            q = shapes[cell, j, 1]+xi*xi*(shapes[cell, j, 3]+xi*xi*shapes[cell, j, 5])
            fb[j] = q*ratio/width[cell]
            hb[j] = ((2.0*shapes[cell, j, 3]*xi*xi+4.0*shapes[cell, j, 5]*xi**4)
                     /width[cell]+gap*value[j])/A0
    else:
        gap = _positive_gap(z)
        f = chi/r
        h = (_horner_low(deltas[cell], xi)/z+gap*chi)/A0
        for j in range(local):
            fb[j] = value[j]/r
            hb[j] = (_horner_low(shape_deltas[cell, j], xi)/z+gap*value[j])/A0
    p = k/A0
    f_T = 0.0
    p_T = 0.0
    h_T = 0.0
    if has_background == 1:
        if use_selected == 1:
            chosen = selected[panel]
        else:
            chosen = np.empty(2, dtype=np.int64)
            for j in range(2):
                label = time-z if j == 0 else time+z
                if label > knots[-1] or (quiet == 0 and label < knots[0]):
                    return _DOMAIN
                chosen[j] = _actual_piece(knots, onset_index, label)
        status, limits = _background_limits(
            coeff, lengths, knots, onset_index, quiet, chosen, time, z)
        if status != _OK:
            return status
        fm, fp = limits[0], limits[1]
        pm, pp = limits[2], limits[3]
        sm, sp = limits[4], limits[5]
        bg_chi = fp-fm
        bg_k = pp-pm
        bg_b = sp-sm
        bg_az = sp+sm
        bg_f = bg_chi/r
        bg_p = bg_k/A0
        bg_h = (pp+pm)/A0-bg_f
        f += bg_f
        p += bg_p
        h += bg_h
        f_T = bg_k/r
        p_T = bg_b/A0
        h_T = bg_az/A0-bg_k/r
    if not (math.isfinite(f) and math.isfinite(p) and math.isfinite(h)):
        return _NONFINITE
    jacobian = width[cell]*A0
    partials = out[5:5+2*count+(1 if has_background == 1 else 0)]
    for i in range(partials.shape[0]):
        partials[i] = 0.0
    out[0] = r
    out[1] = jacobian
    saw = 1.0 if (f != 0.0 or p != 0.0 or h != 0.0) else 0.0
    out[-1] = saw
    if r == 0.0:
        out[2] = 1.0
        out[3] = 0.0
        out[4] = 0.0
        return _OK
    status, P, U = _potential(f, alpha)
    if status != _OK:
        return status
    A = A0-2.0*alpha*mass/r
    if not math.isfinite(A) or A <= 0.0:
        return _METRIC
    square = p*p+h*h
    a = alpha*square/(2.0*r)
    rho = r*r*P/2.0+A*square/4.0
    if not (math.isfinite(a) and math.isfinite(rho)):
        return _NONFINITE
    for j in range(local):
        slot = int(indices[cell, j])
        partials[slot] = r*r*U*fb[j]/2.0+A*h*hb[j]/2.0
        partials[count+slot] = A*p*value[j]/(A0*2.0)
    if has_background == 1:
        partials[partials.shape[0]-1] = r*r*U*f_T/2.0+A*(p*p_T+h*h_T)/2.0
    out[2] = A
    out[3] = a
    out[4] = rho
    return _OK


@njit(cache=True, nogil=True, fastmath=False, parallel=True)
def _points(panels, points, masses, alpha, count, cut, has_background, time, onset_index, quiet,
            use_selected, breaks, owners, left, width, n_local, indices, powers, deltas,
            shapes, shape_deltas, knots, selected, coeff, lengths, out, status):
    """Independent nodes. The caller keeps the ordered adjoint sum."""
    for i in prange(panels.shape[0]):
        status[i] = _point(
            int(panels[i]), points[i], masses[i], alpha, count, cut, has_background, time,
            onset_index, quiet, use_selected, breaks, owners, left, width, n_local, indices,
            powers, deltas, shapes, shape_deltas, knots, selected, coeff, lengths, out[i])


class CompiledCapCoefficients:
    """Stage-local polynomials and one compiled coefficient row.

    ``evaluate`` returns None for a near-centre background point. Those points
    stay on the reference implementation. ``evaluate_batch`` rejects that
    region because a mixed batch would hide which equation was used.
    """

    def __init__(self, fields, breaks, owners, cut, alpha, count, *, background=None,
                 time=None, panel_indices=None, budget):
        budget.check(64)
        self.count = int(count)
        self.alpha = float(alpha)
        self.cut = float(cut)
        self.breaks = np.ascontiguousarray(breaks, dtype=np.float64)
        self.owners = np.ascontiguousarray(owners, dtype=np.int64)
        self.n_panels = int(self.owners.shape[0])
        if self.breaks.shape != (self.n_panels+1,):
            raise ValueError('compiled coefficient panels do not match their breaks')
        n_cells = len(fields.left)
        self.left = np.ascontiguousarray(fields.left, dtype=np.float64)
        self.width = np.ascontiguousarray(fields.width, dtype=np.float64)
        self.n_local = np.zeros(n_cells, dtype=np.int64)
        self.indices = np.zeros((n_cells, _MAX_LOCAL), dtype=np.int64)
        self.powers = np.zeros((n_cells, 2, 6))
        self.deltas = np.zeros((n_cells, 6))
        self.shapes = np.zeros((n_cells, _MAX_LOCAL, 6))
        self.shape_deltas = np.zeros((n_cells, _MAX_LOCAL, 6))
        basis = fields.basis
        for cell in range(n_cells):
            idx = np.asarray(basis.indices(cell), dtype=np.int64)
            local = int(idx.shape[0])
            if not 1 <= local <= _MAX_LOCAL:
                raise ValueError('compiled coefficient basis exceeds its local width')
            self.n_local[cell] = local
            self.indices[cell, :local] = idx
            powers = np.asarray(fields.powers[cell], dtype=np.float64)
            self.powers[cell, :, :powers.shape[1]] = powers
            deltas = np.asarray(fields.deltas[cell], dtype=np.float64)
            self.deltas[cell, :deltas.shape[0]] = deltas
            shapes = np.asarray(fields.shapes[cell], dtype=np.float64)
            self.shapes[cell, :shapes.shape[1], :] = shapes.T
            slope = np.asarray(fields.shape_deltas[cell], dtype=np.float64)
            self.shape_deltas[cell, :slope.shape[1], :] = slope.T
            budget.check(32)
        self.has_background = background is not None
        self.time = 0.0 if time is None else float(time)
        self.onset_index = np.int64(-1)
        self.quiet = np.int64(0)
        self.knots = np.zeros(2, dtype=np.float64)
        self.selected = np.full((self.n_panels, 2), -1, dtype=np.int64)
        self.use_selected = np.int64(0)
        self.coeff = np.zeros((3, 1, 1))
        self.lengths = np.ones(3, dtype=np.int64)
        self.degree = 1
        if self.has_background:
            profile = background.profile
            self.knots = np.ascontiguousarray(profile._labels, dtype=np.float64)
            self.quiet = np.int64(1 if profile.quiet_extension else 0)
            onset = profile.onset
            if onset is not None and int(onset['index']) >= 1:
                self.onset_index = np.int64(onset['index'])
            spline = profile._spline
            derivatives = profile._derivatives[:2]
            degree = int(spline.c.shape[0])
            pieces = int(spline.c.shape[1])
            self.degree = degree
            self.coeff = np.zeros((3, degree, pieces))
            self.lengths = np.zeros(3, dtype=np.int64)
            self.coeff[0, :, :] = spline.c
            self.lengths[0] = spline.c.shape[0]
            for number, polynomial in enumerate(derivatives, start=1):
                rows = int(polynomial.c.shape[0])
                self.coeff[number, :rows, :] = polynomial.c
                self.lengths[number] = rows
            if panel_indices is not None:
                if len(panel_indices) != self.n_panels:
                    raise ValueError('compiled panel ownership does not match the radial panels')
                for panel, pair in enumerate(panel_indices):
                    pair = np.asarray(pair, dtype=np.int64)
                    if pair.shape != (2,):
                        raise ValueError('compiled panel ownership requires two piece indices')
                    self.selected[panel] = pair
                self.use_selected = np.int64(1)
            budget.check(pieces*degree)
        self.derivative_count = 2*self.count+(1 if self.has_background else 0)
        self.saw_total_field = False
        self._row = np.empty(5+self.derivative_count+1)
        if not np.all(np.isfinite(self.powers)) or not np.all(np.isfinite(self.shapes)):
            raise ValueError('compiled coefficient cache is nonfinite')

    def _charge(self, repeats, budget):
        local = int(np.max(self.n_local))
        budget.check(repeats*accounted_evaluations(local, self.derivative_count, self.degree))

    def _raise(self, status):
        raise ValueError(_MESSAGES.get(int(status), 'compiled coefficient point is unresolved'))

    def _call(self, panel, point, mass, out):
        return int(_point(
            int(panel), float(point), float(mass), self.alpha, self.count, self.cut,
            np.int64(1 if self.has_background else 0), self.time, self.onset_index, self.quiet,
            self.use_selected, self.breaks, self.owners, self.left, self.width, self.n_local,
            self.indices, self.powers, self.deltas, self.shapes, self.shape_deltas, self.knots,
            self.selected, self.coeff, self.lengths, out))

    def evaluate(self, panel, point, mass, budget):
        """One ``coefficients`` row, or None when the reference path still owns it."""
        if isinstance(panel, bool) or not isinstance(panel, (int, np.integer)):
            raise ValueError('fixed radial evaluation escaped its physical cell')
        panel = int(panel)
        if not 0 <= panel < self.n_panels:
            raise ValueError('fixed radial evaluation escaped its physical cell')
        point = float(point)
        mass = float(mass)
        if not math.isfinite(mass):
            raise ValueError('compiled coefficient mass is nonfinite')
        left, right = self.breaks[panel:panel+2]
        if not np.nextafter(left, -np.inf) <= point <= np.nextafter(right, np.inf):
            raise ValueError('fixed radial evaluation escaped its physical cell')
        point = min(float(right), max(float(left), point))
        cell = int(self.owners[panel])
        xi = self.cut if point == self.breaks[-1] else point-cell
        z = float(self.left[cell]+self.width[cell]*xi)
        if self.has_background and z <= SMALL_Z:
            return None
        self._charge(1, budget)
        status = self._call(panel, point, mass, self._row)
        if status != _OK:
            self._raise(status)
        self.saw_total_field = bool(self._row[-1] != 0.0)
        partials = self._row[5:5+self.derivative_count].copy()
        return (float(self._row[0]), float(self._row[1]), float(self._row[2]),
                float(self._row[3]), float(self._row[4]), partials)

    def evaluate_batch(self, panels, points, masses, budget):
        """Per-node coefficient rows. The ordered sum remains outside this call."""
        if np.asarray(panels).dtype.kind not in 'iu':
            raise ValueError('compiled coefficient panels must be integer indices')
        panels = np.ascontiguousarray(panels, dtype=np.int64)
        points = np.ascontiguousarray(points, dtype=np.float64)
        masses = np.ascontiguousarray(masses, dtype=np.float64)
        if panels.ndim != 1 or points.shape != panels.shape or masses.shape != panels.shape:
            raise ValueError('compiled coefficient batch requires three equal vectors')
        if not np.all(np.isfinite(points)) or not np.all(np.isfinite(masses)):
            raise ValueError('compiled coefficient batch is nonfinite')
        n = int(panels.shape[0])
        out = np.empty((n, 5+self.derivative_count+1))
        if n == 0:
            return out[:, :-1]
        clamped_points = np.empty(n, dtype=np.float64)
        for index, (panel, point) in enumerate(zip(panels, points)):
            if not 0 <= int(panel) < self.n_panels:
                raise ValueError('fixed radial evaluation escaped its physical cell')
            left, right = self.breaks[int(panel):int(panel)+2]
            if not np.nextafter(left, -np.inf) <= float(point) <= np.nextafter(right, np.inf):
                raise ValueError('fixed radial evaluation escaped its physical cell')
            clamped = min(float(right), max(float(left), float(point)))
            cell = int(self.owners[int(panel)])
            xi = self.cut if clamped == self.breaks[-1] else clamped-cell
            z = float(self.left[cell]+self.width[cell]*xi)
            if self.has_background and z <= SMALL_Z:
                raise ValueError('compiled coefficient batch excludes the near-centre background path')
            clamped_points[index] = clamped
        self._charge(n, budget)
        status = np.empty(n, dtype=np.int64)
        _points(panels, clamped_points, masses, self.alpha, self.count, self.cut,
                np.int64(1 if self.has_background else 0), self.time, self.onset_index,
                self.quiet, self.use_selected, self.breaks, self.owners, self.left, self.width,
                self.n_local, self.indices, self.powers, self.deltas, self.shapes,
                self.shape_deltas, self.knots, self.selected, self.coeff, self.lengths,
                out, status)
        if np.any(status != _OK):
            self._raise(int(status[int(np.flatnonzero(status != _OK)[0])]))
        return out[:, :5+self.derivative_count]
