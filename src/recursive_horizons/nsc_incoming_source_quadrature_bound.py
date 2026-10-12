"""Directed quadrature of the existing order24 incoming source windows.

No radial propagation or physical-mode bound is performed. The complex
extension conjugates Riccati coefficients, never the complex energy.
"""
from functools import lru_cache
import math
import types

import mpmath as mp
import numpy as np

from .nsc_incoming_source_tail import _mp_bloch
from .nsc_incoming_projector_energy_bound import incoming_riccati_interval_coefficients
from .nsc_incoming_middle_bound import _parameters
from .nsc_incoming_state_moments import incoming_group_factor
from .nsc_incoming_vacuum_tail_bound import _precision, _lo, _hi, _range


KINDS = ('low_vacuum', 'middle_correction')


def pack(value):
    """Lossless binary endpoints; decimal formatting never supplies proof."""
    return [list(endpoint) for endpoint in value._mpi_]


def unpack(value):
    return mp.iv.make_mpf(tuple(tuple(int(x) for x in endpoint) for endpoint in value))


def float_enclosure(value):
    return [float(np.nextafter(float(_lo(value)), -np.inf)),
            float(np.nextafter(float(_hi(value)), np.inf))]


def _conj_coeff(value):
    return mp.iv.mpc(value.real, -value.imag)


def _sqrt(value):
    # ctx_iv.sqrt has no complex implementation; its general power does.
    return value**mp.iv.mpf('0.5')


@lru_cache(maxsize=1)
def interval_bloch():
    """Reuse the generated ad4 code with directed arithmetic and constants.

    All generated floating constants must be exact integers or half-integers.
    Refuse an owner change that could silently round a symbolic constant.
    """
    source = _mp_bloch()
    if any(isinstance(c, float) and (not math.isfinite(c) or 2*c != int(2*c))
           for c in source.__code__.co_consts):
        raise ValueError('generated Bloch constant is not an exact half-integer')
    env = dict(source.__globals__)
    env.update(sin=mp.iv.sin, cos=mp.iv.cos, sqrt=_sqrt, pi=mp.iv.pi,
               array=lambda value: np.asarray(value, dtype=object))
    return types.FunctionType(source.__code__, env)


def _series(energy, coefficients):
    x = 1/(2*energy)
    result = mp.iv.mpc(0)
    for c in reversed(coefficients):
        result = (result+c)*x
    return result


class IntervalWindowKernel:
    """Existing source law on one signed channel; coefficients cached once."""

    def __init__(self, channel, sign=1):
        if sign not in (-1, 1) or channel['index'] not in range(1, 33):
            raise ValueError('one signed non-LLL retained channel required')
        self.channel, self.sign = dict(channel), sign
        self.a, self.r, self.mass, magnitude, self.factor = _parameters(channel)
        self.angular = sign*magnitude
        signs = 1 if channel['angular_eigenvalue'] == 0 else 2
        self.factor /= signs
        stored = mp.mpf(incoming_group_factor(channel, (1,) if signs == 1 else (1, -1)))
        self.factor = _range(min(_lo(self.factor), stored), max(_hi(self.factor), stored))
        # The reused coefficient owner accepts the nonnegative inventory and
        # encloses exact labels plus stored rounded labels. Negative angular
        # coefficients follow its exact c_n(-ell)=(-1)^n conjugate(c_n(ell).
        c = incoming_riccati_interval_coefficients(channel, order=24)
        self.coefficients = c if sign == 1 else [(-1)**(j+1)*_conj_coeff(v) for j, v in enumerate(c)]
        self.bar_coefficients = [_conj_coeff(c) for c in self.coefficients]
        self.order = 24

    def normalization_upper(self, minimum_energy, order):
        """Triangle majorant proves |a² S(E) Sbar(E)|<1 on the disk."""
        e = mp.iv.mpf(minimum_energy)
        if _lo(e) <= 0 or order not in (16, 24):
            raise ValueError('positive energy distance and owned Riccati order required')
        magnitude = sum((abs(c)/(2*e)**(j+1)
                         for j, c in enumerate(self.coefficients[:order])), mp.iv.mpf(0))
        return self.a*self.a*magnitude*magnitude

    def __call__(self, energy, kind):
        if kind not in KINDS:
            raise ValueError('existing low vacuum or middle correction required')
        a, r, mass, ell = self.a, self.r, self.mass, self.angular
        E = energy
        if kind == 'middle_correction':
            S = _series(E, self.coefficients[:16])
            Sb = _series(E, self.bar_coefficients[:16])
            dS = _series(E, self.coefficients[16:])*(1/(2*E))**16
            dSb = _series(E, self.bar_coefficients[16:])*(1/(2*E))**16
            D = 1+a*a*S*Sb
            du = a*a*(S*dSb+Sb*dS+dS*dSb)
            denominator = D*(D+du)
            reS, imS = (S+Sb)/2, (S-Sb)/(2*mp.iv.j)
            reD, imD = (dS+dSb)/2, (dS-dSb)/(2*mp.iv.j)
            delta = [2*a*(D*imD-imS*du)/denominator,
                     -2*a*(D*reD-reS*du)/denominator, -2*du/denominator]
        else:
            S, Sb = _series(E, self.coefficients), _series(E, self.bar_coefficients)
            u = a*a*S*Sb
            y = a*a*(mass*mass+ell*ell/(r*r))/(E*E)
            root = _sqrt(1+y)
            delta = [a*(S-Sb)/(mp.iv.j*(1+u))-mass*a/(E*root),
                     -a*(S+Sb)/(1+u)+ell*a/(r*E*root),
                     y/(root*(1+root))-2*u/(1+u)]
            terms = interval_bloch()(3*mp.iv.pi/4, mass, ell, E)
            for i in range(3):
                delta[i] -= sum((terms[j][i, 0] for j in range(1, 5)), mp.iv.mpf(0))
        parallel = -E/a*delta[2]
        return [-mass*delta[0]+ell/r*delta[1]+parallel,
                parallel, mp.iv.mpc(0), ell/(2*r)*delta[1]]


def legendre_pair(n, x):
    """Exact three-term recurrence for P_n and P_(n-1), interval arithmetic."""
    if n < 1:
        raise ValueError('positive Legendre degree required')
    before, value = mp.iv.mpf(1), x
    for j in range(2, n+1):
        before, value = value, ((2*j-1)*x*value-(j-1)*before)/j
    return value, before


def certify_gauss_brackets(brackets, points):
    """n disjoint sign-changing brackets account for all n exact roots."""
    if points < 2 or len(brackets) != points:
        raise ValueError('one bracket per Legendre root required')
    output, previous = [], mp.mpf(-1)
    for bracket in brackets:
        x = unpack(bracket)
        left, right = _lo(x), _hi(x)
        if not previous < left < right < 1:
            raise ValueError('strictly disjoint ordered root brackets inside (-1,1) required')
        pl, _ = legendre_pair(points, mp.iv.mpf(left))
        pr, _ = legendre_pair(points, mp.iv.mpf(right))
        if not ((_hi(pl) < 0 < _lo(pr)) or (_hi(pr) < 0 < _lo(pl))):
            raise ArithmeticError('interval Legendre endpoint signs do not certify this root')
        pn, previous_p = legendre_pair(points, x)
        derivative = points*(previous_p-x*pn)/(1-x*x)
        if _lo(derivative) <= 0 <= _hi(derivative):
            raise ArithmeticError('root derivative enclosure contains zero')
        weight = 2/((1-x*x)*derivative*derivative)
        if _lo(weight) <= 0:
            raise ArithmeticError('Gauss weight not strictly positive')
        output.append((x, weight)); previous = right
    return output


def gauss_rule(points=48, precision=80):
    """Approximate roots only locate brackets; interval signs certify them."""
    if points not in (8, 48, 64) or precision < 50:
        raise ValueError('owned 8-control/48/64 rule and precision>=50 required')
    with _precision(precision):
        nodes, _ = mp.gauss_quadrature(points, 'legendre')
        width = mp.mpf(10)**(-(precision-20))
        brackets = [pack(_range(x-width, x+width)) for x in nodes]
        certified = certify_gauss_brackets(brackets, points)
        return {'points': points, 'precision': precision, 'roots': brackets,
                'weights': [pack(w) for _, w in certified],
                'proof': 'n disjoint sign-changing Legendre brackets; derivative weight identity; all weights positive'}


def disk_bound(kernel, center, kind, *, radius=8, halfwidth=4, points=48):
    """Rectangle contains the disk; analytic branch/denominator margins proved."""
    if kind not in KINDS or not 0 < halfwidth < radius or center <= 2*radius:
        raise ValueError('disk with Re(E)>|Im(E)| and halfwidth<radius required')
    real = _range(mp.mpf(center)-radius, mp.mpf(center)+radius)
    E = mp.iv.mpc(real, _range(-mp.mpf(radius), mp.mpf(radius)))
    normalization = {str(n): kernel.normalization_upper(center-radius, n)
                     for n in ((16, 24) if kind == 'middle_correction' else (24,))}
    if any(_hi(value) >= 1 for value in normalization.values()):
        raise ArithmeticError('normalization pole exclusion unresolved on selected disk')
    # E² lies in Re>0; adding the positive real gap preserves Re>0.
    # Thus both ad4 gap roots and sqrt(1+a² M²/E²) use their positive real
    # continuation without crossing the negative real axis or zero.
    # The full rectangle proves analyticity, but high powers of wide complex
    # rectangles can lose their nonzero modulus. Cover the boundary circle by
    # small directed angular arcs; Cauchy's estimate only needs this circle.
    maxima = [mp.mpf(0)]*4
    for arc in range(32):
        angle = 2*mp.iv.pi*_range(mp.mpf(arc)/32, mp.mpf(arc+1)/32)
        boundary = mp.iv.mpc(center+radius*mp.iv.cos(angle), radius*mp.iv.sin(angle))
        for j, value in enumerate(kernel(boundary, kind)):
            maxima[j] = max(maxima[j], _hi(abs(value)))
    if not all(mp.isfinite(v) for v in maxima):
        raise ArithmeticError('circle source enclosure is not finite')
    M = [mp.iv.mpf(value) for value in maxima]
    ratio = mp.iv.mpf(halfwidth)/radius
    errors = [4*halfwidth*m*ratio**(2*points)/(1-ratio) for m in M]
    return {'center': center, 'radius': radius, 'halfwidth': halfwidth,
            'energy_squared_real_lower': pack(mp.iv.mpf(center-radius)**2-radius**2),
            'normalization_product_upper': {k: pack(v) for k, v in normalization.items()},
            'kernel_absolute_upper': [pack(v) for v in M],
            'quadrature_error_upper': [pack(v) for v in errors],
            'circle_arc_count': 32,
            'bound': '4*h*M*(h/R)^(2*n)/(1-h/R)'}


def certify_source_window(channel, kind, left, right, *, rule=None, points=48, precision=80, progress=None):
    """Enclose one existing order24 window; no radial/mode solves.

    Return source-integral intervals for the kernel itself, physical signed
    factors included. Physical projector/thermal errors remain separate.
    """
    if kind not in KINDS or channel['index'] not in range(1, 33) or left <= 0 or (right-left) % 8:
        raise ValueError('positive eight-unit cells on one retained source window required')
    with _precision(precision):
        rule = gauss_rule(points, precision) if rule is None else rule
        if rule['points'] != points or rule['precision'] != precision:
            raise ValueError('rule point count and directed precision must match')
        nodes = certify_gauss_brackets(rule['roots'], points)
        signs = (1,) if channel['angular_eigenvalue'] == 0 else (1, -1)
        cells = []
        for sign in signs:
            kernel = IntervalWindowKernel(channel, sign)
            for start in range(int(left), int(right), 8):
                center = start+4
                disk = disk_bound(kernel, center, kind, points=points)
                samples = []
                for x, w in nodes:
                    E = center+4*x
                    values = kernel(E, kind)
                    if any(not (_lo(v.imag) <= 0 <= _hi(v.imag)) for v in values):
                        raise ArithmeticError('real-axis source enclosure lost its real value')
                    samples.append([pack(v.real) for v in values])
                cells.append({'sign': sign, 'left': start, 'right': start+8,
                              'physical_factor': pack(kernel.factor), 'disk': disk, 'samples': samples})
                if progress:
                    progress({'kind': kind, 'sign': sign, 'cell': [start, start+8]})
        return {'group': channel['index'], 'channel': dict(channel), 'kind': kind,
                'interval': [left, right], 'points': points, 'precision': precision,
                'source_order': 24, 'rule': rule, 'cells': cells,
                'parameter_policy': 'exact defining labels and stored binary labels enclosed together; no fitted parameters',
                'complex_policy': 'conjugate coefficients at unchanged E; never conjugate complex E'}


def replay_window(payload):
    """Authenticate exact roots and contract saved directed kernel enclosures."""
    with _precision(payload['precision']):
        rule = payload['rule']; nodes = certify_gauss_brackets(rule['roots'], payload['points'])
        for (_, w), saved in zip(nodes, rule['weights']):
            if pack(w) != saved:
                raise ValueError('saved Gauss weight does not match certified root')
        integral = [mp.iv.mpf(0) for _ in range(4)]
        errors = [mp.iv.mpf(0) for _ in range(4)]
        for cell in payload['cells']:
            if len(cell['samples']) != len(nodes) or cell['right']-cell['left'] != 8:
                raise ValueError('complete owned Gauss cell required')
            factor = unpack(cell['physical_factor'])
            disk = cell['disk']; ratio = mp.iv.mpf(disk['halfwidth'])/disk['radius']
            if _lo(unpack(disk['energy_squared_real_lower'])) <= 0 or any(
                    _hi(unpack(v)) >= 1 for v in disk['normalization_product_upper'].values()):
                raise ArithmeticError('saved analytic-domain margin failed')
            for j in range(4):
                quadrature = sum((4*w*unpack(row[j]) for (_, w), row in zip(nodes, cell['samples'])), mp.iv.mpf(0))
                bound = 4*disk['halfwidth']*unpack(disk['kernel_absolute_upper'][j])*ratio**(2*payload['points'])/(1-ratio)
                if pack(bound) != disk['quadrature_error_upper'][j]:
                    raise ValueError('analytic Gauss error replay differs')
                integral[j] += factor*quadrature
                errors[j] += factor*bound
        result = []
        for q, error in zip(integral, errors):
            upper = _hi(error)
            result.append(q+_range(-upper, upper))
        approximant = [float((_lo(v)+_hi(v))/2) for v in integral]
        midpoint = [mp.iv.mpf(value) for value in approximant]
        radius = [mp.iv.mpf(_hi(abs(v-m))) for v, m in zip(integral, midpoint)]
        a, r, _, _, _ = _parameters(payload['channel'])
        lapse = 4*mp.iv.pi*a*r*r*(radius[0]+errors[0])
        return {'quadrature_sum_interval': [pack(v) for v in integral],
                'source_integral_interval': [pack(v) for v in result],
                'source_integral_float_enclosure': [float_enclosure(v) for v in result],
                'certified_approximant': approximant,
                'quadrature_error_upper': [float_enclosure(v)[1] for v in errors],
                'arithmetic_and_parameter_radius_upper': [float_enclosure(v)[1] for v in radius],
                'lapse_error_upper': float_enclosure(lapse)[1],
                'stationarity_tolerance': 3e-11,
                'numerical_accuracy_pass': _hi(lapse) <= mp.mpf('3e-11'),
                'scope': 'specified order24 approximation only; physical projector/thermal errors separate'}
