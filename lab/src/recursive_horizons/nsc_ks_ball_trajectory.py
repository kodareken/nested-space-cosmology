"""Ball-arithmetic Fourier/polynomial representation of saved KS trajectories.

Every input binary endpoint and dense coefficient is enclosed as supplied.
No mode or history equation is solved here. The represented polynomial is
arbitrary numerical data until its continuous PDE residual is bounded.
"""
from math import comb, factorial

from flint import acb, arb, ctx
import numpy as np

from .nsc_ks_trajectory import TrajectorySegment


def complex_ball(value):
    value = complex(value)
    if not np.isfinite([value.real, value.imag]).all():
        raise ValueError('finite binary complex coefficient required')
    return acb(arb(value.real), arb(value.imag))


def exact_upper(value):
    """Serializable exact dyadic upper bound, without decimal reinterpretation."""
    upper = value.abs_upper() if isinstance(value, acb) else value.upper()
    if not upper.is_finite():
        raise ArithmeticError('finite upper enclosure required')
    mantissa, exponent = upper.man_exp()
    return {'mantissa': str(mantissa), 'exponent': int(exponent)}


def restored_upper(item):
    return arb(int(item['mantissa'])) * arb(2)**int(item['exponent'])


def polynomial_at(coefficients, x):
    result = acb(0)
    for coefficient in reversed(coefficients):
        result = result*x + coefficient
    return result


def affine_restriction(coefficients, left, right):
    """p(left+(right-left)x) in power basis, with interval arithmetic."""
    left, right = arb(left), arb(right)
    if not left >= 0 or not right <= 1 or not right > left:
        raise ValueError('ordered subinterval of [0,1] required')
    width = right-left
    return [sum((coefficients[n]*comb(n, j)*left**(n-j)*width**j
                 for n in range(j, len(coefficients))), acb(0))
            for j in range(len(coefficients))]


def power_to_bernstein(coefficients, degree=None):
    degree = len(coefficients)-1 if degree is None else int(degree)
    if degree < len(coefficients)-1:
        raise ValueError('Bernstein degree cannot truncate the polynomial')
    return [sum((coefficients[j]*arb(comb(k, j))/comb(degree, j)
                 for j in range(min(k, len(coefficients)-1)+1)), acb(0))
            for k in range(degree+1)]


class BallFourierSegment:
    """Exact-input enclosure of weighted X=A+D on the numerical period.

    coefficients[time_power][spin][source][Fourier_mode] uses modes in NumPy
    FFT order. The homogeneous A belongs only to the zero Fourier mode.
    """
    def __init__(self, segment, weights, spatial_count, length, *, bits=120):
        if not isinstance(segment, TrajectorySegment):
            raise TypeError('saved anchored TrajectorySegment required')
        weights = np.asarray(weights, float)
        if weights.ndim != 1 or not len(weights) or not np.isfinite(weights).all() or np.any(weights <= 0):
            raise ValueError('original positive source column weights required')
        if isinstance(spatial_count, bool) or int(spatial_count) != spatial_count or spatial_count < 8:
            raise ValueError('owned spatial mode count required')
        self.spatial_count, self.source_count = int(spatial_count), len(weights)
        self.length = arb(length)
        if not self.length > 0:
            raise ValueError('positive exact numerical period required')
        size = 2*len(weights)
        if len(segment.start) < size*(1+self.spatial_count):
            raise ValueError('trajectory does not contain reference and full envelope difference')
        self.rho_start, self.rho_end = segment.rho_start, segment.rho_end
        self.modes = tuple(j if j < (self.spatial_count+1)//2 else j-self.spatial_count
                           for j in range(self.spatial_count))
        self.bits = int(bits)
        with ctx.workprec(self.bits):
            self.length = arb(length)
            self.waves = tuple(2*arb.pi()*k/self.length for k in self.modes)
            rows = [segment.start, segment.end, *segment.coefficients]
            anchor_modes = []
            for row in rows:
                A, D = row[:size].reshape(2, len(weights)), row[size:size*(1+self.spatial_count)].reshape(2, len(weights), self.spatial_count)
                value = []
                for spin in range(2):
                    sources = []
                    for source, weight in enumerate(weights):
                        modes = [v/self.spatial_count for v in acb.dft([complex_ball(x) for x in D[spin, source]])]
                        modes[0] += complex_ball(A[spin, source])
                        sources.append([v*arb(float(weight)) for v in modes])
                    value.append(sources)
                anchor_modes.append(value)
            self.coefficients = [[[[acb(0) for _ in self.modes] for _ in weights] for _ in range(2)] for _ in range(8)]
            for spin in range(2):
                for source in range(len(weights)):
                    for mode in range(self.spatial_count):
                        c = [acb(0) for _ in range(8)]
                        c[0] = anchor_modes[0][spin][source][mode]
                        c[1] = anchor_modes[1][spin][source][mode]-c[0]
                        for j in range(1, 7):
                            p, q = (j+2)//2, (j+1)//2
                            value = anchor_modes[j+1][spin][source][mode]
                            for k in range(q+1):
                                c[p+k] += (-1)**k*comb(q, k)*value
                        for power in range(8):
                            self.coefficients[power][spin][source][mode] = c[power]

    def value(self, fraction, relative_z, *, time_derivative=0, axial_derivative=0):
        """Weighted envelope derivative; time derivative is in step fraction."""
        if min(time_derivative, axial_derivative) < 0:
            raise ValueError('nonnegative derivative orders required')
        with ctx.workprec(self.bits):
            x, z = arb(fraction), arb(relative_z)
            if not x >= 0 or not x <= 1:
                raise ValueError('fraction enclosure must lie in [0,1]')
            phases = [(acb(0, wave*z)).exp() * acb(0, wave)**axial_derivative for wave in self.waves]
            answer = []
            for spin in range(2):
                row = []
                for source in range(self.source_count):
                    terms = []
                    for mode, phase in enumerate(phases):
                        coeff = [self.coefficients[p][spin][source][mode] * (factorial(p)//factorial(p-time_derivative))
                                 for p in range(time_derivative, 8)]
                        terms.append(polynomial_at(coeff, x)*phase)
                    row.append(sum(terms, acb(0)))
                answer.append(row)
            return answer

    def restricted(self, left, right):
        """Exact polynomial restriction; substep endpoints stay dyadic balls."""
        with ctx.workprec(self.bits):
            left, right = arb(left), arb(right)
            if not left >= 0 or not right <= 1 or not right > left:
                raise ValueError('ordered subinterval of [0,1] required')
            result = object.__new__(type(self))
            result.__dict__ = dict(self.__dict__)
            width = arb(self.rho_end)-arb(self.rho_start)
            result.rho_start = arb(self.rho_start)+width*left
            result.rho_end = arb(self.rho_start)+width*right
            result.coefficients = [[[[acb(0) for _ in self.modes] for _ in range(self.source_count)]
                                     for _ in range(2)] for _ in range(8)]
            for spin in range(2):
                for source in range(self.source_count):
                    for mode in range(self.spatial_count):
                        coeff = affine_restriction([self.coefficients[p][spin][source][mode]
                                                    for p in range(8)], left, right)
                        for p in range(8):
                            result.coefficients[p][spin][source][mode] = coeff[p]
            return result

    def derivative_moment_bounds(self, max_order):
        """Uniform spatial derivative bounds for each time-power coefficient."""
        with ctx.workprec(self.bits):
            result = []
            for coefficient in self.coefficients:
                orders = []
                for order in range(max_order+1):
                    spin_bounds = []
                    for spin in range(2):
                        value = arb(0)
                        for k, wave in enumerate(self.waves):
                            norm = sum((coefficient[spin][s][k].abs_upper()**2 for s in range(self.source_count)), arb(0)).upper().sqrt()
                            value += abs(wave)**order * norm
                        spin_bounds.append(value.upper())
                    orders.append(spin_bounds)
                result.append(orders)
            return result
