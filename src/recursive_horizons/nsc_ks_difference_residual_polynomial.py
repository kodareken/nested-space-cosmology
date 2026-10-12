"""Polynomial residual of the owned reference/difference split.

Production interpolates D=X-A with A'=G_ref A and D'=L_g D+M A. A is
z-independent. The identity R_D=R_X-R_A is formed on Bernstein
coefficients. Already bounded norms of R_X and R_A are not subtracted.
"""
from flint import acb, arb, ctx
import numpy as np

from .nsc_ks_ball_trajectory import BallFourierSegment
from .nsc_ks_fourier_residual_bound import operator_remainder_bounds, polynomial_residual_bounds
from .nsc_ks_residual_polynomial import ResidualPolynomial, row_norm
from .nsc_ks_trajectory import TrajectorySegment


def monomial_keys(order=4):
    """w^p U^q keys of the owned reciprocal-radius expansion, p+q<=order."""
    if isinstance(order, bool) or not isinstance(order, int) or order < 1:
        raise ValueError('positive reciprocal-radius order required')
    return tuple((n - q, q) for n in range(1, order + 1) for q in range(n + 1))


def _state_mask(row, source_count, spatial_count, *, keep_reference, keep_difference):
    row = np.array(row, complex, copy=True)
    size = 2 * int(source_count)
    end = size * (1 + int(spatial_count))
    if row.ndim != 1 or len(row) < end:
        raise ValueError('captured state does not contain reference and envelope difference')
    if not keep_reference:
        row[:size] = 0
    if not keep_difference:
        row[size:end] = 0
    return row


def split_trajectory_segment(segment, source_count, spatial_count):
    """Return the captured X segment together with A-only and D-only copies."""
    if not isinstance(segment, TrajectorySegment):
        raise TypeError('saved anchored TrajectorySegment required')
    nsrc, nz = int(source_count), int(spatial_count)

    def build(keep_reference, keep_difference):
        return TrajectorySegment(
            segment.rho_start, segment.rho_end,
            _state_mask(segment.start, nsrc, nz, keep_reference=keep_reference,
                        keep_difference=keep_difference),
            _state_mask(segment.end, nsrc, nz, keep_reference=keep_reference,
                        keep_difference=keep_difference),
            np.stack([_state_mask(row, nsrc, nz, keep_reference=keep_reference,
                                  keep_difference=keep_difference)
                      for row in segment.coefficients]))

    return segment, build(True, False), build(False, True)


def ball_split(segment, weights, spatial_count, length, *, bits=120):
    """Weighted Fourier balls of X, A and D from one captured reconstruction."""
    weights = np.asarray(weights, float)
    x_seg, a_seg, d_seg = split_trajectory_segment(segment, len(weights), spatial_count)
    return tuple(BallFourierSegment(item, weights, spatial_count, length, bits=bits)
                 for item in (x_seg, a_seg, d_seg))


def _matching_zero_potential(potential):
    return {key: [arb(0)] * len(coeff) for key, coeff in potential.items()}


def _sub_cube(left, right):
    if len(left) != len(right):
        raise ValueError('Bernstein cubes of equal degree required')
    degree, sources, modes = len(left) - 1, len(left[0][0]), len(left[0][0][0])
    result = [[[[acb(0) for _ in range(modes)] for _ in range(sources)] for _ in range(2)]
              for _ in range(degree + 1)]
    for p in range(degree + 1):
        for spin in range(2):
            for source in range(sources):
                for mode in range(modes):
                    result[p][spin][source][mode] = left[p][spin][source][mode] - right[p][spin][source][mode]
    return result


class DifferenceResidualPolynomial:
    """S_D=D_xi-h L_g D-h M A on the same Bernstein lattice as ResidualPolynomial.

    G_ref is the background operator with empty physical potential. Matching
    zero potential keys force the same Bernstein degree without putting M A
    into the reference residual.
    """

    def __init__(self, field, reference_field, inv_a2, inv_a, inv_ar, potential,
                 mass, angular, energies):
        if not isinstance(field, BallFourierSegment) or not isinstance(reference_field, BallFourierSegment):
            raise TypeError('combined X and homogeneous A BallFourierSegment fields required')
        if (field.source_count != reference_field.source_count
                or field.spatial_count != reference_field.spatial_count
                or field.bits != reference_field.bits):
            raise ValueError('reference and combined fields must share source, mode and precision layout')
        self.field = field
        self.reference_field = reference_field
        self.bits = field.bits
        self.potential_keys = tuple(potential)
        with ctx.workprec(self.bits):
            self.combined = ResidualPolynomial(
                field, inv_a2, inv_a, inv_ar, potential, mass, angular, energies)
            self.reference = ResidualPolynomial(
                reference_field, inv_a2, inv_a, inv_ar, _matching_zero_potential(potential),
                mass, angular, energies)
            if self.combined.degree != self.reference.degree:
                raise ValueError('reference and combined residual Bernstein degrees must match')
            self.degree = self.combined.degree
            self.free = _sub_cube(self.combined.free, self.reference.free)
            self.terms = self.combined.terms
            self._moment_cache = {}

    def _moment(self, spectra, order):
        key = (id(spectra), order)
        if key not in self._moment_cache:
            self._moment_cache[key] = tuple(
                tuple(sum((abs(self.field.waves[k]) ** order * row_norm(
                    [coefficient[spin][s][k] for s in range(self.field.source_count)])
                    for k in range(self.field.spatial_count)), arb(0)).upper()
                      for spin in range(2)) for coefficient in spectra)
        return [list(row) for row in self._moment_cache[key]]

    def free_contains_zero(self):
        """True when the background residual of D vanishes identically."""
        with ctx.workprec(self.bits):
            for coefficient in self.free:
                for spin in range(2):
                    for row in coefficient[spin]:
                        for value in row:
                            if not value.contains(0):
                                return False
            return True


def difference_operator_remainder_bounds(field_d, field_a, time_errors, profile_sup,
                                         radius_tail, mass, angular, energies):
    """Time-Taylor and radius remainders of L_g on D and of M on A.

    Background G_ref remainders belong to R_A and are omitted here. The two
    remainder balls are added; they act on different fields, so no cancelled
    A background term is recovered by subtracting residual norms.
    """
    with ctx.workprec(max(field_d.bits, field_a.bits)):
        background = {name: time_errors[name] for name in ('inv_a2', 'inv_a', 'inv_ar')}
        potential = {key: value for key, value in time_errors.items() if isinstance(key, tuple)}
        rem_d = operator_remainder_bounds(
            field_d, {**background, **potential}, profile_sup, radius_tail, mass, angular, energies)
        rem_a = operator_remainder_bounds(
            field_a, {**{name: arb(0) for name in background}, **potential},
            profile_sup, radius_tail, mass, angular, energies)
        return tuple((left + right).upper() for left, right in zip(rem_d, rem_a))


def difference_polynomial_bounds(residual, profiles):
    """Row-l2 residual bounds of S_D, including omitted profile bands."""
    if not isinstance(residual, DifferenceResidualPolynomial):
        raise TypeError('DifferenceResidualPolynomial required')
    return polynomial_residual_bounds(residual, profiles)
