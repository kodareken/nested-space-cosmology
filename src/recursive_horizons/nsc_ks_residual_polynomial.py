"""Fourier/Bernstein polynomial part of the owned KS equation residual.

Scalar time polynomials must come from the geometry enclosure owner. Their
remainders and the reciprocal-radius remainder are additional errors; this
algebra alone cannot certify the full continuous equation.
"""
from math import comb, factorial

from flint import arb, acb, ctx

from .nsc_ks_ball_trajectory import power_to_bernstein


def convolution(left, right):
    result = [acb(0) for _ in range(len(left)+len(right)-1)]
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            result[i+j] += x*y
    return result


def row_norm(row):
    return sum((value.abs_upper()**2 for value in row), arb(0)).upper().sqrt().upper()


def _cube(degree, sources, modes):
    return [[[[acb(0) for _ in range(modes)] for _ in range(sources)] for _ in range(2)] for _ in range(degree+1)]


class ResidualPolynomial:
    """S=X_xi-h L_g X with polynomial coefficient approximations.

    inv_a2, inv_a, inv_ar and potential[(p,q)] are power coefficients in
    the same step fraction xi as the captured field. Each potential term is
    the real scalar multiplying i*sigma2*w(z)^p*U(z)^q.
    """
    def __init__(self, field, inv_a2, inv_a, inv_ar, potential, mass, angular, energies):
        if len(energies) != field.source_count:
            raise ValueError('one unchanged energy label per source column required')
        self.field, self.potential_keys = field, tuple(potential)
        self.degree = 7+max(len(inv_a2), len(inv_a), len(inv_ar), *(len(v) for v in potential.values()), 1)-1
        self.bits = field.bits
        self._moment_cache = {}
        with ctx.workprec(self.bits):
            h = arb(field.rho_end)-arb(field.rho_start)
            m, ell = arb(float(mass)), arb(float(angular))
            self.free = _cube(self.degree, field.source_count, field.spatial_count)
            self.terms = {key: _cube(self.degree, field.source_count, field.spatial_count) for key in potential}
            for source, energy in enumerate(energies):
                E = arb(float(energy))
                for mode, wave in enumerate(field.waves):
                    X = [[field.coefficients[p][spin][source][mode] for p in range(8)] for spin in range(2)]
                    for spin, sign in ((0, 1), (1, -1)):
                        derivative = [j*X[spin][j] for j in range(1, 8)]
                        own = convolution(inv_a2, X[spin])
                        mass_mix = convolution(inv_a, X[1-spin])
                        angular_mix = convolution(inv_ar, X[1-spin])
                        value = [acb(0) for _ in range(self.degree+1)]
                        for p, v in enumerate(derivative):
                            value[p] += v
                        for p in range(self.degree+1):
                            diagonal = acb(0, sign*(wave-E))*own[p] if p < len(own) else acb(0)
                            off = sign*ell*angular_mix[p] if p < len(angular_mix) else acb(0)
                            off -= acb(0, m)*mass_mix[p] if p < len(mass_mix) else acb(0)
                            value[p] -= h*(diagonal+off)
                        bernstein = power_to_bernstein(value, self.degree)
                        for p, v in enumerate(bernstein):
                            self.free[p][spin][source][mode] = v
                        for key, scalar in potential.items():
                            term = [-h*sign*v for v in convolution(scalar, X[1-spin])]
                            coefficients = power_to_bernstein(term, self.degree)
                            for p, v in enumerate(coefficients):
                                self.terms[key][p][spin][source][mode] = v

    def _values(self, spectra, relative_z, max_derivative):
        z = arb(relative_z)
        phase = [acb(0, wave*z).exp() for wave in self.field.waves]
        answer = []
        for derivative in range(max_derivative+1):
            factor = [p*acb(0, k)**derivative for p, k in zip(phase, self.field.waves)]
            answer.append([[[sum((row[k]*factor[k] for k in range(self.field.spatial_count)), acb(0))
                             for row in coefficient[spin]] for spin in range(2)] for coefficient in spectra])
        return answer

    def point_coefficients(self, relative_z, profile_jets, max_derivative):
        """Spatial derivatives of each Bernstein coefficient at one z ball.

        profile_jets[key][j] encloses the ordinary j-th derivative of w^p U^q.
        """
        if set(profile_jets) != set(self.terms):
            raise ValueError('all spatial potential profiles required exactly once')
        with ctx.workprec(self.bits):
            values = self._values(self.free, relative_z, max_derivative)
            for key, spectra in self.terms.items():
                profile = profile_jets[key]
                if len(profile) <= max_derivative:
                    raise ValueError('insufficient spatial profile derivatives')
                field = self._values(spectra, relative_z, max_derivative)
                for j in range(max_derivative+1):
                    for q in range(self.degree+1):
                        for spin in range(2):
                            for source in range(self.field.source_count):
                                values[j][q][spin][source] += sum((comb(j, k)*profile[k]*field[j-k][q][spin][source]
                                                                for k in range(j+1)), acb(0))
            return values

    def _moment(self, spectra, order):
        key = (id(spectra), order)
        if key not in self._moment_cache:
            self._moment_cache[key] = tuple(tuple(sum((abs(self.field.waves[k])**order*row_norm([coefficient[spin][s][k]
                      for s in range(self.field.source_count)]) for k in range(self.field.spatial_count)), arb(0)).upper()
                 for spin in range(2)) for coefficient in spectra)
        return [list(row) for row in self._moment_cache[key]]

    def spatial_cell_bounds(self, relative_center, halfwidth, profile_jets_center,
                            profile_derivative_bounds, *, axial_orders=(0, 1), taylor_order=4):
        """Continuous xi/z bound for the polynomial residual part on a cell.

        Bernstein convexity bounds all xi in [0,1]; a spatial Taylor bound
        includes a derivative remainder on the full z cell.
        """
        highest = max(axial_orders)+taylor_order+1
        with ctx.workprec(self.bits):
            half = arb(halfwidth)
            if not half >= 0:
                raise ValueError('nonnegative spatial halfwidth required')
            center = self.point_coefficients(relative_center, profile_jets_center, highest-1)
            outputs = []
            for target in axial_orders:
                derivative = target+taylor_order+1
                remainder = self._moment(self.free, derivative)
                for key, spectra in self.terms.items():
                    if key not in profile_derivative_bounds or len(profile_derivative_bounds[key]) <= derivative:
                        raise ValueError('full-cell spatial profile derivative bounds required')
                    for k in range(derivative+1):
                        bound = arb(profile_derivative_bounds[key][k])
                        if not bound >= 0:
                            raise ValueError('nonnegative spatial derivative bounds required')
                        moments = self._moment(spectra, derivative-k)
                        for q in range(self.degree+1):
                            for spin in range(2):
                                remainder[q][spin] += comb(derivative, k)*bound*moments[q][spin]
                maximum = arb(0)
                for q in range(self.degree+1):
                    for spin in range(2):
                        bound = sum((half**j/factorial(j)*row_norm(center[target+j][q][spin])
                                     for j in range(taylor_order+1)), arb(0))
                        bound += half**(taylor_order+1)/factorial(taylor_order+1)*remainder[q][spin]
                        maximum = arb.max(maximum, bound.upper())
                outputs.append(maximum.upper())
            return outputs
