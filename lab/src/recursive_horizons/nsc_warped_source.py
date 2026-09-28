"""Finite proper-time Dirac modulus on the declared conformal interval.

Two complementary five-dimensional Dirac copies; no link mass or torsion.
The compact weight composes with the existing four-dimensional D4 spectrum.
It is not the complete quantum functional or a Lorentzian vacuum choice.
"""
from math import pi, sqrt

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad_vec
from scipy.linalg import block_diag, cholesky, eigh, solve_triangular
from scipy.special import erfc, exp1


class WarpedCompactWeight:
    """Galerkin quadratic form for D_E^dagger D_E, with its actual measure.

    x=(Y+ell/2)/ell in [0,1]; sigma=-a*(2*x-1)^2. The path parameter s
    multiplies sigma. In the canonical free field chi, the eigenproblem is
    integral exp(-s*sigma)|D0 chi|² = epsilon integral exp(s*sigma)|chi|².
    Allowed values use cosine modes; the complementary component is Dirichlet.
    No derivative condition is imposed as an additional first-order wall.
    The phase basis chi=(cosine*f, i*sine*g) makes the complex Hermitian
    quadratic form real. Its weighted cross term is retained, not discarded
    by restricting the physical spinor to real components.
    """

    def __init__(self, cutoff=2., interval=2., modes=12, points=96,
                 amplitude=18/1015, path=1.):
        if not all(np.isfinite(x) for x in (cutoff, interval, amplitude, path)):
            raise ValueError("finite scales and warp required")
        if cutoff <= 0 or interval <= 0 or amplitude < 0 or not 0 <= path <= 1:
            raise ValueError("positive scales, nonnegative warp, and path in [0,1] required")
        if not isinstance(modes, int) or isinstance(modes, bool) or modes < 2:
            raise ValueError("at least two compact modes required")
        if not isinstance(points, int) or points < 2*modes+4:
            raise ValueError("resolve the compact overlap quadrature")
        self.cutoff, self.interval = cutoff, interval
        self.modes, self.points = modes, points
        self.amplitude, self.path = amplitude, path
        x, w = leggauss(points)
        sigma = -amplitude*x*x
        x, w = (x+1)/2, w/2
        n = np.arange(1, modes+1)
        k = n*pi/interval
        c = np.column_stack((np.ones(points), sqrt(2)*np.cos(pi*x[:, None]*n)))
        sn = sqrt(2)*np.sin(pi*x[:, None]*n)
        dc = np.column_stack((np.zeros(points), -sn*k))
        ds = c[:, 1:]*k

        def gram(left, right, weight):
            return left.T @ (weight[:, None]*right)

        def kinetic(weight):
            k0 = block_diag(gram(dc, dc, weight), gram(ds, ds, weight))
            k2 = block_diag(gram(c, c, weight), gram(sn, sn, weight))
            # |D0(f,i*g)|²=(lambda*g-f')²+(lambda*f-g')².
            # The cross derivative integrates to zero only for constant weight.
            cross = -gram(dc, sn, weight)-gram(c, ds, weight)
            k1 = np.zeros_like(k0)
            k1[:modes+1, modes+1:] = cross
            k1[modes+1:, :modes+1] = cross.T
            return k0, k1, k2

        plus, minus = w*np.exp(path*sigma), w*np.exp(-path*sigma)
        self.mass = block_diag(gram(c, c, plus), gram(sn, sn, plus))
        self.mass_s = block_diag(gram(c, c, sigma*plus), gram(sn, sn, sigma*plus))
        self.k0, self.k1, self.k2 = kinetic(minus)
        self.ks0, self.ks1, self.ks2 = kinetic(-sigma*minus)
        chol = cholesky(self.mass, lower=True)
        self.inverse_chol = solve_triangular(chol, np.eye(2*modes+1), lower=True)

    def response(self, fourD_squared, derivatives=True):
        """h(y)=Tr_compact E1(epsilon(y)/Lambda²), y>0; no extra copy factor.

        The full two-copy functional is (1/2) Tr_4 h(D4²). A discrete y=0
        eigenvalue needs the separate zero-mode/IR prescription, not clipping.
        """
        y = float(fourD_squared)
        if not np.isfinite(y) or y <= 0:
            raise ValueError("positive four-dimensional squared eigenvalue required")
        lam = sqrt(y)
        kinetic = self.k0+lam*self.k1+y*self.k2
        transform = self.inverse_chol
        matrix = transform @ kinetic @ transform.T
        values, u = eigh((matrix+matrix.T)/2, check_finite=False)
        if values[0] <= 0:
            raise ArithmeticError("unresolved positive singular value; refine the representation")
        result = {"weight": float(np.sum(exp1(values/self.cutoff**2)))}
        if not derivatives:
            return result
        v = transform.T @ u
        diagonal = lambda matrix: np.einsum("ij,ij->j", v, matrix @ v)
        exponential = np.exp(-values/self.cutoff**2)
        derivative = -exponential/values
        eigen_y = diagonal(self.k2+self.k1/(2*lam))
        eigen_ell = diagonal(-(2*self.k0+lam*self.k1)/self.interval)
        eigen_s = diagonal(self.ks0+lam*self.ks1+y*self.ks2)-values*diagonal(self.mass_s)
        result.update({
            "d_squared": float(derivative @ eigen_y),
            "d_interval": float(derivative @ eigen_ell),
            "d_path": float(derivative @ eigen_s),
            "d_cutoff": float(2*np.sum(exponential)/self.cutoff),
            "dual_heat_path_derivative": float(2*exponential @ diagonal(self.mass_s)),
        })
        return result


def angular_modes(flux, radius, maximum):
    """BKW integer flux convention: nonzero degeneracies include both signs."""
    if not isinstance(flux, int) or isinstance(flux, bool):
        raise ValueError("integer magnetic flux required")
    if not np.isfinite(radius) or radius <= 0 or maximum < 0:
        raise ValueError("positive radius and nonnegative angular maximum required")
    q = abs(flux)
    return [(n*(n+q)/radius**2, q if n == 0 else 2*(q+2*n))
            for n in range(maximum+1) if n or q]


def integrated_exp1(mass_squared, cutoff):
    """Integral_0^infinity dz E1((z+mass_squared)/cutoff²)."""
    y = np.asarray(mass_squared, dtype=float)
    nonzero = y > 0
    result = np.asarray(cutoff**2*np.exp(-y/cutoff**2))
    result[nonzero] -= y[nonzero]*exp1(y[nonzero]/cutoff**2)
    return result


def omitted_tail_bounds(flux, radius, interval, cutoff, angular_max,
                        compact_max, momentum_extent):
    """Flat comparison bounds for omitted angular/KK/momentum traces.

    The nonpositive warp raises singular values, so these also bound omitted
    continuum sectors there. They do NOT bound the Galerkin errors of retained
    singular values. Bounds are evaluated in ordinary floating point.
    """
    q = abs(flux)
    if (2*angular_max+q)**2 < 2*(cutoff*radius)**2 or compact_max < 1:
        raise ValueError("tail must start in the decreasing part of the spectrum")
    sectors = angular_modes(flux, radius, angular_max)
    angular = sum(d*np.exp(-m/cutoff**2) for m, d in sectors)
    angular_tail = 2*(cutoff*radius)**2*np.exp(-angular_max*(angular_max+q)/(cutoff*radius)**2)
    compact = 1+2*np.sum(np.exp(-(pi*np.arange(1, compact_max+1)/(interval*cutoff))**2))
    compact_tail = interval*cutoff/sqrt(pi)*erfc(pi*compact_max/(interval*cutoff))
    total_a, total_y = angular+angular_tail, compact+compact_tail
    factor = cutoff**2/(4*pi)
    return {
        "angular_potential_bound": float(factor*angular_tail*total_y),
        "omitted_compact_level_potential_bound": float(factor*total_a*compact_tail),
        "momentum_potential_bound": float(factor*np.exp(-momentum_extent**2/cutoff**2)*total_a*total_y),
    }


def product_potential(flux, radius, interval, cutoff, angular_max, compact_max):
    """Analytic unwarped momentum integral, at stated angular/KK truncations."""
    sectors = angular_modes(flux, radius, angular_max)
    a2, d = np.array(sectors).T
    b2 = (pi*np.arange(compact_max+1)/interval)**2
    k = np.ones(compact_max+1)*2
    k[0] = 1
    mass = a2[:, None]+b2[None, :]
    return float(np.sum(d[:, None]*k*integrated_exp1(mass, cutoff))/(4*pi))


def source_potential(weight, flux=1, radius=1., angular_max=10,
                     momentum_extent=12., quadrature_atol=2e-9):
    """Euclidean density per reference R2 area on R2 x S2(r) x I, with warp.

    Returns the same functional's radius, interval, cutoff and warp-path
    variations. These are homogeneous projections, not the full stress on
    a varying throat. Angular/KK truncation errors are assessed separately.
    """
    if momentum_extent <= 0 or quadrature_atol <= 0:
        raise ValueError("positive momentum extent and quadrature tolerance required")
    sectors = angular_modes(flux, radius, angular_max)
    fields = ("weight", "d_interval", "d_path", "d_cutoff", "dual_heat_path_derivative")

    def integrand(momentum):
        result = np.zeros(len(fields))
        for mass, degeneracy in sectors:
            response = weight.response(momentum**2+mass)
            result += degeneracy*np.array([response[name] for name in fields])
        return momentum*result/(2*pi)

    values, error = quad_vec(integrand, 0., momentum_extent,
                            epsabs=quadrature_atol, epsrel=quadrature_atol, norm="max")
    dr = 0.
    for mass, degeneracy in sectors:
        if mass:
            # Fundamental theorem after integrating the R2 momentum:
            # lower spectral endpoint mass(r), not a numerical metric difference.
            start = weight.response(mass, derivatives=False)["weight"]
            end = weight.response(mass+momentum_extent**2, derivatives=False)["weight"]
            dr += degeneracy*mass*(start-end)/(2*pi*radius)
    result = dict(zip(("potential", "d_interval", "d_path", "d_cutoff",
                       "dual_heat_path_derivative"), map(float, values)))
    result["d_radius"] = float(dr)
    result["quadrature_error_estimate"] = float(error)
    result["scaling_residual"] = float(radius*dr+weight.interval*values[1]
                                       -weight.cutoff*values[3]+2*values[0])
    result["path_variation_residual"] = float(values[2]-values[4])
    return result
