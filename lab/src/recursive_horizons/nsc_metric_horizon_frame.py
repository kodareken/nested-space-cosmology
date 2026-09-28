"""Full-profile horizon Frobenius frame with a convergent-series tail bound.

The frame uses the exact NSC compact metric, not its Rindler approximation.
It does not replace the original prepared source or certify its sewing error.
"""
from dataclasses import dataclass

from flint import acb, arb, arb_series, ctx

from .nsc_massive_jost_phase_bound import _series_context


def _profile(q):
    return 3*(arb.pi()-q) + arb(3)/2*(2*q).sin() - q.sin()**2


def horizon_q(hint, *, bits=192):
    """Bracket the unique local zero using directed signs and derivative."""
    with ctx.workprec(bits):
        h = arb(hint)
        if not h.is_finite() or not h > 1:
            raise ValueError('finite positive exterior horizon hint required')
        center = arb.pi()/2 + h.atan()
        left, right = (center-arb('1e-10')).lower(), (center+arb('1e-10')).upper()
        if not _profile(left) > 0 or not _profile(right) < 0:
            raise ValueError('hint does not bracket the fixed metric horizon')
        for _ in range(bits-32):
            mid = (left+right)/2
            value = _profile(mid)
            if value > 0:
                left = mid.lower()
            elif value < 0:
                right = mid.upper()
            else:
                break
        q = left.union(right)
        derivative = -3 + 3*(2*q).cos() - (2*q).sin()
        if not derivative < 0:
            raise ArithmeticError('simple horizon not enclosed')
        return q


@dataclass(frozen=True)
class HorizonFrame:
    q: arb
    slope: arb
    energy: arb
    coefficients: tuple
    analytic_radius: arb
    generator_majorant: arb
    bits: int

    def evaluate(self, distance, *, interior=False):
        """Return component balls, including the infinite-series remainder.

        Rows are spin and columns are the two affine horizon branches. On the
        interior use t=-i sqrt(distance), retaining |u|^lambda separately;
        source thermal factors are not inserted a second time into the frame.
        """
        with ctx.workprec(self.bits):
            distance = arb(distance)
            if not distance.is_finite() or not distance > 0:
                raise ValueError('strictly positive compact distance required')
            radius = distance.sqrt()
            ratio = radius/self.analytic_radius
            if not ratio < 1:
                raise ValueError('evaluation must be strictly inside the analytic disk')
            alpha = 2*self.generator_majorant
            order = len(self.coefficients)-1
            # Cauchy + ||(n/2+lambda-K0)^-1|| <= 2/n majorizes
            # either column by (1-t/R)^(-2M).
            term = arb(1)
            for n in range(1, order+2):
                term *= (alpha+n-1)*ratio/n
            factor = ratio*max(arb(1), ((order+1+alpha)/(order+2)).upper())
            if not factor < 1:
                raise ArithmeticError('geometric remainder ratio is not below one')
            tail = (term/(1-factor)).upper()
            t = acb(0, -radius) if interior else acb(radius)
            answer = [[acb(0), acb(0)] for _ in range(2)]
            for column, sign in enumerate((1, -1)):
                phase = acb(0, sign*self.energy/self.slope*distance.log()).exp()
                for row in range(2):
                    value = acb(0)
                    for c in reversed(self.coefficients):
                        value = value*t + c[row][column]
                    # The l1 tail of each column is <= tail. Component balls
                    # enclose its real and imaginary parts independently.
                    answer[row][column] = phase*(value+acb(arb(0, tail), arb(0, tail)))
            return answer, tail


def metric_horizon_frame(energy, mass, angular, horizon_hint, *, order=16,
                         analytic_radius='0.1', bits=192):
    """Enclose full-metric Frobenius coefficients and a Cauchy majorant.

    All channel parameters are real; Arb intervals are allowed. The bound is
    uniform over them. The stored source surface gravity is deliberately not
    used as a substitute for the exact metric slope.
    """
    if (type(order) is not int or order < 2 or type(bits) is not int or bits < 80):
        raise ValueError('integer order >=2 and precision >=80 required')
    with _series_context(bits, order+2):
        E, m, ell, R = map(arb, (energy, mass, angular, analytic_radius))
        if not all(v.is_finite() for v in (E, m, ell, R)) or not m >= 0 or not R > 0:
            raise ValueError('finite real labels, nonnegative mass and positive disk required')
        q = horizon_q(horizon_hint, bits=bits)
        u = arb_series([q, 1], prec=order+3)
        W = _profile(u)
        # W(q_h)=0 exactly. Derivatives are enclosed on the root bracket.
        H = arb_series([-W[k+1] for k in range(order+2)], prec=order+2)
        h0 = H[0]
        U = R*R
        # |W^(n)| <= 2*2^n for n>=2, valid on the real root;
        # entire Taylor series gives this bound on the complex u disk.
        delta = 2*((2*U).exp()-1-2*U)/U
        lower_h = h0-delta
        lower_sin = q.sin()-q.sin().abs_upper()*(U.cosh()-1)-q.cos().abs_upper()*U.sinh()
        if not lower_h > 0 or not lower_sin > 0:
            raise ArithmeticError('metric or areal-radius analytic disk not enclosed')
        M = (E.abs_upper()*delta/(h0*lower_h)
             + R*(m.abs_upper()/lower_sin+ell.abs_upper())/lower_h.sqrt()).upper()
        inv = H.inv()
        off_mass = -m/u.sin()/H.sqrt()
        off_ell = ell/H.sqrt()
        K = []
        for n in range(order+1):
            if n % 2 == 0:
                z = acb(0, E*inv[n//2])
                K.append(((z, acb(0)), (acb(0), -z)))
            else:
                k = n//2
                K.append(((acb(0), acb(off_mass[k], -off_ell[k])),
                          (acb(off_mass[k], off_ell[k]), acb(0))))
        coefficients = [((acb(1), acb(0)), (acb(0), acb(1)))]
        for n in range(1, order+1):
            entry = [[acb(0), acb(0)] for _ in range(2)]
            for column, sign in enumerate((1, -1)):
                lam = acb(0, sign*E/h0)
                for row in range(2):
                    rhs = sum((K[j][row][spin]*coefficients[n-j][spin][column]
                               for j in range(1, n+1) for spin in range(2)), acb(0))
                    entry[row][column] = rhs/(arb(n)/2+lam-K[0][row][row])
            coefficients.append(tuple(tuple(row) for row in entry))
        return HorizonFrame(q, h0, E, tuple(coefficients), R, M, bits)


def reflection_from_phase(frame, radial_coordinate, phase):
    """Map an enclosed decaying exterior phase to affine horizon reflection.

    The exterior vector is (1, exp(i*phase)); its arbitrary common amplitude
    cancels in c_minus/c_plus. The exact metric frame includes the entire
    collar from this point to the horizon. No finite-collar Rindler matching
    or normalization of the returned reflection is performed.
    """
    if not isinstance(frame, HorizonFrame):
        raise TypeError('a validated full-profile frame is required')
    with ctx.workprec(frame.bits):
        rho, theta = arb(radial_coordinate), arb(phase)
        if not rho.is_finite() or not theta.is_finite():
            raise ValueError('finite radial and phase enclosures required')
        distance = arb.pi()/2+rho.atan()-frame.q
        F, tail = frame.evaluate(distance)
        z = acb(0,theta).exp()
        # Cofactors: the common determinant cancels from the ratio.
        denominator = F[1][1]-F[0][1]*z
        numerator = -F[1][0]+F[0][0]*z
        if denominator.contains(0):
            raise ArithmeticError('outgoing horizon coefficient may vanish')
        return numerator/denominator, distance, tail
