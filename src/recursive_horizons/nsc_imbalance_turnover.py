"""Cesaro circulation of one incoherent regional imbalance.

The reusable kernel projects a state onto polynomials in a Hermitian
generator with the Frobenius product. On the frozen six-mode window that
projection is the time average of unitary evolution when the generator is
nondegenerate. This does not close the local incoming gate.
"""
from __future__ import annotations

import mpmath as mp
import numpy as np
import sympy as sp

from .nsc_finite_turnover import (
    _as_numpy,
    _current,
    _filtered,
    _frobenius,
    _onsite_projectors,
    _pair,
    _rational,
)
from .nsc_nested_qualities import finite_window

SCHEMA = "NSC-IMBALANCE-TURNOVER-v1"
PASS = "PASS_IMBALANCE_MEAN_CIRCULATION"
ZERO = "ZERO_MEAN_CIRCULATION"
MEAN_ONLY = "MEAN_CIRCULATION_RECORDED_FINITE_T_OPEN"
FAIL = "FAIL_IMBALANCE_TURNOVER"

_SQRT29 = sp.sqrt(29)


def _powers(operator, highest):
    """Return [I, J, ..., J**highest] by exact multiplication.

    Each product is reduced to Gaussian-rational `a+b*I` before the next
    multiply. Leaving that product unexpanded makes `J**10` blow up.
    """
    powers = [sp.eye(operator.rows)]
    current = powers[0]
    for _ in range(highest):
        current = _cancel_matrix(current * operator)
        powers.append(current)
    return powers


def _real_rational(value, label):
    value = sp.expand(value)
    real, imag = sp.together(sp.re(value)), sp.together(sp.im(value))
    if imag != 0:
        raise ValueError(f"{label} is not real: {value}")
    if not real.is_rational:
        raise ValueError(f"{label} is not rational: {real}")
    return sp.Rational(real)


def _power_traces(powers):
    return [_real_rational(sp.trace(power), f"trace power {index}") for index, power in enumerate(powers)]


def _gram_matrix(traces, degree):
    return sp.Matrix(degree, degree, lambda row, col: traces[row + col])


def _moment_column(powers, state, degree):
    """Real moments in `Q(sqrt(29))`. Power traces that build the Gram stay rational."""
    column = sp.zeros(degree, 1)
    for index in range(degree):
        column[index, 0] = _real_q_sqrt29(sp.trace(powers[index] * state), f"moment {index}")
    return column


def _cancel_matrix(matrix):
    cancelled = sp.zeros(matrix.rows, matrix.cols)
    for row in range(matrix.rows):
        for col in range(matrix.cols):
            entry = sp.expand(matrix[row, col])
            real, imag = sp.together(sp.re(entry)), sp.together(sp.im(entry))
            cancelled[row, col] = real + sp.I * imag
    return cancelled


def commutant_projection(operator, state, powers=None):
    """Project `state` onto span{I, J, ..., J**(n-1)}.

    `G_ab = Tr(J**(a+b))` and `b_a = Tr(J**a state)`. A nonzero Gram
    determinant means those powers are a basis. Gram entries stay rational.
    Moments and coefficients may be exact reals `a+b*sqrt(29)`, which is
    required when `state` is an onsite mode projector. For nondegenerate
    Hermitian `J` the powers are the commutant, so the solution is the Cesaro
    mean of `exp(-i J t) state exp(i J t)`. Pass `powers` to reuse one tower.
    """
    degree = operator.rows
    if operator.cols != degree:
        raise ValueError("commutant projection expects a square generator")
    if powers is None:
        powers = _powers(operator, 2 * degree - 2)
    if len(powers) < 2 * degree - 1:
        raise ValueError("power tower stops before J**(2n-2)")
    traces = _power_traces(powers)
    gram = _gram_matrix(traces, degree)
    determinant = gram.det()
    moments = _moment_column(powers, state, degree)
    if determinant == 0:
        return {
            "mean": None,
            "alpha": None,
            "gram_determinant": determinant,
            "traces": traces,
            "powers": powers,
            "moments_matched": False,
            "commutes": False,
        }
    alpha = gram.solve(moments)
    mean = sp.zeros(degree)
    for index in range(degree):
        coefficient = _real_q_sqrt29(alpha[index, 0], f"projection coefficient {index}")
        alpha[index, 0] = coefficient
        if coefficient != 0:
            mean += coefficient * powers[index]
    mean = _cancel_matrix(mean)
    moment_error = _moment_column(powers, mean - state, degree)
    commutes = _matrix_is_zero(operator * mean - mean * operator)
    return {
        "mean": mean,
        "alpha": alpha,
        "gram_determinant": sp.together(determinant),
        "traces": traces,
        "powers": powers,
        "moments_matched": all(moment_error[index, 0] == 0 for index in range(degree)),
        "commutes": commutes,
    }


def _matrix_is_zero(matrix):
    for row in range(matrix.rows):
        for col in range(matrix.cols):
            entry = matrix[row, col]
            if entry != 0 and sp.expand(entry) != 0:
                return False
    return True


def _regional_shift(depth=3, block=2):
    size = depth * block
    shift = sp.zeros(size)
    shift[:block, :block] = sp.eye(block)
    shift[-block:, -block:] = -sp.eye(block)
    return shift


def regional_preparation(delta, depth=3, block=2):
    """`I/2 + delta (P_region0 - P_region_last)`, with no off-block coherence."""
    size = depth * block
    return sp.eye(size) / 2 + sp.sympify(delta) * _regional_shift(depth, block)


def _embed(block_matrix, region, depth=3, block=2):
    embedded = sp.zeros(depth * block)
    start = region * block
    embedded[start:start + block, start:start + block] = block_matrix
    return embedded


def _as_q_sqrt29(value):
    value = sp.expand(value)
    coefficient = sp.together(value.coeff(_SQRT29))
    remainder = sp.together(sp.expand(value - coefficient * _SQRT29))
    if remainder.has(_SQRT29) or coefficient.has(_SQRT29):
        return None
    if not remainder.is_rational or not coefficient.is_rational:
        return None
    return sp.Rational(remainder), sp.Rational(coefficient)


def _real_q_sqrt29(value, label):
    """Return an exact real `a+b*sqrt(29)`, or a rational when `b` is zero."""
    value = sp.expand(value)
    real = sp.together(sp.expand(sp.re(value)))
    imag = sp.together(sp.expand(sp.im(value)))
    if imag != 0:
        imag = sp.simplify(imag)
    if imag != 0:
        raise ValueError(f"{label} is not real: {value}")
    split = _as_q_sqrt29(real)
    if split is None:
        raise ValueError(f"{label} is outside Q(sqrt(29)): {real}")
    rational, radical = split
    if radical == 0:
        return sp.Rational(rational)
    return rational + radical * _SQRT29


def _exact_scalar(value):
    split = _as_q_sqrt29(value)
    if split is not None:
        rational, radical = split
        if radical == 0:
            return _rational(rational)
        return str(sp.simplify(rational + radical * _SQRT29))
    value = sp.simplify(value)
    if value.is_rational:
        return _rational(value)
    return str(value)


def _sqrt29_bounds(bits=64):
    scale = 1 << bits
    root = int(sp.integer_nthroot(29 * scale * scale, 2)[0])
    lower, upper = sp.Rational(root, scale), sp.Rational(root + 1, scale)
    if not (lower > 0 and lower ** 2 <= 29 <= upper ** 2):
        raise RuntimeError("sqrt(29) bounds failed")
    return lower, upper


def _bound_linear_sqrt29(rational, radical, lower_sqrt, upper_sqrt):
    samples = (radical * lower_sqrt, radical * upper_sqrt)
    return rational + min(samples), rational + max(samples)


def _abs_lower(value, bits=64):
    """Return `(magnitude, positive rational lower bound)` or `(0, 0)`."""
    split = _as_q_sqrt29(value)
    if split is None:
        simplified = sp.simplify(value)
        if simplified == 0:
            return sp.Integer(0), sp.Integer(0)
        if simplified.is_rational:
            magnitude = abs(sp.Rational(simplified))
            return magnitude, magnitude
        raise ValueError(f"current is outside Q(sqrt(29)): {simplified}")
    rational, radical = split
    if rational == 0 and radical == 0:
        return sp.Integer(0), sp.Integer(0)
    lower_sqrt, upper_sqrt = _sqrt29_bounds(bits)
    lower, upper = _bound_linear_sqrt29(rational, radical, lower_sqrt, upper_sqrt)
    if lower > 0:
        return rational + radical * _SQRT29, lower
    if upper < 0:
        return -(rational + radical * _SQRT29), -upper
    if bits >= 4096:
        raise ValueError("current sign was not separated by rational sqrt(29) bounds")
    return _abs_lower(value, bits * 2)


def _rational_sqrt_upper(value):
    rational = sp.Rational(sp.together(value))
    if rational < 0:
        raise ValueError("square root of a negative value")
    if rational == 0:
        return sp.Integer(0)
    product = rational.p * rational.q
    root, exact = sp.integer_nthroot(product, 2)
    if not exact:
        root += 1
    return sp.Rational(int(root), rational.q)


def _sqrt29_expression_upper(value, bits=64):
    """Rational upper bound of a value in Q(sqrt(29)), not of its square root."""
    value = sp.simplify(sp.together(value))
    split = _as_q_sqrt29(value)
    if split is None:
        if value.is_rational:
            return sp.Rational(value)
        raise ValueError(f"cannot bound {value}")
    rational, radical = split
    if radical == 0:
        return sp.Rational(rational)
    lower_sqrt, upper_sqrt = _sqrt29_bounds(bits)
    _lower, upper = _bound_linear_sqrt29(rational, radical, lower_sqrt, upper_sqrt)
    if upper < 0:
        raise ValueError("negative Frobenius square")
    return upper


def _channel_key(left, right):
    return f"{left[0]}:{left[1]}->{right[0]}:{right[1]}"


def _mode_currents(operator, state, projectors, pairs):
    currents = {}
    for left, right in pairs:
        currents[(left, right)] = sp.simplify(_current(
            operator, state, left[0], projectors[left[1]], right[0], projectors[right[1]]
        ))
    return currents


def _offblocks_zero(state, depth=3, block=2):
    for row in range(depth):
        for col in range(depth):
            if row == col:
                continue
            row_slice = slice(row * block, (row + 1) * block)
            col_slice = slice(col * block, (col + 1) * block)
            if state[row_slice, col_slice] != sp.zeros(block):
                return False
    return True


def _population(state, region, block=2):
    start = region * block
    return _real_rational(sp.trace(state[start:start + block, start:start + block]), f"population {region}")


def _population_rate(operator, state, region, block=2):
    start = region * block
    commutator = operator * state - state * operator
    return _real_rational(-sp.I * sp.trace(commutator[start:start + block, start:start + block]), f"rate {region}")


def _polynomial_annihilates(operator, polynomial):
    """Horner evaluation. Monic degree-n annihilation certifies the characteristic polynomial."""
    coefficients = polynomial.all_coeffs()
    accumulator = sp.eye(operator.rows) * coefficients[0]
    for coefficient in coefficients[1:]:
        accumulator = accumulator * operator + coefficient * sp.eye(operator.rows)
    return _matrix_is_zero(accumulator)


def _real_characteristic_polynomial(operator):
    symbol = sp.Symbol("x")
    polynomial = sp.Poly(sp.expand(operator.charpoly(symbol).as_expr()), symbol)
    converted = []
    for coefficient in polynomial.all_coeffs():
        converted.append(_real_rational(coefficient, "characteristic coefficient"))
    rational = sp.Poly.from_list(converted, symbol, domain=sp.QQ)
    if rational.degree() != operator.rows or rational.LC() != 1:
        raise ValueError("characteristic polynomial is not monic of full degree")
    return rational


def _peval(polynomial, point):
    accumulator = 0
    for coefficient in polynomial.all_coeffs():
        accumulator = accumulator * point + coefficient
    return sp.together(accumulator)


def _shrink_root_bracket(polynomial, left, right):
    left, right = sp.Rational(left), sp.Rational(right)
    if left == right:
        return left, right
    left_value, right_value = _peval(polynomial, left), _peval(polynomial, right)
    if left_value == 0:
        return left, left
    if right_value == 0:
        return right, right
    midpoint = (left + right) / 2
    middle_value = _peval(polynomial, midpoint)
    if middle_value == 0:
        return midpoint, midpoint
    left_count = int(polynomial.count_roots(left, midpoint))
    right_count = int(polynomial.count_roots(midpoint, right))
    if left_count == 1 and right_count == 0:
        return left, midpoint
    if left_count == 0 and right_count == 1:
        return midpoint, right
    if middle_value == 0 or (left_count == 1 and right_count == 1):
        return midpoint, midpoint
    raise ValueError("isolating bracket did not contain one real root")


def isolating_gap(operator):
    """Positive rational lower bound on the eigenvalue gap, from isolating intervals."""
    polynomial = _real_characteristic_polynomial(operator)
    if not _polynomial_annihilates(operator, polynomial):
        raise ValueError("characteristic polynomial does not annihilate J")
    if int(polynomial.count_roots()) != polynomial.degree():
        raise ValueError("J has fewer real roots than its dimension")
    raw_intervals = polynomial.intervals()
    if len(raw_intervals) != polynomial.degree() or any(multiplicity != 1 for _, multiplicity in raw_intervals):
        raise ValueError("eigenvalue isolating intervals are not six simple roots")
    brackets = [(sp.Rational(left), sp.Rational(right)) for (left, right), _ in raw_intervals]
    brackets.sort()
    for _ in range(polynomial.degree() * 8):
        if all(brackets[index][1] < brackets[index + 1][0] for index in range(len(brackets) - 1)):
            break
        tightened = []
        for index, (left, right) in enumerate(brackets):
            overlaps = (
                (index > 0 and brackets[index - 1][1] >= left)
                or (index + 1 < len(brackets) and right >= brackets[index + 1][0])
            )
            tightened.append(_shrink_root_bracket(polynomial, left, right) if overlaps else (left, right))
        brackets = tightened
    else:
        raise ValueError("eigenvalue intervals did not separate")
    if any(int(polynomial.count_roots(left, right)) != 1 for left, right in brackets):
        raise ValueError("a separated interval does not contain exactly one root")
    gaps = [brackets[index + 1][0] - brackets[index][1] for index in range(len(brackets) - 1)]
    gap = min(gaps)
    if gap <= 0:
        raise ValueError("eigenvalue gap lower bound is not positive")
    return {
        "gap": gap,
        "intervals": brackets,
        "polynomial": [str(coefficient) for coefficient in polynomial.all_coeffs()],
    }


def _current_frobenius_square(operator, region_a, projector_a, region_b, projector_b, block=2):
    """`||K||_F^2` for `K = -i(Pa J Pb - Pb J Pa)`.

    Off-diagonal region blocks have disjoint support, so the two directed
    blocks add in the Frobenius square. Projector operator norms are not required.
    """
    a0, b0 = region_a * block, region_b * block
    jab = operator[a0:a0 + block, b0:b0 + block]
    jba = operator[b0:b0 + block, a0:a0 + block]
    forward = sp.trace(projector_b * jba * projector_a * jab)
    backward = sp.trace(projector_a * jab * projector_b * jba)
    return sp.simplify(forward + backward)


def _bitlength(matrix):
    bits = 1
    for row in range(matrix.rows):
        for col in range(matrix.cols):
            entry = sp.together(matrix[row, col])
            for part in (sp.re(entry), sp.im(entry)):
                rational = sp.Rational(sp.together(part))
                bits = max(bits, rational.p.bit_length(), rational.q.bit_length())
    return bits


def _mp_matrix(matrix, precision):
    mp.mp.dps = precision
    converted = mp.matrix(matrix.rows, matrix.cols)
    for row in range(matrix.rows):
        for col in range(matrix.cols):
            entry = sp.together(matrix[row, col])
            real, imag = sp.Rational(sp.re(entry)), sp.Rational(sp.im(entry))
            converted[row, col] = mp.mpc(mp.mpf(real.p) / mp.mpf(real.q), mp.mpf(imag.p) / mp.mpf(imag.q))
    return converted


def _spectral_dephasing_error(operator, state, mean):
    """Compare the exact commutant projection with a high-precision eigenbasis pinch."""
    bits = max(_bitlength(operator), _bitlength(state), _bitlength(mean))
    precision = min(240, max(40, bits // 2 + 20))
    generated = _mp_matrix(operator, precision)
    initial = _mp_matrix(state, precision)
    exact = _mp_matrix(mean, precision)
    eigenvalues, vectors = mp.eigh(generated)
    rotated = vectors.H * initial * vectors
    diagonal = mp.diag([rotated[index, index] for index in range(rotated.rows)])
    pinched = vectors * diagonal * vectors.H
    residual = exact - pinched
    squared = mp.mpf(0)
    for row in range(residual.rows):
        for col in range(residual.cols):
            squared += abs(residual[row, col]) ** 2
    error = float(mp.sqrt(squared))
    gaps = [float(eigenvalues[index + 1] - eigenvalues[index]) for index in range(len(eigenvalues) - 1)]
    return {
        "frobenius_error": error,
        "precision_digits": precision,
        "entry_bitlength": bits,
        "numeric_min_gap": min(gaps),
    }


def _response_pair(operator, state, spectral):
    numeric_operator = _as_numpy(operator)
    numeric_state = _as_numpy(state)
    reduced, _, coupling, _ = _filtered(numeric_operator, numeric_state, spectral)
    resolvent = np.linalg.inv(spectral * np.eye(numeric_operator.shape[0], dtype=np.complex128) - numeric_operator)
    full_rows = resolvent[:2, :]
    full = full_rows @ numeric_state @ full_rows.conj().T
    exterior = numeric_operator[2:, 2:]
    exterior_resolvent = np.linalg.inv(spectral * np.eye(exterior.shape[0], dtype=np.complex128) - exterior)
    self_energy = coupling @ exterior_resolvent @ coupling.conj().T
    return {
        "reduced": reduced,
        "full": full,
        "agreement": _frobenius(reduced - full),
        "self_energy_versus_covariance": _frobenius(self_energy - numeric_state[:2, :2]),
        "response_frobenius": _frobenius(reduced),
    }


def _finite_time_bound(operator, state, mean, current, channel, projectors):
    """Remainder bound `||C_T - Cbar||_F <= 2 ||Y||_F / (gamma T)`.

    `Y` is off-diagonal in the eigenbasis because `Cbar` is the commutant
    projection. Each oscillation is at most `2/(|lambda_m-lambda_n| T)`.
    """
    spectrum = isolating_gap(operator)
    deviation = state - mean
    if not _matrix_is_zero(deviation - deviation.H):
        raise ValueError("deviation is not Hermitian")
    frobenius_square = _real_rational(sp.trace(deviation * deviation), "deviation Frobenius square")
    deviation_upper = _rational_sqrt_upper(frobenius_square)
    (region_a, mode_a), (region_b, mode_b) = channel
    current_square = _current_frobenius_square(
        operator, region_a, projectors[mode_a], region_b, projectors[mode_b]
    )
    embedded_a = _as_numpy(_embed(projectors[mode_a], region_a))
    embedded_b = _as_numpy(_embed(projectors[mode_b], region_b))
    generated = _as_numpy(operator)
    amplitude = embedded_a @ generated @ embedded_b
    numeric_square = float(np.linalg.norm(-1j * (amplitude - amplitude.conj().T), ord="fro") ** 2)
    algebraic_square = complex(sp.N(current_square, 40))
    if abs(algebraic_square.imag) > 1e-8 or abs(algebraic_square.real - numeric_square) > max(1e-6, 1e-6 * abs(numeric_square)):
        raise ValueError(
            f"current-operator norm disagrees with K: algebraic {algebraic_square}, numeric {numeric_square}"
        )
    current_square_upper = _sqrt29_expression_upper(current_square)
    if current_square_upper < 0:
        raise ValueError("current-operator square upper bound is negative")
    operator_upper = _rational_sqrt_upper(current_square_upper)
    magnitude, lower = _abs_lower(current)
    if lower <= 0:
        raise ValueError("selected mean current has no positive lower bound")
    gap = spectrum["gap"]
    threshold = sp.together(2 * deviation_upper * operator_upper / (gap * lower))
    return {
        "obtained": True,
        "gamma_lower_exact": _rational(gap),
        "isolating_intervals": [[_rational(left), _rational(right)] for left, right in spectrum["intervals"]],
        "characteristic_polynomial_high_to_low": spectrum["polynomial"],
        "deviation_frobenius_square_exact": _rational(frobenius_square),
        "deviation_frobenius_upper_exact": _rational(deviation_upper),
        "current_channel": _channel_key(*channel),
        "current_exact": _exact_scalar(current),
        "current_abs_lower_exact": _rational(lower),
        "current_operator_frobenius_square_exact": _exact_scalar(current_square),
        "current_operator_frobenius_upper_exact": _rational(operator_upper),
        "state_error_at_most": "2*deviation_upper/(gamma*T)",
        "mean_current_error_at_most": "2*deviation_upper*current_operator_upper/(gamma*T)",
        "sign_separation_if_T_greater_than": _rational(threshold),
        "separator_magnitude_used": _exact_scalar(magnitude),
        "reading": "finite-window averaging error bound, not a physical cooling law",
    }


def _projector_checks(local, larger, smaller, plus, minus):
    return {
        "onsite_projectors_resolve_identity": bool(sp.simplify(plus + minus - sp.eye(2)) == sp.zeros(2)),
        "onsite_projectors_hermitian_idempotent": bool(
            plus == plus.H and minus == minus.H
            and sp.simplify(plus * plus - plus) == sp.zeros(2)
            and sp.simplify(minus * minus - minus) == sp.zeros(2)
            and sp.simplify(plus * minus) == sp.zeros(2)
        ),
        "onsite_projectors_match_eigenvalues": bool(
            sp.simplify(plus * local - larger * plus) == sp.zeros(2)
            and sp.simplify(minus * local - smaller * minus) == sp.zeros(2)
        ),
    }


def _alpha_list(alpha):
    return [_rational(alpha[index, 0]) for index in range(alpha.rows)]


def _fill_currents(operator, projectors, modes, channels, degree, powers):
    """Linear map from one local mode filling `P_{r,s}` to adjacent mean currents."""
    columns = {}
    projections = {}
    for mode in modes:
        basis = _embed(projectors[mode[1]], mode[0])
        projected = commutant_projection(operator, basis, powers)
        if projected["mean"] is None or not projected["moments_matched"] or not projected["commutes"]:
            raise ValueError(f"mode filling {mode} did not project onto the commutant")
        currents = _mode_currents(operator, projected["mean"], projectors, channels)
        columns[_channel_key(mode, mode).split("->")[0]] = {
            _channel_key(left, right): _exact_scalar(value) for (left, right), value in currents.items()
        }
        projections[mode] = currents
        if projected["mean"].rows != degree:
            raise ValueError("unexpected projection shape")
    return columns, projections


def compute_result():
    """Time-averaged circulation of `C0 = I/2 + (1/4)(P0 - P2)` on the frozen window."""
    local = sp.Matrix([[1, sp.I / 5], [-sp.I / 5, 2]])
    link = sp.Matrix([[sp.Rational(1, 4), sp.I / 7], [sp.Rational(1, 9), sp.Rational(1, 6)]])
    omega = sp.Rational(3, 2)
    delta = sp.Rational(1, 4)
    depth = 3
    spectral_symbol = 1 + 2 * sp.I
    spectral_numeric = complex(spectral_symbol)
    probe = sp.Rational(1, 8)
    window = finite_window(local, link, omega, 0, depth)
    degree = window.rows
    larger, smaller, plus, minus = _onsite_projectors(local)
    projectors = {"plus": plus, "minus": minus}
    modes = [(region, name) for region in range(depth) for name in ("plus", "minus")]
    ordered_pairs = [(left, right) for left in modes for right in modes]
    adjacent = [
        (left, right)
        for left in modes
        for right in modes
        if abs(left[0] - right[0]) == 1
    ]
    cut_pairs = [
        ((0, left_mode), (region, right_mode))
        for left_mode in ("plus", "minus")
        for region in (1, 2)
        for right_mode in ("plus", "minus")
    ]
    far_pairs = [
        ((1, left_mode), (2, right_mode))
        for left_mode in ("plus", "minus")
        for right_mode in ("plus", "minus")
    ]

    initial = regional_preparation(delta, depth)
    reverse = regional_preparation(-delta, depth)
    uniform = regional_preparation(0, depth)
    imbalance = _regional_shift(depth)
    powers = _powers(window, 2 * degree - 2)
    base = commutant_projection(window, initial, powers)
    reverse_mean = commutant_projection(window, reverse, powers)
    uniform_mean = commutant_projection(window, uniform, powers)
    imbalance_mean = commutant_projection(window, imbalance, powers)
    mean = base["mean"]

    derivative = -sp.I * (window * initial - initial * window)
    initial_currents = _mode_currents(window, initial, projectors, ordered_pairs)
    initial_rates = _mode_currents(window, derivative, projectors, adjacent)
    mean_currents = _mode_currents(window, mean, projectors, ordered_pairs)
    reverse_currents = _mode_currents(window, reverse_mean["mean"], projectors, adjacent)
    uniform_currents = _mode_currents(window, uniform_mean["mean"], projectors, adjacent)
    imbalance_currents = _mode_currents(window, imbalance_mean["mean"], projectors, adjacent)

    def _exchange(currents, pairs):
        total = sp.Integer(0)
        lower = sp.Integer(0)
        signed = []
        for pair in pairs:
            magnitude, floor = _abs_lower(currents[pair])
            total += magnitude
            lower += floor
            signed.append(currents[pair])
        return sp.simplify(total), lower, signed

    cut_exchange, cut_lower, cut_signed = _exchange(mean_currents, cut_pairs)
    far_exchange, far_lower, _ = _exchange(mean_currents, far_pairs)
    net_cut = sp.simplify(sum(cut_signed, sp.Integer(0)))
    filling_map = {}
    map_imbalance = {}
    filling_error = None
    try:
        filling_map, filling_currents = _fill_currents(window, projectors, modes, adjacent, degree, powers)
        for pair in adjacent:
            combined = (
                filling_currents[(0, "plus")][pair] + filling_currents[(0, "minus")][pair]
                - filling_currents[(2, "plus")][pair] - filling_currents[(2, "minus")][pair]
            )
            map_imbalance[pair] = sp.simplify(combined - imbalance_currents[pair])
    except Exception as exc:
        filling_error = f"{type(exc).__name__}: {exc}"
        filling_map = {}
        map_imbalance = {}

    populations = [_population(initial, region) for region in range(depth)]
    mean_populations = [_population(mean, region) for region in range(depth)]
    initial_population_rates = [_population_rate(window, initial, region) for region in range(depth)]
    mean_population_rates = [_population_rate(window, mean, region) for region in range(depth)]
    content = _real_rational(sp.trace(initial), "content")
    energy = _real_rational(sp.trace(window * initial), "energy")
    mean_content = _real_rational(sp.trace(mean), "mean content")
    mean_energy = _real_rational(sp.trace(window * mean), "mean energy")

    dephasing = _spectral_dephasing_error(window, initial, mean)
    response = _response_pair(window, mean, spectral_numeric)
    perturbed = sp.Matrix(window)
    perturbed[2, 2] += probe
    perturbed_projection = commutant_projection(perturbed, initial)
    perturbed_response = _response_pair(perturbed, perturbed_projection["mean"], spectral_numeric)
    frozen_response = _response_pair(perturbed, mean, spectral_numeric)
    blocks_fixed = bool(
        perturbed[:2, :2] == window[:2, :2]
        and perturbed[:2, 2:] == window[:2, 2:]
        and sp.simplify(perturbed[2, 2] - window[2, 2] - probe) == 0
    )

    circulation_pairs = [pair for pair in adjacent if mean_currents[pair] != 0]
    selected = None
    selected_lower = None
    for pair in adjacent:
        _, floor = _abs_lower(mean_currents[pair])
        if floor > 0 and (selected_lower is None or floor > selected_lower):
            selected, selected_lower = pair, floor

    bound = None
    bound_error = None
    if selected is not None and base["gram_determinant"] != 0 and base["commutes"] and base["moments_matched"]:
        try:
            bound = _finite_time_bound(window, initial, mean, mean_currents[selected], selected, projectors)
        except Exception as exc:
            bound_error = f"{type(exc).__name__}: {exc}"

    any_initial_rate = any(value != 0 for value in initial_rates.values())
    cut_positive = bool(cut_lower > 0)
    persistent = bool(selected is not None)
    identities = {
        "hermitian_window": bool(window == window.H),
        "dimension_six": bool(degree == 6 and window.cols == 6),
        "no_corner_block": bool(window[:2, 4:] == sp.zeros(2) and window[4:, :2] == sp.zeros(2)),
        "link_j02_is_one_quarter": bool(sp.simplify(window[0, 2] - sp.Rational(1, 4)) == 0),
        "gram_determinant_nonzero": bool(base["gram_determinant"] != 0),
        "mean_commutes_with_generator": bool(base["commutes"]),
        "power_moments_match_through_degree_five": bool(base["moments_matched"]),
        "initial_generator_commutator_nonzero": bool(not _matrix_is_zero(window * initial - initial * window)),
        "initial_state_has_no_interregion_coherence": bool(_offblocks_zero(initial)),
        "initial_admissible": bool(set(populations) == {sp.Rational(3, 2), sp.Integer(1), sp.Rational(1, 2)} and content == 3),
        "initial_interregion_currents_zero": bool(all(
            initial_currents[((left_region, left_mode), (right_region, right_mode))] == 0
            for left_region, left_mode in modes
            for right_region, right_mode in modes
            if left_region != right_region
        )),
        "initial_current_derivative_nonzero": bool(any_initial_rate),
        "content_and_energy_conserved": bool(
            _real_rational(sp.trace(derivative), "content derivative") == 0
            and _real_rational(sp.trace(window * derivative), "energy derivative") == 0
            and mean_content == content
            and mean_energy == energy
        ),
        "mean_inventory_rates_zero": bool(all(rate == 0 for rate in mean_population_rates)),
        "uniform_preparation_is_stationary_and_current_free": bool(
            uniform_mean["commutes"]
            and uniform_mean["mean"] == uniform
            and all(value == 0 for value in uniform_currents.values())
            and _alpha_list(uniform_mean["alpha"]) == ["1/2", "0", "0", "0", "0", "0"]
        ),
        "sign_reversal_negates_adjacent_currents": bool(all(
            sp.simplify(reverse_currents[pair] + mean_currents[pair]) == 0 for pair in adjacent
        )),
        "currents_scale_exactly_with_delta": bool(all(
            sp.simplify(mean_currents[pair] - delta * imbalance_currents[pair]) == 0 for pair in adjacent
        )),
        "reverse_mean_matches_linear_preparation": bool(
            _matrix_is_zero(reverse_mean["mean"] - (uniform + (-delta) * imbalance_mean["mean"]))
        ),
        "filling_map_reproduces_imbalance_currents": bool(
            filling_error is None
            and len(map_imbalance) == len(adjacent)
            and all(value == 0 for value in map_imbalance.values())
        ),
        "currents_antisymmetric": bool(all(
            sp.simplify(mean_currents[(left, right)] + mean_currents[(right, left)]) == 0
            for left in modes for right in modes
        )),
        "mean_mode_divergence_zero": bool(all(
            sp.simplify(sum(mean_currents[(mode, other)] for other in modes)) == 0
            for mode in modes
        )),
        "corner_mean_currents_zero": bool(all(
            mean_currents[((0, left_mode), (2, right_mode))] == 0
            for left_mode in ("plus", "minus")
            for right_mode in ("plus", "minus")
        )),
        "intraregion_mean_currents_zero": bool(all(
            mean_currents[((region, "plus"), (region, "minus"))] == 0 for region in range(depth)
        )),
        "spectral_dephasing_matches_projection": bool(dephasing["frobenius_error"] < 1e-8),
        "full_filter_matches_reduced": bool(response["agreement"] < 1e-8),
        "perturbed_full_filter_matches_reduced": bool(perturbed_response["agreement"] < 1e-8),
        "perturbed_generator_span_remains_full": bool(
            perturbed_projection["gram_determinant"] != 0
            and perturbed_projection["commutes"]
            and perturbed_projection["moments_matched"]
        ),
        "exterior_probe_keeps_initial_state_and_retained_coupling": bool(blocks_fixed),
        "self_energy_differs_from_covariance_block": bool(response["self_energy_versus_covariance"] > 1e-8),
        "admissible_mean_by_commutant_pinch": bool(
            base["gram_determinant"] != 0 and base["commutes"] and base["moments_matched"]
        ),
    }
    identities.update(_projector_checks(local, larger, smaller, plus, minus))
    identities_hold = all(identities.values())
    bound_obtained = bool(bound and bound.get("obtained"))
    if not identities_hold:
        verdict = FAIL
    elif not persistent:
        verdict = ZERO
    elif not bound_obtained:
        verdict = MEAN_ONLY
    else:
        verdict = PASS

    inputs = {
        "H": [_pair(local[row, col]) for row in range(2) for col in range(2)],
        "B": [_pair(link[row, col]) for row in range(2) for col in range(2)],
        "omega": str(omega),
        "depth": depth,
        "first": 0,
        "dimension": degree,
        "delta": _rational(delta),
        "reverse_delta": _rational(-delta),
        "uniform_delta": "0",
        "preparation": "I/2 + delta*(P_region0 - P_region2)",
        "interregion_coherence": "none",
        "spectral_parameter": "1+2*I",
        "retained_region": 0,
        "exterior_probe": {"index": [2, 2], "delta": _rational(probe)},
        "gram": "G_ab=Tr(J**(a+b)), b_a=Tr(J**a*C0), a,b=0..5",
        "current": "j_ab=2*Im*Tr(P_a*J*P_b*C)",
    }
    required = ("H", "B", "omega", "depth", "delta", "preparation", "spectral_parameter", "exterior_probe")
    if any(key not in inputs or inputs[key] in (None, "", []) for key in required):
        verdict = FAIL
        identities["declared_inputs_present"] = False
    else:
        identities["declared_inputs_present"] = True
        identities_hold = identities_hold and True

    metrics = {
        "Gbar_region0_cut_exact": _exact_scalar(cut_exchange),
        "Gbar_region0_cut_lower_exact": _rational(cut_lower),
        "Gbar_region0_cut_positive": cut_positive,
        "net_region0_cut_exact": _exact_scalar(net_cut),
        "Gbar_region1_to_region2_exact": _exact_scalar(far_exchange),
        "Gbar_region1_to_region2_lower_exact": _rational(far_lower),
        "persistent_mean_mode_circulation": persistent,
        "nonzero_mean_channels": [_channel_key(*pair) for pair in circulation_pairs],
        "mean_adjacent_currents_exact": {
            _channel_key(*pair): _exact_scalar(mean_currents[pair]) for pair in adjacent
        },
        "initial_adjacent_current_derivatives_exact": {
            _channel_key(*pair): _exact_scalar(initial_rates[pair]) for pair in adjacent
        },
        "projection_coefficients_exact": _alpha_list(base["alpha"]),
        "gram_determinant_exact": str(base["gram_determinant"]),
        "total_content_exact": _rational(content),
        "total_energy_exact": _rational(energy),
        "initial_region_populations_exact": [_rational(value) for value in populations],
        "mean_region_populations_exact": [_rational(value) for value in mean_populations],
        "initial_region_population_rates_exact": [_rational(value) for value in initial_population_rates],
        "mean_region_population_rates_exact": [_rational(value) for value in mean_population_rates],
        "uniform_adjacent_currents_max_abs_exact": _exact_scalar(max(uniform_currents.values(), key=lambda item: abs(sp.N(item)))),
        "delta_scale_factor_exact": _rational(delta),
        "spectral_dephasing_frobenius_error": dephasing["frobenius_error"],
        "spectral_dephasing_digits": dephasing["precision_digits"],
        "projection_entry_bitlength": dephasing["entry_bitlength"],
        "numeric_eigenvalue_gap_control": dephasing["numeric_min_gap"],
        "filtered_mean_response_frobenius": response["response_frobenius"],
        "full_versus_reduced_frobenius": response["agreement"],
        "self_energy_versus_covariance_frobenius": response["self_energy_versus_covariance"],
        "perturbed_filtered_response_frobenius": perturbed_response["response_frobenius"],
        "perturbed_full_versus_reduced_frobenius": perturbed_response["agreement"],
        "perturbed_versus_baseline_response_frobenius": _frobenius(perturbed_response["reduced"] - response["reduced"]),
        "frozen_cbar_shortcut_versus_recomputed_frobenius": _frobenius(
            frozen_response["reduced"] - perturbed_response["reduced"]
        ),
        "frozen_cbar_shortcut_versus_baseline_frobenius": _frobenius(frozen_response["reduced"] - response["reduced"]),
        "mode_filling_to_mean_current_map": filling_map,
    }
    assumptions = [
        "The generator is the frozen finite_window(H, B, 3/2, 0, 3); geometry is an input.",
        "C0 is the incoherent regional population I/2 + delta (P_region0 - P_region2), not the v1 cubic stationary covariance.",
        "The recorded mean is the Cesaro average of unitary conjugates, obtained as the Frobenius projection onto polynomials in J. It is not a claim that C(t) itself relaxes pointwise.",
        "A positive cut exchange is balanced mode current on this finite window. It is not particle production, heat, entropy production, a clock, antimatter, or an eternity proof.",
        "F Cbar F^dagger is a resolvent-weighted response at z=1+2i, not an equal-time occupation and not a stress tensor.",
        "The exterior self-energy Sigma is not the covariance C or Cbar.",
        "The exterior probe changes J[2,2] only. The initial C0, retained block, and coupling block stay fixed; the perturbed mean is recomputed.",
        "The finite-T estimate, when present, bounds the averaging remainder. It is not a physical cooling law.",
        "The local incoming gate remains OPEN. Its campaign is not used.",
    ]
    if persistent:
        remaining = (
            "Nonzero Cesaro mode current on this finite window is not connected to a self-consistent "
            "radius, a matter-antimatter identification, or a continuing cosmological history."
        )
    else:
        remaining = (
            "This incoherent regional preparation has zero persistent mode circulation on the frozen "
            "window. A different admissible preparation can still be tested; no search was run. The "
            "result is not a radius, matter-antimatter, or cosmological identification."
        )
    if bound_error:
        remaining = f"{remaining} Finite-T bound unfinished: {bound_error}"
    return {
        "schema": SCHEMA,
        "verdict": verdict,
        "verdict_scope": "this finite six-mode window and this incoherent regional preparation only",
        "question": (
            "Does unitary evolution of C0=I/2+(1/4)(P_region0-P_region2), with no inter-region "
            "coherence, produce nonzero persistent time-averaged mode circulation?"
        ),
        "answer": "yes" if persistent and identities_hold else ("no" if identities_hold else "unresolved"),
        "physical_local_gate": "OPEN",
        "inputs": inputs,
        "checks": identities,
        "metrics": metrics,
        "finite_time_bound": bound,
        "finite_time_bound_error": bound_error,
        "filling_control_error": filling_error,
        "preparation_conditions": {
            "family": "incoherent regional imbalance",
            "formula": "C0(delta)=I/2+delta*(P_region0-P_region2)",
            "delta": _rational(delta),
            "controls": ["delta=-1/4 reverses every adjacent mean current", "delta=0 remains I/2 and carries no mode current"],
            "coherence": "regional blocks are scalar and off-blocks are zero",
            "admissible_eigenvalues": ["3/4", "3/4", "1/2", "1/2", "1/4", "1/4"],
            "evolution": "C(t)=exp(-i*J*t)*C0*exp(i*J*t)",
            "conserved": ["Tr(C)", "Tr(J*C)"],
        },
        "assumptions": assumptions,
        "remaining_connection": remaining,
    }
