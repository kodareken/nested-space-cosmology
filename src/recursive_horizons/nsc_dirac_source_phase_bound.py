"""Directed Taylor enclosure of the unit-angular source-phase derivative f_s,z.

This is geometry quadrature of

    f_z / ell^2 = -int_1^{rho_up} r_z / r^3 (rho, z - s D(rho)) d rho

with interval Taylor remainder. It does not evolve a source, fit a stress, or
close the physical local gate. Missing source-tail and other error components
remain unset.
"""
from dataclasses import dataclass, field
from time import process_time

from flint import acb, arb, arb_series, ctx

from . import nsc_ks_ball_geometry as G
from .nsc_ks_ball_operator import AnalyticRadiusFamily
from .nsc_ks_profile_identity import profile_identity


SOURCE_SIGNS = (1, -1)
TAYLOR_DEGREE = 13
REMAINDER_ORDER = 14
MAX_DEPTH = 24
MAX_CELLS = 8000
BITS = 120
CPU_CAP = 60.0
ANGULAR_SQUARE_WEIGHT = 107952
MIN_DYADIC_SCALE = 52
D_DERIVATIVE_BOUND = arb(25) / 16


def unit_angular_target():
    """2e-16 as a directed working-precision ball, not a decimal string."""
    return arb(2) / arb(10)**16
PHYSICAL_LOCAL_GATE = (
    'OPEN: one-point directed f_s,z enclosure only; source-tail, field and '
    'between-node errors remain parent-owned; no assembly correction')


class BudgetExceeded(RuntimeError):
    """A declared cell, depth or CPU limit stopped the enclosure."""


def _inv_a2_integrand(z, analytic):
    value = G.inv_a_integrand(z, analytic)
    return value * value


def distance_series(rho_ball, order):
    """D=int_1^rho a^{-2} d rho, with interval value from |D'|<=25/16."""
    G.restricted_axial_lower_bound()
    with G._JetWork(order) as n:
        rho = G._rho_ball(rho_ball)
        higher = G._finite_series(
            G.background_series(rho, n).inv_a2.integral(), 'D tail', n)
        integral = acb.integral(_inv_a2_integrand, acb(1), acb(rho.mid()))
        if not integral.is_finite() or not integral.imag.contains(0):
            raise ArithmeticError('acb.integral of 1/a^2 did not return a real enclosure')
        value = integral.real + D_DERIVATIVE_BOUND * rho.rad() * arb(0, 1)
        coeffs = G._coeffs(higher, n)
        coeffs[0] = value
        return G._finite_series(arb_series(coeffs, prec=n + 1), 'D', n)


def _compose(series, increment, order):
    """Composition at the same basepoint; increment has exact zero constant."""
    if increment[0] != 0:
        raise ValueError('formal increment must have exact zero constant')
    result = arb_series([0], prec=order + 1)
    for coefficient in reversed(G._coeffs(series, order)):
        result = result * increment + coefficient
    return result


def _dyadic_parts(value, name):
    number = float(value)
    numerator, denominator = number.as_integer_ratio()
    if denominator & (denominator - 1):
        raise ValueError(f'{name} must be a binary dyadic')
    return numerator, denominator.bit_length() - 1


def _at_scale(numerator, native_scale, scale):
    if scale < native_scale:
        raise ValueError('insufficient dyadic scale')
    return numerator << (scale - native_scale)


def dyadic_rho_cover(rho_left, rho_right):
    """Greedy power-of-two cover of an ordered binary-dyadic rho interval."""
    left_n, left_s = _dyadic_parts(rho_left, 'rho_left')
    right_n, right_s = _dyadic_parts(rho_right, 'rho_right')
    scale = max(left_s, right_s, MIN_DYADIC_SCALE)
    start = _at_scale(left_n, left_s, scale)
    stop = _at_scale(right_n, right_s, scale)
    if stop <= start:
        raise ValueError('ordered positive-length dyadic rho interval required')
    cells = []
    pos = start
    while pos < stop:
        align = pos & -pos
        width = 1
        remaining = stop - pos
        while width * 2 <= remaining and width * 2 <= align:
            width *= 2
        cells.append((pos, pos + width))
        pos += width
    return scale, tuple(cells)


def _cell_ball(start, stop, scale):
    ball = arb((start + stop, -(scale + 1)), (stop - start, -(scale + 1)))
    if not (ball >= 1 and ball <= arb(33) / 32):
        raise G.SubdivisionNeeded('dyadic cell ball left the owned slab')
    return ball


def _cell_center(start, stop, scale):
    return arb((start + stop, -(scale + 1)))


def _cell_half(start, stop, scale):
    return arb((stop - start, -(scale + 1)))


def integrand_rho_series(model, rho_ball, z_target, sign, order):
    """Rho-jet of -r_z/r^3 along z_char=z-s D, with exact-zero D increment."""
    if sign not in SOURCE_SIGNS:
        raise ValueError('source sign must be +1 or -1')
    profile_order = order + 1
    if profile_order > G.MAX_JET_ORDER:
        raise ValueError('profile composition exceeds MAX_JET_ORDER')
    with G._JetWork(order) as n:
        background = G.background_series(rho_ball, n)
        shift = G.shift_series(rho_ball, n)
        distance = distance_series(rho_ball, n)
        window = G.plateau_series(-shift, model.inner, model.outer, n)
        z_char0 = arb(z_target) - sign * distance[0]
        increment = G._zero_const(-sign * distance, n)
        if increment[0] != 0:
            raise ArithmeticError('characteristic increment lost its exact zero constant')
    with G._JetWork(profile_order):
        w_zjet, u_zjet = model.profile_series(z_char0, profile_order)
        w_z_zjet = w_zjet.derivative()
        u_z_zjet = u_zjet.derivative()
    with G._JetWork(order) as n:
        w = _compose(w_zjet, increment, n)
        u = _compose(u_zjet, increment, n)
        w_z = _compose(w_z_zjet, increment, n)
        u_z = _compose(u_z_zjet, increment, n)
        first = window * shift
        third = window * shift ** 3 / 6
        radius = background.r + first * w + third * u
        radius_z = first * w_z + third * u_z
        if not (radius[0] > 0):
            raise G.SubdivisionNeeded('radius positivity requires a smaller rho cell')
        # Interval inverse of the proved-positive radius; the upper endpoint
        # is never substituted as a lower denominator.
        integrand = -radius_z * (radius.inv() ** 3)
        return G._finite_series(integrand, 'phase integrand', n)


def integrate_dyadic_cell(model, start, stop, scale, z_target, sign, *,
                          taylor_degree=TAYLOR_DEGREE, remainder_order=REMAINDER_ORDER):
    """Even-power center polynomial plus uniform degree-(remainder) remainder."""
    if remainder_order != taylor_degree + 1:
        raise ValueError('remainder order must be one more than the Taylor degree')
    if not 0 <= taylor_degree < remainder_order <= G.MAX_JET_ORDER - 1:
        raise ValueError('Taylor degree and remainder exceed the owned jet limit')
    ball = _cell_ball(start, stop, scale)
    half = _cell_half(start, stop, scale)
    center = _cell_center(start, stop, scale)
    polynomial = integrand_rho_series(
        model, center, z_target, sign, taylor_degree)
    remainder_jet = integrand_rho_series(
        model, ball, z_target, sign, remainder_order)
    value = arb(0)
    coeffs = G._coeffs(polynomial, taylor_degree)
    for degree in range(0, taylor_degree + 1, 2):
        value += coeffs[degree] * (2 * half ** (degree + 1)) / arb(degree + 1)
    span = (2 * half ** (remainder_order + 1)) / arb(remainder_order + 1)
    magnitude = abs(G._coeffs(remainder_jet, remainder_order)[remainder_order]) * span
    extra = magnitude.upper()
    enclosed = value + arb((0, 0), extra.man_exp())
    return enclosed, extra


def pack_real_ball(value):
    """Exact dyadic midpoint/radius; binary64 is display only."""
    if not value.is_finite():
        raise ArithmeticError('finite real ball required')
    mid = value.mid().man_exp()
    rad = value.rad().man_exp()
    return {
        'mid_mantissa': str(mid[0]),
        'mid_exponent': int(mid[1]),
        'rad_mantissa': str(rad[0]),
        'rad_exponent': int(rad[1]),
        'binary64_mid': float(value.mid()),
        'binary64_rad': float(value.rad()),
    }


def unpack_real_ball(item):
    return arb(
        (int(item['mid_mantissa']), int(item['mid_exponent'])),
        (int(item['rad_mantissa']), int(item['rad_exponent'])),
    )


def pack_upper(value):
    upper = value.upper()
    if not upper.is_finite():
        raise ArithmeticError('finite directed upper required')
    mantissa, exponent = upper.man_exp()
    return {'mantissa': str(mantissa), 'exponent': int(exponent),
            'binary64_upper': float(upper)}


def unpack_upper(item):
    return arb(int(item['mantissa'])) * arb(2) ** int(item['exponent'])


def aggregate_phase_error_uppers(plus_ball, minus_ball, weight=ANGULAR_SQUARE_WEIGHT):
    """Directed N,beta value-error uppers from unit-angular f_z radii.

    The contraction is the owned edge map with a>=4/5 in the N denominator.
    """
    if weight != ANGULAR_SQUARE_WEIGHT:
        raise ValueError('original aggregate angular-square weight per sign required')
    weight = arb(int(weight))
    two_pi = 2 * arb.pi()
    axial_lower = arb(G.AXIAL_LOWER.numerator) / arb(G.AXIAL_LOWER.denominator)
    radii = plus_ball.rad() + minus_ball.rad()
    beta = weight * radii / two_pi
    newton = beta / axial_lower
    return newton.upper(), beta.upper()


@dataclass
class DirectedSourcePhaseZEnclosure:
    signs: tuple
    z: float
    rho_up: float
    angular: float
    profile_identity: str
    balls: tuple
    remainder_uppers: tuple
    cells: tuple
    max_depth: tuple
    splits: tuple
    cover_cells: int
    dyadic_scale: int
    taylor_degree: int
    remainder_order: int
    target_radius: object
    status: str
    obstruction: str | None
    cpu_seconds: float
    coefficient_error_bound: object
    remainder_error_bound: object
    physical_error_bound: object
    source_tail_error_bound: object
    physical_local_gate: str
    assembly_correction: bool
    diagnostics: dict = field(default_factory=dict)


def _deadline(cpu_limit, started):
    if cpu_limit is None:
        return None
    limit = float(cpu_limit)
    if not (limit > 0):
        raise ValueError('positive CPU limit required')
    return started + limit


def _enclose_sign(model, z_target, sign, scale, cover, *, target, taylor_degree,
                  remainder_order, max_depth, max_cells, deadline):
    length = arb(cover[-1][1] - cover[0][0])
    stack = [(start, stop, 0) for start, stop in reversed(cover)]
    total = arb(0)
    remainder_upper = arb(0)
    cells = 0
    splits = 0
    depth_seen = 0
    while stack:
        if deadline is not None and process_time() >= deadline:
            raise BudgetExceeded('cpu_budget')
        start, stop, depth = stack.pop()
        depth_seen = max(depth_seen, depth)
        width = arb(stop - start)
        local = target * width / length
        try:
            value, extra = integrate_dyadic_cell(
                model, start, stop, scale, z_target, sign,
                taylor_degree=taylor_degree, remainder_order=remainder_order)
            can_split = (
                stop - start > 1 and depth < max_depth
                and cells + len(stack) + 2 <= max_cells)
            if not (extra <= local) and can_split:
                raise G.SubdivisionNeeded('Taylor remainder exceeds the cell budget')
        except G.SubdivisionNeeded as error:
            if stop - start <= 1:
                raise BudgetExceeded('terminal_dyadic_cell') from error
            if depth >= max_depth:
                raise BudgetExceeded('max_depth') from error
            if cells + len(stack) + 2 > max_cells:
                raise BudgetExceeded('max_cells') from error
            mid = (start + stop) // 2
            stack.append((mid, stop, depth + 1))
            stack.append((start, mid, depth + 1))
            splits += 1
            continue
        total += value + arb((0, 0), extra.man_exp())
        remainder_upper += extra
        cells += 1
        if cells > max_cells:
            raise BudgetExceeded('max_cells')
    return total, remainder_upper, cells, depth_seen, splits


def directed_source_phase_z_enclosure(
        family, z, *, rho_up, angular=1.0, bits=BITS,
        taylor_degree=TAYLOR_DEGREE, remainder_order=REMAINDER_ORDER,
        target_radius=None, max_depth=MAX_DEPTH, max_cells=MAX_CELLS,
        cpu_limit=CPU_CAP):
    """Enclose unit-angular f_s,z at one z for both original source signs."""
    if isinstance(z, (bool, list, tuple)) or getattr(z, 'ndim', 0) != 0:
        raise ValueError('one real control z required')
    z = float(z)
    rho_up = float(rho_up)
    angular = float(angular)
    if angular != 1.0:
        raise ValueError('this owner encloses the unit-angular coefficient')
    if isinstance(max_depth, bool) or isinstance(max_cells, bool):
        raise ValueError('integer cell and depth limits required')
    max_depth, max_cells = int(max_depth), int(max_cells)
    if max_depth < 0 or max_cells < 1:
        raise ValueError('nonnegative depth and positive cell limit required')
    model = AnalyticRadiusFamily(family)
    identity = profile_identity(family, include_normal_window=True)
    if model.w_zero and model.U_zero:
        with ctx.workprec(int(bits)):
            zeros = (arb(0), arb(0))
            packed = tuple(pack_real_ball(v) for v in zeros)
            packed_target = pack_upper(unit_angular_target())
        return DirectedSourcePhaseZEnclosure(
            signs=SOURCE_SIGNS, z=z, rho_up=rho_up, angular=angular,
            profile_identity=identity, balls=zeros, remainder_uppers=(arb(0), arb(0)),
            cells=(0, 0), max_depth=(0, 0), splits=(0, 0), cover_cells=0,
            dyadic_scale=MIN_DYADIC_SCALE, taylor_degree=taylor_degree,
            remainder_order=remainder_order, target_radius=unit_angular_target(),
            status=('PASS: directed Taylor enclosure of unit-angular f_s,z at one '
                    'incoming node; physical local gate OPEN'),
            obstruction=None, cpu_seconds=0.0, coefficient_error_bound=packed,
            remainder_error_bound=None, physical_error_bound=None,
            source_tail_error_bound=None, physical_local_gate=PHYSICAL_LOCAL_GATE,
            assembly_correction=False,
            diagnostics={
                'rule': 'identically vanishing r_z; no rho quadrature',
                'gauss_nodes_are_bounds': False,
                'local_embedding_reused': True,
                'global_parent_child_matching_required': False,
                'new_field_or_source_runs': 0,
                'remainder_uppers': (pack_upper(arb(0)), pack_upper(arb(0))),
                'target_radius_upper': packed_target,
            },
        )
    scale, cover = dyadic_rho_cover(1.0, rho_up)
    started = process_time()
    deadline = _deadline(cpu_limit, started)
    balls = []
    remainder_uppers = []
    cells = []
    depths = []
    splits = []
    obstruction = None
    packed = None
    packed_rems = ()
    packed_target = None
    with ctx.workprec(int(bits)):
        saved_cap, saved_prec = ctx.cap, ctx.prec
        target = unit_angular_target() if target_radius is None else arb(target_radius)
        try:
            for sign in SOURCE_SIGNS:
                try:
                    ball, rem, n_cells, depth, n_splits = _enclose_sign(
                        model, z, sign, scale, cover, target=target,
                        taylor_degree=taylor_degree, remainder_order=remainder_order,
                        max_depth=max_depth, max_cells=max_cells, deadline=deadline)
                except BudgetExceeded as error:
                    obstruction = str(error)
                    break
                balls.append(ball)
                remainder_uppers.append(rem)
                cells.append(n_cells)
                depths.append(depth)
                splits.append(n_splits)
            if obstruction is None:
                if len(balls) != len(SOURCE_SIGNS):
                    obstruction = 'incomplete_source_signs'
                else:
                    for ball in balls:
                        if not (ball.rad() <= target):
                            obstruction = 'taylor_remainder_exceeds_target'
                            break
            packed = None if obstruction else tuple(pack_real_ball(v) for v in balls)
            packed_rems = tuple(pack_upper(v) for v in remainder_uppers)
            packed_target = pack_upper(target)
        finally:
            ctx.cap, ctx.prec = saved_cap, saved_prec
    cpu = process_time() - started
    if obstruction is None:
        status = (
            'PASS: directed Taylor enclosure of unit-angular f_s,z at one '
            'incoming node; physical local gate OPEN')
    else:
        status = f'OPEN: {obstruction}; physical local gate OPEN'
    return DirectedSourcePhaseZEnclosure(
        signs=SOURCE_SIGNS,
        z=z,
        rho_up=rho_up,
        angular=angular,
        profile_identity=identity,
        balls=tuple(balls),
        remainder_uppers=tuple(remainder_uppers),
        cells=tuple(cells),
        max_depth=tuple(depths),
        splits=tuple(splits),
        cover_cells=len(cover),
        dyadic_scale=scale,
        taylor_degree=taylor_degree,
        remainder_order=remainder_order,
        target_radius=target,
        status=status,
        obstruction=obstruction,
        cpu_seconds=cpu,
        coefficient_error_bound=packed,
        remainder_error_bound=None,
        physical_error_bound=None,
        source_tail_error_bound=None,
        physical_local_gate=PHYSICAL_LOCAL_GATE,
        assembly_correction=False,
        diagnostics={
            'rule': 'even-power center Taylor polynomial plus uniform remainder_order coefficient',
            'gauss_nodes_are_bounds': False,
            'local_embedding_reused': True,
            'global_parent_child_matching_required': False,
            'new_field_or_source_runs': 0,
            'remainder_uppers': packed_rems,
            'target_radius_upper': packed_target,
        },
    )
