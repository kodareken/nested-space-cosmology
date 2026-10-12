"""Shared cubic UV channels on the directed massless integrator.

One massless ell=1 forcing jet is integrated per rho box, incoming z and
characteristic sign. Any later (m, ell) coefficient is an algebraic assembly
of those channels. The factorization proof remains
``massive_cubic_identities`` / ``_factor_blocks``; those 74 residuals are not
recomputed here. Unit angular momentum is an algebraic basis, not a source
species. Nothing here covers I, builds C_M, or changes an infinity or Lambda.
"""
import time
import numpy as np

from flint import arb, arb_series, ctx

from . import nsc_ks_ball_geometry as G
from .nsc_ks_ball_operator import slab_ball
from .nsc_ks_cubic_uv_enclosure import (
    CubicForcing,
    CubicGeometry,
    taylor_volterra_cell,
    volterra_cell,
)
from .nsc_ks_profile_identity import profile_identity


class CubicChannelBudgetExceeded(ArithmeticError):
    """The directed path hit its CPU deadline before completion.

    No coefficient is returned. A partial cell sum is not an enclosure.
    """

    def __init__(self, used):
        self.cpu_seconds = float(used)
        super().__init__(
            'cubic channel basis exceeded its CPU budget after %.3f process '
            'seconds; no enclosure returned' % self.cpu_seconds)


class _PreservingContext:
    """Restore flint precision and series cap even when a jet helper fails."""

    def __enter__(self):
        self._before = (ctx.prec, ctx.cap)
        return self

    def __exit__(self, *exc):
        ctx.prec, ctx.cap = self._before
        return False


_SCOPE = (
    'shared massless unit-angular cubic channels; history-minus-reference; '
    'upstream Qz difference is the homogeneous zero representative; '
    'unit angular is not a source species; one z box is not the incoming '
    'interval I; C_M is not constructed')


def _require_sign(sign):
    if type(sign) is not int or sign not in (-1, 1):
        raise ValueError('explicit characteristic sign +/-1 required')


def _require_order(order):
    if type(order) is not int or not 0 <= order <= G.MAX_JET_ORDER - 3:
        raise ValueError('admissible forcing Taylor order required')


def _require_unit_basis(model):
    """Massive current A double-counts m^2 a^2 f, once in J2 and once in JM."""
    if type(model) is not CubicGeometry:
        raise TypeError(
            'massless CubicGeometry required; massive current A double-counts '
            'm^2 a^2 f')
    if not model.angular.is_exact() or not model.angular == 1:
        raise ValueError(
            'exact unit angular basis required; this is not a source species')


def _required_ball(value, name):
    if value is None or isinstance(value, (bool, np.bool_)):
        raise ValueError('explicit ' + name + ' required')
    ball = arb(value)
    if not ball.is_finite():
        raise ArithmeticError('finite ' + name + ' required')
    return ball


def _finite_series(series, order, name):
    if series is None or isinstance(series, (bool, np.bool_)):
        raise ValueError('explicit ' + name + ' jet required')
    values = G._coeffs(series, order)
    if any(not coefficient.is_finite() for coefficient in values):
        raise ArithmeticError('nonfinite %s channel jet' % name)
    return arb_series(values, prec=order + 1)


def _check_budget(deadline, started):
    if deadline is not None and time.process_time() > deadline:
        raise CubicChannelBudgetExceeded(time.process_time() - started)


def _finite_labels(angular, mass):
    if angular is None or mass is None or isinstance(angular, (bool, np.bool_)) or isinstance(mass, (bool, np.bool_)):
        raise ValueError('explicit mass and angular labels required')
    ell = angular if isinstance(angular, arb) else G._binary_arb(angular, 'angular')
    mass_ball = mass if isinstance(mass, arb) else G._binary_arb(mass, 'mass')
    if not ell.is_finite() or ell.contains(0):
        raise ValueError('finite nonzero angular label required')
    if not mass_ball.is_finite():
        raise ValueError('finite mass required')
    return ell, mass_ball


def _deadline_of(value):
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError('CPU deadline must be a process_time timestamp or None')
    if value != value or value in (float('inf'), float('-inf')):
        raise ValueError('finite CPU deadline required')
    return float(value)


def unit_cubic_geometry(family, *, bits=160):
    """Massless ell=1 geometry. The label 1 is an algebraic basis only."""
    with _PreservingContext():
        return CubicGeometry(family, 1, bits=bits)


def split_massless_channels(qz, qzz, current, coupling, axial, sign, order=0):
    """Split one massless ell=1 jet into the J2, J4 and JM forcings.

    f, g, c and h are the qz forcing, qzz forcing, coupling and massless
    current. B = s a^2 c f and A = h - B. h must already be massless: the
    massive current is h + m^2 a^2 f, and JM' carries a^2 f again.
    """
    with _PreservingContext():
        _require_sign(sign)
        if isinstance(qz, arb_series):
            if type(order) is not int or order < 1:
                raise ValueError('positive Taylor order required for series channels')
            with G._JetWork(order):
                source = _finite_series(qz, order, 'qz')
                second = _finite_series(qzz, order, 'qzz')
                rate = _finite_series(coupling, order, 'coupling')
                height = _finite_series(current, order, 'current')
                scale = _finite_series(axial, order, 'axial')
                scale2 = _finite_series(scale * scale, order, 'a^2')
                quartic = _finite_series(
                    arb(sign) * scale2 * rate * source, order, 'B')
                quadratic = _finite_series(height - quartic, order, 'A')
                carried = _finite_series(scale2 * source, order, 'a^2 f')
                zero = _finite_series(
                    arb_series([arb(0)], prec=order + 1), order, 'zero coupling')
                characteristic = _finite_series(
                    arb_series([arb(sign)], prec=order + 1), order, 'sign coupling')

                def pack(cur, coup):
                    return {
                        'qz': source, 'qzz': second, 'current': cur, 'coupling': coup}

                return {
                    'A': quadratic, 'B': quartic, 'f': source, 'g': second,
                    'c': rate, 'h': height, 'axial': scale,
                    'J2': pack(quadratic, zero),
                    'J4': pack(quartic, rate),
                    'JM': pack(carried, characteristic),
                }
        source, second, height, rate, scale = (_required_ball(v, n) for v, n in zip(
            (qz, qzz, current, coupling, axial), ('qz', 'qzz', 'current', 'coupling', 'axial')))
        if not all(value.is_finite() for value in (source, second, height, rate, scale)):
            raise ArithmeticError('nonfinite massless channel forcing')
        quartic = arb(sign) * scale**2 * rate * source
        quadratic = height - quartic
        if not quartic.is_finite() or not quadratic.is_finite():
            raise ArithmeticError('nonfinite massless channel split')
        return {
            'A': quadratic, 'B': quartic, 'f': source, 'g': second,
            'c': rate, 'h': height, 'axial': scale,
            'J2': CubicForcing(source, second, quadratic, arb(0)),
            'J4': CubicForcing(source, second, quartic, rate),
            'JM': CubicForcing(source, second, scale**2 * source, arb(sign)),
        }


def integrate_shared_cell(state, point, ranges, start, end, order):
    """Advance Qz, Qzz, J2, J4 and JM on one decreasing cell.

    The three Volterra updates all read the same pre-cell Qz and Qzz. Chaining
    an updated Qzz into the next channel would integrate that channel twice.
    Positive order calls ``taylor_volterra_cell`` three times. Order 0 uses the
    owner's uniform ``volterra_cell`` the same way; Taylor order starts at 1.
    """
    with _PreservingContext():
        if not isinstance(state, tuple) or len(state) != 5:
            raise TypeError('five channel states Qz, Qzz, J2, J4, JM required')
        if type(order) is not int or order < 0:
            raise ValueError('nonnegative Taylor order required')
        first, second, current2, current4, current_m = (_required_ball(value, 'channel state') for value in state)
        start, end = _required_ball(start, 'cell start'), _required_ball(end, 'cell end')
        fresh = (first, second)
        updated = []
        for name, current in (('J2', current2), ('J4', current4), ('JM', current_m)):
            if order:
                if ranges is None:
                    raise ValueError('cell-wide jets required for a positive Taylor order')
                stepped = taylor_volterra_cell(
                    (fresh[0], fresh[1], current), point[name], ranges[name],
                    start, end, order)
            else:
                stepped = volterra_cell(
                    (fresh[0], fresh[1], current), point[name], start, end)
            updated.append(stepped)
        for index in (0, 2):
            if (not updated[1][0].overlaps(updated[index][0])
                    or not updated[1][1].overlaps(updated[index][1])):
                raise ArithmeticError('shared Qz/Qzz forcing diverged inside one cell')
        return (
            updated[1][0], updated[1][1], updated[0][2], updated[1][2], updated[2][2])


def _forcing_channels(model, rho, z_in, sign, order):
    """One massless ell=1 jet, then the owned axial series. Geometry is not rebuilt."""
    if order:
        series = CubicGeometry.forcing_series(model, rho, z_in, sign, order)
        axial = G.background_series(G._rho_ball(rho), order).a
        return split_massless_channels(
            series['qz'], series['qzz'], series['current'], series['coupling'],
            axial, sign, order)
    values = CubicGeometry.forcing(model, rho, z_in, sign)
    axial = G.background_series(G._rho_ball(rho), 0).a[0]
    return split_massless_channels(
        values.qz, values.qzz, values.current, values.coupling, axial, sign, 0)


def reconstructed_massive_forcing(model, rho, z_in, sign, order, angular, mass):
    """Massive forcing jet assembled from the massless unit basis.

    qz and qzz scale by ell^2. The current is ell^2 A + ell^4 B + m^2 ell^2 a^2 f
    and the coupling is ell^2 c + s m^2. This is the forcing identity, not an
    integrated coefficient and not a tail.
    """
    with _PreservingContext():
        _require_unit_basis(model)
        _require_sign(sign)
        _require_order(order)
        with ctx.workprec(model.bits):
            ell, mass_ball = _finite_labels(angular, mass)
            bundle = _forcing_channels(model, rho, z_in, sign, order)
            ell2 = ell**2
            ell4 = ell2**2
            mass2 = mass_ball**2
            if order:
                with G._JetWork(order):
                    carried = _finite_series(
                        bundle['axial']**2 * bundle['f'], order, 'a^2 f')
                    current = _finite_series(
                        ell2 * bundle['A'] + ell4 * bundle['B']
                        + (mass2 * ell2) * carried, order, 'massive current')
                    coupling = _finite_series(
                        ell2 * bundle['c'] + arb(sign) * mass2, order, 'massive coupling')
                    return {
                        'qz': _finite_series(ell2 * bundle['f'], order, 'massive qz'),
                        'qzz': _finite_series(ell2 * bundle['g'], order, 'massive qzz'),
                        'current': current,
                        'coupling': coupling,
                    }
            carried = bundle['axial']**2 * bundle['f']
            return CubicForcing(
                ell2 * bundle['f'], ell2 * bundle['g'],
                ell2 * bundle['A'] + ell4 * bundle['B'] + mass2 * ell2 * carried,
                ell2 * bundle['c'] + sign * mass2)


def enclose_cubic_channel_basis(model, z_in, sign, rho_up, *, cells=128,
                                 max_depth=16, order=0, cpu_deadline=None):
    """Full-path enclosure of Qz, Qzz, J2, J4, JM and the Sigma N channels.

    Directed balls, the exact upstream rho, flat preparation and the old
    subdivision failure are the historical rules. An unresolved wall raises
    ``SubdivisionNeeded`` and emits nothing. A z ball is covered only as that
    box. ``entire_incoming_interval`` stays false. C_M stays null.

    ``cpu_deadline`` is an absolute ``time.process_time`` stamp. If it passes
    before the path is complete, ``CubicChannelBudgetExceeded`` is raised and
    no partial state is returned.
    """
    started = time.process_time()
    with _PreservingContext():
        _require_unit_basis(model)
        _require_sign(sign)
        if (type(cells) is not int or cells < 1
                or type(max_depth) is not int or max_depth < 0):
            raise ValueError(
                'positive cell count and nonnegative subdivision depth required')
        _require_order(order)
        deadline = _deadline_of(cpu_deadline)
        _check_budget(deadline, started)
        with ctx.workprec(model.bits):
            upper = G._rho_ball(rho_up)
            if not upper > 1 or not upper.is_exact():
                raise ValueError('exact upstream point in the declared slab required')
            if not -G.shift_series(upper, 0)[0] >= arb(model.family.normal_outer):
                raise ValueError(
                    'history must be flat in the unchanged upstream preparation')
            incoming = G._binary_arb(z_in, 'incoming z')
            state = (arb(0), arb(0), arb(0), arb(0), arb(0))
            visited = 0

            def advance(state, start, end, depth):
                nonlocal visited
                _check_budget(deadline, started)
                try:
                    if order:
                        ranges = _forcing_channels(
                            model, slab_ball(end, start), incoming, sign, order)
                        point = _forcing_channels(model, start, incoming, sign, order)
                    else:
                        ranges = _forcing_channels(
                            model, slab_ball(end, start), incoming, sign, 0)
                        point = None
                except G.SubdivisionNeeded:
                    if depth >= max_depth:
                        raise G.SubdivisionNeeded(
                            'maximum cell subdivision reached; no enclosure emitted')
                    mid = (start + end) / 2
                    state = advance(state, start, mid, depth + 1)
                    return advance(state, mid, end, depth + 1)
                visited += 1
                if order:
                    return integrate_shared_cell(
                        state, point, ranges, start, end, order)
                return integrate_shared_cell(state, ranges, None, start, end, 0)

            for index in range(cells):
                start = upper - (upper - 1) * index / cells
                end = upper - (upper - 1) * (index + 1) / cells
                state = advance(state, start, end, 0)
            _check_budget(deadline, started)
            jets = model.jets(arb(1), incoming)
            scale = jets['a']
            base = jets['b']
            base_rho = jets['b_rho']
            base_ref = jets['b_reference_rho']
            base_rhoz = jets['b_rhoz']
            if not scale > 0:
                raise ArithmeticError('surface axial scale is not positive')
            first, second, current2, current4, current_m = state
            quadratic_n = (
                sign * current2 / scale
                + scale**3 / 2 * (base_rho**2 - base_ref**2)
                + sign * scale * base * base_rhoz)
            quartic_n = sign * current4 / scale - 4 * sign * base**2 * first / scale
            mixed_n = sign * current_m / scale - sign * scale * first
            values = (
                first, second, current2, current4, current_m,
                quadratic_n, quartic_n, mixed_n, scale)
            if not all(value.is_finite() for value in values):
                raise ArithmeticError('nonfinite cubic channel basis')
            return {
                'Qz': first, 'Qzz': second, 'J2': current2, 'J4': current4,
                'JM': current_m, 'N2': quadratic_n, 'N4': quartic_n, 'NM': mixed_n,
                'surface_a': scale, 'accepted_cells': visited, 'z_domain': incoming,
                'rho_up': upper, 'sign': sign, 'order': order, 'cells': cells,
                'max_depth': max_depth, 'bits': model.bits,
                'profile_identity': profile_identity(model.family),
                'angular_basis': 1, 'unit_angular_is_source_species': False,
                'massless_A': True, 'coefficient_scope': _SCOPE,
                'C_M': None, 'higher_uv_remainder_C_M': None,
                'higher_uv_remainder_bound': None, 'uniform_C4_on_I': None,
                'complete_UV_tail': None, 'entire_incoming_interval': False,
                'z_box_is_full_incoming_interval': False,
                'physical_local_gate': 'OPEN',
                'cpu_seconds': time.process_time() - started,
            }


def assemble_cubic_channels(basis, angular, mass):
    """Assemble qz, qzz, J3 and N from one integrated characteristic basis.

    qz = ell^2 Qz, qzz = ell^2 Qzz, and
    J3 = ell^2 J2 + ell^4 J4 + m^2 ell^2 JM. The same powers multiply N2, N4
    and NM. Evenness in ell is algebraic. It does not identify the two
    characteristic signs. This does not modify a residual assembly.
    """
    with _PreservingContext():
        if basis.get('massless_A') is not True:
            raise ValueError('A must come from the massless unit-angular basis')
        if basis.get('C_M') is not None or basis.get('higher_uv_remainder_C_M') is not None:
            raise ValueError('channel basis must leave C_M unconstructed')
        if basis.get('entire_incoming_interval'):
            raise ValueError('one channel box is not the incoming interval')
        bits = basis.get('bits')
        if type(bits) is not int or bits < 80:
            raise ValueError('channel basis precision is missing')
        sign = basis.get('sign')
        _require_sign(sign)
        with ctx.workprec(bits):
            basis = dict(basis)
            for key in ('Qz', 'Qzz', 'J2', 'J4', 'JM', 'N2', 'N4', 'NM', 'surface_a'):
                basis[key] = _required_ball(basis.get(key), key)
            ell, mass_ball = _finite_labels(angular, mass)
            ell2 = ell**2
            ell4 = ell2**2
            mass2 = mass_ball**2
            qz = ell2 * basis['Qz']
            qzz = ell2 * basis['Qzz']
            current = (
                ell2 * basis['J2'] + ell4 * basis['J4'] + mass2 * ell2 * basis['JM'])
            bracket = (
                ell2 * basis['N2'] + ell4 * basis['N4'] + mass2 * ell2 * basis['NM'])
            correction = -sign * basis['surface_a'] * mass2 * qz
            values = (qz, qzz, current, bracket, correction)
            if not all(value.is_finite() for value in values):
                raise ArithmeticError('nonfinite assembled cubic coefficient')
            return {
                'q_z': qz, 'q_zz': qzz, 'J3': current, 'N_bracket3': bracket,
                'surface_mass_correction': correction, 'angular': ell,
                'mass': mass_ball, 'sign': sign,
                'C_M': None, 'higher_uv_remainder_C_M': None,
                'higher_uv_remainder_bound': None, 'uniform_C4_on_I': None,
                'complete_UV_tail': None, 'entire_incoming_interval': False,
                'physical_local_gate': 'OPEN',
                'massless_A': True, 'coefficient_scope': _SCOPE,
            }
