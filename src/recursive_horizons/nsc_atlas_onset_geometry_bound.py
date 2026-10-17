"""Whole-panel enclosures of first geometric coefficients on the quiet onset.

The caller authenticates the source/reference receipts. This is integration
of continuum comparison bounds, not an evolution of sampled field values.
All quantities are dimensionless; x is the original affine source label.
"""
from fractions import Fraction as Q

from flint import arb, ctx

from .nsc_atlas_onset_endpoint_bound import prefix_kernel_integral
from .nsc_atlas_source_cell_bound import _ball, _coordinate, _pulse_jets


def dyadic(value):
    if not value.is_finite():
        raise ValueError('nonfinite geometric enclosure')
    mantissa, exponent = value.man_exp()
    return Q(int(mantissa))*Q(2)**int(exponent)


def bounds(value):
    return dict(lower=dyadic(value.lower()), upper=dyadic(value.upper()))


def _stored(row):
    def end(name):
        value = row[name]
        return Q(value['mantissa'])*Q(2)**value['exponent']
    return _ball(end('lower')).union(_ball(end('upper')))


def _range(lo, hi):
    if lo > hi:
        raise ValueError('reversed interval')
    return lo.union(hi)


def _positive(value, name):
    if not value.is_finite() or not value > 0:
        raise ValueError('positive finite '+name+' required')
    return value


def _absmax(value):
    return abs(value).upper()


def moment_panel(first, second, density, h):
    """Enclose a primitive and its integral, including every interior point.

    For a constant density interval on a panel, F'=density and G'=F.
    The h^2/2 term bounds the triangular integral, not a trapezoid estimate.
    """
    duration = _range(arb(0), h.upper())
    triangle = _range(arb(0), (h*h/2).upper())
    return (first+duration*density, second+duration*first+triangle*density,
            first+h*density, second+h*first+h*h*density/2)


def mass_terms(z, v, bp, sx, w, wx):
    p, j, d = w/z, wx/sx, v-bp
    fp = 1-_ball(Q(4, 5))/z
    lp = -d*j*j/2-fp*j*p/2+v*fp*p*p/4
    ld = -z*j*j+2*z*(bp+z)*j*p-v*z*(bp+z)*p*p-z*z*v*p*p/2
    cp, cd = 1/(6*z)-1/(10*z*z), z/4-1/(10*z*z)
    qp, qd = d*w*j/6+cp*w*w, z*w*j/3+cd*w*w
    return lp, ld, qp, qd


def radius_rate(c, ip, idiff, w):
    """Derivative of E=e/v; frame endpoint w*j has cancelled exactly."""
    z, v, bp, sx = (c[k] for k in ('z', 'v', 'bp', 'sx'))
    return sx*((ip/z+bp*idiff/(2*z*z))/(v*v)+c['Ce']*w*w)


def cap_clock_rate(c, e, rho, difference_prefix, w, jcap):
    """Both cap mass/lapse and local w*j cancellations are explicit."""
    z, bm, sx, az = (c[k] for k in ('z', 'bm', 'sx', 'Az'))
    cd = z/4-1/(10*z*z)
    dbm_reduced = c['bm_z']*rho+(difference_prefix+cd*w*w)/(2*z*z)
    return sx*(dbm_reduced/az+bm*w*w/(12*z*z*az)
               +2*bm*z*e/(az*az)+bm*jcap/(az*az))


def geometry_constants(gauge, left, right):
    """Use only the proved common reference hull; no new gauge inversion."""
    inp, wall, wedge = (gauge[k] for k in ('inputs', 'wall', 'vacuum_wedge'))
    zl, zh = (_ball(Q(inp[k])) for k in ('radius_lower', 'radius_upper'))
    if not 0 < left < right < Q(inp['radius_lower']) <= Q(inp['radius_upper']) < Q(1, 2):
        raise ValueError('require the admitted expanding quiet-collar radius strip')
    z, a, ell = zl.union(zh), _ball(Q(4, 5)), _ball(right-left)
    sx, v, bp, bm, qx = (_stored(wall[k]) for k in
                         ('s_x', 'v', 'beta_plus', 'beta_minus', 'T_x_plus_Z_x'))
    d, az = v-bp, 1-z*z
    for name, value in dict(sx=sx, v=v, bp=bp, bm=bm, d=d, Az=az, Qx=qx).items():
        _positive(value, name)
    bp_z = -2/(5*z**3)-_ball(Q(5, 4))
    bm_z = bp_z+2
    vz = bp*bp_z/v-a/(2*v*z*z)
    dz = vz-bp_z
    cp, cd = 1/(6*z)-1/(10*z*z), z/4-1/(10*z*z)
    ce = vz/(12*z*v)+(cp/z+bp*cd/(2*z*z))/(v*v)+1/(12*z*z)
    bstar = _ball(prefix_kernel_integral(left, Q(inp['radius_upper'])))
    q = ell*bstar
    if not bstar < 1 or not q < 1:
        raise ValueError('require positive-source comparison with Bstar<1 and q<1')
    rstar = _positive(zl-qx.upper()*ell, 'quiet collar radius')
    amin, amax, hmax = 1-zh*zh, 1-rstar*rstar, 1-3*rstar*rstar
    _positive(1-3*zh*zh, 'cap kinetic sign margin')
    ckin = qx.upper()**2*zh*hmax*amax/(2*zl)
    cap21 = (ckin*2/(qx.lower()**2*rstar*rstar*amin*amin)).upper()
    cap01 = (ckin/rstar**4).upper()
    cap00 = (amax*qx.upper()/zl).upper()
    r, ru, W = _stored(wedge['r']), _stored(wedge['r_U']), _stored(wall['W'])
    _positive(r, 'exterior radius'); _positive(ru, 'exterior r_U')
    _positive(W, 'exterior width')
    kmax = _stored(wedge['K']).upper()
    _positive(kmax, 'kernel upper bound')
    cu = kmax*ell/r.lower()+ru.upper()/r.lower()**2
    bulk = (W.upper()*r.upper()*cu*cu/(2*ru.lower())).upper()
    # Separate appearances of current/past geometry share this hull, not a
    # value. Interval operations safely enclose their independent products.
    fp = 1-a/z
    hjp = -z*fp-2*d*z*(bp+z)
    hcp = z*v*fp/2+d*(v*z*(bp+z)+z*z*v/2)
    hlip = (zh*zh*_absmax((z*dz-d)/(z*z))*(v*sx).upper()).upper()
    eta21 = (hlip/sx.lower()).upper()
    etag2 = (_absmax(hjp)/(2*zl)+_absmax(2*z*cp-d*cd)).upper()
    eta00 = (_absmax(hcp)*sx.upper()/(zl*zl)).upper()
    return dict(z=z, v=v, bp=bp, bm=bm, sx=sx, d=d, Az=az, a=a,
        vz=vz, dz=dz, bm_z=bm_z, Ce=ce, Bstar=bstar, q=q,
        rstar=rstar, Qx=qx, bulk=bulk, cap21=cap21, cap01=cap01, cap00=cap00,
        eta21=eta21, etag2=etag2, eta00=eta00)


STATE_KEYS = ('Iplus', 'ID', 'I0', 'I1', 'I2', 'M00', 'M01', 'M20', 'M21', 'E', 'eta', 'T1')


def _step(state, c, x, h, g, gx):
    """A triangular panel: moments, mass prefixes, radius, then clocks."""
    span = _range(arb(0), h.upper())
    # The analytic monotonicity proof allows intersection with nonnegative
    # reals even when rounded ball arithmetic gives a tiny negative floor.
    w = _range(max(arb(0), ((1-c['q'])*g).lower()), max(arb(0), g.upper()))
    wx = _range(max(arb(0), (gx-c['Bstar']*g).lower()), max(arb(0), gx.upper()))
    lp, ld, qp, qd = mass_terms(c['z'], c['v'], c['bp'], c['sx'], w, wx)
    ip = state['Iplus']+span*c['sx']*lp
    idiff = state['ID']+span*c['sx']*ld
    m00, m01, n00, n01 = moment_panel(state['M00'], state['M01'], g*g, h)
    m20, m21, n20, n21 = moment_panel(state['M20'], state['M21'], gx*gx, h)
    ps, psx = g/x, gx/x-g/(x*x)
    i0 = state['I0']+span*ps*ps
    i1, i2, ni1, ni2 = moment_panel(state['I1'], state['I2'], x*psx*psx, h)
    erate = radius_rate(c, ip, idiff, w)
    E = state['E']+span*erate
    e = c['v']*E
    rho = e+w*w/(12*c['z'])
    jcap = _range(-(c['cap21']*m21+c['cap01']*m01).upper(),
                  (c['cap00']*m00).upper())
    # Cancellations are bounded before forming the clock rates.
    jeta_abs = (c['eta21']*m21+c['etag2']*g*g+c['eta00']*m00).upper()
    jeta = _range(-jeta_abs, jeta_abs)
    source_zeta = (-i2/2-i0/6)/x-ps*ps/6-i0/(6*(c['a']-x))
    bulk = _range(arb(0), (c['bulk']*g*g).upper())
    zeta = source_zeta+bulk+w*w/(4*c['z']**2)
    etarate = c['sx']/(2*c['a']*c['d'])*(zeta+c['dz']*rho/c['d']
                   +jeta/(2*c['z']**2*c['v']*c['d']))
    trate = cap_clock_rate(c, e, rho, idiff, w, jcap)
    return dict(Iplus=state['Iplus']+h*c['sx']*lp, ID=state['ID']+h*c['sx']*ld,
        I0=state['I0']+h*ps*ps, I1=ni1, I2=ni2,
        M00=n00, M01=n01, M20=n20, M21=n21,
        E=state['E']+h*erate, eta=state['eta']+h*etarate, T1=state['T1']+h*trate)


def onset_geometry_bound(gauge, geometry, endpoint, *, budget, panels=128):
    if isinstance(panels, bool) or not isinstance(panels, int) or not 1 <= panels <= 128:
        raise ValueError('require one fixed partition of 1 to 128 panels')
    if not callable(getattr(budget, 'check', None)):
        raise ValueError('require a shared cumulative budget')
    budget.check(8192)
    left, right = (Q(gauge['inputs'][k]) for k in ('left', 'right'))
    if not Q(geometry['U0']) < 0 or not Q(geometry['Va']) < Q(geometry['Vb']) < 0:
        raise ValueError('original source signs required')
    with ctx.workprec(160):
        L, R = _ball(left), _ball(right)
        zl, zxl, _ = _coordinate(L, geometry)
        zr, _, _ = _coordinate(R, geometry)
        if not zl < -1 or not -1 < zr or not zr < -_ball(Q(1, 2)) or not zxl > 1:
            raise ValueError('require quiet left edge and strong monotone first-half source')
        c = geometry_constants(gauge, left, right)
        state = {k:arb(0) for k in STATE_KEYS}
        h = _ball((right-left)/panels)
        for index in range(panels):
            budget.check(1024)
            lo, hi = (_ball(left+Q(j, panels)*(right-left)) for j in (index, index+1))
            g, gx, _ = _pulse_jets(lo, hi, geometry)
            state = _step(state, c, lo.union(hi), h, g, gx)
            if not all(v.is_finite() for v in state.values()):
                raise ValueError('nonfinite panel enclosure')
        budget.check(1024)
        def value(name):
            row = endpoint[name]
            return _ball(Q(row['lower'])).union(_ball(Q(row['upper'])))
        w, wx = value('wall_value'), value('wall_derivative')
        _, _, qp, qd = mass_terms(c['z'], c['v'], c['bp'], c['sx'], w, wx)
        e = c['v']*state['E']
        rho = e+w*w/(12*c['z'])
        shift = 2*c['a']*c['d']*state['eta']
        outputs = dict(mu=state['Iplus']+qp, mass_difference=state['ID']+qd,
            e_proper_time=e, rho_proper_time=rho, eta=state['eta'], T1_proper_time=state['T1'],
            radius_fixed_source=e+c['v']*shift,
            clock_fixed_source=state['T1']+c['bm']/c['Az']*shift)
        if not all(v.is_finite() for v in outputs.values()):
            raise ValueError('nonfinite endpoint enclosure')
        budget.check()
        return dict(outputs={k:bounds(v) for k, v in outputs.items()},
            constants={k:bounds(v) for k, v in c.items()},
            moments={k:bounds(state[k]) for k in ('Iplus', 'ID', 'I0', 'I1', 'I2', 'M00', 'M01', 'M20', 'M21', 'E')},
            source_conditions=dict(compact_coordinate_right=bounds(zr), coordinate_derivative_left=bounds(zxl)),
            panels=panels, precision_bits=160,
            scope='first amplitude coefficients at the original fixed source label; '
                  'dimensionless Einstein radius and centre-anchored polar clock; '
                  'no evaluated quadratic Taylor constant or central arrival')
