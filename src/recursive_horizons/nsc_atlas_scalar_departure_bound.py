"""Uniform canonical scalar differences from certified geometric differences.

The exterior equation is mixed: X_Ux=-K X. The cap comparison follows
actual characteristics on one common z collar. All inputs are uniform
enclosures; no field samples, evolution, or amplitude derivatives enter.
"""
from fractions import Fraction as Q

from flint import arb, ctx

from .nsc_atlas_geometric_departure_bound import (
    frozen_context, _stored, _u, _l, _round_upper,
)
from .nsc_atlas_coupled_tube_bound import wall_difference_bounds
from .nsc_atlas_onset_geometry_bound import _stored as interval
from .nsc_atlas_onset_endpoint_bound import prefix_kernel_integral
from .nsc_atlas_reference_gauge_bound import _ball


def exterior_comparison(*, duration, width, kernel, source_X, source_Q,
                        wall_P, kernel_pointwise, kernel_integrated, reference_X,
                        budget):
    """Eliminate DX<=source_X+width*DP in the mixed-equation comparison."""
    budget.check(512)
    raw = (duration, width, kernel, source_X, source_Q, wall_P,
           kernel_pointwise, kernel_integrated, reference_X)
    if any(isinstance(v, (float, bool)) for v in raw):
        raise ValueError('require exact rational comparison inputs')
    dt, w, k, sx, sq, pw, ks, ki, x0 = map(Q, raw)
    if any(v < 0 for v in (dt, w, k, sx, sq, pw, ks, ki, x0)):
        raise ValueError('require nonnegative comparison inputs')
    loop = dt*w*k
    if loop >= 1:
        raise ValueError('exterior comparison loop reaches one')
    p = (pw+dt*(k*sx+ks*x0))/(1-loop)
    x = sx+w*p
    q = sq+w*k*x+ki*x0
    budget.check()
    return dict(loop=loop, X=x, P=p, Q=q)


def reference_jets(parameters, gauge, z_star):
    """Reference spatial jets from the mixed equation and the quiet cap."""
    wall, wedge, inp = (gauge[k] for k in ('wall', 'vacuum_wedge', 'inputs'))
    g0, g1, g2 = (_ball(Q(parameters[k])) for k in ('g', 'g_x', 'g_xx'))
    xl, xr = (Q(inp[k]) for k in ('left', 'right'))
    rl, rh = (Q(inp[k]) for k in ('radius_lower', 'radius_upper'))
    x, z, a = _ball(xl).union(_ball(xr)), _ball(rl).union(_ball(rh)), _ball(Q(4, 5))
    r, rx = interval(wedge['r']), interval(wedge['r_x'])
    k = _stored(wedge['K'])
    # C0=(x/r) exp((x-r)/a), K0=C0*a/(4*r^3) in the fixed Ubar gauge.
    kx = _u(k*abs(1/x+1/a-rx*(1/a+4/r)))
    b = _ball(prefix_kernel_integral(xl, rh))
    q0 = _u(g1+b*g0)
    w2 = _u(g2+_stored(wall['W_x'])*k*g0
            +_stored(wall['W'])*(kx*g0+k*q0))
    sx, v, bp, bm = (interval(wall[key]) for key in ('s_x', 'v', 'beta_plus', 'beta_minus'))
    az, d = 1-z*z, v-bp
    bpz = -2/(5*z**3)-_ball(Q(5, 4))
    bmz = bpz+2
    vz = bp*bpz/v-a/(2*v*z*z)
    zx = v*sx
    sxx = 2*a*d/(a-x)**2+2*x*(vz-bpz)*zx/(a-x)
    qproper = (bm+v)/az
    qproper_z = (bmz+vz)/az+2*z*(bm+v)/az**2
    qxx = _u(abs(sxx*qproper+sx*qproper_z*zx))
    qmin = _stored(wall['T_x_plus_Z_x'], 'lower')
    length = _u(_ball(rh).atanh()-_ball(Q(z_star)))
    if not qmin > 0 or not length > 0:
        raise ValueError('positive reference arrival rate/collar required')
    fp, fpp = _u(q0/qmin), _u(w2/qmin**2+q0*qxx/qmin**3)
    return dict(Q=q0, wall_xx=w2, K_x=kx, arrival_xx=qxx,
                F_prime=fp, F_second=fpp, chi_eta=_u(length*fp),
                E_plus_eta=_u(2*length*fpp), E_plus=_u(2*fp))


def scalar_departure_bound(parameters, source, gauge, tube, geometry, *, budget):
    if not callable(getattr(budget, 'check', None)):
        raise ValueError('require a shared cumulative budget')
    budget.check(8192)
    with ctx.workprec(160):
        c = frozen_context(source, gauge, tube)
        if (Q(parameters['alpha_max']) != c['al']
                or Q(geometry['alpha_max']) != c['al']
                or Q(geometry['duration']) != Q(tube['duration'])):
            raise ValueError('geometry/source amplitude or duration mismatch')
        for key in ('left', 'right'):
            if Q(parameters[key]) != Q(gauge['inputs'][key]):
                raise ValueError('source/reference domain mismatch')
        g = {k:_ball(Q(v)) for k,v in geometry['coefficients'].items()}
        t = {k:_stored(v) for k,v in tube['constants'].items()}
        ref = reference_jets(parameters, gauge, tube['z_star'])
        g0, g1 = (_ball(Q(parameters[k])) for k in ('g', 'g_x'))
        al, dt, W, k = _ball(c['al']), c['dt'], c['W'], c['k']
        rlo, rhi, Rlo, rstar = c['rlo'], c['rhi'], c['Rlo'], t['cap_inner_radius']
        A0min, A0max, A0wall = _l(1-rhi*rhi), _u(1-rstar*rstar), _u(1-Rlo*Rlo)
        length = _l(Rlo.atanh()-_ball(Q(tube['z_star'])))
        if not A0min > 0 or not length > 0:
            raise ValueError('positive sourced cap/collar required')
        # Existing positive formulas have explicit alpha factors; their
        # coefficient bounds use the frozen alpha_max tube throughout.
        gamma_c = _u(t['cap_beta_z_upper']/al)
        reaction_c = _u(t['cap_reaction_upper']/al)
        beta_c = _u(t['cap_mu_upper']*(t['cap_lapse_upper']/al
                    +2*c['N']/(rstar*A0min)))
        tx, beta = t['wall_clock_rate_upper'], t['cap_beta_upper']
        minus_loop = _u(dt*tx*t['cap_beta_z_upper'])
        if not minus_loop < 1:
            raise ValueError('outgoing cap comparison loop reaches one')
        em = _u(dt*tx*A0max*reaction_c/(1-minus_loop))
        fw = _u(t['X_upper']/Rlo)
        jc = _u(c['Sw']*rhi*fw**3/18)
        pw = _u((A0wall*em/(c['Hmin']*A0min)+jc)/(2*c['Upmin']))
        sx = _u(g0**3/(36*rlo**2))
        sq = _u((g0/rlo)**2*g1/12+(g0/rlo)**3*c['vmax']/18)
        def kernel(jet, span):
            return _u(c['vmax']*jet/rlo**2+c['pabs']*span*g['r_x']/rlo**2
                +c['Cmax']*span*g['log_C']/(4*rlo**2)+2*k*span*g['radius']/rlo)
        ks, ki = kernel(g['delta_r_U'], arb(1)), kernel(g['J_r'], W)
        feedback = dict(duration=Q(tube['duration']), width=_round_upper(W),
            kernel=_round_upper(k), source_X=_round_upper(sx), source_Q=_round_upper(sq),
            wall_P=_round_upper(pw), kernel_pointwise=_round_upper(ks),
            kernel_integrated=_round_upper(ki), reference_X=Q(parameters['g']))
        ext = exterior_comparison(**feedback, budget=budget)
        dx, dp, dq = (_ball(ext[key]) for key in ('X', 'P', 'Q'))
        budget.check(8192)
        wd = wall_difference_bounds(eps_r=g['radius'], eps_L=g['log_C'],
            Verror=g['r_x'], Lxerror=g['L_x'], alpha=arb(1), N=c['N'],
            mu_error=c['mu_coefficient'], Rlo=Rlo, Rhi=rhi,
            R0lo=c['R0lo'], R0hi=c['R0hi'], Amin=c['Amin'], Hmin=c['Hmin'],
            Upmin=c['Upmin'], Upmax=c['Upmax'], Vpmin=c['Vpmin'], Vpmax=c['Vpmax'],
            Wxmin=c['Wxlo'], rx0max=c['rx0hi'], Sw=c['Sw'], root_change=c['drt'],
            Sw_change=c['dSw'], u_wall=c['u_wall'], nc=c['nc'], beta_bar_max=c['bbar'])
        # |Delta(2*A0*Vplus/H)|/alpha, retaining actual wall radius.
        aE = _u(2*A0wall*c['Vpmax']/c['Hmin'])
        daE = _u(2*(2*rhi*g['radius']*c['Vpmax']/c['Hmin']
                  +A0wall*wd['Vdot']/c['Hmin']
                  +A0wall*c['Vpmax']*wd['H']/c['Hmin']**2))
        ew = _u(aE*dq+daE*ref['Q']+A0wall*jc/c['Hmin'])
        nu_c = _u(g['delta_T']/dt)
        rx_c = _u(g['r_x']+c['Wxhi']*g['delta_r_U'])
        ref_rx = _stored(gauge['wall']['R_x'])
        ref_zx = _stored(gauge['wall']['Z_x'])
        ref_tx = _stored(gauge['wall']['T_x'])
        z_c = _u(g['radius']/A0min)
        zx_c = _u(rx_c/A0min+ref_rx*2*rhi*g['radius']/A0min**2)
        velocity = _u(nu_c+tx*beta_c)
        chi_speed = _u(zx_c/length+ref_zx*z_c/length**2)
        plus_speed = _u(chi_speed+velocity/length+ref_tx*z_c/length**2)
        ep = _u(ew+dt*(plus_speed*ref['E_plus_eta']
                    +tx*(gamma_c*_ball(Q(tube['box']['E_plus']))+A0max*reaction_c)))
        chi = _u(dx+dt*(chi_speed*ref['chi_eta']
                    +velocity*ref['E_plus']/2+tx*beta*(ep+em)/2))
        coefficients = {key:_round_upper(value) for key,value in
            dict(X=dx,P=dp,Q=dq,E_minus=em,E_plus=ep,chi=chi,wall_P=pw,wall_E_plus=ew).items()}
        comparison = {key:_round_upper(value) for key,value in dict(
            minus_loop=minus_loop,beta=beta_c,beta_z=gamma_c,reaction=reaction_c,J=jc,
            incoming_multiplier=aE,incoming_multiplier_difference=daE,
            clock_rate=nu_c,wall_radius_rate=rx_c,wall_z=z_c,wall_z_rate=zx_c,
            cap_velocity=velocity,chi_speed=chi_speed,plus_speed=plus_speed).items()}
        budget.check(2048)
        return dict(alpha_max=c['al'],duration=Q(tube['duration']),z_star=Q(tube['z_star']),
            coefficients=coefficients,maximum_amplitude_bounds={k:c['al']*v for k,v in coefficients.items()},
            exterior_inputs=feedback,exterior_loop=ext['loop'],comparison=comparison,
            reference={k:_round_upper(v) for k,v in ref.items()},precision_bits=160,
            scope='uniform canonical normalized scalar differences from alpha=0 on the original '
                  'source wedge and common z collar; no quadratic remainder, physical-field '
                  'alpha derivative, later arrival, or stored-trajectory error')
