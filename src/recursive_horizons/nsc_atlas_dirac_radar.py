"""Exact wall-clock null/radar chart of the quiet cap--wall reference.

Units ell=1. Each optical null coordinate is labelled by the proper time
of its wall intersection. No spatial cutoff profile or box boundary is
introduced. The finite null square covers a causal diamond, not a whole
global Cauchy surface or a quantum state.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import lambertw

from .nsc_atlas_cap_wall_response import reference_coefficients
from .nsc_atlas_causal_reference import _integral

RMIN, RMAX = .45, .55
MASS = .4


def tortoise(radius, region):
    r = np.asarray(radius, dtype=float)
    if not np.all(np.isfinite(r)) or np.any(r <= 0):
        raise ValueError('positive finite areal radius required')
    if region == 'cap':
        if np.any(r >= 1):
            raise ValueError('this cap optical chart lies below its static horizon')
        return np.arctanh(r)
    if region == 'exterior':
        if np.any(r >= 2*MASS):
            raise ValueError('this exterior optical chart is on the white-hole interior branch')
        return r+2*MASS*np.log1p(-r/(2*MASS))
    raise ValueError('region must be cap or exterior')


def causal_domain_bounds():
    """Conservative analytic inequalities on the entire null square.

    beta_minus and beta_plus decrease on this fixed wall strip;
    sqrt(-f_plus)<=v<beta_minus there. Bounds use those inequalities,
    not minima found from a sampled chart.
    """
    lo, hi = reference_coefficients(RMIN), reference_coefficients(RMAX)
    bm_max, bm_min = float(lo['beta_minus']), float(hi['beta_minus'])
    bp_abs = max(abs(float(lo['beta_plus'])), abs(float(hi['beta_plus'])))
    fext_abs_min = 2*MASS/RMAX-1
    fcap_wall_min = 1-RMAX**2
    duration_upper = (RMAX-RMIN)/math.sqrt(fext_abs_min)
    cap_z_min = float(tortoise(RMIN,'cap'))-duration_upper/(2*bm_min)
    cap_z_max = float(tortoise(RMIN,'cap'))+duration_upper*bm_max/fcap_wall_min
    cap_r_min, cap_r_max = math.tanh(cap_z_min), math.tanh(cap_z_max)
    cap_C_min = (1-cap_r_max**2)/(2*bm_max**2)
    ext_C_min = fext_abs_min/(bp_abs+bm_max)**2
    if min(cap_r_min, cap_C_min, ext_C_min) <= 0:
        raise ValueError('the declared source strip has no positive radar-chart bound')
    return dict(wall_duration_upper=duration_upper, cap_radius_lower=cap_r_min,
        cap_radius_upper=cap_r_max, exterior_radius_lower=RMIN,
        exterior_radius_upper=RMAX, conformal_factor_lower=min(cap_C_min,ext_C_min),
        lapse_lower=math.sqrt(min(cap_C_min,ext_C_min)))


def _radii(value):
    raw = np.asarray(value)
    if raw.dtype.kind not in 'iuf':
        raise ValueError('wall-radius labels must be real data')
    r = np.asarray(raw,dtype=float)
    if not np.all(np.isfinite(r)) or np.any((r < RMIN)|(r > RMAX)):
        raise ValueError('wall-radius labels must lie in [.45,.55]')
    return r


def direct_wall_maps(radius, order=64):
    """Independent existing Gauss quadrature for (tau,u_cap,v_cap,u_ext,v_ext)."""
    r = float(_radii(radius))
    initial = np.array([0.,-tortoise(RMIN,'cap'),tortoise(RMIN,'cap'),
                        -tortoise(RMIN,'exterior'),tortoise(RMIN,'exterior')])
    def rate(x, index):
        g = reference_coefficients(x)
        v = g['v']
        if index == 0:
            return 1/v
        side = 'minus' if index in (1,2) else 'plus'
        direction = 1 if index in (1,3) else -1
        return 1/(v*(g['beta_'+side]+direction*v))
    return initial+np.array([_integral(lambda x:rate(x,index),RMIN,r,order)
                            for index in range(5)])


class WallRadarChart:
    """U,V are proper wall null labels; T=(U+V)/2,x=(V-U)/2.

    The physical chart uses cap V<=U, exterior V>=U. Explicit single-region
    overrides are geometric continuation controls, not another attachment.
    """
    def __init__(self):
        self.bounds = causal_domain_bounds()
        initial = np.array([0.,-tortoise(RMIN,'cap'),tortoise(RMIN,'cap'),
                            -tortoise(RMIN,'exterior'),tortoise(RMIN,'exterior')])
        def rhs(r, values):
            g = {k:float(v) for k,v in reference_coefficients(r).items()}
            v, bm, bp = g['v'],g['beta_minus'],g['beta_plus']
            return [1/v,1/(v*(bm+v)),1/(v*(bm-v)),1/(v*(bp+v)),1/(v*(bp-v))]
        solution = solve_ivp(rhs,(RMIN,RMAX),initial,method='DOP853',rtol=2e-13,
                             atol=3e-15,max_step=(RMAX-RMIN)/4,dense_output=True)
        if not solution.success:
            raise RuntimeError(solution.message)
        self.solution = solution.sol
        self.duration = float(solution.y[0,-1])
        self.wall_rhs_calls = solution.nfev

    def maps(self, radius):
        r = _radii(radius)
        return self.solution(r.ravel()).reshape((5,)+r.shape)

    def radius_at_clock(self, clock):
        t = np.asarray(clock,dtype=float)
        if not np.all(np.isfinite(t)) or np.any((t < 0)|(t > self.duration)):
            raise ValueError('null clock labels must lie in [0,wall duration]')
        unique, inverse = np.unique(t.ravel(),return_inverse=True)
        radii = []
        for value in unique:
            if value == 0:
                r = RMIN
            elif value == self.duration:
                r = RMAX
            else:
                r = brentq(lambda y:float(self.solution(y)[0])-value,RMIN,RMAX,
                            xtol=5e-15,rtol=8*np.finfo(float).eps)
            radii.append(r)
        return np.asarray(radii)[inverse].reshape(t.shape)

    def fields(self, clock_u, clock_v, *, region='nested', angular_k=1):
        U,V = np.broadcast_arrays(np.asarray(clock_u,dtype=float),np.asarray(clock_v,dtype=float))
        clocks = np.concatenate((U.ravel(),V.ravel()))
        r = self.radius_at_clock(clocks)
        ru,rv = r[:U.size].reshape(U.shape),r[U.size:].reshape(V.shape)
        return self.from_wall_radii(ru,rv,region=region,angular_k=angular_k)

    def from_wall_radii(self, radius_u, radius_v, *, region='nested', angular_k=1):
        ru,rv = np.broadcast_arrays(_radii(radius_u),_radii(radius_v))
        if isinstance(angular_k,bool) or not isinstance(angular_k,int) or not angular_k:
            raise ValueError('a nonzero integer angular Dirac eigenvalue is required')
        if region not in ('nested','cap','exterior'):
            raise ValueError('region must be nested, cap or exterior')
        cap = ru >= rv if region == 'nested' else np.full(ru.shape,region == 'cap')
        mu,mv = self.maps(ru),self.maps(rv)
        gu,gv = reference_coefficients(ru),reference_coefficients(rv)
        bu,bv = [np.where(cap,g['beta_minus'],g['beta_plus']) for g in (gu,gv)]
        bru,brv = [np.where(cap,g['beta_minus_x'],g['beta_plus_x']) for g in (gu,gv)]
        vu,vv = gu['v'],gv['v']
        pu,qv = 1/(bu+vu),1/(bv-vv)
        optical_u = np.where(cap,mu[1],mu[3])
        optical_v = np.where(cap,mv[2],mv[4])
        z = (optical_v-optical_u)/2
        r = np.empty_like(z)
        r[cap] = np.tanh(z[cap])
        if np.any(~cap):
            inverse = lambertw(-np.exp(z[~cap]/(2*MASS)-1),k=0)
            if np.max(np.abs(inverse.imag)) > 1e-13:
                raise ValueError('exterior tortoise inverse left its real white-hole branch')
            r[~cap] = 2*MASS*(1+inverse.real)
        if np.any(r <= 0):
            raise ValueError('radar chart reached nonpositive areal radius')
        f = np.where(cap,1-r*r,1-2*MASS/r)
        fp = np.where(cap,-2*r,2*MASS/r**2)
        fpp = np.where(cap,-2.,-4*MASS/r**3)
        C = f*pu*qv
        if np.any(C <= 0) or not np.all(np.isfinite(C)):
            raise ValueError('radar chart lost its future conformal lapse')
        rU,rV = -f*pu/2,f*qv/2
        p_log = -(vu*bru+gu['radial_acceleration'])/(bu+vu)
        q_log = -(vv*brv-gv['radial_acceleration'])/(bv-vv)
        LU,LV = -fp*pu/2+p_log,fp*qv/2+q_log
        N = np.sqrt(C)
        rt,rx = rU+rV,rV-rU
        Nt,Nx = N*(LU+LV)/2,N*(LV-LU)/2
        W = angular_k*N/r
        return dict(clock_u=mu[0],clock_v=mv[0],time=(mu[0]+mv[0])/2,
            radar_position=(mv[0]-mu[0])/2,is_cap=cap,radius=r,radius_t=rt,radius_x=rx,
            conformal_factor=C,lapse=N,radial_length=N,lapse_t=Nt,lapse_x=Nx,
            log_C_u=LU,log_C_v=LV,radius_uv=-C*fp/4,log_C_uv=-C*fpp/4,
            static_f=f,radial_curvature=-fpp,potential=W,
            potential_t=angular_k*(Nt/r-N*rt/r**2),
            potential_x=angular_k*(Nx/r-N*rx/r**2))


def reference_radar_report(chart):
    """One predeclared finite geometry check and wall-radar clock protocol."""
    labels = np.linspace(0.,chart.duration,9)
    U,V = np.meshgrid(labels,labels,indexing='ij')
    row = chart.fields(U,V)
    r,C = row['radius'],row['conformal_factor']
    mass = r*(1-(row['radius_x']**2-row['radius_t']**2)/C)/2
    expected_mass = np.where(row['is_cap'],r**3/2,MASS)
    R2 = 4*row['log_C_uv']/C
    box_r = -4*row['radius_uv']/C
    R4 = R2+2*(1-row['static_f'])/r**2-4*box_r/r
    errors = dict(geometric_mass=float(np.max(np.abs(mass-expected_mass))),
                  four_dimensional_scalar_curvature=float(np.max(np.abs(R4-np.where(row['is_cap'],12.,0.)))))
    wall_r = chart.radius_at_clock(labels)
    cap = chart.from_wall_radii(wall_r,wall_r,region='cap')
    ext = chart.from_wall_radii(wall_r,wall_r,region='exterior')
    errors.update(wall_lapse=float(max(np.max(np.abs(cap['lapse']-1)),np.max(np.abs(ext['lapse']-1)))),
        angular_Israel=float(np.max(np.abs(ext['radius_x']-cap['radius_x']+2*wall_r))),
        temporal_Israel=float(np.max(np.abs(ext['lapse_x']-cap['lapse_x']+2))),
        Dirac_potential_normal_jump=float(np.max(np.abs(ext['potential_x']-cap['potential_x']))))
    maps = chart.maps(RMAX)
    quad32,quad64 = direct_wall_maps(RMAX,32),direct_wall_maps(RMAX,64)
    errors.update(wall_maps_vs_independent_quadrature=float(np.max(np.abs(maps-quad64))),
                  independent_quadrature_32_64=float(np.max(np.abs(quad32-quad64))))
    if np.min(row['lapse']) < chart.bounds['lapse_lower']:
        raise RuntimeError('sample violated the conservative analytic lapse bound')
    if max(errors.values()) > 1e-9:
        raise RuntimeError('quiet wall-radar metric/matching check exceeded 1e-9')
    probe = {}
    for region in ('nested','cap'):
        event = chart.from_wall_radii(.475,.525,region=region)
        probe[region] = {key:float(event[key]) for key in
            ('clock_u','clock_v','time','radar_position','radius','lapse','potential','potential_x')}
        probe[region]['received_over_wall_frequency'] = 1/float(event['lapse'])
    probe['interpretation'] = ('Photon emitted at wall tau(.475), received by a fixed radar-position '
        'observer at (U,V)=(tau(.475),tau(.525)); cap control continues the same cap '
        'and uses the same observer radar clock protocol. Not a Dirac quantum correction.')
    return dict(inputs=dict(H_ell=1.,k_wall_ell=2.,M_over_ell=MASS,
            wall_radius=[RMIN,RMAX],angular_k=1),
        domain=dict(null_label_range=[0.,chart.duration],radar_time=[0.,chart.duration],
            radar_position_max=chart.duration/2,
            diamond='0<=T-x<=duration and 0<=T+x<=duration',
            cap_triangle='V<=U',exterior_triangle='V>=U',
            is_a_complete_global_Cauchy_surface=False),
        analytic_bounds=chart.bounds,
        sampled=dict(metric_points=81,min_lapse=float(np.min(row['lapse'])),
            min_radius=float(np.min(r)),max_radius=float(np.max(r))),
        checks=errors,radar_probe=probe,
        counts=dict(wall_ODE_rhs=chart.wall_rhs_calls,new_P001_grids=0,
                    Dirac_time_steps=0,P001_time_steps=0,loop_polarizations=0),
        nested_Psi_over_Phi=None,
        scope='Exact quiet finite causal chart and background radar protocol; no state selection, Dirac propagation, renormalized metric response or nested slip.')
