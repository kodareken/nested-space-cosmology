"""Conserved central-mass tangent of the actual finite cap--wall geometry.

The source is present before the proper-wall clock window. Its interior is
outside this annulus and is not modeled. The initial physical wall radius
and speed, H and tension are fixed; the parent Schwarzschild mass changes
as required by both Israel constraints. This is a parametric source
susceptibility, not a retarded switched source or a choice of Dirac vacuum.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from .nsc_atlas_cap_wall_response import reference_coefficients
from .nsc_atlas_dirac_backreaction import MeanBackreaction
from .nsc_atlas_dirac_mean import RadarPointEvaluator

R0 = .45
K_WALL = 2.


def parent_mass_derivative():
    """d M_parent/db at b=0, holding R(0) and proper Rdot(0) fixed."""
    g = reference_coefficients(R0)
    return float(g['beta_plus']/g['beta_minus'])


def anchored_tortoise_mass_derivative(radius, *, region):
    """d r_star/d M at fixed H; the derivative vanishes at the corner R0."""
    r = np.asarray(radius,dtype=float)
    if not np.isfinite(r).all() or np.any(r<=0):
        raise ValueError('positive finite radius required')
    if region=='cap':
        if np.any(r>=1):
            raise ValueError('cap query exits its static optical branch')
        def primitive(x):
            return 2*np.log(x)-np.log1p(-x*x)+1/(1-x*x)
    elif region=='exterior':
        if np.any(r>=.8):
            raise ValueError('exterior query exits the white-hole optical branch')
        def primitive(x):
            return 2*np.log(.8-x)+1.6/(.8-x)
    else:
        raise ValueError('region must be cap or exterior')
    return primitive(r)-primitive(R0)


class _VacuumAnnulus:
    """Constant corner masses, no additional bulk stress in this annulus."""
    def __init__(self,chart):
        self.chart,self.region,self.start = chart,'nested',0.
        self.breakpoints = np.array([])
        self.points = RadarPointEvaluator(chart)

    def evaluate(self,U,V,*,region=None):
        fields = self.points(U,V,region=region)
        zero = np.zeros_like(fields['radius'])
        return dict(fields=fields,
            orthonormal=dict(energy_density=zero,radial_pressure=zero,
                angular_pressure=zero,outward_energy_flux=zero),
            null=dict(uu=zero,uv=zero,vv=zero))


class CompactMassTangent:
    """h_B and delta W at fixed physical proper-wall null labels.

    Reuse the existing GR wall response. Only the anchored optical tangent
    is integrated additionally. All returned changes are d/db at b=0 in
    ell=1, and are bounded to the original quiet causal diamond.
    """
    def __init__(self,chart):
        self.chart = chart
        self.parent_derivative = parent_mass_derivative()
        self.source = _VacuumAnnulus(chart)
        self.mean = MeanBackreaction(self.source,chart.duration,
            initial_masses=(self.parent_derivative,1.),
            wall_rtol=2e-11,wall_atol=2e-13)
        def rhs(t,unused):
            wall = self.mean.wall_state(t)
            g = reference_coefficients(wall['radius'])
            v,rate = float(g['v']),wall['displacement_rate']
            result = []
            for index,side in enumerate(('minus','plus')):
                beta = float(g['beta_'+side])
                db = wall['delta_beta'][index]
                result.extend((-(db+rate)/(beta+v)**2,
                               -(db-rate)/(beta-v)**2))
            return result
        solution = solve_ivp(rhs,(0.,chart.duration),np.zeros(4),method='DOP853',
            rtol=2e-11,atol=2e-13,max_step=chart.duration/40,dense_output=True)
        if not solution.success:
            raise RuntimeError(solution.message)
        self.optical_solution = solution.sol
        self.optical_rhs_calls = solution.nfev

    def _wall(self,clock,*,cap):
        shape = np.asarray(clock).shape
        labels,inverse = np.unique(np.asarray(clock,dtype=float).ravel(),return_inverse=True)
        states = [self.mean.wall_state(float(t)) for t in labels]
        r = np.array([s['radius'] for s in states])[inverse]
        rho = np.array([s['displacement'] for s in states])[inverse]
        rate = np.array([s['displacement_rate'] for s in states])[inverse]
        g = reference_coefficients(r)
        cap = np.broadcast_to(cap,shape).ravel()
        beta = np.where(cap,g['beta_minus'],g['beta_plus'])
        beta_r = np.where(cap,g['beta_minus_x'],g['beta_plus_x'])
        beta_rr = g['beta_second']
        v,acc = g['v'],g['radial_acceleration']
        difference = self.parent_derivative-1
        db = beta_r*rho+difference/(K_WALL*r*r)
        db_dot = beta_rr*v*rho+beta_r*rate-2*difference*v/(K_WALL*r**3)
        acc_r = g['beta_minus_x']**2+g['beta_minus']*beta_rr+1
        acceleration = (acc_r*rho+difference/K_WALL*(g['beta_plus_x']/r**2
            -2*g['beta_plus']/r**3)-self.parent_derivative/r**2)
        optical = self.optical_solution(labels)[:,inverse]
        index = np.where(cap,0,2)
        data = dict(rho=rho,rate=rate,beta=beta,v=v,
            u=optical[index,np.arange(inverse.size)],
            v_opt=optical[index+1,np.arange(inverse.size)],
            u_rate=-(db+rate)/(beta+v)**2,
            v_rate=-(db-rate)/(beta-v)**2,
            log_u=-(db+rate)/(beta+v),
            log_v=-(db-rate)/(beta-v),
            log_u_rate=-(db_dot+acceleration)/(beta+v)
                +(db+rate)*(beta_r*v+acc)/(beta+v)**2,
            log_v_rate=-(db_dot-acceleration)/(beta-v)
                +(db-rate)*(beta_r*v-acc)/(beta-v)**2)
        return {k:x.reshape(shape) for k,x in data.items()}

    def fields(self,U,V,*,region='nested',angular_k=1):
        if isinstance(angular_k,bool) or not isinstance(angular_k,int) or not angular_k:
            raise ValueError('a nonzero integer angular Dirac eigenvalue is required')
        U,V = np.broadcast_arrays(np.asarray(U,dtype=float),np.asarray(V,dtype=float))
        # The existing chart enforces the finite domain and one-sided choice.
        f = self.source.points(U,V,region=region)
        if angular_k!=1:
            f = self.chart.fields(U,V,region=region,angular_k=angular_k)
        cap = f['is_cap']
        wu,wv = self._wall(U,cap=cap),self._wall(V,cap=cap)
        r,C,F = f['radius'],f['conformal_factor'],f['static_f']
        ru,rv = (f['radius_t']-f['radius_x'])/2,(f['radius_t']+f['radius_x'])/2
        fp,fpp = -4*f['radius_uv']/C,-4*f['log_C_uv']/C
        dm = np.where(cap,1.,self.parent_derivative)
        K = np.empty_like(r)
        for mask,side in ((cap,'cap'),(~cap,'exterior')):
            if np.any(mask):
                K[mask] = anchored_tortoise_mass_derivative(r[mask],region=side)
        s = F*((wv['v_opt']-wu['u'])/2-dm*K)
        su = fp*ru*s/F-F*wu['u_rate']/2-2*dm*ru/(r*F)
        sv = fp*rv*s/F+F*wv['v_rate']/2-2*dm*rv/(r*F)
        df = fp*s-2*dm/r
        q = df/F+wu['log_u']+wv['log_v']
        qu = ((fpp*ru*s+fp*su+2*dm*ru/r**2)/F
              -df*fp*ru/F**2+wu['log_u_rate'])
        qv = ((fpp*rv*s+fp*sv+2*dm*rv/r**2)/F
              -df*fp*rv/F**2+wv['log_v_rate'])
        W = f['potential']
        Wu,Wv = (f['potential_t']-f['potential_x'])/2,(
            f['potential_t']+f['potential_x'])/2
        factor = q/2-s/r
        return dict(background=f,radius_change=s,radius_u_change=su,radius_v_change=sv,
            log_C_change=q,log_C_u_change=qu,log_C_v_change=qv,
            potential_change=W*factor,
            potential_u_change=Wu*factor+W*(qu/2-su/r+s*ru/r**2),
            potential_v_change=Wv*factor+W*(qv/2-sv/r+s*rv/r**2))
