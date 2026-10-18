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
from .nsc_atlas_dirac_radar import WallRadarChart,tortoise

R0 = .45
K_WALL = 2.


class CompactMassChart(WallRadarChart):
    """Exact real small-b vacuum family in the same proper-clock window.

    The body/inner-horizon interior is excluded. Positive b is the physical
    compact-source branch; zero recovers the quiet reference. The fixed
    annulus .3<r<.65 and 0<=b<=2e-4 are computational domains, not a global
    formation or finite-amplitude semiclassical validity theorem.
    """
    def __init__(self,reference,central_mass):
        b=float(central_mass)
        if not np.isfinite(b) or not 0<=b<=2e-4:
            raise ValueError('compact mass must lie in [0,2e-4] in ell=1')
        self.central_mass,self.duration = b,reference.duration
        g=reference_coefficients(R0)
        self.initial_speed=float(g['v'])
        bm0=np.sqrt(self.initial_speed**2+1-2*b/R0-R0*R0)
        self.parent_mass=R0*(1+self.initial_speed**2-(bm0-K_WALL*R0)**2)/2
        self._radius_to_clock={}
        self._nodes,self._weights=np.polynomial.legendre.leggauss(12)
        self._anchors=(float(tortoise(R0,'cap')),float(tortoise(R0,'exterior')))
        initial=[R0,-self._anchors[0],self._anchors[0],
                 -self._anchors[1],self._anchors[1]]
        def rhs(t,y):
            geom=self._geometry(y[0])
            v,bm,bp=geom['v'],geom['beta_minus'],geom['beta_plus']
            return [v,1/(bm+v),1/(bm-v),1/(bp+v),1/(bp-v)]
        solution=solve_ivp(rhs,(0.,self.duration),initial,method='DOP853',
            rtol=2e-12,atol=2e-14,max_step=self.duration/40,dense_output=True)
        if not solution.success:
            raise RuntimeError(solution.message)
        self.solution=solution.sol
        self.wall_rhs_calls=solution.nfev
        self.wall_end_radius=float(solution.y[0,-1])

    def _geometry(self,r):
        r=np.asarray(r)
        d=self.parent_mass-self.central_mass
        bm=d/(K_WALL*r*r)+3*r/4
        bp=bm-K_WALL*r
        bmr=-2*d/(K_WALL*r**3)+3/4
        bpr=bmr-K_WALL
        f=1-2*self.central_mass/r-r*r
        fr=2*self.central_mass/r**2-2*r
        v=np.sqrt(bm*bm-f)
        return dict(x=r,v=v,beta_minus=bm,beta_plus=bp,
            beta_minus_x=bmr,beta_plus_x=bpr,beta_second=6*d/(K_WALL*r**4),
            radial_acceleration=bm*bmr-fr/2,cap_f_prime=fr,k_wall=K_WALL)

    def wall_geometry(self,clock):
        return self._geometry(self.radius_at_clock(clock))

    def radius_at_clock(self,clock):
        t=np.asarray(clock,dtype=float)
        if not np.isfinite(t).all() or np.any((t<0)|(t>self.duration)):
            raise ValueError('null clock labels must lie in [0,wall duration]')
        r=self.solution(t.ravel())[0].reshape(t.shape)
        self._radius_to_clock.update((float(x),float(y)) for x,y in zip(r.ravel(),t.ravel()))
        return r

    def maps(self,radius):
        r=np.asarray(radius,dtype=float)
        if not np.isfinite(r).all() or np.any((r<R0)|(r>self.wall_end_radius)):
            raise ValueError('wall-radius query exits this family clock window')
        from scipy.optimize import brentq
        clocks=[]
        for x in r.ravel():
            key=float(x)
            if key not in self._radius_to_clock:
                self._radius_to_clock[key]=brentq(
                    lambda t:float(self.solution(t)[0])-key,0.,self.duration,
                    xtol=5e-15,rtol=8*np.finfo(float).eps)
            clocks.append(self._radius_to_clock[key])
        clocks=np.asarray(clocks)
        out=self.solution(clocks)
        out[0]=clocks
        return out.reshape((5,)+r.shape)

    def _optical_radius(self,z,*,cap):
        z=np.asarray(z,dtype=float)
        anchor=self._anchors[0 if cap else 1]
        mass=self.central_mass if cap else self.parent_mass
        def F(r):
            return 1-2*mass/r-r*r if cap else 1-2*mass/r
        r=np.clip(R0+F(R0)*(z-anchor),.301,.649)
        for unused in range(12):
            half=(r-R0)/2
            query=(r+R0)[...,None]/2+half[...,None]*self._nodes
            value=anchor+half*np.sum(self._weights/F(query),axis=-1)
            residual=value-z
            if np.max(np.abs(residual),initial=0.)<8e-15:
                return r
            r=r-residual*F(r)
            if np.any((r<=.3)|(r>=.65)):
                raise ValueError('optical inversion exits the admitted compact-source annulus')
        raise RuntimeError('compact-source optical inversion did not converge')

    def from_wall_radii(self,radius_u,radius_v,*,region='nested',angular_k=1):
        if isinstance(angular_k,bool) or not isinstance(angular_k,int) or not angular_k:
            raise ValueError('a nonzero integer angular Dirac eigenvalue is required')
        if region not in ('nested','cap','exterior'):
            raise ValueError('region must be nested, cap or exterior')
        ru,rv=np.broadcast_arrays(np.asarray(radius_u,dtype=float),np.asarray(radius_v,dtype=float))
        mu,mv=self.maps(ru),self.maps(rv)
        cap=ru>=rv if region=='nested' else np.full(ru.shape,region=='cap')
        gu,gv=self._geometry(ru),self._geometry(rv)
        bu,bv=[np.where(cap,g['beta_minus'],g['beta_plus']) for g in (gu,gv)]
        bru,brv=[np.where(cap,g['beta_minus_x'],g['beta_plus_x']) for g in (gu,gv)]
        vu,vv=gu['v'],gv['v']
        pu,qv=1/(bu+vu),1/(bv-vv)
        optical_u=np.where(cap,mu[1],mu[3])
        optical_v=np.where(cap,mv[2],mv[4])
        z=(optical_v-optical_u)/2
        r=np.empty_like(z)
        for mask,is_cap in ((cap,True),(~cap,False)):
            if np.any(mask):
                r[mask]=self._optical_radius(z[mask],cap=is_cap)
        mass=np.where(cap,self.central_mass,self.parent_mass)
        F=1-2*mass/r-np.where(cap,r*r,0.)
        fp=2*mass/r**2-np.where(cap,2*r,0.)
        fpp=-4*mass/r**3-np.where(cap,2.,0.)
        C=F*pu*qv
        if not np.isfinite(C).all() or np.any(C<=0):
            raise ValueError('compact-source chart lost its positive radar lapse')
        rU,rV=-F*pu/2,F*qv/2
        LU=-fp*pu/2-(vu*bru+gu['radial_acceleration'])/(bu+vu)
        LV=fp*qv/2-(vv*brv-gv['radial_acceleration'])/(bv-vv)
        N=np.sqrt(C)
        rt,rx=rU+rV,rV-rU
        Nt,Nx=N*(LU+LV)/2,N*(LV-LU)/2
        W=angular_k*N/r
        return dict(clock_u=mu[0],clock_v=mv[0],time=(mu[0]+mv[0])/2,
            radar_position=(mv[0]-mu[0])/2,is_cap=cap,radius=r,
            radius_t=rt,radius_x=rx,conformal_factor=C,lapse=N,radial_length=N,
            lapse_t=Nt,lapse_x=Nx,log_C_u=LU,log_C_v=LV,
            radius_uv=-C*fp/4,log_C_uv=-C*fpp/4,static_f=F,
            radial_curvature=-fpp,potential=W,
            potential_t=angular_k*(Nt/r-N*rt/r**2),
            potential_x=angular_k*(Nx/r-N*rx/r**2))


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
