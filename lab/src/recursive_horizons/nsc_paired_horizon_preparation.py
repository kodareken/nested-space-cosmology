"""Signed-angular horizon preparation with explicit finite-collar matching.

The existing horizon/incoming state and scattering owners supply all inputs.
The retained seed covariance is never an input. Full spatial PG covariance
still requires a globally normalized mode resolution on both sides of the
horizon; a finite horizon-to-seed map is not that resolution.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import atan2, exp, log, pi, sqrt

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.special import expit

from .nsc_unruh_state import ParentDirac, horizon_frame, horizon_covariance
from .nsc_compact_ctp_neck import _massive_reflection
from .nsc_transmitting_dirac_domain import S1, S2, S3, I2, MODE_TO_CURRENT, TransmittingDiracSeamDomain
from .nsc_lorentzian import geometry
from .nsc_common_time_bulk_split import CommonTimeBulkSplit


def source_covariance(E,kappa,omega,mass):
    """Owned affine horizon pair plus the inherited open incoming occupation."""
    if E==0 or min(kappa,omega)<=0 or mass<0:raise ValueError('nonzero signed energy and fixed source scales required')
    f=float(expit(-2*pi*E/kappa));n=float(expit(-2*pi*E/(omega*kappa)))
    if mass>0 and abs(E)<omega*mass:n=0. if E>0 else 1.
    x=pi*E/kappa;s=exp(-abs(x))/(1+exp(-2*abs(x)))
    return np.array([[f,-1j*s,0.],[1j*s,1-f,0.],[0.,0.,n]],complex)


@dataclass(frozen=True)
class PairedHorizonSeedMap:
    horizon_rho: float
    surface_gravity: float
    inheritance_ratio: float
    collar_delta_q: float
    reflection_tolerance: float
    outer_floor: float=60.

    @property
    def background(self):return ParentDirac(self.horizon_rho,self.surface_gravity)

    @property
    def horizon_radius(self):return sqrt(1+self.horizon_rho**2)

    def working_frame(self,y,E,mass,angular,*,legacy_compact_distance=False):
        delta=exp(y);rh=self.horizon_radius
        signed_offset=-rh*rh*delta+self.horizon_rho*rh*rh*delta**2
        effective=np.hypot(mass*rh,angular)
        phase=atan2(mass*rh,angular)
        rotation=np.diag(np.exp(np.array([-1j,1j])*(pi/4+phase/2)))
        distance=sqrt(2*delta/self.surface_gravity)
        if legacy_compact_distance:distance/=rh
        return rotation@horizon_frame(E,self.surface_gravity,effective,distance,
                                      self.background.near_tortoise(signed_offset),timelike=True)

    def generator(self,y,E,mass,angular):
        delta=exp(y);b=self.background;w=b.interior_W(delta)
        radius=1/np.sin(b.horizon_q-delta)
        return delta/sqrt(w)*(-mass*radius*S1+angular*S2)-E*delta/w*S3

    def frame_residual(self,E,mass,angular,*,legacy_compact_distance=False):
        y=log(self.collar_delta_q);h=2e-4
        f=lambda value:self.working_frame(value,E,mass,angular,legacy_compact_distance=legacy_compact_distance)
        derivative=(-f(y+2*h)+8*f(y+h)-8*f(y-h)+f(y-2*h))/(12*h)
        residual=derivative+1j*self.generator(y,E,mass,angular)@f(y)
        return {'differential_norm':float(np.linalg.norm(residual)),
                'off_diagonal_maximum':float(max(abs(residual[0,1]),abs(residual[1,0]))),
                'finite_difference_log_step':h,
                'distance_factor':1/self.horizon_radius if legacy_compact_distance else 1.}

    @lru_cache(maxsize=128)
    def scattering(self,E,mass,angular):
        if mass==0:
            # The signed angular frame includes the constant phase rotation.
            # In that frame the m=0 eta partners have the same R coefficient.
            if angular<0:return self.scattering(E,mass,-angular)
            value=self.background.reflection(E,angular,tolerance=self.reflection_tolerance)
        else:
            value=_massive_reflection(self.background,E,angular,mass,
                inner_offset=self.collar_delta_q,tolerance=self.reflection_tolerance,outer_floor=self.outer_floor)
        return value

    def at_seed(self,E,mass,angular,*,legacy_compact_distance=False):
        return self.at_radius(0.,E,mass,angular,legacy_compact_distance=legacy_compact_distance)

    def pg_generator_residual(self,rho,E,mass,angular):
        beta,bp,_=geometry(rho);radius=sqrt(1+rho*rho);a=sqrt(beta*beta-1)
        seam=TransmittingDiracSeamDomain(1.,1.,float(beta),radius)
        trace,_,_=seam.trace_map(np.ones(1));S=radius*trace[0]
        H=-mass*S1+angular/radius*S2-E/a*S3
        converted=(-beta*bp/(2*a*a)*I2-bp/(2*a*a)*S2-1j*E*beta/(a*a)*I2
                   +S@(1j/a*H)@np.linalg.inv(S))
        v=CommonTimeBulkSplit.principal(1.,1.,float(beta))
        potential=CommonTimeBulkSplit.local_potential(1.,radius,mass,angular)
        direct=np.linalg.solve(v,1j*(E*I2-potential)+.5*bp*I2)
        return float(np.linalg.norm(converted-direct))

    def at_radius(self,rho,E,mass,angular,*,legacy_compact_distance=False):
        """Physical input -> finite modal map, never stored X -> guessed Y.

        The legacy flag is a controlled reproduction of the frozen collar
        convention, not an alternative physical normalization to fit.
        """
        if not -1.<=rho<=1.:
            raise ValueError('this evaluated mode chart is the trapped probe interval [-1,1]')
        if E<0:
            # The angular charge-bundle map supplies eta reversal and its
            # known unit-modulus phase. The radial part is the locked sigma3.
            # Compute from mode columns and the boundary state, not from a
            # complemented stored seed matrix.
            positive=self.at_radius(rho,-E,mass,-angular,legacy_compact_distance=legacy_compact_distance)
            source=source_covariance(E,self.surface_gravity,self.inheritance_ratio,mass)
            modes=S3@positive['seed_mode_map'].conj()
            trace_modes=S3@positive['PG_trace_mode_map'].conj()
            seed=modes@source@modes.conj().T
            trace=trace_modes@source@trace_modes.conj().T
            beta=float(geometry(rho)[0]);radius=sqrt(1+rho*rho)
            _,inverse,_=TransmittingDiracSeamDomain(1.,1.,beta,radius).trace_map(np.ones(1))
            eig=np.linalg.eigvalsh(seed)
            return {
                'source_covariance':source,'seed_mode_map':modes,'seed_covariance':seed,
                'PG_trace_mode_map':trace_modes,'PG_trace_covariance':trace,
                'PG_spatial_half_density_mode_map':radius*trace_modes,
                'eigenvalues':eig,'signed_energy':E,'rho':rho,
                'legacy_compact_distance':legacy_compact_distance,'spatial_PG_covariance':None,
                'positive_partner_scattering':positive['scattering'],
                'residuals':{
                    'modal_coisometry':float(np.linalg.norm(modes@modes.conj().T-np.eye(2))),
                    'T_trace_pullback':float(np.linalg.norm(inverse[0]@trace@inverse[0].conj().T-seed)),
                    'KS_to_PG_Dirac_generator':self.pg_generator_residual(rho,E,mass,angular),
                    'seed_CAR_lower_violation':max(0.,-float(eig.min())),
                    'seed_CAR_upper_violation':max(0.,float(eig.max())-1.),
                    'boundary_state_signed_relation':float(np.linalg.norm(source-(np.eye(3)-positive['source_covariance'].conj()))),
                    'signed_modal_covariance':float(np.linalg.norm(seed-(np.eye(2)-S3@positive['seed_covariance'].conj()@S3))),
                },
            }
        scatter=self.scattering(E,mass,angular)
        R=complex(scatter['reflection']);T=float(scatter['transmission'])
        source=source_covariance(E,self.surface_gravity,self.inheritance_ratio,mass)
        sewing=np.array([[0.,1.,0.],[R,0.,sqrt(max(0.,T))]],complex)
        compressed=sewing@source@sewing.conj().T
        direct=horizon_covariance(E,self.surface_gravity,R,float(source[2,2].real))
        target_q=pi/2+np.arctan(rho)
        begin=log(self.collar_delta_q);end=log(self.background.horizon_q-target_q)
        initial=self.working_frame(begin,E,mass,angular,legacy_compact_distance=legacy_compact_distance)
        run=solve_ivp(lambda y,u:(-1j*self.generator(y,E,mass,angular)@u.reshape(2,2)).ravel(),
                      (begin,end),initial.ravel(),method='DOP853',rtol=2e-12,atol=2e-14,max_step=.1)
        if not run.success:raise ArithmeticError(run.message)
        transported=run.y[:,-1].reshape(2,2)
        mode_map=transported@sewing
        seed=mode_map@source@mode_map.conj().T
        beta,beta_prime,_=geometry(rho);radius=sqrt(1+rho*rho)
        seam=TransmittingDiracSeamDomain(1.,1.,float(beta),radius)
        trace,inverse,gram=seam.trace_map(np.ones(1))
        # static t = tau - F(rho), F'=beta/A; the common origin is at rho=0.
        phase_integral=quad(lambda x:float(geometry(x)[0]/geometry(x)[2]),0.,rho,
                            epsabs=2e-13,epsrel=2e-13)[0] if rho else 0.
        pg_trace_map=np.exp(1j*E*phase_integral)*trace[0]@mode_map
        pg_trace=pg_trace_map@source@pg_trace_map.conj().T
        pulled=inverse[0]@pg_trace@inverse[0].conj().T
        eig=np.linalg.eigvalsh(seed)
        return {
            'source_covariance':source,'sewing':sewing,'seed_mode_map':mode_map,
            'seed_covariance':seed,'PG_trace_mode_map':pg_trace_map,
            'PG_trace_covariance':pg_trace,
            'residuals':{
                'scattering_current':float(scatter['current_defect']),
                'source_compression_vs_owned_horizon':float(np.linalg.norm(compressed-direct)),
                'modal_coisometry':float(np.linalg.norm(mode_map@mode_map.conj().T-np.eye(2))),
                'T_trace_pullback':float(np.linalg.norm(pulled-seed)),
                'KS_to_PG_Dirac_generator':self.pg_generator_residual(rho,E,mass,angular),
                'seed_CAR_lower_violation':max(0.,-float(eig.min())),
                'seed_CAR_upper_violation':max(0.,float(eig.max())-1.),
            },
            'scattering':scatter,'eigenvalues':eig,
            'legacy_compact_distance':legacy_compact_distance,
            'signed_energy':E,
            'rho':rho,'PG_spatial_half_density_mode_map':radius*pg_trace_map,
            'spatial_PG_covariance':None,
        }

    def require_global_PG_covariance(self):
        raise ValueError('finite modal seed/trace map lacks full global PG spectral/current completeness; it is not C_PG')
