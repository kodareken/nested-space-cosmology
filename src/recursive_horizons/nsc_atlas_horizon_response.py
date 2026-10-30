"""Scalar-sourced cosmological horizon area, dynamical entropy and flux directions.

The retained canonical two-mode parent supplies all sources. This is the
quadratic null Einstein response with a bounded asymptotically de Sitter
future, not nonlinear collapse or a daughter-universe construction.
"""
from dataclasses import dataclass
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad_vec, solve_ivp

from .nsc_atlas_parent_circulation import (
    Parent, evolve, horizon_flux, static_patch_energy)

ROOT=Path(__file__).resolve().parents[2]
EIGENVALUES=np.array([3.,8.])
NORMALIZATIONS=np.sqrt([2.,12.])/math.pi


@dataclass
class HorizonSource:
    parent: Parent
    modes: object
    endpoint_u: float

    @property
    def nu_squared(self):
        return self.parent.mass_squared/self.parent.H**2-2.25

    def state(self,u):
        """Recover the original canonical field on the closed parent."""
        y=self.modes.sol(u)
        cosh=math.cosh(u);tanh=math.tanh(u);scale=cosh**-1.5
        z=self.parent.initial()
        z[:2]=[cosh/self.parent.H,math.sinh(u)]
        z[2:6]=[scale*y[0],self.parent.H*scale*(y[1]-1.5*tanh*y[0]),
                scale*y[2],self.parent.H*scale*(y[3]-1.5*tanh*y[2])]
        z[6]=z[3]*z[4]-z[2]*z[5]
        z[8]=u/self.parent.H
        return z

    def null_amplitudes(self,u):
        y=self.modes.sol(u);tanh=math.tanh(u)
        factor=self.parent.H*math.cosh(u)**-2.5
        return factor*NORMALIZATIONS*np.array([
            y[1]-2.5*tanh*y[0],
            tanh*y[3]+(1-3.5*tanh*tanh)*y[2]])

    def integrands(self,u):
        a,b=self.null_amplitudes(u)
        tanh=math.tanh(u)
        F=4*math.pi/(3*self.parent.H**2)*(a*a+b*b)
        return np.array([tanh*F,F/math.cosh(u),tanh*a*a,tanh*a*b,tanh*b*b])

    def tail_bounds(self,U=None):
        """Analytic tails using the nonincreasing rescaled oscillator energies.

        Bounds apply to exact mode solutions at the declared initial constants;
        numerical integration error is checked separately, not certified here.
        """
        U=self.endpoint_u if U is None else U
        q=self.parent.initial()[[2,4]]
        twice_E0=(self.nu_squared+EIGENVALUES+.75)*q*q
        D2=twice_E0*(1+2.5/math.sqrt(self.nu_squared))**2
        BF=float(4*math.pi/3*np.sum(NORMALIZATIONS**2*D2))
        P=32*BF*math.exp(-5*U)/5
        S=64*BF*math.exp(-6*U)/6
        return dict(endpoint_Ht=U,F_sech5_coefficient=BF,
                    P_tail_upper_bound=P,S_tail_upper_bound=S,
                    energy_tail_upper_bound=P/self.parent.H,
                    normalized_area_deficit_tail_upper_bound=self.parent.kappa*P/(4*math.pi),
                    angular_matrix_spectral_norm_tail_upper_bound=3*self.parent.H**2*P/(4*math.pi))


def solve_source(*,parent=None,endpoint_u=8.,rtol=2e-12,atol=2e-14):
    p=parent or Parent.retained()
    if not 2<=endpoint_u<=10:
        raise ValueError('finite future endpoint 2 <= Ht <= 10 required')
    nu2=p.mass_squared/p.H**2-2.25
    if nu2<=0:
        raise ValueError('positive rescaled oscillator floor required')
    q=p.initial()[[2,4]]
    def rhs(u,y):
        w2=nu2+(EIGENVALUES+.75)/math.cosh(u)**2
        return [y[1],-w2[0]*y[0],y[3],-w2[1]*y[2]]
    sol=solve_ivp(rhs,(0.,endpoint_u),[q[0],0.,q[1],0.],method='DOP853',
                  rtol=rtol,atol=atol,dense_output=True,max_step=.015)
    if not sol.success or sol.nfev>140000:
        raise RuntimeError('finite rescaled scalar solve failed or budget exceeded')
    return HorizonSource(p,sol,endpoint_u)


def solve_response(source,*,rtol=2e-11,atol=1e-13):
    # P and S are future integrals. Q is the angular Killing-energy matrix.
    sol=solve_ivp(lambda u,_:-source.integrands(u),(source.endpoint_u,0.),np.zeros(5),
                  method='DOP853',rtol=rtol,atol=atol,dense_output=True,max_step=.015)
    if not sol.success or sol.nfev>140000:
        raise RuntimeError('finite null response failed or budget exceeded')
    return sol


def initial_area_coefficient(parent=None):
    p=parent or Parent.retained();m2=p.mass_squared/p.H**2
    return 2*(m2+4)/(5*(m2+3))+12*(2*m2+13)/(35*(m2+8))


def area_record(source,response,u):
    """All area changes are coefficients of epsilon^2 relative to 4*pi/H^2."""
    p=source.parent;P,S=response.sol(u)[:2]
    deficit=p.kappa/(4*math.pi)*(P-math.sinh(u)*S)
    growth_rate=p.kappa/(4*math.pi)*math.cosh(u)*S
    dynamical_deficit=deficit+math.tanh(u)*growth_rate
    coefficient=initial_area_coefficient(p)
    EK=static_patch_energy(u/p.H,source.state(u),p)
    expected=p.kappa*p.H*EK/(4*math.pi)
    return dict(Ht=u,proper_readout=u/p.H,
                ordinary_area_deficit_over_epsilon_squared=deficit,
                ordinary_area_growth_over_epsilon_squared=coefficient-deficit,
                normalized_area_growth_rate_over_epsilon_squared=growth_rate,
                dynamical_area_deficit_over_epsilon_squared=dynamical_deficit,
                dynamical_area_growth_over_epsilon_squared=coefficient-dynamical_deficit,
                finite_endpoint_correction_over_epsilon_squared=math.tanh(u)*growth_rate,
                patch_Killing_energy_over_epsilon_squared=EK,
                comparison_first_law_absolute_residual=abs(dynamical_deficit-expected),
                other_future_background_expansion_over_H=4*math.tanh(u))


def angular_record(source,response):
    p=source.parent;tail=source.tail_bounds()
    q11,q12,q22=response.sol(0.)[2:]
    Q=np.array([[q11,q12,0.],[q12,q22,0.],[0.,0.,0.]])
    eigenvalues,vectors=np.linalg.eigh(Q)
    axis=vectors[:,-1]
    if axis[0]<0:axis=-axis
    trace=float(np.trace(Q));ratio=float(eigenvalues[-1]/trace)
    c=math.cos(math.pi/6)
    cap_fraction=ratio*(1-c**3)+(1-ratio)/2*(2-3*c+c**3)
    return dict(matrix=Q.tolist(),eigenvalues=eigenvalues.tolist(),
                axis_in_n_coordinates=axis.tolist(),
                axis_angle_from_n1_degrees=math.degrees(math.atan2(axis[1],axis[0])),
                largest_eigenvalue_fraction=ratio,maximum_over_sphere_average=3*ratio,
                two_30_degree_caps_energy_fraction=cap_fraction,
                two_caps_area_fraction=1-c,
                total_energy_from_trace=4*math.pi*trace/(3*p.H**3),
                spectral_norm_tail_bound=tail['angular_matrix_spectral_norm_tail_upper_bound'],
                antipodal_energy_equal=True,
                interpretation='Directions of energy transfer on the same smooth horizon. These are not formed children or a physical cap-size selection.')


def calculation():
    source=solve_source();response=solve_response(source);p=source.parent
    prior=json.loads((ROOT/'results/summaries/nsc-atlas-parent-circulation-v1.json').read_bytes())
    endpoint=prior['background']['finite_proper_endpoint'];uend=p.H*endpoint
    points=[0.,uend,3.]
    rows=[area_record(source,response,u) for u in points]
    optical=[]
    for light in prior['light_comparisons']:
        u=p.H*light['return_proper'];row=area_record(source,response,u)
        optical.append(dict(radius=light['radius'],Ht=u,
            fractional_light_delay_over_dynamical_area_growth=
                light['fractional_delay_over_initial_energy_fraction']/row['dynamical_area_growth_over_epsilon_squared'],
            fractional_light_delay_over_ordinary_area_growth=
                light['fractional_delay_over_initial_energy_fraction']/row['ordinary_area_growth_over_epsilon_squared'],
            dynamical_area_growth_over_epsilon_squared=row['dynamical_area_growth_over_epsilon_squared']))
    reference=evolve(parent=p,clock='proper',endpoint=endpoint)
    times=np.linspace(0.,endpoint,401)
    difference=max(np.max(abs(source.state(p.H*t)[:6]-reference.sol(t)[:6])) for t in times)
    quadrature,estimate=quad_vec(source.integrands,0.,source.endpoint_u,epsabs=1e-12,epsrel=2e-12,limit=1200)
    tighter_source=solve_source(rtol=3e-13,atol=3e-15)
    tighter_response=solve_response(tighter_source,rtol=3e-12,atol=1e-14)
    spread=max(abs(area_record(tighter_source,tighter_response,u)['ordinary_area_deficit_over_epsilon_squared']-
                   row['ordinary_area_deficit_over_epsilon_squared']) for u,row in zip(points,rows))
    controls=dict(mode_rhs_evaluations=source.modes.nfev,null_rhs_evaluations=response.nfev,
                  original_proper_solution_max_state_difference=float(difference),
                  independent_future_integral_max_difference=float(np.max(abs(quadrature-response.sol(0.)))),
                  quadrature_error_estimate=float(estimate),
                  tighter_solve_area_coefficient_max_absolute_spread=spread,
                  initial_area_coefficient_absolute_difference=abs(rows[0]['ordinary_area_deficit_over_epsilon_squared']-initial_area_coefficient(p)),
                  comparison_first_law_max_absolute_residual=max(r['comparison_first_law_absolute_residual'] for r in rows))
    return dict(schema='nsc-atlas-horizon-response-v1',
        inherited=dict(vacuum=p.vacuum,H=p.H,mass_squared=p.mass_squared,
                       m_squared=p.mass_squared/p.H**2,kappa=p.kappa,
                       initial_energy_fraction='epsilon^2',unique_preparation=False),
        exact=dict(initial_normalized_area_deficit=initial_area_coefficient(p),
                   coefficient_formula='2(m^2+4)/(5(m^2+3))+12(2m^2+13)/(35(m^2+8))',
                   source_dependent_future_eigenvalues=[3,8],
                   relation='H*(A_vac-A_dyn)=kappa*epsilon^2*E_K; A_dyn=A-(tanh(Ht)/H)*A_dot',
                   entropy_convention='S_dyn=2*pi*A_dyn/kappa; T=H/(2*pi), hbar=c=kB=1',
                   perturbative_order='epsilon^2 with higher-order backreaction not bounded',
                   future_boundary='Bounded null expansion and asymptotically retained de Sitter area'),
        area_history=rows,light_area_relations=optical,angular_transfer=angular_record(source,response),
        infinite_future_tail=source.tail_bounds(),controls=controls,
        scope=['Explicit scalar-sourced area/entropy response using the published dynamical horizon first law',
               'A monotone geometric record on the chosen expanding branch; not an ultimate origin of time',
               'A cosmological horizon: its other background future expansion is positive for t>0',
               'No SI Lambda, matter-antimatter imbalance, new trapped surface or daughter formation derived'],
        no_observational_fit=True)
