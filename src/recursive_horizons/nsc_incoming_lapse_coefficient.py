"""Closed incoming r_TTT response of the already declared local/reference action."""
import mpmath as mp
import numpy as np
import sympy as sp

from .nsc_incoming_cauchy_jets import incoming_cauchy_jets
from .nsc_incoming_local_constraints import local_action_densities,scalar_taylor_coefficients,taylor_two_jets
from .nsc_light_restoration_action import LightRestorationAction
from .nsc_magnetic_light_reference import magnetic_light_spectrum
from .nsc_spherical_local_history import LockedSphericalLocalAction
from .nsc_incoming_vacuum_tail_bound import _precision,_lo,_hi,_range


def reference_identity():
    m,L,p,Ha,Hr,r,u=sp.symbols('m L p H_a H_r r u',real=True,nonzero=True)
    h=sp.Matrix([-m,L,p]);hd=sp.Matrix([0,-L*Hr,-p*Ha]);variation=sp.Matrix([0,-L*u/r,0])
    expected=L*L*u/r*(m*m*Hr+p*p*(Hr-Ha))
    residual=sp.simplify(h.cross(hd).dot(h.cross(variation))-expected)
    M=sp.symbols('M',positive=True)
    I0=sp.simplify(sp.beta(sp.Rational(1,2),3)/M**6-sp.Rational(16,15)/M**6)
    I2=sp.simplify(sp.expand_func(sp.beta(sp.Rational(3,2),2))/M**4-sp.Rational(4,15)/M**4)
    # Gamma form removes representation-dependent unevaluated beta functions.
    I0=sp.simplify(I0.rewrite(sp.gamma));I2=sp.simplify(I2.rewrite(sp.gamma))
    if [residual,I0,I2]!=[0,0,0]:raise ArithmeticError('closed reference moment identity failed')
    return {'cross_product_residual':str(residual),'integral_residuals':[str(I0),str(I2)],
            'I0':'16/(15*M^6)','I2':'4/(15*M^4)',
            'reference_group_coefficient':'-d*a*L^2/(120*pi*r)*(4*m^2*H_r/M^6+(H_r-H_a)/M^4)',
            'multiplicity':'d=copy_count*degeneracy; full signed group and full momentum line',
            'LLL_and_zero_angular_response':0}


def _hull(expression,stored):
    value=mp.mpf(float(stored))
    return _range(min(_lo(expression),value),max(_hi(expression),value))


def _interval(value):
    # Round the public binary endpoints OUTWARD as well.
    return [float(np.nextafter(float(_lo(value)),-np.inf)),float(np.nextafter(float(_hi(value)),np.inf))]


def closed_lapse_coefficient(channels,ledger,*,precision=50):
    domain=incoming_cauchy_jets();normal=domain.normal_geometry()
    if ledger['magnetic_flux']!=4 or [c['index'] for c in channels]!=list(range(33)):
        raise ValueError('unchanged q4 retained33 inventory required')
    with _precision(precision):
        pi=_hull(mp.iv.pi,np.pi)
        a=_hull(mp.iv.sqrt(3*mp.iv.pi/2-4),normal['intrinsic_a'])
        r=_hull(mp.iv.sqrt(2),normal['intrinsic_r'])
        Ha=_hull((6-3*mp.iv.pi/2)/(2*a),normal['k_parallel'])
        Hr=_hull(-a/2,normal['k_perp'])
        hq=_hull(mp.iv.euler/2-mp.iv.mpf(25)/48,float(magnetic_light_spectrum(4).fourth_order_harmonic))
        CW=mp.iv.mpf(float(ledger['C_Weyl']))
        pieces=[];reference=mp.iv.mpf(0)
        for c in channels[1:]:
            m=_hull(c['compact_level']*mp.iv.pi/2,c['compact_mass'])
            ell=_hull(mp.iv.sqrt(c['angular_level']*(c['angular_level']+4)),c['angular_eigenvalue'])
            L=ell/r;M2=m*m+L*L
            if _lo(M2)<=0:raise ValueError('gapped non-LLL channel required')
            multiplicity=c['copy_count']*c['degeneracy']
            value=-multiplicity*a*L*L/(120*pi*r)*(4*m*m*Hr/M2**3+(Hr-Ha)/M2**2)
            reference+=value;pieces.append({'group':c['index'],'coefficient_interval':_interval(value)})
        local={
            'light/bar_radial_R2_squared':-hq*a*r*(Ha-Hr)/(30*pi),
            'light/WZ_Weyl':-mp.iv.ln(r)*a*r*(Ha-Hr)/(30*pi),
            'light/WZ_boxR':-a*r*(Ha+Hr)/(60*pi),
            'compact/weyl_bulk':32*pi*CW*a*r*(Ha-Hr)/3,
        }
        local_sum=sum(local.values(),mp.iv.mpf(0));total=reference+local_sum
        return {'reference_coefficient_interval':_interval(reference),
                'per_group_reference':pieces,'local_coefficient_intervals':{k:_interval(v) for k,v in local.items()},
                'local_sum_interval':_interval(local_sum),'total_coefficient_interval':_interval(total),
                'strictly_positive':_lo(total)>0,'precision':precision,
                'zero_channels':['light/cylinder','light/LLL_geometry_bar','light/WZ_Euler','light/WZ_gauge',
                                 'compact/einstein_bulk','compact/maxwell_bulk','compact/euler_bulk_diagnostic'],
                'input_interpretation':{'geometry_pi_harmonic':'hulls include exact owner expressions and stored floating evaluations',
                                        'C_Weyl':'the unchanged declared binary ledger coefficient, not an error bound for an unprovided exact coefficient',
                                        'mass_angular':'hulls include existing definitions and stored labels'}}


def local_density_stencil(ledger,step):
    """Exact mixed polynomial stencil; only arithmetic is approximate.

    Homogeneous curvature has no N_TT. The r_TTT Euler coefficient is
    -partial_N_T partial_r_TT L. Degree<=4 excludes odd mixed powers other
    than N_T*r_TT, so four corners extract the coefficient without truncation.
    """
    if step not in(.5,1.):raise ValueError('the two declared numerical stencil controls required')
    domain=incoming_cauchy_jets();raw=taylor_two_jets(scalar_taylor_coefficients(domain.fields),0.,0.)[0]
    light=LightRestorationAction(magnetic_light_spectrum(4));compact=LockedSphericalLocalAction(ledger)
    values={}
    for sn in(-1,1):
        for sr in(-1,1):
            perturbed=raw.copy();perturbed[0,1]+=sn*step;perturbed[3,3]+=sr*step
            for name,value in local_action_densities(perturbed,light,compact).items():
                values[name]=values.get(name,0)-sn*sr*value/(4*step*step)
    return {k:float(v.real) for k,v in values.items()}
