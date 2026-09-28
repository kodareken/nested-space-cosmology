"""State matching for the already adopted finite proper-time light factor.

Apply the known fermionic thermal heat-image formula. This is a flat KMS
compatibility control, not a thermal description of the expanding child.
"""
from math import exp, log, log1p, pi

import numpy as np
from scipy.integrate import quad
from scipy.special import erfc


def canonical_thermal(temperature):
    """One ordinary four-component massless Dirac field, vacuum subtracted."""
    if not np.isfinite(temperature) or temperature<=0:
        raise ValueError('positive finite temperature required')
    pressure=7*pi*pi*temperature**4/180
    return {'free_energy':-pressure,'pressure':pressure,'rho':3*pressure,
            'entropy':4*pressure/temperature,'trace':0.}


def endpoint_images(temperature,cutoff=1.,terms=64):
    """Canonical thermal source plus the finite omitted-winding endpoint.

    Only the exponentially convergent endpoint difference is summed; the
    ordinary Stefan-Boltzmann result is an imported analytic input.
    """
    if not np.isfinite(cutoff) or cutoff<=0 or terms<1:
        raise ValueError('positive cutoff and image count required')
    standard=canonical_thermal(temperature)
    n=np.arange(1,terms+1,dtype=float);x=(n*cutoff/(2*temperature))**2
    signs=np.where(n.astype(int)%2,-1.,1.)
    factor=4*temperature**4/pi**2
    correction_F=float(-factor*np.sum(signs*(1+x)*np.exp(-x)/n**4))
    correction_rho=float(-3*correction_F+2*factor*np.sum(signs*x*x*np.exp(-x)/n**4))
    f=standard['free_energy']+correction_F;rho=standard['rho']+correction_rho
    pressure=-f
    correction={'free_energy':correction_F,'pressure':-correction_F,'rho':correction_rho,
                'entropy':(correction_rho-correction_F)/temperature,
                'trace':correction_rho+3*correction_F}
    # Magnitudes decrease for both alternating endpoint series. The next
    # term bounds the F remainder; combine the two series for rho.
    j=terms+1;y=(j*cutoff/(2*temperature))**2
    log_f=log(factor)+log1p(y)-y-4*log(j)
    log_rho=float(np.logaddexp(log(3)+log_f,log(2*factor)+2*log(y)-y-4*log(j)))
    floor=float(np.nextafter(0.,1.))
    bound_f=exp(max(log_f,log(floor)))
    bound_rho=exp(max(log_rho,log(floor)))
    return {'temperature':temperature,'cutoff':cutoff,'vacuum_free_energy':cutoff**4/(16*pi*pi),
            'canonical':standard,'finite_endpoint':{'free_energy':f,'pressure':pressure,'rho':rho,
                'entropy':(rho+pressure)/temperature,'trace':rho-3*pressure},
            'endpoint_minus_canonical':correction,
            'alternating_remainder_bound':{'free_energy':bound_f,'rho':bound_rho,
                'log10_free_energy':log_f/log(10),'log10_rho':log_rho/log(10),
                'scope':'analytic image truncation only; positive underflow floor; excludes floating-point roundoff'},'terms':terms}


def matsubara_source(temperature,cutoff=1.,terms=32):
    """Independent direct Matsubara heat integral and lapse/temperature variation.

    Frequencies are antiperiodic. The full E1 determinant includes the
    same vacuum term, which is subtracted after each metric variation.
    This route is efficient when T/cutoff is not very small.
    """
    if temperature<=0 or cutoff<=0 or terms<1:
        raise ValueError('positive scales and frequency count required')
    free=0.;energy=0.
    for j in range(terms):
        omega=(2*j+1)*pi*temperature
        a=(omega/cutoff)**2
        if a>740:continue
        # t=s*cutoff², so the lower endpoint is always 1.
        integrals=[quad(lambda t:exp(-a*t)*t**power,1,np.inf,
                        epsabs=1e-13,epsrel=3e-12)[0] for power in (-2.5,-1.5)]
        free+=temperature*cutoff**3*integrals[0]/(2*pi**1.5)
        energy+=temperature*omega**2*cutoff*integrals[1]/pi**1.5
    vacuum=cutoff**4/(16*pi*pi)
    return {'free_energy':free-vacuum,'pressure':vacuum-free,'rho':energy-vacuum,
            'entropy':(energy-free)/temperature,'trace':energy+3*free-4*vacuum,
            'full_free_energy':free,'full_rho':energy,'vacuum_free_energy':vacuum}


def canonical_integral(temperature):
    """Physical occupation integral, independent of the heat-image construction."""
    value=quad(lambda x:x**3*exp(-x)/(1+exp(-x)),0,np.inf,epsabs=2e-12)[0]
    rho=2*temperature**4*value/pi**2
    return {'rho':rho,'pressure':rho/3}


def mode_energy_derivative(energy,temperature,cutoff=1.,images=64):
    """Derivative wrt E of one +/-E Dirac pair's Euclidean free energy.

    The vacuum derivative is -erfc(E/cutoff), rather than -1. The
    canonical thermal increment is 2/(exp(E/T)+1). Their finite-endpoint
    discrepancy contains the off-shell heat contour contribution.
    """
    if energy<=0 or temperature<=0 or cutoff<=0:
        raise ValueError('positive energy, temperature and cutoff required')
    n=np.arange(1,images+1,dtype=float)
    signs=np.where(n.astype(int)%2,-1.,1.)
    a=energy/cutoff;b=n*cutoff/(2*temperature)
    # Gusev-Zelnikov's image heat kernel, integrated before differentiating.
    # exp(n E/T)*erfc(a+b) is evaluated without overflow.
    from scipy.special import erfcx
    thermal=float(np.sum(signs*(np.exp(-n*energy/temperature)*erfc(a-b)
                               +np.exp(-a*a-b*b)*erfcx(a+b))))
    finite_thermal=-thermal
    canonical=2*exp(-energy/temperature)/(1+exp(-energy/temperature))
    return {'vacuum_derivative':-erfc(a),'finite_thermal_derivative':finite_thermal,
            'canonical_thermal_derivative':canonical,
            'difference':finite_thermal-canonical}
