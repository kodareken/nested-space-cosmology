"""Exterior massive mode amplitudes, beyond a stored reflection probability.

The NSC metric's existing asymptotic series and exterior Dirac equation are
reused. This owner retains the complex Jost amplitude needed for spatial mode
fields and assesses the outer boundary error. No seed state is an input.
"""
from dataclasses import dataclass
from math import atan2, exp, log, sqrt

import numpy as np
from scipy.integrate import solve_ivp

from .nsc_unruh_state import horizon_frame


def outgoing_ratio(energy, mass, angular, radius, order=8):
    """NSC-specific massive extension of the owned exterior ratio series.

    Coefficients solve the same Riccati equation as _massive_reflection; no
    new Hamiltonian or interaction enters. The returned remainder is checked
    against the actual metric by the caller, not assumed to be a proof bound.
    """
    if energy <= 0 or mass < 0 or energy == mass or radius <= 0 or order < 1:
        raise ValueError('positive energy/radius away from threshold required')
    momentum = np.sqrt(complex(energy*energy-mass*mass))
    metric = np.zeros(order+2); metric[0] = 1.
    for j in range((order+2)//2):
        power = 2*j+1
        if power < len(metric): metric[power] = -6*(-1)**j/((2*j+1)*(2*j+3))
    root = np.zeros_like(metric); root[0] = 1.
    for n in range(1,len(root)):
        root[n] = (metric[n]-sum(root[j]*root[n-j] for j in range(1,n)))/2
    sphere = np.zeros_like(metric); sphere[0] = 1.; coefficient = 1.
    for n in range(1,len(sphere)//2):
        coefficient *= (-.5-(n-1))/n
        sphere[2*n] = coefficient
    first = np.zeros_like(metric)
    first[1:] = angular*np.convolve(root,sphere)[:len(metric)-1]
    second = mass*root
    minus = first-1j*second; plus = first+1j*second
    ratio = np.zeros(order+1,complex)
    ratio[0] = 1j*mass/(energy+momentum) if mass else 0.
    for n in range(1,order+1):
        derivative = -sum(metric[n-j-1]*j*ratio[j] for j in range(1,n))
        square = np.convolve(ratio,ratio)
        forcing = 1j*(plus[n]+np.convolve(minus,square)[n])
        ratio[n] = (forcing-derivative)/(2j*momentum)
    powers = radius**(-np.arange(order+1,dtype=float))
    value = ratio@powers
    derivative = -(np.arange(order+1)*ratio@powers)/radius
    return complex(value), complex(derivative), ratio


@dataclass
class ExteriorJostMode:
    background: object
    energy: float
    mass: float
    angular: float
    radial_collar: float
    outer_radius: float
    run: object
    horizon_coefficients: np.ndarray
    reflection: complex
    transmission: float
    residuals: dict

    @property
    def transmission_phase(self):
        """Jost phase in the current-normalized outer outgoing coordinate."""
        return complex(np.exp(-1j*self.run.y[1,-1].imag)*np.conj(self.horizon_coefficients[0])/abs(self.horizon_coefficients[0]))

    @property
    def complex_transmission(self):
        return sqrt(self.transmission)*self.transmission_phase

    def field(self, rho):
        """Static two-spinor with unit outgoing-horizon amplitude."""
        if not self.background.horizon_rho+self.radial_collar <= rho <= self.outer_radius:
            raise ValueError('Jost field evaluation outside its controlled radial interval')
        ratio, log_amplitude = self.run.sol(log(rho-self.background.horizon_rho))
        log_inner = self.run.y[1,-1]
        amplitude = np.exp(log_amplitude-log_inner)/self.horizon_coefficients[0]
        return amplitude*np.array([1.,ratio],complex)


def solve_jost(background, energy, mass, angular, *, radial_collar=1e-10,
               outer_radius=None, order=8, rtol=2e-13, atol=2e-15):
    momentum = np.sqrt(complex(energy*energy-mass*mass))
    if energy <= 0 or mass < 0 or energy == mass:
        raise ValueError('positive frequency away from massive threshold required')
    end = max(60.,30*max(1.,abs(angular),mass)/max(abs(momentum),.2)) if outer_radius is None else outer_radius
    for outer_expansions in range(16):
        ratio, derivative, _ = outgoing_ratio(energy,mass,angular,end,order)
        metric_a = background.A_from_offset(end-background.horizon_rho)
        v1 = angular*sqrt(metric_a)/sqrt(1+end*end);v2 = mass*sqrt(metric_a)
        asymptotic_residual = abs(metric_a*derivative-1j*(v1+1j*v2)+2j*energy*ratio-1j*(v1-1j*v2)*ratio*ratio)
        initial_current = float(1-abs(ratio)**2)
        if outer_radius is not None or (asymptotic_residual < 3e-13 and (energy>mass or abs(initial_current)<3e-13)):
            break
        end *= 2
    else: raise ArithmeticError('outer mode expansion failed its algebraic accuracy gate')
    if energy > mass and initial_current <= 0:
        raise ArithmeticError('outer outgoing current is unresolved')
    if energy < mass and abs(initial_current) > 1e-8:
        raise ArithmeticError('decaying outer expansion is not current-null to its declared accuracy')
    def rhs(y,state):
        offset = exp(y);rho = background.horizon_rho+offset
        a = background.A_from_offset(offset)
        first = angular*sqrt(a)/sqrt(1+rho*rho); second = mass*sqrt(a)
        z = state[0]
        return [offset/a*(1j*(first+1j*second)-2j*energy*z+1j*(first-1j*second)*z*z),
                1j*offset/a*(energy-(first-1j*second)*z)]
    if energy < mass:
        # The decaying solution has exactly zero current. Evolving its phase
        # preserves this fact without subtracting exponentially large modes.
        # The finite series' norm error is retained in closed_outer_current.
        def phase_rhs(y,state):
            theta = float(state[0].real); z = np.exp(1j*theta)
            dr,dlog = rhs(y,np.array([z,state[1]],complex))
            return [float((dr/(1j*z)).real),dlog]
        run = solve_ivp(phase_rhs,(log(end-background.horizon_rho),log(radial_collar)),
                        np.array([np.angle(ratio),0j]),method='DOP853',rtol=rtol,atol=atol,max_step=.15,dense_output=True)
        dense = run.sol
        run.sol = lambda y: np.array([np.exp(1j*dense(y)[0].real),dense(y)[1]])
        run.y[0] = np.exp(1j*run.y[0].real)
    else:
        run = solve_ivp(rhs,(log(end-background.horizon_rho),log(radial_collar)),
                        np.array([ratio,0j]),method='DOP853',rtol=rtol,atol=atol,max_step=.15,dense_output=True)
    if not run.success: raise ArithmeticError(run.message)
    rh = sqrt(1+background.horizon_rho**2);phase = atan2(mass*rh,angular)
    orient = np.diag(np.exp(np.array([-1j,1j])*phase/2))
    frame = orient@horizon_frame(energy,background.surface_gravity,np.hypot(angular,mass*rh),
                                sqrt(2*radial_collar/background.surface_gravity)/rh,
                                background.near_tortoise(radial_collar))
    coefficients = np.linalg.solve(frame,np.array([1.,run.y[0,-1]],complex))
    reflection = coefficients[1]/coefficients[0]
    transmission = (initial_current*exp(-2*float(run.y[1,-1].real))/abs(coefficients[0])**2) if energy>mass else 0.
    defect = abs(abs(reflection)**2+transmission-1)
    return ExteriorJostMode(background,energy,mass,angular,radial_collar,end,run,coefficients,
                            complex(reflection),float(transmission),{
                                'horizon_current':float(defect),
                                'outer_Riccati_residual':float(asymptotic_residual),
                                'closed_outer_current':abs(initial_current) if energy<mass else 0.,
                                'outer_radius':float(end),
                                'outer_expansions':outer_expansions,
                            })
