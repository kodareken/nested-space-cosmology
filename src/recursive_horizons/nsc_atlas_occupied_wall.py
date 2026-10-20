"""A declared worldvolume Dirac state, its Israel wall and a shared dust parent.

This is a new surface-material candidate. The transparent bulk Dirac probe
is unchanged. R_ref=1, c=1, and masses below are geometric G*M_phys/c^2.
The finite occupied-state term uses ell_P^2=G*hbar/c^3 explicitly.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from scipy.integrate import quad
from scipy.optimize import brentq

from .nsc_atlas_parent_observables import DustParent


def filled_shell(level, species=2):
    """Exact positive spectrum n/R, multiplicity 2*n for each two-spinor."""
    if type(level) is not int or level < 0 or type(species) is not int or species < 1:
        raise ValueError("nonnegative integer level and positive integer species required")
    return dict(level=level, species=species,
                number=species*level*(level+1),
                energy_coefficient=species*level*(level+1)*(2*level+1)//3)


@dataclass(frozen=True)
class OccupiedWall:
    """Fixed action, with the state changed only through its spectral level."""

    reference_level: int = 1000
    species: int = 2
    bare_k: float = 1/225
    H_squared: float = 49/625
    reference_gas: float = 7/45
    parent_H0: float = .075

    def __post_init__(self):
        reference = filled_shell(self.reference_level, self.species)
        if (reference["energy_coefficient"] == 0 or
                any(not math.isfinite(x) or x <= 0 for x in (
                    self.bare_k, self.H_squared, self.reference_gas, self.parent_H0))):
            raise ValueError("positive finite action data and occupied reference required")

    @property
    def planck_length_squared(self):
        return self.reference_gas/filled_shell(
            self.reference_level, self.species)["energy_coefficient"]

    def gas_coefficient(self, level):
        return self.planck_length_squared*filled_shell(level, self.species)["energy_coefficient"]

    def matter(self, radius, level):
        if not math.isfinite(radius) or radius <= 0:
            raise ValueError("finite positive wall radius required")
        d = self.gas_coefficient(level)
        density = self.bare_k+d/radius**3
        pressure = -self.bare_k+d/(2*radius**3)
        return dict(four_pi_G_sigma=density, four_pi_G_pressure=pressure,
                    shell_mass=radius**2*density,
                    gas_coefficient=d, s=radius*density,
                    s_prime=self.bare_k-2*d/radius**3,
                    s_second=6*d/radius**4)

    def mass_jet(self, radius, level):
        """Static geometric mass and its first two R derivatives at fixed state."""
        if not math.isfinite(radius) or not 0 < radius < 1/math.sqrt(self.H_squared):
            raise ValueError("static core radius must satisfy 0<HR<1")
        R, h, t, d = radius, self.H_squared, self.bare_k, self.gas_coefficient(level)
        a = math.sqrt(1-h*R*R)
        ap, app = -h*R/a, -h/a**3
        m, mp, mpp = t*R*R+d/R, 2*t*R-d/R**2, 2*t+2*d/R**3
        M = h*R**3/2+a*m-m*m/(2*R)
        Mr = 1.5*h*R*R+ap*m+a*mp-m*mp/R+m*m/(2*R*R)
        Mrr = (3*h*R+app*m+2*ap*mp+a*mpp
               -(mp*mp+m*mpp)/R+2*m*mp/R**2-m*m/R**3)
        return dict(mass=M, derivative=Mr, second_derivative=Mrr,
                    beta_minus=a, beta_plus=a-m/R, shell_mass=m)

    def potential(self, radius, mass, level):
        m = self.matter(radius, level)
        fm, fp = 1-self.H_squared*radius**2, 1-2*mass/radius
        beta_minus = (fm-fp+m["s"]**2)/(2*m["s"])
        return fm-beta_minus**2

    def static_state(self, level):
        """The admitted neighboring filled shells, without a parameter scan."""
        if level not in (self.reference_level-1, self.reference_level, self.reference_level+1):
            raise ValueError("only the three declared neighboring states are admitted")
        R = brentq(lambda r: self.mass_jet(r, level)["derivative"],
                   .97, 1.03, xtol=3e-14)
        j, stress = self.mass_jet(R, level), self.matter(R, level)
        M, a, b = j["mass"], j["beta_minus"], j["beta_plus"]
        if not 0 < b < a:
            raise ValueError("candidate did not retain the ordinary exterior orientation")
        # An independent derivative of the signed rest junction F(R,z)=0.
        frr = (-self.H_squared/a**3
               +2*M/(b*R**3)+M*M/(b**3*R**4)-stress["s_second"])
        fz = 1/(2*a)-1/(2*b)
        Vrr = frr/fz
        Vrr_mass = 2*a*j["second_derivative"]/j["shell_mass"]
        temporal = M/(R*R*b)+self.H_squared*R/a
        temporal_residual = temporal+(b-a)/R-2*stress["four_pi_G_pressure"]
        return dict(**filled_shell(level, self.species), radius=R, mass=M,
                    beta_minus=a, beta_plus=b, **stress,
                    mass_radial_second_derivative=j["second_derivative"],
                    V_second=Vrr, V_second_mass_check=Vrr_mass,
                    radial_omega=math.sqrt(Vrr/2),
                    surface_level_spacing=1/R,
                    radial_to_spectral_frequency_ratio=R*math.sqrt(Vrr/2),
                    angular_junction_residual=a-b-stress["s"],
                    temporal_junction_residual=temporal_residual,
                    potential_residual=self.potential(R, M, level),
                    mass_stationarity_residual=j["derivative"],
                    mass_derivative_wrt_spectral_coefficient=self.planck_length_squared*b/R,
                    central_to_static_infinity_redshift=a/b)

    def central_parent_signal(self, state, *, observer_chi=8., emission_time=0.):
        """Central beacon -> static wall -> compensated join -> fixed FLRW observer."""
        R, M, a, b = (state[k] for k in ("radius", "mass", "beta_minus", "beta_plus"))
        rb0 = (2*M/self.parent_H0**2)**(1/3)
        if not math.isfinite(observer_chi) or observer_chi <= rb0:
            raise ValueError("receiver must lie outside the comoving vacuole")
        parent = DustParent(M, rb0)
        H = math.sqrt(self.H_squared)
        cap_time = math.atanh(H*R)/H
        cap_T_delay = a/b*cap_time
        distance = cap_T_delay+parent.tortoise(rb0)-parent.tortoise(R)+a/b*emission_time
        if distance <= 0:
            raise ValueError("this reception chart requires positive outgoing phase")
        def phase(t):
            return quad(lambda s: 1/(1+parent.state(s)["radial_velocity"]),
                        0., t, epsabs=2e-12, epsrel=2e-12)[0]-distance
        upper = distance*(1+parent.state(0.)["radial_velocity"])
        tb = brentq(phase, 0., upper, xtol=3e-13)
        boundary = parent.state(tb)
        eta_b = 2/self.parent_H0*((1+1.5*self.parent_H0*tb)**(1/3)-1)
        eta_r = eta_b+observer_chi-rb0
        tr = 2/(3*self.parent_H0)*((1+self.parent_H0*eta_r/2)**3-1)
        ar = (1+1.5*self.parent_H0*tr)**(2/3)
        redshift = a/b*(1+boundary["radial_velocity"])*ar/boundary["scale_factor"]
        return dict(initial_join_radius=rb0, parent_H0=self.parent_H0,
                    parent_density=3*self.parent_H0**2/(8*math.pi),
                    observer_chi=observer_chi, central_emission_time=emission_time,
                    cap_Killing_delay=cap_T_delay, join_reception_time=tb,
                    reception_time=tr, one_plus_redshift=redshift,
                    join_radius_at_reception=boundary["radius"],
                    join_scale_factor=boundary["scale_factor"], observer_scale_factor=ar,
                    static_limit_redshift=a/b)


def nested_branch_curvature(radius, mass, H_squared, gas_coefficient):
    """Exact negative rest-junction curvature for beta_minus>0>beta_plus.

    At a static solution F_R=0 this equals V_RR. It is only a necessary
    stability diagnostic; arbitrary inputs need not be equilibria.
    """
    if not (radius > 2*mass > 0 and 0 <= H_squared*radius**2 < 1
            and gas_coefficient >= 0):
        raise ValueError("positive-mass timelike static lapses and nonnegative gas required")
    a = math.sqrt(1-H_squared*radius**2)
    b = math.sqrt(1-2*mass/radius)
    Frr = (-H_squared/a**3-2*mass/(b*radius**3)
           -mass**2/(b**3*radius**4)-6*gas_coefficient/radius**4)
    Fz = 1/(2*a)+1/(2*b)
    return dict(F_RR=Frr, F_z=Fz, necessary_V_second=Frr/Fz)
