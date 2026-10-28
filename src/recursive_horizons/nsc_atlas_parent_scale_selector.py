"""A physical volume mode coupled to the inherited parent source.

This is a separate, homogeneous effective-action candidate. The ancestor
vacuum sheet is an explicit averaged source, not a constructed localized
higher-dimensional solution. Neither its bulk offset nor its UV charge is
selected by inheritance. No observed cosmological scale enters this module.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from .nsc_atlas_inherited_matter import exact_initial_geometry


def matched_coefficients(d: int, ancestor_steps: int, matter=F(3, 8)):
    """Match the old density and an averaged initial volume equilibrium once.

    A and B are requirements on the new parent action, not UV predictions.
    The integer d and ancestor sheet are declared candidate choices.
    """
    if not isinstance(d, int) or d < 1 or not isinstance(ancestor_steps, int) or ancestor_steps < 1:
        raise ValueError("positive integer dimension and ancestor steps required")
    rho = F(9, 32)
    q = F(3, 2)+ancestor_steps*F(1, 2)
    tension = q*q/8
    B = d*(rho+tension+matter/2)/2
    A = B+rho-tension
    stiffness = -d*(d+2)*rho+d*(d-2)*tension-d*(3*d+4)*matter/4
    return dict(d=d, ancestor_steps=ancestor_steps, ancestor_flux=q,
                vacuum=rho, matter=matter, A=A, B=B, T=tension,
                kinetic=F(d*(d+2), 2), initial_stiffness=stiffness,
                required_UV_bare=A-rho)


@dataclass(frozen=True)
class VolumeParent:
    """Minimal two-ancestor stable initial example, kappa_4=L_unit(0)=1."""
    d: int = 12
    ancestor_steps: int = 2

    @property
    def coefficients(self):
        return matched_coefficients(self.d, self.ancestor_steps)

    @property
    def K(self):
        return float(self.coefficients['kinetic'])

    def potential(self, s, derivative=0):
        c = self.coefficients
        return sum(sign*float(c[key])*(-p)**derivative*math.exp(-p*s)
                   for sign, key, p in ((1, 'A', self.d), (-1, 'B', self.d+2),
                                        (1, 'T', 2*self.d)))

    def matter(self, N, s):
        return .375*math.exp(-3*N-self.d*s/2)

    def hubble(self, N, s, v, E):
        C = float(exact_initial_geometry()['curvature'])
        square = (self.K*v*v/2+self.potential(s)+E)/3-C*math.exp(-2*N)
        if square <= 0:
            raise ValueError("positive expanding homogeneous branch ended")
        return math.sqrt(square)

    def averaged_rhs(self, t, y):
        N, s, v = y
        E = self.matter(N, s)
        H = self.hubble(N, s, v, E)
        force = -self.potential(s, 1)+self.d*E/2
        return np.array([H, v, -3*H*v+force/self.K])

    def scalar_rhs(self, t, y, omega0):
        N, s, v, chi, pi = y
        w2 = omega0**2*math.exp(-self.d*s)
        E = (pi*pi+w2*chi*chi)/2
        H = self.hubble(N, s, v, E)
        force = -self.potential(s, 1)+self.d*w2*chi*chi/2
        return np.array([H, v, -3*H*v+force/self.K, pi, -3*H*pi-w2*chi])

    def vacuum_roots(self):
        """All positive extrema via the single turning point of a polynomial.

        For d=12 set x=exp(-2s); V=x^6(A-B*x+T*x^6).
        The derivative polynomial has one positive turning point, so the
        two brackets exhaust positive roots. No spatial or parameter scan.
        """
        if self.d != 12:
            raise ValueError("analytic two-root brackets currently own d=12")
        c = self.coefficients
        A, B, T = (float(c[k]) for k in ('A', 'B', 'T'))
        f = lambda x: 6*A-7*B*x+12*T*x**6
        turn = (7*B/(72*T))**.2
        if f(turn) >= 0:
            return []
        upper = 2*turn
        while f(upper) <= 0:
            upper *= 2
        roots = (brentq(f, 0., turn, xtol=5e-15),
                 brentq(f, turn, upper, xtol=5e-15))
        rows = []
        for x in roots:
            s = -math.log(x)/2
            V, ss = self.potential(s), self.potential(s, 2)
            rows.append(dict(s=s, vacuum=V, stiffness=ss,
                canonical_mass_squared=ss/self.K,
                mass_squared_over_H_vacuum_squared=3*ss/(self.K*V),
                volume_ratio=math.exp(self.d*s),
                charge_ratio=math.exp(-self.d*s/2),
                kind='minimum' if ss > 0 else 'maximum'))
        return rows


def reduction_record():
    return dict(metric='ds_D^2=exp(-d*s) ds_4E^2+exp(2*s) gamma_ab dy^a dy^b',
        assumptions=['D=4+d Einstein gravity, canonical bulk scalar and a D-form',
                     'charged (D-2)-branes wrap the complete internal Einstein space',
                     'one homogeneous volume mode; no claim of shape or KK stability',
                     'ancestor 4D vacuum sheet is an additional averaged source'],
        kappa_4='kappa_D/V0', kinetic='d*(d+2)/(2*kappa_4)',
        three_form_kinetic='Z3=exp(d*s)/V0',
        canonical_charge='e_D*sqrt(V0)*exp(-d*s/2)',
        wrapped_tension='sigma_D*V0*exp(-d*s/2)',
        scalar_mass_and_mixing='(m_D,mu_D)*exp(-d*s/2)',
        inherited_constant_s_invariants=['q/e', 'epsilon', 'sigma/e', 'omega/e'],
        scalar_minimum='-mu0*q0/omega0^2, independent of s',
        gauge_rule='Keep the topological charge and unrescaled A3 in the action; a dynamic canonical field redefinition carries derivative terms.',
        UV_dimensionless_charge='eta_D=e_D*kappa_D^(D/(2*(D-2)))',
        physical_internal_volume='v=V0*exp(d*s)/kappa_D^(d/(D-2))',
        Planck_charge='eta_4=eta_D/sqrt(v)',
        SI_scale_selected=False)


def evolve(*, model=None, endpoint=4., clock='expansion', omega=None,
           rtol=2e-10, atol=2e-12):
    """One finite homogeneous evolution with energy/action and clock controls.

    omega=None is the averaged source; finite omega resolves the exact
    scalar. Full scalar controls are short, not 10^9-cycle integrations.
    """
    model = model or VolumeParent()
    if endpoint <= 0 or endpoint > 6 or clock not in ('expansion', 'proper'):
        raise ValueError("finite 0<N<=6 and a declared clock required")
    if omega is not None and (omega <= 0 or endpoint > .15):
        raise ValueError("resolved scalar controls require positive omega and N<=.15")
    rhs = model.averaged_rhs if omega is None else lambda t, y: model.scalar_rhs(t, y, omega)
    y0 = [0., 0., 0.] if omega is None else [0., 0., 0., math.sqrt(.75)/omega, 0.]
    # Integral ledger components: expansion work; corrected-action remainder;
    # integral |c_dot|/(1-|c|) for a finite adiabatic enclosure.
    def augmented(t, z):
        y = z[:-3]
        dy = rhs(t, y)
        N, s, v = y[:3]
        H = dy[0]
        if omega is None:
            work = 3*H*(model.K*v*v+model.matter(N, s))
            extra = [work, 0., 0.]
        else:
            chi, pi = y[3:]
            w = omega*math.exp(-model.d*s/2)
            Hdot = float(exact_initial_geometry()['curvature'])*math.exp(-2*N)-(model.K*v*v+pi*pi)/2
            c = (3*H-model.d*v/2)/(2*w)
            cdot = (3*Hdot-model.d*dy[2]/2)/(2*w)+c*model.d*v/2
            if abs(c) >= 1:
                raise ValueError("corrected oscillator action not positive")
            S = math.exp(3*N)*chi*pi
            extra = [3*H*(model.K*v*v+pi*pi), cdot*S, abs(cdot)/(1-abs(c))]
        return np.r_[dy, extra]
    initial = np.r_[y0, [0., 0., 0.]]
    if clock == 'expansion':
        def f(N, z):
            dz = augmented(0., z)
            return dz/dz[0]
        solution = solve_ivp(f, (0., endpoint), initial, method='DOP853',
                             rtol=rtol, atol=atol, dense_output=True)
    else:
        def event(t, z):
            return z[0]-endpoint
        event.terminal, event.direction = True, 1
        solution = solve_ivp(augmented, (0., 100.), initial, method='DOP853',
                             rtol=rtol, atol=atol, events=event, dense_output=True)
        if not len(solution.t_events[0]):
            raise RuntimeError("finite proper-clock endpoint not reached")
    if not solution.success or solution.nfev > 180000:
        raise RuntimeError("finite solve failed or RHS budget exhausted")
    y = solution.y[:, -1]
    N, s, v = y[:3]
    E = model.matter(N, s) if omega is None else (y[4]**2+omega**2*math.exp(-model.d*s)*y[3]**2)/2
    energy = model.K*v*v/2+model.potential(s)+E
    record = dict(clock=clock, endpoint_N=N, s=s, v=v,
        internal_volume_ratio=math.exp(model.d*s),
        charge_tension_frequency_ratio=math.exp(-model.d*s/2),
        vacuum_density=model.potential(s), matter_density=E,
        volume_kinetic_density=model.K*v*v/2,
        H=model.hubble(N, s, v, E), energy_ledger_residual=energy-21/32+y[-3],
        nfev=solution.nfev)
    if omega is not None:
        w = omega*math.exp(-model.d*s/2)
        D = math.exp(3*N)*E/w
        c = (3*record['H']-model.d*v/2)/(2*w)
        I = D+c*math.exp(3*N)*y[3]*y[4]
        I0 = .375/omega
        record.update(omega0=omega, corrected_action_ledger_relative=(I-I0-y[-2])/I0,
            action_ratio=D/I0, corrected_action_remainder_bound=math.expm1(y[-1]),
            action_ratio_enclosure=[math.exp(-y[-1])/(1+abs(c)),
                                    math.exp(y[-1])/(1-abs(c))])
    return record, solution


def calculation():
    model = VolumeParent()
    c = model.coefficients
    averaged, _ = evolve(model=model)
    proper, _ = evolve(model=model, clock='proper')
    controls = []
    for omega in (200., 400.):
        exact, _ = evolve(model=model, endpoint=.1, omega=omega, atol=2e-14)
        mean, _ = evolve(model=model, endpoint=.1)
        exact['absolute_volume_field_difference_from_average'] = abs(exact['s']-mean['s'])
        controls.append(exact)
    vacua = model.vacuum_roots()
    stable = next(row for row in vacua if row['kind'] == 'minimum')
    barrier = next(row for row in vacua if row['kind'] == 'maximum')
    C = float(exact_initial_geometry()['curvature'])
    energy = sum(averaged[k] for k in ('vacuum_density','matter_density','volume_kinetic_density'))
    return dict(schema='nsc-atlas-parent-volume-selector-v1', reduction=reduction_record(),
        minimal_vacuum=dict(mass_squared_over_H_vacuum_squared=-6,
            obstruction='Bulk vacuum plus internal curvature cannot stabilize its matched positive vacuum.'),
        extra_stiffness='S_d[W]=W_ss+2*(d+1)*W_s+d*(d+2)*W; a power p contributes (p-d)*(p-d-2)*W.',
        candidate={key: str(val) if isinstance(val, F) else val for key, val in c.items()},
        one_ancestor_initial_obstruction='32*Vtot_ss/d=-2*d-62<0 for every d>0',
        two_ancestor_initial_condition='32*Vtot_ss/d=7*d-80; first positive integer d=12',
        required_wall_volume_normal_derivative_jump=-model.d*.2/(2*model.K),
        wall_matching_obstruction='K*[normal derivative s]=d(sigma)/ds=-d*sigma/2; a globally uniform s with zero momentum cannot satisfy this new junction.',
        source_matching='Same initial charge, flux, scalar state, tension and total Friedmann density as 485b3b9 in the 4D homogeneous sector; volume is a new evolving field.',
        future_is_a_new_candidate=True, vacuum_extrema=vacua,
        averaged_evolution=averaged, independent_proper_clock=proper,
        volume_field_clock_difference=abs(averaged['s']-proper['s']),
        resolved_scalar_controls=controls,
        averaged_basin_entry=dict(total_density=energy,
            barrier_energy_margin=barrier['vacuum']-energy,
            expanding_branch_H_squared_floor=stable['vacuum']/3-C*math.exp(-8),
            left_of_barrier=bool(averaged['s'] < barrier['s']),
            implication='Conditional on this finite numerical entry, the exact averaged energy law confines s and ensures future approach to the positive minimum. Not a full-scalar or localized uplift theorem.'),
        late_vacuum_Lref0_squared=stable['vacuum']/C,
        initial_vacuum_Lref0_squared=float(c['vacuum'])/C,
        remaining_requirements=['derive UV charge eta_D and reference internal volume, not set them to one',
            'derive required bulk offset A-rho and internal curvature, not confuse matching with selection',
            'construct localized ancestor source and test other internal modes',
            'solve the new volume junction and initial constraints, then establish the shared finite domain',
            'new observer comparison: old P004/P005 and likelihoods do not validate this changed dynamics'],
        SI_scale_selected=False, Lambda_obs_derived=False)
