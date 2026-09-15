"""Massive horizon/infinity modes on the common PG slice.

This owner extends the matched *mode columns*, not the stored covariance,
onto both radial charts.  The radial integrations select no evolution time.
The geometry-only flow primitives supply the clock change; the massless LLL
mode transform and covariance are never used.

The continuous-energy definition and its evaluated finite sections are kept
distinct from a proof of global spectral completeness or a spatial covariance.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import atan2, exp, log, pi, sqrt

import numpy as np
from scipy.integrate import solve_ivp

from .nsc_common_time_bulk_split import CommonTimeBulkSplit
from .nsc_lorentzian import geometry
from .nsc_paired_horizon_preparation import PairedHorizonSeedMap, source_covariance
from .nsc_pg_lll_preparation import LLLFlowAtlas
from .nsc_transmitting_dirac_domain import I2, S1, S2, S3, MODE_TO_CURRENT
from .nsc_unruh_state import horizon_frame


def open_source_projector(energy, mass):
    """Infinity carries flux only above |E|=m; horizon channels always exist."""
    if energy == 0 or mass < 0 or abs(energy) == mass:
        raise ValueError('nonzero energy away from the massive threshold required')
    return np.diag([1., 1., float(abs(energy) > mass)])


def physical_sewing(energy, mass, reflection, transmission):
    """Current sewing from the owned scattering coefficients, not a selector.

    The last row describes a normalized future infinity observable, up to its
    Jost phase. It certifies current balance only; it is not substituted for
    an exterior mode field or promoted to a physical Cauchy map.
    """
    projection = open_source_projector(energy, mass)
    if not 0 <= transmission <= 1+3e-11 or abs(abs(reflection)**2 + transmission - 1) > 3e-11:
        raise ValueError('owned current-normalized scattering data required')
    if not projection[2, 2] and sqrt(transmission) > 3e-11:
        raise ValueError('closed infinity channel has resolved spurious transmission')
    t = sqrt(transmission)
    inner = np.array([[0, 1, 0], [reflection, 0, t]], complex) @ projection
    outer = np.array([[1, 0, 0], [reflection, 0, t]], complex)
    future = np.vstack((inner, [t, 0, -np.conj(reflection)])) if projection[2, 2] else inner
    future = future @ projection
    return inner, outer @ projection, future, projection


@dataclass
class MassiveModeSection:
    energy: float
    compact_mass: float
    angular: float
    rho_interior: np.ndarray
    interior: np.ndarray
    rho_exterior: np.ndarray
    exterior: np.ndarray
    source_covariance: np.ndarray
    source_projector: np.ndarray
    residuals: dict
    function_evaluations: int
    interior_fundamental: np.ndarray
    exterior_ingoing: np.ndarray

    def spatial_covariance(self):
        raise ValueError('finite sections do not supply global spectral completeness or the continuum contact term')

    def retarded_boundary_value(self, side_x, index_x, side_y, index_y):
        """(H-E-i0)^-1 from the physical outgoing/incoming solutions.

        This is the resolvent convention R, not G_R=-R.  It imports the
        first-order Green jump and applies it to these NSC mode fields.
        """
        if self.energy < 0: raise ValueError('this control uses the directly evaluated positive-frequency fields')
        if side_y == 'interior':
            if side_x == 'exterior': return np.zeros((2,2),complex)
            x = self.rho_interior[index_x];y = self.rho_interior[index_y]
            theta = 1. if x < y else (.5 if x == y else 0.)
            return 1j*theta*self.interior_fundamental[index_x]@self.interior_fundamental[index_y].conj().T
        y = self.rho_exterior[index_y]
        up_y = self.exterior[index_y,:,0];in_y = self.exterior_ingoing[index_y]
        velocity = S2-float(geometry(y)[0])*I2
        coefficients = np.linalg.solve(np.column_stack((up_y,-in_y)), 1j*np.linalg.inv(velocity))
        if side_x == 'interior':
            return np.outer(self.interior_fundamental[index_x,:,1],coefficients[1])
        x = self.rho_exterior[index_x]
        if x > y: return np.outer(self.exterior[index_x,:,0],coefficients[0])
        if x < y: return np.outer(self.exterior_ingoing[index_x],coefficients[1])
        return .5*(np.outer(self.exterior[index_x,:,0],coefficients[0])+np.outer(self.exterior_ingoing[index_x],coefficients[1]))

    def spectral_jump_residual(self):
        """Local Stone normalization; neither finite-grid CAR nor completeness."""
        entries = [('interior',i,self.interior[i]) for i in range(len(self.interior))]
        entries += [('exterior',i,self.exterior[i]) for i in range(len(self.exterior))]
        absolute = 0.;relative = 0.
        for sx,ix,fx in entries:
            for sy,iy,fy in entries:
                discontinuity = (self.retarded_boundary_value(sx,ix,sy,iy)-self.retarded_boundary_value(sy,iy,sx,ix).conj().T)/1j
                modes = fx@fy.conj().T
                error = float(np.linalg.norm(discontinuity-modes))
                absolute = max(absolute,error)
                relative = max(relative,error/max(1.,float(np.linalg.norm(modes))))
        return {'spectral_jump_absolute':absolute,'spectral_jump_relative':relative}


class MassivePGModeResolution:
    """Same H_D and fixed horizon, with a signed massive spectral parameter.

    `preparation` owns C_H, n_in, the matched collar, and the radial operator.
    Reflection/transmission may be read from its authenticated scattering
    artifact; no seed covariance or extra particle sector is an input.
    """

    def __init__(self, preparation: PairedHorizonSeedMap, *, rtol=2e-13, atol=2e-15):
        self.preparation = preparation
        self.atlas = LLLFlowAtlas(preparation.horizon_rho, preparation.surface_gravity)
        self.rtol = rtol
        self.atol = atol

    def clock(self, rho):
        # Only coordinate primitives of the already-owned metric are reused.
        return .5 * (self.atlas.coordinate(1, rho) + self.atlas.coordinate(-1, rho))

    def frame(self, rho):
        beta, bp, metric_a = map(float, geometry(rho))
        if rho == self.preparation.horizon_rho:
            raise ValueError('horizon eigenmodes are distributions; use their two chart limits')
        if metric_a < 0:
            a = sqrt(-metric_a)
            xi = np.arctanh(1 / beta)
            frame = (np.cosh(xi / 2)*I2 + np.sinh(xi / 2)*S2) @ MODE_TO_CURRENT / sqrt(a)
            connection = (-beta*bp*I2 - bp*S2) / (2*a*a)
        else:
            xi = np.arctanh(beta)
            phase = np.diag(np.exp(np.array([-1j, 1j])*pi/4))
            frame = (np.cosh(xi / 2)*I2 + np.sinh(xi / 2)*S2) @ MODE_TO_CURRENT @ phase / metric_a**.25
            connection = (beta*bp*I2 + bp*S2) / (2*metric_a)
        return frame, connection

    def pg_generator(self, rho, energy, mass, angular):
        beta, bp, _ = geometry(rho)
        v = S2 - beta*I2
        potential = CommonTimeBulkSplit.local_potential(1., sqrt(1+rho*rho), mass, angular)
        return np.linalg.solve(v, 1j*(energy*I2-potential) + .5*bp*I2)

    def exterior_generator(self, log_offset, energy, mass, angular):
        offset = exp(log_offset)
        p = self.preparation
        rho = p.horizon_rho + offset
        metric_a = p.background.A_from_offset(offset)
        potential = sqrt(metric_a)*(angular/sqrt(1+rho*rho)*S1 + mass*S2)
        return 1j*offset/metric_a*S3 @ (energy*I2-potential)

    def frame_residual(self, rho, energy, mass, angular):
        frame, connection = self.frame(rho)
        beta, _, metric_a = map(float, geometry(rho))
        if metric_a < 0:
            a = sqrt(-metric_a)
            old = 1j/a*(-mass*S1 + angular/sqrt(1+rho*rho)*S2 - energy/a*S3)
            norm_target = -I2
        else:
            old = self.exterior_generator(log(rho-self.preparation.horizon_rho), energy, mass, angular)/(rho-self.preparation.horizon_rho)
            norm_target = S3
        transformed = connection + 1j*energy*beta/metric_a*I2 + frame@old@np.linalg.inv(frame)
        v = S2-beta*I2
        return {
            'Dirac_frame': float(np.linalg.norm(transformed-self.pg_generator(rho, energy, mass, angular))),
            'frame_current': float(np.linalg.norm(frame.conj().T@v@frame-norm_target)),
        }

    def section(self, energy, mass, angular, interior_points, exterior_points, *, scattering=None, jost=None):
        """Evaluate physical mode columns on both charts; seed C is not accepted.

        Negative energy uses the locked paired sigma3-conjugation map and
        recomputes the source law, never complements a stored seed matrix.
        """
        if energy < 0:
            other = self.section(-energy, mass, -angular, interior_points, exterior_points, scattering=scattering, jost=jost)
            other.energy = energy
            other.angular = angular
            other.interior = np.einsum('ij,njk->nik', S3, other.interior.conj())
            other.exterior = np.einsum('ij,njk->nik', S3, other.exterior.conj())
            other.interior_fundamental = np.einsum('ij,njk->nik', S3, other.interior_fundamental.conj())
            other.exterior_ingoing = np.einsum('ij,nj->ni', S3, other.exterior_ingoing.conj())
            other.source_covariance = source_covariance(energy, self.preparation.surface_gravity, self.preparation.inheritance_ratio, mass)
            return other
        p = self.preparation
        inner_rho = np.asarray(interior_points, float)
        outer_rho = np.asarray(exterior_points, float)
        if (inner_rho.ndim != 1 or outer_rho.ndim != 1 or not len(inner_rho) or not len(outer_rho)
                or not np.isfinite(np.r_[inner_rho, outer_rho]).all()
                or np.any(inner_rho >= p.horizon_rho) or np.any(outer_rho <= p.horizon_rho)):
            raise ValueError('finite samples on both sides of the fixed horizon required')
        if jost is not None:
            if (jost.energy,jost.mass,jost.angular) != (energy,mass,angular):
                raise ValueError('Jost mode labels differ from the requested positive partner')
            if scattering is not None: raise ValueError('choose the mode-bearing Jost owner or an authenticated scattering control, not both')
            scattering = {'reflection':jost.reflection,'transmission':jost.transmission}
        scattering = p.scattering(energy, mass, angular) if scattering is None else scattering
        R = complex(scattering['reflection']); transmission = float(scattering['transmission'])
        sewing, exterior_sewing, future, projection = physical_sewing(energy, mass, R, transmission)
        begin = log(p.collar_delta_q)
        inner_y = np.log(p.background.horizon_q - (pi/2 + np.arctan(inner_rho)))
        if np.min(inner_y) < begin:
            raise ValueError('interior samples must lie beyond the declared matched numerical collar')
        frame_h = p.working_frame(begin, energy, mass, angular)
        run_i = solve_ivp(lambda y,u: (-1j*p.generator(y, energy, mass, angular)@u.reshape(2,2)).ravel(),
                          (begin, float(inner_y.max())), frame_h.ravel(), method='DOP853',
                          rtol=self.rtol, atol=self.atol, max_step=.1, dense_output=True)
        epsilon = p.collar_delta_q  # exterior owner uses radial rho offset
        outer_y = np.log(outer_rho-p.horizon_rho)
        if np.min(outer_y) < log(epsilon):
            raise ValueError('exterior samples must lie beyond the declared radial collar')
        rh = p.horizon_radius; phase = atan2(mass*rh, angular)
        orient = np.diag(np.exp(np.array([-1j,1j])*phase/2))
        ext_h = orient@horizon_frame(energy, p.surface_gravity, np.hypot(angular,mass*rh),
                                    sqrt(2*epsilon/p.surface_gravity)/rh, p.background.near_tortoise(epsilon))
        # Evolve physical combinations directly, avoiding cancellation of two
        # exponentially large basis columns after crossing an evanescent barrier.
        initial_exterior = ext_h@exterior_sewing
        # With a Jost field the up column is evaluated by stable inward
        # propagation, rather than a cancellation in an outward solution.
        if jost is not None:
            initial_exterior[:,0] = 0.
            initial_exterior[:,2] = ext_h[:,1]
        run_e = solve_ivp(lambda y,u: (self.exterior_generator(y,energy,mass,angular)@u.reshape(2,3)).ravel(),
                          (log(epsilon), float(outer_y.max())), initial_exterior.ravel(), method='DOP853',
                          rtol=self.rtol, atol=self.atol, max_step=.1, dense_output=True)
        if not run_i.success or not run_e.success:
            raise ArithmeticError('massive radial mode integration failed')
        inner = []; outer = []; frame_checks = []; interior_current = []; exterior_current = []
        inner_frames = []; outer_incoming = []
        for rho,y in zip(inner_rho,inner_y):
            f,_ = self.frame(rho)
            fundamental = np.exp(1j*energy*self.clock(np.array(rho)))*f@run_i.sol(y).reshape(2,2)
            modes = fundamental@sewing
            inner_frames.append(fundamental)
            inner.append(modes)
            current = modes.conj().T@(S2-float(geometry(rho)[0])*I2)@modes
            interior_current.append(float(np.linalg.norm(current+sewing.conj().T@sewing)))
            frame_checks.append(self.frame_residual(rho,energy,mass,angular))
        target = exterior_sewing.conj().T@S3@exterior_sewing
        for rho,y in zip(outer_rho,outer_y):
            f,_ = self.frame(rho)
            static = run_e.sol(y).reshape(2,3).copy()
            if jost is not None:
                static[:,0] = jost.field(rho)
                incoming = np.exp(1j*energy*self.clock(np.array(rho)))*f@static[:,2]
                static[:,2] *= sqrt(transmission)
            else:
                incoming = (np.exp(1j*energy*self.clock(np.array(rho)))*f@static[:,2]/sqrt(transmission)) if transmission else np.full(2,np.nan+0j)
            outer_incoming.append(incoming)
            modes = np.exp(1j*energy*self.clock(np.array(rho)))*f@static
            outer.append(modes)
            current = modes.conj().T@(S2-float(geometry(rho)[0])*I2)@modes
            exterior_current.append(float(np.linalg.norm(current-target)))
            frame_checks.append(self.frame_residual(rho,energy,mass,angular))
        residuals = {
            'open_source_current_sewing': float(np.linalg.norm(future.conj().T@future-projection)),
            'interior_current': max(interior_current),
            'exterior_current': max(exterior_current),
            'Dirac_frame': max(v['Dirac_frame'] for v in frame_checks),
            'frame_current': max(v['frame_current'] for v in frame_checks),
            'closed_infinity_column': float(max(np.linalg.norm(a[:,2]) for a in inner+outer)) if not projection[2,2] else 0.,
            'closed_channel_roundoff_amplitude': sqrt(transmission) if not projection[2,2] else 0.,
        }
        return MassiveModeSection(energy,mass,angular,inner_rho,np.array(inner),outer_rho,np.array(outer),
                                  source_covariance(energy,p.surface_gravity,p.inheritance_ratio,mass),
                                  projection,residuals,run_i.nfev+run_e.nfev,np.array(inner_frames),np.array(outer_incoming))
