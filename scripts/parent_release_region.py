"""Fixed source-lobe region and passive full-parent energy/work readout.

No field, metric, incoming waveform, radius or receiving state is selected
by the evolution. The region is a frozen reference-coordinate measurement
window derived from the ORIGINAL initial carrier's nearest nodal surface.
Coordinate energy-ledger rates use GH time; the reference geodesic clock only
labels the sampled leaf. Tube-proper-unit area flux is reported separately.
No initial-data or evolution function is called.
"""
import math

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq

from recursive_horizons.nsc_atlas_gh_matter_work import matter_work
from recursive_horizons.nsc_atlas_parent_gh import FRAME_CHAR, PAIRS
from recursive_horizons.nsc_atlas_s3_frame import arbitrary_points_values
from recursive_horizons.nsc_atlas_trapped_parent import _first_geometry
from recursive_horizons.nsc_atlas_trapped_surfaces import (
    _angular_basis, _embedding, _unit_frame, round_surface)


def source_region(parent):
    """Largest round-reference ball in the initial positive carrier lobe.

    Set q=a/b and p=(X2+X4)/sqrt2,m=(X2-X4)/sqrt2. On f=0,
    m²=p²+2qx, 2p²=1-x²-2qx-X3². Thus anchor.X is bounded by
    q(x+q)+sqrt((1-q²)/2)*sqrt(1+q²-(x+q)²)-q².
    Cauchy--Schwarz gives (1+q²)/sqrt2-q² globally; equality is feasible
    at x=(sqrt2-1)q,X3=0,p=sqrt((1-q²)/2). This fixes the radius,
    rather than choosing one to obtain positive influx or compactness.
    """
    z = parent.mass_squared/parent.H**2
    a = math.sqrt(12/(parent.kappa*(z+3)))
    b = math.sqrt(72/(parent.kappa*(z+8)))
    q = a/b
    if not 0 < q < 1:
        raise ValueError('Nearest-node formula requires the original0<q<1 source')
    c = math.sqrt((1-q*q)/2)
    anchor = np.array([q, c, 0., c])
    maximum_dot = (1+q*q)/math.sqrt(2)-q*q
    radius = math.acos(maximum_dot)
    x = (math.sqrt(2)-1)*q
    m = math.sqrt(c*c+2*q*x)
    y, w = (c+m)/math.sqrt(2), (c-m)/math.sqrt(2)
    return dict(sampler_arguments=dict(anchor=anchor.tolist(), radius=radius),
        carrier_a=a, carrier_b=b, q=q, closest_node_points=[[x, y, 0., w], [x, w, 0., y]],
        maximum_anchor_dot_on_node=maximum_dot,
        defining_source='f=aX1+bX2X4 at the original prepared slice',
        derivation='Global Cauchy--Schwarz upper bound with an explicitly feasible equality point; rho=acos(max anchor.X on f=0).',
        boundary_role='Fixed INITIAL reference-coordinate ball; not a material surface, horizon or newly prepared receiver.')


def nodal_radius_control(parent):
    """Independent scalar stationary-root check, not a radius-selection fit."""
    definition = source_region(parent); q = definition['q']
    R = math.sqrt(1+q*q); lower, upper = q-R, R-q
    coefficient = math.sqrt((1-q*q)/2)
    def derivative(x):
        return q-coefficient*(x+q)/math.sqrt(max(1e-300, 1+q*q-(x+q)**2))
    root = brentq(derivative, lower+1e-12, upper-1e-12, xtol=2e-14)
    value = q*root+coefficient*math.sqrt(1+q*q-(root+q)**2)
    points = np.asarray(definition['closest_node_points'])
    f = definition['carrier_a']*points[:, 0]+definition['carrier_b']*points[:, 1]*points[:, 3]
    return dict(stationary_root=root, feasible_interval=[lower, upper],
        root_error=abs(root-(math.sqrt(2)-1)*q),
        maximum_dot_error=abs(value-definition['maximum_anchor_dot_on_node']),
        nodal_residual=float(max(abs(f))), unit_sphere_residual=float(max(abs(np.sum(points*points, axis=1)-1))),
        global_bound_reason='The reduced objective has strictlynegative second derivative throughout its feasible interval.')


def _orders(order):
    if isinstance(order, (int, np.integer)) and not isinstance(order, bool):
        result = (int(order), int(order), 2*int(order))
    else:
        result = tuple(order)
    if len(result) != 3 or any(isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 2 for n in result):
        raise ValueError('order is integer>=2 or(radial,polar,azimuth) integer tuple')
    return tuple(map(int, result))


def region_quadrature(anchor, radius, order=12):
    """Unit-S3 volume weights and exact reference boundary embedding jets."""
    anchor = np.asarray(anchor, dtype=float)
    if anchor.shape != (4,) or np.any(~np.isfinite(anchor)) or abs(anchor@anchor-1) > 2e-12:
        raise ValueError('Unit source-anchor quaternion required')
    radius = float(radius)
    if not math.isfinite(radius) or not 0 < radius < math.pi:
        raise ValueError('Reference angular radius must lie in(0,pi)')
    nr, nt, np_ = _orders(order)
    angular = _angular_basis(0, (nt, np_))
    surface = round_surface(anchor, radius)
    rho, boundary, dX, _, radial = _embedding(surface, angular)
    x, w = leggauss(nr); rr = (x+1)*radius/2; wr = w*radius/2
    triad = _unit_frame(anchor[None, :])[0]
    directions = angular['directions']@triad
    volume = (np.cos(rr)[:, None, None]*anchor[None, None, :]
              +np.sin(rr)[:, None, None]*directions[None, :, :]).reshape(-1, 4)
    weights = (wr[:, None]*np.sin(rr[:, None])**2*angular['weights'][None, :]).ravel()
    return dict(order=(nr, nt, np_), volume_points=volume, volume_weights=weights,
                boundary_points=boundary, boundary_dX=dX, boundary_radial=radial,
                angular=angular, actual_boundary_rho=rho)


def analytic_geometry_controls(anchor, radius, frame_scale=1., order=12):
    """Independent flat-ball and round-S3 integrals; no physical state/PDE."""
    rule = region_quadrature(anchor, radius, order)
    r = float(frame_scale)
    if not math.isfinite(r) or r <= 0:
        raise ValueError('Positive reference inverse radius required')
    volume = float(sum(rule['volume_weights']))/r**3
    X = rule['boundary_points']; frame = _unit_frame(X)
    V = np.einsum('nik,nak->nai', frame, rule['boundary_dX'])/r
    gram = np.einsum('nai,nbi->nab', V, V)
    area = float(rule['angular']['weights']@(np.sqrt(np.linalg.det(gram))/np.sin(rule['angular']['theta'])))
    exact_volume = 2*math.pi*(radius-math.sin(2*radius)/2)/r**3
    exact_area = 4*math.pi*math.sin(radius)**2/r**2
    nr = _orders(order)[0]; x, w = leggauss(nr); rr = (x+1)*radius/2
    flat_volume = 4*math.pi*float((w*radius/2)@(rr*rr))
    return dict(round_volume_error=abs(volume-exact_volume), round_area_error=abs(area-exact_area),
                flat_volume_error=abs(flat_volume-4*math.pi*radius**3/3),
                flat_area_error=abs(float(sum(rule['angular']['weights']))*radius**2-4*math.pi*radius**2),
                round_exact_volume=exact_volume, round_exact_area=exact_area)


def normal_boost_flux(rho, normal_pressure, outward_flux, relative_speed):
    """Outward stress flux measured by a tube moving at v relative to n."""
    v = np.asarray(relative_speed, dtype=float)
    if np.any(~np.isfinite(v)) or np.any(abs(v) >= 1):
        raise ValueError('Timelike relative normal speed required')
    return ((1+v*v)*outward_flux-v*(rho+normal_pressure))/(1-v*v)


class RegionSampler:
    """Cached fixed-node first-jet readout; never evaluates a physical RHS."""
    def __init__(self, system, anchor, radius, order=12):
        self.system = system
        self.anchor = np.asarray(anchor, dtype=float).copy(); self.radius = float(radius)
        self.rule = region_quadrature(self.anchor, self.radius, order)
        self.nv = len(self.rule['volume_points']); self.nb = len(self.rule['boundary_points'])
        self.boundary = slice(self.nv, self.nv+self.nb)
        points = np.vstack([self.rule['volume_points'], self.rule['boundary_points'], self.anchor[None, :]])
        values = arbitrary_points_values(points, system.basis.labels)
        self.blocks = {}
        for char in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
            indices = system.sectors[char][0]
            self.blocks[char] = (indices, np.ascontiguousarray(values[:, indices]))
        self.count = len(points)
        self.selected = np.r_[np.arange(10), [20, 21]]
        self.characters = np.vstack([system.chars]+[system.chars[self.selected]*FRAME_CHAR[i] for i in range(3)])
        frame = _unit_frame(self.rule['boundary_points'])
        self.tangents = np.einsum('nik,nak->nai', frame, self.rule['boundary_dX'])/system.r
        self.dF = system.r*np.einsum('nik,nk->ni', frame, self.rule['boundary_radial'])

    def _jets(self, state):
        s = self.system; C = s.unpack(state)
        coefficients = np.vstack([C]+[s.r*(D@C[self.selected].T).T for D in s.basis.generators])
        value = np.empty((self.count, len(coefficients)))
        for char, (indices, Y) in self.blocks.items():
            fields = np.flatnonzero(np.all(self.characters == char, axis=1))
            if len(fields):
                value[:, fields] = Y@coefficients[np.ix_(fields, indices)].T
        g = np.empty((self.count, 4, 4)); v = np.empty_like(g)
        dg = np.empty((self.count, 3, 4, 4)); dq = np.empty((self.count, 3, 2))
        for k, (a, b) in enumerate(PAIRS):
            g[:, a, b] = g[:, b, a] = value[:, k]
            v[:, a, b] = v[:, b, a] = value[:, 10+k]
            for i in range(3):
                dg[:, i, a, b] = dg[:, i, b, a] = value[:, 24+12*i+k]
        for i in range(3):
            dq[:, i] = value[:, 24+12*i+10:24+12*i+12]
        return dict(g=g, v=v, dg=dg, q=value[:, 20:22], p=value[:, 22:24], dq=dq)

    def sample(self, state):
        """Scalar budget row; caller owns coordinate t and time integration."""
        s = self.system; state = np.asarray(state, dtype=float)
        if state.shape != (s.dof+10,) or np.any(~np.isfinite(state)):
            raise ValueError('Actual finite saved ParentGH state required')
        J = self._jets(state); geometry = _first_geometry(J['g'], J['v'], J['dg'], s.r)
        W = s.potential(J['q'])
        budget = matter_work(**J, r=s.r, potential=W)
        h = geometry['h']; hi = np.linalg.inv(h); N = geometry['lapse']; beta = geometry['shift']
        Pi = (J['p']-np.einsum('ni,nia->na', beta, J['dq']))/N[:, None]
        jcov = -np.einsum('na,nia->ni', Pi, J['dq']); jup = np.einsum('nij,nj->ni', hi, jcov)
        scalar_s = s.sstar+J['q'][:, 0]/math.sqrt(s.K)
        normal_s = Pi[:, 0]/math.sqrt(s.K)
        normal_expansion = -budget['extrinsic_trace']
        coordinate_volume_rate = .5*np.einsum('nij,nji->n', hi, J['v'][:, 1:, 1:])
        vol = slice(0, self.nv); boundary = self.boundary
        weight = self.rule['volume_weights']*budget['volume_weight'][vol]/s.r**3
        volume = float(sum(weight)); average = lambda a: float(weight@a[vol]/volume)
        energy = float(weight@budget['rho'][vol]); work = float(weight@budget['work_integrand'][vol])
        vb = np.einsum('nij,nj->ni', hi[boundary], self.dF)
        spatial_norm = np.sqrt(np.einsum('ni,ni->n', self.dF, vb))
        unit_cov = self.dF/spatial_norm[:, None]
        induced = np.einsum('nai,nij,nbj->nab', self.tangents, h[boundary], self.tangents)
        area_weight = self.rule['angular']['weights']*np.sqrt(np.linalg.det(induced))/np.sin(self.rule['angular']['theta'])
        area = float(sum(area_weight)); area_radius = math.sqrt(area/(4*math.pi))
        # Fixed-coordinate tube timelikeness, independent of an incoming mask.
        gb_inverse = np.linalg.inv(J['g'][boundary]); spacetime_dF = np.column_stack([np.zeros(self.nb), self.dF])
        tube_norm = np.einsum('ni,nij,nj->n', spacetime_dF, gb_inverse, spacetime_dF)
        if min(tube_norm) <= 0:
            raise ValueError('Frozen diagnostic tube is no longer timelike; not a physical singularity')
        material = N[boundary]*np.einsum('ni,ni->n', jup[boundary], unit_cov)
        shift = -budget['rho'][boundary]*np.einsum('ni,ni->n', beta[boundary], unit_cov)
        speed = np.einsum('ni,ni->n', beta[boundary], unit_cov)/N[boundary]
        if np.any(abs(speed) >= 1):
            raise ValueError('Frozen tube normal speed is not timelike; no flux classification')
        gamma = 1/np.sqrt(1-speed*speed)
        spatial_unit = vb/spatial_norm[:, None]
        normal4 = np.column_stack([1/N[boundary], -beta[boundary]/N[boundary, None]])
        spatial4 = np.column_stack([np.zeros(self.nb), spatial_unit])
        U = gamma[:, None]*(normal4+speed[:, None]*spatial4)
        nB = gamma[:, None]*(spatial4+speed[:, None]*normal4)
        derivatives = np.concatenate([J['p'][boundary, None, :], J['dq'][boundary]], axis=1)
        Uphi = np.einsum('ni,nia->na', U, derivatives)
        nBphi = np.einsum('ni,nia->na', nB, derivatives)
        tube_flux = -np.sum(Uphi*nBphi, axis=1)
        normal_gradient = np.einsum('ni,nia->na', spatial_unit, J['dq'][boundary])
        normal_pressure = np.sum(Pi[boundary]**2+normal_gradient**2, axis=1)-budget['rho'][boundary]
        boosted = normal_boost_flux(budget['rho'][boundary], normal_pressure, material/N[boundary], speed)
        tube_weighted = N[boundary]/gamma*tube_flux
        supplied = np.einsum('ni,ni->n', budget['energy_flux'][boundary]/budget['volume_weight'][boundary, None], unit_cov)
        outside = float(area_weight@(material+shift))
        vacuum = s.parent.vacuum
        excitation_flux = float(area_weight@(material-(budget['rho'][boundary]-vacuum)*
                                             np.einsum('ni,ni->n', beta[boundary], unit_cov)))
        excitation_work = work+float(weight@(N[vol]*vacuum*budget['extrinsic_trace'][vol]))
        area_rate_density = .5*np.einsum('nab,nba->n', np.linalg.inv(induced),
            np.einsum('nai,nij,nbj->nab', self.tangents, J['v'][boundary, 1:, 1:], self.tangents))
        area_rate = float(area_weight@area_rate_density)
        external_weight = weight*np.exp(-18*scalar_s[vol]); external_volume = float(sum(external_weight))
        external_area_weight = area_weight*np.exp(-12*scalar_s[boundary]); external_area = float(sum(external_area_weight))
        external_area_rate = float(external_area_weight@(area_rate_density-12*J['p'][boundary, 0]/math.sqrt(s.K)))
        gradient = np.einsum('nij,nia,nja->na', hi, J['dq'], J['dq'])
        pressure = .5*np.sum(Pi*Pi, axis=1)-np.sum(gradient, axis=1)/6-W
        mu2 = s.parent.mass_squared*np.exp(-12*J['q'][:, 0]/math.sqrt(s.K))
        loading = .5*(Pi[:, 1]**2+mu2*J['q'][:, 1]**2)*np.exp(6*scalar_s)
        return dict(reference_observer_proper_duration=float(state[s.dof+8]),
            energy=energy, excitation_energy=energy-vacuum*volume, volume=volume,
            rho_mean=energy/volume, pressure_mean=average(pressure), w=average(pressure)/(energy/volume),
            volume_coordinate_rate=float(weight@coordinate_volume_rate[vol]),
            normal_expansion_mean=average(normal_expansion),
            s_mean=average(scalar_s), s_min=float(min(scalar_s)), s_max=float(max(scalar_s)),
            normal_s_rate_mean=average(normal_s), normal_s_rate_min=float(min(normal_s)), normal_s_rate_max=float(max(normal_s)),
            gradient_energy_mean_radion=.5*average(gradient[:, 0]), gradient_energy_mean_carrier=.5*average(gradient[:, 1]),
            normal_kinetic_energy_mean_radion=.5*average(Pi[:, 0]**2), normal_kinetic_energy_mean_carrier=.5*average(Pi[:, 1]**2),
            oscillator_loading_mean=average(loading), oscillator_loading_min=float(min(loading)), oscillator_loading_max=float(max(loading)),
            work=work, excitation_work=excitation_work, outward_GH_flux=outside,
            outward_excitation_GH_flux=excitation_flux,
            outward_material_flux=float(area_weight@material), shift_flux=float(area_weight@shift),
            inward_material_flux_positive=float(area_weight@np.maximum(-material, 0)),
            outward_material_flux_positive=float(area_weight@np.maximum(material, 0)),
            tube_outward_area_flux=float(area_weight@tube_flux),
            tube_outward_area_positive=float(area_weight@np.maximum(tube_flux, 0)),
            tube_outward_area_negative=float(area_weight@np.minimum(tube_flux, 0)),
            tube_outward_GH_weighted_flux=float(area_weight@tube_weighted),
            tube_outward_GH_weighted_positive=float(area_weight@np.maximum(tube_weighted, 0)),
            tube_outward_GH_weighted_negative=float(area_weight@np.minimum(tube_weighted, 0)),
            tube_flux_boost_identity_max=float(max(abs(tube_flux-boosted))),
            boundary_normal_speed_min=float(min(speed)), boundary_normal_speed_max=float(max(speed)),
            tube_flux_role='Outward jB=-sum(Uphi*nBphi); area integral is localtube-proper-unit flux, GHweighted usesN/gamma. Neither is the Eulerian normal-energy ledger or a child-feeding claim.',
            instantaneous_energy_coordinate_rate_from_balance=-outside+work,
            area=area, area_radius=area_radius, area_coordinate_rate=area_rate,
            area_radius_coordinate_rate=area_rate/(8*math.pi*area_radius),
            external_volume=external_volume,
            external_volume_coordinate_rate=float(external_weight@(coordinate_volume_rate[vol]-18*J['p'][vol, 0]/math.sqrt(s.K))),
            external_normal_expansion_mean=float(external_weight@(np.exp(6*scalar_s[vol])*(normal_expansion[vol]-18*normal_s[vol]))/external_volume),
            external_area=external_area, external_area_radius=math.sqrt(external_area/(4*math.pi)),
            external_area_coordinate_rate=external_area_rate,
            minimum_tube_normal_squared=float(min(tube_norm)),
            energy_flux_formula_error=float(max(abs(supplied-material-shift))),
            rates_role='Eulerian energy/volume/work and GHweighted flux rates use GH time; unweighted tube-area flux has its separate role. No anchor-clock energy derivative.',
            geometry_role='Full actual first jets and physical measure in a frozen reference ball; external=exp(-12s)Einstein metric, not a new state or child.')
