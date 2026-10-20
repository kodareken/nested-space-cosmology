"""Compact two-centre completion of the checked quiet nested parent.

The shared Schwarzschild bridge has one mass, not two independently chosen
end masses. Right-side reflection changes radial orientation, not time.
All transport uses the existing normal geodesic and Jacobi engine.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.integrate import solve_ivp

from . import nsc_atlas_nested_slice as ns
from .nsc_atlas_geometric_average import KAPPA, partitioned_scalar_average
from .nsc_atlas_parent_evolution import (
    NormalParent, ParentBudget, normal_bulk, normal_moments,
)


class _RegularCapSide(NormalParent):
    """The original one-level body, retained only from centre to throat."""

    def __init__(self, order, *, budget, rtol=2e-12):
        if order not in (32, 64):
            raise ValueError("admitted Gauss orders are 32 and 64")
        self.order, self.rtol, self.budget = order, rtol, budget
        self.nodes, self.weights = np.polynomial.legendre.leggauss(order)
        self.initial_walls = {"outer": self.wall("outer", 1.5)}
        self.cross_cache = {}

        def rhs(_, y):
            budget.add()
            w = self.wall("outer", y[0])
            return [w["Rdot"], w["beta_plus"]/w["f_plus"]]

        sol = solve_ivp(rhs, (0., .04), [1.5, 0.], method="DOP853",
                        dense_output=True, rtol=rtol, atol=rtol*.02)
        if not sol.success:
            raise RuntimeError(sol.message)
        self.wall_solutions = {"outer": sol}

    @staticmethod
    def wall(name, R):
        if name != "outer":
            raise ValueError("one-level side has only its outer membrane")
        return ns.tension_wall(R, 0., 1., .4, 0., 2.)


def reflected_initial_wall(wall):
    """Right wall in the left-to-right global radial orientation.

    Its local minus=cap and plus=connector become global plus=cap and
    minus=connector. U is future in both charts; n and b change sign.
    """
    R, k, gamma, b = (wall[x] for x in ("radius", "tension", "gamma", "gamma_v"))
    return dict(
        Rprime_minus=-wall["Rprime_plus"], Rprime_plus=-wall["Rprime_minus"],
        V_minus=wall["V_plus"], V_plus=wall["V_minus"],
        beta_minus=-wall["beta_plus"], beta_plus=-wall["beta_minus"],
        gamma=gamma, b=-b,
        radial_jump_residual=-wall["Rprime_minus"]+wall["Rprime_plus"]+k*R*gamma,
        velocity_jump_residual=wall["V_minus"]-wall["V_plus"]+k*R*b,
        angular_Israel_residual=(-wall["beta_minus"]+wall["beta_plus"])/R+k,
        temporal_Israel_residual=-wall["worldtube_Ktau_minus"]+wall["worldtube_Ktau_plus"]+k,
    )


def quiet_throat_continuation(mass, hubble, *, order=64, budget=None):
    """Finite lifetime of the E=0 quiet throat and a same-radius bounce test.

    The source bound is necessary at a temporal turn with E=0; it does not
    construct that stress, a bounce solution, or a continuation through R=0.
    """
    if order not in (32, 64):
        raise ValueError("admitted Gauss orders are 32 and 64")
    rb = ns.kottler_horizons(mass, hubble, budget)["rb"]
    a = (hubble*rb)**2
    nodes, weights = np.polynomial.legendre.leggauss(order)
    chi = (nodes+1)*math.pi/4
    u = np.sin(chi)**2
    if budget is not None:
        budget.add(order)
    lifetime = float(math.pi/4*weights @ (2*rb*u/np.sqrt(1-a*(1+u+u*u))))
    lower = math.pi*rb/(2*math.sqrt(1-a))
    upper = math.pi*rb/(2*math.sqrt(1-3*a))
    return dict(throat_radius=rb, normal_lifetime=lifetime,
                analytic_lower_bound=lower, analytic_upper_bound=upper,
                kappa_radial_pressure_bounce_threshold=-1/rb**2,
                supplied_kappa_radial_pressure=-3*hubble*hubble,
                necessary_added_kappa_radial_pressure=3*hubble*hubble-1/rb**2,
                condition="Strict inequalities for a regular temporal bounce at this radius with a strict spatial minimum; radial rho+p must be negative in the Einstein model.",
                scope="Upper limit for a smooth whole normal foliation; not proof that this is its first breakdown or that topology physically changes.")


def _side_snapshot(model, t, *, nested):
    if not 0 <= t <= model.endpoint:
        raise ValueError("compact parent time outside [0,.01] ell")
    names = ("inner", "outer") if nested else ("outer",)
    cuts = {name: model.crossing(name, t) for name in names}
    controls = dict(max_bulk_constraint=0., min_Jacobi_scale=math.inf,
                    max_energy_identity_gap=0.)
    pieces = {}

    def accumulate_controls(rays, constraint, M, H):
        controls["max_bulk_constraint"] = max(controls["max_bulk_constraint"],
                                              float(max(abs(constraint))))
        controls["min_Jacobi_scale"] = min(controls["min_Jacobi_scale"],
                                           float(min(rays[:, 2])))
        gap = rays[:, 1]**2+1-2*M/rays[:, 0]-H*H*rays[:, 0]**2-rays[:, 4]**2
        controls["max_energy_identity_gap"] = max(
            controls["max_energy_identity_gap"], float(max(abs(gap))))

    def original_piece(name, lower, upper):
        def integrand(x):
            engine_name = "host_bridge" if name.startswith("host_bridge_") else name
            initial, jac, M, H = model.initial_rays(engine_name, x)
            rays = normal_bulk(initial, t, M, H, budget=model.budget, rtol=model.rtol)
            values, constraint = normal_moments(rays, M, H)
            accumulate_controls(rays, constraint, M, H)
            return values*jac[:, None]
        pieces[name] = np.asarray(model.quad(integrand, lower, upper))

    original_piece("cap", 0., cuts["inner"]["original_radius"] if nested else 1.5)
    if nested:
        original_piece("host_left", .8, .1)
        original_piece("host_bridge_left", -math.sqrt(.1-model.rb), 0.)
        original_piece("host_bridge_right", 0., math.sqrt(.5-model.rb))
        original_piece("host_transition", .5, .9)
        original_piece("host_sheet", .9, 1.5)
    original_piece("outer_taper", cuts["outer"]["original_radius"], 1.)
    original_piece("outer_bridge", -math.sqrt(.2), 0.)

    for name in names:
        if t == 0:
            pieces[name+"_crossed"] = np.zeros(7)
            continue
        times = (model.nodes+1)*t/2
        crossings = [model.crossing(name, s) for s in times]
        initial = np.asarray([c["outgoing"] for c in crossings])
        M, H = (.01, 1.) if nested else (0., 1.)
        rays = normal_bulk(initial, t-times, M, H, budget=model.budget, rtol=model.rtol)
        values, constraint = normal_moments(rays, M, H)
        accumulate_controls(rays, constraint, M, H)
        label_weight = np.asarray([abs(c["label_rate"]) for c in crossings])
        pieces[name+"_crossed"] = t/2*model.weights @ (values*label_weight[:, None])

    total = sum(pieces.values())
    walls = {}
    for name, c in cuts.items():
        R, gamma = c["radius"], c["gamma"]
        energy = 2*R*R*gamma
        stress = 2*R*R*(gamma-3/gamma)
        curvature = 2*KAPPA*energy
        total[4:] += [curvature, energy, stress]
        walls[name] = dict(radius=R, proper_time=c["tau"], gamma=gamma,
                           local_v=c["v"], energy=energy, stress_trace=stress,
                           curvature=curvature)
    # The vacuum region is the shared Schwarzschild connector, counted
    # once as the sum of its left and right halves.
    connector = pieces["outer_taper"][0]+pieces["outer_bridge"][0]
    return dict(moments=total, pieces=pieces, walls=walls, controls=controls,
                vacuum_region_volume=float(connector),
                curved_region_volume=float(total[0]-connector))


class CompactParent:
    """A nested side and a regular/nested side joined at the same throat."""

    def __init__(self, order=64, *, budget=None, rtol=2e-12):
        self.budget = budget if budget is not None else ParentBudget()
        self.nested = NormalParent(order, budget=self.budget, rtol=rtol)
        self.regular = _RegularCapSide(order, budget=self.budget, rtol=rtol)
        self.cache = {}
        self.clock_origins = {}

    def side(self, name, t):
        if name not in ("nested", "regular"):
            raise ValueError("side must be nested or regular")
        key = name, float(t)
        if key not in self.cache:
            model = getattr(self, name)
            self.cache[key] = _side_snapshot(model, t, nested=name == "nested")
        return self.cache[key]

    def snapshot(self, t, *, right="regular"):
        left, other = self.side("nested", t), self.side(right, t)
        combined = partitioned_scalar_average(
            {"left": left["moments"], "right": other["moments"]})
        volume = combined["moment_vector"][0]
        curved = left["curved_region_volume"]+other["curved_region_volume"]
        wall_trace = KAPPA/(4*volume)*sum(
            w["energy"]-w["stress_trace"] for side in (left, other)
            for w in side["walls"].values())
        return dict(time=t, right_side=right, total_volume=volume,
                    moments=combined["moment_vector"], sources=combined["sources"],
                    side_volume_fractions=combined["volume_fractions"],
                    curved_region_volume_fraction=curved/volume,
                    zero_vacuum_connector_fraction=1-curved/volume,
                    supplied_bulk_trace=3*curved/volume,
                    material_wall_trace=wall_trace,
                    within_side_Q=combined["within_region_Q"],
                    between_side_Q=combined["between_region_Q"],
                    Q_composition_residual=combined["Q_composition_residual"],
                    wall_count=3 if right == "regular" else 4,
                    left_walls=left["walls"], right_walls=other["walls"],
                    controls=dict(
                        max_bulk_constraint=max(s["controls"]["max_bulk_constraint"] for s in (left, other)),
                        min_Jacobi_scale=min(s["controls"]["min_Jacobi_scale"] for s in (left, other))),
                    left_flash=self.nested.central_flash(t),
                    domain="Complete compact S3 slice, two regular centres and no outer cut; inherited normal clock.",
                    limits="H values remain supplied. Global scalar average is not a locally inferred Lambda.")

    def connector_event(self, t, *, side, right="regular"):
        """Both outer walls in ONE bridge-anchored Schwarzschild chart.

        Standard Kruskal normalization: -UV=(R/2M-1)*exp(R/2M).
        The central bridge has T=0 initially. Reflection reverses Killing
        T as well as the local radial label, preserving future orientation.
        """
        if side not in ("left", "right") or right not in ("regular", "nested"):
            raise ValueError("left/right side and regular/nested completion required")
        kind = "nested" if side == "left" else right
        model = getattr(self, kind)
        if kind not in self.clock_origins:
            def integrand(R):
                initial, _, _, _ = model.initial_rays("outer_taper", R)
                return initial[:, 1]/((1-.8/R)*initial[:, 4])
            self.clock_origins[kind] = float(model.quad(integrand, 1., 1.5))
        c = model.crossing("outer", t)
        orientation = 1 if side == "left" else -1
        T = orientation*(self.clock_origins[kind]+model.wall_solutions["outer"].sol(c["tau"])[1])
        R = c["radius"]
        w = c["wall"]
        f = 1-.8/R
        rstar = R+.8*math.log(R/.8-1)
        U = orientation*math.exp(-(T-rstar)/1.6)
        V = -orientation*math.exp((T+rstar)/1.6)
        Tdot, Rdot = orientation*w["beta_plus"]/f, w["Rdot"]
        Udot = -U*(Tdot-Rdot/f)/1.6
        Vdot = V*(Tdot+Rdot/f)/1.6
        return dict(normal_time=t, wall_proper_time=c["tau"], R=R, T=float(T),
                    U=U, V=V, U_tau=Udot, V_tau=Vdot,
                    radius_identity=-U*V-(R/.8-1)*math.exp(R/.8),
                    tangent_norm=-f*Tdot*Tdot+Rdot*Rdot/f)
