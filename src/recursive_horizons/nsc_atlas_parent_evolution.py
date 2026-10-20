"""Normal-clock evolution of the existing quiet two-wall parent.

All regional masses, vacuum values, walls and the initial slice are imported.
The only new choice is the freely falling normal clock, through t/ell=.01.
Wall-crossing Jacobi impulses are retained; no scalar field grid is evolved.
"""
from __future__ import annotations

import math
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root

from . import nsc_atlas_nested_slice as ns
from .nsc_atlas_geometric_average import KAPPA, scalar_average_sources


class ParentBudget:
    def __init__(self, limit=3_000_000, cpu_seconds=120.):
        self.limit, self.cpu_seconds = limit, cpu_seconds
        self.count, self.started = 0, time.process_time()

    def add(self, count=1):
        self.count += count
        if self.count > self.limit or time.process_time()-self.started > self.cpu_seconds:
            raise RuntimeError("normal-parent calculation allowance exhausted")


def normal_bulk(initial, elapsed, mass, hubble, *, budget=None, rtol=2e-12):
    """R,V,A,B and a static-clock increment; elapsed may differ per ray.

    Columns: R,V,A,B,E,E_l. E is SIGNED Rprime, not initial normal velocity.
    A=dl(t)/dl(0), B=A_t; E_l is the original-label derivative of E.
    The static clock is evaluated only where explicitly requested by a
    nonzero seventh input column; the evolution itself is horizon regular.
    """
    initial = np.atleast_2d(np.asarray(initial, dtype=float))
    if initial.shape[1] not in (6, 7) or not np.isfinite(initial).all():
        raise ValueError("six or seven finite ray columns required")
    elapsed = np.broadcast_to(np.asarray(elapsed, dtype=float), (len(initial),))
    if np.any(elapsed < 0) or np.any(initial[:, [0, 2]] <= 0):
        raise ValueError("positive radius/Jacobi scale and nonnegative time required")
    E = initial[:, 4]
    clock = initial[:, 6] if initial.shape[1] == 7 else np.zeros(len(initial))
    y0 = np.column_stack((initial[:, :4], np.zeros(len(initial)))).T

    def rhs(_, flat):
        if budget is not None:
            budget.add(len(initial))
        R, V, A, B, _ = flat.reshape(5, -1)
        if np.any(R <= 0):
            raise ValueError("normal ray reached the curvature singularity")
        f = 1-2*mass/R-hubble*hubble*R*R
        Tdot = np.divide(clock*E, f, out=np.zeros_like(E), where=clock != 0)
        return (np.array([V, -mass/R**2+hubble*hubble*R, B,
                          (2*mass/R**3+hubble*hubble)*A, Tdot])*elapsed).ravel()

    solution = solve_ivp(rhs, (0., 1.), y0.ravel(), method="DOP853",
                         rtol=rtol, atol=rtol*.02)
    if not solution.success:
        raise RuntimeError(solution.message)
    final = solution.y[:, -1].reshape(5, -1).T
    if np.any(final[:, [0, 2]] <= 0):
        raise ValueError("normal congruence left its positive-radius/Jacobi domain")
    return np.column_stack((final[:, :4], E, initial[:, 5], final[:, 4]))


def normal_moments(rays, mass, hubble):
    """Independent intrinsic curvature from conserved E_l/A, plus sources."""
    R, V, A, B, E, El = np.asarray(rays).T[:6]
    radial, angular = B/A, V/R
    theta = radial+2*angular
    shear = (radial-angular)**2/3
    curvature = -4*El/(A*R)+2*(1-E*E)/R**2
    measure = 4*math.pi*A*R*R
    rho = 3*hubble*hubble/KAPPA
    constraint = curvature+4*radial*angular+2*angular*angular-6*hubble*hubble
    return np.column_stack((measure, measure*theta, measure*theta*theta,
                            measure*shear, measure*curvature,
                            measure*rho, -3*measure*rho)), constraint


class NormalParent:
    """The full centre-to-cut slice, evolved on one normal proper clock.

    Computations use the original ell as unit. The outer boundary follows
    its original normal ray; it is not held at R=4 after time zero.
    """
    endpoint = .01

    def __init__(self, order=64, *, budget=None, rtol=2e-12):
        if order not in (32, 64):
            raise ValueError("admitted Gauss orders are 32 and 64")
        self.order, self.rtol = order, rtol
        self.budget = budget if budget is not None else ParentBudget()
        self.nodes, self.weights = np.polynomial.legendre.leggauss(order)
        self.rb = ns.kottler_horizons(.01, 1.)["rb"]
        self.initial_walls = {name: self.wall(name, R) for name, R
                              in (("inner", .8), ("outer", 1.5))}
        self.wall_solutions = {}
        self.cross_cache = {}
        for name, R0 in (("inner", .8), ("outer", 1.5)):
            def rhs(_, y, name=name):
                self.budget.add()
                w = self.wall(name, y[0])
                clock = (w["gamma"] if name == "inner"
                         else w["beta_plus"]/w["f_plus"])
                return [w["Rdot"], clock]
            sol = solve_ivp(rhs, (0., .04), [R0, 0.], method="DOP853",
                            dense_output=True, rtol=rtol, atol=rtol*.02)
            if not sol.success:
                raise RuntimeError(sol.message)
            self.wall_solutions[name] = sol

    @staticmethod
    def wall(name, R):
        if name == "inner":
            return ns.tension_wall(R, 0., 1., .01, 1., 2.)
        if name == "outer":
            return ns.tension_wall(R, .01, 1., .4, 0., 2.)
        raise ValueError("inner or outer wall required")

    def quad(self, function, low, high):
        if low == high:
            return 0.
        x = (low+high)/2+(high-low)*self.nodes/2
        self.budget.add(len(x))
        return (high-low)/2*self.weights @ function(x)

    def initial_rays(self, name, x):
        x = np.asarray(x, dtype=float)
        if name == "cap":
            R, M, H = x, 0., 1.
            V, Vr, E, jac = R, np.ones_like(R), np.ones_like(R), np.ones_like(R)
        elif name in ("host_bridge", "outer_bridge"):
            M, H, rb = (.01, 1., self.rb) if name == "host_bridge" else (.4, 0., .8)
            s = ns.bridge_chart(x, rb, H, M)
            R, V, Vr, E, jac = (s[k] for k in
                               ("radius", "velocity", "velocity_radius", "Rprime", "dl_dy"))
        else:
            R = x
            if name == "host_left":
                M, H, orientation = .01, 1., -1
                V, Vr = ns.taper_velocity(R, .1, .8, self.initial_walls["inner"]["V_plus"])
            elif name == "host_transition":
                M, H, orientation = .01, 1., 1
                V, Vr = ns.transition_velocity(R, .5, .9, M, H)
            elif name == "host_sheet":
                M, H, orientation = .01, 1., 1
                V, Vr = ns.positive_sheet_velocity(R, M, H)
            elif name == "outer_taper":
                M, H, orientation = .4, 0., -1
                V, Vr = ns.taper_velocity(R, 1., 1.5, self.initial_walls["outer"]["V_plus"])
            else:
                raise ValueError("unknown original region")
            E = ns.signed_slope(R, M, H, V, orientation)
            jac = 1/E
        El = M/R**2-H*H*R+V*Vr
        return np.column_stack((R, V, np.ones_like(R), Vr, E, El)), jac, M, H

    def crossing(self, name, t):
        """Incoming ray and transmitted Jacobi data at a common-time wall."""
        t = float(t)
        if t < 0 or t > self.endpoint*(1+1e-12):
            raise ValueError("parent time outside [0,.01] ell")
        key = name, t
        if key in self.cross_cache:
            return self.cross_cache[key]
        sol = self.wall_solutions[name]
        if name == "inner":
            tau = 0. if t == 0 else brentq(lambda z: sol.sol(z)[1]-t, 0., .04,
                                         xtol=1e-15)
            R = float(sol.sol(tau)[0])
            origin = R*math.exp(-t)
            incoming = np.array([R, R, math.exp(t), math.exp(t), 1., 0., 0.])
        else:
            def incoming_ray(origin):
                data, _, _, _ = self.initial_rays("outer_taper", np.array([origin]))
                data = np.column_stack((data, [1.]))
                ray = normal_bulk(data, t, .4, 0., budget=self.budget, rtol=self.rtol)[0]
                def clock_integrand(x):
                    a, _, _, _ = self.initial_rays("outer_taper", x)
                    return a[:, 1]/(a[:, 4]*(1-.8/x))
                T0 = self.quad(clock_integrand, 1.5, origin)
                return ray, float(T0)

            def residual(z):
                origin, tau = z
                ray, T0 = incoming_ray(origin)
                Rw, Tw = sol.sol(tau)
                return [ray[0]-Rw, T0+ray[6]-Tw]

            w0 = self.initial_walls["outer"]
            if t == 0:
                origin, tau = 1.5, 0.
            else:
                guess = [1.5+w0["Rprime_plus"]*w0["v"]*t, t/w0["gamma"]]
                found = root(residual, guess, tol=2e-11)
                if np.max(np.abs(residual(found.x))) > 3e-12:
                    raise RuntimeError("outer wall/normal clock matching failed")
                origin, tau = map(float, found.x)
            incoming, _ = incoming_ray(origin)
            R = float(sol.sol(tau)[0])
        w = self.wall(name, R)
        side = "minus" if name == "inner" else "plus"
        beta, f = w["beta_"+side], w["f_"+side]
        gamma = (beta*incoming[4]-w["Rdot"]*incoming[1])/f
        b = (w["Rdot"]*incoming[4]-beta*incoming[1])/f
        v = b/gamma
        if not 0 <= tau <= .04 or gamma < 1-1e-12 or (b >= 0 if name == "inner" else b <= 0):
            raise ValueError("crossing left its admitted incoming branch")
        eta_rate = w["worldtube_Ktau_"+side]/gamma-v*incoming[3]/incoming[2]
        gamma_rate = b*eta_rate
        direction = -1 if name == "inner" else 1
        outgoing = incoming[:6].copy()
        outgoing[1] -= 2*R*abs(b)
        outgoing[3] += 2*incoming[2]/abs(b)
        outgoing[4] += direction*2*R*gamma
        outgoing[5] += direction*2*incoming[2]/v*(w["Rdot"]+R*gamma_rate)
        result = dict(t=t, tau=tau, radius=R, original_radius=origin,
                      incoming=incoming, outgoing=outgoing, gamma=gamma, b=b, v=v,
                      label_rate=v/incoming[2], wall=w,
                      normalization_residual=gamma*gamma-b*b-1)
        self.cross_cache[key] = result
        return result

    def snapshot(self, t):
        if not 0 <= t <= self.endpoint:
            raise ValueError("parent time outside the proved finite domain")
        crossings = {name: self.crossing(name, t) for name in ("inner", "outer")}
        worst_constraint, min_A, min_R, max_energy_gap = 0., math.inf, math.inf, 0.

        def integrated(name, lower, upper):
            def integrand(x):
                nonlocal worst_constraint, min_A, min_R, max_energy_gap
                initial, jac, M, H = self.initial_rays(name, x)
                rays = normal_bulk(initial, t, M, H, budget=self.budget, rtol=self.rtol)
                values, constraint = normal_moments(rays, M, H)
                worst_constraint = max(worst_constraint, float(max(abs(constraint))))
                min_A, min_R = min(min_A, min(rays[:, 2])), min(min_R, min(rays[:, 0]))
                gap = rays[:, 1]**2+1-2*M/rays[:, 0]-H*H*rays[:, 0]**2-rays[:, 4]**2
                max_energy_gap = max(max_energy_gap, float(max(abs(gap))))
                return values*jac[:, None]
            return np.asarray(self.quad(integrand, lower, upper))

        pieces = {}
        pieces["cap"] = integrated("cap", 0., crossings["inner"]["original_radius"])
        pieces["host_left"] = integrated("host_left", .8, .1)
        pieces["host_bridge_left"] = integrated("host_bridge", -math.sqrt(.1-self.rb), 0.)
        pieces["host_bridge_right"] = integrated("host_bridge", 0., math.sqrt(.5-self.rb))
        pieces["host_transition"] = integrated("host_transition", .5, .9)
        pieces["host_sheet"] = integrated("host_sheet", .9, 1.5)
        pieces["outer_taper"] = integrated("outer_taper", crossings["outer"]["original_radius"], 1.)
        pieces["outer_bridge_left"] = integrated("outer_bridge", -math.sqrt(.2), 0.)
        pieces["outer_bridge_right"] = integrated("outer_bridge", 0., math.sqrt(3.2))
        for name in ("inner", "outer"):
            if t == 0:
                pieces[name+"_crossed"] = np.zeros(7)
                continue
            times = (self.nodes+1)*t/2
            cuts = [self.crossing(name, z) for z in times]
            initial = np.asarray([c["outgoing"] for c in cuts])
            rays = normal_bulk(initial, t-times, .01, 1., budget=self.budget, rtol=self.rtol)
            values, constraint = normal_moments(rays, .01, 1.)
            worst_constraint = max(worst_constraint, float(max(abs(constraint))))
            min_A = min(min_A, min(rays[:, 2]))
            min_R = min(min_R, min(rays[:, 0]))
            gap = rays[:, 1]**2+1-.02/rays[:, 0]-rays[:, 0]**2-rays[:, 4]**2
            max_energy_gap = max(max_energy_gap, float(max(abs(gap))))
            jac = np.asarray([abs(c["label_rate"]) for c in cuts])
            pieces[name+"_crossed"] = t/2*self.weights @ (values*jac[:, None])
        total = sum(pieces.values())
        walls = {}
        for name, c in crossings.items():
            R, gamma = c["radius"], c["gamma"]
            energy = 2*R*R*gamma
            stress = 2*R*R*(gamma-3/gamma)
            curvature = 2*KAPPA*energy
            total[4:] += [curvature, energy, stress]
            walls[name] = dict(radius=R, proper_time=c["tau"], gamma=gamma, v=c["v"],
                               original_radius=c["original_radius"],
                               energy=energy, stress_trace=stress, curvature=curvature)
        averages = scalar_average_sources(*total)
        return dict(time=t, moments=total.tolist(), sources=averages, walls=walls,
                    pieces={k: v.tolist() for k, v in pieces.items()},
                    controls=dict(max_bulk_constraint=worst_constraint, min_Jacobi_scale=min_A,
                                  min_sample_radius=min_R, max_energy_identity_gap=max_energy_gap),
                    scope="Normal-clock finite scalar response of supplied quiet two-wall geometry.")

    def central_flash(self, t):
        """Pair a simultaneous parent cut with an actual delayed cap signal.

        The original quiet cap's established continuation supplies reception.
        No continuation of the whole normal foliation to reception is assumed.
        """
        from .nsc_atlas_nested_spacetime import central_inward_flash
        c = self.crossing("inner", t)
        w = c["wall"]
        flash = central_inward_flash(radius=w["radius"], H=1.,
                                    beta_minus=w["beta_minus"], Rdot=w["Rdot"])
        return dict(emission_normal_time=t, emission_wall_proper_time=c["tau"],
                    central_reception_time=t+flash["central_proper_time"],
                    one_plus_redshift=flash["one_plus_redshift"],
                    protocol="Negligible inward flash with calibrated wall proper frequency; quiet cap continuation.")
