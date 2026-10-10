"""Finite local vacuum equivalence for a retained right Schwarzschild exterior."""

import math

import numpy as np

SCHEMA = "nsc-atlas-connection-observer-v1"


def _parameters(mass, observer_radius, cut_radius):
    if not all(math.isfinite(x) for x in (mass, observer_radius, cut_radius)):
        raise ValueError("observer parameters must be finite")
    if not 0 < 2*mass < observer_radius < cut_radius:
        raise ValueError("require 0 < 2M < R_obs < R_cut on the right exterior")


def tortoise(radius, mass):
    if not math.isfinite(mass) or mass <= 0 or not math.isfinite(radius) or radius <= 2*mass:
        raise ValueError("tortoise coordinate requires positive mass and R > 2M")
    return radius+2*mass*math.log(radius/(2*mass)-1)


def protected_window(mass, observer_radius, cut_radius):
    """Open future certificate, conditional on identical vacuum initial data."""
    _parameters(mass, observer_radius, cut_radius)
    # Stable difference avoids cancellation of the arbitrary tortoise origin.
    coordinate = cut_radius-observer_radius+2*mass*math.log1p(
        (cut_radius-observer_radius)/(observer_radius-2*mass))
    lapse = math.sqrt(1-2*mass/observer_radius)
    return {"coordinate_duration": coordinate, "proper_duration": lapse*coordinate,
            "observer_lapse": lapse, "upper_endpoint_included": False}


def require_protected_clock(proper_time, mass, observer_radius, cut_radius):
    window = protected_window(mass, observer_radius, cut_radius)
    if not math.isfinite(proper_time) or not 0 <= proper_time < window["proper_duration"]:
        raise ValueError("clock event is outside the certified open vacuum window")
    return proper_time/window["observer_lapse"]


def require_compatible_slice(record, mass_over_ell=.4, cut_over_ell=4.):
    """Check the supplied body, not the separately conditional gluing operation."""
    if (not math.isfinite(mass_over_ell) or not math.isfinite(cut_over_ell)
            or mass_over_ell <= 0 or cut_over_ell <= 2*mass_over_ell):
        raise ValueError("selected mass and retained cut must be finite and outside 2M")
    end, bridge = record["parent_end"], record["bridge"]
    if (not end["exact_vacuum_overlap_available"] or not bridge["K_is_zero"]
            or end["geometry"] != "exact Schwarzschild AE end"
            or not math.isclose(end["ADM_mass_over_ell"], mass_over_ell, rel_tol=1e-13)
            or not math.isclose(bridge["throat_R_over_ell"], 2*mass_over_ell, rel_tol=1e-13)
            or cut_over_ell <= 2*mass_over_ell):
        raise ValueError("supplied slice does not provide the selected exact vacuum exterior")


def causal_window_record(slice_record, length_scale=1.):
    ell = float(length_scale)
    if not math.isfinite(ell) or ell <= 0:
        raise ValueError("length scale must be finite and positive")
    require_compatible_slice(slice_record)
    M, Robs, Rcut = .4*ell, 1.6*ell, 4*ell
    window = protected_window(M, Robs, Rcut)
    times = []
    for order in (32, 64):
        nodes, weights = np.polynomial.legendre.leggauss(order)
        R = (Rcut+Robs)/2+(Rcut-Robs)/2*nodes
        times.append(float((Rcut-Robs)/2*np.dot(weights, 1/(1-2*M/R))))
    # A constant ingoing Kruskal V=exp((t+r*)/(4M)) reaches the observer
    # at the limiting coordinate time. Compare logarithms to avoid overflow.
    kruskal_error = abs((window["coordinate_duration"]+tortoise(Robs, M))/(4*M)
                        -tortoise(Rcut, M)/(4*M))
    return {"schema": SCHEMA, "status": "FINITE_VACUUM_GEOMETRIC_DEGENERACY",
        "admission": {"owner": "lab/docs/nsc-atlas-connection-observer.md",
            "borrowed": ["gr-qc/9403003v1 Section2.1 Eq11-13", "1309.7591v2 Theorems2.6-2.8"],
            "missing": "finite observer domain in the retained right parent exterior",
            "decides": "local vacuum metric cannot distinguish retained-data-equivalent far assemblies in this window",
            "unlocks": "identify observables or later domains that can break this degeneracy",
            "finite_domain": "static R_obs=4M; initial 2M<R<R_cut=10M; open proper-time window",
            "cpu_seconds_max": 60, "quadrature_nodes": [32, 64],
            "stop": "incompatible vacuum initial data or failed propagation check", "priority_claim": False},
        "inputs": {"M_over_ell": M/ell, "R_obs_over_ell": Robs/ell, "R_cut_over_ell": Rcut/ell,
            "initial_slice": "right Schwarzschild t=0, K=0, 2M<R<R_cut",
            "gluing_requirement": "retain closure through R_cut and an identical open vacuum overlap beyond it",
            "separation_threshold_supplied": False},
        "window": {"coordinate_duration_over_ell": window["coordinate_duration"]/ell,
            "observer_proper_duration_over_ell": window["proper_duration"]/ell,
            "observer_lapse": window["observer_lapse"], "upper_endpoint_included": False,
            "inner_optical_endpoint": "minus infinity", "time_origin": "observer crossing of the initial slice"},
        "independent_checks": {"quadrature_analytic_max_error_over_ell": max(abs(t-window["coordinate_duration"]) for t in times)/ell,
            "quadrature_32_64_error_over_ell": abs(times[0]-times[1])/ell,
            "ingoing_kruskal_log_coordinate_error": kruskal_error,
            "quadrature_comparison_is_rigorous_bound": False},
        "consequence": {"equivalence": "same local vacuum metric in the open protected domain",
            "equivalent_inputs": "identical retained initial h,K modulo isometry",
            "different_inputs_allowed": "outside the retained vacuum initial region",
            "new_physical_mechanism": False},
        "limits": ["conditional on retained data and a smooth local vacuum development",
            "no global shell evolution theorem used", "no predicted arrival at the window boundary",
            "field and quantum state data are not equated", "no sibling signal or cross-universe age"]}
