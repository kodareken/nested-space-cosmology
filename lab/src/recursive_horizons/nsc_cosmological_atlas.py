"""Finite exterior-first atlas: imported geometry, declared hypotheses and gaps.

An atlas record is not a solved spacetime. Geometric attachment, descent,
causal access and correlation are independent relations. The optional finite
Hamiltonian adapter exercises an existing reducer; it does not source a metric.
All numerical examples below are dimensionless compatibility/software controls.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import heapq
import json
import math

import numpy as np
import sympy as sp

from .nsc_evolving_reduction import evolve_retained_region

SCHEMA = "nsc-cosmological-atlas-v1"
LINK_KINDS = frozenset(("descent", "attachment", "causal", "correlation"))
SOURCES = {
    "barrabes_frolov": {
        "url": "https://arxiv.org/html/hep-th/9511136v1",
        "read": "Introduction and sections 2-4, equations 3, 6, 7, 12-22",
        "use": "Multiple daughters in ONE hole; null-shell matching and collisions",
        "limit": "Kinematic toy model; l>a condition uses near-singularity approximation, equal nucleation radius and a<<R0<<2m; nucleation probability absent",
    },
    "smolin": {
        "url": "https://arxiv.org/html/gr-qc/9404011v1",
        "read": "Introduction and dimensionless-parameter hypothesis",
        "use": "Descent with potentially inherited or changed dimensionless parameters",
        "limit": "Hypothesis; no physical transition kernel imported",
    },
    "poplawski": {
        "url": "https://arxiv.org/html/2007.11556v2",
        "read": "Sections 1-8 through equation 33; v2, 23 June 2026",
        "use": "Spin-fluid source, shear and conditional bounce",
        "limit": "Exact Kantowski-Sachs geometry is not reached; S3 continuation missing",
    },
    "multimouth": {
        "url": "https://arxiv.org/html/2012.07821v2",
        "read": "Sections 1-3.4 and discussion; Appendix A not fully checked",
        "use": "Controlled three-mouth candidate, energy/size limits and signaling",
        "limit": "Charged near-extremal semiclassical domain with many fermion species",
    },
    "causal_review": {
        "url": "https://arxiv.org/html/1801.04539v1",
        "read": "Introduction and causal-structure discussion",
        "use": "Causal diagrams distinguish connection from signal access",
        "limit": "Review, not proof that astrophysical black holes host children",
    },
    "finite_reduction": {
        "path": "lab/src/recursive_horizons/nsc_evolving_reduction.py",
        "read": "evolve_retained_region and fixed-observer interface",
        "use": "Imported finite retarded elimination with an explicitly supplied H(t)",
        "limit": "Finite matrix control; no inferred geometry or relativistic signal cone",
    },
    "nested_qualities": {
        "path": "docs/nsc-nested-qualities.md",
        "read": "Sections 1-4",
        "use": "Known Schur reduction; inherited law differs from regional state",
        "limit": "Finite operators; no determination of physical constants",
    },
    "israel_note": {
        "path": "lab/docs/child-neck/israel-matching-note.md",
        "read": "Angular junction formula and ball versus small-neck orientation",
        "use": "Imported angular thin-shell matching, not complete surface dynamics",
        "limit": "Pressure, conservation and a physical shell source still needed",
    },
    "fixed_theory": {
        "url": "https://arxiv.org/html/2604.27897v1",
        "read": "Sections II.2.1, III matching/domain passages, IV",
        "use": "Hold the effective theory fixed when perturbing a state",
        "limit": "Idealized null shells; not a nonlinear stability proof",
    },
}


def canonical_bytes(value):
    return (json.dumps(value, sort_keys=True, indent=2, ensure_ascii=False,
                       allow_nan=False) + "\n").encode("utf-8")


def fingerprint(value):
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def _positive(value, label):
    if isinstance(value, bool) or not math.isfinite(float(value)) or value <= 0:
        raise ValueError(f"{label} must be finite and positive")
    return float(value)


def schwarzschild_f(radius, geometric_mass):
    """Imported f=1-2mu/R, in G=c=1 units; mu is a length, not kilograms."""
    return 1.0 - 2.0 * _positive(geometric_mass, "geometric mass") / _positive(radius, "radius")


def angular_junction(radius_minus, radius_plus, f_minus, f_plus, radial_speed,
                     *, normal_minus=1, normal_plus=1):
    """Angular Israel constraint only: 4*pi*R*sigma = beta_minus-beta_plus.

    Both R and dR/dtau use the same induced proper-time metric. This does
    not supply the other extrinsic-curvature equation or shell matter.
    """
    radius = _positive(radius_minus, "minus radius")
    other = _positive(radius_plus, "plus radius")
    if not math.isclose(radius, other, rel_tol=1e-12, abs_tol=0):
        raise ValueError("incompatible induced S2 radii at the junction")
    if normal_minus not in (-1, 1) or normal_plus not in (-1, 1):
        raise ValueError("normal orientations must be signed unit choices")
    if not all(math.isfinite(float(x)) for x in (f_minus, f_plus, radial_speed)):
        raise ValueError("junction data must be finite")
    if radial_speed == 0 and (f_minus == 0 or f_plus == 0):
        raise ValueError("stationary horizon is null, not a timelike junction")
    squares = (radial_speed**2 + f_minus, radial_speed**2 + f_plus)
    if min(squares) < 0:
        raise ValueError("junction has no real timelike normal for supplied data")
    beta_minus, beta_plus = (normal_minus * math.sqrt(squares[0]),
                             normal_plus * math.sqrt(squares[1]))
    return {
        "beta_minus": beta_minus, "beta_plus": beta_plus,
        "four_pi_R_sigma": beta_minus - beta_plus,
        "sigma_geometric": (beta_minus - beta_plus) / (4 * math.pi * radius),
        "full_junction_solved": False,
        "missing": ["temporal extrinsic curvature", "shell source and conservation"],
    }


def null_collision_residual(f_a, f_b, f_c, f_d):
    """Borrowed null-shell condition at the collision sphere, not everywhere."""
    if not all(math.isfinite(float(x)) for x in (f_a, f_b, f_c, f_d)):
        raise ValueError("null-shell metric data must be finite")
    return f_a * f_b - f_c * f_d


def static_clock(parent_time, radius, geometric_mass, valid_interval):
    """Proper elapsed time of a static Schwarzschild observer from interval start.

    A static clock ends at its domain boundary. There is no continuation of
    this observer through a horizon and no automatic parent/child clock map.
    """
    start, stop = valid_interval
    if not all(math.isfinite(float(x)) for x in (parent_time, start, stop)) or stop <= start:
        raise ValueError("clock interval must be finite and ordered")
    if not start <= parent_time <= stop:
        raise ValueError("observer has exited its declared clock domain")
    lapse_squared = schwarzschild_f(radius, geometric_mass)
    if lapse_squared <= 0:
        raise ValueError("static observer is unavailable at or inside the horizon")
    return math.sqrt(lapse_squared) * (parent_time - start)


def _region(identifier, title, role, geometry, source, state, clock, parameters,
            *, generation=0, contained_in=None):
    return {"id": identifier, "title": title, "role": role, "geometry": geometry,
            "source": source, "state": state, "clock": clock,
            "dimensionless_parameters": parameters,
            "generation": generation, "contained_in": contained_in,
            "validity": "finite candidate patch; no global spacetime completion asserted"}


def _link(identifier, left, right, kind, source, evidence, *, status="unresolved",
          missing=None, delay=None, interval=None, clock=None):
    return {"id": identifier, "from": left, "to": right, "kind": kind,
            "source_refs": list(source), "evidence": evidence, "status": status,
            "missing": list(missing or []), "minimum_delay": delay,
            "valid_interval": interval, "clock": clock}


def build_atlas():
    """Return four frozen finite candidate families and their explicit gaps."""
    sphere = {"symmetry_orbits": "S2", "spatial_topology": "unspecified",
              "embedding": "none; drawings are patch diagrams"}
    parent = _region("parent", "Unresolved parent cosmos", "parent",
        {**sphere, "chart": None, "metric": "g_parent(x), not supplied",
         "globally_spherical": False},
        "parent stress tensor and field equation not supplied", "not supplied",
        {"kind": "no common clock supplied"}, {})
    regions = [parent]
    for key, mass, position in (("a", 1.0, [-3.0, 0.0]), ("b", 1.5, [3.0, 0.0])):
        regions.append(_region(f"hole_{key}", f"Black-hole interface {key.upper()}", "interface",
            {**sphere, "chart": "static Schwarzschild exterior, separate local chart",
             "metric": "ds2=-f dt2+f^-1 dR2+R2 dOmega2; f=1-2mu/R",
             "domain": "R>2mu", "schematic_parent_position": position,
             "placement_is_solution": False},
            "vacuum Einstein exterior only", "quantum state not supplied",
            {"kind": "local exterior Killing coordinate; no parent synchronization"},
            {"mu_over_reference_length": mass}, contained_in="parent"))
        regions.append(_region(f"child_{key}", f"Potential child {key.upper()}", "child",
            {**sphere, "chart": "closed FLRW ansatz, not solved",
             "spatial_topology": "S3 candidate; not the horizon S2",
             "metric": "ds2=-dτ2+a(τ)^2[dχ2+sin(χ)^2 dOmega2]",
             "areal_radius": "R=a(τ) sin(χ)"},
            "transition and child source not supplied", "may differ from sibling",
            {"kind": "child comoving proper time; interpatch conversion missing"},
            {"inherited_law": "unspecified", "state_is_shared": False}, generation=1))
    regions.append(_region("grandchild", "Potential grandchild", "grandchild",
        {**sphere, "chart": None, "metric": "unresolved; finite depth stops here"},
        "formation source not supplied", "not supplied", {"kind": "not supplied"}, {}, generation=2))
    links = []
    for key in ("a", "b"):
        links.extend([
            _link(f"parent-hole-{key}", "parent", f"hole_{key}", "attachment",
                  ["causal_review"], "A separate local Schwarzschild chart is available",
                  missing=["shared nonspherical parent metric and patch matching"]),
            _link(f"descent-{key}", f"hole_{key}", f"child_{key}", "descent",
                  ["smolin", "poplawski"], "Published black-hole offspring hypothesis class",
                  missing=["formation mechanism", "state and dimensionless-parameter transfer"]),
            _link(f"neck-{key}", f"hole_{key}", f"child_{key}", "attachment",
                  ["israel_note"], "Known junction equations can constrain a proposed neck",
                  missing=["induced metric", "physical layer source", "complete matching"]),
            _link(f"infall-{key}", "parent", f"child_{key}", "causal",
                  ["causal_review"], "Inward causal access is possible in appropriate collapse geometries",
                  missing=["actual collapse geometry", "arrival kernel", "clock conversion"]),
        ])
    links.extend([
        _link("siblings", "child_a", "child_b", "correlation",
              ["nested_qualities"], "Common history can correlate unequal regional states",
              missing=["joint initial state"], status="hypothesis"),
        _link("grandchild-descent", "child_a", "grandchild", "descent",
              ["smolin"], "One further possible black-hole formation; finite cutoff",
              missing=["child A black-hole interface", "transition and state transfer"]),
        _link("grandchild-neck", "child_a", "grandchild", "attachment",
              ["israel_note"], "A further interface requires its own junction",
              missing=["child A black-hole interface", "geometry and layer source"]),
    ])
    default = {
        "id": "separate_holes", "title": "One parent, two different black holes",
        "theory": {"family": "unselected source theory", "couplings": {},
                   "rule": "no coupling inferred from a drawing or retuned between states"},
        "assumptions": ["each distinct black hole may host a child",
                        "one further generation only", "no sibling signaling assumed",
                        "local spherical patches do not make the common parent spherical"],
        "regions": regions, "links": links, "max_descent_depth": 2,
        "solution_status": "candidate atlas with unresolved physical interfaces",
    }
    bubble_regions = [_region("hole", "One Schwarzschild interior", "interface",
        {**sphere, "chart": "Schwarzschild trapped region", "domain": "0<R<2m"},
        "limiting-curvature thin-shell model", "not a specified quantum state",
        {"kind": "interior time; exterior t becomes spatial"}, {"m_over_a": 1000.0})]
    for key in ("a", "b"):
        bubble_regions.append(_region(f"bubble_{key}", f"de Sitter daughter {key.upper()}", "child",
            {**sphere, "chart": "de Sitter patch", "metric": "f=1-R2/a2"},
            "false-vacuum interior plus negative-energy null shells", "not supplied",
            {"kind": "local de Sitter chart"}, {"a_over_reference_length": 1.0}, generation=1))
    bubble_links = []
    for key in ("a", "b"):
        bubble_links += [
            _link(f"bubble-descent-{key}", "hole", f"bubble_{key}", "descent",
                  ["barrabes_frolov"], "Independent nucleation sites in a single interior",
                  status="imported_conditional", missing=["nucleation microphysics and probability"]),
            _link(f"bubble-shell-{key}", "hole", f"bubble_{key}", "attachment",
                  ["barrabes_frolov"], "Equations 3, 6, 7: null-shell angular matching",
                  status="imported_conditional", missing=["physical transition source"]),
        ]
    bubble_links.append(_link("bubble-collision", "bubble_a", "bubble_b", "causal",
        ["barrabes_frolov"], "Section 3, equations 9-12: near-singularity approximation with equal nucleation radius and a<<R0<<2m permits intersections for l<a",
        status="imported_conditional", missing=["choose collision continuation and common clock"]))
    bubbles = {"id": "one_hole_multiple_daughters", "title": "Several daughters inside one hole",
        "theory": {"family": "Barrabes-Frolov limiting-curvature toy model", "couplings": {"a": 1.0}},
        "assumptions": ["fixed limiting de Sitter radius a", "thin null shells", "spherical toy model",
                        "l>a disconnection uses near-singularity approximation and equal nucleation radius with a<<R0<<2m"],
        "regions": bubble_regions, "links": bubble_links, "max_descent_depth": 1,
        "fixed_parameter_bindings": [{"region": f"bubble_{key}", "parameter": "a_over_reference_length",
                                      "theory_key": "a"} for key in ("a", "b")],
        "constraints": {"nucleation": "R0^3=2*m*a^2", "collision": "f_A*f_B=f_C*f_D",
                        "disconnection": "l>a prevents shell intersections in the near-singularity, equal-R0 approximation, a<<R0<<2m"},
        "solution_status": "borrowed kinematic candidate; no NSC formation claim"}
    mouths = []
    for name in ("A", "B", "C"):
        mouths.append(_region(name, f"Mouth {name}", "mouth",
            {**sphere, "chart": "matched near-extremal throat patch",
             "metric": "near-throat AdS2 x S2 with controlled perturbations"},
            "charged fermion Casimir energy; magnetic flux and supporting structure",
            "state of imported construction, not the NSC child state",
            {"kind": "coordinate conversion follows matched geometry"}, {}))
    connected_links = [_link(f"mouth-{left}-{right}", left, right, "causal", ["multimouth"],
        "Section 3 establishes conditional traversability; numeric transit time not assigned",
        status="imported_conditional", missing=["NSC coefficient/state identification", "numerical clock map"])
        for left, right in (("A", "B"), ("B", "A"), ("A", "C"), ("C", "A"), ("B", "C"), ("C", "B"))]
    connected = {"id": "connected_mouths", "title": "A connected three-mouth candidate",
        "theory": {"family": "MMP plus Emparan et al.", "couplings": {"g": "fixed per comparison"}},
        "assumptions": ["near-extremal magnetic mouths", "large enough fermion species count",
                        "small third mouth", "supporting structure and controlled backreaction"],
        "regions": mouths, "links": connected_links, "max_descent_depth": 0,
        "constraints": {"energy": "E_small < |E_binding|", "size": "G*m_small < r_e",
                        "source": "m_small < N_f*Q/(8*r_e)", "topology": "F2, not three independent pairwise handles F3"},
        "solution_status": "published conditional causal candidate; no astrophysical identification"}
    equivalence = {
        "id": "response_equivalence", "title": "Different exteriors, the same local response",
        "theory": {"family": "finite Hermitian operator control", "couplings": {"visible_link": "1/3"}},
        "assumptions": ["hidden mode is exactly uncoupled", "same visible state", "finite dimension 3"],
        "regions": [_region(k, title, "operator_block", {"chart": None, "spatial_topology": "none: finite Hilbert space"},
                    "explicitly supplied matrix", "Gaussian covariance", {"kind": "dimensionless model time"}, {})
                    for k, title in (("local", "Local mode"), ("visible", "Visible exterior"), ("hidden", "Hidden exterior"))],
        "links": [_link("operator-link", "local", "visible", "attachment", ["nested_qualities"],
                        "Declared matrix coupling, not geometric gluing", status="algebra_control")],
        "max_descent_depth": 0, "solution_status": "known algebraic control; no spacetime claim"}
    atlas = {"schema": SCHEMA, "purpose": "exterior-first finite cosmological candidate map",
        "admission": {"category": "operator/domain/state compatibility map",
                      "new_physics_result": False, "borrowed": "matching, causal geometry and finite elimination",
                      "missing_arrow": "physically sourced exterior-to-child interfaces and response"},
        "sources": deepcopy(SOURCES), "scenarios": [default, bubbles, connected, equivalence],
        "local_questions_later": ["dark-sector response after a causal/source map exists",
                                 "matter-antimatter symmetry selected by an actual attachment",
                                 "clock and entropy behavior in a named observer domain"],
        "fitted_observations": []}
    validate_atlas(atlas)
    return atlas


def validate_atlas(atlas):
    if atlas.get("schema") != SCHEMA:
        raise ValueError("unsupported atlas schema")
    for scenario in atlas["scenarios"]:
        ids = [r["id"] for r in scenario["regions"]]
        if not ids or len(ids) != len(set(ids)):
            raise ValueError("region identifiers must be unique and finite")
        depth = scenario["max_descent_depth"]
        if isinstance(depth, bool) or not isinstance(depth, int) or not 0 <= depth <= 8:
            raise ValueError("finite descent cutoff must be an integer between 0 and 8")
        children = {key: [] for key in ids}
        regions = {r["id"]: r for r in scenario["regions"]}
        for region in scenario["regions"]:
            generation = region["generation"]
            if isinstance(generation, bool) or not isinstance(generation, int) or not 0 <= generation <= depth:
                raise ValueError("region generation exceeds its finite cutoff")
            host = region["contained_in"]
            if host is not None:
                if host not in regions or regions[host]["generation"] != generation:
                    raise ValueError("a contained interface must share its host generation")
                children[host].append((region["id"], 0))
        _validate_fixed_bindings(scenario)
        link_ids = set()
        for link in scenario["links"]:
            if link["id"] in link_ids:
                raise ValueError("link identifiers must be unique")
            link_ids.add(link["id"])
            if link["from"] not in children or link["to"] not in children or link["kind"] not in LINK_KINDS:
                raise ValueError("invalid link endpoint or kind")
            if not link["evidence"] or not link["source_refs"] or any(s not in atlas["sources"] for s in link["source_refs"]):
                raise ValueError("every link needs evidence and a known source")
            if link["kind"] == "descent":
                children[link["from"]].append((link["to"], 1))
            if link["minimum_delay"] is not None:
                if link["kind"] != "causal" or not link["clock"]:
                    raise ValueError("a causal delay requires a named common clock")
                _positive(link["minimum_delay"], "causal delay")
                if not link["valid_interval"] or len(link["valid_interval"]) != 2:
                    raise ValueError("a delay requires a finite validity interval")
                start, stop = link["valid_interval"]
                if not all(math.isfinite(float(x)) for x in (start, stop)) or stop <= start:
                    raise ValueError("causal validity interval is invalid")
        def descend(node, ancestors):
            if node in ancestors:
                raise ValueError("descent cannot contain a cycle")
            return max((weight + descend(c, ancestors | {node}) for c, weight in children[node]), default=0)
        if max(descend(node, set()) for node in ids) > depth:
            raise ValueError("descent graph exceeds its finite cutoff")
        for link in scenario["links"]:
            if link["kind"] == "descent" and regions[link["to"]]["generation"] != regions[link["from"]]["generation"] + 1:
                raise ValueError("descent must advance exactly one declared generation")
    return True


def earliest_arrival(scenario, origin, destination, departure, *, clock):
    """Conditional graph bound on declared positive-delay edges, not geodesic solving.

    Correlation/descent/attachment edges never become signal edges. Symbolic
    or unresolved routes have no numerical arrival bound and are skipped.
    Edge intervals cover the whole traversal; waiting is permitted within them.
    """
    ids = {r["id"] for r in scenario["regions"]}
    if origin not in ids or destination not in ids or not math.isfinite(departure):
        raise ValueError("invalid observer or departure time")
    queue, best = [(float(departure), origin)], {origin: float(departure)}
    while queue:
        time, node = heapq.heappop(queue)
        if time != best[node]:
            continue
        if node == destination:
            return time
        for edge in scenario["links"]:
            if edge["from"] != node or edge["kind"] != "causal" or edge["minimum_delay"] is None:
                continue
            if edge["status"] not in ("imported_conditional", "declared_control") or edge["clock"] != clock:
                continue
            if edge["to"] not in ids:
                raise ValueError("causal endpoint is not a known region")
            delay = _positive(edge["minimum_delay"], "causal delay")
            interval = edge["valid_interval"]
            if not isinstance(interval, (list, tuple)) or len(interval) != 2:
                raise ValueError("causal interval must contain two finite endpoints")
            start, stop = interval
            if not all(math.isfinite(float(x)) for x in (start, stop)) or stop <= start:
                raise ValueError("causal interval must be finite and ordered")
            arrival = max(time, start) + delay
            if arrival <= stop and arrival < best.get(edge["to"], math.inf):
                best[edge["to"]] = arrival
                heapq.heappush(queue, (arrival, edge["to"]))
    return None


def _validate_fixed_bindings(scenario):
    regions = {r["id"]: r for r in scenario["regions"]}
    for binding in scenario.get("fixed_parameter_bindings", []):
        value = regions[binding["region"]]["dimensionless_parameters"][binding["parameter"]]
        fixed = scenario["theory"]["couplings"][binding["theory_key"]]
        if value != fixed:
            raise ValueError("regional parameter disagrees with its fixed theory coupling")


def require_fixed_theory(reference, perturbed):
    """Changing a state is permitted; secretly retuning the theory is not."""
    if fingerprint(reference["theory"]) != fingerprint(perturbed["theory"]):
        raise ValueError("theory changed during a fixed-theory comparison")
    if reference.get("fixed_parameter_bindings", []) != perturbed.get("fixed_parameter_bindings", []):
        raise ValueError("fixed theory parameter bindings changed")
    _validate_fixed_bindings(reference)
    _validate_fixed_bindings(perturbed)
    return True


def evolve_supplied_hamiltonian(matrix, times, retained_indices, *, covariance, origin):
    """Finite static/callable H(t) adapter; origin exposes who supplied H and C.

    Uses the existing retarded reducer without a geometry inference, added
    evolution algorithm, or a claim of relativistic finite propagation speed.
    """
    if not isinstance(origin, str) or not origin.strip():
        raise ValueError("explicit Hamiltonian origin is required")
    grid = np.asarray(times, dtype=float)
    if grid.ndim != 1 or not 2 <= len(grid) <= 4097 or not np.isfinite(grid).all():
        raise ValueError("finite time grid must contain 2 to 4097 finite nodes")
    if np.any(np.diff(grid) <= 0):
        raise ValueError("finite time grid must increase")
    def checked(value, dimension=None):
        h = np.array(value, dtype=complex, copy=True)
        if h.ndim != 2 or h.shape[0] != h.shape[1] or not 2 <= h.shape[0] <= 64:
            raise ValueError("finite Hamiltonian dimension must be between 2 and 64")
        if dimension is not None and h.shape != (dimension, dimension):
            raise ValueError("supplied Hamiltonian changed dimension")
        if not np.isfinite(h).all() or not np.allclose(h, h.conj().T, rtol=0, atol=1e-12):
            raise ValueError("supplied Hamiltonian must be finite and Hermitian")
        return h
    h = checked(matrix(float(grid[0])) if callable(matrix) else matrix)
    dimension = len(h)
    hamiltonian = (lambda t: checked(matrix(t), dimension)) if callable(matrix) else (lambda _t: h)
    indices = list(retained_indices)
    if not indices or len(set(indices)) != len(indices) or len(indices) >= len(h):
        raise ValueError("retained indices must be a nonempty proper subset")
    if any(isinstance(i, bool) or not isinstance(i, int) or not 0 <= i < len(h) for i in indices):
        raise ValueError("retained index outside Hamiltonian")
    basis = np.eye(len(h), dtype=complex)[:, indices]
    result = evolve_retained_region(hamiltonian, basis, grid, covariance=covariance)
    result["atlas_scope"] = {"origin": origin, "physical_hamiltonian_derived": False,
                             "relativistic_causality_proved": False,
                             "reused": "nsc_evolving_reduction.evolve_retained_region"}
    return result


def response_equivalence_control():
    """Exact known Schur control: a dark exterior mode is locally invisible."""
    z = sp.Symbol("z")
    local, coupling = sp.Matrix([[1]]), sp.Matrix([[sp.Rational(1, 3), 0]])
    responses = []
    for hidden in (3, 5):
        exterior = sp.diag(2, hidden)
        responses.append(sp.simplify((z * sp.eye(1) - local - coupling *
            (z * sp.eye(2) - exterior).inv() * coupling.T).inv()[0]))
    difference = sp.simplify(responses[0] - responses[1])
    return {"exterior_eigenvalues": [[2, 3], [2, 5]], "local_resolvent": str(responses[0]),
            "exact_difference": str(difference), "equivalent_for_declared_local_probe": difference == 0,
            "assumptions": ["no coupling to the hidden mode", "identical visible state"],
            "resolvent_domain": "Im(z)>0; boundary values exclude full-system eigenvalue poles",
            "additional_discriminator": "A direct hidden-mode spectral probe distinguishes energy 3 from 5; adding a link is a separate model",
            "claim": "imported algebra; distinct exteriors need not be identifiable locally"}


def control_admission(question, missing, unlock, budget, stop):
    """Six-field per-calculation contract; adapter controls claim no new physics."""
    return {"question": question, "borrowed": "published matching or existing finite reduction; see source records",
            "missing_connection": missing, "decision_and_unlock": unlock,
            "finite_budget": budget, "stop": stop}


def compatibility_controls():
    """Cheap imported-formula controls; no cosmological fit or PDE campaign."""
    # Rational data make the null matching exactly checkable without a solver.
    m, a = sp.Integer(4), sp.Integer(1)
    r0 = sp.Integer(2)
    matching = sp.simplify(1 - 2 * m / r0 - (1 - r0**2 / a**2))
    neck = angular_junction(4., 4., .75, .5, 0., normal_minus=-1)
    atlas = build_atlas()
    baseline = atlas["scenarios"][0]
    perturbed = deepcopy(baseline)
    perturbed["regions"][1]["dimensionless_parameters"]["mu_over_reference_length"] = 1.1
    return {"admission": control_admission("Do the imported patch and graph constraints agree with the proposed motif?",
                "Physical interface source, joint state and complete geometry remain absent",
                "Reject incompatible attachments; expose which response connections still need sourcing",
                {"scenarios": 4, "maximum_descent_depth": 2, "symbolic_matrix_dimension": 3},
                "Stop on invalid radius, causal clock, lineage, source binding or changed theory; no PDE continuation"),
            "bubble_matching": {"m": 4, "a": 1, "R0": 2, "exact_f_jump": str(matching),
                                "deep_interior_hierarchy_asserted": False},
            "small_neck_angular_constraint": neck,
            "default_sibling_numeric_arrival": earliest_arrival(baseline, "child_a", "child_b", 0., clock="none"),
            "same_theory_after_state_perturbation": require_fixed_theory(baseline, perturbed),
            "response_equivalence": response_equivalence_control()}


def integration_control():
    """A supplied four-mode motif tests wiring, not a cosmological Hamiltonian."""
    h = np.array([[1., .2, .1, 0.], [.2, 1.5, 0., .15],
                  [.1, 0., 2., 0.], [0., .15, 0., 2.5]])
    covariance = np.diag([.2, .7, .4, .1])
    times = np.linspace(0., .25, 33)
    result = evolve_supplied_hamiltonian(h, times, [1, 2], covariance=covariance,
        origin="explicit dimensionless four-mode adapter fixture; no geometry-derived coefficients")
    return {"admission": control_admission("Does the atlas accept an explicitly supplied H through the existing reducer?",
                "No physically sourced cosmological H has been supplied",
                "Validate the adapter for a future sourced finite response; do not admit matrix fixtures as cosmology",
                {"dimension": 4, "output_nodes": 33, "dimensionless_time_end": .25},
                "Stop after this one bounded fixture; reject invalid H, covariance or clock grid"),
            "scope": result["atlas_scope"], "H": h.tolist(), "C_initial": covariance.tolist(),
            "time_grid": times.tolist(), "retained_indices": [1, 2],
            "final_occupations": np.diag(result["covariance_total"][-1]).real.tolist(),
            "covariance_assembly_residual_max": result["covariance_assembly_residual_max"],
            "trapezoid_residual_max": result["trapezoid_residual_max"]}
