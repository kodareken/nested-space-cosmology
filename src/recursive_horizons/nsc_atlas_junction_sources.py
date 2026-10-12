"""Finite source-compatibility screen using imported spherical junction laws.

This is an application of Israel matching and shell conservation, not a new
junction law or a child-formation solution.  The exact closed-radiation world
tube is prescribed independently of the tested constant-tension source.
See docs/nsc-atlas-junction-sources.md for the analytic obstruction.
"""

from __future__ import annotations

import math


def _positive(name: str, value: float) -> None:
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")


def _stresses(radius: float, velocity: float, acceleration: float,
              mass: float, beta_minus: float) -> dict[str, float]:
    """Comoving FLRW interior, outward Schwarzschild normal, G=c=1.

    The explicit static Schwarzschild chart is restricted to R>2M.  Its
    exterior proper-time tangent is future directed.  Ktau_minus=0 because
    the comoving FLRW fluid world line is geodesic.
    """
    f = 1.0 - 2.0 * mass / radius
    if not math.isfinite(f) or f <= 0:
        raise ValueError("outside static Schwarzschild chart: require R > 2M")
    beta_plus = math.sqrt(f + velocity * velocity)
    ktau_plus = (acceleration + mass / radius**2) / beta_plus
    sigma = (beta_minus - beta_plus) / (4.0 * math.pi * radius)
    pressure = (ktau_plus + (beta_plus - beta_minus) / radius) / (8.0 * math.pi)
    beta_dot = velocity * ktau_plus
    sigma_dot = (-beta_dot / radius - (beta_minus - beta_plus) * velocity / radius**2) / (4.0 * math.pi)
    return {
        "R": radius, "Rdot": velocity, "Rddot": acceleration,
        "f_plus": f, "beta_minus": beta_minus, "beta_plus": beta_plus,
        "Ktau_minus": 0.0, "Ktau_plus": ktau_plus,
        "Ktheta_minus": beta_minus / radius,
        "Ktheta_plus": beta_plus / radius,
        "sigma": sigma, "pressure": pressure, "sigma_dot": sigma_dot,
        "bulk_energy_flux": 0.0,
        "codazzi_residual": sigma_dot + 2.0 * velocity * (sigma + pressure) / radius,
    }


def radiation_junction(eta: float, *, a_max: float = 1.0,
                       mass: float = 0.4, chi: float = math.pi / 3) -> dict[str, float]:
    """Required stress, not a supplied physical source, on an exact world tube.

    a=A sin(eta), tau=A(1-cos(eta)), R=a sin(chi), 0<eta,chi<pi.
    The retained FLRW region is 0<=chi_child<=chi; chi>pi/2 is a large
    spherical region with decreasing areal radius at its boundary.  The
    normal points toward increasing chi and toward increasing exterior R.
    This local identification does not establish a black-hole child or a
    nonspherical common-parent embedding.
    """
    _positive("a_max", a_max)
    if not math.isfinite(mass) or mass < 0:
        raise ValueError("mass must be finite and nonnegative")
    if not math.isfinite(eta) or not 0 < eta < math.pi:
        raise ValueError("eta must lie strictly between 0 and pi")
    if not math.isfinite(chi) or not 0 < chi < math.pi:
        raise ValueError("chi must lie strictly between 0 and pi")
    a = a_max * math.sin(eta)
    s = math.sin(chi)
    radius = a * s
    velocity = s / math.tan(eta)
    acceleration = -a_max**2 * s / a**3
    result = _stresses(radius, velocity, acceleration, mass, math.cos(chi))
    result.update({
        "eta": eta, "a": a, "tau": a_max * (1.0 - math.cos(eta)),
        "rho_child": 3.0 * a_max**2 / (8.0 * math.pi * a**4),
        "p_child": a_max**2 / (8.0 * math.pi * a**4),
        "m_child_boundary": a_max**2 * s**3 / (2.0 * a),
    })
    return result


def constant_wall_polynomial(radius: float, *, a_max: float, mass: float,
                             chi: float, sigma_wall: float) -> float:
    """Necessary angular condition after squaring; roots alone are not enough.

    With k=4*pi*sigma_wall, the condition is
    k^2 R^4 - 2 cos(chi) k R^3 + 2 M R - A^2 sin(chi)^4 = 0.
    It cannot vanish on any nonzero radius interval, for either sign of k,
    because its constant coefficient is nonzero for A>0 and 0<chi<pi.
    """
    _positive("radius", radius)
    _positive("a_max", a_max)
    if not math.isfinite(mass) or mass < 0:
        raise ValueError("mass must be finite and nonnegative")
    if not math.isfinite(chi) or not 0 < chi < math.pi:
        raise ValueError("chi must lie strictly between 0 and pi")
    if not math.isfinite(sigma_wall):
        raise ValueError("sigma_wall must be finite")
    k = 4.0 * math.pi * sigma_wall
    return k*k*radius**4 - 2.0*math.cos(chi)*k*radius**3 + 2.0*mass*radius - a_max*a_max*math.sin(chi)**4


def _metric_acceleration_check(eta: float, *, a_max: float, mass: float,
                             chi: float, step: float = 1e-5) -> float:
    """Independent Ktau via metric Christoffels and a differentiated tangent.

    Uses tau as the difference coordinate, not d(sqrt(f+Rdot^2))/dtau
    algebra from the production formula.  The five-point stencil is a
    finite integration check, not the proof of incompatibility.
    """
    tau = a_max * (1.0 - math.cos(eta))
    h = step * a_max

    def tangent(t: float) -> tuple[float, float, float]:
        a = math.sqrt(2.0*a_max*t - t*t)
        r = a * math.sin(chi)
        v = (a_max-t) * math.sin(chi) / a
        f = 1.0 - 2.0*mass/r
        return math.sqrt(f+v*v)/f, v, r

    ut, ur, r = tangent(tau)
    vals = [tangent(tau + j*h) for j in (-2, -1, 1, 2)]
    dut = (vals[0][0]-8*vals[1][0]+8*vals[2][0]-vals[3][0])/(12*h)
    dur = (vals[0][1]-8*vals[1][1]+8*vals[2][1]-vals[3][1])/(12*h)
    f = 1.0 - 2.0*mass/r
    fp = 2.0*mass/r**2
    at = dut + fp*ut*ur/f
    ar = dur + f*fp*ut*ut/2.0 - fp*ur*ur/(2.0*f)
    # Covariant unit normal (-Rdot, Tdot), and K^tau_tau=n_mu a^mu.
    return -ur*at + ut*ar


def _independent_checks() -> dict[str, float | str | int]:
    acceleration_errors, conservation_errors, gauss_errors = [], [], []
    for chi in (math.pi/3, 2*math.pi/3):
        for eta in (1.27, 1.35, 1.43):
            q = radiation_junction(eta, chi=chi)
            acceleration_errors.append(abs(_metric_acceleration_check(
                eta, a_max=1.0, mass=0.4, chi=chi) - q["Ktau_plus"]))
            h = 1e-5
            t = q["tau"]
            sigmas = [radiation_junction(math.acos(1-(t+j*h)), chi=chi)["sigma"]
                      for j in (-2, -1, 1, 2)]
            derivative = (sigmas[0]-8*sigmas[1]+8*sigmas[2]-sigmas[3])/(12*h)
            conservation_errors.append(abs(derivative+2*q["Rdot"]*(q["sigma"]+q["pressure"])/q["R"]))
            # Difference of normal Einstein constraints: intrinsic curvature
            # cancels, so [G_nn]/8pi must equal -p_child in these conventions.
            def k_constraint(theta: float, temporal: float) -> float:
                return 2*theta**2 + 4*theta*temporal
            jump = (k_constraint(q["Ktheta_plus"], q["Ktau_plus"])
                    - k_constraint(q["Ktheta_minus"], q["Ktau_minus"]))/ (16*math.pi)
            gauss_errors.append(abs(jump+q["p_child"]))

    dust_sigma, dust_pressure = [], []
    for chi in (0.3, 0.7, 1.0):
        for eta in (2.2, 2.5, 2.8):
            a_max = 2.0
            a = a_max*(1.0-math.cos(eta))/2.0
            s = math.sin(chi)
            mass = a_max*s**3/2.0
            q = _stresses(a*s, s/math.tan(eta/2),
                          -a_max*s/(2*a*a), mass, math.cos(chi))
            dust_sigma.append(abs(q["sigma"]))
            dust_pressure.append(abs(q["pressure"]))
    return {
        "metric_acceleration_max_abs_error": max(acceleration_errors),
        "finite_difference_codazzi_max_abs_error": max(conservation_errors),
        "normal_einstein_jump_max_abs_error": max(gauss_errors),
        "oppenheimer_snyder_max_abs_sigma": max(dust_sigma),
        "oppenheimer_snyder_max_abs_pressure": max(dust_pressure),
        "finite_difference_step_over_A": 1e-5,
        "tested_cap_orientations": 2,
        "scope": "Metric/constraint and dust controls; no stability calculation",
    }


def junction_record() -> dict:
    """Deterministic, finite, JSON-serializable evidence for the atlas producer."""
    a_max, mass, chi = 1.0, 0.4, math.pi/3
    sigma_wall = 1.0/(64.0*math.pi*a_max)
    samples = []
    for i in range(33):
        eta = 1.25 + 0.20*i/32
        q = radiation_junction(eta, a_max=a_max, mass=mass, chi=chi)
        samples.append({
            "eta": eta, "tau_over_A": q["tau"]/a_max,
            "R_over_A": q["R"]/a_max, "Rdot": q["Rdot"],
            "beta_minus": q["beta_minus"], "beta_plus": q["beta_plus"],
            "f_plus": q["f_plus"],
            "A_sigma_required": a_max*q["sigma"],
            "A_pressure_required": a_max*q["pressure"],
            "A_sigma_wall": a_max*sigma_wall,
            "A_pressure_wall": -a_max*sigma_wall,
            "A_eos_residual": a_max*(q["sigma"]+q["pressure"]),
            "angular_residual": q["beta_minus"]-q["beta_plus"]-4*math.pi*sigma_wall*q["R"],
            "polynomial_over_A2": constant_wall_polynomial(q["R"], a_max=a_max,
                mass=mass, chi=chi, sigma_wall=sigma_wall)/a_max**2,
            "A2_codazzi_residual": a_max*a_max*q["codazzi_residual"],
            "bulk_energy_flux": 0.0,
        })
    return {
        "schema_version": 1,
        "status": "INCOMPATIBLE_FIXED_SOURCE",
        "claim_type": "Conditional compatibility obstruction; prior-law application, no novelty claim",
        "admission": {
            "borrowed_result": "Israel surface stress and Codazzi laws: gr-qc/0507063v2 (3),(4),(9),(10),(24); fixed vacuum-wall EOS in 0711.2198v1 section 3 and 0803.0268v2 (4)-(6); existing child-neck Israel note and implementation.",
            "missing_connection": "Whether an independently fixed vacuum wall can support a comoving exact radiation child joined to Schwarzschild over a nonzero time interval.",
            "decision_and_unlock": "Exclude only this source/world-tube pairing; a moving boundary, another bulk, or independently justified different source requires a new admitted candidate.",
            "finite_domain": "33 points eta in [1.25,1.45], A=1, M/A=.4, chi=pi/3; R>2M; shell proper time equals child FLRW proper time; no centre, horizon or turnaround extrapolation.",
            "budget_and_stop": "At most 33 witness samples, six differential checks, nine dust controls and one polynomial argument; less than 60 CPU seconds; no solver, parameter fit or scan. Stop this pairing when incompatible.",
            "verification_and_attribution": "Christoffel acceleration, independent finite-difference conservation, normal Einstein jump, shell-free dust control and symbolic polynomial in tests; bounded source search, no worldwide priority assertion.",
        },
        "theory": {
            "bulk": "4D general relativity; closed radiation FLRW and Schwarzschild vacuum",
            "G": 1, "c": 1, "A": a_max, "M_over_A": mass/a_max,
            "chi": chi, "kA": 1/16,
            "source": "sigma_wall=1/(64*pi*A), P_wall=-sigma_wall, fixed positive tension",
            "source_selection": "Prescribed before evaluation, no calibration point or fitted stress function",
            "normal": "Increasing child chi; beta_minus=cos(chi); exterior beta_plus>0",
            "world_tube": "chi fixed; a=A sin(eta); tau=A(1-cos(eta)); 0<=child_chi<=chi",
        },
        "domain": {"eta": [1.25, 1.45], "sample_count": 33,
                   "tau_over_A": [samples[0]["tau_over_A"], samples[-1]["tau_over_A"]],
                   "R_over_A": [samples[0]["R_over_A"], samples[-1]["R_over_A"]],
                   "exterior_horizon_R_over_A": 0.8},
        "analytic_obstruction": {
            "polynomial": "k^2 R^4 - 2 cos(chi) k R^3 + 2 M R - A^2 sin(chi)^4",
            "constant_coefficient": "-A^2 sin(chi)^4 != 0",
            "argument": "A polynomial zero throughout a varying-radius interval has all coefficients zero; the constant coefficient forbids this for A>0, 0<chi<pi.",
            "scope": "Any finite real constant sigma of either sign and either sign cos(chi), either exterior normal sign; exact comoving radiation boundary and fixed Schwarzschild mass. Squaring is only a necessary condition.",
            "not_covered": "An isolated instant, noncomoving wall, de Sitter/dust child, different source or exterior, physical production, stability, nonspherical embedding",
        },
        "samples": samples,
        "independent_checks": _independent_checks(),
        "source_reading": [
            {"url": "https://arxiv.org/html/gr-qc/0507063v2", "read": "Sections 2, 3.1, 3.2; use (3),(4),(9),(10),(24), general metric definition rather than importing a static interior"},
            {"url": "https://arxiv.org/html/0711.2198v1", "read": "Sections 2-3 and overview of formation constraints; BGG is prior art via this review, original BGG not independently reread"},
            {"url": "https://arxiv.org/html/0803.0268v2", "read": "Sections II, IV, V; domain-wall plus dust source assumptions and proposed origins/exchange; no numerical breathing-bubble result imported", "erratum": "https://journals.aps.org/prd/pdf/10.1103/PhysRevD.80.049901", "erratum_read": "Full one-page 2009 correction: effective-potential equations (9),(15) corrected; neither is used here. Authors state their figures/discussion used corrected expressions."},
        ],
        "limits": [
            "Not a solution for producing a child universe or gluing this child into Brill-Lindquist two-hole data.",
            "Required surface stress is a diagnostic, not a discovered physical substance.",
            "No perturbation stability result; the unsatisfied background is stopped before such a calculation.",
            "No inference about constants, dark matter, dark energy, antimatter or cross-universe ages.",
            "Existing balanced torsion/vacuum layer still lacks its physical exchange mechanism; this screen does not certify it.",
        ],
    }
