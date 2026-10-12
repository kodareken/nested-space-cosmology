"""Classical Einstein–Maxwell reference for the ingoing RN radial chart.

The action is the laboratory bulk term

    S = ∫ √|g| [-A R_L - C_F F^2],

signature +---, with no Weyl-squared pole and no cosmological term. The
executable flux is the q of the spherical reduction, so

    G_N = 1/(16 π A),    magnetic_r2 = C_F flux^2 / (4 A).

Neutral massless Dirac matter is one closed angular shell, κ = 1, with
multiplicity 4κ inserted once. The angular coupling is κ/r. It is not a
4D mass, and the shell carries no electric current. The magnetic flux is
fixed.

The vacuum chart is the ingoing flat-spatial Painlevé–Gullstrand line element

    ds² = N² dt² - (dr + β dt)² - r² dΩ².

For the source-free magnetic solution the Einstein equation plus unit
asymptotic Killing normalization gives N = 1 and

    β = √(2 mass/r - magnetic_r2/r²),    f = N² - β².

N = 1 is that vacuum result. Coupled evolution does not keep N = 1.
Curvature below uses

    R^a_bcd = ∂_c Γ^a_db - ∂_d Γ^a_cb + Γ^a_ce Γ^e_db - Γ^a_de Γ^e_cb,
    Ricci_bd = R^a_bad,

the same convention as the laboratory tidal consumer. Direct contraction
of this RN metric gives R = 0, Ricci² = 4 magnetic_r2²/r^8 and

    K = 48 mass²/r^6 - 96 mass magnetic_r2/r^7 + 56 magnetic_r2²/r^8.

Constructors and diagnostics do not write. ``write_assessment`` can create
only ``results/development/nsc-rn-reference-assessment-v1.json`` after a
frozen commit, and this slice does not call it.
"""
from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass
from pathlib import Path

for _thread_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_var, "1")

import numpy as np

SCHEMA = "NSC-RN-REFERENCE-v1"
ASSESSMENT_SCHEMA = "NSC-RN-REFERENCE-ASSESSMENT-v1"
LAB = Path(__file__).resolve().parents[2]
HISTORICAL_PERIOD = 8.0
KAPPA = 1
MULTIPLICITY = 4 * KAPPA
ELECTRIC_CURRENT = 0.0
AUDITED_A = 0.04501936182826115
AUDITED_C_F = 0.01125484045706527
AUDITED_FLUX = 4.0
AUDITED_V4 = 0.0
CANONICAL_PI = "canonical_pi"
MAX_ASSESSMENT_BYTES = 64 * 1024 * 1024
ASSESSMENT_OUTPUT = LAB / "results/development/nsc-rn-reference-assessment-v1.json"
SEALED_INPUT_DIRECTORIES = (
    LAB / "results/development/nsc-discovery-parent-strong-empty-k-v1",
    LAB / "results/development/nsc-discovery-parent-empty-source-v1",
    LAB / "results/development/nsc-discovery-parent-cut-confirmation-v1",
    LAB / "results/development/nsc-discovery-vacuum-control-v1",
)
CUT_CONFIRMATION_COMMIT = "5cb775a751add5f72eb06cb681ca16c24a3527d2"
RIEMANN_CONVENTION = (
    "R^a_bcd = d_c Gamma^a_db - d_d Gamma^a_cb "
    "+ Gamma^a_ce Gamma^e_db - Gamma^a_de Gamma^e_cb; "
    "Ricci_bd = R^a_bad; signature (+---)"
)
SIGMA2 = np.diag([1.0, -1.0])
SIGMA_ANGULAR = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)

ACTION = {
    "schema": SCHEMA,
    "name": "classical-leading-Einstein-Maxwell-plus-neutral-massless-Dirac",
    "density": "-A R_L - C_F F^2",
    "signature": "+---",
    "weyl_pole": False,
    "C_W_included": False,
    "cosmological_term": False,
    "V4": 0.0,
    "branch": "fixed-charge",
    "G_N": "1/(16*pi*A)",
    "P2": "C_F*flux**2/(4*A)",
    "magnetic_r2": "C_F*flux**2/(4*A)",
    "maxwell_coupling": "g4**2 = 1/(4*C_F)",
    "flux_convention": (
        "executable flux is q in r_Q^2 = C_F q^2/(4A); "
        "the prose label 2*pi*q is not this parameter"
    ),
    "chart": "ingoing flat-spatial Painleve-Gullstrand",
    "line_element": "N**2 dt**2 - (dr + beta dt)**2 - r**2 dOmega**2",
    "periodic_autonomous_exchange": False,
    "outer_boundary": "nonperiodic absorbing",
    "excision": "pure outflow on an open interval between r_minus and r_plus",
    "parent": "large-r exterior",
    "child": "near the outer horizon, different areal interval",
    "lapse": (
        "vacuum N=1 is derived by asymptotic Killing normalization; "
        "coupled evolution must not keep N=1"
    ),
    "riemann": RIEMANN_CONVENTION,
    "dirac": {
        "neutral": True,
        "massless": True,
        "electric_current": ELECTRIC_CURRENT,
        "magnetic_flux": "fixed",
        "kappa": KAPPA,
        "multiplicity": MULTIPLICITY,
        "multiplicity_counted_once": True,
        "angular_operator": "kappa/r",
        "angular_is_4d_mass": False,
        "exterior_killing_frequency": "supplied by the scattering-packet preparation, not by this reference",
        "phi_shape": ["spinor", "points", "rank"],
        "spinor_dimension": 2,
        "characteristic_sigma2": "diagonal",
        "quadrature": "radial diagonal-norm SBP weights",
        "manufactured_quadrature_control": "positivity_unestablished",
    },
    "k_is_mass": False,
    "bare_misner_sharp_difference": "magnetic_r2/(2*r)",
}


def _positive(value, name):
    number = float(value)
    if not np.isfinite(number) or number <= 0.0:
        raise ValueError(name + " must be finite and positive")
    return number


def _finite(value, name):
    number = float(value)
    if not np.isfinite(number):
        raise ValueError(name + " must be finite")
    return number


def newton_constant(A):
    """G_N = 1/(16 π A) from S = ∫ √|g| [-A R_L + …]."""
    return 1.0 / (16.0 * np.pi * _positive(A, "A"))


def magnetic_radius_square(A, C_F, flux):
    """Geometric charge square C_F flux²/(4A). The audited coefficients give 1."""
    return _positive(C_F, "C_F") * _finite(flux, "flux") ** 2 / (4.0 * _positive(A, "A"))


def angular_coupling(radius, kappa=KAPPA):
    """Sphere Dirac coupling κ/r. This is not a 4D mass."""
    radius = np.asarray(radius, dtype=float)
    if radius.size == 0 or not np.isfinite(radius).all() or np.min(radius) <= 0.0:
        raise ValueError("areal radius must be finite and positive")
    if int(kappa) != KAPPA:
        raise ValueError("the first shell is the closed κ=1 multiplet")
    return kappa / radius


def radial_sbp_weights(r_inner, r_outer, points):
    """Diagonal-norm SBP weights on a uniform nonperiodic radial grid.

    Interior weights equal the spacing and the two endpoints carry half that
    spacing, so the weights sum to the areal interval length.
    """
    count = int(points)
    if count < 2:
        raise ValueError("an SBP grid needs at least two radii")
    inner, outer = float(r_inner), float(r_outer)
    if not np.isfinite(inner) or not np.isfinite(outer) or not outer > inner:
        raise ValueError("radial SBP interval must be finite and ordered")
    spacing = (outer - inner) / (count - 1)
    weights = np.full(count, spacing, dtype=float)
    weights[0] = weights[-1] = 0.5 * spacing
    return weights


def weighted_gram(phi, weights):
    """⟨φ_i, φ_j⟩ = Σ_{a,p} w_p φ*_{a p i} φ_{a p j}."""
    field = np.asarray(phi, dtype=complex)
    quadrature = np.asarray(weights, dtype=float)
    if field.ndim != 3 or field.shape[0] != 2:
        raise ValueError("canonical phi shape is (2, points, rank)")
    if quadrature.shape != (field.shape[1],) or not np.isfinite(quadrature).all():
        raise ValueError("SBP weights must be finite and match the radial points")
    if np.any(quadrature <= 0.0):
        raise ValueError("SBP weights must be positive")
    return np.einsum("apr,p,aps->rs", np.conjugate(field), quadrature, field, optimize=True)


def manufactured_quadrature_control(*, center, width, r_inner, r_outer, points, occupations):
    """Quadrature smoke shape only. It is not a prepared positive-energy source.

    No Killing phase or radial mode is solved here, so a positive frequency is
    not established. The scattering-packet owner supplies that preparation.
    """
    weights_c = np.asarray(occupations, dtype=float)
    if weights_c.ndim != 1 or weights_c.size < 1:
        raise ValueError("Gaussian occupations need one weight per column")
    if not np.isfinite(weights_c).all() or np.any(weights_c < 0.0) or np.any(weights_c > 1.0):
        raise ValueError("Gaussian occupations must lie in [0, 1]")
    radius = np.linspace(float(r_inner), float(r_outer), int(points), dtype=float)
    quadrature = radial_sbp_weights(r_inner, r_outer, int(points))
    centre, scale = _finite(center, "center"), _positive(width, "width")
    envelope = np.exp(-0.5 * ((radius - centre) / scale) ** 2)
    phi = np.zeros((2, radius.size, weights_c.size), dtype=complex)
    for column in range(weights_c.size):
        row = column % 2
        phi[row, :, column] = envelope
        norm = np.sqrt(np.real(np.sum(quadrature * np.abs(phi[:, :, column]) ** 2)))
        if norm <= 0.0:
            raise ValueError("Gaussian shell vanished on this interval")
        phi[:, :, column] /= norm
    return {
        "phi": phi,
        "radius": radius,
        "weights": quadrature,
        "occupations": weights_c.copy(),
        "gram": weighted_gram(phi, quadrature),
        "sigma2": SIGMA2.copy(),
        "sigma_angular": SIGMA_ANGULAR.copy(),
        "angular_coupling": angular_coupling(radius),
        "kappa": KAPPA,
        "multiplicity": MULTIPLICITY,
        "multiplicity_counted_once": True,
        "electric_current": ELECTRIC_CURRENT,
        "angular_is_4d_mass": False,
        "role": "manufactured_quadrature_control",
        "positivity_unestablished": True,
        "prepared_positive_energy_source": False,
        "exterior_killing_frequency": None,
    }


@dataclass(frozen=True)
class RNReference:
    """Source-free magnetic RN data fixed by the action and one mass parameter."""

    A: float
    C_F: float
    flux: float
    mass: float
    magnetic_r2: float
    G_N: float
    V4: float = AUDITED_V4

    @classmethod
    def from_action(cls, A, C_F, flux, mass):
        """Fixed-charge branch. V4 stays 0; P² is magnetic_r2."""
        newton = newton_constant(A)
        charge = magnetic_radius_square(A, C_F, flux)
        return cls(
            A=_positive(A, "A"),
            C_F=_positive(C_F, "C_F"),
            flux=_finite(flux, "flux"),
            mass=_positive(mass, "mass"),
            magnetic_r2=charge,
            G_N=newton,
            V4=AUDITED_V4,
        )

    @property
    def discriminant(self):
        return self.mass ** 2 - self.magnetic_r2

    @property
    def extremal(self):
        scale = max(self.mass ** 2, 1.0)
        return abs(self.discriminant) <= 1e-12 * scale

    @property
    def horizons(self):
        """Roots of r² - 2 mass r + magnetic_r2 = 0. Empty when the discriminant is negative."""
        if self.discriminant < 0.0 and not self.extremal:
            return {"r_minus": None, "r_plus": None, "extremal": False, "black_hole": False}
        root = 0.0 if self.extremal else float(np.sqrt(self.discriminant))
        return {
            "r_minus": self.mass - root,
            "r_plus": self.mass + root,
            "extremal": bool(self.extremal or root == 0.0),
            "black_hole": True,
        }

    @property
    def shift_domain_lower(self):
        """Smallest radius with a real ingoing shift, magnetic_r2/(2 mass)."""
        return self.magnetic_r2 / (2.0 * self.mass)

    def f(self, radius):
        """Metric factor 1 - 2 mass/r + magnetic_r2/r²."""
        r = np.asarray(radius, dtype=float)
        values = 1.0 - 2.0 * self.mass / r + self.magnetic_r2 / r ** 2
        return float(values) if values.ndim == 0 else values

    def beta_argument(self, radius):
        r = np.asarray(radius, dtype=float)
        return 2.0 * self.mass / r - self.magnetic_r2 / r ** 2

    def beta_rn(self, radius):
        """Positive ingoing shift. Real only for r >= magnetic_r2/(2 mass)."""
        argument = self.beta_argument(radius)
        if np.any(argument < -1e-12):
            raise ValueError("real ingoing shift requires r >= magnetic_r2/(2*mass)")
        values = np.sqrt(np.maximum(np.asarray(argument, dtype=float), 0.0))
        return float(values) if values.ndim == 0 else values

    def beta_derivative(self, radius):
        beta = np.asarray(self.beta_rn(radius), dtype=float)
        r = np.asarray(radius, dtype=float)
        slope = (-self.mass / r ** 2 + self.magnetic_r2 / r ** 3) / beta
        return float(slope) if slope.ndim == 0 else slope

    def beta_second_derivative(self, radius):
        beta = np.asarray(self.beta_rn(radius), dtype=float)
        r = np.asarray(radius, dtype=float)
        numerator = -self.mass / r ** 2 + self.magnetic_r2 / r ** 3
        numerator_r = 2.0 * self.mass / r ** 3 - 3.0 * self.magnetic_r2 / r ** 4
        second = numerator_r / beta - numerator ** 2 / beta ** 3
        return float(second) if second.ndim == 0 else second

    def vacuum_lapse(self, radius):
        """Derived vacuum normalization. The value is 1; the radius only checks the chart."""
        r = np.asarray(radius, dtype=float)
        if r.size == 0 or not np.isfinite(r).all() or np.min(r) <= 0.0:
            raise ValueError("lapse radius must be finite and positive")
        ones = np.ones(r.shape, dtype=float)
        return float(ones) if ones.ndim == 0 else ones

    def characteristic_speeds(self, radius):
        """Radial null speeds dr/dt of ds² = dt² - (dr + β dt)². Outgoing, then ingoing."""
        beta = np.asarray(self.beta_rn(radius), dtype=float)
        outgoing, ingoing = 1.0 - beta, -(1.0 + beta)
        if outgoing.ndim == 0:
            return float(outgoing), float(ingoing)
        return outgoing, ingoing

    def pure_outflow(self, radius):
        """Both radial null speeds leave a domain whose inner edge is this radius."""
        outgoing, ingoing = self.characteristic_speeds(radius)
        return np.asarray(outgoing) < 0.0 and np.asarray(ingoing) < 0.0 if np.ndim(outgoing) == 0 else (
            (outgoing < 0.0) & (ingoing < 0.0)
        )

    def pg_metric_jets(self, radius):
        """Covariant (t, r) jets of the vacuum chart. N = 1 is derived, and ∂_r N = 0."""
        r = np.asarray(radius, dtype=float)
        beta = np.asarray(self.beta_rn(r), dtype=float)
        beta_r = np.asarray(self.beta_derivative(r), dtype=float)
        factor = np.asarray(self.f(r), dtype=float)
        factor_r = 2.0 * self.mass / r ** 2 - 2.0 * self.magnetic_r2 / r ** 3
        return {
            "r": r,
            "N": np.ones(r.shape, dtype=float),
            "N_r": np.zeros(r.shape, dtype=float),
            "beta": beta,
            "beta_r": beta_r,
            "beta_rr": np.asarray(self.beta_second_derivative(r), dtype=float),
            "f": factor,
            "f_r": factor_r,
            "g_tt": factor,
            "g_tr": -beta,
            "g_rr": -np.ones(r.shape, dtype=float),
            "g_tt_r": factor_r,
            "g_tr_r": -beta_r,
            "g_rr_r": np.zeros(r.shape, dtype=float),
            "lapse_derived": True,
            "coupled_lapse_kept_at_one": False,
        }

    def invariants(self, radius):
        """Scalar curvature, Ricci square, Kretschmann and Weyl square."""
        r = np.asarray(radius, dtype=float)
        if r.size == 0 or not np.isfinite(r).all() or np.min(r) <= 0.0:
            raise ValueError("curvature radius must be finite and positive")
        charge = self.magnetic_r2
        ricci2 = 4.0 * charge ** 2 / r ** 8
        kretschmann = (
            48.0 * self.mass ** 2 / r ** 6
            - 96.0 * self.mass * charge / r ** 7
            + 56.0 * charge ** 2 / r ** 8
        )
        values = {
            "R4": np.zeros(r.shape, dtype=float),
            "Ricci2": ricci2,
            "K": kretschmann,
            "Weyl2": kretschmann - 2.0 * ricci2,
        }
        if r.ndim == 0:
            return {name: float(value) for name, value in values.items()}
        return values

    def tetrad_riemann(self, radius):
        """Orthonormal PG components. R_0101 is the radial tide in the ingoing frame."""
        self.beta_rn(radius)
        r = np.asarray(radius, dtype=float)
        charge = self.magnetic_r2
        radial = (2.0 * self.mass * r - 3.0 * charge) / r ** 4
        angular = -(self.mass * r - charge) / r ** 4
        values = {
            "R_0101": radial,
            "R_0202": angular,
            "R_0303": angular,
            "R_1212": -angular,
            "R_0121": np.zeros(r.shape, dtype=float),
            "R_0212": np.zeros(r.shape, dtype=float),
        }
        if r.ndim == 0:
            return {name: float(value) for name, value in values.items()}
        return values

    def misner_sharp_pg(self, radius):
        """Bare Misner–Sharp mass (r/2) β² on the vacuum chart."""
        beta = np.asarray(self.beta_rn(radius), dtype=float)
        r = np.asarray(radius, dtype=float)
        values = 0.5 * r * beta ** 2
        return float(values) if values.ndim == 0 else values

    def charged_mass_pg(self, radius):
        """Bare Misner–Sharp mass plus magnetic_r2/(2r). Constant on exact RN."""
        r = np.asarray(radius, dtype=float)
        values = np.asarray(self.misner_sharp_pg(r), dtype=float) + self.magnetic_r2 / (2.0 * r)
        return float(values) if values.ndim == 0 else values


def misner_sharp(radius, grad_r_squared):
    """(+---) Misner–Sharp mass (r/2) (1 + ∇^a r ∇_a r)."""
    r = np.asarray(radius, dtype=float)
    grad2 = np.asarray(grad_r_squared, dtype=float)
    return 0.5 * r * (1.0 + grad2)


def charged_mass_from_gradient(radius, grad_r_squared, magnetic_r2):
    """Charged mass. The bare Misner–Sharp value is smaller by magnetic_r2/(2r)."""
    r = np.asarray(radius, dtype=float)
    charge = float(magnetic_r2)
    return misner_sharp(r, grad_r_squared) + charge / (2.0 * r)


def conformal_grad_r_squared(Q, radius, radius_x, p_Q, A):
    """∇r·∇r on the historical conformal chart N = r Q, β = 0.

    r_t = Q p_Q / (2 a r) with a = -8 π A, so
    (r_t² - r_x²) / N² = p_Q² / (4 a² r⁴) - r_x² / (r² Q²).
    The kinetic anchor is not an input.
    """
    lapse_density = np.asarray(Q, dtype=float)
    r = np.asarray(radius, dtype=float)
    slope = np.asarray(radius_x, dtype=float)
    momentum = np.asarray(p_Q, dtype=float)
    coefficient = -8.0 * np.pi * _positive(A, "A")
    kinetic = momentum ** 2 / (4.0 * coefficient ** 2 * r ** 4)
    return kinetic - slope ** 2 / (r ** 2 * lapse_density ** 2)


def _stats(values):
    array = np.asarray(values, dtype=float)
    return {
        "mean": float(np.mean(array)),
        "min": float(np.min(array)),
        "max": float(np.max(array)),
        "max_abs": float(np.max(np.abs(array))),
    }



def sourcefree_static_rn(radius, mass, *, A=AUDITED_A, C_F=AUDITED_C_F, flux=AUDITED_FLUX):
    """Static source-free RN sample on the ingoing PG chart.

    Adapter contract, with no field state and no campaign:

    sourcefree_static_rn(radius, mass, *, A=AUDITED_A, C_F=AUDITED_C_F, flux=AUDITED_FLUX)
        -> dict with V4=0, G_N, P2=magnetic_r2, horizons, metric_jets,
           invariants, tetrad, misner_sharp, charged_mass.
    """
    reference = RNReference.from_action(A, C_F, flux, mass)
    if reference.V4 != 0.0:
        raise ValueError("the fixed-charge branch keeps V4 at 0")
    jets = reference.pg_metric_jets(radius)
    return {
        "interface": "nsc_rn_reference.sourcefree_static_rn",
        "branch": "fixed-charge",
        "V4": reference.V4,
        "G_N": reference.G_N,
        "P2": reference.magnetic_r2,
        "magnetic_r2": reference.magnetic_r2,
        "A": reference.A,
        "C_F": reference.C_F,
        "flux": reference.flux,
        "mass": reference.mass,
        "horizons": reference.horizons,
        "shift_domain_lower": reference.shift_domain_lower,
        "metric_jets": jets,
        "invariants": reference.invariants(radius),
        "tetrad": reference.tetrad_riemann(radius),
        "misner_sharp": reference.misner_sharp_pg(radius),
        "charged_mass": reference.charged_mass_pg(radius),
        "source_free": True,
        "static": True,
        "field_state": None,
        "campaign": None,
        "coupled_calibration": False,
        "br_radius_rescaled": False,
        "coupled_lapse_kept_at_one": False,
    }


def placement(reference, *, parent, child, excision, outer="absorbing"):
    """Parent exterior and a different near-horizon child interval."""
    horizons = reference.horizons
    reasons = []
    if not horizons["black_hole"]:
        return {
            "supported": False,
            "blocker": "no real horizon, so pure-outflow excision is undefined",
            "horizons": horizons,
        }
    r_minus, r_plus = horizons["r_minus"], horizons["r_plus"]

    def _interval(pair, name):
        if len(pair) != 2:
            raise ValueError(name + " interval needs two radii")
        inner, outer_radius = float(pair[0]), float(pair[1])
        if not np.isfinite(inner) or not np.isfinite(outer_radius) or not outer_radius > inner:
            raise ValueError(name + " interval must be finite and ordered")
        return inner, outer_radius

    parent_inner, parent_outer = _interval(parent, "parent")
    child_inner, child_outer = _interval(child, "child")
    excision_radius = float(excision)
    disjoint = parent_inner >= child_outer or child_inner >= parent_outer
    near = child_inner > r_minus and child_outer > r_plus and child_inner < r_plus + max(r_plus - r_minus, r_plus)
    if not parent_inner > r_plus:
        reasons.append("parent inner radius is not outside the outer horizon")
    if not disjoint or (parent_inner, parent_outer) == (child_inner, child_outer):
        reasons.append("parent and child areal intervals are not different")
    if not near:
        reasons.append("child interval is not a near-outer-horizon collar")
    if horizons["extremal"]:
        reasons.append("extremal horizons leave no open pure-outflow excision interval")
    elif not (r_minus < excision_radius < r_plus):
        reasons.append("excision is not strictly between the horizons")
    elif not reference.pure_outflow(excision_radius):
        reasons.append("excision radius is not pure outflow")
    if outer != "absorbing":
        reasons.append("outer boundary must be nonperiodic and absorbing")
    return {
        "supported": not reasons,
        "reasons": reasons,
        "parent": [parent_inner, parent_outer],
        "child": [child_inner, child_outer],
        "excision": excision_radius,
        "outer_boundary": outer,
        "periodic": False,
        "horizons": horizons,
        "coupled_lapse_kept_at_one": False,
    }


def _sha256_file(path):
    hasher = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            hasher.update(block)
    return hasher.hexdigest()


def _load_json(path):
    return json.loads(Path(path).read_text())


def decode_canonical_pi(record, arrays):
    """Nodal momenta p = W @ (pi / dx_g) for an explicit canonical_pi record.

    The kinetic anchor is not read. A missing representation is rejected.
    dx_g comes from the arrays, the record, or nf matched to stored W.
    """
    if not isinstance(record, dict):
        raise ValueError("checkpoint record must be a mapping")
    if "momentum_representation" not in record or record.get("momentum_representation") in (None, ""):
        raise ValueError("momentum representation is missing; refusing to guess")
    representation = record["momentum_representation"]
    if representation != CANONICAL_PI:
        raise ValueError("refusing momentum representation " + str(representation))
    for name in ("pi_Q", "pi_r", "W"):
        if name not in arrays:
            raise ValueError("canonical_pi checkpoint is missing " + name)
    frame = np.asarray(arrays["W"], dtype=float)
    if frame.ndim != 2 or frame.shape[0] != frame.shape[1]:
        raise ValueError("stored W must be square")
    count = int(frame.shape[0])
    pi_q = np.asarray(arrays["pi_Q"], dtype=float)
    pi_r = np.asarray(arrays["pi_r"], dtype=float)
    if pi_q.shape != (count,) or pi_r.shape != (count,):
        raise ValueError("canonical pi does not match W")
    spacing = _canonical_dx(record, arrays, count)
    return {
        "representation": CANONICAL_PI,
        "dx_g": spacing,
        "n": count,
        "p_Q": frame @ (pi_q / spacing),
        "p_r": frame @ (pi_r / spacing),
        "kinetic_anchor_ignored": True,
    }


def _canonical_dx(record, arrays, count):
    if "dx_g" in arrays:
        return float(np.asarray(arrays["dx_g"]).reshape(-1)[0])
    if record.get("dx_g") is not None:
        return float(record["dx_g"])
    fermions = record.get("nf")
    if fermions is None or int(fermions) - 1 != int(count):
        raise ValueError("canonical_pi decode needs dx_g or nf matching W")
    return HISTORICAL_PERIOD / float(count)


def _flat_radius_gradient(radius):
    values = np.asarray(radius, dtype=float)
    span = float(np.max(values) - np.min(values))
    scale = max(float(np.max(np.abs(values))), 1.0)
    if span > 1e-12 * scale:
        return None
    return np.zeros(values.shape, dtype=float)


def _mass_from_decoded(lapse, radius, momentum, charge):
    slope = _flat_radius_gradient(radius)
    if slope is None:
        return {
            "status": "unresolved",
            "reason": "areal radius is not flat and this slice does not estimate its gradient",
            "charged_mass": None,
            "misner_sharp": None,
            "kinetic_anchor_from_decoded_p": None,
        }
    grad2 = conformal_grad_r_squared(lapse, radius, slope, momentum, AUDITED_A)
    bare = misner_sharp(radius, grad2)
    charged = charged_mass_from_gradient(radius, grad2, charge)
    return {
        "status": "decoded_gradient",
        "radius_gradient": "flat sealed radius; spatial derivative not estimated",
        "charged_mass": _stats(charged),
        "misner_sharp": _stats(bare),
        "bare_difference": _stats(charged - bare),
        "charge_over_2r": _stats(charge / (2.0 * np.asarray(radius, dtype=float))),
        "kinetic_anchor_from_decoded_p": float(np.mean(np.asarray(momentum, dtype=float) ** 2 / np.asarray(radius, dtype=float) ** 3)),
        "fitted_to_horizons": False,
        "fitted_to_kinetic_anchor": False,
    }


CURVATURE_UNRESOLVED = {
    "status": "unresolved",
    "physical_ground_truth": False,
    "coupled_calibration": False,
    "method_used": None,
    "reason": (
        "leading.metric_jets was not called. An owned Fourier estimate of this "
        "slice previously reached |R4| about 0.757, while prior authenticated "
        "metric jets on the same family are about 0.167 to 0.160. Neither "
        "number is adopted."
    ),
}


def _source_free_record(directory, case_name, ordinal=0):
    from .nsc_discovery_episode import load_checkpoint
    folder = Path(directory)
    record, arrays = load_checkpoint(folder, case_name, ordinal)
    decoded = decode_canonical_pi(record, arrays)
    charge = magnetic_radius_square(AUDITED_A, AUDITED_C_F, AUDITED_FLUX)
    mass = _mass_from_decoded(arrays["Q"], arrays["r"], decoded["p_Q"], charge)
    json_path = folder / f"{case_name}-{int(ordinal):06d}.json"
    npz_path = folder / f"{case_name}-{int(ordinal):06d}.npz"
    phi_max = max(float(np.max(np.abs(arrays["phi0"]))), float(np.max(np.abs(arrays["phi1"]))))
    return {
        "json": str(json_path.relative_to(LAB)),
        "npz": str(npz_path.relative_to(LAB)),
        "authenticated": True,
        "ordinal": int(ordinal),
        "case_id": case_name,
        "momentum_representation": record["momentum_representation"],
        "dx_g": decoded["dx_g"],
        "n": decoded["n"],
        "nf": record.get("nf"),
        "declared_k": None if record.get("common_k") is None else float(record["common_k"]),
        "phi_max_abs": phi_max,
        "areal_radius": _stats(arrays["r"]),
        "stored_constraints": record.get("control_preparation", {}).get("leading_constraints"),
        "mass": mass,
        "curvature": dict(CURVATURE_UNRESOLVED),
        "input_hashes": {
            str(json_path.relative_to(LAB)): _sha256_file(json_path),
            str(npz_path.relative_to(LAB)): _sha256_file(npz_path),
        },
    }


def _walk(node, key, found):
    if isinstance(node, dict):
        if key in node and isinstance(node[key], dict):
            found.append(node[key])
        for value in node.values():
            _walk(value, key, found)
    elif isinstance(node, list):
        for value in node:
            _walk(value, key, found)


def _vacuum_control(lab):
    folder = lab / "results/development/nsc-discovery-vacuum-control-v1"
    record_path = folder / "record.json"
    record = _load_json(record_path)
    payload_path = folder / record["payload"]
    payload_hash = _sha256_file(payload_path)
    coefficients = record["locked_coefficients"]
    coefficient_match = (
        coefficients["A"] == AUDITED_A
        and coefficients["C_F"] == AUDITED_C_F
        and coefficients["flux"] == AUDITED_FLUX
        and float(coefficients.get("V_rel", 0.0)) == 0.0
    )
    stored = []
    _walk(record, "metric_invariants", stored)
    def _errors(name):
        return [float(item[name]["max_error_from_exact"]) for item in stored if name in item]
    return {
        "authenticated": payload_hash == record["payload_sha256"] and coefficient_match,
        "payload_sha256_matches": payload_hash == record["payload_sha256"],
        "coefficients_match_audited_action": coefficient_match,
        "used_as_rn_test": False,
        "br_radius_rescaled": False,
        "stored_station_count": len(stored),
        "stored_R4_error_max": max(_errors("R4")) if stored else None,
        "stored_Ricci2_error_max": max(_errors("Ricci2")) if stored else None,
        "stored_K_error_max": max(_errors("K")) if stored else None,
        "note": "Stored conformal Bertotti-Robinson invariants are not an RN calibration.",
        "rerun": False,
        "input_hashes": {
            str(record_path.relative_to(LAB)): _sha256_file(record_path),
            str(payload_path.relative_to(LAB)): payload_hash,
        },
    }


def _cut_confirmation(lab):
    folder = lab / "results/development/nsc-discovery-parent-cut-confirmation-v1"
    record_path = folder / "confirmation.json"
    payload_path = folder / "confirmation.npz"
    record = _load_json(record_path)
    payload_hash = _sha256_file(payload_path)
    commit = str(record.get("producing_commit", ""))
    arms = record.get("arms", [])
    matched = [bool(arm.get("matched")) for arm in arms]
    constraints = []
    for arm in arms:
        constraints.append({
            "preparation": arm.get("preparation", {}).get("actual_constraints"),
            "endpoint": arm.get("endpoint", {}).get("constraints"),
        })
    representation = record.get("momentum_representation")
    return {
        "authenticated": (
            payload_hash == record.get("payload_sha256")
            and commit == CUT_CONFIRMATION_COMMIT
            and bool(arms)
            and all(matched)
        ),
        "producing_commit": commit,
        "commit_is_5cb775a7": commit.startswith("5cb775a7"),
        "payload_sha256_matches": payload_hash == record.get("payload_sha256"),
        "arms_matched": matched,
        "complete": bool(arms) and all(matched) and commit == CUT_CONFIRMATION_COMMIT,
        "constraints": constraints,
        "momentum_decoded": False,
        "decode_rejected": representation != CANONICAL_PI,
        "curvature": dict(CURVATURE_UNRESOLVED),
        "unresolved_jets": [
            "cut-confirmation arrays are not a canonical_pi checkpoint, so momenta were not decoded",
            "occupied accelerations and matter-coupled curvature remain unresolved",
        ],
        "rerun": False,
        "input_hashes": {
            str(record_path.relative_to(LAB)): _sha256_file(record_path),
            str(payload_path.relative_to(LAB)): payload_hash,
        },
    }


def assess_historical(lab=None):
    """Read sealed checkpoints and report decoded masses. Curvature stays unresolved.

    No trajectory is started and no file is written.
    """
    root = LAB if lab is None else Path(lab)
    charge = magnetic_radius_square(AUDITED_A, AUDITED_C_F, AUDITED_FLUX)
    strong = root / "results/development/nsc-discovery-parent-strong-empty-k-v1"
    weak = root / "results/development/nsc-discovery-parent-empty-source-v1"
    strong_plus = _source_free_record(strong, "nf128_plus_balanced_source_free_dt0.001", 0)
    strong_minus = _source_free_record(strong, "nf128_minus_balanced_source_free_dt0.001", 0)
    weak_plus = _source_free_record(weak, "nf128_plus_balanced_source_free_dt0.001", 0)
    vacuum = _vacuum_control(root)
    confirmation = _cut_confirmation(root)
    decoded_masses = []
    for label, case in (
        ("strong-empty-k-plus", strong_plus),
        ("strong-empty-k-minus", strong_minus),
        ("old-weak-control", weak_plus),
    ):
        mass = case["mass"]
        decoded_masses.append({
            "case": label,
            "declared_k": case["declared_k"],
            "kinetic_anchor_from_decoded_p": mass.get("kinetic_anchor_from_decoded_p"),
            "charged_mass": None if mass["charged_mass"] is None else mass["charged_mass"]["mean"],
            "misner_sharp": None if mass["misner_sharp"] is None else mass["misner_sharp"]["mean"],
            "areal_radius": case["areal_radius"]["mean"],
            "magnetic_r2": charge,
            "dx_g": case["dx_g"],
            "n": case["n"],
            "fitted_to_kinetic_anchor": False,
            "fitted_to_horizons": False,
        })
    hashes = {}
    for case in (strong_plus, strong_minus, weak_plus, vacuum, confirmation):
        hashes.update(case["input_hashes"])
    unresolved = [
        "k=0.02 has no sealed canonical_pi checkpoint",
        CURVATURE_UNRESOLVED["reason"],
        "coupled lapse on the new radial chart has not been integrated",
        *confirmation["unresolved_jets"],
    ]
    authenticated = all((
        strong_plus["authenticated"],
        strong_minus["authenticated"],
        weak_plus["authenticated"],
        vacuum["authenticated"],
        confirmation["authenticated"],
    ))
    return {
        "schema": ASSESSMENT_SCHEMA,
        "authenticated": authenticated,
        "blocker": None if authenticated else "sealed historical authentication failed",
        "records_mutated": False,
        "production_trajectory": False,
        "rerun": False,
        "final_record_written": False,
        "audited_coefficients": {
            "A": AUDITED_A,
            "C_F": AUDITED_C_F,
            "flux": AUDITED_FLUX,
            "V4": AUDITED_V4,
            "magnetic_r2": charge,
            "P2": charge,
            "G_N": newton_constant(AUDITED_A),
        },
        "decoded_masses": decoded_masses,
        "magnetic_charge_rescaled": False,
        "br_radius_rescaled_with_k": False,
        "coupled_calibration": False,
        "campaigns": {
            "strong-empty-k-v1": {"plus": strong_plus, "minus": strong_minus},
            "old-weak-control": weak_plus,
            "cut-confirmation-v1": confirmation,
            "vacuum-br-control": vacuum,
        },
        "supported_domains": [
            {
                "domain": "exact static source-free RN",
                "statement": "sourcefree_static_rn at fixed G_N and P2 with V4=0",
            },
            {
                "domain": "decoded conformal mass",
                "statement": "charged mass from canonical_pi decode and the conformal gradient, not from k",
            },
        ],
        "unresolved_jets": unresolved,
        "input_hashes": hashes,
        "k_is_mass": False,
    }


def write_assessment(assessment, path, *, producer_commit):
    """Create only results/development/nsc-rn-reference-assessment-v1.json.

    Requires a frozen commit and the assessment input hashes. Refuses sealed
    input directories, any other path, an existing file, and payloads over 64 MiB.
    This function is not called by the reference tests.
    """
    commit = str(producer_commit or "").lower()
    if len(commit) != 40 or any(character not in "0123456789abcdef" for character in commit):
        raise RuntimeError("final evidence writer requires a frozen 40-character root commit")
    target = Path(path).resolve()
    allowed = ASSESSMENT_OUTPUT.resolve()
    if target != allowed:
        raise RuntimeError(
            "writer accepts only results/development/nsc-rn-reference-assessment-v1.json"
        )
    for sealed in SEALED_INPUT_DIRECTORIES:
        root = sealed.resolve()
        if target == root or root in target.parents:
            raise RuntimeError("refusing to write inside a sealed input directory")
    hashes = assessment.get("input_hashes") if isinstance(assessment, dict) else None
    if not isinstance(hashes, dict) or not hashes:
        raise RuntimeError("assessment input_hashes are required before writing")
    if target.exists():
        raise FileExistsError(target)
    payload = dict(assessment)
    payload["producer_commit"] = commit
    payload["written_by"] = "nsc_rn_reference.write_assessment"
    payload["input_hashes"] = dict(hashes)
    body = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()
    if len(body) > int(MAX_ASSESSMENT_BYTES):
        raise RuntimeError("assessment exceeds the 64 MiB creation limit")
    temporary = target.with_name(".nsc-rn-reference-assessment-v1.json.tmp")
    flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY
    descriptor = os.open(temporary, flags, 0o644)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(body)
        os.chmod(temporary, 0o444)
        os.link(temporary, target)
    except FileExistsError:
        os.unlink(temporary)
        raise FileExistsError(target) from None
    finally:
        if temporary.exists():
            os.chmod(temporary, 0o644)
            temporary.unlink()
    os.chmod(target, 0o444)
    return target
