"""Branch-owned SGB-L adapter for the universal covariant principal symbol.

The tensor-level principal algebra is shared ACT1 mathematics. Its existing
spherical builder is FGC-QR-specific, so SGB-L supplies a separate background:
``F=Mpl^2``, ``F'=0``, ``f'=alpha_gb`` and
``nabla nabla f=alpha_gb*nabla nabla phi``. Point coefficient and directional
records retain tensors, hashes and sampled spectra. They do not classify
health, enclose a neighborhood, run a source, or authorize a trajectory.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from math import isfinite
from typing import Sequence

import numpy as np

from .evolution.covariant_principal_health import (
    CovariantPrincipalBackground,
    FIELD_ORDER,
    PrincipalCoefficientTensors,
    _orthonormal_frame,
    _spherical_scalar_hessian,
    coefficient_tensor_symbol,
    complete_modified_harmonic_principal_symbol,
    frobenius_normalized_principal_symbol,
    point_spectrum_diagnostics,
    principal_coefficient_tensors,
    second_order_coefficient_matrices,
)
from .evolution.floating_jet import scalar_primal
from .evolution.multidirectional_health import (
    CoefficientDeformationBounds,
    coefficient_deformation_bounds,
)
from .modified_harmonic import modified_harmonic_symbol
from .sgb1_ctl1_source import HAT_NORMAL_FACTOR, TILDE_NORMAL_FACTOR
from .spherical_reduction import (
    BASE_FIELD_ORDER,
    SphericalState,
    direct_4d_curvature,
)


Q = Fraction
PINNED_UNIVERSAL_GEOMETRY_DEPENDENCIES = (
    "covariant_principal_health._orthonormal_frame",
    "covariant_principal_health._spherical_scalar_hessian",
)
FIELD_ORDER_12 = FIELD_ORDER
DECLARED_TILDE_NORMAL_FACTOR = Q(TILDE_NORMAL_FACTOR)
DECLARED_HAT_NORMAL_FACTOR = Q(HAT_NORMAL_FACTOR)
DECLARED_DIRECTION_COVECTORS = (
    (1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.6, 0.8, 0.0),
)
RADIAL_DIRECTION = DECLARED_DIRECTION_COVECTORS[0]
ANGULAR_DIRECTION = DECLARED_DIRECTION_COVECTORS[1]
MIXED_DIRECTION = DECLARED_DIRECTION_COVECTORS[2]
POINT_SPECTRUM_CLASSIFICATION = (
    "pointwise_binary64_all_polarization_spectrum_not_a_uniform_HYP2_certificate"
)
SYMBOL_RECONSTRUCTION_ATOL = 2.0e-12
_MIXED_FOUR_COVECTOR = (0.27, 0.3, -0.4, 0.75**0.5)


def _finite_array(name: str, value: object, shape: tuple[int, ...]) -> np.ndarray:
    array = np.asarray(value, dtype=np.float64)
    if array.shape != shape or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a finite array of shape {shape}")
    return array.copy()


def _companion(a_matrix, b_matrix, c_matrix) -> np.ndarray:
    size = a_matrix.shape[0]
    return np.block(
        [
            [np.zeros((size, size)), np.eye(size)],
            [
                -np.linalg.solve(a_matrix, c_matrix),
                -np.linalg.solve(a_matrix, b_matrix),
            ],
        ]
    )


def _spherical_embedding() -> np.ndarray:
    embedding = np.zeros((12, 6), dtype=np.float64)
    for column, row in enumerate((0, 1, 4)):
        embedding[row, column] = 1.0
    embedding[7, 3] = embedding[9, 3] = 1.0
    embedding[10, 4] = embedding[11, 5] = 1.0
    return embedding


def sgbl_covariant_principal_background(
    state: SphericalState,
) -> CovariantPrincipalBackground:
    """Map one complete SGB-L spherical two-jet to branch principal data.

    Reusing the universal curvature/frame arithmetic is not reuse of an
    FGC-QR health result. The returned point background is not an interval or
    strong-hyperbolicity certificate.
    """

    if not isinstance(state, SphericalState):
        raise TypeError("state must be a SphericalState")
    if state.branch != "SGB-L":
        raise ValueError("SGB-L principal background requires the linear branch")
    if state.beta != 0 or state.eta != 0 or state.alpha == 0:
        raise ValueError(
            "SGB-L principal background requires alpha*phi with beta=eta=0"
        )
    if tuple(BASE_FIELD_ORDER) != (
        "h_tt",
        "h_tr",
        "h_rr",
        "areal_radius",
        "phi",
        "chi",
    ):
        raise RuntimeError("spherical field order changed")

    riemann_coordinate, _inverse, connection = direct_4d_curvature(state)
    riemann_coordinate_array = np.vectorize(scalar_primal)(
        np.asarray(riemann_coordinate, dtype=object)
    ).astype(np.float64)
    hessian_coordinate = _spherical_scalar_hessian(state, connection)
    frame = _orthonormal_frame(state)
    riemann_orthonormal = np.einsum(
        "ma,nb,pc,qd,mnpq->abcd",
        frame,
        frame,
        frame,
        frame,
        riemann_coordinate_array,
    )
    hessian_phi = frame.T @ hessian_coordinate @ frame
    alpha_gb = float(state.alpha)
    planck_mass = float(state.planck_mass)
    if not np.isfinite(alpha_gb) or not np.isfinite(planck_mass) or planck_mass <= 0:
        raise ValueError("SGB-L action parameters must remain finite with Mpl>0")
    return CovariantPrincipalBackground(
        effective_planck_squared=planck_mass * planck_mass,
        effective_planck_prime=0.0,
        gb_coupling_prime=alpha_gb,
        riemann_lower=riemann_orthonormal,
        hessian_gb_lower=alpha_gb * hessian_phi,
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLDirectionSpectrum:
    """Binary64 spectrum at one orthonormal spatial direction.

    The classification is retained so a finite imaginary part cannot be
    relabelled as a uniform cone certificate.
    """

    spatial_covector: tuple[float, float, float]
    kinetic_determinant: float
    kinetic_smallest_singular_value: float
    maximum_abs_imaginary_part: float
    eigenvalues: tuple[complex, ...]
    spatial_symbol: np.ndarray
    classification: str = POINT_SPECTRUM_CLASSIFICATION

    def __post_init__(self) -> None:
        covector = tuple(float(entry) for entry in self.spatial_covector)
        if len(covector) != 3 or not all(isfinite(entry) for entry in covector):
            raise ValueError("spatial_covector must be three finite reals")
        object.__setattr__(self, "spatial_covector", covector)
        for name in (
            "kinetic_determinant",
            "kinetic_smallest_singular_value",
            "maximum_abs_imaginary_part",
        ):
            value = float(getattr(self, name))
            if not isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)
        eigenvalues = tuple(complex(value) for value in self.eigenvalues)
        if len(eigenvalues) != 24 or not all(
            isfinite(value.real) and isfinite(value.imag) for value in eigenvalues
        ):
            raise ValueError("direction spectrum must retain24 finite eigenvalues")
        object.__setattr__(self, "eigenvalues", eigenvalues)
        symbol = _finite_array("spatial_symbol", self.spatial_symbol, (12, 12))
        symbol.setflags(write=False)
        object.__setattr__(self, "spatial_symbol", symbol)
        if self.classification != POINT_SPECTRUM_CLASSIFICATION:
            raise ValueError("spectrum classification is not a health certificate")


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLPrincipalPointFacts:
    """Immutable point principal facts with health claims closed.

    Coefficient tensors, hashes, radial/angular symbols and spectra are
    retained. Cone, interval and pass fields are derived and false.
    """

    background: CovariantPrincipalBackground
    tilde_normal_factor: Fraction
    hat_normal_factor: Fraction
    coefficient_tensors: PrincipalCoefficientTensors
    coefficient_tensor_sha256: str
    deformation: CoefficientDeformationBounds
    radial_spectrum: SGBLDirectionSpectrum
    angular_spectrum: SGBLDirectionSpectrum
    mixed_spectrum: SGBLDirectionSpectrum
    radial_restriction_agrees_with_spherical_symbol: bool
    angular_spatial_symbol_differs_from_radial: bool
    reconstructed_symbol_agrees_at_mixed_covector: bool
    action_hessian_symmetry_defect: float
    gauge_identity_defect: float
    finite_coefficient_tensors: bool
    binary64_uncertainty_visible: bool = True

    def __post_init__(self) -> None:
        if type(self.background) is not CovariantPrincipalBackground:
            raise TypeError("background must be CovariantPrincipalBackground")
        if type(self.coefficient_tensors) is not PrincipalCoefficientTensors:
            raise TypeError("coefficient_tensors must be PrincipalCoefficientTensors")
        if type(self.deformation) is not CoefficientDeformationBounds:
            raise TypeError("deformation must be CoefficientDeformationBounds")
        object.__setattr__(
            self, "tilde_normal_factor", Fraction(self.tilde_normal_factor)
        )
        object.__setattr__(
            self, "hat_normal_factor", Fraction(self.hat_normal_factor)
        )
        if not 1 < self.tilde_normal_factor < self.hat_normal_factor:
            raise ValueError("point facts require 1 < tilde factor < hat factor")
        expected = _tensor_sha256(self.coefficient_tensors)
        if self.coefficient_tensor_sha256 != expected:
            raise ValueError("coefficient_tensor_sha256 does not match the tensors")
        if not self.binary64_uncertainty_visible:
            raise ValueError("binary64 directional facts must keep uncertainty visible")
        if not self.finite_coefficient_tensors:
            raise ValueError("coefficient tensors must be finite")
        if tuple(item.spatial_covector for item in self.spectra) != (
            RADIAL_DIRECTION, ANGULAR_DIRECTION, MIXED_DIRECTION
        ):
            raise ValueError("point spectra do not follow the declared direction order")
        for array in (
            self.background.riemann_lower,
            self.background.hessian_gb_lower,
            self.coefficient_tensors.time_time,
            self.coefficient_tensors.time_space,
            self.coefficient_tensors.space_space,
        ):
            array.setflags(write=False)

    @property
    def spectra(self) -> tuple[SGBLDirectionSpectrum, ...]:
        return (self.radial_spectrum, self.angular_spectrum, self.mixed_spectrum)

    @property
    def linear_branch_background_has_declared_derivatives(self) -> bool:
        background = self.background
        return (
            background.effective_planck_prime == 0.0
            and background.effective_planck_squared > 0.0
            and background.gb_coupling_prime != 0.0
        )

    @property
    def ten_metric_polarizations_plus_phi_chi_retained(self) -> bool:
        return FIELD_ORDER_12 == FIELD_ORDER and len(FIELD_ORDER_12) == 12

    @property
    def kr_sufficient_envelope_is_not_this_certificate(self) -> bool:
        return True

    @property
    def quantitative_all_covector_weak_coupling_health_envelope_passed(self) -> bool:
        return False

    @property
    def strongly_hyperbolic(self) -> bool:
        return False

    @property
    def interval_invertible(self) -> bool:
        return False

    @property
    def SGBL_branch_owned_and_healthy(self) -> bool:
        return False

    @property
    def cone_certificate(self) -> str:
        return "unqualified"


def _tensor_sha256(tensors: PrincipalCoefficientTensors) -> str:
    payload = (
        np.ascontiguousarray(tensors.time_time, dtype=np.float64).tobytes()
        + np.ascontiguousarray(tensors.time_space, dtype=np.float64).tobytes()
        + np.ascontiguousarray(tensors.space_space, dtype=np.float64).tobytes()
    )
    return sha256(payload).hexdigest()


def _direction_spectrum(
    background: CovariantPrincipalBackground,
    spatial_covector: Sequence[float],
    *,
    tilde: Fraction,
    hat: Fraction,
) -> SGBLDirectionSpectrum:
    tilde_f = float(tilde)
    hat_f = float(hat)
    diagnostics = point_spectrum_diagnostics(
        background,
        spatial_covector,
        tilde_normal_factor=tilde_f,
        hat_normal_factor=hat_f,
    )
    _a_matrix, _b_matrix, spatial = second_order_coefficient_matrices(
        background,
        spatial_covector,
        tilde_normal_factor=tilde_f,
        hat_normal_factor=hat_f,
    )
    kinetic = diagnostics["kinetic"]
    return SGBLDirectionSpectrum(
        spatial_covector=tuple(spatial_covector),  # type: ignore[arg-type]
        kinetic_determinant=float(kinetic["kinetic_determinant"]),
        kinetic_smallest_singular_value=float(
            kinetic["kinetic_smallest_singular_value"]
        ),
        maximum_abs_imaginary_part=float(diagnostics["maximum_abs_imaginary_part"]),
        eigenvalues=tuple(diagnostics["eigenvalues"]),
        spatial_symbol=spatial,
        classification=str(diagnostics["classification"]),
    )


def _radial_restriction_agrees(
    state: SphericalState,
    background: CovariantPrincipalBackground,
    *,
    tilde: Fraction,
    hat: Fraction,
) -> bool:
    a_matrix, b_matrix, c_matrix = second_order_coefficient_matrices(
        background,
        RADIAL_DIRECTION,
        tilde_normal_factor=float(tilde),
        hat_normal_factor=float(hat),
    )
    embedding = _spherical_embedding()
    covariant_roots = np.linalg.eigvals(
        _companion(
            embedding.T @ a_matrix @ embedding,
            embedding.T @ b_matrix @ embedding,
            embedding.T @ c_matrix @ embedding,
        )
    )
    exact_symbol = modified_harmonic_symbol(
        state,
        tilde_normal_factor=int(tilde),
        hat_normal_factor=int(hat),
    )
    exact_a = np.asarray(
        [
            [float(entry[2] if len(entry) > 2 else 0) for entry in row]
            for row in exact_symbol
        ]
    )
    exact_b = np.asarray(
        [
            [float(entry[1] if len(entry) > 1 else 0) for entry in row]
            for row in exact_symbol
        ]
    )
    exact_c = np.asarray(
        [
            [float(entry[0] if entry else 0) for entry in row]
            for row in exact_symbol
        ]
    )
    coordinate_speeds = np.linalg.eigvals(_companion(exact_a, exact_b, exact_c))
    orbit = np.asarray(
        [
            [float(state.h_tt.value), float(state.h_tr.value)],
            [float(state.h_tr.value), float(state.h_rr.value)],
        ]
    )
    inverse = np.linalg.inv(orbit)
    lapse = 1.0 / (-inverse[0, 0]) ** 0.5
    normal = -lapse * inverse[:, 0]
    radial_metric = float(state.h_rr.value) ** 0.5
    expected = radial_metric * (normal[1] - normal[0] * coordinate_speeds)
    return bool(
        np.allclose(
            np.sort_complex(covariant_roots),
            np.sort_complex(expected),
            rtol=0.0,
            atol=SYMBOL_RECONSTRUCTION_ATOL,
        )
    )


def sgbl_principal_point_facts(state: SphericalState) -> SGBLPrincipalPointFacts:
    """Return closed point facts for one complete SGB-L two-jet.

    Binary64 directional spectra remain diagnostics. They do not freeze a
    cone, an interval, or branch-owned health.
    """

    background = sgbl_covariant_principal_background(state)
    tilde = DECLARED_TILDE_NORMAL_FACTOR
    hat = DECLARED_HAT_NORMAL_FACTOR
    tensors = principal_coefficient_tensors(
        background,
        tilde_normal_factor=float(tilde),
        hat_normal_factor=float(hat),
    )
    reconstructed = coefficient_tensor_symbol(tensors, _MIXED_FOUR_COVECTOR)
    direct = frobenius_normalized_principal_symbol(
        complete_modified_harmonic_principal_symbol(
            background,
            _MIXED_FOUR_COVECTOR,
            tilde_normal_factor=float(tilde),
            hat_normal_factor=float(hat),
        )
    )
    reconstruction_agrees = bool(
        np.allclose(reconstructed, direct, rtol=0.0, atol=2.0e-14)
    )
    if not reconstruction_agrees:
        raise ValueError("extracted coefficient tensors do not reassemble the symbol")
    deformation = coefficient_deformation_bounds(background)
    radial = _direction_spectrum(
        background, RADIAL_DIRECTION, tilde=tilde, hat=hat
    )
    angular = _direction_spectrum(
        background, ANGULAR_DIRECTION, tilde=tilde, hat=hat
    )
    mixed = _direction_spectrum(
        background, MIXED_DIRECTION, tilde=tilde, hat=hat
    )
    angular_differs = float(
        np.linalg.norm(angular.spatial_symbol - radial.spatial_symbol)
    ) > 0.0
    finite = all(
        np.all(np.isfinite(array))
        for array in (
            tensors.time_time,
            tensors.time_space,
            tensors.space_space,
            radial.spatial_symbol,
            angular.spatial_symbol,
        )
    )
    return SGBLPrincipalPointFacts(
        background=background,
        tilde_normal_factor=tilde,
        hat_normal_factor=hat,
        coefficient_tensors=tensors,
        coefficient_tensor_sha256=_tensor_sha256(tensors),
        deformation=deformation,
        radial_spectrum=radial,
        angular_spectrum=angular,
        mixed_spectrum=mixed,
        radial_restriction_agrees_with_spherical_symbol=_radial_restriction_agrees(
            state, background, tilde=tilde, hat=hat
        ),
        angular_spatial_symbol_differs_from_radial=angular_differs,
        reconstructed_symbol_agrees_at_mixed_covector=reconstruction_agrees,
        action_hessian_symmetry_defect=float(
            deformation.action_hessian_symmetry_defect
        ),
        gauge_identity_defect=float(deformation.gauge_identity_defect),
        finite_coefficient_tensors=finite,
        binary64_uncertainty_visible=True,
    )


__all__ = [
    "ANGULAR_DIRECTION",
    "DECLARED_DIRECTION_COVECTORS",
    "DECLARED_HAT_NORMAL_FACTOR",
    "DECLARED_TILDE_NORMAL_FACTOR",
    "FIELD_ORDER_12",
    "MIXED_DIRECTION",
    "PINNED_UNIVERSAL_GEOMETRY_DEPENDENCIES",
    "POINT_SPECTRUM_CLASSIFICATION",
    "RADIAL_DIRECTION",
    "SGBLDirectionSpectrum",
    "SGBLPrincipalPointFacts",
    "sgbl_covariant_principal_background",
    "sgbl_principal_point_facts",
]
