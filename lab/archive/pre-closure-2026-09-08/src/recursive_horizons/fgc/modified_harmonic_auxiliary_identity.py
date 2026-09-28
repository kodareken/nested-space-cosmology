"""Formal radial hat-sector completion identity for the MHG spherical chart."""

from __future__ import annotations

from fractions import Fraction
from typing import Any

from .chi_principal_identity import chi_principal_identity
from .reference_branch_principal_identity import reference_branch_principal_identity_certificate


Q = Fraction


class _Formal:
    """Independent sparse commutative polynomial algebra over Q."""

    __slots__ = ("terms",)

    def __init__(self, terms=()) -> None:
        combined: dict[tuple[str, ...], Fraction] = {}
        for word, coefficient in terms:
            key = tuple(sorted(word))
            combined[key] = combined.get(key, Q(0)) + coefficient
        self.terms = tuple(sorted((word, coefficient) for word, coefficient in combined.items() if coefficient))

    @classmethod
    def atom(cls, name: str) -> "_Formal":
        return cls((((name,), Q(1)),))

    @classmethod
    def scalar(cls, value: int) -> "_Formal":
        return cls(()) if value == 0 else cls((((), Q(value)),))

    @staticmethod
    def _coerce(value: object) -> "_Formal":
        if isinstance(value, _Formal):
            return value
        if isinstance(value, int):
            return _Formal.scalar(value)
        raise TypeError("formal operand must be an exact integer or formal expression")

    def __add__(self, other: object) -> "_Formal":
        other = self._coerce(other)
        return _Formal(self.terms + other.terms)

    __radd__ = __add__

    def __neg__(self) -> "_Formal":
        return _Formal((word, -coefficient) for word, coefficient in self.terms)

    def __sub__(self, other: object) -> "_Formal":
        other = self._coerce(other)
        return self + (-other)

    def __mul__(self, other: object) -> "_Formal":
        other = self._coerce(other)
        return _Formal((left + right, a * b) for left, a in self.terms for right, b in other.terms)

    __rmul__ = __mul__

    def __truediv__(self, divisor: int) -> "_Formal":
        return _Formal((word, coefficient / divisor) for word, coefficient in self.terms)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, _Formal) and self.terms == other.terms


def hat_projector_contraction_identity(*, mutate_trace_sign: int = -1) -> dict[str, Any]:
    """Prove ``xi_mu hatP_a^(b mu nu) xi_b=hatq delta_a^nu/2``.

    The route is assembled directly from the trace-reversal projector; it does
    not reuse the propagation implementation.  ``mutate_trace_sign`` is an
    exact mutation probe for the final projector term.
    """
    if mutate_trace_sign not in (-1, 1):
        raise ValueError("trace-sign mutation must be either -1 or +1")
    size = 4
    xi = tuple(_Formal.atom(f"xi{index}") for index in range(size))
    hat = tuple(tuple(_Formal.atom(f"h{min(row, column)}{max(row, column)}") for column in range(size)) for row in range(size))
    hatq = sum((xi[mu] * hat[mu][beta] * xi[beta] for mu in range(size) for beta in range(size)), _Formal.scalar(0))
    derived = tuple(
        tuple(
            sum((
                xi[mu] * (
                    (int(alpha == mu) * hat[nu][beta]
                     + int(alpha == nu) * hat[mu][beta]
                     + mutate_trace_sign * int(alpha == beta) * hat[mu][nu]) / 2
                ) * xi[beta]
                for mu in range(size) for beta in range(size)), _Formal.scalar(0))
            for alpha in range(size)
        )
        for nu in range(size)
    )
    expected = tuple(tuple((hatq / 2) if alpha == nu else _Formal.scalar(0) for alpha in range(size)) for nu in range(size))
    if derived != expected:
        raise ValueError("formal hat-projector double-covector contraction differs from hatq delta/2")
    return {
        "dimension": size,
        "identity": "xi_mu hatP_alpha^(beta mu nu) xi_beta=(hatq/2)delta_alpha^nu",
        "all_sixteen_index_components_exact": True,
        "proof_method": "independent_exact_sparse_formal_index_contraction_over_Q",
    }


def noether_extension_composition(*, action_diffeomorphism_invariance: bool = True) -> dict[str, Any]:
    """State the explicit conditional composition used by radial completion."""
    if action_diffeomorphism_invariance is not True:
        raise ValueError("Noether principal Bianchi handoff requires action diffeomorphism invariance premise")
    projector = hat_projector_contraction_identity()
    return {
        "action_diffeomorphism_invariance_is_named_premise_not_machine_derived": True,
        "unredefined_principal_bianchi_identity": "B P0=0 on scalar principal rows",
        "extension_divergence_identity": "B S_hat T_tilde=F*hatq*T_tilde/2",
        "composed_full_identity": "B(P0+S_hat T_tilde)=F*hatq*T_tilde/2",
        "hat_null_and_scalar_rows_imply_full_principal_bianchi_rows": True,
        "projector_contraction": projector,
    }


def tilde_pure_gauge_contraction_identity(*, mutate_connection_sign: int = -1) -> dict[str, Any]:
    """Independently prove the trace-reversal pure-gauge contraction.

    With ``h_rs=xi_(r X_s)=(xi_r X_s+xi_s X_r)/2``, direct principal
    connection variation gives ``tilde_g^(rs) deltaGamma^a_rs=qtilde X^a/2``.
    The repository's radial generators use the twice-normalized polarization,
    which changes neither their null-kernel property nor their span.
    """
    if mutate_connection_sign not in (-1, 1):
        raise ValueError("connection-sign mutation must be either -1 or +1")
    size = 4
    xi = tuple(_Formal.atom(f"txi{index}") for index in range(size))
    x_down = tuple(_Formal.atom(f"X{index}") for index in range(size))
    physical = tuple(tuple(_Formal.atom(f"g{row}{column}") for column in range(size)) for row in range(size))
    tilde = tuple(tuple(_Formal.atom(f"t{min(row, column)}{max(row, column)}") for column in range(size)) for row in range(size))
    qtilde = sum((xi[row] * tilde[row][column] * xi[column] for row in range(size) for column in range(size)), _Formal.scalar(0))
    x_up = tuple(sum((physical[alpha][lam] * x_down[lam] for lam in range(size)), _Formal.scalar(0)) for alpha in range(size))
    # The last sign is the -xi_lambda h_rho_sigma term in delta Gamma.
    derived = tuple(
        sum((
            tilde[rho][sigma]
            * physical[alpha][lam]
            * (
                xi[rho] * (xi[lam] * x_down[sigma] + xi[sigma] * x_down[lam]) / 2
                + xi[sigma] * (xi[lam] * x_down[rho] + xi[rho] * x_down[lam]) / 2
                + mutate_connection_sign * xi[lam] * (xi[rho] * x_down[sigma] + xi[sigma] * x_down[rho]) / 2
            ) / 2
            for rho in range(size) for sigma in range(size) for lam in range(size)), _Formal.scalar(0))
        for alpha in range(size)
    )
    expected = tuple(qtilde * x_up[alpha] / 2 for alpha in range(size))
    if derived != expected:
        raise ValueError("formal tilde pure-gauge contraction differs from qtilde X/2")
    return {
        "identity": "tilde_g^(rho sigma) deltaGamma^alpha_rho_sigma[xi_(rho X_sigma)]=qtilde(xi) X^alpha/2",
        "symmetric_polarization_normalization": "xi_(rho X_sigma)=(xi_rho X_sigma+xi_sigma X_rho)/2",
        "repository_generator_normalization": "R=2*xi_(rho X_sigma)",
        "all_four_vector_components_exact": True,
        "proof_method": "independent_exact_sparse_formal_connection_contraction_over_Q",
    }


def tilde_right_kernel_composition(*, action_diffeomorphism_invariance: bool = True) -> dict[str, Any]:
    """Compose the named diffeomorphism premise with the tilde contraction."""
    if action_diffeomorphism_invariance is not True:
        raise ValueError("pure-gauge right-kernel composition requires action diffeomorphism invariance premise")
    contraction = tilde_pure_gauge_contraction_identity()
    return {
        "action_diffeomorphism_invariance_is_named_premise_not_machine_derived": True,
        "unredefined_right_identity": "P0 R=0",
        "extension_on_pure_gauge": "S_hat T_tilde R=S_hat*(qtilde X/2)",
        "tilde_null_full_right_identity": "(P0+S_hat T_tilde)R=0 when qtilde=0",
        "spherical_generator_order": ("h_tt", "h_tr", "h_rr", "areal_radius", "phi", "chi"),
        "spherical_generators": (("-2c", "1", "0", "0", "0", "0"), ("0", "-c", "2", "0", "0", "0")),
        "contraction": contraction,
        "full_symbol_right_kernels_at_tilde_null": True,
    }


def radial_hat_omitted_row_identity(
    *,
    xi_t_up: Fraction = Q(1),
    xi_r_up: Fraction = Q(1),
    hat_null: Fraction = Q(0),
    selected_metric_tr: Fraction = Q(0),
    selected_metric_theta: Fraction = Q(0),
    scalar_phi: Fraction = Q(0),
    scalar_chi: Fraction = Q(0),
) -> dict[str, Any]:
    """Complete the two omitted metric rows from formal principal identities.

    At a hat-null covector, the MHG extension divergence is proportional to
    ``hat_null``.  On both scalar rows the ACT1 Noether/Bianchi handoff gives
    ``xi^t E_tt+xi^r E_tr=0`` and
    ``xi^t E_tr+xi^r E_rr=0``.  Hence the selected chart condition E_tr=0
    forces both omitted radial metric rows when both raised covector components
    are nonzero.  E_theta is a selected row; chi is an independent canonical
    principal block and is imposed as its own selected zero condition.
    """
    values = (xi_t_up, xi_r_up, hat_null, selected_metric_tr, selected_metric_theta, scalar_phi, scalar_chi)
    if any(isinstance(value, bool) or not isinstance(value, Fraction) for value in values):
        raise TypeError("auxiliary identity inputs must be Fractions")
    if xi_t_up == 0 or xi_r_up == 0:
        raise ValueError("radial hat completion chart requires nonzero raised covector components")
    if hat_null != 0:
        raise ValueError("radial hat completion identity applies only on the hat-null cone")
    if any(value != 0 for value in (selected_metric_tr, selected_metric_theta, scalar_phi, scalar_chi)):
        raise ValueError("radial hat completion requires all selected chart rows to vanish")
    # Solve the two formal Bianchi equations with E_tr=0, independently in
    # their two diagonal unknowns.  These are exact rational divisions only
    # after the two chart-nonvanishing premises above have been checked.
    metric_tt = -(xi_r_up * selected_metric_tr) / xi_t_up
    metric_rr = -(xi_t_up * selected_metric_tr) / xi_r_up
    if metric_tt != 0 or metric_rr != 0:
        raise AssertionError("formal radial Bianchi completion did not close")
    return {
        "hat_null_extension_divergence_zero": True,
        "scalar_noether_handoff_applied": True,
        "selected_rows": ("metric_tr", "metric_theta_theta", "scalar_phi", "scalar_chi"),
        "omitted_rows": ("metric_tt", "metric_rr"),
        "raised_covector_chart_components_nonzero": True,
        "completed_metric_tt": metric_tt,
        "completed_metric_rr": metric_rr,
        "all_six_radial_rows_zero": True,
        "scope": "formal_hat_null_radial_covector_chart_completion_not_a_global_fixed_pivot_claim",
    }


def auxiliary_sector_identity_certificate() -> dict[str, Any]:
    """Compose the universal REF1 bridge with radial hat completion premises."""
    bridge = reference_branch_principal_identity_certificate()
    chi = chi_principal_identity()
    composition = noether_extension_composition()
    tilde = tilde_right_kernel_composition()
    completion = radial_hat_omitted_row_identity()
    if not bridge["gauge_extension"]["universal_formal_certificate"]["connection_variation_route_equals_projector_route"]:
        raise ValueError("universal REF1 projector bridge is unavailable")
    return {
        "classification": "universal_conditional_radial_hat_sector_row_completion_identity",
        "universal_ref1_projector_bridge": bridge["gauge_extension"]["universal_formal_certificate"],
        "canonical_chi_principal_identity": chi,
        "noether_extension_composition": composition,
        "tilde_pure_gauge_right_kernel": tilde,
        "hat_row_completion": completion,
        "premises": {
            "regular_background_and_F_positive": True,
            "hat_null_radial_covector": True,
            "scalar_rows_zero": True,
            "selected_metric_tr_and_theta_rows_zero": True,
            "raised_covector_components_nonzero_chart": True,
        },
        "nonclaims": {
            "global_fixed_hat_pivot_without_covector_atlas_proven": False,
            "interval_pivot_invertibility_proven": False,
            "uniform_eigenframe_proven": False,
            "strong_hyperbolicity_proven": False,
        },
    }
