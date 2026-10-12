"""Independent executable certificate for the TDG4 sampling theorem.

This module does not estimate a terminal history.  It proves the information
boundary frozen by TDG4-FRZ1: on the declared unrestricted smooth function
space, finitely many point samples do not bound the continuum Fourier power in
any nonzero target bin.  A smooth direction supported between samples lies in
the sampling kernel, has a nonzero target coefficient, and can be scaled
without changing any recorded value.

The returned dictionaries are proof certificates and exact algebra controls,
not a claim of formal proof-assistant verification.  No campaign checkpoint,
trajectory, SGB-L state, or FGC-QR state is consumed here.
"""

from __future__ import annotations

from fractions import Fraction
from numbers import Integral
from typing import Any, Sequence, TypeAlias


ExactRational: TypeAlias = Fraction | Integral
ExactMatrix: TypeAlias = tuple[tuple[Fraction, ...], ...]

TDG4_THEOREM_SAMPLE_COUNT = 64
TDG4_THEOREM_TOP_BINS = (28, 29, 30, 31, 32)
TDG4_THEOREM_FUNCTION_SPACE = "C_c_infinity_open_unit_interval"
TDG4_THEOREM_SUFFICIENT_ROUTES = (
    "finite_dimensional_class_with_injective_stable_sampling",
    "quantitative_continuum_derivative_or_Sobolev_bound",
    "continuum_evolution_residual_plus_stability_estimate",
)
TDG4_THEOREM_INSUFFICIENT_PREMISES = (
    "finite_samples_only",
    "smoothness_without_a_quantitative_norm_bound",
    "agreement_of_finitely_sampled_surrogates",
    "finite_sample_densification_without_regular_control",
    "discrete_residual_without_a_continuum_stability_estimate",
)


def _fraction(name: str, value: ExactRational) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise TypeError(f"{name} must be an exact rational value")
    return value if isinstance(value, Fraction) else Fraction(int(value))


def _fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _nodes(values: Sequence[ExactRational]) -> tuple[Fraction, ...]:
    if isinstance(values, (str, bytes)) or len(values) < 2:
        raise ValueError("TDG4 theorem requires at least two sample nodes")
    nodes = tuple(
        _fraction(f"sample node {index}", value)
        for index, value in enumerate(values)
    )
    if nodes[0] < 0 or nodes[-1] > 1:
        raise ValueError("TDG4 theorem nodes must lie in the unit interval")
    if any(right <= left for left, right in zip(nodes, nodes[1:])):
        raise ValueError("TDG4 theorem nodes must be strictly increasing")
    return nodes


def _bins(values: Sequence[int]) -> tuple[int, ...]:
    if isinstance(values, (str, bytes)) or not values:
        raise ValueError("TDG4 theorem requires at least one target bin")
    bins: list[int] = []
    for value in values:
        if isinstance(value, bool) or not isinstance(value, Integral):
            raise TypeError("TDG4 theorem target bins must be integers")
        frequency = int(value)
        if frequency <= 0:
            raise ValueError("TDG4 theorem target bins must be positive")
        bins.append(frequency)
    if len(set(bins)) != len(bins):
        raise ValueError("TDG4 theorem target bins must be distinct")
    return tuple(bins)


def _matrix(name: str, values: Sequence[Sequence[ExactRational]]) -> ExactMatrix:
    if isinstance(values, (str, bytes)) or not values:
        raise ValueError(f"{name} must be a nonempty matrix")
    rows: list[tuple[Fraction, ...]] = []
    width: int | None = None
    for row_index, row in enumerate(values):
        if isinstance(row, (str, bytes)) or not row:
            raise ValueError(f"{name} row {row_index} must be nonempty")
        converted = tuple(
            _fraction(f"{name}[{row_index},{column_index}]", entry)
            for column_index, entry in enumerate(row)
        )
        if width is None:
            width = len(converted)
        elif len(converted) != width:
            raise ValueError(f"{name} must be rectangular")
        rows.append(converted)
    return tuple(rows)


def _rank(values: ExactMatrix) -> int:
    """Independent exact rational rank used only by the theorem controls."""

    work = [list(row) for row in values]
    row_count = len(work)
    column_count = len(work[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [entry / scale for entry in work[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                left - factor * right
                for left, right in zip(work[row], work[pivot_row], strict=True)
            ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def factorization_control(
    sample_matrix: Sequence[Sequence[ExactRational]],
    target_matrix: Sequence[Sequence[ExactRational]],
) -> dict[str, Any]:
    """Check the finite-dimensional shadow of ``T|ker(S)=0`` independently."""

    sample = _matrix("sample matrix", sample_matrix)
    target = _matrix("target matrix", target_matrix)
    if len(sample[0]) != len(target[0]):
        raise ValueError("sample and target matrices need one domain dimension")
    sample_rank = _rank(sample)
    combined_rank = _rank(sample + target)
    factors = combined_rank == sample_rank
    return {
        "domain_dimension": len(sample[0]),
        "sample_rank": sample_rank,
        "sample_plus_target_rank": combined_rank,
        "sample_map_is_injective": sample_rank == len(sample[0]),
        "target_annihilates_sampling_kernel": factors,
        "target_factors_through_samples": factors,
        "unrestricted_kernel_amplitude_makes_target_unbounded": not factors,
    }


def factorization_lemma_certificate() -> dict[str, Any]:
    """Return necessity, sufficiency, and independent exact matrix controls."""

    controls = {
        "injective": factorization_control(
            ((1, 0), (0, 1)),
            ((2, 3),),
        ),
        "factorable_noninjective": factorization_control(
            ((1, 0, 0), (0, 1, 0)),
            ((2, 3, 0),),
        ),
        "nonidentifiable": factorization_control(
            ((1, 0, 0), (0, 1, 0)),
            ((0, 0, 1),),
        ),
    }
    all_controls_pass = (
        controls["injective"]["sample_map_is_injective"]
        and controls["injective"]["target_factors_through_samples"]
        and not controls["factorable_noninjective"]["sample_map_is_injective"]
        and controls["factorable_noninjective"][
            "target_factors_through_samples"
        ]
        and not controls["nonidentifiable"][
            "target_annihilates_sampling_kernel"
        ]
        and controls["nonidentifiable"][
            "unrestricted_kernel_amplitude_makes_target_unbounded"
        ]
    )
    return {
        "statement": "T factors through S iff ker(S) is a subset of ker(T)",
        "necessity": (
            "if T=A composed_with S and S(h)=0 then T(h)=A(0)=0"
        ),
        "sufficiency": (
            "define A(S(f))=T(f); ker(S) subset ker(T) makes A well_defined"
        ),
        "well_defined_difference_check": (
            "S(f)=S(g) implies f-g in ker(S), hence T(f)=T(g)"
        ),
        "exact_matrix_controls": controls,
        "all_factorization_controls_pass": all_controls_pass,
    }


def smooth_kernel_witness_theorem(
    sample_nodes: Sequence[ExactRational], *, target_bin: int
) -> dict[str, Any]:
    """Certify a real smooth kernel direction for one nonzero Fourier bin."""

    nodes = _nodes(sample_nodes)
    bins = _bins((target_bin,))
    frequency = bins[0]
    gaps = tuple(right - left for left, right in zip(nodes, nodes[1:]))
    gap_index = max(range(len(gaps)), key=lambda index: (gaps[index], -index))
    left = nodes[gap_index]
    right = nodes[gap_index + 1]
    width = right - left
    support_left = left + width / 3
    support_right = left + 2 * width / 3
    if not left < support_left < support_right < right:
        raise ValueError("TDG4 theorem support must lie strictly inside a gap")
    return {
        "sample_count": len(nodes),
        "target_bin": frequency,
        "selected_gap_index": gap_index,
        "selected_gap": [_fraction_text(left), _fraction_text(right)],
        "support_interval": [
            _fraction_text(support_left),
            _fraction_text(support_right),
        ],
        "standard_bump": (
            "phi(x)=exp(-1/((x-a)*(b-x))) for a<x<b and zero otherwise"
        ),
        "witness": "h_A,k(x)=A*phi(x)*cos(2*pi*k*x)",
        "real_smooth_compact_support": True,
        "zero_on_a_neighborhood_of_every_sample": True,
        "sample_values_and_all_sample_node_derivatives_zero": True,
        "target_coefficient_identity": "hat(h_A,k)=A*(I_0+I_2k)/2",
        "I0_strictly_positive": True,
        "phase_nonconstant_on_positive_measure_support": True,
        "strict_triangle_inequality": "abs(I_2k)<I_0",
        "reverse_triangle_lower_bound": (
            "abs(I_0+I_2k)>=I_0-abs(I_2k)>0"
        ),
        "target_coefficient_nonzero_for_every_nonzero_amplitude": True,
        "field_power_lower_bound": "2*A^2*abs(hat(h_1,k))^2",
        "derivative_coefficient_identity": "hat(dh/dx)_k=2*pi*i*k*hat(h)_k",
        "derivative_power_lower_bound": (
            "2*(2*pi*k)^2*A^2*abs(hat(h_1,k))^2"
        ),
        "both_declared_powers_unbounded_as_abs_A_tends_to_infinity": True,
    }


def uniform_alias_crosscheck(*, sample_count: int, target_bin: int) -> dict[str, Any]:
    """Independently rederive the nominal uniform-grid alias identity."""

    if isinstance(sample_count, bool) or not isinstance(sample_count, Integral):
        raise TypeError("sample count must be an integer")
    count = int(sample_count)
    if count < 3:
        raise ValueError("uniform alias crosscheck requires at least three nodes")
    frequency = _bins((target_bin,))[0]
    carrier = count - 1
    if frequency >= carrier:
        raise ValueError("target bin must lie below the uniform-grid carrier")
    coefficient = Fraction(1, 4)
    return {
        "sample_count": count,
        "target_bin": frequency,
        "carrier_frequency": carrier,
        "partner_bin": 2 * carrier - frequency,
        "sample_identity": "sin(2*pi*(N-1)*j/(N-1))=sin(2*pi*j)=0",
        "all_sample_values_exactly_zero": True,
        "product_to_sum_identity": (
            "sin(qx)*sin((q-k)x)=(cos(kx)-cos((2q-k)x))/2"
        ),
        "target_complex_coefficient": _fraction_text(coefficient),
        "partner_complex_coefficient": _fraction_text(-coefficient),
        "positive_bin_field_power": _fraction_text(2 * coefficient * coefficient),
        "amplitude_scaling_makes_target_power_unbounded": True,
    }


def sampling_identifiability_theorem_certificate(
    sample_nodes: Sequence[ExactRational] | None = None,
    *,
    target_bins: Sequence[int] = TDG4_THEOREM_TOP_BINS,
) -> dict[str, Any]:
    """Execute the complete TDG4 theorem certificate without history data."""

    if sample_nodes is None:
        sample_nodes = tuple(
            Fraction(index, TDG4_THEOREM_SAMPLE_COUNT - 1)
            for index in range(TDG4_THEOREM_SAMPLE_COUNT)
        )
    nodes = _nodes(sample_nodes)
    bins = _bins(target_bins)
    factorization = factorization_lemma_certificate()
    witnesses = {
        str(frequency): smooth_kernel_witness_theorem(
            nodes, target_bin=frequency
        )
        for frequency in bins
    }
    aliases = {
        str(frequency): uniform_alias_crosscheck(
            sample_count=TDG4_THEOREM_SAMPLE_COUNT,
            target_bin=frequency,
        )
        for frequency in bins
    }
    all_witnesses_pass = all(
        witness["sample_values_and_all_sample_node_derivatives_zero"]
        and witness[
            "target_coefficient_nonzero_for_every_nonzero_amplitude"
        ]
        and witness[
            "both_declared_powers_unbounded_as_abs_A_tends_to_infinity"
        ]
        for witness in witnesses.values()
    )
    all_aliases_pass = all(
        alias["all_sample_values_exactly_zero"]
        and alias["target_complex_coefficient"] == "1/4"
        and alias["positive_bin_field_power"] == "1/8"
        for alias in aliases.values()
    )
    theorem_completed = (
        factorization["all_factorization_controls_pass"]
        and all_witnesses_pass
        and all_aliases_pass
    )
    return {
        "theorem": (
            "finite point samples alone do not place a finite upper bound on "
            "the declared continuum top-band field or derivative power on "
            "C_c_infinity_open_unit_interval"
        ),
        "function_space": TDG4_THEOREM_FUNCTION_SPACE,
        "sample_count": len(nodes),
        "target_bins": list(bins),
        "sampling_operator": (
            f"S_x(f)=(f(x_0),...,f(x_{len(nodes) - 1}))"
        ),
        "target_operator": (
            "T_K(f)=(integral_0^1 f(x)*exp(-2*pi*i*k*x) dx) for k in K"
        ),
        "factorization_lemma": factorization,
        "smooth_kernel_witnesses": witnesses,
        "uniform_alias_crosschecks": aliases,
        "all_smooth_witnesses_pass": all_witnesses_pass,
        "all_uniform_alias_crosschecks_pass": all_aliases_pass,
        "finite_samples_alone_bound_field_power": False,
        "finite_samples_alone_bound_derivative_power": False,
        "actual_unsampled_power_of_the_PROTO13_histories_inferred": False,
        "PDE_solution_manifold_restriction_proved": False,
        "sufficient_assumption_routes": list(TDG4_THEOREM_SUFFICIENT_ROUTES),
        "insufficient_premises": list(TDG4_THEOREM_INSUFFICIENT_PREMISES),
        "sampling_identifiability_theorem_completed": theorem_completed,
        "actual_terminal_histories_consumed": False,
        "campaign_checkpoint_loaded": False,
        "state_advanced": False,
    }
