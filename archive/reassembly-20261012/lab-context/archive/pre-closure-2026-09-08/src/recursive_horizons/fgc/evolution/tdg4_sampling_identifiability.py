"""Exact sampling-identifiability controls for the prospective TDG4 gate.

TDG3 showed that two native-grid estimators can agree while still measuring
only sampled surrogates.  This module isolates the underlying mathematical
question.  Let ``S`` evaluate a continuum history at finitely many times and
let ``T`` return its declared top-band Fourier coefficients.  Samples alone
identify ``T`` exactly if and only if ``T`` annihilates ``ker(S)``.  If a
single smooth kernel direction has nonzero target coefficient, then its
amplitude can be scaled without changing a sample and the target power is
unbounded on the sample fibre.

The universal witness used here is supported strictly between adjacent sample
times.  An independent, exact trigonometric witness covers the nominal
64-point uniform phase grid.  Neither construction reads a campaign history,
defines a replacement admission, or says anything about FGC-QR dynamics.
"""

from __future__ import annotations

from fractions import Fraction
from numbers import Integral
from typing import Any, Sequence, TypeAlias

from ..exact_linear_algebra import rank


ExactRational: TypeAlias = Fraction | Integral
ExactMatrix: TypeAlias = tuple[tuple[Fraction, ...], ...]

TDG4_SAMPLE_COUNT = 64
TDG4_TOP_BINS = (28, 29, 30, 31, 32)
TDG4_FUNCTION_SPACE = "C_c_infinity_open_unit_interval"
TDG4_SUFFICIENT_ASSUMPTION_ROUTES = (
    "finite_dimensional_class_with_injective_stable_sampling",
    "quantitative_continuum_derivative_or_Sobolev_bound",
    "continuum_evolution_residual_plus_stability_estimate",
)
TDG4_INSUFFICIENT_PREMISES = (
    "finite_samples_only",
    "smoothness_without_a_quantitative_norm_bound",
    "agreement_of_finitely_sampled_surrogates",
    "finite_sample_densification_without_regular_control",
    "discrete_residual_without_a_continuum_stability_estimate",
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _fraction(name: str, value: ExactRational) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise TypeError(f"{name} must be an exact rational value")
    return value if isinstance(value, Fraction) else Fraction(int(value))


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _exact_nodes(values: Sequence[ExactRational]) -> tuple[Fraction, ...]:
    if isinstance(values, (str, bytes)) or len(values) < 2:
        raise ValueError("TDG4 requires at least two sample nodes")
    nodes = tuple(_fraction(f"sample node {index}", value) for index, value in enumerate(values))
    if nodes[0] < 0 or nodes[-1] > 1:
        raise ValueError("TDG4 sample nodes must lie in the unit interval")
    if any(right <= left for left, right in zip(nodes, nodes[1:])):
        raise ValueError("TDG4 sample nodes must be strictly increasing")
    return nodes


def _exact_matrix(name: str, value: Sequence[Sequence[ExactRational]]) -> ExactMatrix:
    if isinstance(value, (str, bytes)) or not value:
        raise ValueError(f"{name} must be a nonempty matrix")
    rows: list[tuple[Fraction, ...]] = []
    width: int | None = None
    for row_index, row in enumerate(value):
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


def generic_smooth_kernel_witness(
    sample_nodes: Sequence[ExactRational], *, target_bin: int
) -> dict[str, Any]:
    """Construct the exact support geometry for a smooth kernel witness.

    For the widest adjacent-node gap ``(a,b)``, choose an interior support
    interval ``(a + (b-a)/3, a + 2(b-a)/3)`` and any nonnegative, nonzero
    ``phi in C_c^infinity`` on that interval.  Then

    ``h_A,k(x) = A phi(x) cos(2 pi k x)``

    vanishes in a neighborhood of every sample.  Its coefficient at ``k`` is
    ``A (I_0 + I_2k)/2``.  Strict triangle inequality gives
    ``|I_2k| < I_0`` because the phase is nonconstant on a positive-measure
    interval, so the coefficient is nonzero and scales linearly with ``A``.
    """

    nodes = _exact_nodes(sample_nodes)
    if isinstance(target_bin, bool) or not isinstance(target_bin, Integral):
        raise TypeError("TDG4 target bin must be an integer")
    frequency = int(target_bin)
    if frequency <= 0:
        raise ValueError("TDG4 target bin must be positive")
    gaps = tuple(right - left for left, right in zip(nodes, nodes[1:]))
    gap_index = max(range(len(gaps)), key=lambda index: (gaps[index], -index))
    left = nodes[gap_index]
    right = nodes[gap_index + 1]
    width = right - left
    support_left = left + width / 3
    support_right = left + 2 * width / 3
    _require(
        left < support_left < support_right < right,
        "TDG4 witness support is not strictly inside its sample gap",
    )
    return {
        "sample_count": len(nodes),
        "target_bin": frequency,
        "function_space": TDG4_FUNCTION_SPACE,
        "sample_gap_index": gap_index,
        "sample_gap": [_fraction_text(left), _fraction_text(right)],
        "sample_gap_width": _fraction_text(width),
        "support_interval": [
            _fraction_text(support_left),
            _fraction_text(support_right),
        ],
        "witness_family": "h_A,k(x)=A*phi(x)*cos(2*pi*k*x)",
        "bump_contract": (
            "phi is nonnegative, nonzero, smooth, and compactly supported "
            "inside the declared open interval"
        ),
        "all_sample_values_are_exactly_zero": True,
        "all_sample_neighborhoods_are_zero": True,
        "all_derivatives_at_sample_nodes_are_zero": True,
        "target_coefficient_identity": "hat(h_A,k)=A*(I_0+I_2k)/2",
        "strict_triangle_inequality": "abs(I_2k)<I_0",
        "strictness_reason": (
            "the Fourier phase is nonconstant on positive-measure bump support"
        ),
        "target_coefficient_nonzero_for_nonzero_amplitude": True,
        "field_top_band_power_scales_as_amplitude_squared": True,
        "derivative_coefficient_identity": "hat(dh/dx)_k=2*pi*i*k*hat(h)_k",
        "derivative_top_band_power_scales_as_amplitude_squared": True,
        "sample_fibre_target_power_is_unbounded": True,
        "finite_densification_alone_remains_nonidentifying": True,
    }


def uniform_alias_witness(*, sample_count: int, target_bin: int) -> dict[str, Any]:
    """Return an exact trigonometric alias witness on ``j/(N-1)`` nodes.

    With ``q=N-1``, the function

    ``sin(2*pi*q*x) sin(2*pi*(q-k)*x)``

    vanishes at every sample.  Product-to-sum gives coefficients ``1/4`` at
    ``+/-k`` and ``-1/4`` at ``+/-(2q-k)``.  The declared positive-bin field
    power at unit amplitude is therefore exactly ``2*(1/4)^2 = 1/8``.
    """

    if isinstance(sample_count, bool) or not isinstance(sample_count, Integral):
        raise TypeError("TDG4 sample count must be an integer")
    count = int(sample_count)
    if count < 3:
        raise ValueError("TDG4 uniform alias control requires at least three samples")
    if isinstance(target_bin, bool) or not isinstance(target_bin, Integral):
        raise TypeError("TDG4 target bin must be an integer")
    frequency = int(target_bin)
    carrier = count - 1
    if frequency <= 0 or frequency >= carrier:
        raise ValueError("TDG4 target bin must lie strictly below the carrier")
    partner = 2 * carrier - frequency
    coefficient = Fraction(1, 4)
    return {
        "sample_count": count,
        "sample_node_formula": "x_j=j/(N-1), j=0,...,N-1",
        "carrier_frequency": carrier,
        "target_bin": frequency,
        "partner_bin": partner,
        "witness": "sin(2*pi*(N-1)*x)*sin(2*pi*((N-1)-k)*x)",
        "sample_annihilation_identity": "sin(2*pi*j)=0",
        "all_sample_values_are_exactly_zero": True,
        "target_complex_coefficient": _fraction_text(coefficient),
        "negative_target_complex_coefficient": _fraction_text(coefficient),
        "partner_complex_coefficient": _fraction_text(-coefficient),
        "negative_partner_complex_coefficient": _fraction_text(-coefficient),
        "target_coefficient_squared": _fraction_text(coefficient * coefficient),
        "unit_amplitude_positive_bin_field_power": _fraction_text(
            2 * coefficient * coefficient
        ),
        "unit_amplitude_positive_bin_derivative_power": (
            f"(2*pi*{frequency})^2/8"
        ),
        "same_sample_vector_for_every_amplitude": True,
        "target_power_unbounded_under_amplitude_scaling": True,
    }


def finite_dimensional_sampling_assessment(
    sample_matrix: Sequence[Sequence[ExactRational]],
    target_matrix: Sequence[Sequence[ExactRational]],
) -> dict[str, Any]:
    """Assess exact linear identifiability with rational row-space ranks.

    ``T`` annihilates ``ker(S)`` exactly when every target row belongs to the
    row space of ``S``.  Equivalently, stacking the two matrices does not
    increase rank.  In that case ``T=A*S`` on the declared finite-dimensional
    class.  Otherwise an unrestricted kernel amplitude makes the target
    unbounded on a sample fibre.
    """

    sample = _exact_matrix("sample matrix", sample_matrix)
    target = _exact_matrix("target matrix", target_matrix)
    if len(sample[0]) != len(target[0]):
        raise ValueError("TDG4 sample and target matrices need one domain dimension")
    sample_rank = rank(sample)
    augmented_rank = rank(sample + target)
    factors = augmented_rank == sample_rank
    domain_dimension = len(sample[0])
    return {
        "domain_dimension": domain_dimension,
        "sample_row_count": len(sample),
        "target_row_count": len(target),
        "sample_rank": sample_rank,
        "sample_plus_target_rank": augmented_rank,
        "sample_map_is_injective": sample_rank == domain_dimension,
        "target_annihilates_sampling_kernel": factors,
        "target_factors_through_samples": factors,
        "unrestricted_sample_fibre_target_is_unbounded": not factors,
        "quantitative_stability_constant_still_required_for_noisy_data": True,
    }


def lipschitz_coefficient_error_bound(
    sample_nodes: Sequence[ExactRational], *, derivative_bound: ExactRational
) -> dict[str, Any]:
    """Bound Fourier-coefficient error from a declared continuum ``C^1`` norm.

    If ``|f'| <= L`` and ``f_PL`` interpolates the exact samples, then on a
    gap of width ``Delta`` the interpolation error has integral at most
    ``L*Delta^2/2``.  Consequently every Fourier coefficient obeys

    ``|hat(f)_k-hat(f_PL)_k| <= L/2 * sum_j Delta_j^2 <= L*h_max/2``.

    The result is useful only when ``L`` is supplied independently of the
    samples being bounded.
    """

    nodes = _exact_nodes(sample_nodes)
    if nodes[0] != 0 or nodes[-1] != 1:
        raise ValueError("TDG4 coefficient bound requires nodes spanning [0,1]")
    bound = _fraction("derivative bound", derivative_bound)
    if bound < 0:
        raise ValueError("TDG4 derivative bound must be nonnegative")
    gaps = tuple(right - left for left, right in zip(nodes, nodes[1:]))
    sum_squares = sum((gap * gap for gap in gaps), Fraction(0))
    maximum_gap = max(gaps)
    sharp = bound * sum_squares / 2
    coarse = bound * maximum_gap / 2
    _require(sharp <= coarse, "TDG4 interpolation bounds are inconsistent")
    return {
        "sample_count": len(nodes),
        "derivative_bound": _fraction_text(bound),
        "sum_gap_squares": _fraction_text(sum_squares),
        "maximum_gap": _fraction_text(maximum_gap),
        "coefficient_error_bound": _fraction_text(sharp),
        "coarser_max_gap_bound": _fraction_text(coarse),
        "bound_is_frequency_uniform": True,
        "same_coefficient_bound_feeds_field_and_phase_derivative_power": True,
        "derivative_bound_is_an_external_quantitative_premise": True,
        "samples_alone_do_not_supply_the_derivative_bound": True,
    }


def sampling_identifiability_preflight() -> dict[str, Any]:
    """Run exact no-history controls for the prospective TDG4 theorem."""

    uniform_nodes = tuple(Fraction(index, TDG4_SAMPLE_COUNT - 1) for index in range(TDG4_SAMPLE_COUNT))
    smooth = {
        str(frequency): generic_smooth_kernel_witness(
            uniform_nodes, target_bin=frequency
        )
        for frequency in TDG4_TOP_BINS
    }
    aliases = {
        str(frequency): uniform_alias_witness(
            sample_count=TDG4_SAMPLE_COUNT, target_bin=frequency
        )
        for frequency in TDG4_TOP_BINS
    }
    nonidentifiable = finite_dimensional_sampling_assessment(
        ((1, 0, 0), (0, 1, 0)),
        ((0, 0, 1),),
    )
    factorable = finite_dimensional_sampling_assessment(
        ((1, 0, 0), (0, 1, 0)),
        ((2, 3, 0),),
    )
    injective = finite_dimensional_sampling_assessment(
        ((1, 0), (0, 1)),
        ((2, 3),),
    )
    lipschitz = lipschitz_coefficient_error_bound(
        (Fraction(0), Fraction(1, 3), Fraction(1)),
        derivative_bound=Fraction(6),
    )
    all_exact_controls_pass = (
        all(record["sample_fibre_target_power_is_unbounded"] for record in smooth.values())
        and all(record["unit_amplitude_positive_bin_field_power"] == "1/8" for record in aliases.values())
        and nonidentifiable["unrestricted_sample_fibre_target_is_unbounded"]
        and factorable["target_factors_through_samples"]
        and injective["sample_map_is_injective"]
        and lipschitz["coefficient_error_bound"] == "5/3"
        and lipschitz["coarser_max_gap_bound"] == "2"
    )
    _require(all_exact_controls_pass, "TDG4 exact preflight differs")
    return {
        "sample_count": TDG4_SAMPLE_COUNT,
        "top_bins": list(TDG4_TOP_BINS),
        "function_space": TDG4_FUNCTION_SPACE,
        "generic_smooth_kernel_witnesses": smooth,
        "exact_uniform_alias_witnesses": aliases,
        "finite_dimensional_controls": {
            "nonidentifiable": nonidentifiable,
            "factorable_noninjective": factorable,
            "injective": injective,
        },
        "quantitative_regular_bound_control": lipschitz,
        "sufficient_assumption_routes": list(TDG4_SUFFICIENT_ASSUMPTION_ROUTES),
        "insufficient_premises": list(TDG4_INSUFFICIENT_PREMISES),
        "all_exact_controls_pass": True,
        "actual_terminal_histories_consumed": False,
        "campaign_checkpoint_loaded": False,
        "state_advanced": False,
        "replacement_temporal_admission_defined": False,
        "PROTO14_frozen": False,
    }
