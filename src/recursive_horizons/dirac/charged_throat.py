"""Bind the NSC charged subsector to the imported MMP throat equations.

This module does not solve the Maldacena--Milekhin--Popov geometry.  It maps
the already calculated NSC field content, coefficients and state data into
their normalization and reports the first unresolved common-action input.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from math import pi, sqrt
from typing import Any, Mapping

import numpy as np


@dataclass(frozen=True)
class CommonCoefficients:
    """Local Lorentzian coefficients in integral sqrt(-g)(A R-C F^2-V)."""

    einstein: float
    gauge: float
    vacuum: float

    def __post_init__(self) -> None:
        values = (self.einstein, self.gauge, self.vacuum)
        if not all(np.isfinite(values)) or self.einstein <= 0 or self.gauge <= 0:
            raise ValueError("finite positive Einstein and gauge coefficients required")

    def mmp_map(self, flux: int = 1) -> dict[str, float]:
        if isinstance(flux, bool) or not isinstance(flux, int) or flux == 0:
            raise ValueError("nonzero integer magnetic flux required")
        q = abs(flux)
        newton = 1.0 / (16.0 * pi * self.einstein)
        gauge_squared = 1.0 / (4.0 * self.gauge)
        charge_radius_squared = q * q * self.gauge / (4.0 * self.einstein)
        charge_radius = sqrt(charge_radius_squared)
        throat_length = 16.0 * charge_radius**3 / (q * newton)
        cosmological_curvature = self.vacuum / (2.0 * self.einstein)
        xi = cosmological_curvature * charge_radius_squared
        return {
            "G_N": newton,
            "g_squared": gauge_squared,
            "r_e_squared": charge_radius_squared,
            "ell_MMP": throat_length,
            "lambda_4": cosmological_curvature,
            "Xi": xi,
            "Xi_per_q_squared": xi / (q * q),
        }


@dataclass(frozen=True)
class MMPEmbeddingResult:
    """Result of the algebraic/domain embedding gate."""

    field_content_passes: bool
    anomaly_and_link_parity_passes: bool
    antiperiodic_state_passes: bool
    coefficient_normalization_passes: bool
    retained_partial_seed_passes: bool
    full_embedding_decidable: bool
    first_unresolved_owner: str | None
    coefficient_map: Mapping[str, float]
    domain_map: Mapping[str, Any]
    state_map: Mapping[str, Any]
    residuals: Mapping[str, Any]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _matching_row(compact_matching: Mapping[str, Any], matching_cutoff: float) -> Mapping[str, Any]:
    rows = compact_matching.get("matched_coefficients", ())
    selected = [row for row in rows if abs(row["matching_cutoff"] - matching_cutoff) < 1e-14]
    if len(selected) != 1:
        raise ValueError("one compact-matching coefficient row required")
    return selected[0]


def evaluate_mmp_embedding(
    charged_sector: Mapping[str, Any],
    compact_matching: Mapping[str, Any],
    compact_casimir: Mapping[str, Any],
    vacuum_matching: Mapping[str, Any],
    *,
    matching_cutoff: float = 1.0,
) -> MMPEmbeddingResult:
    """Evaluate the narrow imported-solution embedding without a gravity run."""
    schemas = {
        charged_sector.get("schema"), compact_matching.get("schema"),
        compact_casimir.get("schema"), vacuum_matching.get("schema"),
    }
    expected = {
        "NSC-CHARGED-SECTOR-v1", "NSC-COMPACT-MATCHING-v1",
        "NSC-COMPACT-CASIMIR-v1", "NSC-VACUUM-CHARGE-MATCHING-v2",
    }
    if schemas != expected:
        raise ValueError("embedding inputs have unexpected schemas")

    multiplicity = charged_sector["multiplicity"]
    field_pass = (
        multiplicity["bulk_complex_Dirac_fields"] == 2
        and multiplicity["massless_4D_Dirac_zero_fields_before_link_mass"] == 1
        and multiplicity["massive_4D_Dirac_fields_per_nonzero_compact_level"] == 2
    )
    admissible = [
        row for row in charged_sector["parity_and_anomaly_rows"]
        if row["anomaly_free_zero_sector"] and row["neutral_even_scalar_link_allowed"]
    ]
    parity_pass = len(admissible) == 2 and all(
        row["left_Weyl_zero_modes"] == row["right_Weyl_zero_modes"] == 1
        for row in admissible
    )
    holonomy = compact_casimir["holonomy"]
    ap_pass = (
        "alpha=1/2" in holonomy["minimum"]
        and holonomy["second_variation_control"]["matrix_hessian"] > 0
    )

    row = _matching_row(compact_matching, matching_cutoff)
    coefficients = CommonCoefficients(
        einstein=row["A_Dirac"], gauge=row["C_gauge_Dirac"], vacuum=row["V_Dirac"]
    )
    mapped = coefficients.mmp_map()
    coefficient_definitions = vacuum_matching["coefficient_definitions"]
    coefficient_pass = (
        coefficient_definitions["GN"] == "1/(16*pi*A)"
        and coefficient_definitions["gauge_coupling_squared"] == "1/(4*C)"
        and coefficient_definitions["rQ_squared"] == "q^2*C/(4*A)"
    )
    retained_seed_pass = mapped["Xi_per_q_squared"] <= 0.25
    scope = vacuum_matching["scope"]
    full_decidable = bool(
        scope["complete_quantum_measure_or_vacuum_coefficient_fixed"]
        and compact_casimir["holonomy"]["complete_state_or_recursive_return_path_derived"]
    )
    unresolved = None
    if not retained_seed_pass and not scope["complete_quantum_measure_or_vacuum_coefficient_fixed"]:
        unresolved = "full common-action vacuum coefficient V_full"
    elif not compact_casimir["holonomy"]["complete_state_or_recursive_return_path_derived"]:
        unresolved = "physical closed magnetic return domain and state"
    elif not full_decidable:
        unresolved = "complete common-action embedding"

    return MMPEmbeddingResult(
        field_content_passes=field_pass,
        anomaly_and_link_parity_passes=parity_pass,
        antiperiodic_state_passes=ap_pass,
        coefficient_normalization_passes=coefficient_pass,
        retained_partial_seed_passes=retained_seed_pass,
        full_embedding_decidable=full_decidable,
        first_unresolved_owner=unresolved,
        coefficient_map=mapped,
        domain_map={
            "reference": "MMP closed magnetic field-line return path",
            "nsc_unwrapped_radial_line_identified_with_reference": False,
            "admissible_compact_parities": [
                [row["parent_parity_sign"], row["child_parity_sign"]]
                for row in admissible
            ],
            "magnetic_flux": "nonzero integer q retained symbolically",
            "lowest_Landau_channels": "abs(q) complex two-dimensional Dirac fields",
        },
        state_map={
            "phase": "alpha=1/2 modulo one",
            "selection": holonomy["minimum"],
            "AP_minus_P_energy_control": holonomy["separate_resolution_controls"][-1]["AP_minus_P_energy"],
            "complete_physical_return_state": holonomy["complete_state_or_recursive_return_path_derived"],
        },
        residuals={
            "coefficient_map": 0,
            "field_multiplicity": 0 if field_pass else 1,
            "parity_and_anomaly": 0 if parity_pass else 1,
            "AP_stationarity": holonomy["cutoff_controls"][-1]["phases"][4]["phase_derivative"],
            "retained_partial_Xi_minus_RNdS_bound_per_q_squared": mapped["Xi_per_q_squared"] - 0.25,
            "asymptotically_flat_required_lambda_4": 0.0,
            "retained_partial_lambda_4": mapped["lambda_4"],
        },
    )
