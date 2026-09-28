"""Prospective successor product-box continuous-no-trap resource ladder.

This owner does not edit the depth-8/10/12 instrument.  It binds that
predecessor's nominal contract hash and then evaluates a separately frozen
depth-14/16/18 ladder with practical resource caps declared before any
Picard run.  Aggregate health, ``FRZ1``, ``PREF1``, execution, and holdout
remain false.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from types import MappingProxyType
from typing import Any, Mapping

from .sgb1_ctl1_initial_health import (
    CERTIFICATE_CHART_KIND,
    DECLARED_BASE_CELLS,
    DECLARED_C_DOMAIN,
    DECLARED_K_DOMAIN,
    DECLARED_LAMBDA_DOMAIN,
    DECLARED_MAX_PICARD_ITERATIONS,
    SGBLExactInitialSlice,
)
from .sgb1_ctl1_trap_refinement import (
    DECLARED_COMPACTNESS_STRICT_UPPER,
    DECLARED_LADDER_DEPTHS as PREDECESSOR_DECLARED_DEPTHS,
    INSTRUMENT_ID as PREDECESSOR_INSTRUMENT_ID,
    NOMINAL_CHI_AMPLITUDE,
    PRODUCT_BOX_CHART_KIND,
    SMALL_AMPLITUDE_CONTROL,
    SGBLProductBoxLadderLevel,
    SGBLProductBoxLadderLevelRecord,
    sgbl_evaluate_product_box_ladder_level,
    sgbl_product_box_shared_coverage_monotone,
    sgbl_trap_refinement_contract_sha256,
)


Q = Fraction
INSTRUMENT_ID = "FGC-1-SGB1-CTL1-TRAP-REFINEMENT2"
PREDECESSOR_NOMINAL_CONTRACT_SHA256 = (
    "a82363a29b641c538a0e15f14ec2f4c99a49182a9a01c84ae7b1d14a953834c8"
)
DECLARED_LADDER2_DEPTHS = (14, 16, 18)
DECLARED_LADDER2_CELL_CAPS = (4096, 8192, 16384)
DECLARED_LADDER2_RHS_CAPS = (131072, 262144, 524288)
DECLARED_LADDER2_BIT_CAP = 65536
SUCCESSOR_STOP_REASONS = frozenset(
    {
        "wrong_predecessor",
        "changed_ladder",
        "reordered_levels",
        "skipped_level",
        "tampered_contract",
    }
)

if PREDECESSOR_INSTRUMENT_ID != "FGC-1-SGB1-CTL1-TRAP-REFINEMENT":
    raise RuntimeError("successor requires the preserved TRAP-REFINEMENT predecessor")
if PREDECESSOR_DECLARED_DEPTHS != (8, 10, 12):
    raise RuntimeError("successor requires the preserved depth-8/10/12 nonpass")
if set(DECLARED_LADDER2_DEPTHS) & set(PREDECESSOR_DECLARED_DEPTHS):
    raise RuntimeError("successor depths must not reuse the predecessor ladder")


def _fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


class SGBLTrapRefinement2Stop(ValueError):
    """Typed contract stop of the successor product-box ladder."""

    def __init__(
        self,
        reason: str,
        message: str,
        payload: Mapping[str, Any] | None = None,
    ) -> None:
        if reason not in SUCCESSOR_STOP_REASONS:
            raise ValueError("unknown trap-refinement-2 stop reason")
        self.reason = reason
        self.payload = dict(payload or {})
        super().__init__(message)


def _declare_product_box_resource_ladder2() -> tuple[SGBLProductBoxLadderLevel, ...]:
    levels = []
    for depth, cells, rhs in zip(
        DECLARED_LADDER2_DEPTHS,
        DECLARED_LADDER2_CELL_CAPS,
        DECLARED_LADDER2_RHS_CAPS,
    ):
        levels.append(
            SGBLProductBoxLadderLevel(
                name=f"depth_{depth}",
                max_bisection_depth=depth,
                max_cells=cells,
                max_rhs_evaluations=rhs,
                max_rational_bit_length=DECLARED_LADDER2_BIT_CAP,
            )
        )
    return tuple(levels)


DECLARED_PRODUCT_BOX_RESOURCE_LADDER2 = _declare_product_box_resource_ladder2()


def _require_declared_ladder2(
    ladder: object,
) -> tuple[SGBLProductBoxLadderLevel, ...]:
    if type(ladder) is not tuple:
        raise SGBLTrapRefinement2Stop(
            "changed_ladder",
            "successor resource ladder must be the frozen declared tuple",
        )
    if not ladder:
        raise SGBLTrapRefinement2Stop(
            "skipped_level",
            "declared successor ladder levels were skipped",
        )
    if not all(type(level) is SGBLProductBoxLadderLevel for level in ladder):
        raise SGBLTrapRefinement2Stop(
            "changed_ladder",
            "successor ladder levels must be SGBLProductBoxLadderLevel",
        )
    depths = tuple(level.max_bisection_depth for level in ladder)
    declared_set = set(DECLARED_LADDER2_DEPTHS)
    if len(ladder) < len(DECLARED_LADDER2_DEPTHS) and set(depths) <= declared_set:
        raise SGBLTrapRefinement2Stop(
            "skipped_level",
            "a declared successor ladder level was skipped",
            {"depths": depths},
        )
    if (
        len(ladder) == len(DECLARED_LADDER2_DEPTHS)
        and set(depths) == declared_set
        and depths != DECLARED_LADDER2_DEPTHS
    ):
        raise SGBLTrapRefinement2Stop(
            "reordered_levels",
            "declared successor ladder levels were reordered",
            {"depths": depths},
        )
    if ladder != DECLARED_PRODUCT_BOX_RESOURCE_LADDER2:
        raise SGBLTrapRefinement2Stop(
            "changed_ladder",
            "successor resource ladder is not the frozen declaration",
            {"depths": depths},
        )
    return ladder


def _require_predecessor_hash(value: object) -> str:
    if type(value) is not str or not value:
        raise SGBLTrapRefinement2Stop(
            "wrong_predecessor",
            "predecessor nominal contract hash must be the frozen digest",
        )
    if value != PREDECESSOR_NOMINAL_CONTRACT_SHA256:
        raise SGBLTrapRefinement2Stop(
            "wrong_predecessor",
            "successor does not bind the preserved depth-8/10/12 nominal hash",
            {"supplied": value, "required": PREDECESSOR_NOMINAL_CONTRACT_SHA256},
        )
    return value


def _shared_coverage_monotone(
    levels: tuple[SGBLProductBoxLadderLevelRecord, ...],
) -> bool:
    for earlier, later in zip(levels, levels[1:]):
        earlier_cells = earlier.ode.cells if earlier.ode is not None else ()
        later_cells = later.ode.cells if later.ode is not None else ()
        if not sgbl_product_box_shared_coverage_monotone(earlier_cells, later_cells):
            return False
    return True


def _compactness_upper_nonincreasing(
    levels: tuple[SGBLProductBoxLadderLevelRecord, ...],
) -> bool:
    previous: Fraction | None = None
    for record in levels:
        if record.ode is None or not record.ode.cells:
            continue
        if previous is not None and record.compactness_upper > previous:
            return False
        previous = record.compactness_upper
    return True


def _tiles_and_c_strictly_below_one(record: SGBLProductBoxLadderLevelRecord) -> bool:
    return (
        record.complete_support_coverage
        and record.compactness_upper < DECLARED_COMPACTNESS_STRICT_UPPER
        and record.exterior_compactness_upper < DECLARED_COMPACTNESS_STRICT_UPPER
        and record.compactness_margin > 0
    )


def _contract_payload2(
    spec: SGBLExactInitialSlice,
    ladder: tuple[SGBLProductBoxLadderLevel, ...],
    levels: tuple[SGBLProductBoxLadderLevelRecord, ...],
    *,
    shared_coverage_monotone: bool,
    compactness_upper_nonincreasing: bool,
    predecessor_nominal_contract_sha256: str,
) -> dict[str, Any]:
    proved_names = tuple(
        record.level.name
        for record in levels
        if record.continuous_no_initial_trapped_sphere
    )
    tiled_names = tuple(
        record.level.name for record in levels if record.complete_support_coverage
    )
    tiled_c_lt_1 = tuple(
        record.level.name for record in levels if _tiles_and_c_strictly_below_one(record)
    )
    return {
        "INSTRUMENT_ID": INSTRUMENT_ID,
        "predecessor_instrument_id": PREDECESSOR_INSTRUMENT_ID,
        "predecessor_declared_depths": list(PREDECESSOR_DECLARED_DEPTHS),
        "predecessor_nominal_contract_sha256": predecessor_nominal_contract_sha256,
        "chart_kind": PRODUCT_BOX_CHART_KIND,
        "certificate_chart_kind": CERTIFICATE_CHART_KIND,
        "product_box_is_not_the_certificate_chart": True,
        "chi_amplitude": _fraction_text(spec.chi_amplitude),
        "phi_amplitude": _fraction_text(spec.phi_amplitude),
        "center": _fraction_text(spec.center),
        "half_width": _fraction_text(spec.half_width),
        "compactness_gate": "strict_C_lt_1",
        "compactness_strict_upper": _fraction_text(Q(DECLARED_COMPACTNESS_STRICT_UPPER)),
        "base_cells": DECLARED_BASE_CELLS,
        "max_picard_iterations": DECLARED_MAX_PICARD_ITERATIONS,
        "lambda_domain": (
            _fraction_text(DECLARED_LAMBDA_DOMAIN.lower),
            _fraction_text(DECLARED_LAMBDA_DOMAIN.upper),
        ),
        "k_domain": (
            _fraction_text(DECLARED_K_DOMAIN.lower),
            _fraction_text(DECLARED_K_DOMAIN.upper),
        ),
        "compactness_domain": (
            _fraction_text(DECLARED_C_DOMAIN.lower),
            _fraction_text(DECLARED_C_DOMAIN.upper),
        ),
        "declared_depths": list(DECLARED_LADDER2_DEPTHS),
        "ladder": [level.as_contract_mapping() for level in ladder],
        "levels": [record.as_contract_mapping() for record in levels],
        "evaluated_every_declared_level": len(levels) == len(DECLARED_LADDER2_DEPTHS),
        "shared_coverage_monotone": shared_coverage_monotone,
        "compactness_upper_nonincreasing": compactness_upper_nonincreasing,
        "proved_level_names": list(proved_names),
        "tiled_support_level_names": list(tiled_names),
        "tiled_support_with_c_lt_1_level_names": list(tiled_c_lt_1),
        "any_declared_level_proves_continuous_no_initial_trap": bool(proved_names),
        "any_declared_level_tiles_compact_support": bool(tiled_names),
        "any_declared_level_tiles_and_c_strictly_below_one": bool(tiled_c_lt_1),
        "SGBL_branch_owned_and_healthy": False,
        "FRZ1": False,
        "PREF1": False,
        "holdout_authorized": False,
        "execution_authorized": False,
        "sampled_nodes_are_not_the_certificate": True,
    }


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLProductBoxLadder2Record:
    """Complete evaluation of the frozen successor product-box resource ladder."""

    slice: SGBLExactInitialSlice
    ladder: tuple[SGBLProductBoxLadderLevel, ...]
    levels: tuple[SGBLProductBoxLadderLevelRecord, ...]
    predecessor_nominal_contract_sha256: str
    shared_coverage_monotone: bool
    compactness_upper_nonincreasing: bool
    proved_level_names: tuple[str, ...]
    tiled_support_level_names: tuple[str, ...]
    tiled_support_with_c_lt_1_level_names: tuple[str, ...]
    any_declared_level_proves_continuous_no_initial_trap: bool
    any_declared_level_tiles_compact_support: bool
    any_declared_level_tiles_and_c_strictly_below_one: bool
    contract_payload: Mapping[str, Any]
    contract_sha256: str

    def __post_init__(self) -> None:
        if type(self.slice) is not SGBLExactInitialSlice:
            raise TypeError("slice must be SGBLExactInitialSlice")
        _require_predecessor_hash(self.predecessor_nominal_contract_sha256)
        if self.ladder != DECLARED_PRODUCT_BOX_RESOURCE_LADDER2:
            raise SGBLTrapRefinement2Stop(
                "changed_ladder",
                "successor aggregate record must keep the frozen resource ladder",
            )
        if len(self.levels) != len(DECLARED_LADDER2_DEPTHS):
            raise SGBLTrapRefinement2Stop(
                "skipped_level",
                "successor aggregate record must evaluate every declared level",
            )
        if tuple(record.level for record in self.levels) != DECLARED_PRODUCT_BOX_RESOURCE_LADDER2:
            raise SGBLTrapRefinement2Stop(
                "reordered_levels",
                "successor aggregate record levels must follow the frozen ladder order",
            )
        expected_monotone = _shared_coverage_monotone(self.levels)
        expected_nonincreasing = _compactness_upper_nonincreasing(self.levels)
        if self.shared_coverage_monotone != expected_monotone:
            raise SGBLTrapRefinement2Stop(
                "tampered_contract",
                "shared-coverage monotone bit does not match the returned cells",
            )
        if self.compactness_upper_nonincreasing != expected_nonincreasing:
            raise SGBLTrapRefinement2Stop(
                "tampered_contract",
                "compactness-upper monotone bit does not match the returned levels",
            )
        proved = tuple(
            record.level.name
            for record in self.levels
            if record.continuous_no_initial_trapped_sphere
        )
        tiled = tuple(
            record.level.name
            for record in self.levels
            if record.complete_support_coverage
        )
        tiled_c = tuple(
            record.level.name
            for record in self.levels
            if _tiles_and_c_strictly_below_one(record)
        )
        if self.proved_level_names != proved:
            raise SGBLTrapRefinement2Stop(
                "tampered_contract",
                "proved level names do not match the level records",
            )
        if self.tiled_support_level_names != tiled:
            raise SGBLTrapRefinement2Stop(
                "tampered_contract",
                "tiled-support level names do not match the level records",
            )
        if self.tiled_support_with_c_lt_1_level_names != tiled_c:
            raise SGBLTrapRefinement2Stop(
                "tampered_contract",
                "tiled C<1 level names do not match the level records",
            )
        if self.any_declared_level_proves_continuous_no_initial_trap != bool(proved):
            raise SGBLTrapRefinement2Stop(
                "tampered_contract",
                "aggregate pass bit does not match proved levels",
            )
        if self.any_declared_level_tiles_compact_support != bool(tiled):
            raise SGBLTrapRefinement2Stop(
                "tampered_contract",
                "aggregate tiling bit does not match tiled levels",
            )
        if self.any_declared_level_tiles_and_c_strictly_below_one != bool(tiled_c):
            raise SGBLTrapRefinement2Stop(
                "tampered_contract",
                "aggregate tiled-C<1 bit does not match tiled C<1 levels",
            )
        if proved and not self.shared_coverage_monotone:
            raise SGBLTrapRefinement2Stop(
                "tampered_contract",
                "a product-box pass cannot be claimed on nonmonotone shared coverage",
            )
        expected_payload = _contract_payload2(
            self.slice,
            self.ladder,
            self.levels,
            shared_coverage_monotone=self.shared_coverage_monotone,
            compactness_upper_nonincreasing=self.compactness_upper_nonincreasing,
            predecessor_nominal_contract_sha256=self.predecessor_nominal_contract_sha256,
        )
        if dict(self.contract_payload) != expected_payload:
            raise SGBLTrapRefinement2Stop(
                "tampered_contract",
                "contract payload does not match the frozen successor ladder record",
            )
        expected_hash = sgbl_trap_refinement_contract_sha256(expected_payload)
        if self.contract_sha256 != expected_hash:
            raise SGBLTrapRefinement2Stop(
                "tampered_contract",
                "contract hash does not match the nonpromoting payload",
            )
        if any(
            dict(self.contract_payload)[name] is not False
            for name in (
                "SGBL_branch_owned_and_healthy",
                "FRZ1",
                "PREF1",
                "holdout_authorized",
                "execution_authorized",
            )
        ):
            raise SGBLTrapRefinement2Stop(
                "tampered_contract",
                "nonpromoting contract payload cannot set aggregate health true",
            )

    @property
    def sampled_nodes_are_not_the_certificate(self) -> bool:
        return True

    @property
    def product_box_is_not_the_certificate_chart(self) -> bool:
        return True

    @property
    def SGBL_branch_owned_and_healthy(self) -> bool:
        return False

    @property
    def FRZ1(self) -> bool:
        return False

    @property
    def PREF1(self) -> bool:
        return False

    @property
    def holdout_authorized(self) -> bool:
        return False

    @property
    def execution_authorized(self) -> bool:
        return False

    @property
    def evaluated_every_declared_level(self) -> bool:
        return len(self.levels) == len(DECLARED_LADDER2_DEPTHS)

    @property
    def diagnostic_wall_seconds(self) -> tuple[float, ...]:
        return tuple(record.diagnostic_wall_seconds for record in self.levels)


def sgbl_evaluate_product_box_resource_ladder2(
    spec: SGBLExactInitialSlice | None = None,
    *,
    ladder: tuple[SGBLProductBoxLadderLevel, ...] | None = None,
    predecessor_nominal_contract_sha256: str | None = None,
) -> SGBLProductBoxLadder2Record:
    """Evaluate every frozen successor level.  A pass does not stop later levels."""

    if spec is None:
        spec = SGBLExactInitialSlice()
    elif type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    if spec.chi_amplitude != NOMINAL_CHI_AMPLITUDE and spec.chi_amplitude != SMALL_AMPLITUDE_CONTROL:
        raise ValueError(
            "successor ladder slice must be nominal A_chi=3 or the named small-amplitude control"
        )
    predecessor = _require_predecessor_hash(
        PREDECESSOR_NOMINAL_CONTRACT_SHA256
        if predecessor_nominal_contract_sha256 is None
        else predecessor_nominal_contract_sha256
    )
    frozen = _require_declared_ladder2(
        DECLARED_PRODUCT_BOX_RESOURCE_LADDER2 if ladder is None else ladder
    )
    records = []
    for level in frozen:
        records.append(sgbl_evaluate_product_box_ladder_level(spec, level))
    levels = tuple(records)
    monotone = _shared_coverage_monotone(levels)
    nonincreasing = _compactness_upper_nonincreasing(levels)
    proved = tuple(
        record.level.name for record in levels if record.continuous_no_initial_trapped_sphere
    )
    tiled = tuple(
        record.level.name for record in levels if record.complete_support_coverage
    )
    tiled_c = tuple(
        record.level.name for record in levels if _tiles_and_c_strictly_below_one(record)
    )
    payload = _contract_payload2(
        spec,
        frozen,
        levels,
        shared_coverage_monotone=monotone,
        compactness_upper_nonincreasing=nonincreasing,
        predecessor_nominal_contract_sha256=predecessor,
    )
    return SGBLProductBoxLadder2Record(
        slice=spec,
        ladder=frozen,
        levels=levels,
        predecessor_nominal_contract_sha256=predecessor,
        shared_coverage_monotone=monotone,
        compactness_upper_nonincreasing=nonincreasing,
        proved_level_names=proved,
        tiled_support_level_names=tiled,
        tiled_support_with_c_lt_1_level_names=tiled_c,
        any_declared_level_proves_continuous_no_initial_trap=bool(proved),
        any_declared_level_tiles_compact_support=bool(tiled),
        any_declared_level_tiles_and_c_strictly_below_one=bool(tiled_c),
        contract_payload=MappingProxyType(payload),
        contract_sha256=sgbl_trap_refinement_contract_sha256(payload),
    )


def sgbl_trap_refinement2_health_gate(
    record: SGBLProductBoxLadder2Record,
) -> dict[str, Any]:
    """Aggregate flags stay false even if a declared successor level locally passes."""

    if type(record) is not SGBLProductBoxLadder2Record:
        raise TypeError("record must be SGBLProductBoxLadder2Record")
    return {
        "any_declared_level_proves_continuous_no_initial_trap": (
            record.any_declared_level_proves_continuous_no_initial_trap
        ),
        "any_declared_level_tiles_compact_support": (
            record.any_declared_level_tiles_compact_support
        ),
        "any_declared_level_tiles_and_c_strictly_below_one": (
            record.any_declared_level_tiles_and_c_strictly_below_one
        ),
        "proved_level_names": record.proved_level_names,
        "evaluated_every_declared_level": record.evaluated_every_declared_level,
        "shared_coverage_monotone": record.shared_coverage_monotone,
        "predecessor_nominal_contract_sha256": record.predecessor_nominal_contract_sha256,
        "sampled_nodes_are_not_the_certificate": True,
        "product_box_is_not_the_certificate_chart": True,
        "SGBL_branch_owned_and_healthy": False,
        "FRZ1": False,
        "PREF1": False,
        "holdout_authorized": False,
        "execution_authorized": False,
        "copied_gr0_or_fgcqr_health_evidence": False,
    }


__all__ = [
    "DECLARED_LADDER2_BIT_CAP",
    "DECLARED_LADDER2_CELL_CAPS",
    "DECLARED_LADDER2_DEPTHS",
    "DECLARED_LADDER2_RHS_CAPS",
    "DECLARED_PRODUCT_BOX_RESOURCE_LADDER2",
    "INSTRUMENT_ID",
    "PREDECESSOR_INSTRUMENT_ID",
    "PREDECESSOR_NOMINAL_CONTRACT_SHA256",
    "SGBLProductBoxLadder2Record",
    "SGBLTrapRefinement2Stop",
    "sgbl_evaluate_product_box_resource_ladder2",
    "sgbl_trap_refinement2_health_gate",
]
