"""Remaining assembly arithmetic: exact rationals, rounding-last and nulls."""
from fractions import Fraction

import pytest

from recursive_horizons.nsc_ks_assembly_arithmetic import (
    DECLARED_SIGNED_FAMILIES,
    Interval,
    OperationInventory,
    abs_matrix_product,
    add_residual_parts_abs,
    apply_gamma_n_last,
    diagnostic_control_arithmetic,
    family_sum_abs,
    frobenius_abs_contraction,
    production_assembly_arithmetic,
    signed_family_coverage,
    stepwise_gamma_sum,
    trace_abs_contraction,
    validate_assembly_arithmetic_record,
)
from recursive_horizons.nsc_ks_continuous_constraint_bound import gamma_n


def test_interval_product_and_abs_matrix_are_exact():
    left = Interval(1, 2)
    right = Interval(3, 4)
    assert (left * right).lo == 3
    assert (left * right).hi == 8
    assert (left + right).lo == 4
    assert (left + right).hi == 6
    product, operations = abs_matrix_product(
        ((1, 2), (0, 1)), ((3, 0), (1, 4)))
    assert product == ((5, 8), (1, 4))
    assert operations == 1 + 2 + 1 + 2 + 1 + 2 + 1 + 2
    contraction, count = frobenius_abs_contraction(product, ((1, 0), (0, 1)))
    assert contraction == 5 + 4
    assert count == 4 + 3


def test_trace_contraction_and_family_sum_keep_N_beta_independent():
    value, operations = trace_abs_contraction(
        ((1, 2), (3, 4)), ((5, 6), (7, 8)))
    assert value == 1 * 5 + 2 * 7 + 3 * 6 + 4 * 8
    assert operations == 4 + 3
    total, inventory, opened = family_sum_abs((
        (Fraction(1, 2), Fraction(1, 3)),
        (Fraction(1, 4), Fraction(1, 5)),
        (Fraction(1, 8), Fraction(1, 7)),
    ))
    assert opened is None
    assert total[0] == Fraction(1, 2) + Fraction(1, 4) + Fraction(1, 8)
    assert total[1] == Fraction(1, 3) + Fraction(1, 5) + Fraction(1, 7)
    assert inventory.family_summations == 4
    assert inventory.matrix_products == 0
    null, _inventory, status = family_sum_abs((
        (1, 1), None, (1, 1)))
    assert null is None
    assert "missing" in status


def test_rounding_is_applied_last_not_per_summand():
    inventory = OperationInventory(residual_additions=2)
    exact = (Fraction(1), Fraction(2))
    rounded, growth = apply_gamma_n_last(exact, inventory, exact, unit_roundoff=Fraction(1, 1000))
    assert growth == gamma_n(2, Fraction(1, 1000))
    assert rounded[0] == 1 + growth * 1
    last = apply_gamma_n_last(
        (3, 0), OperationInventory(residual_additions=2), (3, 0),
        unit_roundoff=Fraction(1, 1000))[0][0]
    stepwise = stepwise_gamma_sum((1, 1, 1), unit_roundoff=Fraction(1, 1000))
    assert last == 3 + gamma_n(2, Fraction(1, 1000)) * 3
    assert stepwise != last
    with pytest.raises(ValueError, match="phase arithmetic"):
        OperationInventory(phase_arithmetic_operations=1)


def test_edge_and_observed_drift_cannot_fill_the_component():
    with pytest.raises(ValueError, match="phase_value"):
        add_residual_parts_abs((1, 1), (1, 1), (1, 1), edge=(1, 1))
    with pytest.raises(ValueError, match="drift"):
        production_assembly_arithmetic(observed_drift=(1e-16, 1e-16))
    opened = production_assembly_arithmetic()
    assert opened["bound"] is None
    assert opened["certificate_use"] is False
    assert opened["missing_inputs"]
    assert opened["gamma_n_applied_last"]
    assert opened["phase_value_excluded"]
    validate_assembly_arithmetic_record(opened)


def _synthetic_arrays(*, family_terms):
    return dict(
        delta_abs=((1, 0), (0, 1)),
        source_abs=((1, 0), (0, 1)),
        reference_abs=((1, 0), (0, 1)),
        vertex_n_abs=((1, 0), (0, 0)),
        vertex_beta_abs=((0, 0), (0, 1)),
        family_terms=family_terms,
        baseline=(Fraction(1, 2), Fraction(1, 4)),
        geometry=(Fraction(1, 8), Fraction(1, 16)),
        matter=(Fraction(1, 8), Fraction(1, 16)),
    )


def test_complete_directed_inputs_enclose_and_nulls_do_not_pass():
    one_family = production_assembly_arithmetic(
        **_synthetic_arrays(family_terms=((Fraction(1, 8), Fraction(1, 16)),)))
    assert one_family["status"].startswith("OPEN")
    assert one_family["certificate_use"] is False
    assert one_family["bound"] is None
    assert one_family["source_inventory_authenticated"] is False
    assert one_family["covered_signed_families"] is None
    assert one_family["declared_signed_families"] == DECLARED_SIGNED_FAMILIES
    synthetic_120 = production_assembly_arithmetic(
        **_synthetic_arrays(family_terms=tuple(
            (Fraction(1, 8), Fraction(1, 16))
            for _ in range(DECLARED_SIGNED_FAMILIES))))
    assert synthetic_120["status"].startswith("OPEN")
    assert synthetic_120["certificate_use"] is False
    assert synthetic_120["bound"] is None
    assert synthetic_120["source_inventory_authenticated"] is False
    coverage = signed_family_coverage({"signed_families": 120})
    assert coverage["authenticated"] is False
    assert coverage["covered_signed_families"] is None
    assert one_family["operation_inventory"]["phase_arithmetic_operations"] == 0
    assert one_family["gamma_n_applied_last"]
    forged = dict(one_family)
    forged["physical_EXISTENCE_certificate"] = True
    with pytest.raises(ValueError, match="EXISTENCE"):
        validate_assembly_arithmetic_record(forged)
    null = dict(production_assembly_arithmetic())
    null["status"] = "PASS: invented"
    null["certificate_use"] = True
    with pytest.raises(ValueError, match="null|authenticated"):
        validate_assembly_arithmetic_record(null)


def test_diagnostic_control_is_not_a_certificate():
    diagnostic = dict(diagnostic_control_arithmetic())
    assert diagnostic["certificate_use"] is False
    assert diagnostic["gamma_n_applied_last"]
    assert diagnostic["stepwise_gamma_is_not_production"]
    validate_assembly_arithmetic_record(diagnostic)
    diagnostic["status"] = "PASS: control"
    with pytest.raises(ValueError, match="certificate"):
        validate_assembly_arithmetic_record(diagnostic)
