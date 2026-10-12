"""Directed interval arithmetic for remaining local-gate residual assembly.

The remaining arithmetic scope is matrix products, contractions, signed-family
summation and the baseline+geometry+matter additions. N and beta stay
independent. Phase contraction arithmetic is already enclosed by
``phase_value`` and is not counted again. The exact-rational ``gamma_n``
term is applied last. Missing interval primitives or final arrays keep the
production component null; observed float drift is not a bound.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from types import MappingProxyType

from .nsc_ks_continuous_constraint_bound import gamma_n
from .nsc_ks_continuous_constraint_enclosure import (
    CONSTRAINTS,
    serialize_fraction,
    outward_float,
    _nonneg_fraction,
)


SCHEMA = "NSC-KS-ASSEMBLY-ARITHMETIC-v1"
UNIT_ROUNDOFF = Fraction(1, 2 ** 53)
DECLARED_POSITIVE_FAMILIES = 60
DECLARED_SIGNED_FAMILIES = 120
PHASE_VALUE_SCOPE = (
    "phase contraction arithmetic is already enclosed by phase_value; "
    "the remaining arithmetic component excludes that work")
REQUIRED_PRODUCTION_ARRAYS = (
    "interval/directed abs majorants of delta, C_src and reference columns "
    "at the assembly nodes for every contracted family",
    "interval/directed abs majorants of the N,beta vertex tensors",
    "family-summand N,beta abs enclosures for all 120 signed families",
    "interval/directed abs enclosures of baseline, geometry and matter; "
    "edge/phase terms are excluded from this component",
    "authenticated original weights, signs, degeneracies and coherences "
    "for all 120 signed families from the source inventory",
    "final assembled residual arrays only as a replay identity, never as a "
    "substitute for the directed operation inventory",
)


def signed_family_coverage(inventory=None):
    """Measured 60/120 coverage from an authenticated source inventory.

    Caller-supplied family arrays and flags are not coverage. This module
    does not authenticate RetainedSourceInventory bytes, so production
    coverage stays not covered.
    """
    _ = inventory
    return MappingProxyType({
        "declared_positive_families": DECLARED_POSITIVE_FAMILIES,
        "declared_signed_families": DECLARED_SIGNED_FAMILIES,
        "covered_positive_families": None,
        "covered_signed_families": None,
        "all_declared_source_paths_closed": False,
        "authenticated": False,
        "weights_bound": False,
        "signs_bound": False,
        "degeneracies_bound": False,
        "coherences_bound": False,
    })


def _inventory_closed(coverage):
    return (
        coverage.get("authenticated") is True
        and coverage.get("covered_positive_families") == DECLARED_POSITIVE_FAMILIES
        and coverage.get("covered_signed_families") == DECLARED_SIGNED_FAMILIES
        and coverage.get("all_declared_source_paths_closed") is True
        and coverage.get("weights_bound") is True
        and coverage.get("signs_bound") is True
        and coverage.get("degeneracies_bound") is True
        and coverage.get("coherences_bound") is True
    )


def _matrix(values, name):
    if values is None:
        raise ValueError("explicit " + name + " required")
    if isinstance(values, (bool, str)):
        raise ValueError(name + " must be a rectangular nonnegative matrix")
    rows = tuple(tuple(row) for row in values)
    if not rows or any(len(row) != len(rows[0]) for row in rows):
        raise ValueError(name + " must be a rectangular matrix")
    return tuple(tuple(_nonneg_fraction(item, name) for item in row) for row in rows)


def _pair(value, name):
    if value is None:
        return None
    if isinstance(value, bool) or len(value) != 2:
        raise ValueError(name + " must contain N and beta")
    return (_nonneg_fraction(value[0], name + " N"),
            _nonneg_fraction(value[1], name + " beta"))


@dataclass(frozen=True, slots=True)
class Interval:
    """Closed rational interval. Endpoints are exact; no binary64 widening."""

    lo: Fraction
    hi: Fraction

    def __post_init__(self):
        lo = self.lo if isinstance(self.lo, Fraction) else Fraction(self.lo)
        hi = self.hi if isinstance(self.hi, Fraction) else Fraction(self.hi)
        if lo > hi:
            raise ValueError("interval endpoints must be ordered")
        object.__setattr__(self, "lo", lo)
        object.__setattr__(self, "hi", hi)

    @property
    def abs_upper(self):
        return max(abs(self.lo), abs(self.hi))

    def __add__(self, other):
        other = other if isinstance(other, Interval) else Interval(other, other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    def __mul__(self, other):
        other = other if isinstance(other, Interval) else Interval(other, other)
        products = (self.lo * other.lo, self.lo * other.hi,
                    self.hi * other.lo, self.hi * other.hi)
        return Interval(min(products), max(products))

    def as_record(self):
        return {
            "lo": {"numerator": str(self.lo.numerator),
                   "denominator": str(self.lo.denominator)},
            "hi": {"numerator": str(self.hi.numerator),
                   "denominator": str(self.hi.denominator)},
            "abs_upper": serialize_fraction(self.abs_upper),
        }


@dataclass(frozen=True, slots=True)
class OperationInventory:
    matrix_products: int = 0
    contractions: int = 0
    family_summations: int = 0
    residual_additions: int = 0
    phase_arithmetic_operations: int = 0

    def __post_init__(self):
        for name in self.__dataclass_fields__:
            value = getattr(self, name)
            if isinstance(value, bool) or int(value) != value or value < 0:
                raise ValueError("nonnegative integer " + name + " required")
            object.__setattr__(self, name, int(value))
        if self.phase_arithmetic_operations != 0:
            raise ValueError(
                "phase arithmetic is already enclosed by phase_value and "
                "must not enter the remaining arithmetic inventory")

    @property
    def total(self):
        return (self.matrix_products + self.contractions
                + self.family_summations + self.residual_additions)

    def added(self, **counts):
        data = {name: getattr(self, name) for name in self.__dataclass_fields__}
        for name, value in counts.items():
            data[name] = data[name] + int(value)
        return OperationInventory(**data)

    def as_record(self):
        return {
            "matrix_products": self.matrix_products,
            "contractions": self.contractions,
            "family_summations": self.family_summations,
            "residual_additions": self.residual_additions,
            "phase_arithmetic_operations": self.phase_arithmetic_operations,
            "total": self.total,
            "gamma_n_applied_last": True,
            "phase_value_excluded": True,
            "scope": PHASE_VALUE_SCOPE,
        }


def interval_matrix_product(left, right):
    """Exact interval product of real matrices. Counts every multiply-add."""
    a = tuple(tuple(item if isinstance(item, Interval) else Interval(item, item)
                    for item in row) for row in left)
    b = tuple(tuple(item if isinstance(item, Interval) else Interval(item, item)
                    for item in row) for row in right)
    if not a or not b or len(a[0]) != len(b):
        raise ValueError("compatible matrix shapes required")
    inner = len(b)
    result = []
    operations = 0
    for row in a:
        out = []
        for column in range(len(b[0])):
            total = row[0] * b[0][column]
            operations += 1
            for k in range(1, inner):
                total = total + row[k] * b[k][column]
                operations += 2
            out.append(total)
        result.append(tuple(out))
    return tuple(result), operations


def abs_matrix_product(left, right):
    """Directed |AB| majorant from entrywise nonnegative |A|, |B|."""
    a = _matrix(left, "left abs matrix")
    b = _matrix(right, "right abs matrix")
    if len(a[0]) != len(b):
        raise ValueError("compatible abs-matrix shapes required")
    inner = len(b)
    result = []
    operations = 0
    for row in a:
        out = []
        for column in range(len(b[0])):
            total = row[0] * b[0][column]
            operations += 1
            for k in range(1, inner):
                total += row[k] * b[k][column]
                operations += 2
            out.append(total)
        result.append(tuple(out))
    return tuple(result), operations


def frobenius_abs_contraction(left, right):
    """sum_ij |L_ij| |R_ij|. One multiply per entry; additions after the first."""
    a = _matrix(left, "left contraction abs")
    b = _matrix(right, "right contraction abs")
    if len(a) != len(b) or len(a[0]) != len(b[0]):
        raise ValueError("contraction factors must share shape")
    total = Fraction(0)
    operations = 0
    first = True
    for i, row in enumerate(a):
        for j, value in enumerate(row):
            term = value * b[i][j]
            operations += 1
            if first:
                total = term
                first = False
            else:
                total += term
                operations += 1
    return total, operations


def trace_abs_contraction(left, right):
    """sum_ij |L_ij| |R_ji|, the vertex-density pairing of the matter owner."""
    a = _matrix(left, "vertex abs")
    b = _matrix(right, "density abs")
    if len(a) != len(b[0]) or len(a[0]) != len(b):
        raise ValueError("trace contraction requires matching transpose shape")
    total = Fraction(0)
    operations = 0
    first = True
    for i, row in enumerate(a):
        for j, value in enumerate(row):
            term = value * b[j][i]
            operations += 1
            if first:
                total = term
                first = False
            else:
                total += term
                operations += 1
    return total, operations


def componentwise_vertex_contraction(vertex_n, vertex_beta, density):
    """Independent N,beta contractions against one density abs majorant."""
    n_value, n_ops = trace_abs_contraction(vertex_n, density)
    beta_value, beta_ops = trace_abs_contraction(vertex_beta, density)
    return (n_value, beta_value), n_ops + beta_ops


def family_sum_abs(terms):
    """Directed N,beta sum. Any missing family keeps the production sum null."""
    if terms is None:
        return None, OperationInventory(), "OPEN: family summands missing"
    sequence = tuple(terms)
    if not sequence:
        raise ValueError("at least one family summand required")
    pairs = []
    for index, term in enumerate(sequence):
        if term is None:
            return None, OperationInventory(), (
                "OPEN: family summand " + str(index) + " missing")
        pairs.append(_pair(term, "family summand"))
    total = list(pairs[0])
    operations = 0
    for term in pairs[1:]:
        total[0] += term[0]
        total[1] += term[1]
        operations += 2
    return (total[0], total[1]), OperationInventory(family_summations=operations), None


def add_residual_parts_abs(baseline, geometry, matter, edge=None, *,
                           include_edge=False):
    """Baseline+geometry+matter, and edge only when phase is not already enclosed."""
    if include_edge:
        raise ValueError(
            "edge/phase arithmetic is already enclosed by phase_value; "
            "do not add it in the remaining arithmetic component")
    if edge is not None:
        raise ValueError(
            "remaining arithmetic excludes the source-cutoff edge; "
            "phase_value already owns that contraction arithmetic")
    parts = []
    missing = []
    for name, value in (("baseline", baseline), ("geometry", geometry),
                        ("matter", matter)):
        if value is None:
            missing.append(name)
        else:
            parts.append(_pair(value, name))
    if missing:
        return None, OperationInventory(), (
            "OPEN: residual parts missing: " + ", ".join(missing))
    total = [parts[0][0], parts[0][1]]
    operations = 0
    for term in parts[1:]:
        total[0] += term[0]
        total[1] += term[1]
        operations += 2
    return (total[0], total[1]), OperationInventory(residual_additions=operations), None


def apply_gamma_n_last(exact_pair, inventory, term_abs_sum, *,
                       unit_roundoff=UNIT_ROUNDOFF):
    """Exact directed sum first, then one ``gamma_n`` term on the abs mass."""
    exact = _pair(exact_pair, "exact arithmetic pair")
    mass = _pair(term_abs_sum, "term absolute sum")
    if not isinstance(inventory, OperationInventory):
        raise TypeError("OperationInventory required")
    if inventory.total == 0:
        return exact, Fraction(0)
    growth = gamma_n(inventory.total, unit_roundoff)
    return (exact[0] + growth * mass[0], exact[1] + growth * mass[1]), growth


def stepwise_gamma_sum(terms, *, unit_roundoff=UNIT_ROUNDOFF):
    """Roundoff applied after every addition. Not the production rule."""
    values = tuple(_nonneg_fraction(item, "stepwise term") for item in terms)
    if not values:
        raise ValueError("stepwise terms required")
    total = values[0]
    for term in values[1:]:
        growth = gamma_n(1, unit_roundoff)
        total = total + term + growth * (abs(total) + abs(term))
    return total


def _open_arithmetic(missing, coverage, *, status=None):
    unique_missing = tuple(dict.fromkeys(missing))
    return {
        "schema": SCHEMA,
        "bound": None,
        "status": status or (
            "OPEN: remaining assembly arithmetic lacks directed interval "
            "inputs and authenticated 120-family inventory; "
            "production component stays null"),
        "missing_inputs": unique_missing,
        "operation_inventory": OperationInventory().as_record(),
        "gamma_n_applied_last": True,
        "phase_value_excluded": True,
        "observed_drift_used_as_bound": False,
        "certificate_use": False,
        "physical_EXISTENCE_certificate": False,
        "physical_NONEXISTENCE_certificate": False,
        "declared_positive_families": DECLARED_POSITIVE_FAMILIES,
        "declared_signed_families": DECLARED_SIGNED_FAMILIES,
        "covered_positive_families": coverage["covered_positive_families"],
        "covered_signed_families": coverage["covered_signed_families"],
        "source_inventory_authenticated": coverage["authenticated"] is True,
        "all_declared_source_paths_closed": coverage[
            "all_declared_source_paths_closed"] is True,
    }


def production_assembly_arithmetic(*,
                                   delta_abs=None, source_abs=None, reference_abs=None,
                                   vertex_n_abs=None, vertex_beta_abs=None,
                                   family_terms=None, baseline=None, geometry=None,
                                   matter=None, edge=None, observed_drift=None,
                                   include_edge=False, family_inventory=None):
    """Return a finite bound only from complete directed inputs.

    ``observed_drift`` is accepted only as a rejected diagnostic. It never
    becomes the production arithmetic component. Family coverage is not
    inferred from supplied arrays or caller-written flags.
    """
    if observed_drift is not None:
        raise ValueError(
            "observed float drift is not a directed arithmetic bound")
    coverage = signed_family_coverage(family_inventory)
    missing = []
    if any(item is None for item in (delta_abs, source_abs, reference_abs,
                                     vertex_n_abs, vertex_beta_abs)):
        missing.extend(REQUIRED_PRODUCTION_ARRAYS[:2])
    terms = None if family_terms is None else tuple(family_terms)
    if terms is None or len(terms) != DECLARED_SIGNED_FAMILIES:
        missing.append(REQUIRED_PRODUCTION_ARRAYS[2])
    if any(item is None for item in (baseline, geometry, matter)):
        missing.append(REQUIRED_PRODUCTION_ARRAYS[3])
    if not _inventory_closed(coverage):
        missing.append(REQUIRED_PRODUCTION_ARRAYS[4])
    if include_edge or edge is not None:
        raise ValueError(
            "remaining arithmetic excludes edge/phase work already in phase_value")
    if missing:
        return _open_arithmetic(missing, coverage)
    product_dc, ops_dc = abs_matrix_product(delta_abs, source_abs)
    product_dcr, ops_dcr = abs_matrix_product(product_dc, _transpose(reference_abs))
    (n_term, beta_term), ops_contract = componentwise_vertex_contraction(
        vertex_n_abs, vertex_beta_abs, product_dcr)
    family_total, family_inv, family_open = family_sum_abs(terms)
    if family_open:
        return _open_arithmetic(
            (REQUIRED_PRODUCTION_ARRAYS[2],), coverage, status=family_open)
    residual, residual_inv, residual_open = add_residual_parts_abs(
        baseline, geometry, matter)
    if residual_open:
        return _open_arithmetic(
            (REQUIRED_PRODUCTION_ARRAYS[3],), coverage, status=residual_open)
    inventory = OperationInventory(
        matrix_products=ops_dc + ops_dcr,
        contractions=ops_contract,
        family_summations=family_inv.family_summations,
        residual_additions=residual_inv.residual_additions,
    )
    exact = (n_term + family_total[0] + residual[0],
             beta_term + family_total[1] + residual[1])
    mass = (n_term + family_total[0] + residual[0],
            beta_term + family_total[1] + residual[1])
    rounded, growth = apply_gamma_n_last(exact, inventory, mass)
    return {
        "schema": SCHEMA,
        "bound": [outward_float(rounded[0]), outward_float(rounded[1])],
        "exact_upper_evidence": {
            "N": serialize_fraction(rounded[0]),
            "beta": serialize_fraction(rounded[1]),
        },
        "exact_before_rounding": {
            "N": serialize_fraction(exact[0]),
            "beta": serialize_fraction(exact[1]),
        },
        "gamma_n": serialize_fraction(growth),
        "status": "PASS: remaining assembly arithmetic with gamma_n last",
        "missing_inputs": (),
        "operation_inventory": inventory.as_record(),
        "gamma_n_applied_last": True,
        "phase_value_excluded": True,
        "observed_drift_used_as_bound": False,
        "certificate_use": True,
        "physical_EXISTENCE_certificate": False,
        "physical_NONEXISTENCE_certificate": False,
        "contraction_abs": {
            "N": serialize_fraction(n_term),
            "beta": serialize_fraction(beta_term),
        },
        "family_sum_abs": {
            "N": serialize_fraction(family_total[0]),
            "beta": serialize_fraction(family_total[1]),
        },
        "residual_parts_abs": {
            "N": serialize_fraction(residual[0]),
            "beta": serialize_fraction(residual[1]),
        },
        "declared_positive_families": DECLARED_POSITIVE_FAMILIES,
        "declared_signed_families": DECLARED_SIGNED_FAMILIES,
        "covered_positive_families": coverage["covered_positive_families"],
        "covered_signed_families": coverage["covered_signed_families"],
        "source_inventory_authenticated": coverage["authenticated"] is True,
        "all_declared_source_paths_closed": coverage[
            "all_declared_source_paths_closed"] is True,
    }


def _transpose(matrix):
    rows = _matrix(matrix, "transpose abs")
    return tuple(tuple(row[index] for row in rows) for index in range(len(rows[0])))


def diagnostic_control_arithmetic():
    """Exact rational control of the remaining operations. certificate_use is false."""
    left = ((Fraction(1), Fraction(2)), (Fraction(0), Fraction(1)))
    right = ((Fraction(3), Fraction(0)), (Fraction(1), Fraction(4)))
    product, product_ops = abs_matrix_product(left, right)
    contraction, contraction_ops = frobenius_abs_contraction(product, product)
    family, family_inv, _open = family_sum_abs((
        (Fraction(1, 8), Fraction(1, 16)),
        (Fraction(1, 8), Fraction(1, 16)),
        (Fraction(1, 4), Fraction(1, 8)),
    ))
    residual, residual_inv, _residual_open = add_residual_parts_abs(
        (Fraction(1, 2), Fraction(1, 3)),
        (Fraction(1, 4), Fraction(1, 5)),
        family)
    inventory = OperationInventory(
        matrix_products=product_ops,
        contractions=contraction_ops,
        family_summations=family_inv.family_summations,
        residual_additions=residual_inv.residual_additions,
    )
    exact = (contraction + residual[0], contraction + residual[1])
    rounded, growth = apply_gamma_n_last(exact, inventory, exact)
    stepwise = stepwise_gamma_sum(
        (contraction, residual[0]), unit_roundoff=Fraction(1, 1000))
    return MappingProxyType({
        "schema": SCHEMA,
        "scope": "synthetic exact-rational remaining-assembly control",
        "certificate_use": False,
        "status": "DIAGNOSTIC: exact rational matrix/contraction/family/residual control",
        "product": [[serialize_fraction(item) for item in row] for row in product],
        "contraction": serialize_fraction(contraction),
        "family_sum": {
            "N": serialize_fraction(family[0]),
            "beta": serialize_fraction(family[1]),
        },
        "residual_parts": {
            "N": serialize_fraction(residual[0]),
            "beta": serialize_fraction(residual[1]),
        },
        "operation_inventory": inventory.as_record(),
        "gamma_n": serialize_fraction(growth),
        "rounded": {
            "N": serialize_fraction(rounded[0]),
            "beta": serialize_fraction(rounded[1]),
        },
        "stepwise_gamma_is_not_production": True,
        "stepwise_control_unit_1e3_N": serialize_fraction(stepwise),
        "gamma_n_applied_last": True,
        "phase_value_excluded": True,
        "observed_drift_used_as_bound": False,
        "physical_EXISTENCE_certificate": False,
        "physical_NONEXISTENCE_certificate": False,
        "constraints": list(CONSTRAINTS),
    })


def validate_assembly_arithmetic_record(record):
    if not isinstance(record, dict) or record.get("schema") != SCHEMA:
        raise ValueError("NSC-KS-ASSEMBLY-ARITHMETIC-v1 required")
    if record.get("phase_value_excluded") is not True:
        raise ValueError("remaining arithmetic must exclude phase_value work")
    if record.get("gamma_n_applied_last") is not True:
        raise ValueError("gamma_n must be applied last")
    if record.get("observed_drift_used_as_bound"):
        raise ValueError("observed drift cannot be the arithmetic bound")
    inventory = record.get("operation_inventory")
    if not isinstance(inventory, dict) or inventory.get("phase_arithmetic_operations") not in (0, None):
        raise ValueError("phase arithmetic operations must stay zero")
    if record.get("physical_EXISTENCE_certificate"):
        raise ValueError("OPEN arithmetic may not be reported as EXISTENCE")
    if record.get("physical_NONEXISTENCE_certificate"):
        raise ValueError("OPEN arithmetic may not be reported as NON_EXISTENCE")
    token = str(record.get("status", "")).split(":", 1)[0]
    if token in ("PASS", "EXISTENCE") and record.get("certificate_use") is not True:
        raise ValueError("PASS arithmetic requires certificate_use")
    if token == "DIAGNOSTIC":
        if record.get("certificate_use") is not False:
            raise ValueError("diagnostic arithmetic cannot be used as a certificate")
        return record
    if record.get("source_inventory_authenticated") is not True:
        if record.get("covered_signed_families") not in (None,):
            raise ValueError(
                "unauthenticated arithmetic cannot report signed-family coverage")
        if record.get("covered_positive_families") not in (None,):
            raise ValueError(
                "unauthenticated arithmetic cannot report positive-family coverage")
        if token in ("PASS", "EXISTENCE"):
            raise ValueError(
                "PASS arithmetic requires authenticated 120-family inventory")
    if record.get("bound") is None:
        if record.get("certificate_use") is True:
            raise ValueError("null arithmetic bound cannot claim certificate_use")
        if not record.get("missing_inputs"):
            raise ValueError("null arithmetic bound must list missing inputs")
        if token != "OPEN":
            raise ValueError("null arithmetic bound remains OPEN")
    elif record.get("certificate_use") is not True:
        if token == "PASS":
            raise ValueError("PASS arithmetic requires certificate_use")
    if token == "PASS":
        if record.get("covered_signed_families") != DECLARED_SIGNED_FAMILIES:
            raise ValueError("PASS arithmetic requires 120 signed families")
        if record.get("source_inventory_authenticated") is not True:
            raise ValueError(
                "PASS arithmetic requires authenticated 120-family inventory")
    return record
