"""Conditional between-node and floating assembly bounds for the local gate."""
from fractions import Fraction

import numpy as np

VERIFICATION_NODE_COUNT = 257
CELL_COUNT = 256


def _fraction(value, name):
    if isinstance(value, bool) or value is None:
        raise ValueError("explicit numerical " + name + " required")
    value = Fraction(str(value)) if isinstance(value, float) else Fraction(value)
    if value < 0:
        raise ValueError(name + " must be nonnegative")
    return value


def lipschitz_continuous_max(nodes, value_uppers, derivative_uppers):
    """Enclose sup |f| from nodal |f| uppers and one derivative bound per cell.

    On [z_i,z_{i+1}], both endpoint cones hold. The simpler uniform
    intersection ``min(a,b)+D*h`` is used. A valid D necessarily encloses the
    endpoint difference, so this is never below either endpoint.
    """
    nodes = tuple(Fraction(str(value)) if isinstance(value, float) else Fraction(value)
                  for value in nodes)
    values = tuple(_fraction(value, "nodal value upper") for value in value_uppers)
    derivatives = tuple(_fraction(value, "cell derivative upper")
                        for value in derivative_uppers)
    if len(nodes) != len(values) or len(nodes) < 2 or len(derivatives) != len(nodes) - 1:
        raise ValueError("matching nodes, nodal values and cell derivatives required")
    if any(right <= left for left, right in zip(nodes, nodes[1:])):
        raise ValueError("strictly increasing nodes required")
    cells = []
    for index, (left, right, derivative) in enumerate(
            zip(nodes, nodes[1:], derivatives)):
        width = right - left
        endpoint_gap = abs(values[index + 1] - values[index])
        if endpoint_gap > derivative * width:
            raise ValueError("derivative bound does not contain endpoint movement")
        upper = min(values[index], values[index + 1]) + derivative * width
        cells.append(upper)
    nodal = max(values)
    continuous = max(cells)
    return {
        "nodal_maximum_upper": nodal,
        "continuous_maximum_upper": continuous,
        "between_node_remainder_upper": max(Fraction(), continuous - nodal),
        "cell_maximum_uppers": tuple(cells),
        "derivative_bound_required": True,
        "dense_sampling_used_as_bound": False,
    }


def componentwise_continuous_max(nodes, values, derivative_uppers):
    values = np.asarray(values, object)
    derivative_uppers = np.asarray(derivative_uppers, object)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError("nodal values must have shape (nodes,2)")
    if derivative_uppers.shape != (len(values) - 1, 2):
        raise ValueError("derivative uppers must have shape (nodes-1,2)")
    return {
        name: lipschitz_continuous_max(nodes, values[:, index],
                                       derivative_uppers[:, index])
        for index, name in enumerate(("N", "beta"))
    }


def gamma_n(operations, unit_roundoff=Fraction(1, 2**53)):
    if isinstance(operations, bool) or int(operations) != operations or operations < 0:
        raise ValueError("nonnegative integer operation count required")
    operations = int(operations)
    unit = _fraction(unit_roundoff, "unit roundoff")
    if operations * unit >= 1:
        raise ValueError("roundoff model requires n*u < 1")
    return operations * unit / (1 - operations * unit)


def rounded_sum_error(term_absolute_uppers, *, unit_roundoff=Fraction(1, 2**53)):
    terms = tuple(_fraction(value, "term absolute upper")
                  for value in term_absolute_uppers)
    if not terms:
        return Fraction()
    return gamma_n(len(terms) - 1, unit_roundoff) * sum(terms)


def verification_componentwise_continuous_max(nodes, values, derivative_uppers):
    """Componentwise enclosure on the owned 257-node / 256-cell verification grid.

    This is the same Lipschitz intersection as ``componentwise_continuous_max``,
    restricted to the certification grid. It does not accept a shorter sample.
    """
    nodes = tuple(nodes)
    values = np.asarray(values, object)
    derivative_uppers = np.asarray(derivative_uppers, object)
    if len(nodes) != VERIFICATION_NODE_COUNT:
        raise ValueError("exactly 257 verification nodes required")
    if values.shape != (VERIFICATION_NODE_COUNT, 2):
        raise ValueError("nodal values must have shape (257, 2)")
    if derivative_uppers.shape != (CELL_COUNT, 2):
        raise ValueError("derivative uppers must have shape (256, 2)")
    return componentwise_continuous_max(nodes, values, derivative_uppers)


def residual_assembly_roundoff(component_term_uppers, *, extra_operations=0,
                                unit_roundoff=Fraction(1, 2**53)):
    """Roundoff of final baseline+geometry+matter+edge additions.

    Matrix products and source contractions need their own directed owner;
    this primitive covers only assembly of their already-enclosed outputs.
    """
    if set(component_term_uppers) != {"N", "beta"}:
        raise ValueError("N and beta term lists required")
    result = {}
    for name, terms in component_term_uppers.items():
        terms = tuple(_fraction(value, name + " term") for value in terms)
        operations = max(0, len(terms) - 1) + int(extra_operations)
        result[name] = gamma_n(operations, unit_roundoff) * sum(terms)
    return {
        "N": result["N"],
        "beta": result["beta"],
        "scope": "final residual additions only; field/contraction arithmetic separate",
    }
