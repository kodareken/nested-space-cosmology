"""One shared exterior departure derivative for prefix jets and ALE rates.

For Z_x = a(eta,x) Z_eta + ..., a >= 0, characteristics travel towards
decreasing eta.  The physical wall at eta=1 supplies inflow.  The forward
quadratic stencil closes with a centred penultimate row and a backward
wall row.  No source or wall value is changed and no required wall jet is
substituted.  This polynomial consistency statement is not a stability
theorem for the coupled Einstein/scalar discretization.
"""
import math
from numbers import Real
import time

import numpy as np


MAX_CELLS = 512
MAX_CHANNELS = 3


def _array(value, name):
    # Reject booleans before list coercion can turn a mixed bool/float list
    # into an apparently real array.  Existing numeric arrays need no copy.
    def contains_bool(item):
        if isinstance(item, (bool, np.bool_)):
            return True
        if isinstance(item, (list, tuple)):
            return any(contains_bool(part) for part in item)
        return False

    if contains_bool(value):
        raise ValueError(f"{name} must contain finite real numbers, not booleans")
    try:
        raw = np.asarray(value)
        if raw.dtype.kind not in "fiu":
            raise ValueError(f"{name} must contain finite real numbers")
        result = raw.astype(float, copy=False)
    except (TypeError, OverflowError) as error:
        raise ValueError(f"{name} must contain finite real numbers") from error
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain finite real numbers")
    return result


def _inputs(eta, values, budget, minimum_cells):
    eta, values = _array(eta, "eta"), _array(values, "values")
    if (eta.ndim != 1 or not minimum_cells+1 <= eta.size <= MAX_CELLS+1
            or eta[0] != 0. or eta[-1] != 1.):
        raise ValueError(f"eta must have {minimum_cells}..{MAX_CELLS} cells on [0,1]")
    h = float(eta[1]-eta[0])
    if (h <= 0 or np.any(np.diff(eta) <= 0)
            or not np.allclose(np.diff(eta), h, rtol=5e-13, atol=0)):
        raise ValueError("eta must be a uniform increasing closed grid")
    if (values.ndim not in (1, 2) or values.shape[-1] != eta.size
            or (values.ndim == 2 and not 1 <= values.shape[0] <= MAX_CHANNELS)):
        raise ValueError("values must have shape (nodes,) or (1..3,nodes)")
    if not callable(getattr(budget, "check", None)):
        raise ValueError("a live shared cumulative budget is required")
    # Use the caller's allowance; this operator does not restart its clock or
    # replace a physical caller's larger, separately admitted budget.
    for name in ("maximum_evaluations", "evaluations"):
        item = getattr(budget, name, None)
        if (isinstance(item, (bool, np.bool_)) or not isinstance(item, (int, np.integer))
                or item < (1 if name == "maximum_evaluations" else 0)):
            raise ValueError("budget counts must be nonnegative integers with a positive limit")
    for name in ("cpu_seconds", "started"):
        item = getattr(budget, name, None)
        if (isinstance(item, (bool, np.bool_)) or not isinstance(item, Real)
                or not math.isfinite(float(item)) or item < 0
                or (name == "cpu_seconds" and item == 0)):
            raise ValueError("budget clock and positive CPU limit must be finite")
    if budget.started > time.process_time():
        raise ValueError("budget cannot start in future process time")
    budget.check(int(values.size))
    return values, h


def quadratic_departure_gradient(eta, values, *, budget):
    """Differentiate 4..512 uniform cells, along the final array axis.

    ``values`` is one channel (N+1,) or a stack of 1..3 channels (k,N+1).
    Each returned derivative value costs one shared-budget point.  The
    returned array is new; inputs remain unchanged.  The caller owns source
    anchoring and the wall value and must use this SAME result in its
    transverse prefix jets and moving-grid term.

    Rows through N-2 use (-3Z_i+4Z_{i+1}-Z_{i+2})/(2h); N-1 uses
    (Z_N-Z_{N-2})/(2h), and N uses (3Z_N-4Z_{N-1}+Z_{N-2})/(2h).
    All rows differentiate quadratics exactly, subject to roundoff.
    """
    values, h = _inputs(eta, values, budget, 4)
    result = np.empty_like(values)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            result[..., :-2] = (-3*values[..., :-2]+4*values[..., 1:-1]-values[..., 2:])/(2*h)
            result[..., -2] = (values[..., -1]-values[..., -3])/(2*h)
            result[..., -1] = (3*values[..., -1]-4*values[..., -2]+values[..., -3])/(2*h)
    except FloatingPointError as error:
        raise ValueError("nonfinite departure-gradient arithmetic") from error
    if not np.all(np.isfinite(result)):
        raise ValueError("nonfinite departure-gradient result")
    budget.check()
    return result


def first_order_departure_gradient(eta, values, *, budget):
    """The frozen forward/secant comparator, with the same shape and charge.

    This permits 2..512 cells to match the previous stage's grid contract.
    Its final wall row repeats the last interval secant.  It changes neither
    the input values nor any boundary data.
    """
    values, h = _inputs(eta, values, budget, 2)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            result = np.concatenate((np.diff(values, axis=-1),
                                     (values[..., -1]-values[..., -2])[..., None]), axis=-1)/h
    except FloatingPointError as error:
        raise ValueError("nonfinite departure-gradient arithmetic") from error
    if not np.all(np.isfinite(result)):
        raise ValueError("nonfinite departure-gradient result")
    budget.check()
    return result
