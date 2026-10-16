"""Ordered characteristic integrals on one positive-radius exterior slice.

The caller supplies X=(r,log(C),c), its independent Ubar jets, and the
longitudinal source jets. This stage reconstructs X_x; it neither chooses
spatial differences nor advances time, the source, the wall or the cap.
"""
from collections.abc import Mapping

import numpy as np
from scipy.integrate import cumulative_simpson, cumulative_trapezoid


MAX_NODES = 513
STATE_KEYS = ("r", "L", "c")
JET_KEYS = ("r_Ubar", "log_C_Ubar", "c_Ubar")
RATE_KEYS = ("r_x", "log_C_x", "c_x")


def _array(value, name, shape=None):
    try:
        raw = np.asarray(value)
        if raw.dtype.kind not in "iuf" or raw.size == 0:
            raise ValueError
        if isinstance(value, (list, tuple)) and any(
                isinstance(item, (bool, np.bool_)) for item in np.asarray(value, dtype=object).flat):
            raise ValueError
        result = np.array(raw, dtype=float, copy=True)
        if not np.all(np.isfinite(result)) or (shape is not None and result.shape != shape):
            raise ValueError
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(name + " must contain finite real values of the required shape") from error
    return result


def _mapping(value, keys, name, shape):
    if not isinstance(value, Mapping) or set(value) != set(keys):
        raise ValueError(name + " must contain exactly " + ", ".join(keys))
    return {key: _array(value[key], name + "." + key, shape) for key in keys}


def _check_tube(states, tube):
    if not isinstance(tube, Mapping) or set(tube) != set(STATE_KEYS):
        raise ValueError("tube must contain exactly r, L, c bounds")
    bounds = {key: _array(tube[key], "tube." + key, (2,)) for key in STATE_KEYS}
    for key, (low, high) in bounds.items():
        if not low < high:
            raise ValueError("tube bounds must be strictly increasing")
        if np.any(states[key] < low) or np.any(states[key] > high):
            raise ValueError("slice state exits the declared tube: " + key)
    if bounds["r"][0] <= 0:
        raise ValueError("tube radius floor must be positive")
    with np.errstate(over="ignore", under="ignore"):
        lapses = np.exp(bounds["L"])
    if np.any(lapses <= 0) or not np.all(np.isfinite(lapses)):
        raise ValueError("tube lapse bounds must be positive and representable")


def exterior_characteristic_slice(Ubar, states, first_jets, source_rates, *,
                                  alpha, quadrature_order=2, tube, budget):
    """Reconstruct longitudinal first jets by three deterministic prefixes.

    Ubar contains 3..513 uniform increasing nodes starting at the source 0.
    states has exactly r,L,c arrays; first_jets has r_Ubar,log_C_Ubar,c_Ubar
    arrays at the same nodes. The lapse Ubar jet is retained for interface
    compatibility but does not enter these exterior mixed equations.
    source_rates has three finite scalars: r_x,log_C_x,c_x at Ubar=0.

    The canonical normalized scalar f=Psi/sqrt(alpha) is used through its
    continuous jets f_a=u_a/sqrt(1+alpha*u**2/6), u=c/r. Ordered integrals
    solve (r*r_x)_Ubar=-C/4, (r*f_x)_Ubar=-r_x*f_Ubar, and the lapse mixed
    equation. alpha=0 keeps the normalized scalar present. Quadrature order
    2 is cumulative trapezoid; order 4 is cumulative Simpson, using adjacent
    slice nodes also for odd prefixes. No spatial field jets are inferred
    from samples. The caller must compare quadrature/spatial discretizations.

    tube supplies finite inclusive engineering bounds for r,L,c; it is not
    a proved invariant region. A supplied budget.check(count) owns the live
    cumulative ceiling. Four N-node passes are charged: normalized inputs
    and each of the three ordered prefixes. There is no internal budget
    reset, source evaluation, time integration, adaptive retry or fallback.

    Returned mixed_rhs values are pointwise mixed-equation evaluations at
    the reconstructed jets, not independently measured constraint defects.
    This stage imposes no pure null constraint or wall boundary condition.
    """
    Ubar = _array(Ubar, "Ubar")
    if (Ubar.ndim != 1 or not 3 <= Ubar.size <= MAX_NODES or Ubar[0] != 0
            or np.any(np.diff(Ubar) <= 0)):
        raise ValueError("Ubar must be 3..513 increasing nodes starting at zero")
    steps = np.diff(Ubar)
    if not np.allclose(steps, Ubar[-1]/(Ubar.size-1), rtol=5e-13, atol=0):
        raise ValueError("Ubar grid must be uniform")
    if type(quadrature_order) is not int or quadrature_order not in (2, 4):
        raise ValueError("quadrature_order must be 2 or 4")
    alpha = float(_array(alpha, "alpha", ()))
    if alpha < 0:
        raise ValueError("alpha must be nonnegative")
    states = _mapping(states, STATE_KEYS, "states", Ubar.shape)
    jets = _mapping(first_jets, JET_KEYS, "first_jets", Ubar.shape)
    source = _mapping(source_rates, RATE_KEYS, "source_rates", ())
    _check_tube(states, tube)
    if not callable(getattr(budget, "check", None)):
        raise ValueError("budget must expose a live cumulative check(count)")
    budget.check()
    count = int(Ubar.size)

    def prefix(values):
        budget.check(count)
        if quadrature_order == 2:
            return cumulative_trapezoid(values, x=Ubar, initial=0.)
        return cumulative_simpson(values, x=Ubar, initial=0.)

    try:
        with np.errstate(over="raise", divide="raise", invalid="raise", under="ignore"):
            budget.check(count)
            r, L, c = (states[key] for key in STATE_KEYS)
            rU, cU = jets["r_Ubar"], jets["c_Ubar"]
            C = np.exp(L)
            u = c/r
            D = 1+alpha*u*u/6
            root_D = np.sqrt(D)
            uU = (cU-u*rU)/r
            fU = uU/root_D
            source_fx = (source["c_x"]-u[0]*source["r_x"])/(r[0]*root_D[0])
            integral_C = prefix(C/4)
            rx = (r[0]*source["r_x"]-integral_C)/r
            integral_scalar = prefix(rx*fU)
            fx = (r[0]*source_fx-integral_scalar)/r
            cx = u*rx+r*root_D*fx
            # Preserve the supplied first jet exactly at the common source.
            rx[0], cx[0] = source["r_x"], source["c_x"]
            lapse_integrand = 2*rU*rx/(r*r)+C/(2*r*r)-alpha*fU*fx
            integral_lapse = prefix(lapse_integrand)
            Lx = source["log_C_x"]+integral_lapse
            radial_rhs = -(rU*rx+C/4)/r
            kernel = (rU*rx+C/4)/(r*r)-alpha*fU*fx/6
            scalar_rhs = -kernel*c
            arrays = {
                "Ubar": Ubar, "r_x": rx, "log_C_x": Lx, "c_x": cx,
                "D": D, "f_Ubar": fU, "f_x": fx,
                "integral_C_over_four": integral_C,
                "integral_r_x_f_Ubar": integral_scalar,
                "integral_lapse": integral_lapse,
                "r_Ubar_x": radial_rhs, "log_C_Ubar_x": lapse_integrand,
                "c_Ubar_x": scalar_rhs,
            }
    except (FloatingPointError, OverflowError) as error:
        raise ValueError("nonfinite exterior characteristic slice arithmetic") from error
    if not all(np.all(np.isfinite(value)) for value in arrays.values()):
        raise ValueError("nonfinite exterior characteristic slice output")
    budget.check()
    return {**arrays, "settings": {"nodes": count, "quadrature_order": quadrature_order,
                                  "charged_evaluations": 4*count},
            "source_rates": {key: float(value) for key, value in source.items()},
            "scope": {"source_integrations": 0, "time_steps": 0,
                      "pure_constraints_evaluated": False,
                      "quadrature_error_is_enclosed": False}}
