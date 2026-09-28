"""Decision gate for the joint NSC metric/state history boundary problem."""
from __future__ import annotations

import numpy as np


def nonuniqueness_witness(control_arrays):
    required = {"U_short", "U_long", "C_short", "C_long"}
    if not required <= set(control_arrays):
        raise ValueError("two conditional U/C controls required")
    u_short = np.asarray(control_arrays["U_short"])
    u_long = np.asarray(control_arrays["U_long"])
    c_short = np.asarray(control_arrays["C_short"])
    c_long = np.asarray(control_arrays["C_long"])
    if not (u_short.shape == u_long.shape == c_short.shape == c_long.shape):
        raise ValueError("control block shapes differ")
    return {
        "maximum_absolute_U_difference": float(np.max(abs(u_short-u_long))),
        "maximum_absolute_C_difference": float(np.max(abs(c_short-c_long))),
        "frobenius_U_difference": float(np.linalg.norm((u_short-u_long).ravel())),
        "frobenius_C_difference": float(np.linalg.norm((c_short-c_long).ravel())),
        "distinct_admissible_control_count_lower_bound": 2,
    }
