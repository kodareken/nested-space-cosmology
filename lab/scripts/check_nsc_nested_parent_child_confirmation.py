#!/usr/bin/env python3
"""Independent resident-data check with frozen bases and matched coordinate clocks."""
import json
import numpy as np
import derive_nsc_nested_parent_child_replay_basis as basis


def matched_times(cases):
    times = [np.array([row["time"] for row in case["rows"]]) for case in cases]
    if any(t.shape != times[0].shape or not np.allclose(t, times[0], rtol=0, atol=1e-12) for t in times[1:]):
        raise ValueError("cross-band confirmation uses mismatched coordinate clocks")


def physical_curve(case, observable):
    rows = case["rows"]
    if observable == "child_probability":
        values = np.array([row["windows"]["child"]["probability"] for row in rows])
    elif observable == "parent_annulus_probability":
        values = np.array([row["windows"]["left_annulus"]["probability"] +
                           row["windows"]["right_annulus"]["probability"] for row in rows])
    else:
        values = np.array([row["metrics"][observable] for row in rows])
    return values - values[0]


def verify():
    basis_report = basis.verify()
    witness, confirmation = basis.records()
    pairs = [witness["results"][f"nf256_{c}_dt0.0005"] for c in ("baseline", "parent_only")]
    pairs += [confirmation["results"][f"nf512_{c}_dt0.0005"] for c in ("baseline", "parent_only")]
    matched_times(pairs)
    fine = physical_curve(pairs[3], "child_probability") - physical_curve(pairs[2], "child_probability")
    coarse = physical_curve(pairs[1], "child_probability") - physical_curve(pairs[0], "child_probability")
    effect = float(np.max(abs(fine)))
    gap = float(np.max(abs(fine - coarse)))
    if effect <= 0 or gap >= .01 * effect:
        raise ValueError("parent-to-child physical state is not space resolved")
    for name in ("child_to_parent_state", "parent_to_child_geometry", "child_to_parent_annulus_geometry"):
        row = witness["two_way_effects"][name]
        if max(row["time_indicator"], row["space_indicator"]) >= .01 * row["effect_max"]:
            raise ValueError("inherited primary physical effect is unresolved: " + name)
    with np.load(basis.INPUTS[3], allow_pickle=False) as stored:
        child = confirmation["local_response"]["child"]
        array = lambda value: stored[value["array"]]
        reduced, full, refined = (array(child[name]) for name in
                                  ("occupation", "occupation_full", "occupation_full_refined"))
        local_effect = np.max(abs(full - full[0]))
        if local_effect <= 0 or max(np.max(abs(reduced - full)), np.max(abs(refined - full))) >= .01 * local_effect:
            raise ValueError("same-source local reduction is unresolved")
        control_fractions = {}
        for name in child["controls"]:
            control, reference, ref_refined = (array(child[name + suffix]) for suffix in ("", "_full", "_full_refined"))
            signal = reference - full
            magnitude = np.max(abs(signal))
            error = np.max(abs((control - reduced) - signal))
            temporal = np.max(abs((ref_refined - refined) - signal))
            if magnitude > 128 * np.finfo(float).eps and max(error, temporal) >= .01 * magnitude:
                raise ValueError("local control's own effect is unresolved: " + name)
            control_fractions[name] = float(error / magnitude) if magnitude else None
    return {"whole_finite_goal": "PASS", "basis": basis_report,
            "parent_child_state_effect": effect, "space_fraction": gap / effect,
            "local_reduction_fraction": float(np.max(abs(reduced - full)) / local_effect),
            "local_control_fractions": control_fractions, "trajectory_reexecuted": False}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
