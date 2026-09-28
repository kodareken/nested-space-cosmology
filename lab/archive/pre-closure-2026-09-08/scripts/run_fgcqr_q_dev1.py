#!/usr/bin/env python3
"""Measure complete metric-null Raychaudhuri Q on the activated FGC-QR trajectory."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from recursive_horizons.fgc.evolution.affine_measurement import (  # noqa: E402
    ADMHistory,
    analyze_affine_null,
)
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    radial_null_observables,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    EvolutionState,
    PRIMARY_METHOD,
)


ARTIFACT_ID = "FGC-1-FGCQR-Q-A6042-DEV1"
SCHEMA = "FGC-1-FGCQR-Q-A6042-DEV1-v1"
SOURCES = {
    1: ROOT / "runs/fgc-2-sf1/fgcqr-activation-a6042-dev1-history.npz",
    2: ROOT / "runs/fgc-2-sf1/fgcqr-activation-a6042-dev2-history.npz",
    3: ROOT / "runs/fgc-2-sf1/fgcqr-activation-a6042-dev3-history.npz",
}
OUTPUTS = {
    1: ROOT / "runs/fgc-2-sf1/fgcqr-q-a6042-dev1.json",
    2: ROOT / "runs/fgc-2-sf1/fgcqr-q-a6042-dev2.json",
    3: ROOT / "runs/fgc-2-sf1/fgcqr-q-a6042-dev3.json",
}


def _result_arrays(result):
    interior = result.interior
    return {
        "coordinate_times_hex": [float(value).hex() for value in result.coordinate_times],
        "affine_parameters_hex": [float(value).hex() for value in result.affine_parameters],
        "coordinate_radii_hex": [float(value).hex() for value in result.coordinate_radii],
        "theta_hex": [float(value).hex() for value in result.theta],
        "direct_dtheta_dlambda_hex": [
            float(value).hex() for value in result.direct_dtheta_dlambda
        ],
        "complete_Q_hex": [
            float(value).hex() for value in result.complete_raychaudhuri_rhs
        ],
        "ricci_null_hex": [float(value).hex() for value in result.ricci_null],
        "interior_sample_count": int(result.complete_raychaudhuri_rhs[interior].size),
        "minimum_complete_Q": float(np.min(result.complete_raychaudhuri_rhs[interior])),
        "maximum_complete_Q": float(np.max(result.complete_raychaudhuri_rhs[interior])),
        "positive_complete_Q_sample_count": int(
            np.sum(result.complete_raychaudhuri_rhs[interior] > 0.0)
        ),
        "maximum_route_disagreement": result.maximum_route_disagreement,
        "maximum_null_residual": result.maximum_null_residual,
        "maximum_affine_residual": result.maximum_affine_residual,
    }


def run(
    *,
    development: int,
    source_tag: str | None,
    launch_sample: int,
    affine_step: float,
) -> dict[str, object]:
    if source_tag is None:
        source_path = SOURCES[development]
        output_path = OUTPUTS[development]
        artifact_id = ARTIFACT_ID.replace("DEV1", f"DEV{development}")
    else:
        source_path = ROOT / f"runs/fgc-2-sf1/fgcqr-activation-a6042-{source_tag}-history.npz"
        output_path = ROOT / f"runs/fgc-2-sf1/fgcqr-q-a6042-{source_tag}.json"
        artifact_id = f"FGC-1-FGCQR-Q-A6042-{source_tag.upper()}"
    schema = f"{artifact_id}-v1"
    if output_path.exists():
        raise RuntimeError("FGC-QR Q DEV1 one-shot output already exists")
    with np.load(source_path, allow_pickle=False) as data:
        arrays = {name: np.array(data[name], copy=True) for name in data.files}
    raw_times = arrays["times"]
    uniform_times = np.linspace(float(raw_times[0]), float(raw_times[-1]), raw_times.size)

    def interpolate_time(values):
        upper = np.searchsorted(raw_times, uniform_times, side="left")
        upper = np.clip(upper, 1, raw_times.size - 1)
        lower = upper - 1
        denominator = raw_times[upper] - raw_times[lower]
        weight = (uniform_times - raw_times[lower]) / denominator
        shape = (weight.size,) + (1,) * (values.ndim - 1)
        return (1.0 - weight.reshape(shape)) * values[lower] + weight.reshape(shape) * values[upper]

    for name in ("u", "p", "q", "a", "p_r", "q_r"):
        arrays[name] = interpolate_time(arrays[name])
    arrays["times"] = uniform_times
    history = ADMHistory(
        times=arrays["times"],
        radii=arrays["radii"],
        u=arrays["u"],
        p=arrays["p"],
        q=arrays["q"],
        a=arrays["a"],
        p_r=arrays["p_r"],
        q_r=arrays["q_r"],
    )
    if not 0 <= launch_sample < history.times.size - 1:
        raise ValueError("launch sample must precede the final history slice")
    state = EvolutionState(
        history.u[launch_sample],
        history.p[launch_sample],
        history.q[launch_sample],
    )
    observables = radial_null_observables(state)
    measured = (history.radii > 0.0) & (history.radii <= 24.0)
    scores = np.minimum(-observables.theta_plus, -observables.theta_minus)
    trapped = measured & (scores > 0.0)
    if not np.any(trapped):
        raise RuntimeError("selected launch slice has no trapped cell")
    launch_index = int(np.argmax(np.where(trapped, scores, -np.inf)))
    launch_time = float(history.times[launch_sample])
    launch_radius = float(history.radii[launch_index])
    final_time = float(history.times[-1])
    results = {}
    for label, method in (("RK4", PRIMARY_METHOD), ("SSPRK3", COMPARATOR_METHOD)):
        results[label] = _result_arrays(
            analyze_affine_null(
                history,
                method=method,
                launch_time=launch_time,
                launch_radius=launch_radius,
                final_time=final_time,
                step_size=affine_step,
            )
        )
    positive_both = all(item["maximum_complete_Q"] > 0.0 for item in results.values())
    record = {
        "artifact_id": artifact_id,
        "schema": schema,
        "source_history": str(source_path.relative_to(ROOT)),
        "time_interpolation": {
            "method": "piecewise_linear_to_uniform_saved_slice_grid",
            "source_sample_count": int(raw_times.size),
            "uniform_spacing": float(uniform_times[1] - uniform_times[0]),
            "maximum_source_spacing_deviation": float(
                np.max(np.abs(np.diff(raw_times) - np.mean(np.diff(raw_times))))
            ),
        },
        "launch": {
            "sample_index": launch_sample,
            "time": launch_time,
            "radius": launch_radius,
            "grid_index": launch_index,
            "trapped_score": float(scores[launch_index]),
            "theta_plus": float(observables.theta_plus[launch_index]),
            "theta_minus": float(observables.theta_minus[launch_index]),
        },
        "final_time": final_time,
        "affine_step": affine_step,
        "results": results,
        "classification": (
            "activated_trapped_trajectory_contains_positive_complete_Q_samples"
            if positive_both
            else "activated_trapped_trajectory_complete_Q_nonpositive"
        ),
        "terminal": True,
        "nonclaims": {
            "positive_Q_finite_interval_with_error_margin": False,
            "two_evolution_method_three_resolution_confirmation": False,
            "metric_null_defocusing_demonstrated": False,
        },
    }
    output_path.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--development", type=int, choices=(1, 2, 3), default=1)
    parser.add_argument("--source-tag")
    parser.add_argument("--launch-sample", type=int, default=12)
    parser.add_argument("--affine-step", type=float, default=0.0005)
    args = parser.parse_args()
    record = run(
        development=args.development,
        source_tag=args.source_tag,
        launch_sample=args.launch_sample,
        affine_step=args.affine_step,
    )
    print(
        json.dumps(
            {
                "classification": record["classification"],
                "launch": record["launch"],
                "results": {
                    label: {
                        key: value
                        for key, value in result.items()
                        if key
                        in {
                            "minimum_complete_Q",
                            "maximum_complete_Q",
                            "positive_complete_Q_sample_count",
                            "maximum_route_disagreement",
                            "maximum_null_residual",
                            "maximum_affine_residual",
                        }
                    }
                    for label, result in record["results"].items()
                },
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
