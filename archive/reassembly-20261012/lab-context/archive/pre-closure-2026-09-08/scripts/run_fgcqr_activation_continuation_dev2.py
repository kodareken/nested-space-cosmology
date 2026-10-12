#!/usr/bin/env python3
"""Continue the activated FGC-QR trajectory until backreaction is measurable."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import multiprocessing
from pathlib import Path
from types import SimpleNamespace
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
    PRIMARY_METHOD,
    accept_step,
    propose_step,
)
from scripts.run_fgcqr_activation_dev1 import (  # noqa: E402
    ACTION,
    FGCQROperator,
    _projector,
    _sample,
    _sha,
    _verify_fast_batch_gb,
)


ARTIFACT_ID = "FGC-1-FGCQR-ACTIVATION-A6042-DEV2"
SCHEMA = "FGC-1-FGCQR-ACTIVATION-A6042-DEV2-v1"
def run(
    *, source_development: int, target_time: float, stop_factor: float, workers: int
) -> dict[str, object]:
    next_development = source_development + 1
    source_json = ROOT / f"runs/fgc-2-sf1/fgcqr-activation-a6042-dev{source_development}.json"
    source_history = ROOT / f"runs/fgc-2-sf1/fgcqr-activation-a6042-dev{source_development}-history.npz"
    output = ROOT / f"runs/fgc-2-sf1/fgcqr-activation-a6042-dev{next_development}.json"
    history_path = ROOT / f"runs/fgc-2-sf1/fgcqr-activation-a6042-dev{next_development}-history.npz"
    artifact_id = ARTIFACT_ID.replace("DEV2", f"DEV{next_development}")
    schema = SCHEMA.replace("DEV2", f"DEV{next_development}")
    if output.exists() or history_path.exists():
        raise RuntimeError(f"activation DEV{next_development} one-shot output already exists")
    source = json.loads(source_json.read_text())
    with np.load(source_history, allow_pickle=False) as data:
        arrays = {name: np.array(data[name], copy=True) for name in data.files}
    state = EvolutionState(arrays["u"][-1], arrays["p"][-1], arrays["q"][-1])
    reference_state = EvolutionState(arrays["u"][0], arrays["p"][0], arrays["q"][0])
    grid = SimpleNamespace(
        coordinates=arrays["radii"],
        point_count=arrays["radii"].size,
        spacing=float(arrays["radii"][1] - arrays["radii"][0]),
        minimum=float(arrays["radii"][0]),
        maximum=float(arrays["radii"][-1]),
    )
    # SBPFirstDerivative requires its concrete grid type; reconstruct it exactly.
    from recursive_horizons.fgc.evolution.numerical_engine import UniformRadialGrid

    grid = UniformRadialGrid(grid.minimum, grid.maximum, grid.point_count)
    initial = SimpleNamespace(state=reference_state)
    projector = _projector(initial)
    seed_peak = float(np.max(np.abs(arrays["u"][0, :, 4])))
    time = float(arrays["times"][-1])
    step_index = int(
        source["execution"].get(
            "accepted_steps_total", source["execution"].get("accepted_steps")
        )
    )
    serial = 5 * step_index
    retries = 0
    samples = list(source["samples"])
    history = {
        name: list(arrays[name])
        for name in ("times", "u", "p", "q", "a", "p_r", "q_r")
    }
    terminal_kind = "target_time_before_backreaction_target"
    error = None
    comparison_error = _verify_fast_batch_gb()
    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=multiprocessing.get_context("fork"),
    ) as pool:
        operator = FGCQROperator(
            grid,
            spatial_order=4,
            ko_dissipation=1.0 / 64.0,
            pool=pool,
            workers=workers,
        )
        _sample(time, state, operator, seed_peak)
        while time < target_time - 1.0e-14:
            speed = max(float(operator.last_diagnostics["coordinate_speed_upper"]), 1.0e-12)
            width = min(target_time - time, 0.125 * grid.spacing / speed)
            local_retries = 0
            while True:
                try:
                    proposal = propose_step(
                        method=PRIMARY_METHOD,
                        time=time,
                        step_size=width,
                        state=state,
                        rhs=operator,
                        projector=projector,
                    )
                    candidate = proposal.candidate_state
                    effective_planck = ACTION.planck_mass**2 + ACTION.beta * candidate.u[:, 4] ** 2
                    if (
                        np.any(candidate.u[:, 0] <= 0.0)
                        or np.any(candidate.u[:, 2] <= 0.0)
                        or np.any(effective_planck <= 0.0)
                    ):
                        raise ValueError("candidate crossed a positive metric or F boundary")
                    accepted = accept_step(
                        proposal,
                        previous_step_index=step_index,
                        previous_transaction_serial=serial,
                    )
                    break
                except (ArithmeticError, ValueError) as exc:
                    local_retries += 1
                    retries += 1
                    width *= 0.5
                    if local_retries > 20 or width < 2.0**-30:
                        error = {"type": type(exc).__name__, "message": str(exc)}
                        terminal_kind = "typed_numerical_or_health_stop"
                        break
            if error is not None:
                break
            state = accepted.state
            time = accepted.time
            step_index = accepted.step_index
            serial = accepted.transaction_serial
            item = _sample(time, state, operator, seed_peak)
            samples.append(item)
            history["times"].append(time)
            history["u"].append(state.u.copy())
            history["p"].append(state.p.copy())
            history["q"].append(state.q.copy())
            history["a"].append(operator.last_physical_acceleration.copy())
            history["p_r"].append(operator.last_p_r.copy())
            history["q_r"].append(operator.last_q_r.copy())
            if item["activation_factor"] >= stop_factor:
                terminal_kind = "backreaction_activation_target_reached"
                break
    np.savez_compressed(
        history_path,
        **{name: np.asarray(values) for name, values in history.items()},
        radii=arrays["radii"],
    )
    record = {
        "artifact_id": artifact_id,
        "schema": schema,
        "source_terminal": str(source_json.relative_to(ROOT)),
        "source_history": str(source_history.relative_to(ROOT)),
        "classification": terminal_kind,
        "action": source["action"],
        "execution": {
            "method": "RK4",
            "spatial_order": 4,
            "source_workers": workers,
            "source_final_time": float(arrays["times"][-1]),
            "target_time": target_time,
            "stop_activation_factor": stop_factor,
            "final_time": time,
            "accepted_steps_total": step_index,
            "continuation_retries": retries,
            "error": error,
            "fast_batch_GB_RED1_comparison_error": comparison_error,
        },
        "samples": samples,
        "history": {
            "path": str(history_path.relative_to(ROOT)),
            "sha256": _sha(history_path),
            "sample_count": len(history["times"]),
        },
        "terminal": True,
        "nonclaims": {
            "positive_complete_Q": False,
            "two_method_three_resolution_confirmation": False,
            "metric_null_defocusing_demonstrated": False,
        },
    }
    output.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-development", type=int, choices=(1, 2), default=1)
    parser.add_argument("--target-time", type=float, default=1.0)
    parser.add_argument("--stop-factor", type=float, default=128.0)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    result = run(
        source_development=args.source_development,
        target_time=args.target_time,
        stop_factor=args.stop_factor,
        workers=args.workers,
    )
    print(
        json.dumps(
            {
                "classification": result["classification"],
                "execution": result["execution"],
                "final_sample": result["samples"][-1],
                "history": result["history"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
