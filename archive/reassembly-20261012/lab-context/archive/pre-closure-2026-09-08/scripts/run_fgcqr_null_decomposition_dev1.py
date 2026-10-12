#!/usr/bin/env python3
"""Decompose null Ricci focusing into the complete FGC-QR action terms."""

from __future__ import annotations

import json
import argparse
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from recursive_horizons.fgc.evolution.floating_jet import (  # noqa: E402
    FloatJet2,
    scalar_primal,
)
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    radial_null_observables,
)
from recursive_horizons.fgc.evolution.numerical_engine import EvolutionState  # noqa: E402
from recursive_horizons.fgc.spherical_reduction import (  # noqa: E402
    _state_from_adm_pg_jets,
    residuals,
)


ARTIFACT_ID = "FGC-1-FGCQR-NULL-DECOMPOSITION-A6042-DEV1"
SCHEMA = "FGC-1-FGCQR-NULL-DECOMPOSITION-A6042-DEV1-v1"
SOURCE = ROOT / "runs/fgc-2-sf1/fgcqr-activation-a6042-dev3-history.npz"
OUTPUT = ROOT / "runs/fgc-2-sf1/fgcqr-null-decomposition-a6042-dev1.json"
ACTION = {
    "planck_mass": 2,
    "scalar_mass": 0.02,
    "quartic_coupling": 256,
    "beta": -128,
    "alpha_gb": 0,
    "eta": 16384,
}
FIELDS = ("alpha", "shift", "lambda", "areal_radius", "phi", "chi")


def _scalar(value) -> float:
    return float(scalar_primal(value))


def _contract(matrix, null) -> float:
    return float(
        sum(
            null[a] * null[b] * _scalar(matrix[a][b])
            for a in range(2)
            for b in range(2)
        )
    )


def _point(index, radii, u, p, q, a, p_r, q_r):
    jets = {
        name: FloatJet2(
            u[index, column],
            p[index, column],
            q[index, column],
            a[index, column],
            p_r[index, column],
            q_r[index, column],
        )
        for column, name in enumerate(FIELDS)
    }
    state = _state_from_adm_pg_jets(
        {"model_id": "FGC-QR", "action_parameters": ACTION}, jets
    )
    result = residuals(state, use_warped=True)
    alpha, shift, radial = u[index, :3]
    null = np.asarray((1.0 / alpha, -shift / alpha + 1.0 / radial, 0.0, 0.0))
    metric = np.asarray(
        (
            (-alpha * alpha + radial * radial * shift * shift, radial * radial * shift),
            (radial * radial * shift, radial * radial),
        )
    )
    nullness = float(null[:2] @ metric @ null[:2])
    ricci_null = _contract(result["ricci"], null)
    effective_planck = _scalar(result["F"])
    phi = _contract(result["phi_stress"], null) / effective_planck
    chi = _contract(result["chi_stress"], null) / effective_planck
    nonminimal = -_contract(result["nonminimal_term"], null) / effective_planck
    gauss_bonnet = -_contract(result["gb_residual_term"], null) / effective_planck
    reconstructed = phi + chi + nonminimal + gauss_bonnet
    return {
        "index": int(index),
        "coordinate_radius": float(radii[index]),
        "areal_radius": float(u[index, 3]),
        "phi": float(u[index, 4]),
        "outgoing_null_vector": null.tolist(),
        "null_residual": nullness,
        "effective_planck_coefficient": effective_planck,
        "ricci_null": ricci_null,
        "contributions_to_ricci_null": {
            "phi_stress": phi,
            "chi_stress": chi,
            "negative_nonminimal_term": nonminimal,
            "negative_gauss_bonnet_term": gauss_bonnet,
        },
        "reconstructed_ricci_null": reconstructed,
        "reconstruction_error": abs(reconstructed - ricci_null),
        "complete_Q_at_slice": -0.5
        * float((2.0 * ((p[index, 3] - shift * q[index, 3]) / alpha + q[index, 3] / radial) / u[index, 3]) ** 2)
        - ricci_null,
    }


def run(*, source_tag: str | None = None) -> dict[str, object]:
    source = (
        SOURCE
        if source_tag is None
        else ROOT / f"runs/fgc-2-sf1/fgcqr-activation-a6042-{source_tag}-history.npz"
    )
    output = (
        OUTPUT
        if source_tag is None
        else ROOT / f"runs/fgc-2-sf1/fgcqr-null-decomposition-a6042-{source_tag}.json"
    )
    artifact_id = (
        ARTIFACT_ID
        if source_tag is None
        else f"FGC-1-FGCQR-NULL-DECOMPOSITION-A6042-{source_tag.upper()}"
    )
    schema = f"{artifact_id}-v1"
    if output.exists():
        raise RuntimeError("null decomposition DEV1 one-shot output already exists")
    with np.load(source, allow_pickle=False) as data:
        arrays = {name: np.array(data[name], copy=True) for name in data.files}
    u, p, q = arrays["u"][-1], arrays["p"][-1], arrays["q"][-1]
    state = EvolutionState(u, p, q)
    observables = radial_null_observables(state)
    radii = arrays["radii"]
    trapped = np.flatnonzero(
        (radii > 0.0)
        & (radii <= 24.0)
        & (observables.theta_plus < 0.0)
        & (observables.theta_minus < 0.0)
    )
    points = [
        _point(
            index,
            radii,
            u,
            p,
            q,
            arrays["a"][-1],
            arrays["p_r"][-1],
            arrays["q_r"][-1],
        )
        for index in trapped
    ]
    dominant = max(
        (
            (name, max(abs(point["contributions_to_ricci_null"][name]) for point in points))
            for name in points[0]["contributions_to_ricci_null"]
        ),
        key=lambda item: item[1],
    )
    record = {
        "artifact_id": artifact_id,
        "schema": schema,
        "source_history": str(source.relative_to(ROOT)),
        "source_time": float(arrays["times"][-1]),
        "action": ACTION,
        "trapped_point_count": len(points),
        "points": points,
        "dominant_absolute_contribution": {
            "term": dominant[0],
            "magnitude": dominant[1],
        },
        "classification": "complete_FGCQR_null_Ricci_action_term_decomposition",
        "terminal": True,
        "nonclaims": {"metric_null_defocusing_demonstrated": False},
    }
    output.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-tag")
    args = parser.parse_args()
    record = run(source_tag=args.source_tag)
    print(
        json.dumps(
            {
                "dominant": record["dominant_absolute_contribution"],
                "points": [
                    {
                        "radius": point["coordinate_radius"],
                        "ricci_null": point["ricci_null"],
                        "Q": point["complete_Q_at_slice"],
                        "contributions": point["contributions_to_ricci_null"],
                        "reconstruction_error": point["reconstruction_error"],
                    }
                    for point in record["points"]
                ],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
