#!/usr/bin/env python3
"""Actual geometry-derived vertices and their selected tree exchange."""
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from recursive_horizons.nsc_atlas_cosmohedron_interaction import (
    body_with_child, first_axial_mode, interaction,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT/"results/summaries/nsc-atlas-cosmohedron-interaction-v1.json"
SOURCES = ["src/recursive_horizons/nsc_atlas_cosmohedron_interaction.py",
           "src/recursive_horizons/nsc_atlas_gravitational_response.py",
           "scripts/check_nsc_atlas_cosmohedron_interaction.py",
           "tests/test_nsc_atlas_cosmohedron_interaction.py"]
TOLERANCES = (1e-8, 2e-10, 2e-12)


def serializable(value):
    if isinstance(value, dict):
        return {k: serializable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serializable(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def evaluate(job):
    child_radius, scale, rtol = job
    body = body_with_child(child_radius, scale)
    mode = first_axial_mode(body, rtol=rtol)
    result = interaction(mode)
    result.update(inner_radius=body.inner.radius, mass_scale=scale, rtol=rtol,
                  norm_quadrature=mode.norm_quadrature(64),
                  mode_endpoint_gap=float(mode.end[0]-np.pi),
                  rayleigh_upper=mode.rayleigh_upper,
                  interfaces=[body.inner.record(), body.outer.record()],
                  rhs_calls=mode.rhs_calls, mode_cpu_seconds=mode.cpu_seconds)
    return result


def main():
    started = time.perf_counter()
    jobs = [(r, 1., tol) for r in (.25, .30) for tol in TOLERANCES]
    jobs.append((.25, 2., TOLERANCES[-1]))
    with ProcessPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(evaluate, jobs))
    cpu = sum(r["mode_cpu_seconds"] for r in rows)
    rhs = sum(r["rhs_calls"] for r in rows)
    if cpu > 180 or rhs > 400000:
        raise RuntimeError("declared interaction allowance exhausted")
    families = []
    for radius in (.25, .30):
        levels = [r for r in rows if r["inner_radius"] == radius and r["mass_scale"] == 1.]
        fine = levels[-1]
        indicators = {}
        for key in ("gravity_frequency", "coefficient_shape", "correlation_shape"):
            coarse_gap = abs(levels[0][key]-levels[1][key])
            fine_gap = abs(levels[1][key]-fine[key])
            if fine_gap > coarse_gap:
                raise RuntimeError("interaction refinement failed; diagnose before further calculation")
            indicators[key] = max(2*fine_gap, abs(fine[key])*1e-10)
        if abs(fine["norm_quadrature"]-1) > 2e-8:
            raise RuntimeError("independent mode normalization failed")
        families.append(dict(child_radius=radius, levels=levels, indicators=indicators))
    a, b = [f["levels"][-1] for f in families]
    scaled = rows[-1]
    scale_gap = abs(scaled["coefficient_shape"]/a["coefficient_shape"]-1)
    if scale_gap > 1e-8:
        raise RuntimeError("fixed-G geometric scaling failed")
    record = dict(
        schema="NSC-ATLAS-COSMOHEDRON-INTERACTION-v1", date="2026-10-27",
        preparation=dict(theory="Einstein/four-form/supporting surface scalar + one zero-background canonical surface probe per wall",
                         scope="conditional finite positive Gaussian sector; selected first axial-mode tree exchange only",
                         outer_worldtube_over_M=10., outer_mode_boundary="Dirichlet axial master",
                         gravity_multipole=[2, 0], scalar_harmonics=["real l1m1 cosine", "real l2m1 sine"],
                         scalar_zero_mode="outside retained sector", source_couplings="derived, no cubic fit",
                         state="positive-frequency reference from the common static time",
                         full_Einstein_four_point=False, quantum_backreaction_included=False,
                         observational_validation=False),
        families=families, mass_rescaling_control=scaled,
        comparison=dict(coefficient_relative_change=b["coefficient_shape"]/a["coefficient_shape"]-1,
                        correlation_relative_change=b["correlation_shape"]/a["correlation_shape"]-1,
                        axial_frequency_relative_change=b["gravity_frequency"]/a["gravity_frequency"]-1,
                        same_exterior_and_outer_wall=True,
                        note="separately prepared static states; no instantaneous change or causal signal claim"),
        inheritance=dict(shape="Gamma_selected M^4/G invariant at fixed dimensionless geometry",
                         normalized_selected_moment="C_selected=(G/M^2)*shape in hbar=c=1 units",
                         fixed_G_full_strength_invariant=False, scale_relative_gap=scale_gap),
        attribution=dict(cosmohedra="https://arxiv.org/html/2412.19881v1; Eq49 energy kernel",
                         gravity="https://arxiv.org/html/gr-qc/0502028; Eq66 master relation, Eq84 flux normalization",
                         scope="one graph energy kernel; no full positive-geometry construction for Einstein theory"),
        compute=dict(rhs_calls=rhs, cpu_seconds=cpu, workers=4,
                     elapsed_seconds=time.perf_counter()-started, rhs_limit=400000, cpu_limit_seconds=180,
                     note="mode solver process CPU excludes imports and independent cap-norm quadrature"),
        source_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES},
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(serializable(record), indent=2, allow_nan=False)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps(dict(output=OUTPUT.relative_to(ROOT).as_posix(), rhs=rhs, cpu=cpu,
                          base_shape=a["coefficient_shape"], control_shape=b["coefficient_shape"],
                          contrast=record["comparison"], error_indicators=[f["indicators"] for f in families])))


if __name__ == "__main__":
    main()
