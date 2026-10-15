#!/usr/bin/env python3
"""Record the two-sided mass/displacement identities and their model inputs.

--check authenticates the saved result; --check-algebra also repeats its
bounded symbolic calculation. Neither mode integrates a waveform or metric.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-cap-wall-response-v1"
RELATIVE_OUTPUT = "results/development/" + SCHEMA
DEFAULT = ROOT / RELATIVE_OUTPUT
OWNER = "docs/nsc-atlas-cap-wall-response.md"
SOURCES = (
    OWNER,
    "scripts/record_nsc_atlas_cap_wall_response.py",
    "src/recursive_horizons/nsc_atlas_cap_wall_response.py",
    "src/recursive_horizons/nsc_atlas_causal_reference.py",
    "src/recursive_horizons/nsc_cosmological_atlas.py",
    "tests/test_nsc_atlas_cap_wall_response.py",
    "tests/test_nsc_cosmological_atlas.py",
    "tests/test_record_nsc_atlas_cap_wall_response.py",
)
PREDECESSOR = "results/development/nsc-atlas-wall-mass-closure-v1"
PREDECESSOR_SHA256 = "23159f708b77539d3e47f9b584bed3854ced9d8bb6a71179bb6ac38733f0a4d5"
CPU_SECONDS = 60


def canonical_bytes(value):
    return (json.dumps(value, sort_keys=True, indent=2, ensure_ascii=False,
                       allow_nan=False) + "\n").encode("utf-8")


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def read_relative(relative):
    path = ROOT / relative
    if (Path(relative).is_absolute() or ".." in Path(relative).parts
            or path.is_symlink() or not path.resolve().is_relative_to(ROOT.resolve())
            or not path.is_file()):
        raise ValueError("missing or unsafe source: " + relative)
    return path.read_bytes()


def source_bindings():
    return {path: sha(read_relative(path)) for path in sorted(SOURCES)}


def authenticate_predecessor():
    raw = read_relative(PREDECESSOR + "/manifest.json")
    if sha(raw) != PREDECESSOR_SHA256:
        raise ValueError("exterior predecessor manifest mismatch")
    manifest = json.loads(raw)
    for path, expected in manifest["source_bindings"].items():
        if sha(read_relative(path)) != expected:
            raise ValueError("exterior predecessor source mismatch: " + path)
    for name, binding in manifest["outputs"].items():
        payload = read_relative(PREDECESSOR + "/" + name)
        if sha(payload) != binding["sha256"] or len(payload) != binding["bytes"]:
            raise ValueError("exterior predecessor output mismatch: " + name)
    spec = importlib.util.spec_from_file_location(
        "cap_wall_exterior_predecessor", ROOT / "scripts/record_nsc_atlas_wall_mass_closure.py")
    predecessor = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(predecessor)
    predecessor.check(ROOT / PREDECESSOR)
    return {"manifest": PREDECESSOR + "/manifest.json", "sha256": sha(raw),
            "role": "exterior functional, quiet anchor and original scalar/stress context"}


def declaration():
    return {
        "schema": SCHEMA,
        "status": "CONDITIONAL_TWO_SIDED_MASS_AND_WALL_DISPLACEMENT_FUNCTIONALS",
        "owner": OWNER,
        "clock": "s=tau_wall/ell, shared quiet corner rho(s_a)=0",
        "order": "r=x+alpha*rho, alpha=kappa_4*epsilon**2/ell**2",
        "cap_mass": "m_minus/ell=h**2*r**3/2+alpha*nu_minus+O(alpha**2)",
        "cap_pulled_coefficient": "nu_minus+3*h**2*x**2*rho/2",
        "mass_difference": "D=mu_plus-nu_minus",
        "displacement": "rho_s=v_x*rho+(mu_plus/x+beta_plus*D/(b*x**2))/v",
        "signed_normal_response": "delta_beta_pm=beta_pm_x*rho+D/(b*x**2)",
        "acceleration_response": "delta_a_pm=6*m*rho/(b*x**4)-2*D/(b*x**3)",
        "mean_curvature_response": "delta_mean_K=0 (includes metric and wall displacement)",
        "domain": {"x": [0.45, 0.55], "h": 1, "b": 2, "m": 0.4,
                   "source": "original C3 compact preparation, quiet cap and bottom, no-return"},
        "derived": ["cap vacuum-subtracted mass functional",
                    "quiet-anchored wall displacement functional",
                    "compatible angular, temporal and normal junction conditions"],
        "remaining": ["controlled numerical evaluation", "compatible bulk metric evolution",
                      "corrected observer clock and field response", "finite-amplitude remainder bound",
                      "formation, stability and physical source selection"],
        "execution": {"cpu_seconds_max": CPU_SECONDS, "numerical_integrations": 0,
                      "waveform_solves": 0, "bulk_metric_solves": 0},
    }


def algebra():
    from recursive_horizons.nsc_atlas_cap_wall_response import exact_algebra_checks
    return exact_algebra_checks()


def expected_certificate():
    from recursive_horizons.nsc_atlas_cap_wall_response import expected_certificate as expected
    return expected()


def model():
    from recursive_horizons.nsc_cosmological_atlas import model_specification
    return model_specification()


def create(directory):
    started = time.process_time()
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only evidence: use --check or a fresh directory")
    sources, predecessor = source_bindings(), authenticate_predecessor()
    certificate = algebra()
    if certificate != expected_certificate():
        raise ValueError("cap-wall algebra certificate mismatch")
    record = {**declaration(), "exact_algebra": certificate}
    outputs = {"cap-wall.json": canonical_bytes(record),
               "model-specification.json": canonical_bytes(model())}
    if sources != source_bindings() or predecessor != authenticate_predecessor():
        raise RuntimeError("source or prerequisite changed during recording")
    if time.process_time() - started > CPU_SECONDS:
        raise RuntimeError("cap-wall record exceeded process CPU budget")
    manifest = {"schema": SCHEMA + "-manifest", "source_bindings": sources,
                "predecessor": predecessor,
                "outputs": {name: {"sha256": sha(raw), "bytes": len(raw)}
                            for name, raw in outputs.items()}}
    directory.mkdir(parents=True, exist_ok=False)
    for name, raw in {**outputs, "manifest.json": canonical_bytes(manifest)}.items():
        with (directory / name).open("xb") as stream:
            stream.write(raw)
    return {"created": directory.name, "outputs": len(outputs),
            "cpu_seconds": time.process_time() - started, "numerical_integrations": 0}


def check(directory):
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_bytes())
    if manifest.get("schema") != SCHEMA + "-manifest":
        raise ValueError("cap-wall manifest schema mismatch")
    if manifest.get("source_bindings") != source_bindings():
        raise ValueError("cap-wall source binding mismatch")
    if manifest.get("predecessor") != authenticate_predecessor():
        raise ValueError("cap-wall predecessor mismatch")
    expected = {"cap-wall.json", "model-specification.json"}
    if (set(manifest.get("outputs", {})) != expected
            or {p.name for p in directory.iterdir()} != expected | {"manifest.json"}):
        raise ValueError("cap-wall output inventory mismatch")
    records = {}
    for name, binding in manifest["outputs"].items():
        path = directory / name
        if path.is_symlink():
            raise ValueError("unsafe cap-wall output")
        raw = path.read_bytes()
        if sha(raw) != binding["sha256"] or len(raw) != binding["bytes"]:
            raise ValueError("cap-wall artifact binding mismatch: " + name)
        records[name] = json.loads(raw)
    record = dict(records["cap-wall.json"])
    certificate = record.pop("exact_algebra", None)
    if record != declaration() or certificate != expected_certificate():
        raise ValueError("cap-wall declaration mismatch")
    if records["model-specification.json"] != model():
        raise ValueError("cap-wall model specification mismatch")
    return {"checked": directory.name, "sources": len(SOURCES), "outputs": 2,
            "numerical_integrations": 0}


def check_algebra(directory):
    started = time.process_time()
    check(directory)
    saved = json.loads((Path(directory) / "cap-wall.json").read_bytes())["exact_algebra"]
    if algebra() != saved:
        raise ValueError("cap-wall algebra replay differs from record")
    if time.process_time() - started > CPU_SECONDS:
        raise RuntimeError("cap-wall replay exceeded process CPU budget")
    return {"replayed": Path(directory).name, "cpu_seconds": time.process_time() - started}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--check", action="store_true")
    modes.add_argument("--check-algebra", action="store_true")
    args = parser.parse_args(argv)
    action = check_algebra if args.check_algebra else check if args.check else create
    print(json.dumps(action(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
