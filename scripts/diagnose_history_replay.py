#!/usr/bin/env python3
"""Save original-producer output for diagnosing a failed replay; never accept it."""
from __future__ import annotations

import argparse
from contextlib import redirect_stdout
import hashlib
from importlib.metadata import version
import io
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys

import replay_history as replay


def json_mismatches(stored, regenerated, *, relative=1e-10, absolute=1e-12, limit=20) -> dict:
    """Locate differing leaves with the original atlas comparator's semantics.

    This is a diagnostic report only. It never replaces or changes the original
    checker, its tolerances, or the frozen record.
    """
    differences, count = [], 0

    def add(path, reason, actual, expected):
        nonlocal count
        count += 1
        if len(differences) < limit:
            row = {"path": path, "reason": reason, "stored": actual, "regenerated": expected}
            if (isinstance(actual, (int, float)) and isinstance(expected, (int, float))
                    and math.isfinite(actual) and math.isfinite(expected)):
                row["absolute_difference"] = abs(actual-expected)
                row["allowed_difference"] = max(absolute, relative*max(abs(actual), abs(expected)))
            differences.append(row)

    def visit(actual, expected, path):
        if isinstance(expected, bool) or expected is None or isinstance(expected, str):
            if actual != expected:
                add(path, "value", actual, expected)
        elif isinstance(expected, (int, float)):
            if not isinstance(actual, (int, float)) or not math.isclose(actual, expected, rel_tol=relative, abs_tol=absolute):
                add(path, "number", actual, expected)
        elif isinstance(expected, list):
            if not isinstance(actual, list):
                add(path, "type", type(actual).__name__, "list")
            elif len(actual) != len(expected):
                add(path, "list_length", len(actual), len(expected))
            else:
                for index, (left, right) in enumerate(zip(actual, expected)):
                    visit(left, right, f"{path}[{index}]")
        elif isinstance(expected, dict):
            if not isinstance(actual, dict):
                add(path, "type", type(actual).__name__, "dict")
            elif actual.keys() != expected.keys():
                add(path, "key_inventory", sorted(actual), sorted(expected))
            else:
                for key in sorted(expected):
                    visit(actual[key], expected[key], path + "." + key)
        else:
            add(path, "unsupported_json_type", type(actual).__name__, type(expected).__name__)

    visit(stored, regenerated, "$")
    return {"mismatch_count": count, "reported": differences, "truncated": count > len(differences),
            "comparison": {"relative": relative, "absolute": absolute}}


def target(atlas_version: int) -> dict:
    """Select an original context explicitly; no live-source fallback exists."""
    if atlas_version not in range(1, 9):
        raise ValueError("unsupported historical atlas version")
    if atlas_version >= 7:
        prepare = getattr(replay, f"prepare_v{atlas_version}")
        context = getattr(replay, f"V{atlas_version}_CONTEXT")
        bundle = getattr(replay, f"V{atlas_version}_BUNDLE")
        info = prepare(context)
        owner = context
        env = replay.environment(context)
        env["GIT_ALTERNATE_OBJECT_DIRECTORIES"] = str(replay.ROOT / "archive/source-history/objects")
        env["PYTHONNOUSERSITE"] = "1"
    else:
        info = replay.prepare()
        context, owner = replay.CONTEXT, replay.CONTEXT / "lab"
        bundle = replay.BUNDLES[atlas_version-1]
        env = replay.environment(context, True)
    suffix = "" if atlas_version == 1 else f"_v{atlas_version}"
    return {"source_commit": info["source_commit"], "context": context, "owner": owner,
            "bundle": context / bundle, "environment": env,
            "source": f"scripts/derive_nsc_cosmological_atlas{suffix}.py"}


def runtime_metadata() -> dict:
    import matplotlib
    from matplotlib import font_manager, ft2font
    import numpy as np
    from numpy._core._multiarray_umath import __cpu_features__
    from PIL import features
    import scipy
    # Load its BLAS wrapper without running an operation, so threadpoolctl can
    # report SciPy's dispatcher alongside NumPy's separately linked library.
    import scipy.linalg
    font = Path(font_manager.findfont('DejaVu Sans'))
    configurations = {}
    for name, capture in (("numpy_runtime", np.show_runtime), ("numpy_build", np.show_config),
                          ("scipy_build", scipy.show_config)):
        stream = io.StringIO()
        with redirect_stdout(stream):
            capture()
        configurations[name] = stream.getvalue()
    report = {
        'python': platform.python_version(), 'platform': platform.platform(),
        'python_build': sys.version,
        'machine': platform.machine(), 'processor': platform.processor(),
        'dependencies': {name: version(name) for name in (
            'numpy', 'scipy', 'sympy', 'matplotlib', 'Pillow', 'contourpy',
            'fonttools', 'kiwisolver', 'pyparsing',
        )},
        'freetype': ft2font.__freetype_version__,
        'pillow_zlib': features.version('zlib'),
        'numpy_cpu_features': __cpu_features__,
        'matplotlib_backend': matplotlib.get_backend(),
        'dejavu_sans_sha256': hashlib.sha256(font.read_bytes()).hexdigest(),
        'configurations': configurations,
        'runtime_probe_imports': ['numpy', 'scipy.linalg'],
        'numeric_environment': {name: os.environ[name] for name in (
            'OPENBLAS_CORETYPE', 'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS',
            'NPY_DISABLE_CPU_FEATURES', 'MKL_CBWR', 'MKL_NUM_THREADS',
        ) if name in os.environ},
    }
    try:
        from threadpoolctl import threadpool_info
        report['threadpools'] = threadpool_info()
    except ImportError:
        report['threadpools'] = 'threadpoolctl unavailable; build configuration remains recorded'
    return report


def compare_outputs(original: Path, regenerated: Path) -> dict:
    manifest = json.loads((original / 'manifest.json').read_bytes())
    tolerance = manifest['numerical_replay_tolerance']
    comparison = {}
    for name in manifest['outputs']:
        if Path(name).name != name:
            raise ValueError('unsafe historical output name')
        stored = (original / name).read_bytes()
        current = regenerated / name
        if not current.is_file():
            comparison[name] = {'missing_regenerated_output': True}
            continue
        raw = current.read_bytes()
        row = {'stored_sha256': hashlib.sha256(stored).hexdigest(),
               'regenerated_sha256': hashlib.sha256(raw).hexdigest(),
               'byte_equal': stored == raw}
        if name.endswith('.json'):
            row['json_comparison'] = json_mismatches(
                json.loads(stored), json.loads(raw), relative=tolerance['relative'],
                absolute=tolerance['absolute'])
        comparison[name] = row
    return comparison


def diagnose(atlas_version: int, destination: Path | None = None) -> dict:
    chosen = target(atlas_version)
    destination = (replay.ROOT / '.build/history-diagnostics' / f'v{atlas_version}'
                   if destination is None else Path(destination))
    destination.mkdir(parents=True, exist_ok=False)
    result = subprocess.run(
        [sys.executable, chosen['source'], '--output', str(destination / 'regenerated')],
        cwd=chosen['owner'], env=chosen['environment'], capture_output=True, text=True,
    )
    report = {
        'purpose': 'Diagnostics only; unchanged historical checks remain authoritative.',
        'accepted': False,
        'source_commit': chosen['source_commit'], 'atlas_version': atlas_version,
        'returncode': result.returncode, 'stdout': result.stdout, 'stderr': result.stderr,
    }
    try:
        report.update(runtime_metadata())
    except Exception as error:
        report['runtime_capture_error'] = f'{type(error).__name__}: {error}'
    try:
        report['outputs'] = compare_outputs(chosen['bundle'], destination / 'regenerated')
    except Exception as error:
        report['comparison_capture_error'] = f'{type(error).__name__}: {error}'
    (destination / 'runtime.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    return report


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', type=int, choices=range(1, 9), default=5)
    args = parser.parse_args(argv)
    report = diagnose(args.version)
    print(json.dumps(report, indent=2))
    # Diagnostic success means only that the original producer ran. The failed
    # authoritative replay remains failed and no regenerated bytes are adopted.
    raise SystemExit(report['returncode'])


if __name__ == '__main__':
    main()
