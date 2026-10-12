#!/usr/bin/env python3
"""Save original-producer output for diagnosing a failed replay; never accept it."""
from __future__ import annotations

import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import subprocess
import sys

import replay_history as replay


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', type=int, choices=range(1, 7), default=5)
    args = parser.parse_args()
    replay.prepare()
    destination = replay.ROOT / '.build/history-diagnostics' / f'v{args.version}'
    destination.mkdir(parents=True, exist_ok=False)
    suffix = '' if args.version == 1 else f'_v{args.version}'
    source = f'scripts/derive_nsc_cosmological_atlas{suffix}.py'
    result = subprocess.run(
        [sys.executable, source, '--output', str(destination / 'regenerated')],
        cwd=replay.CONTEXT / 'lab', env=replay.environment(replay.CONTEXT, True),
        capture_output=True, text=True,
    )
    import matplotlib
    from matplotlib import font_manager, ft2font
    from numpy._core._multiarray_umath import __cpu_features__
    from PIL import features
    font = Path(font_manager.findfont('DejaVu Sans'))
    report = {
        'purpose': 'Diagnostics only; unchanged historical checks remain authoritative.',
        'source_commit': replay.BASELINE,
        'python': platform.python_version(), 'platform': platform.platform(),
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
        'returncode': result.returncode, 'stdout': result.stdout, 'stderr': result.stderr,
    }
    (destination / 'runtime.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2))
    raise SystemExit(result.returncode)


if __name__ == '__main__':
    main()
