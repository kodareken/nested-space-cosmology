#!/usr/bin/env python3
"""One original source-cell interval certificate; no physical solve or scan."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

from diagnose_nsc_atlas_boundary_rows import authenticate, producer
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_characteristic_trace import SourceSnapshot
from recursive_horizons.nsc_atlas_source_cell_bound import source_cell_error_bound
from recursive_horizons.nsc_source_binding_graph import original_source_graph

ROOT = Path(__file__).resolve().parents[1]
OWNERS = ('scripts/check_nsc_atlas_source_cell_bound.py',
          'scripts/diagnose_nsc_atlas_boundary_rows.py',
          'src/recursive_horizons/nsc_atlas_characteristic_trace.py',
          'src/recursive_horizons/nsc_atlas_source_cell_bound.py',
          'pyproject.toml', 'scripts/run_python.py')


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('preserve the existing result; choose a new admitted output')
    if not args.output.parent.is_dir():
        raise ValueError('output parent must already exist')
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    owners = set(OWNERS)
    for module in tuple(sys.modules.values()):
        name = getattr(module, '__file__', None)
        if name is None:
            continue
        path = Path(name).resolve()
        if path.is_relative_to(ROOT) and path.suffix == '.py':
            relative = path.relative_to(ROOT)
            if relative.parts[0] in ('src', 'scripts'):
                owners.add(relative.as_posix())
    bindings = original_source_graph(ROOT, tuple(sorted(owners)), commit)
    for name, expected in bindings.items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != expected:
            raise ValueError('commit calculation source before execution: '+name)
    budget = Budget(maximum_evaluations=100000, cpu_seconds=30.)
    report = dict(schema='nsc-atlas-source-cell-bound-v2', source_commit=commit, source_bindings=bindings,
                  runtime=producer.runtime_versions(), status='incomplete',
                  limits=dict(charged_values=100000, cpu_seconds=30., source_cells=1, panels=2))
    try:
        _, _, _, identity = authenticate(ROOT, budget)
        raw = (ROOT/producer.BASE_INPUT).read_bytes()
        if hashlib.sha256(raw).hexdigest() != identity['source_bindings'][producer.BASE_INPUT]:
            raise ValueError('original source input changed')
        original = json.loads(raw); data = original['initial_data']['finite']
        rows = {k: np.asarray(v, dtype=float)[::2] for k, v in data['exact_rows'].items()}
        nodes = rows['x']; geometry = original['geometry']; front = geometry['support_x'][0]
        j = int(np.searchsorted(nodes, front)-1)
        if j != 35 or len(nodes) != 129 or rows['h'][j] != 0 or rows['h_x'][j] != 0:
            raise ValueError('original source cell or quiet jet changed')
        snapshot = SourceSnapshot(geometry, alpha=data['alpha'], x_nodes=nodes, source_rows=rows, budget=budget)
        left, right = Q(float(nodes[j])), Q(float(nodes[j+1])); width = right-left
        # Inspect the actual PPoly coefficients; there is no second fit.
        coefficients = snapshot._h.c[:, j]
        power = tuple(Q(float(coefficients[3-i]))*width**i for i in range(4))
        right_jet = tuple(Q(float(snapshot._h.c[k, j+1])) for k in (3, 2))
        if right_jet != (Q(float(rows['h'][j+1])), Q(float(rows['h_x'][j+1]))):
            raise ValueError('right PPoly owner differs from the original source jet')
        report.update(input=dict(path=producer.BASE_INPUT, sha256=hashlib.sha256(raw).hexdigest(),
                operator_sha256=identity['sha256'], manifest_sha256=identity['manifest_sha256']),
            cell=dict(index=j, left=float(left), right=float(right), front=front,
                alpha=data['alpha'], source_knots=129, normalized_exact_binary_powers=[str(v) for v in power],
                right_owner_value_and_first=[str(v) for v in right_jet]),
            polynomial='actual SourceSnapshot PPoly binary coefficients; original scalar source formula')
        report['certificate'] = source_cell_error_bound(geometry, alpha=data['alpha'], left=left, right=right,
                                                        power=power, right_jet=right_jet, budget=budget)
        budget.check(); report['status'] = 'certified_source_representation_cell'
    except (ValueError, RuntimeError, OverflowError, FloatingPointError, KeyError, TypeError) as error:
        report['status'] = 'stopped'; report['failure'] = str(error)
    report['execution'] = dict(charged_values=budget.evaluations, cpu_seconds=time.process_time()-budget.started,
        source_integrations=0, source_queries=0, time_steps=0, retries=0, raw_history_saved=False)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(dict(status=report['status'], execution=report['execution'], failure=report.get('failure'))))
    return 0 if report['status'] == 'certified_source_representation_cell' else 1


if __name__ == '__main__':
    raise SystemExit(main())
