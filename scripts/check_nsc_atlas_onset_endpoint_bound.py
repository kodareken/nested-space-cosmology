#!/usr/bin/env python3
"""One leading wall-field enclosure from authenticated source/reference bounds."""
import argparse
from fractions import Fraction as Q
import hashlib
from importlib.metadata import version
import json
import math
from pathlib import Path
import subprocess
import sys
import time

from flint import ctx

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_onset_endpoint_bound import monotone_wall_enclosure
from recursive_horizons.nsc_atlas_source_cell_bound import _ball, _coordinate, _pulse_jets
from recursive_horizons.nsc_source_binding_graph import original_source_graph

ROOT = Path(__file__).resolve().parents[1]
INPUTS = {
    'results/summaries/nsc-atlas-source-cell-bound-v2.json':
        '707bb8d76359edaf2a26228657a219933ebd164edcf89eb87742f9e66bec2bdd',
    'results/summaries/nsc-atlas-reference-gauge-bound-v1.json':
        'e8a8400cec0f6ae0bd831723696ebf1d41c38bdf0b052301ffab837465f1d216',
}
OWNERS = (
    'scripts/check_nsc_atlas_onset_endpoint_bound.py',
    'src/recursive_horizons/nsc_atlas_onset_endpoint_bound.py',
    'src/recursive_horizons/nsc_atlas_source_cell_bound.py',
    'src/recursive_horizons/nsc_atlas_reference_gauge_bound.py',
    'src/recursive_horizons/nsc_atlas_characteristic_initialization.py',
    'src/recursive_horizons/nsc_atlas_transmission_bound.py',
    'src/recursive_horizons/nsc_source_binding_graph.py',
    'scripts/run_python.py', 'pyproject.toml',
)


def plain(value):
    if isinstance(value, Q):
        return str(value)
    if isinstance(value, dict):
        return {key:plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [plain(item) for item in value]
    return value


def rational_end(value):
    if not value.is_finite():
        raise ValueError('nonfinite source endpoint')
    mantissa, exponent = value.man_exp()
    return Q(int(mantissa))*Q(2)**int(exponent)


def interval(value):
    return rational_end(value.lower()), rational_end(value.upper())


def display(bounds):
    lo, hi = bounds['lower'], bounds['upper']
    return dict(lower=0. if lo == 0 else math.nextafter(float(lo), -math.inf),
                upper=0. if hi == 0 else math.nextafter(float(hi), math.inf))


def read_bound(name, digest):
    raw = (ROOT/name).read_bytes()
    if hashlib.sha256(raw).hexdigest() != digest:
        raise ValueError('original input certificate changed: '+name)
    return json.loads(raw)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or not args.output.parent.is_dir():
        raise ValueError('require an existing parent and a new output path')
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    bindings = original_source_graph(ROOT, OWNERS, commit)
    for name, digest in bindings.items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
            raise ValueError('commit calculation source before execution: '+name)
    source, gauge = [read_bound(name, digest) for name, digest in INPUTS.items()]
    source_input = source['input']
    original = read_bound(source_input['path'], source_input['sha256'])
    if (source['status'] != 'certified_source_representation_cell'
            or gauge['status'] != 'certified_exact_reference_gauge_cell'
            or source['cell'] != gauge['cell']):
        raise ValueError('require the same original certified onset cell')
    geometry = original['geometry']
    cell = gauge['certificate']['inputs']
    left, right, rl, rh = (Q(cell[key]) for key in ('left', 'right', 'radius_lower', 'radius_upper'))
    if (left != Q(source['cell']['left']) or right != Q(source['cell']['right'])
            or not Q(geometry['U0']) < 0
            or not Q(geometry['Va']) < Q(geometry['Vb']) < 0):
        raise ValueError('source labels or original IV source signs changed')
    budget = Budget(maximum_evaluations=4096, cpu_seconds=5.)
    report = dict(schema='nsc-atlas-onset-endpoint-bound-v1', source_commit=commit,
                  source_bindings=bindings, status='incomplete',
                  runtime=dict(python=sys.version.split()[0], python_flint=version('python-flint'),
                               platform=sys.platform),
                  inputs=[dict(path=name, sha256=digest) for name, digest in INPUTS.items()]
                         + [dict(path=source_input['path'], sha256=source_input['sha256'])],
                  limits=dict(charged_values=4096, cpu_seconds=5., precision_bits=160,
                              source_endpoints=2, physical_events=1),
                  parameters=plain(dict(left=left, right=right, radius_lower=rl, radius_upper=rh,
                      source_map_scale=Q(.8), reference_scale=Q(4, 5),
                      source_geometry={key:Q(geometry[key]) for key in ('U0', 'Va', 'Vb')})))
    try:
        budget.check(1024)
        with ctx.workprec(160):
            L, R = _ball(left), _ball(right)
            zl, zxl, _ = _coordinate(L, geometry)
            zr, zxr, _ = _coordinate(R, geometry)
            if not zl < -1 or not zr > -1 or not zr < 0 or not zxl > 0 or not zxr > 0:
                raise ValueError('require a quiet left endpoint and increasing first-half pulse')
            g, gx, _ = _pulse_jets(R, R, geometry)
            if not g > 0 or not gx > 0:
                raise ValueError('endpoint source value and derivative must be positive')
            values = dict(g=interval(g), g_x=interval(gx))
            report['source'] = plain(dict(**values, compact_coordinate_left=interval(zl),
                                         compact_coordinate_right=interval(zr),
                                         monotone_first_half=True, quiet_left=True))
        certificate = monotone_wall_enclosure(left=left, right=right, radius_lower=rl,
            radius_upper=rh, **values, budget=budget)
        report.update(status='certified_leading_wall_endpoint', certificate=plain(certificate),
                      display={key:display(certificate[key]) for key in
                               ('wall_value', 'wall_derivative', 'source_minus_wall',
                                'source_derivative_minus_wall_derivative')})
        budget.check()
    except (ValueError, RuntimeError, KeyError, TypeError, OverflowError) as error:
        report.update(status='stopped', failure=str(error))
    report['execution'] = dict(charged_values=budget.evaluations,
        cpu_seconds=time.process_time()-budget.started, field_evolutions=0,
        source_integrations=0, time_steps=0, root_searches=0, retries=0, raw_history_saved=False)
    report['regenerate'] = ('python scripts/run_python.py scripts/check_nsc_atlas_onset_endpoint_bound.py'
                            ' --output <new-output-path>')
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+'\n',
                           encoding='utf-8', newline='\n')
    print(json.dumps({key:report.get(key) for key in ('status', 'display', 'execution', 'failure')}))
    return 0 if report['status'] == 'certified_leading_wall_endpoint' else 1


if __name__ == '__main__':
    raise SystemExit(main())
