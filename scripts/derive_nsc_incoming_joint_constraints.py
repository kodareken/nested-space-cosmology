#!/usr/bin/env python3
"""Prepare/replay the missing reference response and joint incoming assembly."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_incoming_cauchy_jets import incoming_cauchy_jets, IncomingNormalJetChange
from recursive_horizons.nsc_incoming_reference_response import incoming_reference_response
from recursive_horizons.nsc_incoming_joint_constraints import baseline_matter_source, source_action_gradient, joint_constraint_approximant
from recursive_horizons.nsc_incoming_local_constraints import incoming_local_constraints
from recursive_horizons.nsc_light_restoration_action import LightRestorationAction
from recursive_horizons.nsc_magnetic_light_reference import magnetic_light_spectrum
from recursive_horizons.nsc_spherical_local_history import LockedSphericalLocalAction
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT = ROOT/'results/development/nsc-incoming-joint-constraints.json'
NAMES = ('nsc-incoming-spectral-source', 'nsc-incoming-subgap-source',
         'nsc-incoming-source-tail', 'nsc-incoming-local-constraints',
         'nsc-reference-band-bulk', 'nsc-mode-resolved-cauchy-state')
SOURCES = ('src/recursive_horizons/nsc_incoming_reference_response.py',
           'src/recursive_horizons/nsc_incoming_joint_constraints.py',
           'scripts/derive_nsc_incoming_joint_constraints.py',
           'tests/test_nsc_incoming_reference_response.py',
           'docs/nsc-incoming-joint-constraints.md')
OWNERS = ('src/recursive_horizons/nsc_incoming_cauchy_jets.py',
          'src/recursive_horizons/nsc_spatial_reference_symbol.py',
          'src/recursive_horizons/nsc_reference_band_action.py',
          'src/recursive_horizons/nsc_mode_resolved_cauchy_state.py')


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def jsonable(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.floating): return float(value)
    if isinstance(value, np.integer): return int(value)
    raise TypeError(type(value).__name__)


def records():
    result = {}
    for name in NAMES:
        data = json.loads((ROOT/f'results/development/{name}.json').read_text())
        for key in ('source_hashes', 'input_hashes', 'numerical_owner_and_input_hashes', 'prepared_dependency_hashes'):
            for path, digest in data.get(key, {}).items():
                if sha(path) != digest: raise ValueError('upstream dependency changed: '+path)
        payloads = list(data.get('input_payloads', []))
        for key in ('payload', 'input_payload'):
            if key in data: payloads.append(data[key])
        for payload in payloads:
            if sha(payload['path']) != payload['sha256']:
                raise ValueError('upstream artifact changed: '+payload['path'])
        result[name] = data
    return result


def signature():
    return {p: sha(p) for p in (*SOURCES, *OWNERS,
            *(f'results/development/{n}.json' for n in NAMES))}


def probe():
    return incoming_cauchy_jets((IncomingNormalJetChange('r', 1, 1, .01),))


def prepare():
    source = records()
    before = signature()
    domain = probe()
    controls, arrays = {}, {}
    for name, nodes, scale in [('n16', 16, 1.), ('n24', 24, 1.), ('scale07', 24, .7)]:
        result = incoming_reference_response(domain, source[NAMES[-1]]['channels'], nodes=nodes, scale_factor=scale)
        samples = result.pop('quadrature_samples')
        for i, row in enumerate(samples):
            prefix = f'{name}/{i}/'
            arrays[prefix+'momenta'] = row.pop('momenta')
            arrays[prefix+'weights'] = row.pop('weights')
            arrays[prefix+'insertions'] = row.pop('insertions')
        controls[name] = {'summary': result, 'signed_inventory': samples}
        print(f'{name}: integrated all 63 signed families; '+str(result['reference_action_gradient_change']), flush=True)
    ledger = source['nsc-incoming-local-constraints']['locked_inputs']
    local = incoming_local_constraints(domain, LightRestorationAction(magnetic_light_spectrum(4)),
                                       LockedSphericalLocalAction(ledger))
    if signature() != before: raise ValueError('dependency changed during reference preparation')
    meta = {'signature': before, 'controls': controls, 'local_probe': local,
            'probe_changes': domain.changed_normal_entries()}
    arrays['metadata_json'] = np.frombuffer(json.dumps(meta, sort_keys=True, default=jsonable).encode(), np.uint8)
    raw = deterministic_npz_bytes(arrays)
    digest = hashlib.sha256(raw).hexdigest()
    path = ROOT/f'results/development/artifacts/nsc-incoming-reference-response.{digest}.npz'
    if path.exists() and path.read_bytes() != raw: raise ValueError('content address collision')
    if not path.exists(): path.write_bytes(raw)
    return path


def replay(path):
    source = records()
    with np.load(path, allow_pickle=False) as f:
        arrays = {k: f[k].copy() for k in f.files}
    meta = json.loads(arrays['metadata_json'].tobytes())
    if meta['signature'] != signature(): raise ValueError('prepared reference dependency changed')
    domain = probe()
    if meta['probe_changes'] != domain.changed_normal_entries(): raise ValueError('probe changed')
    summaries = {}
    replay_error = 0.
    for name, control in meta['controls'].items():
        groups = {'0': np.zeros((4, 4))}
        for i, row in enumerate(control['signed_inventory']):
            prefix = f'{name}/{i}/'
            weight, insertions = arrays[prefix+'weights'], arrays[prefix+'insertions']
            k = arrays[prefix+'momenta']
            n = len(k)//2
            if not np.array_equal(k[:n], -k[n:]) or not np.array_equal(weight[:n], weight[n:]):
                raise ValueError('both reference momentum signs required')
            value = row['multiplicity']/(2*np.pi)*np.einsum('n,onb->ob', weight, insertions)
            group = str(row['group'])
            groups[group] = groups.get(group, np.zeros((4, 4)))+value
        orders = sum(groups.values())
        summary = control['summary']
        replay_error = max(replay_error, float(np.max(abs(orders-np.array(summary['formal_order_gradient_changes'])))))
        for group, value in groups.items():
            replay_error = max(replay_error, float(np.max(abs(value-np.array(summary['group_formal_order_gradient_changes'][group])))))
        summaries[name] = {**summary, 'reference_action_gradient_change': orders.sum(axis=0),
                           'formal_order_gradient_changes': orders, 'group_formal_order_gradient_changes': groups}
    fine = summaries['n24']
    group_error = lambda left, right: np.sum([np.abs(np.asarray(left['group_formal_order_gradient_changes'][g]).sum(axis=0)
        -np.asarray(right['group_formal_order_gradient_changes'][g]).sum(axis=0)) for g in left['group_formal_order_gradient_changes']], axis=0)
    errors = {'nodes16_24_summed_absolute_group_change': group_error(summaries['n16'], fine),
              'map_scale_summed_absolute_group_change': group_error(summaries['scale07'], fine)}
    maxima = {'quadrature_replay': replay_error,
              **{key: float(np.max(value)) for key, value in errors.items()},
              **{key: max(s['maxima'][key] for s in summaries.values()) for key in fine['maxima']},
              'local_probe_numerical': meta['local_probe']['maximum_numerical_indicator'],
              'local_probe_Euler_identity': meta['local_probe']['Euler_bulk_identity_residual']}
    matter = baseline_matter_source(source[NAMES[0]], source[NAMES[1]], source[NAMES[2]])
    zero_reference = {**fine, 'reference_action_gradient_change': np.zeros(4), 'changed_normal_entries': []}
    baseline = joint_constraint_approximant(incoming_cauchy_jets(), matter,
        source['nsc-incoming-local-constraints']['baseline'], zero_reference)
    changed = joint_constraint_approximant(domain, matter, meta['local_probe'], fine)
    panel_indicator = np.abs(source_action_gradient(incoming_cauchy_jets(),
        source['nsc-incoming-spectral-source']['sum_absolute_group_quadrature_changes'])[:2])
    failures = {k: v for k, v in maxima.items() if v > 3e-11}
    return {'schema': 'NSC-INCOMING-JOINT-CONSTRAINT-ASSEMBLY-v1', 'accountable_author': 'Douglas Ek',
            'status': ('FAIL' if failures else 'PASS')+': incoming reference response and joint bulk assembly; physical constraints OPEN',
            'source_hashes': {p: sha(p) for p in SOURCES},
            'input_hashes': {p: d for p, d in signature().items() if p not in SOURCES},
            'payload': {'path': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()},
            'locked_inputs': source['nsc-incoming-local-constraints']['locked_inputs'],
            'reference_controls': summaries, 'local_probe': meta['local_probe'], 'matter': matter,
            'baseline': baseline, 'mixed_normal_probe': changed,
            'maxima': maxima, 'quadrature_indicators': errors, 'numerical_tolerance': 3e-11,
            'inherited_finite_panel_action_indicator_not_error_bound': panel_indicator,
            'failures': failures, 'physical_source_error_bound': None,
            'scope': {'physical_IV_selected': False, 'constraint_root_claimed': False,
                      'metric_evolution': False, 'old_mode_generators_rerun': False,
                      'Gamma_rest_assigned': False, 'LLL_geometric_count': 1,
                      'physical_endpoint_completion': False, 'push_or_PDF': False},
            'reproducer': 'python3 scripts/derive_nsc_incoming_joint_constraints.py --check',
            'comparison': {'fields': 'all', 'float_atol': 3e-13, 'float_rtol': 3e-13}}


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.prepare:
        if OUTPUT.exists(): raise FileExistsError('existing record is not overwritten')
        path = prepare()
    else:
        old = json.loads(OUTPUT.read_text())
        path = ROOT/old['payload']['path']
        if sha(old['payload']['path']) != old['payload']['sha256']: raise ValueError('reference artifact changed')
    result = json.loads(json.dumps(replay(path), default=jsonable))
    if args.prepare:
        OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    else:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(old, result)
    print(json.dumps({k: result[k] for k in ('status', 'maxima', 'baseline', 'mixed_normal_probe')}, indent=2))
    if result['failures']: raise SystemExit(1)


if __name__ == '__main__': main()
