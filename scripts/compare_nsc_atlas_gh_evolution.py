"""Compare compact nonlinear-parent records at one proper observer endpoint.

This consumes run summaries only. Spatial changes are empirical indicators,
not certified error bounds; wrapped phase differences never certify winding.
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

from record_nsc_atlas_gh_campaign import source_binding

ROOT = Path(__file__).resolve().parents[1]
OBSERVABLES = (
    'geodesic_observer_total_density_over_initial_vacuum',
    'observer_volume_mode',
    'slice_proper_volume_over_initial_round_reference',
    'g_KK_over_initial_local_value',
    'g_KK_fractional_change_from_initial_local_value',
)
ERRORS = (
    'ADM_H_RMS_over_H_squared', 'ADM_M_RMS_over_H_squared',
    'GH_C_RMS_over_H', 'scalar_last_two_shell_power_fraction',
    'volume_last_two_shell_power_fraction',
    'scalar_acceleration_unresolved_fraction',
    'volume_acceleration_unresolved_fraction',
)


def principal_difference(a, b):
    """The final-minus-initial angle in [-pi, pi], with no winding claim."""
    return math.atan2(math.sin(a-b), math.cos(a-b))


def _tolerances(record):
    return tuple(record['temporal_tolerances'][name]
                 for name in ('rtol', 'atol', 'max_step'))


def _family(record):
    names = ('epsilon', 'prescribed_GH_gauge', 'gauge_constraint_damping',
             'lapse_feedback_coefficient', 'formulation', 'symmetry', 'source',
             'requested_observer_proper_endpoint', 'model_inputs')
    return json.dumps({name: record[name] for name in names}, sort_keys=True)


def summarize(record, *, input_sha256=None, binding_verified=False):
    first, final = record['snapshots'][0], record['snapshots'][-1]
    endpoint = record['requested_observer_proper_endpoint']
    if record['status'] != 'endpoint_reached' or not math.isclose(
            final['observer_proper_duration'], endpoint, rel_tol=0, abs_tol=1e-9):
        raise ValueError('comparison requires the same reached proper endpoint')
    if not final['metric_valid']:
        raise ValueError('comparison requires a valid endpoint metric')
    delta_s = final['observer_volume_mode']-first['observer_volume_mode']
    local = {name: final[name] for name in OBSERVABLES if name in final}
    local.update(
        observer_phase_wrapped=final['observer_phase_wrapped'],
        g_KK_over_initial_local_value=math.exp(-7*delta_s),
        g_KK_fractional_change_from_initial_local_value=math.expm1(-7*delta_s))
    initial_energy = first.get('excitation_normal_energy_over_initial_round_vacuum_energy')
    residual = final.get('excitation_work_balance_residual')
    work = dict(
        enabled=record.get('passive_matter_work', {}).get('enabled', False),
        valid=record.get('passive_matter_work', {}).get('valid'),
        initial_excitation_normalized=initial_energy,
        final_excitation_normalized=final.get('excitation_normal_energy_over_initial_round_vacuum_energy'),
        integrated_work_normalized=final.get('integrated_excitation_work_normalized'),
        residual_normalized=residual,
        residual_over_initial_excitation=(residual/initial_energy
            if residual is not None and initial_energy not in (None, 0) else None))
    return dict(degree=record['degree'], epsilon=record['epsilon'],
        producer_commit=record.get('producer_commit'), input_sha256=input_sha256,
        source_binding_verified=binding_verified,
        temporal_tolerances=record['temporal_tolerances'],
        proper_endpoint=final['observer_proper_duration'],
        observer_phase_winding_certified=False,
        final=local, initial_errors={name:first[name] for name in ERRORS},
        final_errors={name:final[name] for name in ERRORS}, work=work)


def difference(before, after):
    old, new = before['final'], after['final']
    delta = {name:new[name]-old[name] for name in OBSERVABLES}
    delta['observer_phase_principal_difference'] = principal_difference(
        new['observer_phase_wrapped'], old['observer_phase_wrapped'])
    ratios = {name:(after['final_errors'][name]/before['final_errors'][name]
                   if before['final_errors'][name] != 0 else None)
              for name in ERRORS}
    return dict(from_degree=before['degree'], to_degree=after['degree'],
        from_temporal_tolerances=before['temporal_tolerances'],
        to_temporal_tolerances=after['temporal_tolerances'],
        signed_change=delta,
        absolute_change={name:abs(value) for name,value in delta.items()},
        final_error_ratios=ratios)


def compare(records, *, hashes=None, bindings=None):
    if not records:
        raise ValueError('at least one run is required')
    if len({_family(record) for record in records}) != 1:
        raise ValueError('do not combine different preparations, gauges, model inputs or endpoints')
    keys = [(record['degree'], _tolerances(record)) for record in records]
    if len(set(keys)) != len(keys):
        raise ValueError('duplicate degree and time-control settings')
    hashes = [None]*len(records) if hashes is None else hashes
    bindings = [False]*len(records) if bindings is None else bindings
    if len(hashes) != len(records) or len(bindings) != len(records):
        raise ValueError('each record needs one hash and binding flag')
    cases = [summarize(record, input_sha256=digest, binding_verified=bound)
             for record,digest,bound in zip(records, hashes, bindings)]
    paired = sorted(zip(records, cases), key=lambda item:(item[0]['degree'], _tolerances(item[0])))
    spatial = []
    for tolerance in sorted({_tolerances(record) for record in records}):
        group = [case for record,case in paired if _tolerances(record) == tolerance]
        spatial.extend(difference(a,b) for a,b in zip(group, group[1:]))
    temporal = []
    for (ra,a),(rb,b) in itertools.combinations(paired, 2):
        ta,tb = _tolerances(ra),_tolerances(rb)
        if ra['degree'] == rb['degree'] and ta != tb:
            if all(x >= y for x,y in zip(ta,tb)):
                temporal.append(difference(a,b))
            elif all(y >= x for x,y in zip(ta,tb)):
                temporal.append(difference(b,a))
            else:
                raise ValueError('time controls must be consistently tightened')
    return dict(schema='nsc-atlas-gh-comparison-v1',
        epsilon=records[0]['epsilon'], model_inputs=records[0]['model_inputs'],
        requested_proper_endpoint=records[0]['requested_observer_proper_endpoint'],
        cases=[case for _,case in paired], spatial_changes=spatial,
        temporal_changes=temporal, physical_prediction_admitted=False,
        qualification='Adjacent spatial changes are empirical indicators, not certified bounds. A time control at one degree does not certify another degree. Constraint and observable convergence are both required.',
        phase_scope='Direct phase is circular. Principal endpoint differences do not certify accumulated winding.',
        coupling_scope='For d=12, the retained internal action coefficient gives g_KK(s)/g_KK(s_initial)=exp[-7(s-s_initial)]; no Standard Model identification or absolute coupling is selected.',
        slice_scope='Slice volume and normal integrated energy are slicing diagnostics, distinct from local geodesic density.',
        source_sha256={'scripts/compare_nsc_atlas_gh_evolution.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', action='append', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = (ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():
        raise ValueError('fresh repository output required')
    records, hashes = [], []
    for name in args.input:
        path = (ROOT/name).resolve()
        if not path.is_relative_to(ROOT):
            raise ValueError('repository input required')
        raw = path.read_bytes()
        records.append(json.loads(raw)); hashes.append(hashlib.sha256(raw).hexdigest())
    result = compare(records, hashes=hashes,
                     bindings=[source_binding(record) for record in records])
    output.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n',
                      encoding='utf-8', newline='\n')
    print(json.dumps(dict(cases=len(records), spatial_comparisons=len(result['spatial_changes']),
                         temporal_comparisons=len(result['temporal_changes']),
                         physical_prediction_admitted=False)))


if __name__ == '__main__':
    main()
