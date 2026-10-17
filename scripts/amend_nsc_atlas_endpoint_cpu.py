#!/usr/bin/env python3
"""One recorded CPU-only allowance extension for the three running v4 jobs.

Uses CPython's remote_exec on a specifically identified task-owned process.
No field state, counter, numerical tolerance or acceptance guard is changed.
Fresh reproductions use the runner's --extended-cpu-allowance option instead.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT/'scripts/run_nsc_atlas_fixed_endpoint.py'
ADMISSION = ROOT/'docs/prior-art-and-open-claim.md'
V4_SOURCE = 'ae2675980de38f547efa6b4fb1afbd2d12eee62e'
BASE_CPU_SECONDS = (10000., 20000., 40000.)


def extension_record(budget, report, *, accepted, last_x, provenance):
    """Validate and prepare a CPU-only amendment without changing any state."""
    level = report.get('level')
    if type(level) is not int or level not in range(3):
        raise ValueError('unknown endpoint resolution')
    old = BASE_CPU_SECONDS[level]
    if (report.get('schema') != 'nsc-atlas-fixed-endpoint-v1'
            or report.get('status') != 'incomplete'
            or report.get('source_commit') != V4_SOURCE
            or report.get('cpu_allowance_amendment') is not None
            or report['limits']['cpu_seconds'] != old
            or budget.cpu_seconds != old):
        raise ValueError('not an unamended running v4 CPU allowance')
    elapsed = time.process_time()-budget.started
    if not math.isfinite(elapsed) or not 0 <= elapsed < old:
        raise ValueError('the original CPU allowance has already expired')
    if budget.evaluations > budget.maximum_evaluations:
        raise ValueError('charged-value limit has expired')
    record = dict(schema='nsc-atlas-cpu-allowance-amendment-v1', status='applied',
        level=level, process_id=os.getpid(), utc=datetime.now(timezone.utc).isoformat(),
        old_cpu_seconds=old, new_cpu_seconds=2*old, elapsed_cpu_seconds=elapsed,
        charged_values=budget.evaluations, accepted_steps=accepted, last_accepted_x=last_x,
        progress_is_restart_checkpoint=False,
        source_commit=report['source_commit'], provenance=provenance,
        reason='Measured per-step CPU cost grows as additional retained history panels enter the cap; preserve the full-endpoint execution under a finite extended allowance.',
        counters_reset=False, physical_state_changed=False,
        numerical_settings_changed=False,
        unchanged_limits={k:v for k,v in report['limits'].items() if k!='cpu_seconds'},
        radial_accuracy=dict(report['radial_accuracy']))
    # Prepare the complete serializable record before either assignment.
    json.dumps(record, allow_nan=False)
    return record


def extend_allowance(budget, report, *, accepted, last_x, provenance):
    """Double the original ceiling once, without resetting any accounting."""
    record = extension_record(budget, report, accepted=accepted, last_x=last_x,
                              provenance=provenance)
    budget.cpu_seconds = record['new_cpu_seconds']
    report['limits']['cpu_seconds'] = record['new_cpu_seconds']
    report['cpu_allowance_amendment'] = record
    return record


def apply_to_current_process(*, pid, level, receipt, provenance):
    """Validate the interrupted main frame before touching its CPU allowance."""
    if os.getpid() != pid:
        raise ValueError('target process identity changed')
    frame = sys._getframe()
    while frame is not None:
        if (frame.f_code.co_name == 'main'
                and Path(frame.f_code.co_filename).resolve() == RUNNER):
            break
        frame = frame.f_back
    if frame is None:
        raise ValueError('no active endpoint runner main frame')
    scope = frame.f_locals
    expected_output = ROOT/f'results/summaries/nsc-atlas-fixed-endpoint-level{level}-v4.json'
    if (scope['args'].level != level or scope['args'].output.resolve() != expected_output
            or expected_output.exists()):
        raise ValueError('target is not the selected unfinished v4 evolution')
    path = (ROOT/receipt).resolve()
    if path.parent != ROOT/'results/summaries' or path.exists():
        raise ValueError('require a new compact amendment receipt')
    budget, report = scope['budget'], scope['report']
    old = budget.cpu_seconds
    # Reject a wrong or already amended target before reserving any file.
    record = extension_record(budget, report, accepted=scope['accepted'],
                              last_x=scope['last_x'], provenance=provenance)
    # Reserve this receipt before changing the in-memory allowance. A write
    # failure restores the original allowance; it never leaves a silent change.
    reserved = False
    try:
        with path.open('x', encoding='utf-8', newline='\n') as stream:
            reserved = True
            budget.cpu_seconds = record['new_cpu_seconds']
            report['limits']['cpu_seconds'] = record['new_cpu_seconds']
            report['cpu_allowance_amendment'] = record
            stream.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False)+'\n')
            stream.flush()
            os.fsync(stream.fileno())
    except BaseException:
        if reserved:
            budget.cpu_seconds = old
            report['limits']['cpu_seconds'] = old
            report.pop('cpu_allowance_amendment', None)
            # Only the newly created task-owned receipt is removed, after close.
            path.unlink(missing_ok=True)
        raise
    print(json.dumps(dict(event='cpu_allowance_amended', **record)), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pid', type=int, required=True)
    parser.add_argument('--level', type=int, choices=range(3), required=True)
    args = parser.parse_args()
    if args.pid <= 0:
        raise ValueError('positive task-owned process id required')
    source = Path(__file__).resolve()
    relative = source.relative_to(ROOT).as_posix()
    raw = source.read_bytes()
    commit = subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip()
    committed = subprocess.check_output(['git','show',f'{commit}:{relative}'], cwd=ROOT)
    if committed != raw:
        raise ValueError('commit the allowance controller before execution')
    admission = ADMISSION.read_bytes()
    if subprocess.check_output(['git','show',f'{commit}:docs/prior-art-and-open-claim.md'],cwd=ROOT)!=admission:
        raise ValueError('commit the admission before execution')
    provenance = dict(controller_commit=commit, controller_path=relative,
        controller_sha256=hashlib.sha256(raw).hexdigest(),
        admission_path=ADMISSION.relative_to(ROOT).as_posix(),
        admission_sha256=hashlib.sha256(admission).hexdigest())
    receipt = f'results/summaries/nsc-atlas-fixed-endpoint-level{args.level}-v4-cpu-amendment.json'
    if (ROOT/receipt).exists():
        raise ValueError('an amendment receipt already exists')
    payload = ROOT/f'.build/results/p001-cpu-amendment-level{args.level}-pid{args.pid}.py'
    call = dict(pid=args.pid, level=args.level, receipt=receipt, provenance=provenance)
    body = ('import hashlib, runpy\n'
            f'p = {str(source)!r}\n'
            f'assert hashlib.sha256(open(p,"rb").read()).hexdigest() == {provenance["controller_sha256"]!r}\n'
            'controller = runpy.run_path(p, run_name="nsc_cpu_allowance_controller")\n'
            f'controller["apply_to_current_process"](**{call!r})\n')
    with payload.open('x', encoding='utf-8', newline='\n') as stream:
        stream.write(body)
    sys.remote_exec(args.pid, str(payload))
    print(json.dumps(dict(status='dispatched_not_yet_verified', process_id=args.pid,
                          receipt=receipt, payload=str(payload))))


if __name__ == '__main__':
    main()
