"""Prevent the known duplicate metadata key in the two running v4 summaries.

Only diagnostic serialization is changed. No field, integration, tolerance,
budget, or accumulated numeric diagnostic is changed or replayed.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import runpy
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
V4_SOURCE = 'ae2675980de38f547efa6b4fb1afbd2d12eee62e'


def wrap_transport(transport, report, provenance):
    if 'summary_serialization_patch' in report:
        raise ValueError('summary patch already applied')
    if transport.last is None or transport.last.get('state_corrected') is not False:
        raise ValueError('unexpected legacy diagnostic')
    original = transport.advance
    def advance(*args, **kwargs):
        result = original(*args, **kwargs)
        transport.last.pop('state_corrected')
        return result
    transport.advance = advance
    # The legacy summary constructor supplies this same False explicitly.
    transport.last.pop('state_corrected')
    report['summary_serialization_patch'] = dict(provenance,
        change='remove duplicated state_corrected=False from diagnostic expansion',
        numerical_state_changed=False)


def apply(pid, level, provenance):
    if os.getpid() != pid:
        raise ValueError('wrong process')
    frame = sys._getframe()
    while frame is not None:
        if (frame.f_code.co_name == 'main' and
            Path(frame.f_code.co_filename).resolve() == ROOT/'scripts/run_nsc_atlas_fixed_endpoint.py'):
            break
        frame = frame.f_back
    if frame is None:
        raise ValueError('no live endpoint main frame')
    scope = frame.f_locals
    output = ROOT/f'results/summaries/nsc-atlas-fixed-endpoint-level{level}-v4.json'
    report = scope['report']
    if (scope['args'].level != level or scope['args'].output.resolve() != output or output.exists()
        or report['source_commit'] != V4_SOURCE or report['status'] != 'incomplete'):
        raise ValueError('not the expected unfinished v4 run')
    wrap_transport(scope['transport'], report, provenance)
    print(json.dumps(dict(event='summary_serialization_repaired', level=level,
                         **report['summary_serialization_patch'])), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pid', type=int, required=True)
    parser.add_argument('--level', type=int, choices=(1, 2), required=True)
    args = parser.parse_args()
    source = Path(__file__).resolve()
    sha = hashlib.sha256(source.read_bytes()).hexdigest()
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    pinned = subprocess.check_output(['git', 'show', commit+':scripts/repair_nsc_atlas_live_summary.py'], cwd=ROOT)
    if hashlib.sha256(pinned).hexdigest() != sha:
        raise ValueError('commit the summary repair before applying it')
    provenance = dict(controller_commit=commit, controller_sha256=sha)
    payload = ROOT/f'.build/p001-summary-repair-level{args.level}-pid{args.pid}.py'
    body = ('import hashlib, runpy\n'
            f'p={str(source)!r}\n'
            f'assert hashlib.sha256(open(p,"rb").read()).hexdigest()=={sha!r}\n'
            f'runpy.run_path(p)["apply"]({args.pid!r}, {args.level!r}, {provenance!r})\n')
    with payload.open('x', encoding='utf-8', newline='\n') as stream:
        stream.write(body)
    sys.remote_exec(args.pid, str(payload))
    print('Summary repair dispatched; verify the target event.')
