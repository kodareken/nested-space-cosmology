from pathlib import Path
import json
import subprocess

root = Path(__file__).resolve().parent / 'BlackHoles-Infinity'
paths = subprocess.check_output(['git', 'diff', '--cached', '--name-only'], cwd=root, text=True).splitlines()
forbidden = [p for p in paths if p.startswith('.pi/') or p.endswith('.pdf') or '/nsc-cosmological-atlas-v1/' in p or '/nsc-cosmological-atlas-v2/' in p]
assert not forbidden, forbidden
before = subprocess.check_output(['git', 'show', 'HEAD:docs/prior-art-and-open-claim.md'], cwd=root).decode('utf-8')
after = (root / 'docs/prior-art-and-open-claim.md').read_text(encoding='utf-8')
start, stop = '## Imported results', '## What is established'
assert before[before.index(start):before.index(stop)] == after[after.index(start):after.index(stop)]
staged_byte_mismatches = []
for p in paths:
    staged = subprocess.check_output(['git', 'show', ':' + p], cwd=root)
    if staged != (root / p).read_bytes():
        staged_byte_mismatches.append(p)
assert not staged_byte_mismatches, staged_byte_mismatches
print(json.dumps({'staged_paths': paths, 'count': len(paths), 'imported_results_unchanged': True, 'staged_bytes_match_working_tree': True}, indent=2))
