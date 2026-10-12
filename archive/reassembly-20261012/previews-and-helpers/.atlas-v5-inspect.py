from pathlib import Path
import json

repo = Path(__file__).parent / 'BlackHoles-Infinity'
preview = Path(__file__).parent / '.atlas-v5-preview-20261010'
for name in ['lab/scripts/derive_nsc_cosmological_atlas_v5.py', 'lab/src/recursive_horizons/nsc_atlas_causal_reference.py']:
    text = (repo / name).read_text(encoding='utf-8')
    bad = [(i, ord(c)) for i, c in enumerate(text) if ord(c) < 32 and c not in '\r\n\t']
    print(json.dumps({'file':name, 'bad_controls':bad, 'arrow_lines': [l for l in text.splitlines() if 'From' in l or "row['from']" in l]}, ensure_ascii=True))
record = json.loads((preview/'causal-map.json').read_bytes())['reference']
print(json.dumps({'classification': record['classification'], 'events':record['events'], 'samples_keys':list(record['wall_samples']), 'first':{k:v[0] for k,v in record['wall_samples'].items() if isinstance(v,list)}}, ensure_ascii=True))
for name, needle in [('docs/arxiv-studies.md','### First sourced interior'), ('lab/docs/nsc-atlas-causal-classification.md','| Directed question')]:
    text = (repo/name).read_text(encoding='utf-8')
    print(json.dumps({'file':name, 'portion':text[text.index(needle):text.index(needle)+1700]},ensure_ascii=True))
