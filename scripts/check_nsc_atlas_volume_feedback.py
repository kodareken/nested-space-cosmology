"""Compute the same-source causal volume feedback; save only a compact record."""
import argparse
import hashlib
import json
from pathlib import Path

from recursive_horizons.nsc_atlas_volume_feedback import calculation

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',required=True,type=Path)
    out=(ROOT/p.parse_args().output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('fresh repository output required')
    r=calculation()
    names=['src/recursive_horizons/nsc_atlas_volume_feedback.py',
        'src/recursive_horizons/nsc_atlas_horizon_response.py',
        'src/recursive_horizons/nsc_atlas_parent_circulation.py',
        'src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
        'src/recursive_horizons/nsc_atlas_inherited_charge.py',
        'src/recursive_horizons/nsc_atlas_s3_polynomials.py',
        'scripts/check_nsc_atlas_volume_feedback.py']
    r['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    event=r['samples'][4]
    print(json.dumps(dict(reference_event=event,
        numerical_control=r['numerical_control'],causal_memory=r['causal_memory'],
        Schur_interaction=r['Schur_interaction'],elapsed_seconds=r['elapsed_seconds']),indent=2))


if __name__=='__main__':main()
