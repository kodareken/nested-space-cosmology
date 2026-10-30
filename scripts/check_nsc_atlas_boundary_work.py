"""Compute the averaged work closure from the saved finite flow; no new grid."""
import argparse,hashlib,json
from pathlib import Path
from recursive_horizons.nsc_atlas_boundary_work import calculation
ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    out=(ROOT/p.parse_args().output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('fresh repository output required')
    source=ROOT/'results/summaries/nsc-atlas-flow-evolution-v1.json'
    result=calculation(json.loads(source.read_bytes()))
    result['input_sha256']={source.relative_to(ROOT).as_posix():hashlib.sha256(source.read_bytes()).hexdigest()}
    names=('src/recursive_horizons/nsc_atlas_boundary_work.py','scripts/check_nsc_atlas_boundary_work.py')
    result['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(cases=[{k:r[k] for k in ('n','dt','b2','db2_dslowtime','b2_from_integrated_loss')} for r in result['cases']],refinement=result['refinement']),indent=2))


if __name__=='__main__':main()
