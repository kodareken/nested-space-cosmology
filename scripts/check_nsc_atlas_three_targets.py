"""Evaluate the three dimensionless targets from frozen atlas inputs."""
import argparse,hashlib,json,math
from pathlib import Path
from recursive_horizons.nsc_atlas_three_targets import calculation
ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    out=(ROOT/p.parse_args().output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('new repository output required')
    names=('results/development/scale-binding.json','results/summaries/nsc-atlas-parent-volume-selector-v1.json',
           'results/summaries/nsc-atlas-parke-bubble-v1.json','results/summaries/nsc-atlas-occupied-wall-v1.json')
    scale,volume,parke,occupied=[json.loads((ROOT/n).read_text(encoding='utf-8')) for n in names]
    minimum=next(x for x in volume['vacuum_extrema'] if x['kind']=='minimum')
    r=calculation(scale['development_branch']['omega'],minimum['vacuum'],occupied['reference']['level'])
    assert math.isclose(r['inputs']['Parke_R_over_Lunit'],parke['bubble']['radius'],rel_tol=2e-14)
    r['input_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    sources=('src/recursive_horizons/nsc_atlas_three_targets.py','scripts/check_nsc_atlas_three_targets.py',
             'THEORY.md','docs/rn-pg-wall.md','src/recursive_horizons/unified_action.py')
    r['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in sources}
    out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print('Conditional band rows:',r['band']['inheritance']['conditional_rows'])
    print('Orientation controls:',r['asymmetry']['symmetric_orientation_value'],r['asymmetry']['alternative_same_operator_population_value'])
    print('Outside controls:',r['split']['volume_controls'])
    print(r['stopping_result'])


if __name__=='__main__':main()
