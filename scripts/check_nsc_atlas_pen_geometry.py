"""Compute the declared pen geometry with finite area enclosures."""
import argparse,hashlib,json
from pathlib import Path
from recursive_horizons.nsc_atlas_pen_geometry import calculation,action_phase_clock,phase_residual
ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    out=(ROOT/p.parse_args().output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('new repository output required')
    names=('results/development/scale-binding.json','results/summaries/nsc-atlas-parke-bubble-v1.json',
           'results/summaries/nsc-atlas-parent-volume-selector-v1.json','results/summaries/nsc-atlas-parent-source-v1.json',
           'results/summaries/nsc-atlas-inherited-matter-forecast-v1.json')
    scale,wall,volume,parent,forecast=[json.loads((ROOT/n).read_text(encoding='utf-8')) for n in names]
    assert volume['candidate']['d']==12 and volume['candidate']['ancestor_steps']==2
    minimum=next(x for x in volume['vacuum_extrema'] if x['kind']=='minimum')
    clock=action_phase_clock(minimum,parent['inputs']['omega'])
    r=calculation(scale['development_branch']['omega'],wall['bubble']['radius'],clock=clock)
    r['retained_frequency_member_controls']=[]
    for member in forecast['frequency_members']:
        c=action_phase_clock(minimum,member['omega'])
        c.update(residual_100=phase_residual(100,c['delta']),residual_200=phase_residual(200,c['delta']),
                 log_expansion_100=100*c['log_expansion_per_scalar_cycle'])
        r['retained_frequency_member_controls'].append(c)
    auxiliary=calculation(scale['development_branch']['omega'],wall['bubble']['radius'])
    r['auxiliary_phase_control']={key:auxiliary[key] for key in ('construction','target_comparisons','persistent_coverage')}
    r['auxiliary_phase_control']['probe_ratio_2_over_1']=auxiliary['probe']['ratio_depth_2_over_1']
    r['input_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    sources=('src/recursive_horizons/nsc_atlas_pen_geometry.py','scripts/check_nsc_atlas_pen_geometry.py')
    r['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in sources}
    out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print('delta:',r['construction']['delta'])
    for row in r['target_comparisons']:print('pass:',row)
    for row in r['persistent_coverage']:print('persistent:',row)
    print('coupling ratio:',r['probe']['ratio_depth_2_over_1'],r['probe']['ratio_enclosure'])


if __name__=='__main__':main()
