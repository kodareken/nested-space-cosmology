"""Record the energy/frequency/radius identities using retained endpoint data."""
import argparse
import hashlib
import json
from pathlib import Path
from recursive_horizons.nsc_atlas_energy_frequency import calculation
ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    out=(ROOT/p.parse_args().output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('new repository output required')
    inputs=('results/summaries/nsc-atlas-occupied-wall-v1.json',
            'results/summaries/nsc-atlas-parent-volume-selector-v1.json')
    records=[json.loads((ROOT/n).read_text(encoding='utf-8')) for n in inputs]
    r=calculation(*records)
    r['input_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in inputs}
    sources=('src/recursive_horizons/nsc_atlas_energy_frequency.py',
             'src/recursive_horizons/nsc_atlas_occupied_wall.py',
             'scripts/check_nsc_atlas_energy_frequency.py')
    r['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in sources}
    out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print('Retained gravitating endpoints:',r['retained_gravitating_endpoint_comparison'])
    print('Thermal half-radius examples:',r['thermal_benchmark']['half_radius_half_energy_T_ratio'],r['thermal_benchmark']['half_radius_one_32_energy_T_ratio'])
    print('Retained volume endpoint:',r['retained_volume_example'])


if __name__=='__main__':main()
