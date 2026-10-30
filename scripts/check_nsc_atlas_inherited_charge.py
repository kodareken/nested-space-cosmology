"""Retain the compact inherited harmonic/state calculation, without refitting."""
import argparse,hashlib,json
from pathlib import Path
from recursive_horizons.nsc_atlas_inherited_charge import calculation
ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    out=(ROOT/p.parse_args().output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():
        raise ValueError('fresh repository output required')
    r=calculation()
    names=('src/recursive_horizons/nsc_atlas_inherited_charge.py',
           'src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
           'scripts/check_nsc_atlas_inherited_charge.py')
    r['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    out.write_text(json.dumps(r,sort_keys=True,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('Geometric charged sector:',r['sector'])
    print('Thermal comparison:',r['thermal_control'])
    print('Imbalance:',r['imbalance'])


if __name__=='__main__':main()
