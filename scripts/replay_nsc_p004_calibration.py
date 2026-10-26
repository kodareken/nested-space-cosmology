"""Recompute calibration in a temporary output, preserving the frozen forecast."""
import json
from pathlib import Path
import tempfile

import check_nsc_atlas_dynamics_lensing as calculation


def main():
    frozen_path=calculation.FORECAST;raw=frozen_path.read_bytes()
    frozen=json.loads(raw)
    parent=(calculation.ROOT/'.build').resolve();parent.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='p004-replay-',dir=parent) as folder:
        target=Path(folder).resolve()
        if not target.is_relative_to(parent):raise RuntimeError('temporary replay directory escaped workspace')
        calculation.FORECAST=target/'forecast.json'
        try:replayed=calculation.prepare()
        finally:calculation.FORECAST=frozen_path
        gap=abs(replayed['prediction']['value_radians']-frozen['prediction']['value_radians'])
        if gap>max(frozen['prediction']['numerical_envelope'],2e-15):
            raise RuntimeError('calibration replay differs from frozen forecast: '+str(gap))
    if frozen_path.read_bytes()!=raw:raise RuntimeError('frozen forecast changed')
    print('NSC-P004 calibration replay passed; frozen bytes unchanged; angular difference',gap)


if __name__=='__main__':main()
