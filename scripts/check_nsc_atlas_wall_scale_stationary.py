"""Evaluate the inherited constant-tension stationary-radius test, no evolution."""
import argparse
import hashlib
import json
import math
from pathlib import Path

from recursive_horizons.nsc_atlas_flux_conversion import ScalarFlux
from recursive_horizons.nsc_atlas_shared_scalar_slice import SharedScalarSlice

ROOT=Path(__file__).resolve().parents[1]


def calculation():
    atlas=SharedScalarSlice(scalar=ScalarFlux(omega=10**9));rows=[]
    for i,wall in enumerate(atlas.reference_walls):
        A,D,C=wall.coefficients()
        cube=(D+math.sqrt(D*D+8*C*A*A))/(2*C)
        radius=cube**(1/3)
        slope=4*A*A/radius**5+2*D/radius**2-2*C*radius
        second=-20*A*A/radius**6-4*D/radius**3-2*C
        stationary_identity=-12*A*A/radius**6-6*C
        rows.append(dict(wall=i,stationary_radius_model_units=radius,
            potential_at_stationary_radius=wall.potential(radius),slope=slope,
            second_derivative=second,independent_stationary_identity=stationary_identity,
            original_radius=atlas.wall(i)['radius'],original_velocity=atlas.wall(i)['velocity']))
    return dict(schema='nsc-atlas-wall-scale-stationary-v1',walls=rows,
        scope='Instantaneous frozen-field Kottler representation only; not the coupled scalar-wall future or a new EOS.',
        result='Neither stationary point satisfies V=0; both V_second are negative. No stable equilibrium radius selects an SI scale.',
        source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in (
            'scripts/check_nsc_atlas_wall_scale_stationary.py',
            'src/recursive_horizons/nsc_atlas_flux_transition.py',
            'src/recursive_horizons/nsc_atlas_shared_scalar_slice.py')})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();out=(ROOT/args.output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('new output required')
    result=calculation();out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result))


if __name__=='__main__':main()
