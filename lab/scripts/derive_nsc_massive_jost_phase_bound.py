#!/usr/bin/env python3
"""Record/replay a small directed original-energy exterior phase pilot."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_evaluation_binding import implementation_hashes
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_massive_jost_modes import outgoing_ratio
from recursive_horizons.nsc_massive_jost_phase_bound import subgap_phase_bound

OUTPUT = 'results/development/nsc-massive-jost-phase-bound-v1.json'
INVENTORY = 'results/development/nsc-ks-source-inventory.json'
OWNERS = (
    'src/recursive_horizons/nsc_massive_jost_phase_bound.py',
    'src/recursive_horizons/nsc_massive_jost_modes.py',
    'src/recursive_horizons/nsc_ks_ball_trajectory.py',
    'scripts/derive_nsc_massive_jost_phase_bound.py',
    'tests/test_nsc_massive_jost_phase_bound.py',
    'docs/nsc-massive-jost-phase-bound.md',
)


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def calculate():
    inventory = json.loads((ROOT/INVENTORY).read_text())
    payload = inventory['payload']
    if digest(payload['path']) != payload['sha256']:
        raise ValueError('original source inventory payload changed')
    cases = []
    with np.load(ROOT/payload['path'], allow_pickle=False) as source, ctx.workprec(192):
        metadata = json.loads(source['metadata_json'].tobytes())
        for sign in (-1, 1):
            panel = f'group14/low16_{sign}'
            info = metadata['panels'][panel]
            mass, magnitude = info['mass'], info['angular_magnitude']
            if (info['group'] != 14 or info['angular_sign'] != sign
                    or mass != np.pi/2 or magnitude != np.sqrt(5.)):
                raise ValueError('pilot requires its original group-14 channel')
            mass_ball = arb(mass).union(arb.pi()/2)
            angular_ball = sign*arb(magnitude).union(arb(5).sqrt())
            for row in (0, 7, 15):
                energy = float(source[panel+'/energies'][row])
                for radius in (60., 120.):
                    value, _, coefficients = outgoing_ratio(
                        energy, mass, sign*magnitude, radius, order=8)
                    angle = float(np.angle(value))
                    proof = subgap_phase_bound(energy, mass_ball, angular_ball,
                                               radius, coefficients, angle)
                    cases.append({
                        'source_panel': panel, 'source_row': row,
                        'energy_hex': energy.hex(), 'radius_hex': radius.hex(),
                        'mass_hex': float(mass).hex(),
                        'angular_hex': float(sign*magnitude).hex(),
                        'ratio_polynomial_hex': [[z.real.hex(), z.imag.hex()]
                                                 for z in map(complex, coefficients)],
                        'initializer_phase_hex': angle.hex(),
                        'proof': proof,
                    })
    return {
        'schema': 'NSC-MASSIVE-JOST-PHASE-BOUND-v1',
        'status': 'ENCLOSED: exterior subgap phase initializer for declared pilot inputs',
        'cases': cases,
        'analytic_parameter_cover': 'pi/2 and signed sqrt(5), together with binary source labels',
        'source_inventory_unchanged': True,
        'source_columns_replaced': False,
        'full_family_coverage': False,
        'radial_transport_error': None, 'horizon_sewing_error': None,
        'physical_rho1_source_error': None,
        'physical_local_gate': 'OPEN',
        'source_hashes': implementation_hashes(ROOT, owners=OWNERS),
        'input_hashes': {p: digest(p) for p in (INVENTORY, payload['path'])},
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = calculate()
    if args.record:
        publish_exclusive_file(ROOT, OUTPUT,
            (json.dumps(result, sort_keys=True, indent=2)+'\n').encode())
    elif json.loads((ROOT/OUTPUT).read_text()) != result:
        raise ValueError('exterior phase pilot differs from its immutable record')
    print(json.dumps({
        'status': result['status'], 'cases': len(result['cases']),
        'display_max_phase_error': max(float(restored_upper(
            c['proof']['combined_phase_initial_error_upper'])) for c in result['cases']),
        'physical_rho1_source_error': None, 'physical_local_gate': 'OPEN',
    }, indent=2))
