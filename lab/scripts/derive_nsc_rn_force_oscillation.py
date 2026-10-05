"""Creation-only consumer of the recorded RN family."""
from pathlib import Path
import json
from recursive_horizons.nsc_rn_force_oscillation import analyze_family


def main():
    root = Path(__file__).resolve().parents[1]
    output = root / 'results/development/nsc-rn-force-oscillation-v1.json'
    if output.exists():
        raise FileExistsError(f'Creation-only record already exists: {output}')
    record = analyze_family(root / 'results/development/nsc-rn-neutral-family-v1')
    assert all(v == 'CPUDispatcher' for v in record['numba_dispatchers'].values())
    with output.open('x', encoding='utf-8') as stream:
        json.dump(record, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(output)
    for case in record['cases']:
        print(case['case'], 'stations:', len(case['stations']))
        for station in case['stations']:
            fields = station['profiles']
            print(station['station_label'], station['points'],
                  [(name, fields[name]['sign_changes'], fields[name]['mean_cells'])
                   for name in ('radial_pressure', 'angular_pressure')])


if __name__ == '__main__':
    main()
