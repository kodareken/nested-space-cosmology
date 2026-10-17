"""Plot the recorded radial light-cone observables; no scientific replay."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = json.loads(args.input.read_text())
    if result['status'] != 'completed_finite_lightcone':
        raise ValueError('a completed finite light-cone result is required')
    figure, axes = plt.subplots(2, 1, figsize=(6.4, 6.2), sharex=True,
                                layout='constrained')
    for direction, color, label in ((-1, '#1b7186', r'$x_e=x_o-\chi$'),
                                     (1, '#c36a25', r'$x_e=x_o+\chi$')):
        points = sorted((row['observables'] for row in result['primary_rows']
                         if row['observables']['direction'] == direction),
                        key=lambda point: point['coordinate_lookback'])
        radar = [point['radar_distance_over_Lref'] for point in points]
        axes[0].plot(radar, [1e3*point['redshift'] for point in points],
                     'o-', color=color, label=label)
        axes[1].plot(radar, [100*(point['luminosity_over_radar']-1) for point in points],
                     'o-', color=color)
    axes[0].axhline(0, color='#555555', linewidth=.8, linestyle='--')
    axes[1].axhline(0, color='#555555', linewidth=.8, linestyle='--')
    axes[0].set_ylabel(r'Redshift $10^3 z$')
    axes[1].set_ylabel(r'$100(D_L/D_{\rm radar}-1)$ [%]')
    axes[1].set_xlabel(r'Radar distance $D_{\rm radar}/L_{\rm ref}$')
    axes[0].legend(frameon=False, loc='best')
    for axis in axes:
        axis.grid(alpha=.2)
        axis.spines[['top', 'right']].set_visible(False)
    figure.suptitle('Radial light and clock distances\nNormal probes at $T_o=0.3$, $x_o=1$',
                    fontsize=12)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=180)
    plt.close(figure)
    print('Wrote', args.output)


if __name__ == '__main__':
    main()
