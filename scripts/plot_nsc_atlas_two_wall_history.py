"""Illustrate the admitted shared host history and its incoming characteristics."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from recursive_horizons.nsc_atlas_nested_spacetime import TwoWallHistoryBudget, TwoWallHostHistory


ROOT = Path(__file__).resolve().parents[1]
SUMMARY = 'results/summaries/nsc-atlas-two-wall-history-v1.json'
OUTPUT = ROOT / 'docs/figures/nsc-atlas-two-wall-history'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    record = json.loads((ROOT / SUMMARY).read_bytes())
    for source, expected in record['source_sha256'].items():
        if digest(ROOT / source) != expected:
            raise ValueError('figure input source differs: '+source)
    budget = TwoWallHistoryBudget(limit=50000, cpu_seconds=30.)
    model = TwoWallHostHistory(order=64, budget=budget)
    times = np.linspace(0., .04, 17)
    blue, green, grey = '#1e6285', '#16856c', '#687985'
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11.5,
                         'svg.hashsalt': 'nsc-two-wall-history-v1'})
    fig, axes = plt.subplots(1, 2, figsize=(9., 4.7))
    fig.subplots_adjust(left=.085, right=.97, top=.76, bottom=.26, wspace=.40)
    fig.suptitle('Two walls in one host: distinct causal access',
                 fontsize=15, fontweight='bold', y=.97, color='#183848')
    fig.text(.5, .895, 'Host-side rays reach each wall from the original slice.',
             ha='center', fontsize=11, color='#526776')
    for ax, name in zip(axes, ('inner', 'outer')):
        saved = record['cases'][name]
        sign = -1 if name == 'inner' else 1
        origin = saved['initial_event']['null_coordinate']
        events = [model.wall_event(name, float(tau)) for tau in times]
        radii = np.array([row['radius'] for row in events])
        clock = sign*(np.array([row['null_coordinate'] for row in events])-origin)
        foot = saved['endpoint_source']['radius']
        initial_radius = saved['initial_event']['radius']
        pad = .18*(initial_radius-foot)
        initial_radii = np.linspace(foot-pad, initial_radius, 17)
        initial_clock = sign*(np.array([model.initial_event(name, float(r))['null_coordinate']
                                       for r in initial_radii])-origin)
        ax.plot(initial_radii, initial_clock, color=grey, lw=2.1,
                label='Original host slice')
        ax.plot(radii, clock, color=blue, lw=2.4, label='Material wall')
        ax.scatter([radii[0], radii[-1]], [clock[0], clock[-1]], s=24, color=blue, zorder=4)
        for index in (8, 16):
            event = events[index]
            source = model.characteristic_origin(name, event)
            y = clock[index]
            ax.plot([source['radius'], event['radius']], [y, y], color=green, lw=1.7,
                    label='Host-side null ray' if index == 8 else None)
            midpoint = (source['radius']+event['radius'])/2
            ax.annotate('', xy=(midpoint+.07*(event['radius']-source['radius']), y),
                        xytext=(midpoint-.07*(event['radius']-source['radius']), y),
                        arrowprops={'arrowstyle': '-|>', 'color': green, 'lw': 1.5})
            ax.scatter([source['radius']], [y], s=22, color=green, zorder=4)
            ax.annotate(r'$\tau/\ell='+f'{times[index]:.2f}'+r'$',
                        (event['radius'], y), xytext=(-5, 7), textcoords='offset points',
                        ha='right', fontsize=10.5, color=blue)
        ax.set_title(('Inner wall: left host patch' if name == 'inner' else
                      'Outer wall: expanding host patch'), fontsize=11.5, pad=12)
        ax.set_xlabel(r'Areal radius $R/\ell$')
        ax.set_ylabel(r'$-(p-p_0)/\ell$' if name == 'inner' else r'$(q-q_0)/\ell$')
        ax.set_xlim(foot-1.4*pad, radii[-1]+pad)
        ax.set_ylim(-.06*clock[-1], 1.18*max(clock[-1], initial_clock.max()))
        ax.grid(alpha=.18)
        ax.spines[['top', 'right']].set_visible(False)
        labels = saved['source_proper_label_interval']
        ax.text(.5, -.39, f'Source labels: {labels[0]:.6f} to {labels[1]:.6f} '+r'$\ell$',
                transform=ax.transAxes, ha='center', fontsize=10, color='#526776')
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(.5, .095),
               ncol=3, frameon=False, fontsize=10.5)
    fig.text(.5, .015, 'The two proper-time endpoints do not define a simultaneity slice.',
             ha='center', fontsize=10, color='#526776')
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ('.png', '.svg'):
        kwargs = {'metadata': {'Software': 'NSC two-wall history'}} if suffix == '.png' else {'metadata': {'Date': None}}
        fig.savefig(OUTPUT.with_suffix(suffix), dpi=180, **kwargs)
    plt.close(fig)
    budget.add(0)
    sources = [SUMMARY, 'src/recursive_horizons/nsc_atlas_nested_spacetime.py',
               'scripts/plot_nsc_atlas_two_wall_history.py']
    result = dict(schema='nsc-atlas-two-wall-history-figure-v1',
        inputs={path: digest(ROOT/path) for path in sources},
        outputs={OUTPUT.with_suffix(suffix).relative_to(ROOT).as_posix(): digest(OUTPUT.with_suffix(suffix))
                 for suffix in ('.png', '.svg')},
        geometry_evaluations=budget.count, wall_samples_per_side=17,
        process_cpu_seconds=time.process_time()-budget.started,
        scope='Regenerated display of the admitted finite geometry; no new field or parameter experiment.')
    OUTPUT.with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
