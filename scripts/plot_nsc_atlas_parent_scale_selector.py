"""Render the conditional volume potential and regenerate one cheap trajectory."""
from pathlib import Path
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from recursive_horizons.nsc_atlas_parent_scale_selector import VolumeParent, evolve

ROOT = Path(__file__).resolve().parents[1]


def main():
    record = json.loads((ROOT/'results/summaries/nsc-atlas-parent-volume-selector-v1.json').read_text(encoding='utf-8'))
    model = VolumeParent()
    _, solution = evolve(model=model)
    N = np.linspace(0, 4, 300)
    s = solution.sol(N)[1]
    plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False})
    fig, axes = plt.subplots(1,2,figsize=(10.5,4.2))
    x = np.linspace(-.065,.075,400)
    ax = axes[0]
    ax.plot(x,[model.potential(t) for t in x],color='#176d8a',lw=2,label='Derived vacuum powers')
    for point in record['vacuum_extrema']:
        ax.scatter([point['s']],[point['vacuum']],color='#b46132',zorder=3)
        ax.annotate(point['kind'].capitalize(),(point['s'],point['vacuum']),
                    xytext=(5,12 if point['kind']=='minimum' else -21),textcoords='offset points',fontsize=9)
    ax.scatter([0],[9/32],color='#737a82',s=25,label='Matched initial state')
    ax.set(xlabel='Physical volume coordinate s',ylabel='Vacuum density in initial model units',
           title='A positive minimum with an ancestor source',ylim=(.20,.46))
    ax.legend(frameon=False,fontsize=8)
    ax.grid(alpha=.15)
    ax=axes[1]
    ax.plot(N,np.exp(12*s),color='#176d8a',lw=2,label='Internal volume / initial volume')
    ax.plot(N,np.exp(-6*s),color='#b46132',lw=2,label='Charge, tension and frequency / initial')
    ax.axhline(1,color='#737a82',lw=.7)
    ax.set(xlabel='Expansion N = log(a/a0)',ylabel='Ratio',title='Evolved together, with coefficients held fixed')
    ax.legend(frameon=False,fontsize=8,loc='center right')
    ax.grid(alpha=.15)
    fig.suptitle('Conditional homogeneous candidate: d = 12, ancestor flux = 5/2',fontsize=11)
    fig.tight_layout(rect=(0,0,1,.94))
    output=ROOT/'docs/figures/nsc-atlas-parent-volume-selector.png'
    fig.savefig(output,dpi=180)
    plt.close(fig)
    print(output)


if __name__=='__main__': main()
