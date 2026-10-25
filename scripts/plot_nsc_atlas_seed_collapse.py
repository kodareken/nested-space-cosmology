"""Contrast, fixed-mass contraction and finite source controls."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def render(record):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(1,3,figsize=(12,3.8));curves=record['curves']
    for name,label,color in (('fine','10% seed','#176d8a'),('weak_fine','5% seed','#a56338'),('no_gravity','10%, no self-gravity','#777777')):
        c=curves[name];A=np.array(c['A']);delta=np.array(c['contrast'])
        ax[0].loglog(A,1+delta,color=color,label=label)
        shell=np.array(c['shell_radii'])[:,2]
        ax[1].loglog(A,shell/record['runs'][name]['config']['length'],color=color,label=label)
    ax[0].set(title='A large contrast is not enough',xlabel='Scale factor A',ylabel='Central density / ambient')
    ax[0].legend(frameon=False,fontsize=7)
    ax[1].set(title='The same fixed mass can contract',xlabel='Scale factor A',ylabel='Physical shell radius / supplied L')
    c=curves['fine'];ax[2].semilogy(c['A'],np.array(c['half_height_radius'])/record['runs']['fine']['config']['length'],color='#176d8a')
    ax[2].set(xlim=(.3,1.2),title='Finite concentration and rebound',xlabel='Scale factor A',ylabel='Physical half-height radius / L')
    for a in ax:a.grid(alpha=.15)
    fig.tight_layout();path=Path(__file__).resolve().parents[1]/'docs/figures/nsc-atlas-seed-collapse.png'
    fig.savefig(path,dpi=180);plt.close(fig);print(path)
