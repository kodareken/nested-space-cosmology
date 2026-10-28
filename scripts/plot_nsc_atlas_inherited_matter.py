"""Illustrate the retained inherited-source theorem; no new evolution."""
from pathlib import Path
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def main():
    f=json.loads((ROOT/'results/summaries/nsc-atlas-inherited-matter-forecast-v1.json').read_text(encoding='utf-8'))
    v=json.loads((ROOT/'results/summaries/nsc-atlas-inherited-matter-validation-v1.json').read_text(encoding='utf-8'))
    relation=json.loads((ROOT/'results/summaries/nsc-atlas-inherited-density-relation-v1.json').read_text(encoding='utf-8'))
    from fractions import Fraction
    first=f['frequency_members'][0]
    L,U=[float(Fraction(first['comoving_energy_factors_exact'][k])) for k in ('lower','upper')]
    blue,orange,gray='#176d8a','#b46132','#737a82'
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(11,7.2))
    A=np.linspace(1,2,160);reference=4/(3*A**3)
    ax=axes[0,0]
    ax.fill_between(A,L*reference,U*reference,color=blue,alpha=.2,label='Proved envelope, omega = 20')
    ax.plot(A,reference,'--',color=blue,label='Inherited matter comparison')
    ax.axhline(1,color=gray,lw=.8,label='Equal excitation and vacuum')
    ax.scatter([2],[(4/3)*v['unaveraged_rows'][0]['final_comoving_energy_ratio']/8],color=orange,zorder=3,
               label='Unaveraged endpoint')
    ax.set(title='One prepared field supplies both components',xlabel='Parent scale A = a/a0',ylabel='Excitation / inherited vacuum')
    ax.legend(fontsize=7.5,frameon=False)
    ax.grid(alpha=.15)

    ax=axes[0,1]
    ratio=np.geomspace(.002,4/3,200)
    coefficient=relation['coefficient_approximate']
    curve=lambda c:1/(1+ratio-(c*ratio**2)**(1/3))
    Ov=curve(coefficient);Ochi=ratio*Ov
    ax.plot(Ochi,Ov,color=blue,label='Inherited high-frequency comparison')
    ax.plot([0,.7],[1,.3],'--',color=gray,lw=1,label='Spatially flat comparison')
    eq=relation['dust_comparison_energy_equality']
    ax.scatter([eq['Omega_chi']],[eq['Omega_vacuum']],color=orange,label='Equal-energy comparator')
    ax.set(title='Vacuum, scalar and curvature share one relation',xlabel='Omega_chi: scalar energy fraction',
           ylabel='Omega_vacuum',xlim=(0,.7),ylim=(.3,1.05))
    ax.text(.04,.08,'Epoch and SI ruler cancel.\nC_Omega = 0.008066215066',transform=ax.transAxes,fontsize=9)
    ax.legend(frameon=False,fontsize=7,loc='upper right')
    ax.grid(alpha=.15)

    ax=axes[1,0]
    lo,hi=map(float,first['Lambda_DL_squared_interval'])
    values=[r['Lambda_DL_squared'] for r in v['unaveraged_rows'][:2]]
    ax.axvspan(lo,hi,color=blue,alpha=.18,label='Frozen analytic interval')
    ax.scatter(values,[.65,.35],color=[orange,blue],s=48)
    for value,height,label in zip(values,(.65,.35),('Proper time','Expansion clock')):
        ax.annotate(label,(value,height),xytext=(10,0),textcoords='offset points',va='center',fontsize=9)
    ax.set(title='The unaveraged scalar/gravity check fits the bound',xlabel='Lambda_v D_L^2, A observer = 2, z = 1/4',
           ylim=(0,1),yticks=[],xlim=(lo-.001,hi+.001))
    ax.text(.04,.05,'Q = 0.25169929348\nClock comparison: 2.05e-12',transform=ax.transAxes,fontsize=9)
    ax.legend(frameon=False,fontsize=8,loc='upper left')

    ax=axes[1,1]
    ref=sum(map(float,f['dust_comparison']['Lambda_DL_squared_interval']))/2
    lo,hi=map(float,f['frequency_members'][-1]['Lambda_DL_squared_interval'])
    ax.hlines(.5,(lo-ref)*1e10,(hi-ref)*1e10,color=blue,lw=7)
    ax.scatter([0],[.5],color=orange,s=45,zorder=3,label='Rational dust comparator')
    ax.set(title='High-frequency light-cone enclosure',
           xlabel='(Q - Q_d) x 10^10, omega = 10^9',ylim=(0,1),yticks=[],xlim=(-2.3,2.3))
    ax.text(.5,.83,'Q in [0.251803608594062,\n0.251803608946144]',ha='center',transform=ax.transAxes,fontsize=9)
    ax.text(.5,.12,'Analytic bound on the inherited parent.',ha='center',transform=ax.transAxes,fontsize=8)
    ax.legend(frameon=False,fontsize=8,loc='lower center')
    fig.tight_layout(h_pad=2.8,w_pad=2.5)
    path=ROOT/'docs/figures/nsc-atlas-inherited-matter.png'
    fig.savefig(path,dpi=180);plt.close(fig)
    print(path)


if __name__=='__main__':main()
