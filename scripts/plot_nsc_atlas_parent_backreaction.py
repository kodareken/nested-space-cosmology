"""Regenerate the constrained feedback, clock and matter-work comparison."""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from recursive_horizons.nsc_atlas_parent_backreaction import (
    constraints,energy_and_work,harmonic_grams,observer_record,solve_backreaction)

ROOT=Path(__file__).resolve().parents[1]


def main():
    r=solve_backreaction();u=np.linspace(0,8,801);G=harmonic_grams()
    degrees={l:[] for l in (3,4,5,6)};energies=[];works=[];errors=[]
    E0=energy_and_work(r.volume,r.operators,r.G,r.GP,r.GT,0.,r.solution.y[:,0])[0]
    for v in u:
        _,_,pairs=r.fields(float(v));c=pairs[0][0]+pairs[1][0]
        for l in degrees:degrees[l].append(max(0.,float(c@G[l]@c)*2*np.pi**2))
        y=r.solution.sol(v)
        E,_,_=energy_and_work(r.volume,r.operators,r.G,r.GP,r.GT,float(v),y)
        energies.append(E-E0);works.append(y[54])
        h,m,V=constraints(r,float(v));H=r.source.parent.H
        errors.append([max(np.abs(h))/H**2,max(np.abs(m))/H,abs(V)/H])
    short=np.linspace(0,2,301);clocks=[observer_record(r,float(v)) for v in short]
    metric=np.array([v['metric_feedback']['H_phase_event_proper_shift_over_epsilon_squared'] for v in clocks])
    volume=np.array([v['volume_feedback']['H_phase_event_proper_shift_over_epsilon_squared'] for v in clocks])
    fig,ax=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for l,values in degrees.items():ax[0,0].plot(u,values,label=f'Degree {l}')
    ax[0,0].set(xlim=(0,2.5),xlabel='Offset readout Ht',ylabel=r'Squared field norm / $\epsilon^6$',title='A  The original two modes generate finer harmonics')
    ax[0,0].legend(fontsize=9)
    ax[0,1].plot(short,volume,label='Volume response',color='#b27124')
    ax[0,1].plot(short,metric,label='Metric response',color='#3b8668')
    ax[0,1].plot(short,volume+metric,label='Combined',color='#225b91',ls='--')
    ax[0,1].set(xlabel='Offset readout Ht',ylabel=r'$H\,\delta\tau_{phase}/\epsilon^2$',title='B  The inherited scalar changes its own phase clock')
    ax[0,1].legend(fontsize=9)
    lo=min(min(energies),min(works));hi=max(max(energies),max(works))
    ax[1,0].plot([lo,hi],[lo,hi],ls='--',color='grey',label='Exact work identity')
    ax[1,0].scatter(works,energies,s=6,color='#225b91',label='Evolved finite history')
    ax[1,0].set(xlabel='Integrated matter work',ylabel=r'$E_4(t)-E_4(0)$',title='C  Matter, volume and spatial measure close one balance')
    ax[1,0].legend(fontsize=9)
    errors=np.array(errors)
    for i,label in enumerate(('Hamiltonian / H squared','Scalar momentum / H','Vector momentum / H')):
        ax[1,1].semilogy(u,np.maximum(errors[:,i],1e-17),label=label,lw=1)
    ax[1,1].set(xlabel='Offset readout Ht',ylabel='Dimensionless constraint residual',title='D  The evolving metric satisfies all constraint sectors')
    ax[1,1].legend(fontsize=8.5)
    for a in ax.flat:a.grid(alpha=.18)
    fig.suptitle('Inherited wave → metric and volume → scalar response and clock comparison\nComplete stated perturbative orders; finite-amplitude collapse and SI scale remain open',fontsize=13)
    out=ROOT/'paper/figures/nsc-atlas-parent-backreaction.png'
    fig.savefig(out,dpi=180);plt.close(fig)
    print(out.relative_to(ROOT))


if __name__=='__main__':main()
