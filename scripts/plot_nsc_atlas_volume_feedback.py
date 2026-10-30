"""Show causal volume memory and its linked local spectral/gauge response."""
from pathlib import Path
import math

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from recursive_horizons.nsc_atlas_volume_feedback import (
    gram_matrix,sample,solve_feedback,source_coefficients)

ROOT=Path(__file__).resolve().parents[1]


def main():
    r=solve_feedback();G=gram_matrix();u=np.linspace(0,8,3201)
    rows=[sample(r,float(v),G) for v in u]
    mean=np.array([v['spatial_mean_volume_shift_over_epsilon_squared'] for v in rows])
    center=np.array([v['central_volume_shift_over_epsilon_squared'] for v in rows])
    horizon=np.array([v['horizon_mean_volume_shift_over_epsilon_squared'] for v in rows])
    energy=np.array([v['induced_volume_energy_over_epsilon_fourth'] for v in rows])
    static=np.array([r.d*r.source.parent.mass_squared*source_coefficients(r.source,float(v))[0]/(2*r.stiffness) for v in u])
    b=r.solution.sol(u);nu=math.sqrt(r.nu_squared)
    envelope=np.sqrt(b[0]**2+(b[5]/nu)**2)
    fig,ax=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for vals,label,color in ((mean,'Spatial mean','#225b91'),(center,'Central observer','#8b4a1c'),(horizon,'Horizon mean','#3b8668')):
        ax[0,0].plot(u,vals,label=label,color=color)
    ax[0,0].set(xlim=(0,3),xlabel='Offset readout Ht',ylabel=r'$s_2$ (coefficient of $\epsilon^2$)',title='A  The inherited wave drives the internal volume')
    ax[0,0].axhline(0,color='grey',lw=.7);ax[0,0].legend(fontsize=9)
    ax[0,1].plot(-center,-7*center,label='Central response',color='#225b91',lw=2)
    ax[0,1].scatter(-horizon[::50],-7*horizon[::50],s=10,label='Horizon mean',color='#b27124')
    xx=np.array([min(-center.max(),-horizon.max()),max(-center.min(),-horizon.min())])
    ax[0,1].plot(xx,7*xx,color='grey',ls=':',lw=1,zorder=0)
    ax[0,1].set(xlabel=r'$\delta\ln R_M/\epsilon^2$',ylabel=r'$\delta\ln g_{KK}/\epsilon^2$',title='B  One geometry fixes the slope: 7')
    ax[0,1].legend(fontsize=9)
    ax[1,0].semilogy(u,envelope,label='Retained volume quadrature norm',color='#225b91',lw=2)
    ax[1,0].semilogy(u,np.cosh(u)**1.5*static,label='Instantaneous static response, rescaled',color='#b27124',alpha=.55,lw=.7)
    ax[1,0].set(xlabel='Offset readout Ht',ylabel='Rescaled mean-response magnitude',ylim=(1e-8,.1),title='C  Causal memory outlives the local source')
    ax[1,0].legend(fontsize=8.5)
    ax[1,1].plot(u,energy,color='#225b91',label=r'$E_s/\epsilon^4$')
    ax[1,1].set(xlabel='Offset readout Ht',ylabel='Induced volume energy',title='D  Energy is transferred into an initially quiet mode')
    ax[1,1].legend(fontsize=9)
    for a in ax.flat:a.grid(alpha=.18)
    fig.suptitle('The same inherited field changes geometry, local spectra and coupling\nLeading volume response; no new source, SI ruler or present-day fraction selected',fontsize=13)
    out=ROOT/'paper/figures/nsc-atlas-volume-feedback.png'
    fig.savefig(out,dpi=180);plt.close(fig)
    print(out.relative_to(ROOT))


if __name__=='__main__':main()
