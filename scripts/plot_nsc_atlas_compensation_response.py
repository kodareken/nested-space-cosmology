"""One inherited preparation, three spatial views."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def render(record):
    c=record['curves'];row=record['runs']['fine'];K=np.array(c['K'])
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,3,figsize=(12,3.8))
    axes[0].semilogx(c['radius_over_Rd'],c['enclosed_mass_over_core'],color='#176d8a',lw=2)
    axes[0].axhline(0,color='#888888',lw=.7)
    axes[0].axvline(row['core_radius']/row['reservoir_radius'],color='#a56338',ls=':',label='Core cutoff')
    axes[0].set(title='The core borrows from its surroundings',xlabel=r'$R/R_d$',ylabel=r'$\mathfrak{M}(<R)/\mathfrak{M}_c$')
    axes[0].legend(frameon=False,fontsize=8)
    axes[1].loglog(K,np.array(c['geometry_window'])**2,color='#176d8a',lw=2,label='Core + reservoir')
    evolved=record['evolved_spatial_response']['curves']
    axes[1].loglog(evolved['K'],np.array(evolved['window'])**2,color='#926299',lw=1.5,label='Evolved, same initial units')
    axes[1].loglog(K,np.array(c['core_window'])**2,color='#777777',ls='--',label='Core alone')
    mask=K<1
    axes[1].loglog(K[mask],row['low_K_quadratic_coefficient']**2*K[mask]**4,color='#a56338',ls=':',label=r'$K^4$ limit')
    axes[1].axvline(row['half_power_K'],color='#a56338',ls='--',lw=.9)
    axes[1].set(title='A computed compensation turnover',xlabel=r'$K=kR_d$',ylabel='Squared mass window',ylim=(1e-10,2),xlim=(.01,40))
    axes[1].legend(frameon=False,fontsize=8)
    b=np.array(c['projected_radius_over_Rd']);mask=b>=.08
    axes[2].plot(b[mask],np.array(c['projected_surface'])[mask],color='#a56338',label='Projected contrast')
    axes[2].plot(b[mask],np.array(c['projected_mean_minus_local'])[mask],color='#176d8a',label='Mean minus local')
    axes[2].axhline(0,color='#888888',lw=.7)
    axes[2].set(title='An annular spatial signature',xlabel=r'Projected radius / $R_d$',ylabel=r'Surface units $\mathfrak{M}_c/R_d^2$',ylim=(-1.2,6))
    axes[2].legend(frameon=False,fontsize=8)
    for ax in axes:ax.grid(alpha=.15)
    fig.tight_layout();path=Path(__file__).resolve().parents[1]/'docs/figures/nsc-atlas-compensation-response.png'
    fig.savefig(path,dpi=180);plt.close(fig);print(path)
