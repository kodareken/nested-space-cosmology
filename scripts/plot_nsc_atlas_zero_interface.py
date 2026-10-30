"""A calculated decaying mode, its spiral/sphere representation and spectral tail."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from recursive_horizons.nsc_atlas_horizon_response import solve_source

ROOT=Path(__file__).resolve().parents[1]
r=json.loads((ROOT/'results/summaries/nsc-atlas-zero-interface-v1.json').read_bytes())
source=solve_source();nu=np.sqrt(source.nu_squared)
u=np.linspace(4.,8.,6001);y=source.modes.sol(u);O=1/np.cosh(u)
q=O**1.5*y[0];pquad=O**1.5*y[1]/nu
normalizer=np.hypot(q[0],pquad[0]);X=q/normalizer;Y=pquad/normalizer
turn=(u-u[0])*nu/(2*np.pi)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(9,7.8),constrained_layout=True)
ax=fig.add_subplot(221);mask=turn<=8
ax.axhline(0,color='#122033',lw=1.5)
ax.plot(turn[mask],X[mask],color='#245B9C',lw=1.4)
ax.fill_between(turn[mask],0,X[mask],where=X[mask]>=0,color='#245B9C',alpha=.16)
ax.fill_between(turn[mask],0,X[mask],where=X[mask]<0,color='#BD7B1A',alpha=.18)
envelope=np.cosh(u[0])**1.5*O**1.5
ax.plot(turn[mask],envelope[mask],':',color='#7F8B97',lw=1);ax.plot(turn[mask],-envelope[mask],':',color='#7F8B97',lw=1)
ax.set(title='A  Signed oscillation about zero\nEnergy remains quadratic and positive',xlabel='Late carrier turns',ylabel='Normalized field amplitude');ax.grid(alpha=.15)
ax=fig.add_subplot(222);mask=turn<=30
ax.plot(X[mask],Y[mask],color='#245B9C',lw=.9);ax.scatter([0],[0],s=12,color='#BD7B1A')
ax.set(aspect='equal',title='B  The same mode in phase space\nDamping produces an inward spiral',xlabel='Normalized field',ylabel='Normalized quadrature');ax.grid(alpha=.15)
ax=fig.add_subplot(223,projection='3d');radial=X*X+Y*Y
sx=2*X/(1+radial);sy=2*Y/(1+radial);sz=(1-radial)/(1+radial)
theta=np.linspace(0,np.pi,13);phi=np.linspace(0,2*np.pi,25);T,P=np.meshgrid(theta,phi,indexing='ij')
ax.plot_wireframe(np.sin(T)*np.cos(P),np.sin(T)*np.sin(P),np.cos(T),color='#9AAABD',alpha=.25,lw=.5)
ax.plot(sx,sy,sz,color='#245B9C',lw=1);ax.set_box_aspect((1,1,1),zoom=1.25);ax.view_init(20,35);ax.set_axis_off()
ax.set_title('C  Stereographic sphere representation\nA mode-state map, not a new metric')
ax=fig.add_subplot(224);a=np.cosh(u)/np.cosh(u[0])
ax.loglog(a,a**-3,color='#245B9C',label='Excitation density: a^-3')
ax.loglog(a,a,color='#BD7B1A',label='Spatial wavelength: a')
ax.loglog(a,a**-1,'--',color='#67816C',label='Massless probe frequency: a^-1')
ax.axhline(1,color='#7F8B97',ls=':',label='Massive carrier: finite limit')
ax.set(title='D  Different spectral quantities\nTail powers follow the declared operator',xlabel='Spatial scale / initial scale',ylabel='Normalized measure');ax.grid(which='both',alpha=.15);ax.legend(fontsize=7.7,loc='lower left')
fig.suptitle('Offset-parametrized mode evolution: zero baseline, damping and a geometric spectrum\nThe inherited scalar supplies the decay rate; sphere projection is a representation',fontsize=12)
fig.savefig(ROOT/'paper/figures/nsc-atlas-zero-interface.png',dpi=170,metadata={'Software':'NSC offset thermodynamics'})
plt.close(fig)
print('Wrote zero-interface and mode-thermodynamics figure.')
