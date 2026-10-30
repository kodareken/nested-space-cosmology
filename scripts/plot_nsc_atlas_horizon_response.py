"""Render the same-source horizon response and its angular/light relations."""
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from recursive_horizons.nsc_atlas_horizon_response import solve_source,solve_response,initial_area_coefficient

ROOT=Path(__file__).resolve().parents[1]
r=json.loads((ROOT/'results/summaries/nsc-atlas-horizon-response-v1.json').read_bytes())
source=solve_source();response=solve_response(source);p=source.parent
u=np.linspace(0.,1.6,401);P,S=response.sol(u)[:2]
deficit=p.kappa/(4*np.pi)*(P-np.sinh(u)*S)
dynamic_deficit=p.kappa*P/(4*np.pi)
c=initial_area_coefficient(p)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(9,7.8),constrained_layout=True)
ax=fig.add_subplot(221)
ax.plot(u,c-deficit,color='#245B9C',label='Ordinary area')
ax.plot(u,c-dynamic_deficit,color='#A65B1D',label='Dynamical area',ls='--')
ax.axhline(c,color='#778494',ls=':',lw=1)
ax.set(title='A  Depletion increases horizon area',xlabel=r'$H t$',ylabel=r'$\Delta A/(\varepsilon^2 A_{\rm vac})$');ax.grid(alpha=.2);ax.legend(fontsize=8.5,loc='lower right')
ax=fig.add_subplot(222,projection='polar')
Q=np.array(r['angular_transfer']['matrix']);angle=np.linspace(0,2*np.pi,601)
n=np.c_[np.cos(angle),np.sin(angle),np.zeros(len(angle))]
value=3*np.einsum('ni,ij,nj->n',n,Q,n)/np.trace(Q)
ax.plot(angle,value,color='#245B9C',lw=1.6)
ax.set_title('B  Directions of outgoing energy\nEquatorial slice; antipodal maxima',pad=16)
ax.set_rmax(3);ax.set_rticks([1,2,3]);ax.tick_params(labelsize=8)
ax=fig.add_subplot(223)
row=r['area_history'][1]
values=[row['ordinary_area_growth_over_epsilon_squared'],row['dynamical_area_growth_over_epsilon_squared'],row['finite_endpoint_correction_over_epsilon_squared']]
ax.bar([0,1,2],values,color=['#245B9C','#A65B1D','#9EB8D7'],width=.6)
for i,val in enumerate(values):ax.text(i,val+.035,f'{val:.5f}',ha='center',fontsize=9)
ax.set(title='C  Finite-time boundary term matters',xticks=[0,1,2],xticklabels=['Ordinary','Dynamical','Correction'],ylim=(0,1.22),ylabel='Area change coefficient');ax.grid(axis='y',alpha=.2)
ax=fig.add_subplot(224)
old=json.loads((ROOT/'results/summaries/nsc-atlas-parent-circulation-v1.json').read_bytes())
for light,relation in zip(old['light_comparisons'],r['light_area_relations']):
    x=relation['dynamical_area_growth_over_epsilon_squared'];y=light['fractional_delay_over_initial_energy_fraction']*1e6
    ax.scatter(x,y,color='#245B9C',s=35)
    ax.annotate(f"R={light['radius']:g}",(x,y),xytext=(0,9),textcoords='offset points',ha='center',fontsize=8.5)
ax.set(title='D  Same preparation fixes both records',xlabel=r'$\Delta A_{\rm dyn}/(\varepsilon^2 A_{\rm vac})$',ylabel=r'$10^6\Delta\tau/(\varepsilon^2 T_{\rm loop})$',xlim=(0,.95),ylim=(0,26));ax.grid(alpha=.2)
fig.suptitle('One scalar source → horizon flux → area, entropy and directional transfer\nQuadratic response; bounded asymptotically de Sitter future; no SI fit',fontsize=12)
fig.savefig(ROOT/'paper/figures/nsc-atlas-horizon-response.png',dpi=170,metadata={'Software':'NSC same-source horizon response'})
plt.close(fig)
print('Wrote horizon-response figure.')
