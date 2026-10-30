"""Plot the smooth-parent field, compensated response and finite light record."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from recursive_horizons.nsc_atlas_parent_circulation import Parent,evolve,static_patch_energy

ROOT=Path(__file__).resolve().parents[1]
record=json.loads((ROOT/'results/summaries/nsc-atlas-parent-circulation-v1.json').read_bytes())
p=Parent.retained();sol=evolve(parent=p)
theta=np.linspace(0.,sol.t[-1],2401);z=sol.sol(theta);n=theta/(2*np.pi*p.offset)
plt.rcParams.update({'font.size':10.5,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,2,figsize=(9,7.5),constrained_layout=True)
ax=axes[0,0];ax.plot(n,z[0]*p.H,color='#245B9C',lw=1.8)
ax.set(title='A  Smooth closed parent expands',xlabel='Reference phase events',ylabel='Spatial radius / initial radius');ax.grid(alpha=.2)
ax=axes[0,1];ax.plot(n,p.kappa*z[6]/(5*p.H),color='#A65B1D',lw=1.3)
ax.set(title='B  Regional circulation develops\nAll six total rotation charges remain zero',xlabel='Reference phase events',ylabel=r'$\beta/H$ per initial energy fraction');ax.grid(alpha=.2)
ax=axes[1,0];sel=np.linspace(0,len(theta)-1,151,dtype=int)
E0=record['horizon_exchange']['initial_Killing_energy_over_epsilon_squared']
remaining=np.array([static_patch_energy(z[8,i],z[:,i],p)/E0 for i in sel])
ax.plot(n[sel],remaining,color='#245B9C',lw=1.8,label='Retained in observer patch')
ax.plot(n[sel],1-remaining,color='#A65B1D',lw=1.8,label='Crossed smooth horizon')
ax.set(title='C  Same field crosses a smooth horizon\nKilling-energy balance closes',xlabel='Reference phase events',ylabel='Fraction of initial patch energy',ylim=(0,1));ax.grid(alpha=.2);ax.legend(fontsize=8.5)
ax=axes[1,1]
rows=list(reversed(record['light_comparisons']));values=[r['H_delay_over_epsilon_squared']*1e6 for r in rows]
ax.bar([0,1],values,color=['#9EB8D7','#245B9C'],width=.55)
for i,val in enumerate(values):ax.text(i,val+.3,f'{val:.5f}',ha='center',fontsize=10)
ax.set(title='D  Differential guided-light record\nTwo equal-radius rings; phase averaged',xticks=[0,1],xticklabels=['R = 1/16','R = 1/8'],ylabel=r'$10^6 H\,\Delta\tau/\varepsilon^2$',ylim=(0,21));ax.grid(axis='y',alpha=.2)
fig.suptitle('Generated spatial geometry → inherited scalar → vector response → light record\nConditional two-harmonic preparation; leading quadratic response; no SI ruler',fontsize=12)
fig.savefig(ROOT/'paper/figures/nsc-atlas-parent-circulation.png',dpi=170,metadata={'Software':'NSC smooth-parent circulation'})
plt.close(fig)
print('Wrote smooth-parent circulation figure.')
