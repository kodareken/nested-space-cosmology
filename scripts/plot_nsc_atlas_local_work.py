"""Plot retained local constraint data and independent bulk verification."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
ROOT=Path(__file__).resolve().parents[1]
r=json.loads((ROOT/'results/summaries/nsc-atlas-local-work-v1.json').read_bytes())
w=json.loads((ROOT/'results/summaries/nsc-atlas-boundary-work-v1.json').read_bytes())
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,3,figsize=(10,3.8),constrained_layout=True)
source=np.asarray(r['states'][-1]['source_slice_z0'])
scale=max(abs(source.min()),abs(source.max()))
size=r['states'][-1]['padded_size'];stride=max(1,size//24);spacing=4*stride/size
extent=(-2-spacing/2,-2+(source.shape[0]-.5)*spacing,-2-spacing/2,-2+(source.shape[1]-.5)*spacing)
im=axes[0].imshow(source.T,origin='lower',extent=extent,cmap='RdBu_r',norm=TwoSlopeNorm(0,-scale,scale))
axes[0].set(title='A  Local fourth-order source\nFinite vortex at t = .05',xlabel='Slow x',ylabel='Slow y')
fig.colorbar(im,ax=axes[0],shrink=.75,label='S(x,y,0); downsampled stored slice')
h=w['cases'][-1]['history']
axes[1].plot([v['t'] for v in h],[v['db2_dslowtime'] for v in h],color='#245B9C',label='Earlier energy/work identity')
axes[1].scatter([v['t'] for v in r['states']],[v['affine_rate'] for v in r['states']],facecolors='none',edgecolors='#C27517',s=60,label='Independent local source mean / 3',zorder=5)
axes[1].set(title='B  Local and average\nexpansion rates agree',xlabel='Slow time',ylabel='Required affine expansion coefficient')
axes[1].legend(fontsize=7,loc='lower left');axes[1].grid(alpha=.15)
good=max(max(v['maximum_Ricci_coefficient_by_order']) for v in r['direct_bulk_tensor_controls'])
bad=max(v['omitted_expansion_order4_residual'] for v in r['direct_bulk_tensor_controls'])
axes[2].bar(['Completed\nlocal metric','Expansion\nomitted'],[good,bad],color=['#245B9C','#C27517'],width=.55)
axes[2].set(yscale='log',ylim=(1e-16,2),title='C  Independent bulk check\nAnalytic shear-flow control',ylabel='Maximum Ricci series coefficient')
axes[2].text(.03,.15,'Exact radial check: all 125 entries zero.\nThis control is separate from the vortex.',transform=axes[2].transAxes,fontsize=7.8)
axes[2].grid(axis='y',alpha=.15)
fig.suptitle('From a scalar average to a local constraint and a regular bulk witness',fontsize=12)
fig.savefig(ROOT/'paper/figures/nsc-atlas-local-work.png',dpi=170,metadata={'Software':'NSC local work'});plt.close(fig)
print('Wrote local-work figure.')
