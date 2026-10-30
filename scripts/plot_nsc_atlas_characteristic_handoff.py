"""Show the inherited null data, exterior receiver and fixed-grid convergence."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from recursive_horizons.nsc_atlas_horizon_response import solve_source
from recursive_horizons.nsc_atlas_characteristic_handoff import (
    Diamond,evaluate_case,inherited_second_leg,inherited_reference,geometry)

ROOT=Path(__file__).resolve().parents[1]
r=json.loads((ROOT/'results/summaries/nsc-atlas-characteristic-handoff-v1.json').read_bytes())
source=solve_source(endpoint_u=2.);domain=Diamond.inherited(source.parent)
axis=r['domain']['receiver_angular_axis']
_,(u,v,psi,exact)=evaluate_case(source,domain,128,axis,return_arrays=True)
U,V=np.meshgrid(u,v,indexing='ij');eta,chi,a,_=geometry(source,U,V)
field=np.einsum('ijk,k->ij',psi,np.asarray(axis[:2]))/(a*np.sin(chi))
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,2,figsize=(9,7.5),constrained_layout=True)
ax=axes[0,0];scale=float(np.max(abs(field)))
im=ax.pcolormesh(eta,chi,field,cmap='RdBu_r',vmin=-scale,vmax=scale,shading='auto',rasterized=True)
ax.plot(eta[:,0],chi[:,0],color='#122033',lw=2,label='Parent horizon trace')
ax.plot(eta[0],chi[0],color='#BD7B1A',lw=2,label='Inherited second trace')
ax.plot([domain.eta_corner,domain.eta_corner+domain.width],[domain.chi_receiver]*2,'--',color='#245B9C',lw=1.2,label='Receiving worldline')
ax.set(title='A  Smooth continuation\nOutside the parent horizon',xlabel=r'Conformal $\eta$',ylabel=r'Radial $\chi$')
ax.legend(fontsize=7,loc='upper right');fig.colorbar(im,ax=ax,shrink=.75,label=r'Field / $\varepsilon$')
ax=axes[0,1];g,_=inherited_second_leg(source,domain,v)
direct=inherited_reference(source,np.full_like(v,domain.u0),v)[0]
ax.plot(v-domain.v0,g[:,0],color='#245B9C',label='From signed horizon trace')
ax.plot(v-domain.v0,direct[:,0],'--',color='#BD7B1A',lw=1.2,label='Direct parent restriction')
ax.set(title='B  The second datum\nHas the inherited phase',xlabel=r'$v-v_H$',ylabel=r'First reduced component $\Psi_1$');ax.grid(alpha=.2);ax.legend(fontsize=8)
ax=axes[1,0];n=np.array([x['cells_per_leg'] for x in r['resolutions']])
ax.loglog(n,[x['relative_RMS_field_error'] for x in r['resolutions']],'-o',color='#245B9C',label='Field RMS')
ax.loglog(n,[x['horizon_normal_derivative_relative_RMS'] for x in r['resolutions']],'-s',color='#BD7B1A',label='Horizon normal derivative')
ax.loglog(n,[x['noether_flux']['relative_residual'] for x in r['resolutions']],'-^',color='#67846C',label='Noether flux balance')
ax.set(title='C  Three fixed grids\nResolve the same state',xlabel='Cells per null leg',ylabel='Relative error');ax.grid(which='both',alpha=.2);ax.legend(fontsize=8)
ax=axes[1,1];indices=np.arange(len(u));etadiag=eta[indices,indices]
proper=np.arcsinh(np.tan(etadiag))/source.parent.H
computed=field[indices,indices]
expected=np.einsum('ik,k->i',exact[indices,indices],np.asarray(axis[:2]))/(a[indices,indices]*np.sin(domain.chi_receiver))
ax.plot(proper-proper[0],computed,color='#245B9C',label='Continued field')
ax.plot(proper-proper[0],expected,'--',color='#BD7B1A',lw=1.1,label='Global parent reference')
ax.set(title='D  The receiving observer\nReads the inherited oscillation',xlabel='Proper duration since crossing',ylabel=r'Field / $\varepsilon$');ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('One parent state → two inherited null traces → a continuous receiving field\nFixed-background massive scalar; signed phase, physical flux and composition retained',fontsize=12)
fig.savefig(ROOT/'paper/figures/nsc-atlas-characteristic-handoff.png',dpi=170,metadata={'Software':'NSC inherited characteristic handoff'})
plt.close(fig)
print('Wrote inherited characteristic-handoff figure.')
