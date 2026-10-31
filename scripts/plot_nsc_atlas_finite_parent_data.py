"""Regenerate the finite nonspherical geometry and full-constraint comparison."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_finite_parent_data import solve_initial_data,source_profile

ROOT=Path(__file__).resolve().parents[1]


def main():
    record=json.loads((ROOT/'results/summaries/nsc-atlas-finite-parent-data-v1.json').read_text(encoding='utf-8'))
    with threadpool_limits(limits=4):
        data=solve_initial_data(16)
        theta=np.linspace(.001,np.pi-.001,121);phi=np.linspace(0,2*np.pi,241)
        T,P=np.meshgrid(theta,phi,indexing='ij')
        points=np.c_[np.cos(T).ravel(),(np.sin(T)*np.cos(P)).ravel(),np.zeros(T.size),(np.sin(T)*np.sin(P)).ravel()]
        psi=data.values(points).reshape(T.shape);f,D=source_profile(points,data.parent)
        rho=(data.parent.H**2*psi.ravel()**-4*D+data.parent.mass_squared*f*f)/(2*data.parent.vacuum)
    fig=plt.figure(figsize=(12,8),layout='constrained')
    a=fig.add_subplot(2,2,1,projection='mollweide')
    im=a.pcolormesh(P-np.pi,np.pi/2-T,psi,shading='auto',cmap='viridis')
    a.grid(alpha=.3);a.set_xticklabels([]);a.set_title('A  Positive conformal geometry: section X3=0',fontsize=11)
    fig.colorbar(im,ax=a,shrink=.7,label=r'$\Psi$')
    b=fig.add_subplot(2,2,2,projection='mollweide')
    im=b.pcolormesh(P-np.pi,np.pi/2-T,rho.reshape(T.shape),shading='auto',cmap='magma')
    b.grid(alpha=.3);b.set_xticklabels([]);b.set_title('B  The same finite scalar source: section X3=0',fontsize=11)
    fig.colorbar(im,ax=b,shrink=.7,label=r'$\rho_{exc}/V_*$')
    c=fig.add_subplot(2,2,3)
    c.semilogy([r['degree'] for r in record['levels']],[r['full_conformal_PDE_RMS_residual'] for r in record['levels']],'o-',color='#225b91')
    c.set(xlabel='Maximum harmonic degree',ylabel='Full conformal-constraint RMS residual',title='C  Resolve the full equation, beyond its retained projection')
    c.grid(alpha=.2)
    d=fig.add_subplot(2,2,4);steps=record['levels'][-1]['continuation']
    d.plot([0]+[r['epsilon']**2 for r in steps],[1]+[r['proper_volume_over_reference_round_volume'] for r in steps],'o-',label='Proper-volume ratio')
    d.plot([0]+[r['epsilon']**2 for r in steps],[1]+[r['volume_equivalent_length_over_reference_length'] for r in steps],'s-',label='Volume-equivalent length ratio')
    d.set(xlabel=r'Prepared amplitude $\epsilon^2$',ylabel='Relative to the reference round vacuum',title='D  Different initial preparations; this is not time evolution')
    d.grid(alpha=.2);d.legend(fontsize=9)
    fig.suptitle('Finite-amplitude geometry from the same inherited field\nFull nonlinear initial constraint; future formation requires nonlinear evolution',fontsize=13)
    out=ROOT/'paper/figures/nsc-atlas-finite-parent-data.png'
    fig.savefig(out,dpi=180);plt.close(fig)
    print(out.relative_to(ROOT))


if __name__=='__main__':main()
