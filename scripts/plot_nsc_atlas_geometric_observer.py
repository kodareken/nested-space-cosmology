"""Pressure-work mass growth and the source-eliminated tidal relation."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from recursive_horizons.nsc_atlas_geometric_observer import (
    ScalarWallObserver,regular_wall_geometry,regular_bulk_Ruu)


def main():
    o=ScalarWallObserver()
    sol,_=o.prepared.solve(tolerance=1e-9)
    eta=np.linspace(0,1,101)
    rows=[o.reconstruct(regular_wall_geometry(o.prepared,sol,float(x))) for x in eta]
    R0=sol.y[0,-1]
    tau=R0*eta
    growth=np.array([r['vacuum_subtracted_mass'] for r in rows])-rows[0]['vacuum_subtracted_mass']
    p=sol.y[3,-1]**2/2-(o.prepared.child_vacuum+o.prepared.scalar.omega**2*sol.y[2,-1]**2/2)
    work=-o.prepared.scalar.kappa*(p+o.prepared.child_vacuum)*R0**3*(np.cosh(eta)**3-1)/6
    fig,ax=plt.subplots(1,2,figsize=(10.8,3.7),layout='constrained')
    ax[0].plot(tau,growth,color='#245b8b',label='Change in mass above vacuum')
    ax[0].plot(tau,work,color='#288573',linestyle='--',label='Pressure-work integral')
    ax[0].axhline(0,color='#ba7840',label='Integrated scalar influx: zero')
    ax[0].set(xlabel='Wall proper time / (L_ref/c)',ylabel='Geometric mass change / L_ref',title='Mass growth is not a flux measurement')
    P=np.linspace(.13,.29,151)
    C=o.tidal_quadratic_coefficient
    curve=lambda p:-o.prepared.scalar.kappa*o.prepared.child_vacuum-C*(p-o.ground_pressure)**2
    ax[1].plot(P,curve(P),color='#245b8b',label='Same-action tidal relation')
    ax[1].scatter([.2,.22],curve(np.array([.2,.22])),color='#288573',label='Synthetic calibration samples',zorder=4)
    ax[1].scatter([.24],[curve(.24)],color='#ba7840',marker='^',label='Withheld synthetic sample',zorder=4)
    r=o.reconstruct(regular_wall_geometry(o.prepared,sol,0.),measured_Ruu=regular_bulk_Ruu(sol))
    ax[1].scatter([r['pressure']],[r['measured_tidal_corrected']],color='#753c76',marker='*',s=90,label='Solved coupled reference',zorder=5)
    ax[1].set(xlabel='Pressure jump P',ylabel='Velocity-corrected tidal curvature Y',title='A test independent of source history')
    ax[0].legend(fontsize=9.5,loc='upper left')
    ax[1].legend(fontsize=9.5,loc='lower right')
    for a in ax:
        a.grid(alpha=.2)
    out=Path('docs/figures/nsc-atlas-geometric-observer.png')
    fig.savefig(out,dpi=180)
    plt.close(fig)
    print(out)


if __name__=='__main__':main()
