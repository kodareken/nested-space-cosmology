"""Finite parent preparation and its jointly solved spatial geometry."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from recursive_horizons.nsc_atlas_parent_source import ParentScalarSource


def main():
    base=ParentScalarSource()
    bg,_=base.prepared.solve(tolerance=1e-9)
    Li,Lo=bg.p
    fig,ax=plt.subplots(1,2,figsize=(10.8,3.7),layout='constrained')
    for A,color in ((.01,'#245b8b'),(.02,'#288573')):
        model=ParentScalarSource(amplitude=A)
        sol,_=model.solve(bg,tolerance=1e-9)
        d=np.linspace(0,.5,501)
        f,_=model.packet(Lo-d,Lo)
        g,_=model.collar(Lo-d,Lo,sol.p[0])
        ax[0].plot(d,f+g,color=color,label=f'Supplied amplitude {A:g}')
        x=np.linspace(0,1,401)
        xi=np.r_[Li*x,Li+Lo*(1-x[::-1])]
        psi=np.r_[sol.sol(x)[0],sol.sol(x[::-1])[2]]
        ax[1].plot(xi,psi,color=color,label=f'Amplitude {A:g}')
    ax[0].axvspan(0,base.collar_width,color='#ba7840',alpha=.16,label='Solved wall-layer correction')
    ax[0].set(xlabel='Parent seed distance from wall / L_ref',ylabel='Added scalar field',title='Specified packet and required layer')
    ax[1].axhline(1,color='gray',linestyle=':',label='Regular reference')
    ax[1].axvline(Li,color='#ba7840',linestyle='--',label='Charged wall')
    ax[1].set(xlabel='Seed pole-to-pole distance / L_ref',ylabel='Conformal metric factor',title='Gravity responds to the finite source')
    for a in ax:
        a.legend(fontsize=8)
        a.grid(alpha=.2)
    out=Path('docs/figures/nsc-atlas-parent-source.png')
    fig.savefig(out,dpi=180)
    plt.close(fig)
    print(out)


if __name__=='__main__':main()
