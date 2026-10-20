"""The attributed flux law, its shared evolving geometry and curvature readout."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np

from recursive_horizons.nsc_atlas_flux_transition import reference_geometry, evolve_reference


def main():
    _,outer,inner,parent=reference_geometry()
    solution,_=evolve_reference()
    t=np.linspace(0,.25,129)
    y=solution.sol(t)
    fig,axes=plt.subplots(1,3,figsize=(11.8,3.8),layout="constrained")
    ax=axes[0]
    for radius,color in [(1.2,"#e7edf5"),(.8,"#c0d9ef"),(.4,"#8ac7ba")]:
        ax.add_patch(Circle((0,0),radius,facecolor=color,edgecolor="#243e56",linewidth=1.4))
    ax.text(0,0,"Grandchild\n$\\Lambda L^2=0.125$",ha="center",va="center",fontsize=9)
    ax.text(0,.6,"Child: 0.5",ha="center",fontsize=9)
    ax.text(0,1.,"Parent: 1.125",ha="center",fontsize=9)
    ax.text(0,-1.43,"Same charge and tension at both walls",ha="center",fontsize=9)
    ax.set(xlim=(-1.5,1.5),ylim=(-1.55,1.5),aspect="equal",title="One action, inherited vacuum")
    ax.axis("off")
    ax=axes[1]
    ax.plot(t,y[0],label="Outer charged wall",color="#286396")
    ax.plot(t,y[4],label="Inner charged wall",color="#22836e")
    ax.plot(t,[parent.state(z)["radius"] for z in y[8]],label="Cosmological join",color="#5c667a",linestyle="--")
    ax.set(xlabel="Shared middle clock $T/L$",ylabel="Areal radius $R/L$",title="A common finite evolution")
    ax.legend(fontsize=8,loc="center",bbox_to_anchor=(.5,.35))
    ax.grid(alpha=.2)
    ax=axes[2]
    x=np.arange(3)
    rows=[reg.tidal(r) for reg,r in zip((outer.outside,outer.inside,inner.inside),(1.1,.6,.1))]
    ax.plot(x,[z["radial"] for z in rows],"o--",label="$A_r L^2$",color="#b97d32")
    ax.plot(x,[2*z["transverse"] for z in rows],"s--",label="$2A_t L^2$",color="#758a9a")
    ax.plot(x,[z["trace"] for z in rows],"o-",label="Sum = $\\Lambda L^2$",color="#753c76",linewidth=2)
    ax.set(xticks=x,xticklabels=["Parent","Child","Grandchild"],ylabel="Tidal acceleration / separation",title="Mass cancels from the sum")
    ax.legend(fontsize=8)
    ax.grid(alpha=.2)
    target=Path("docs/figures/nsc-atlas-flux-transition.png")
    fig.savefig(target,dpi=180)
    plt.close(fig)
    print(target)


if __name__=="__main__":
    main()
