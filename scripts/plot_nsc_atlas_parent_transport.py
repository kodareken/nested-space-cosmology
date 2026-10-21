"""Causal parent collar and independently normalized radial child trajectories."""
from pathlib import Path
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from recursive_horizons.nsc_atlas_scalar_preparation import ScalarPreparation
from recursive_horizons.nsc_atlas_parent_transport import RadialCapTransport


def main():
    model = ScalarPreparation()
    solution, _ = model.solve(tolerance=1e-9)
    t = RadialCapTransport.from_preparation(model,solution)
    fig, ax = plt.subplots(1,2,figsize=(10.8,3.7),layout="constrained")
    depth = np.linspace(0,t.parent_length,301)
    aa = solution.sol(1-depth/t.parent_length)[4]
    collar = t.null_collar()["depth"]
    ax[0].plot(depth,aa,color="#245b8b")
    ax[0].axvspan(0,collar,color="#288573",alpha=.16,label="Can signal center from eta = 0")
    ax[0].axvline(.3,color="#ba7840",linestyle="--",label="Specified source")
    ax[0].set(xlabel="Parent distance from wall / L_ref",ylabel="a / L_ref",
              title="Only a finite parent region is accessible",xlim=(0,t.parent_length))
    ax[0].legend(fontsize=8,loc="lower center")
    for gamma,color in ((1.1,"#ba7840"),(2.,"#245b8b"),(5.,"#753c76")):
        row = t.massive(.3,gamma)
        end = row.get("arrival_time",8.)
        tau = np.linspace(.08,end,301)
        z = 1/np.tanh(t.child_h*tau)
        c = 1/(t.child_h*row["radial_charge"])**2-1
        I = math.log(2/t.child_h)-np.log(z+np.sqrt(z*z+c))
        chi = row["log_label"]-I
        ax[1].plot(tau,np.maximum(chi,0),color=color,label=f"Source gamma = {gamma:g}")
        if row["central_arrival"]:
            ax[1].scatter([end],[0],color=color,s=25,zorder=5)
    ax[1].set(xlabel="Child proper time / (L_ref/c)",ylabel="Comoving distance to center",
              title="Arrival depends on source energy",ylim=(0,4),xlim=(0,8))
    ax[1].legend(fontsize=8)
    for a in ax:
        a.grid(alpha=.2)
    out = Path("docs/figures/nsc-atlas-parent-transport.png")
    fig.savefig(out,dpi=180)
    plt.close(fig)
    print(out)


if __name__ == "__main__":
    main()
