"""Original atlas plate: two sourced rooms, their restoring mode and clocks."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np

from recursive_horizons.nsc_atlas_occupied_wall import OccupiedWall

ROOT = Path(__file__).resolve().parents[1]


def main():
    w = OccupiedWall()
    states = [w.static_state(q) for q in (999,1000,1001)]
    signals = [w.central_parent_signal(s) for s in states]
    ink, teal, gold, red = "#172e43", "#247e8d", "#bb8132", "#b44f4c"
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,
                         "text.color":ink,"axes.labelcolor":ink,"axes.edgecolor":"#afb7bf",
                         "xtick.color":ink,"ytick.color":ink})
    fig = plt.figure(figsize=(10.6,7.8), facecolor="#faf9f5")
    grid = fig.add_gridspec(2,2, left=.06,right=.97,bottom=.17,top=.86,
                           width_ratios=[1.22,1], hspace=.5,wspace=.28)
    left = fig.add_subplot(grid[:,0])
    top, bottom = fig.add_subplot(grid[0,1]),fig.add_subplot(grid[1,1])
    fig.text(.06,.95,"FROM A SPHERE SPECTRUM TO TWO SOURCED ROOMS",weight="bold",fontsize=15)
    fig.text(.06,.914,"One material law, different occupied states, one expanding parent",fontsize=11)
    left.set_aspect("equal")
    for x,sig,q in [(-6,signals[0],999),(6,signals[2],1001)]:
        rb=sig["initial_join_radius"]
        left.add_patch(Circle((x,0),rb,facecolor="#e6f0ef",edgecolor=teal,lw=1.8,ls="--"))
        left.add_patch(Circle((x,0),1,facecolor="#eedcbf",edgecolor=gold,lw=2))
        left.scatter([x],[0],color=ink,s=18,zorder=6)
        left.annotate("",xy=(x,8),xytext=(x,0),arrowprops={"arrowstyle":"->","lw":1.8,"color":gold})
        left.scatter([x],[8],color=red,s=28,zorder=6)
        left.text(x,9.1,f'q = {q}',ha="center",weight="bold")
        left.text(x,7.0,f'1+z = {sig["one_plus_redshift"]:.6f}',ha="center",fontsize=9.5,
                  bbox={"facecolor":"#faf9f5","edgecolor":"none","pad":1.5},zorder=7)
        left.text(x,-5.7,"Occupied wall\nand regular cap",ha="center",fontsize=10)
    left.text(0,12,"Shared dust geometry",ha="center",color=teal,weight="bold",fontsize=12)
    left.text(0,-9,"Dashed spheres: compensated joins\nSolid spheres: material walls",ha="center",fontsize=10)
    left.set(xlim=(-11,11),ylim=(-11,14))
    left.axis("off")
    R=np.linspace(.8,1.2,170)
    top.plot(R, [1000*(w.mass_jet(r,1000)["mass"]-.18) for r in R],color=teal,lw=2.2)
    top.scatter([1],[0],color=gold,s=32,zorder=5)
    top.set(xlabel=r"Radius $R/R_{\rm ref}$",ylabel=r"$10^3(M-M_{\rm eq})/R_{\rm ref}$",
            title="A restoring radial mode")
    top.text(.46,.8,r"$\omega R_{\rm ref}=1.562561$",transform=top.transAxes,fontsize=10)
    q=[999,1000,1001]
    shift=[10000*(s["one_plus_redshift"]-signals[1]["one_plus_redshift"]) for s in signals]
    bottom.plot(q,shift,"o-",color=teal,lw=1.8,ms=6)
    bottom.axhline(0,color=ink,ls=":",lw=.8)
    bottom.set(xlabel="Filled positive sphere level q",ylabel=r"$10^4\Delta(1+z)$",
               title="A clock response with fixed couplings",xticks=q)
    for ax in (top,bottom):
        ax.spines[["top","right"]].set_visible(False)
        ax.grid(alpha=.13)
    fig.text(.06,.055,
             "Attachment diagram: parent comoving positions; interiors schematic. Linear radial stability only.\n"
             "Sphere spectrum and Fermi-wall inputs are attributed; the shared-parent signals are calculated here.",
             fontsize=9.5,color="#526673")
    out=ROOT/"docs/figures/nsc-atlas-occupied-wall.png"
    fig.savefig(out,dpi=190,facecolor=fig.get_facecolor())
    plt.close(fig)
    print(out)


if __name__ == "__main__":
    main()
