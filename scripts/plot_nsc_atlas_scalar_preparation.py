"""The scalar transition layer, its shared geometry and Lorentzian wall."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from recursive_horizons.nsc_atlas_scalar_preparation import ScalarPreparation


def main():
    model=ScalarPreparation()
    sol,_=model.solve(tolerance=1e-8)
    Li,Lo=sol.p
    x=np.linspace(0,1,501)
    y=sol.sol(x)
    fig,ax=plt.subplots(1,3,figsize=(11.5,3.7),layout="constrained")
    xi=np.concatenate((Li*x,Li+Lo*(1-x[::-1])))
    aa=np.concatenate((y[0],y[4,::-1]))
    ax[0].plot(xi,aa,color="#245b8b")
    ax[0].axvline(Li,color="#ba7840",linestyle="--")
    ax[0].set(xlabel="Pole-to-pole coordinate",ylabel="Scale factor a",title="Gravity and layer solved together")
    ax[0].text(.06,.1,"Child pole",transform=ax[0].transAxes,fontsize=9)
    ax[0].text(.66,.1,"Parent pole",transform=ax[0].transAxes,fontsize=9)
    z=np.linspace(0,1,1601)
    v=sol.sol(z)
    ax[1].plot(Li*(z-1),model.scalar.minimum(model.child_flux)+v[2],color="#288573",label="Child side")
    ax[1].plot(Lo*(1-z),model.scalar.minimum(model.parent_flux)+v[6],color="#245b8b",label="Parent side")
    ax[1].axvline(0,color="#ba7840",linestyle="--",label="Charged membrane")
    ax[1].set(xlim=(-.3,.3),xlabel="Signed distance from membrane",ylabel="Scalar field",title="A continuous spatial transition")
    ax[1].legend(fontsize=8)
    eta=np.linspace(0,1,101)
    R=sol.y[0,-1]
    ax[2].plot(R*eta,R*np.cosh(eta),color="#753c76")
    ax[2].set(xlabel="Wall proper time",ylabel="Wall areal radius",title="The continued expanding boundary")
    for a in ax:a.grid(alpha=.2)
    target=Path("docs/figures/nsc-atlas-scalar-preparation.png")
    fig.savefig(target,dpi=180)
    plt.close(fig)
    print(target)


if __name__=="__main__":main()
