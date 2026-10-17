"""Draw the actual two-level initial radius profile from its compact record."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from recursive_horizons.nsc_atlas_nested_slice import (
    bridge_chart, kottler_horizons, positive_sheet_velocity, signed_slope,
    taper_velocity, tension_wall, transition_velocity,
)

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args();started=time.process_time()
    payload=json.loads(args.input.read_text(encoding="utf-8"));record=payload["record"]
    name="src/recursive_horizons/nsc_atlas_nested_slice.py"
    if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=payload["source_sha256"][name]:
        raise ValueError("plot helper differs from the recorded method")
    if args.output.exists():
        raise FileExistsError("use a fresh figure path")
    inner=tension_wall(.8,0.,1.,.01,1.,2.)
    outer=tension_wall(1.5,.01,1.,.4,0.,2.)
    host=kottler_horizons(.01,1.);end=kottler_horizons(.4,0.)
    nodes,weights=np.polynomial.legendre.leggauss(8)
    pieces=[];distance=0.;samples=0

    def curve(a,b,radius,density,label,color):
        nonlocal distance,samples
        q=np.linspace(a,b,161)
        grid=(q[:-1,None]+q[1:,None])/2+(q[1:,None]-q[:-1,None])/2*nodes
        dl=(q[1:]-q[:-1])/2*np.sum(weights*density(grid),axis=1)
        samples+=grid.size+q.size
        l=distance+np.r_[0.,np.cumsum(dl)];R=radius(q)
        if np.any(np.diff(l)<=0) or not np.all(np.isfinite(R)):
            raise ValueError("invalid display profile")
        pieces.append((l,R,label,color));distance=float(l[-1])

    ident=lambda r:r
    curve(0.,.8,ident,lambda r:np.ones_like(r),"cap","#2563a6")
    def taper_density(r,lo,hi,v,M,H):
        value,_=taper_velocity(r,lo,hi,v)
        return 1/signed_slope(r,M,H,value,-1)
    curve(.8,.1,ident,lambda r:taper_density(r,.1,.8,inner["V_plus"],.01,1.),"host","#bb7929")
    curve(-np.sqrt(.1-host["rb"]),np.sqrt(.5-host["rb"]),
          lambda y:host["rb"]+y*y,
          lambda y:bridge_chart(y,host["rb"],1.,.01)["dl_dy"],"host","#bb7929")
    def right_density(r):
        value,_=transition_velocity(r,.5,.9,.01,1.)
        return 1/signed_slope(r,.01,1.,value,1)
    curve(.5,.9,ident,right_density,"host","#bb7929")
    curve(.9,1.5,ident,lambda r:np.ones_like(r),"host","#bb7929")
    curve(1.5,1.,ident,lambda r:taper_density(r,1.,1.5,outer["V_plus"],.4,0.),"exterior","#53606c")
    curve(-np.sqrt(1.-end["rb"]),np.sqrt(4.-end["rb"]),
          lambda y:end["rb"]+y*y,
          lambda y:bridge_chart(y,end["rb"],0.,.4)["dl_dy"],"exterior","#53606c")
    points={row["name"]:row for row in record["milestones"]}
    if abs(distance-points["retained_cut"]["l_over_ell"])>1e-6:
        raise ValueError("display integration differs from the compact calculation")
    if samples>100_000:
        raise RuntimeError("display sample limit exceeded")
    fig,(ax,zoom)=plt.subplots(1,2,figsize=(11.6,4.4),gridspec_kw={"width_ratios":[2.2,1]})
    seen=set()
    for l,R,label,color in pieces:
        ax.plot(l,R,color=color,lw=2.2,label=label if label not in seen else None)
        zoom.plot(l,R,color=color,lw=2.2);seen.add(label)
    for name,label in (("inner_wall","inner wall"),("outer_wall","outer wall")):
        row=points[name];x=row["l_over_ell"];y=row["R_over_ell"]
        ax.plot(x,y,"o",color="#b63536",ms=5)
        ax.annotate(label,(x,y),xytext=(7,9),textcoords="offset points",fontsize=9)
    for name,label in (("host_minimum","host throat"),("outer_minimum","outer throat")):
        row=points[name];ax.annotate(label,(row["l_over_ell"],row["R_over_ell"]),
            xytext=(7,17),textcoords="offset points",fontsize=9,
            arrowprops={"arrowstyle":"-","color":"#555"})
    c=points["cosmological_horizon"]
    ax.plot(c["l_over_ell"],c["R_over_ell"],"d",color="#202020",ms=4,label="host cosmological horizon")
    h=points["host_minimum"]
    zoom.set_xlim(h["l_over_ell"]-.16,h["l_over_ell"]+.23)
    zoom.set_ylim(-.002,.18)
    zoom.set_title("Inner throat resolved",fontsize=11)
    zoom.axhline(h["R_over_ell"],ls=":",color="#777",lw=.8)
    for panel in (ax,zoom):
        panel.set_xlabel(r"proper radius $l/\ell$")
        panel.set_ylabel(r"areal radius $R/\ell$")
        panel.grid(alpha=.2)
    ax.set_ylim(0,4.15);ax.legend(frameon=False,fontsize=8,loc="upper left")
    ax.set_title("Cap, host and exterior on one initial slice",fontsize=11)
    fig.suptitle("Two successive sourced boundaries",fontsize=15)
    fig.text(.5,.015,"Radius profile, not a Euclidean embedding or a causal-passage diagram. H, tension and masses are supplied.",
             ha="center",fontsize=8,color="#444")
    fig.tight_layout(rect=(0,.05,1,.94))
    fig.savefig(args.output,dpi=170)
    plt.close(fig)
    elapsed=time.process_time()-started
    if elapsed>30:
        raise RuntimeError("display CPU limit exceeded")
    print(f"Wrote {args.output}; {samples} regenerated points; {elapsed:.3f} CPU seconds")


if __name__=="__main__":
    main()
