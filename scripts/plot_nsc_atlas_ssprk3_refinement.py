#!/usr/bin/env python3
"""Render retained SSPRK3 histories; no source, cache or trajectory replay."""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PRODUCER = "scripts/plot_nsc_atlas_ssprk3_refinement.py"
DEFAULT_OUTPUT = "docs/figures/nsc-atlas-ssprk3-refinement.png"
SCHEMA = "nsc-atlas-ssprk3-refinement-figure-v1"
INPUTS = (("nsc-atlas-ssprk3-evolution-v1","ssprk3.json"),
          ("nsc-atlas-ssprk3-refinement-v1","refinement.json"))


def identity(path):
    data=Path(path).read_bytes()
    return dict(sha256=hashlib.sha256(data).hexdigest(),bytes=len(data))


class Allowance:
    def __init__(self):
        self.started=time.process_time(); self.points=0

    def check(self,points=0):
        self.points+=int(points)
        if self.points>200000 or time.process_time()-self.started>30:
            raise RuntimeError("saved-array figure allowance exhausted")


def retained(directory,filename,allowance):
    folder=ROOT/"results/development"/directory
    if {p.name for p in folder.iterdir()}!={filename,"manifest.json"} or any(p.is_symlink() for p in folder.iterdir()):
        raise ValueError("figure input inventory mismatch")
    manifest=json.loads((folder/"manifest.json").read_bytes())
    item=identity(folder/filename)
    if manifest.get("schema")!=directory+"-manifest" or manifest.get("outputs")!={filename:item}:
        raise ValueError("figure input hash mismatch")
    bindings=manifest["source_bindings"]
    for relative,expected in bindings.items():
        path=(ROOT/relative).resolve()
        if not path.is_relative_to(ROOT) or identity(path)["sha256"]!=expected:
            raise ValueError("figure input source mismatch: "+relative)
    record=json.loads((folder/filename).read_bytes())
    if record.get("schema")!=directory:
        raise ValueError("figure input schema mismatch")
    allowance.check()
    return record,{p.relative_to(ROOT).as_posix():identity(p) for p in (folder/filename,folder/"manifest.json")},bindings


def paired_wall(record,steps,allowance):
    chosen=[]
    for amplitude in ("zero","finite"):
        cases=[c for c in record["cases"] if (c["amplitude"],c["cells"],c["steps"])==(amplitude,32,steps)]
        if len(cases)!=1:
            raise ValueError("figure requires one matching amplitude/time case")
        chosen.append(cases[0])
    rows=chosen[0]["snapshots"]
    x=np.array([s["x"] for s in rows])
    if len(x)!=129 or not np.array_equal(x,[s["x"] for s in chosen[1]["snapshots"]]) or np.any(np.diff(x)<=0):
        raise ValueError("paired figure coordinates differ")
    c=np.array([f["physical"]["wall"]["c"]-z["physical"]["wall"]["c"]
                for z,f in zip(rows,chosen[1]["snapshots"])])
    if not np.all(np.isfinite(c)):
        raise ValueError("figure field contains nonfinite values")
    allowance.check(5*len(x))
    return x,c,chosen[1]["alpha"]


def render(output):
    allowance=Allowance()
    output=Path(output).resolve(); manifest_path=output.with_suffix(".manifest.json")
    if output.exists() or manifest_path.exists():
        raise FileExistsError("figure output is creation-only")
    if output.suffix.lower()!=".png" or not output.is_relative_to(ROOT):
        raise ValueError("figure must be a PNG inside the repository")
    producer=identity(ROOT/PRODUCER)
    old,old_inputs,old_bindings=retained(*INPUTS[0],allowance)
    new,new_inputs,new_bindings=retained(*INPUTS[1],allowance)
    parent=new["parent_input"]
    old_path="results/development/"+INPUTS[0][0]+"/"
    if (parent["sha256"]!=old_inputs[old_path+"ssprk3.json"]["sha256"]
            or parent["manifest_sha256"]!=old_inputs[old_path+"manifest.json"]["sha256"]):
        raise ValueError("refinement has a different retained SSPRK3 parent")
    histories=[paired_wall(old,k,allowance) for k in (128,256)]+[paired_wall(new,512,allowance)]
    x=histories[0][0]
    if any(not np.array_equal(x,h[0]) or h[2]!=histories[0][2] for h in histories[1:]):
        raise ValueError("figure changes coordinates or amplitude")
    curves=[h[1] for h in histories]
    gaps=[curves[1]-curves[0],curves[2]-curves[1]]
    for actual,recorded in ((gaps[0],new["diagnostics"]["retained_128_256"]["paired"]["arrays"]["wall"]["c"]),
                            (gaps[1],new["diagnostics"]["paired_time"]["arrays"]["wall"]["c"])):
        if not np.array_equal(actual,recorded):
            raise ValueError("saved figure difference does not reproduce retained diagnostics")
    order=new["diagnostics"]["profile_case_order"]
    index=order.index(dict(amplitude="finite",steps=512))
    brackets=new["diagnostics"]["profile_intervals"][index]
    if not np.array_equal(brackets["x"],x):
        raise ValueError("profile brackets use different saved cuts")
    crossings=np.asarray(brackets["strictly_interior_knots_per_saved_interval"])>0
    peaks=[int(np.argmax(abs(g))) for g in gaps]
    low=max(0,min(peaks)-8); high=min(len(x)-1,max(peaks)+8)
    allowance.check(14*len(x))
    colors=("#b36b24","#167d8d","#273d75")
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.titlesize":11,
        "axes.labelsize":10,"axes.spines.top":False,"axes.spines.right":False})
    fig=plt.figure(figsize=(11.6,7.7))
    grid=fig.add_gridspec(2,2,left=.085,right=.98,top=.85,bottom=.15,hspace=.48,wspace=.26)
    top=fig.add_subplot(grid[0,:]); lower=fig.add_subplot(grid[1,0]); zoom=fig.add_subplot(grid[1,1])
    fig.suptitle("Paired wall scalar: fixed-input SSPRK3 comparison",fontsize=16,fontweight="bold",y=.965,color="#17283a")
    fig.text(.5,.91,"32 cells  •  unchanged source, computational profile and spatial operators",ha="center",fontsize=11,color="#4b5968")
    for curve,steps,color,style in zip(curves,(128,256,512),colors,("--","-",":")):
        top.plot(x,1e5*curve,style,color=color,linewidth=1.8,label=str(steps)+" time steps")
    top.set(title="A  Finite minus zero at matching source labels",xlabel="Incoming source coordinate x",ylabel="Paired wall c (×0.00001)",xlim=x[[0,-1]])
    top.legend(frameon=False,fontsize=9,ncol=3,loc="upper left")
    labels=("256 − 128 steps","512 − 256 steps")
    for gap,color,label in zip(gaps,colors[1:],labels):
        lower.plot(x,1e6*gap,color=color,linewidth=1.8,label=label)
        zoom.plot(x,1e6*gap,color=color,linewidth=1.8,label=label)
    lower.set(title="B  Signed temporal changes of the paired field",xlabel="Incoming source coordinate x",ylabel="Paired c difference (×0.000001)",xlim=x[[0,-1]])
    lower.legend(frameon=False,fontsize=8.5,loc="best")
    # Mark only intervals known to contain a knot from retained endpoint labels.
    # These are not reconstructed internal stage paths or causal classifications.
    marked=np.where(crossings)[0]+1
    zoom.plot(x[marked],1e6*gaps[1][marked],"o",markersize=3.2,color=colors[2],
              label="512-step saved interval contains a profile knot")
    zoom.set(title="C  Region containing the largest temporal changes",xlabel="Incoming source coordinate x",ylabel="Paired c difference (×0.000001)",xlim=(x[low],x[high]))
    zoom.legend(frameon=False,fontsize=7.5,loc="best")
    for ax in (top,lower,zoom):
        ax.grid(alpha=.18,linewidth=.7); ax.tick_params(labelsize=9)
    fig.text(.085,.075,f"Finite amplitude α = {histories[0][2]:.8g}. Curves use 129 retained cuts; markers refer only to the fine finite history.",fontsize=9,color="#46525e")
    fig.text(.085,.04,"Temporal differences and knot locations are diagnostics, not a measured convergence order, causal attribution or error bound.",fontsize=9,color="#46525e")
    allowance.check(20*(3*len(x)+4*len(x)))
    stream=io.BytesIO(); fig.savefig(stream,dpi=180,facecolor="white",metadata={"Software":PRODUCER}); plt.close(fig)
    image=stream.getvalue()
    inputs={**old_inputs,**new_inputs}
    if identity(ROOT/PRODUCER)!=producer or any(identity(ROOT/p)!=v for p,v in inputs.items()):
        raise RuntimeError("figure producer/input changed during rendering")
    for relative,expected in {**old_bindings,**new_bindings}.items():
        if identity(ROOT/relative)["sha256"]!=expected:
            raise RuntimeError("figure bound source changed during rendering")
    allowance.check()
    metadata=dict(schema=SCHEMA,producer={PRODUCER:producer},inputs=inputs,
        outputs={output.relative_to(ROOT).as_posix():dict(sha256=hashlib.sha256(image).hexdigest(),bytes=len(image))},
        runtime=dict(numpy=np.__version__,matplotlib=matplotlib.__version__),
        execution=dict(cpu_seconds=time.process_time()-allowance.started,arithmetic_points=allowance.points,
            cpu_limit_seconds=30,point_limit=200000,renders=1,source_integrations=0,cache_reconstructions=0,
            time_evolutions=0,recursive_scientific_authentications=0,
            count_scope="array operands and plot vertices; excludes file bytes and rendering pixels"),
        comparison=dict(alpha=histories[0][2],steps=[128,256,512],cells=32,
            peak_paired_c=[float(max(abs(c))) for c in curves],max_temporal_change=[float(max(abs(g))) for g in gaps],
            maximum_indices=peaks,zoom_indices=[low,high],
            knot_markers="strictly interior profile knots between saved fine finite wall labels",
            scope="saved finite-window temporal diagnostic; no convergence-order or error-bound claim"))
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open("xb") as target: target.write(image)
    with manifest_path.open("xb") as target:
        target.write((json.dumps(metadata,indent=2,sort_keys=True,allow_nan=False)+"\n").encode())
    return metadata


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT/DEFAULT_OUTPUT)
    args=parser.parse_args()
    result=render(args.output)
    print(json.dumps({key:result[key] for key in ("outputs","execution","comparison")},sort_keys=True))


if __name__=="__main__": main()
