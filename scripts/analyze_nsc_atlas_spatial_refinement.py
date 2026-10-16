"""Saved-state spatial comparisons for the fixed-input 32/64-cell experiment.

No source, wall, metric, profile-fit or evolution call occurs here. Independent
pure constraints retain source interpolation curvature and diagnostic sampling
sensitivity. Two spatial resolutions do not establish an order or an enclosure.
"""
from __future__ import annotations

import copy
import math
import numpy as np

import record_nsc_atlas_ssprk3_refinement as inherited
from recursive_horizons.nsc_atlas_vacuum_departure import vacuum_slice

FIELDS = ("r", "L", "c")
CHANNELS = ("E_x", "inherited", "residual", "transport", "E_Ubar", "mass")
plain = inherited.plain


def array(value, shape, name, budget):
    raw = np.asarray(value)
    if raw.shape != shape or raw.dtype.kind not in "iuf" or not np.all(np.isfinite(raw)):
        raise ValueError(name + " must have the declared finite real shape")
    budget.check(raw.size)
    return np.asarray(raw, dtype=float)


def same(a, b, name):
    if not np.array_equal(a, b):
        raise ValueError(name + " must match exactly")


def restrict_case(case, budget):
    """Restrict exterior nodes only; do not reinterpret the native cap state."""
    n = case["cells"]
    if type(n) is not int or n < 8 or n > 128 or n % 2:
        raise ValueError("restriction needs an even bounded exterior cell count")
    rows = []
    for sample in case["snapshots"]:
        state, physical = dict(sample["state"]), dict(sample["physical"])
        for owner, key in ((state, "departure"), (physical, "increment")):
            owner[key] = {k: array(owner[key][k], (n-1,), key+"."+k, budget)[1::2].copy()
                          for k in FIELDS}
        rows.append(dict(sample, state=state, physical=physical))
    return dict(case, cells=n//2, snapshots=rows)


def restrict_table(table, budget):
    """Common eta nodes, unchanged wall, SUM adjacent total physical cap integrals."""
    result = dict(exterior={}, wall=table["wall"], cap={})
    for key, values in table["exterior"].items():
        values = np.asarray(values)
        if values.ndim != 2 or values.shape[1] % 2 != 1:
            raise ValueError("exterior restriction requires an even number of cells")
        budget.check(values.size)
        result["exterior"][key] = values[:, ::2]
    for key, values in table["cap"].items():
        values = np.asarray(values)
        if values.ndim != 2 or values.shape[1] % 2:
            raise ValueError("cap restriction requires pairs of total integrals")
        budget.check(values.size)
        result["cap"][key] = values[:, ::2]+values[:, 1::2]
    return result


def source_constraint(source, alpha, budget):
    n = len(source["r"])
    v = {k: array(source[k], (n,), "source."+k, budget)
         for k in ("r", "r_x", "r_xx", "log_C_x", "c", "c_x")}
    if not math.isfinite(alpha) or alpha < 0 or np.any(v["r"] <= 0):
        raise ValueError("invalid source radius or alpha")
    u = v["c"]/v["r"]
    ux = (v["c_x"]-u*v["r_x"])/v["r"]
    budget.check(5*n)
    return v["r_xx"]-v["log_C_x"]*v["r_x"]+alpha*v["r"]*ux**2/(2*(1+alpha*u*u/6))


def field_summaries(fields, masks, x, eta, U, budget):
    result = {}
    for key, value in fields.items():
        result[key] = {}
        for name, mask in masks.items():
            budget.check(value.size)
            ix = np.flatnonzero(mask)
            if not len(ix):
                result[key][name] = dict(count=0, maximum_absolute=None, rms=None, location=None)
                continue
            selected = value.ravel()[ix]
            k = int(ix[np.argmax(abs(selected))]); i,j = np.unravel_index(k,value.shape)
            scale = float(np.max(abs(selected)))
            result[key][name] = dict(count=len(ix), minimum=float(np.min(selected)),
                maximum=float(np.max(selected)), maximum_absolute=scale,
                rms=scale*float(np.sqrt(np.mean((selected/scale)**2))) if scale else 0.,
                location=dict(interior_index=[int(i),int(j)],x=float(x[i]),eta=float(eta[j]),
                              Ubar=float(U[i,j]),signed=float(value[i,j])))
    return result


def constraint_grid(history, x, eta, W, radius, source, alpha, fronts, budget):
    """Subtract the ACTUAL inherited source constraint, never a zero ODE target."""
    nt, ne = len(x), len(eta)
    if not 4 <= nt <= 129 or not 5 <= ne <= 129:
        raise ValueError("constraint grid exceeds bounded saved coordinates")
    x=array(x,(nt,),"x",budget); eta=array(eta,(ne,),"eta",budget)
    W=array(W,(nt,),"W",budget); radius=array(radius,(nt,ne),"radius",budget)
    if (np.any(np.diff(x)<=0) or np.any(np.diff(eta)<=0) or eta[0]!=0 or eta[-1]!=1
            or np.any(np.diff(W)<=0) or W[0]<=0 or np.any(radius<=0)):
        raise ValueError("invalid ordered constraint coordinates or radius")
    for axis in (x,eta):
        if not np.allclose(np.diff(axis),np.diff(axis)[0],rtol=1e-10,atol=1e-14):
            raise ValueError("constraint coordinates must be uniform")
    shape=(nt-2,ne-2)
    Ex=array(history["E_x"],shape,"E_x",budget)
    EU=array(history["E_Ubar"],shape,"E_Ubar",budget)
    mass=array(history["mass"],shape,"mass",budget)
    Es=source_constraint(source,alpha,budget)
    saved=array(history["source_constraint"],(nt,),"history source E",budget)
    if not np.allclose(saved,Es,rtol=2e-13,atol=3e-14):
        raise ValueError("history source constraint disagrees with actual cached jets")
    sr=np.asarray(source["r"]); r=radius[1:-1,1:-1]
    weight=sr*Es; baseline=weight[1:-1,None]/r
    fields=dict(E_x=Ex,inherited=baseline,residual=Ex-baseline,
                transport=r*Ex-weight[1:-1,None],E_Ubar=EU,mass=mass)
    U=W[1:-1,None]*eta[None,1:-1]
    bottom=W[1:-1,None]*eta[None,2:]<=W[0]
    wall=W[1:-1,None]*eta[None,:-2]>W[0]
    near=np.zeros(nt-2,dtype=bool)
    fronts=array(fronts,(2,),"source fronts",budget)
    if not x[0]<=fronts[0]<fronts[1]<=x[-1]:
        raise ValueError("source fronts must lie in the declared interval")
    for front in fronts:
        near |= (x[:-2]<=front)&(front<=x[2:])
    near=np.broadcast_to(near[:,None],shape).copy()
    masks=dict(all_interior=np.ones(shape,dtype=bool),front_near=near,away_from_front=~near,
               bottom_origin=bottom,wall_origin=wall,separator=~(bottom|wall))
    budget.check(8*Ex.size)
    return dict(x=x[1:-1],eta=eta[1:-1],Ubar=U,W=W[1:-1],radius=r,
        source_weighted=weight[1:-1],source_constraint=Es,
        source_arithmetic_difference=saved-Es,arrays=fields,masks=masks,
        summaries=field_summaries(fields,masks,x[1:-1],eta[1:-1],U,budget),
        source_curvature="actual Hermite right trace at interior knots; no averaging correction",
        incoming_wall_baseline="unknown; nonzero E_Ubar is not by itself a junction failure")


def select_grid(grid, rows=slice(None), cols=slice(None)):
    """Select common physical coordinate labels, retaining each stencil's masks."""
    result=dict(grid)
    for key in ("x","W","source_weighted"):
        result[key]=np.asarray(grid[key])[rows]
    result["eta"]=np.asarray(grid["eta"])[cols]
    for key in ("Ubar","radius"):
        result[key]=np.asarray(grid[key])[rows,:][:,cols]
    for key in ("arrays","masks"):
        result[key]={k:np.asarray(v)[rows,:][:,cols] for k,v in grid[key].items()}
    return result


def compare_constraints(left, right, budget, *, same_spatial_stencil=False):
    """Signed left-minus-right; union front reach and conservative origin masks."""
    for key in ("x","eta","W","Ubar","source_weighted"):
        same(left[key],right[key],"common "+key)
    near=left["masks"]["front_near"]|right["masks"]["front_near"]
    bottom=left["masks"]["bottom_origin"]&right["masks"]["bottom_origin"]
    wall=left["masks"]["wall_origin"]&right["masks"]["wall_origin"]
    masks=dict(all_interior=np.ones_like(near),front_near=near,away_from_front=~near,
               bottom_origin=bottom,wall_origin=wall,separator=~(bottom|wall))
    fields={k:np.asarray(left["arrays"][k])-np.asarray(right["arrays"][k]) for k in CHANNELS}
    if same_spatial_stencil and np.any(fields["E_Ubar"]!=0):
        raise ValueError("time coarsening changed an identical pure-U stencil")
    budget.check((len(CHANNELS)+4)*near.size)
    return dict(arrays=fields,masks=masks,x=left["x"],eta=left["eta"],Ubar=left["Ubar"],
        summaries=field_summaries(fields,masks,left["x"],left["eta"],left["Ubar"],budget),
        same_spatial_stencil=same_spatial_stencil,
        interpretation="diagnostic differentiation or actual evolution difference as named; neither is an error enclosure")


def backward_jets(values, *, points, stride, budget):
    """Generic saved-node backward stencils; 3 means a 3/4-node first/second pair."""
    z=np.asarray(values); nt,ne=z.shape if z.ndim==2 else (0,0)
    if (not 2<=nt<=129 or not 9<=ne<=129 or type(points)is not int or points not in (3,5)
            or type(stride)is not int or stride not in (1,2)):
        raise ValueError("invalid bounded backward stencil")
    z=array(z,(nt,ne),"departure",budget); h=stride/(ne-1)
    if points==3:
        first=np.array([3.,-4.,1.])/(2*h); second=np.array([2.,-5.,4.,-1.])/h**2
    else:
        first=np.array([25.,-48.,36.,-16.,3.])/(12*h)
        second=np.array([35.,-104.,114.,-56.,11.])/(12*h*h)
    i=ne-1-stride*np.arange(len(first)); j=ne-1-stride*np.arange(len(second))
    if min(i[-1],j[-1])<0:
        raise ValueError("backward stencil exits grid")
    budget.check(nt*(len(i)+len(j)))
    return dict(first=z[:,i]@first,second=z[:,j]@second,first_indices=i,second_indices=j,eta_spacing=h)


def wall_estimate(Z, wall, vacuum, W, alpha, *, points, stride, budget):
    nt=len(W); W=array(W,(nt,),"wall W",budget)
    if np.any(W<=0) or isinstance(alpha,bool) or not math.isfinite(alpha) or alpha<0:
        raise ValueError("invalid wall width/alpha")
    r=array(wall["r"],(nt,),"wall r",budget); c=array(wall["c"],(nt,),"wall c",budget)
    if np.any(r<=0): raise ValueError("positive wall radius required")
    jets={k:backward_jets(Z[k],points=points,stride=stride,budget=budget) for k in FIELDS}
    rU=np.asarray(vacuum["r_Ubar"])+jets["r"]["first"]/W
    rUU=np.asarray(vacuum["r_Ubar_Ubar"])+jets["r"]["second"]/W**2
    LU=np.asarray(vacuum["log_C_Ubar"])+jets["L"]["first"]/W; cU=jets["c"]["first"]/W
    u=c/r; uU=(cU-u*rU)/r
    terms=dict(curvature=rUU,lapse_product=-LU*rU,scalar=alpha*r*uU*uU/(2*(1+alpha*u*u/6)))
    E=sum(terms.values()); B=r*E
    if not all(np.all(np.isfinite(v)) for v in (*terms.values(),E,B)):
        raise ValueError("nonfinite wall constraint")
    budget.check(9*nt)
    arrays=dict(E_Ubar=E,weighted=B,**terms)
    return dict(points=points,stride=stride,eta_spacing=jets["r"]["eta_spacing"],
        first_indices=jets["r"]["first_indices"],second_indices=jets["r"]["second_indices"],
        jets=dict(r_Ubar=rU,r_Ubar_Ubar=rUU,log_C_Ubar=LU,c_Ubar=cU,u_Ubar=uU),
        arrays=arrays,statistics={k:inherited.statistics(v) for k,v in arrays.items()})


def wall_departures(case, source, vacuum, budget):
    n=case["cells"]; nt=len(case["snapshots"])
    wall={k:array([s["physical"]["wall"][k] for s in case["snapshots"]],(nt,),"wall."+k,budget) for k in FIELDS}
    x=np.asarray([s["x"] for s in case["snapshots"]])
    offsets=dict(r=np.asarray(source["r"])-x,L=np.asarray(source["log_C"]),c=np.asarray(source["c"]))
    Z={}
    for key in FIELDS:
        middle=array([s["state"]["departure"][key] for s in case["snapshots"]],(nt,n-1),"Z."+key,budget)
        bg=np.asarray(vacuum[key]) if key!="c" else np.zeros(nt)
        Z[key]=np.column_stack((np.zeros(nt),middle,wall[key]-bg-offsets[key]))
    budget.check(3*nt)
    return Z,wall


def wall_comparison(left,right,x,budget):
    fields={k:np.asarray(left["arrays"][k])-np.asarray(right["arrays"][k]) for k in left["arrays"]}
    budget.check(sum(v.size for v in fields.values()))
    return dict(arrays=fields,statistics={k:inherited.statistics(v) for k,v in fields.items()},
        locations={k:dict(index=int(np.argmax(abs(v))),x=float(x[np.argmax(abs(v))]),
                         signed=float(v[np.argmax(abs(v))])) for k,v in fields.items()})


def diagnostics(record,parent,base,profile,budget):
    """Two fine64 cases versus fixed32 cases; only saved arithmetic/reconstruction."""
    budget.check(); initial=budget.evaluations
    if len(record["cases"])!=2 or len(parent["cases"])!=2:
        raise ValueError("ordered zero/finite cases required")
    x=np.asarray(parent["cache"]["x_nodes"])[::8]
    if len(x)!=129: raise ValueError("original129 saved cuts required")
    source_rows=base["initial_data"]; gauge=inherited.subset(base["gauge"]["rows"],slice(None,None,4))
    same(gauge["x"],x,"base/refined saved cuts")
    for key in ("x","W","W_x","W_xx"):
        same(gauge[key],np.asarray(parent["cache"]["gauge_rows"][key])[::8],"base/refined gauge "+key)
    W=np.asarray(gauge["W"]); vacuum=vacuum_slice(W,x,budget=budget)
    histories=[]; transports=[]; tables=[]; old_tables=[]; caps=[]; controls=[]; wall_results=[]
    for ci,label in enumerate(("zero","finite")):
        case,old=record["cases"][ci],parent["cases"][ci]
        if ((case["amplitude"],case["cells"],case["steps"])!=(label,64,512)
                or (old["amplitude"],old["cells"],old["steps"])!=(label,32,512)
                or case["alpha"]!=old["alpha"]):
            raise ValueError("fixed64/32 spatial cases and amplitudes required")
        for c in (case,old): same([s["x"] for s in c["snapshots"]],x,"saved cuts")
        source=inherited.subset(source_rows[label]["snapshot_rows"],slice(None,None,4))
        cached=inherited.subset(parent["cache"]["source_rows"][label],slice(None,None,8))
        for k in ("x","r","r_x","r_xx","log_C","log_C_x","c","c_x"):
            same(source[k],cached[k],"old/refined source "+k)
        table=inherited.tables(case,source_rows[label]["snapshot_rows"],budget)
        older=inherited.tables(old,source_rows[label]["snapshot_rows"],budget)
        tables.append(table); old_tables.append(older)
        restricted=restrict_case(case,budget); restricted_table=restrict_table(table,budget)
        history={}; grids={}
        for name,proxy in (("native64",case),("restricted32",restricted)):
            history[name]={}
            for nt,stride in ((129,1),(65,2)):
                raw,dep=inherited._histories(base,proxy,budget,stride)
                history[name]["raw"+str(nt)]=inherited.compact_history(raw)
                history[name]["departure"+str(nt)]=inherited.compact_history(dep)
                radius=(table if name=="native64" else restricted_table)["exterior"]["r"][::stride]
                src={k:np.asarray(v)[::stride] for k,v in source.items()}
                grids[name+"_"+str(nt)]=constraint_grid(dep,x[::stride],np.linspace(0,1,proxy["cells"]+1),
                    W[::stride],radius,src,case["alpha"],base["geometry"]["support_x"],budget)
        for nt,stride in ((129,1),(65,2)):
            saved=parent["diagnostics"]["histories"][ci]["departure"+str(nt)]
            src={k:np.asarray(v)[::stride] for k,v in source.items()}
            grids["old32_"+str(nt)]=constraint_grid(saved,x[::stride],np.linspace(0,1,33),W[::stride],
                older["exterior"]["r"][::stride],src,case["alpha"],base["geometry"]["support_x"],budget)
        sampling={name:compare_constraints(grids[name+"_65"],select_grid(grids[name+"_129"],rows=slice(1,None,2)),
                    budget,same_spatial_stencil=True) for name in ("native64","restricted32","old32")}
        spatial={}
        for nt in (129,65):
            native=select_grid(grids["native64_"+str(nt)],cols=slice(1,None,2))
            restricted_grid=grids["restricted32_"+str(nt)]; old_grid=grids["old32_"+str(nt)]
            spatial[str(nt)]=dict(native_common_minus_restricted=compare_constraints(native,restricted_grid,budget),
                restricted_minus_old=compare_constraints(restricted_grid,old_grid,budget),
                native_common_minus_old=compare_constraints(native,old_grid,budget))
        histories.append(history); transports.append(dict(amplitude=label,grids=grids,sampling=sampling,spatial=spatial))
        caps.append(inherited.cap_history(case,profile,budget))
        quiet=x<=base["geometry"]["support_x"][0]
        channels={g+"."+k:inherited.statistics(np.asarray(table[g][k])[quiet])
                  for g,k in (("exterior","c"),("wall","c"),("cap","Ip"),("cap","Iq"))}
        budget.check(sum(np.asarray(table[g][k])[quiet].size for g,k in
                        (("exterior","c"),("wall","c"),("cap","Ip"),("cap","Iq"))))
        controls.append(dict(quiet_scalar_channels=channels,
            minimum_upper_label_margin=min(caps[-1]["arrays"]["upper_label_margin"]),
            final_wall_kinematic_defect=case["final"]["wall_kinematic_defect"],
            final_canonical_wall_flux_defect=case["final"]["canonical_wall_flux_defect"],
            final_conformal_wall_flux_defect=case["final"]["conformal_wall_flux_defect"],
            constructed_identities_are_accuracy_tests=False))
        estimates={}
        for name,c in (("fine64",case),("old32",old)):
            Z,wall=wall_departures(c,source,vacuum,budget); estimates[name]={}
            for points in (3,5):
                for stride in ((1,2) if name=="fine64" else (1,)):
                    estimates[name][f"p{points}_s{stride}"]=wall_estimate(Z,wall,vacuum,W,case["alpha"],
                        points=points,stride=stride,budget=budget)
        differences={}
        for points in (3,5):
            f=estimates["fine64"]; o=estimates["old32"][f"p{points}_s1"]
            differences[str(points)]=dict(native_minus_old=wall_comparison(f[f"p{points}_s1"],o,x,budget),
                equal_width_minus_old=wall_comparison(f[f"p{points}_s2"],o,x,budget),
                native_minus_equal_width=wall_comparison(f[f"p{points}_s1"],f[f"p{points}_s2"],x,budget))
        wall_results.append(dict(amplitude=label,estimates=estimates,differences=differences))
    paired=inherited._difference(tables[1],tables[0]); oldpaired=inherited._difference(old_tables[1],old_tables[0])
    oracle=copy.deepcopy(record["cases"][0])
    for row in oracle["snapshots"]: row["state"]["cap"]=row["physical"]["cap"]
    output=dict(histories=histories,constraint_transport=transports,cap_histories=caps,controls=controls,
        zero_cap_pullback=inherited.cap_pullback(oracle,budget),
        finite_minus_zero64=inherited.comparison(tables[1],tables[0],x,budget),
        spatial_zero=inherited.comparison(restrict_table(tables[0],budget),old_tables[0],x,budget),
        spatial_finite=inherited.comparison(restrict_table(tables[1],budget),old_tables[1],x,budget),
        paired_spatial=inherited.comparison(restrict_table(paired,budget),oldpaired,x,budget),
        wall_curvature=dict(x=x,W=W,cases=wall_results,analytic_vacuum_points=len(x),
            birth_transport_evaluations=0,scope="one-sided saved curvature; same-width comparison, no prescribed wall acceleration"),
        comparison_coordinates="common x, exterior eta and cap normalized radial cells; physical radii may differ",
        cap_restriction="sum adjacent TOTAL physical dr integrals, not average or departure-only restriction",
        source_constraint_is_zero=False,spatial_order_established=False,uniform_error_enclosure=False)
    budget.check(); output["diagnostic_evaluations"]=budget.evaluations-initial
    return plain(output)
