"""Midpoint completion of the retained source/reference interface, without ODEs.

The source's scalar is analytic. Its geometry and transverse jets retain the
same 129-knot Hermite approximation, including its recorded residuals. Gauge
quadrature sensitivity is not an error enclosure. Original rows are copied,
not replaced by newly evaluated approximations.
"""
from collections.abc import Mapping
import json
import math
import time

import numpy as np

from .nsc_atlas_characteristic_initialization import Budget
from .nsc_atlas_characteristic_trace import SourceSnapshot
from .nsc_atlas_reference_wall_gauge import ReferenceWallGauge

SCHEMA = "nsc-atlas-refined-cache-v1"
GAUGE_TOLERANCE = 5e-12
INVERSE_TOLERANCE = 3e-13
SCOPE = ("original513 rows retained exactly; missing512 midpoint source rows use the same129-knot "
         "Hermite source geometry/transverse jets and analytic compact scalar; midpoint reference "
         "gauge uses orders64/32; no new source integration and no error enclosure")


def _budget(budget):
    if (not isinstance(budget,Budget) or type(budget.maximum_evaluations) is not int
            or not 0<budget.maximum_evaluations<=50000000
            or type(budget.evaluations) is not int or budget.evaluations<0
            or isinstance(budget.cpu_seconds,bool) or not math.isfinite(budget.cpu_seconds)
            or not 0<budget.cpu_seconds<=300. or not math.isfinite(budget.started)
            or budget.started>time.process_time()):
        raise ValueError("one live shared Budget within300CPU/50M is required")
    budget.check()
    return budget


def _plain(value):
    if isinstance(value,np.ndarray):
        return value.tolist()
    if isinstance(value,np.generic):
        return value.item()
    if isinstance(value,Mapping):
        return {key:_plain(item) for key,item in value.items()}
    if isinstance(value,(tuple,list)):
        return [_plain(item) for item in value]
    return value


def _table(rows,size,label,budget):
    if not isinstance(rows,Mapping) or not rows or any(not isinstance(k,str) for k in rows):
        raise ValueError(label+" must be a nonempty scalar-row mapping")
    result={}
    for key,value in rows.items():
        raw=np.asarray(value)
        if (raw.shape!=(size,) or raw.dtype.kind not in "iuf"
                or (isinstance(value,(list,tuple)) and any(isinstance(v,(bool,np.bool_)) for v in value))):
            raise ValueError(label+"."+key+" must have finite real row shape")
        budget.check(size)
        result[key]=np.asarray(raw,dtype=float)
        if not np.all(np.isfinite(result[key])):
            raise ValueError(label+" contains nonfinite values")
    if "x" not in result or np.any(np.diff(result["x"])<=0):
        raise ValueError(label+" requires strictly increasing x")
    return result


def _exact_rows(actual,expected,label,budget):
    if set(actual)!=set(expected):
        raise ValueError(label+" field inventory differs")
    for key in expected:
        a,b=np.asarray(actual[key],dtype=float),np.asarray(expected[key],dtype=float)
        budget.check(b.size)
        if a.shape!=b.shape or not np.array_equal(a.view(np.uint64),b.view(np.uint64)):
            raise ValueError(label+" differs: "+key)


def _midpoints(x,budget):
    if (x.shape!=(513,) or np.any(np.diff(x)<=0)
            or not np.allclose(np.diff(x),np.diff(x)[0],rtol=1e-10,atol=1e-14)):
        raise ValueError("original513 coordinates must be a uniform increasing grid")
    mid=x[:-1]+(x[1:]-x[:-1])/2
    budget.check(512)
    if np.any(mid<=x[:-1]) or np.any(mid>=x[1:]):
        raise ValueError("midpoints must be strictly inside their old cells")
    combined=np.empty(1025)
    combined[::2],combined[1::2]=x,mid
    return mid,combined


def _interleave(old,new,budget):
    if set(old)!=set(new):
        raise ValueError("old/new row inventories differ")
    result={}
    for key in old:
        a,b=_plain(old[key]),_plain(new[key])
        if len(a)!=513 or len(b)!=512:
            raise ValueError("interleave requires513 old and512 new rows")
        value=[None]*1025
        value[::2],value[1::2]=a,b
        result[key]=value
        budget.check(1025)
    return result


def _gauge_checks(rows,budget):
    required={"x","W","W_x","C","r","s","s_x","x_s","Ubar_s","T","T_x","T_s",
              "inverse_iterations","inverse_residual","normalization_defect"}
    if not required<=set(rows):
        raise ValueError("gauge rows omit orientation or normalization data")
    n=len(rows["x"])
    budget.check(8*n)
    if (any(np.any(rows[k]<=0) for k in ("W","W_x","C","s_x","x_s","Ubar_s","T_x","T_s"))
            or np.any(rows["r"]<=0) or np.any(rows["r"]>=1)
            or np.any(np.diff(rows["s"])<=0) or np.any(np.diff(rows["T"])<=0)):
        raise ValueError("gauge leaves its positive future-oriented branch")
    iterations=rows["inverse_iterations"]
    if (np.any(iterations!=np.floor(iterations)) or np.any(iterations<0) or np.any(iterations>32)
            or np.max(abs(rows["inverse_residual"]))>INVERSE_TOLERANCE):
        raise ValueError("gauge exceeds original inverse controls")
    defects=dict(normalization=rows["C"]*rows["W_x"]/rows["s_x"]**2-1,
        inverse_clock=rows["x_s"]*rows["s_x"]-1,
        outgoing_clock=rows["Ubar_s"]-rows["W_x"]/rows["s_x"],
        central_clock=rows["T_x"]-rows["T_s"]*rows["s_x"])
    maxima={key:float(np.max(abs(value))) for key,value in defects.items()}
    maxima["reported_normalization"]=float(np.max(abs(rows["normalization_defect"])))
    maxima["normalization_report_difference"]=float(np.max(abs(
        rows["normalization_defect"]-defects["normalization"])))
    if max(maxima.values())>GAUGE_TOLERANCE:
        raise ValueError("gauge normalization or clock identity exceeds tolerance")
    return {**maxima,"inverse_residual":float(np.max(abs(rows["inverse_residual"]))) }


def _compare_gauges(fine,coarse,budget):
    if set(fine)!=set(coarse):
        raise ValueError("gauge orders returned different field inventories")
    differences={}
    for key in fine:
        if key=="inverse_iterations":
            continue
        budget.check(len(fine[key]))
        differences[key]=float(np.max(abs(fine[key]-coarse[key])))
    if max(differences.values())>GAUGE_TOLERANCE:
        raise ValueError("reference32/64 sensitivity exceeds5e-12")
    return differences


def make_refined_cache(base,budget):
    """Return JSON-compatible1025-row cache; all work uses the caller budget.

    No initial budget count or elapsed time is stored, so rebuilding the same
    cache on a previously used budget yields the same deterministic record.
    The two SourceSnapshots charge2314 points including old-row verification
    queries; validation/comparison/copy arithmetic is charged separately.
    """
    budget=_budget(budget)
    start=budget.evaluations
    original_gauge=base["gauge"]["rows"]
    oldg=_table(original_gauge,513,"old gauge",budget)
    x=oldg["x"]
    mid,combined=_midpoints(x,budget)
    if x[0]!=base["geometry"]["x_corner"] or x[-1]!=base["geometry"]["x_last"]:
        raise ValueError("cache does not cover the original source interval")
    original_checks=_gauge_checks(oldg,budget)
    source_rows,source_summaries,source_checks={},{},{}
    for label in ("zero","finite"):
        data=base["initial_data"][label]
        old=_table(data["snapshot_rows"],513,label+" old source",budget)
        exact=_table(data["exact_rows"],257,label+" retained exact source",budget)
        _exact_rows({"x":old["x"]},{"x":x},"source/gauge keys",budget)
        _exact_rows({"x":exact["x"]},{"x":x[::2]},"exact source keys",budget)
        if (label=="zero" and data["alpha"]!=0.) or (label=="finite" and not data["alpha"]>0.):
            raise ValueError("preserve the zero/finite source amplitude identities")
        snapshot=SourceSnapshot(base["geometry"],alpha=data["alpha"],x_nodes=x[::4],
            source_rows={key:value[::2] for key,value in exact.items()},budget=budget)
        reproduced=snapshot.sample(x)
        _exact_rows(reproduced,old,label+" old source reproduction",budget)
        new=_table(snapshot.sample(mid),512,label+" midpoint source",budget)
        _exact_rows({"x":new["x"]},{"x":mid},"new source keys",budget)
        source_rows[label]=_interleave(data["snapshot_rows"],new,budget)
        source_summaries[label]=snapshot.summary()
        source_checks[label]=dict(old_rows_exact=True,maximum_midpoint_source_constraint=float(
            np.max(abs(new["source_constraint"]))),
            maximum_midpoint_curvature_defect=float(np.max(abs(new["source_h_xx_defect"]))),
            scope="actual Hermite residuals retained; source integration error is not enclosed")
    gauges,summaries,checks={},{},{}
    for order in (64,32):
        reference=ReferenceWallGauge(base["geometry"],order=order,budget=budget)
        rows=_table(reference.sample(mid),512,"midpoint gauge"+str(order),budget)
        _exact_rows({"x":rows["x"]},{"x":mid},"new gauge keys",budget)
        gauges[order],summaries[str(order)]=rows,reference.summary()
        checks[str(order)]=_gauge_checks(rows,budget)
        if (abs(summaries[str(order)]["endpoint_integral_defect"])>INVERSE_TOLERANCE
                or summaries[str(order)]["maximum_inverse_residual"]>INVERSE_TOLERANCE
                or summaries[str(order)]["maximum_inverse_iterations"]>32):
            raise ValueError("reference summary exceeds original inverse controls")
    differences=_compare_gauges(gauges[64],gauges[32],budget)
    gauge_rows=_interleave(original_gauge,gauges[64],budget)
    # Check orientation through interleaved old/new neighbors as well.
    joint_checks=_gauge_checks({key:np.asarray(value) for key,value in gauge_rows.items()},budget)
    source_points=sum(s["evaluations"] for s in source_summaries.values())
    if source_points!=2314:
        raise ValueError("unexpected source snapshot construction/query accounting")
    gauge_points=sum(s["evaluations"] for s in summaries.values())
    budget.check()
    total=budget.evaluations-start
    return _plain(dict(schema=SCHEMA,x_nodes=combined,source_rows=source_rows,gauge_rows=gauge_rows,
        diagnostics=dict(scope=SCOPE,old_rows_preserved_exactly=True,source=source_checks,
            original_gauge=original_checks,midpoint_gauge=checks,combined_gauge=joint_checks,
            gauge64_minus32_maximum_absolute=differences,gauge_comparison_tolerance=GAUGE_TOLERANCE,
            inverse_tolerance=INVERSE_TOLERANCE,comparison_is_rigorous_enclosure=False),
        execution=dict(evaluations=total,source_snapshot_points=source_points,
            gauge_evaluations=gauge_points,adapter_evaluations=total-source_points-gauge_points,
            source_snapshots=source_summaries,gauge_orders=summaries,
            source_integrations=0,regional_marches=0,time_evolutions=0)))


def validate_refined_cache(cache,base,budget):
    """Rebuild the same bounded interface and compare every deterministic field."""
    budget=_budget(budget)
    start=budget.evaluations
    try:
        supplied=json.dumps(cache,sort_keys=True,separators=(",",":"),allow_nan=False)
    except (ValueError,TypeError,OverflowError) as error:
        raise ValueError("refined cache must be finite JSON-compatible data") from error
    reconstructed=make_refined_cache(base,budget)
    if supplied!=json.dumps(reconstructed,sort_keys=True,separators=(",",":"),allow_nan=False):
        raise ValueError("refined cache differs from deterministic reconstruction")
    budget.check()
    return dict(validation_evaluations=budget.evaluations-start,source_integrations=0,
                regional_marches=0,time_evolutions=0)
