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
from .nsc_atlas_validation import source_snapshot_roundoff_tolerances

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


def _replay_rows(actual,expected,label,*,geometry=None,alpha=0.,source_knots=None,budget=None,gauge=False):
    """Validation-only roundoff comparison; copied rows use _exact_rows instead.

    The source-frame allowance propagates null-label roundoff into scalar jets
    and cancellation residuals such as h_xx-Q; it is not a source error bound.
    Coordinates stay exact. Fresh adaptive stopping counts may differ within
    the original inverse bounds; copied historical counts stay exact. These comparisons do not add
    scientific queries; when replacing an existing comparison its charges are
    retained. The caller's CPU clock covers every comparison.
    """
    if set(actual)!=set(expected):
        raise ValueError(label+" field inventory differs")
    conditioned = {} if geometry is None else source_snapshot_roundoff_tolerances(
        expected,geometry,alpha,source_knots=source_knots)
    for key in expected:
        a,b=np.asarray(actual[key],dtype=float),np.asarray(expected[key],dtype=float)
        if budget is not None:
            budget.check(b.size)
        if a.shape!=b.shape or not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
            raise ValueError(label+" differs: "+key)
        if key=="inverse_iterations" and gauge:
            same=bool(np.all(a==np.floor(a)) and np.all(b==np.floor(b))
                and np.all((a>=0)&(a<=32)) and np.all((b>=0)&(b<=32))
                and all(type(x) is type(y) for x,y in zip(_plain(actual[key]),_plain(expected[key]))))
        elif key in {"x","inverse_iterations"}:
            same=np.array_equal(a.view(np.uint64),b.view(np.uint64))
            if key=="inverse_iterations":
                same=bool(same and all(type(x) is type(y) for x,y in zip(_plain(actual[key]),_plain(expected[key]))))
        else:
            atol=GAUGE_TOLERANCE if gauge else 2e-13+conditioned.get(key,0.)
            same=np.allclose(a,b,rtol=2e-12,atol=atol)
            if key in {"c","c_x","c_xx"} and key in conditioned:
                quiet=(b==0)&(conditioned[key]==0)
                same=bool(same and np.all(a[quiet]==0))
        if not same:
            if key=="inverse_iterations" and gauge:
                detail=" (integer stopping counts must retain types and stay within 0..32)"
            elif key in {"x","inverse_iterations"}:
                detail=" (exact coordinate/count comparison)"
            else:
                allowance=np.broadcast_to(atol+2e-12*abs(b),b.shape)
                excess=abs(a-b)-allowance
                worst=int(np.argmax(excess))
                detail=(f" (maximum difference={float(np.max(abs(a-b))):.17g}, "
                        f"maximum allowance={float(np.max(allowance)):.17g}, "
                        f"worst index={worst}, difference={float(abs(a-b).flat[worst]):.17g}, "
                        f"allowance={float(allowance.flat[worst]):.17g})")
            raise ValueError(label+" differs: "+key+detail)


def _same_replayed(actual,expected,absolute=None):
    """Exact structure/counters with bounded floating replay differences."""
    if isinstance(expected,dict):
        return isinstance(actual,dict) and set(actual)==set(expected) and all(
            _same_replayed(actual[k],v,absolute.get(k) if isinstance(absolute,dict) else None)
            for k,v in expected.items())
    if isinstance(expected,list):
        return isinstance(actual,list) and len(actual)==len(expected) and all(
            _same_replayed(a,b) for a,b in zip(actual,expected))
    if type(expected) is float:
        return type(actual) is float and math.isfinite(actual) and math.isfinite(expected) and math.isclose(
            actual,expected,rel_tol=2e-12,abs_tol=2e-13 if absolute is None else absolute)
    return type(actual) is type(expected) and actual==expected


def _source_summary_allowances(reconstructed,base):
    """Only maxima of conditioned source residuals inherit their row allowance."""
    result={}
    for label in ("zero","finite"):
        rows={key:np.asarray(value)[1::2] for key,value in reconstructed["source_rows"][label].items()}
        data=base["initial_data"][label]
        limits=source_snapshot_roundoff_tolerances(rows,base["geometry"],data["alpha"],
                                                  source_knots=_source_knots(data))
        result[label]={key:2e-13+float(np.max(limits.get(field,0.))) for key,field in (
            ("maximum_midpoint_source_constraint","source_constraint"),
            ("maximum_midpoint_curvature_defect","source_h_xx_defect"))}
    return {"diagnostics":{"source":result}}


def _source_knots(data):
    """Original129 knot values, including the transverse mixed-slope inputs."""
    exact=data.get("exact_rows",{})
    if not {"x","h","h_x"}<=set(exact):
        return None
    return {key:np.asarray(value)[::2] for key,value in exact.items()}


def _same_copied_rows(actual,original):
    """Copied JSON rows retain numeric types and signed zero as well as values."""
    if set(actual)!=set(original):
        return False
    copied={key:actual[key][::2] for key in original}
    return json.dumps(copied,sort_keys=True,separators=(",",":"),allow_nan=False)==json.dumps(
        _plain(original),sort_keys=True,separators=(",",":"),allow_nan=False)


_ADAPTIVE_GAUGE_FIELDS=frozenset({"evaluations","quadrature_nodes","maximum_inverse_iterations",
                                "maximum_inverse_residual"})
_STOPPING_DIAGNOSTIC_ALLOWANCES=dict(
    midpoint_gauge={str(order):{"inverse_residual":INVERSE_TOLERANCE} for order in (64,32)},
    combined_gauge={"inverse_residual":INVERSE_TOLERANCE},
    gauge64_minus32_maximum_absolute={"inverse_residual":2*INVERSE_TOLERANCE})


def _gauge_work(summary,order,midpoint_rows=None):
    """Check a historical or fresh inverse-work receipt on its own terms.

    Brent's stopping iteration can change under roundoff. Each of the 512
    midpoint inversions still converges within 32 iterations and may request
    at most one fresh integral per iteration plus a returned-root evaluation.
    Each integral costs exactly `order` nodes; construction adds one endpoint
    integral and three fixed checks. These are work ceilings, not error bounds.
    """
    for key in ("quadrature_order","evaluations","quadrature_nodes","query_points",
                "cache_hits","inverse_calls","maximum_inverse_iterations"):
        if type(summary.get(key)) is not int or summary[key]<0:
            raise ValueError("gauge work receipt requires exact nonnegative integer "+key)
    if (summary["quadrature_order"]!=order or summary["query_points"]!=512
            or summary["cache_hits"]!=0 or summary["inverse_calls"]!=512
            or not 1<=summary["maximum_inverse_iterations"]<=32):
        raise ValueError("gauge work receipt violates fixed queries or inverse bounds")
    nodes=summary["quadrature_nodes"]
    maximum=order*(1+512*(summary["maximum_inverse_iterations"]+1))
    if nodes%order or not order*513<=nodes<=maximum:
        raise ValueError("gauge quadrature work exceeds its bounded inverse envelope")
    if summary["evaluations"]!=3+summary["query_points"]+nodes:
        raise ValueError("gauge evaluations do not equal checks, queries and quadrature nodes")
    for key in ("maximum_inverse_residual","endpoint_integral_defect"):
        if not math.isfinite(summary[key]) or abs(summary[key])>INVERSE_TOLERANCE:
            raise ValueError("gauge receipt inverse residual exceeds its original tolerance")
    if summary["maximum_inverse_residual"]<0:
        raise ValueError("gauge maximum inverse residual must be nonnegative")
    if midpoint_rows is not None:
        iterations=np.asarray(midpoint_rows["inverse_iterations"],dtype=float)
        residual=np.asarray(midpoint_rows["inverse_residual"],dtype=float)
        if (iterations.shape!=(512,) or residual.shape!=(512,)
                or not np.all(np.isfinite(iterations)) or not np.all(np.isfinite(residual))
                or np.any(iterations!=np.floor(iterations)) or np.any(iterations<1)
                or np.any(iterations>32) or np.max(abs(residual))>INVERSE_TOLERANCE
                or summary["maximum_inverse_iterations"]!=int(np.max(iterations))
                or summary["maximum_inverse_residual"]!=float(np.max(abs(residual)))):
            raise ValueError("gauge receipt does not match its retained midpoint stopping data")


def _execution_replay(saved,reconstructed):
    """Keep historical accounting exact internally; fresh work is independent."""
    fixed=[]
    for cache in (saved,reconstructed):
        execution=cache["execution"]
        for key in ("evaluations","source_snapshot_points","gauge_evaluations","adapter_evaluations"):
            if type(execution.get(key)) is not int or execution[key]<0:
                raise ValueError("cache execution requires exact nonnegative integer "+key)
        if execution["source_snapshot_points"]!=2314:
            raise ValueError("cache source snapshot accounting changed")
        summaries=execution["gauge_orders"]
        if set(summaries)!={"64","32"}:
            raise ValueError("cache gauge-order receipt inventory changed")
        midpoint={key:np.asarray(value)[1::2] for key,value in cache["gauge_rows"].items()}
        for order in (64,32):
            _gauge_work(summaries[str(order)],order,midpoint if order==64 else None)
        if "diagnostics" in cache:
            diagnostic=cache["diagnostics"]
            for order in (64,32):
                if diagnostic["midpoint_gauge"][str(order)]["inverse_residual"]!=summaries[str(order)]["maximum_inverse_residual"]:
                    raise ValueError("gauge stopping diagnostic does not match its receipt")
            original=diagnostic["original_gauge"]["inverse_residual"]
            if not 0<=original<=INVERSE_TOLERANCE or diagnostic["combined_gauge"]["inverse_residual"]!=max(
                    original,summaries["64"]["maximum_inverse_residual"]):
                raise ValueError("combined gauge stopping diagnostic is inconsistent")
            difference=diagnostic["gauge64_minus32_maximum_absolute"]["inverse_residual"]
            if not 0<=difference<=2*INVERSE_TOLERANCE:
                raise ValueError("gauge stopping comparison exceeds its residual envelope")
        if execution["gauge_evaluations"]!=sum(s["evaluations"] for s in summaries.values()):
            raise ValueError("cache gauge evaluation total is inconsistent")
        if execution["evaluations"]!=sum(execution[key] for key in
                ("source_snapshot_points","gauge_evaluations","adapter_evaluations")):
            raise ValueError("cache total evaluations are inconsistent")
        comparable={key:value for key,value in execution.items()
                    if key not in {"evaluations","gauge_evaluations","gauge_orders"}}
        comparable["gauge_orders"]={label:{key:value for key,value in summary.items()
            if key not in _ADAPTIVE_GAUGE_FIELDS} for label,summary in summaries.items()}
        fixed.append(comparable)
    residual_allowance={"gauge_orders":{str(order):{"endpoint_integral_defect":2*INVERSE_TOLERANCE}
                                        for order in (64,32)}}
    if not _same_replayed(fixed[0],fixed[1],residual_allowance):
        raise ValueError("cache fixed execution declarations differ")


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


def make_refined_cache(base,budget,*,_validation=False):
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
        if _validation:
            _replay_rows(reproduced,old,label+" old source reproduction",
                         geometry=base["geometry"],alpha=data["alpha"],
                         source_knots=_source_knots(data),budget=budget)
        else:
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
    reconstructed=make_refined_cache(base,budget,_validation=True)
    # The historical513 row bytes and all cache keys remain exact. Only newly
    # recomputed floating arithmetic is compared numerically; nothing is replaced.
    try:
        if not np.array_equal(np.asarray(cache["x_nodes"]).view(np.uint64),
                              np.asarray(reconstructed["x_nodes"]).view(np.uint64)):
            raise ValueError("cache keys differ")
        if set(cache["source_rows"])!=set(reconstructed["source_rows"]):
            raise ValueError("source label inventory differs")
        for label in ("zero","finite"):
            old=base["initial_data"][label]["snapshot_rows"]
            table=cache["source_rows"][label]
            if not _same_copied_rows(table,old):
                raise ValueError("copied source rows differ")
            _replay_rows(table,reconstructed["source_rows"][label],label+" source replay",
                         geometry=base["geometry"],alpha=base["initial_data"][label]["alpha"],
                         source_knots=_source_knots(base["initial_data"][label]))
        old=base["gauge"]["rows"]
        table=cache["gauge_rows"]
        if not _same_copied_rows(table,old):
            raise ValueError("copied gauge rows differ")
        _replay_rows(table,reconstructed["gauge_rows"],"gauge replay",gauge=True)
        _execution_replay(cache,reconstructed)
        keys=set(reconstructed)-{"x_nodes","source_rows","gauge_rows","execution"}
        source_limits=_source_summary_allowances(reconstructed,base)
        for label,fields in source_limits["diagnostics"]["source"].items():
            for field,atol in fields.items():
                actual=cache["diagnostics"]["source"][label][field]
                expected=reconstructed["diagnostics"]["source"][label][field]
                if not _same_replayed(actual,expected,atol):
                    raise ValueError(f"source diagnostic {label}.{field} differs: "
                        f"difference={abs(actual-expected):.17g}, "
                        f"allowance={atol+2e-12*abs(expected):.17g}")
        source_limits["diagnostics"].update(_STOPPING_DIAGNOSTIC_ALLOWANCES)
        if set(cache)!=set(reconstructed) or not _same_replayed(
                {k:cache[k] for k in keys},{k:reconstructed[k] for k in keys},source_limits):
            raise ValueError("cache declarations or diagnostics differ")
    except (KeyError,TypeError,ValueError) as error:
        raise ValueError("refined cache differs from deterministic reconstruction: "+str(error)) from error
    budget.check()
    return dict(validation_evaluations=budget.evaluations-start,
                saved_cache_evaluations=cache["execution"]["evaluations"],
                replayed_cache_evaluations=reconstructed["execution"]["evaluations"],source_integrations=0,
                regional_marches=0,time_evolutions=0)
