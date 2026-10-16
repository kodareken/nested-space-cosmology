"""Build computational cap histories from existing wall or prepared-field jets.

The caller supplies values and material rates evaluated at the same saved
stages. This adapter does not prepare a source, reconstruct a physical stage,
or evolve a field. Its second jets are interpolation curvature; they are not
an independently measured or computed physical second derivative.
The prepared-field adapter instead uses the existing continuum-integral
second jet and preserves the analytic quartic onset explicitly.
"""
from collections.abc import Mapping
from fractions import Fraction
import math

import numpy as np

from .nsc_atlas_characteristic_background import CharacteristicProfile, MAX_POINTS


CURVATURE_METHOD = "numpy.gradient(F_prime, arrival_labels, edge_order=2)"
SCOPE = ("computational C2 quintic background from supplied zero-case wall values "
         "and same-stage first jets; finite boundary discrepancies remain evolved; "
         "interpolation curvature is not physical second-jet evidence")


def _rows(rows, name):
    if not isinstance(rows, Mapping) or not {"T", "r", "c"}.issubset(rows):
        raise ValueError(name+" requires T,r,c arrays")
    result = {}
    for key in ("T", "r", "c"):
        value = rows[key]
        if isinstance(value, (list, tuple)) and any(
                isinstance(v, (bool, np.bool_)) for v in np.asarray(value, dtype=object).flat):
            raise ValueError(name+" contains boolean data")
        raw = np.asarray(value)
        if (raw.ndim != 1 or raw.dtype.kind not in "iuf"
                or not 3 <= raw.size <= MAX_POINTS):
            raise ValueError(name+" requires 3..513 real one-dimensional rows")
        result[key] = np.array(raw, dtype=float, copy=True)
        if not np.all(np.isfinite(result[key])):
            raise ValueError(name+" must be finite")
    if len({a.shape for a in result.values()}) != 1:
        raise ValueError(name+" array shapes differ")
    return result


def build_background_profile(walls, wall_rates, *, alpha, budget):
    """Construct one bounded profile from wall arrays with matching row order.

    Both mappings contain ``T,r,c`` vectors. ``wall_rates`` are derivatives
    with respect to the same common evolution parameter, not derivatives of
    a saved-value interpolant. At alpha=0, before any reflected return,
    zeta=T+atanh(r), F=c and F'=c_x/zeta_x. Second jets are the declared
    second-order gradient of these supplied first jets in zeta. A single
    quintic Hermite interpolant then owns every runtime derivative.

    The first value and both resulting derivatives must be exactly quiet.
    Nothing is extrapolated beyond the last label. Later background queries
    must satisfy T+atanh(r)<=labels[-1]; the returned profile enforces this.
    Lower extension is only the certified zero branch. The caller retains
    responsibility for authenticating the supplied case and its same-stage
    provenance; numerical arrays alone cannot establish that provenance.

    Charges are six input scalar arrays, three derived arrays and one knot
    jet evaluation per row, all on the caller's shared budget.
    """
    if (isinstance(alpha, (bool, np.bool_)) or not np.isscalar(alpha)
            or not isinstance(alpha, (int, float, np.integer, np.floating))
            or not math.isfinite(float(alpha)) or float(alpha) != 0.):
        raise ValueError("background input must be the alpha=0 case")
    if not callable(getattr(budget, "check", None)):
        raise ValueError("one shared cumulative budget.check(count) is required")
    budget.check()
    before = budget.evaluations if hasattr(budget, "evaluations") else None
    walls, rates = _rows(walls, "walls"), _rows(wall_rates, "wall_rates")
    if walls["T"].shape != rates["T"].shape:
        raise ValueError("wall values and rates must have the same row count")
    n = len(walls["T"])
    budget.check(6*n)
    r = walls["r"]
    if np.any(r <= 0) or np.any(r >= 1):
        raise ValueError("wall radii must lie strictly between zero and one")
    with np.errstate(over="raise", divide="raise", invalid="raise", under="ignore"):
        try:
            labels = walls["T"]+np.arctanh(r)
            arrival_rate = rates["T"]+rates["r"]/(1-r*r)
            if (not np.all(np.isfinite(labels)) or not np.all(np.diff(labels) > 0)
                    or not np.all(np.isfinite(arrival_rate)) or np.any(arrival_rate <= 0)):
                raise ValueError("arrival labels must increase and their rates must be positive")
            first = rates["c"]/arrival_rate
            second = np.gradient(first, labels, edge_order=2)
        except (FloatingPointError, OverflowError) as error:
            raise ValueError("nonfinite profile chain or computational curvature") from error
    budget.check(3*n)
    values = walls["c"].copy()
    if not all(np.all(np.isfinite(a)) for a in (first, second)):
        raise ValueError("nonfinite profile derivatives")
    if values[0] != 0 or first[0] != 0 or second[0] != 0:
        raise ValueError("first value, first jet and computational curvature must be exactly quiet")
    no_return_margin = float(labels[0]-np.max(walls["T"]))
    if no_return_margin <= 0:
        raise ValueError("input must remain before the first central arrival")
    profile = CharacteristicProfile(labels, values, first, second)
    knots = profile.jet(labels, budget=budget)
    diagnostics = {
        "no_return_margin": no_return_margin,
        "minimum_arrival_rate": float(np.min(arrival_rate)),
        "minimum_label_spacing": float(np.min(np.diff(labels))),
        "max_abs_material_chain_defect": float(np.max(abs(first*arrival_rate-rates["c"]))),
        "max_abs_knot_value_defect": float(np.max(abs(knots["F"]-values))),
        "max_abs_knot_first_defect": float(np.max(abs(knots["F_prime"]-first))),
        "max_abs_knot_second_defect": float(np.max(abs(knots["F_second"]-second))),
        "quiet_anchor": {"value": float(values[0]), "first": float(first[0]),
                         "second": float(second[0])},
    }
    budget.check()
    return {
        "profile": profile, "labels": labels, "values": values, "first": first,
        "second": second, "walls": walls, "wall_rates": rates,
        "arrival_rate": arrival_rate, "label_range": [float(labels[0]), float(labels[-1])],
        "alpha": 0., "curvature_method": CURVATURE_METHOD, "scope": SCOPE,
        "diagnostics": diagnostics,
        "charged_evaluations": None if before is None else budget.evaluations-before,
    }


def build_prepared_background(cap, *, quiet_tolerances=(0.,0.,0.)):
    """Retain the existing prepared field, with its mathematical quartic front.

    Existing field cuts and geometry knots supply the labels; there is no
    solve, knot scan or waveform differentiation. Values/first jets and the
    continuum-integral second jet come from CapLeadingField.profile. The
    resulting Hermite profile owns its own runtime derivatives. Their
    difference from the prepared field needs a separate error comparison.

    The analytic quiet branch replaces supplied quiet jets only within the
    caller's three declared numerical tolerances, which default to zero.
    Their measured discrepancies remain in the result. The rounded source
    front is recorded relative to the exact root of the binary endpoints.
    cap.budget owns every query and all adapter accounting.
    """
    from .nsc_atlas_cap_field import CapLeadingField
    from .nsc_atlas_cap_background_jets import RetainedCapBackground
    from .nsc_atlas_characteristic_background import _array,_budget
    if not isinstance(cap,CapLeadingField):
        raise ValueError('prepared background requires the existing CapLeadingField')
    budget=_budget(cap.budget)
    before=budget.kernel_values
    tolerances=_array(quiet_tolerances,'quiet tolerances')
    if tolerances.shape!=(3,) or np.any(tolerances<0):
        raise ValueError('three nonnegative quiet tolerances are required')
    if cap.quiet_certificate.get('quiet_cap') is not True:
        raise ValueError('prepared background requires the quiet-cap certificate')
    cuts=_array(cap.prepared_field.dense_data['cuts'],'prepared field cuts')
    geometry=_array(cap.geometry['arrival'].x,'arrival geometry knots')
    if (cuts.ndim!=1 or len(cuts)<2 or cuts[0]!=0. or cuts[-1]!=1.
            or np.any(np.diff(cuts)<=0.) or geometry.ndim!=1
            or np.any(np.diff(geometry)<=0.)):
        raise ValueError('prepared field cuts and arrival knots must be ordered')
    Va,Vb=cap.Va,cap.Vb;dV=cap.dV
    if not np.isfinite([Va,Vb,dV]).all() or not Va<Vb or dV!=Vb-Va:
        raise ValueError('prepared source endpoints are inconsistent')
    edges=Va+dV*np.array([.25,.75])
    native=geometry[(geometry>=Va)&(geometry<=Vb)]
    V=np.unique(np.r_[Va,Vb,native,np.clip(Va+dV*cuts,Va,Vb),edges])
    if not 2<=len(V)<=MAX_POINTS:
        raise ValueError('native prepared-history union exceeds 513 labels')
    labels=np.asarray(cap.geometry['arrival'](V),dtype=float)
    if (labels.shape!=V.shape or not np.all(np.isfinite(labels)) or np.any(np.diff(labels)<=0.)
            or labels[0]!=cap.domain['qa'] or labels[-1]!=cap.domain['qb']):
        raise ValueError('prepared-history arrival labels are unresolved or leave the finite domain')
    index=int(np.searchsorted(V,edges[0]));qV=float(cap.geometry['arrival'](edges[0],1))
    if not np.isfinite(qV) or qV<=0.: raise ValueError('source-front arrival rate must be positive')
    try:
        quartic=float((8/(dV*qV))**4)
    except (OverflowError,ZeroDivisionError) as error:
        raise ValueError('source-front quartic coefficient is unresolved') from error
    if not np.isfinite(quartic) or quartic==0.:
        raise ValueError('source-front quartic coefficient is unresolved')
    budget.check(8*len(V))
    sample=cap.profile(labels,second=True)
    raw=np.array([sample[k] for k in ('f','f_prime','f_second')],dtype=float,copy=True)
    if raw.shape!=(3,len(V)) or not np.all(np.isfinite(raw)):
        raise ValueError('prepared history returned invalid value/first/second jets')
    quiet=np.max(abs(raw[:,:index+1]),axis=1)
    if np.any(quiet>tolerances):
        raise ValueError('prepared quiet-jet discrepancy exceeds the declared tolerances: '
                         f'{quiet.tolist()} versus {tolerances.tolist()}')
    supplied=raw.copy();supplied[:,:index+1]=0.
    background=RetainedCapBackground(labels,*supplied,
        onset=dict(index=index,quartic_coefficient=quartic))
    jet=background.profile.jet(labels,budget=budget)
    defects=np.array([np.max(abs(jet[k]-supplied[i])) for i,k in enumerate(('F','F_prime','F_second'))])
    exact_front=(3*Fraction(float(Va))+Fraction(float(Vb)))/4
    budget.check(6*len(V))
    return dict(background=background,source_null_labels=V.copy(),
        physical_input_domain=dict(cap.domain),
        onset=dict(index=index,label=float(labels[index]),source_V=float(edges[0]),arrival_rate=qV,
            quartic_coefficient=quartic,source_V_roundoff=float(Fraction(float(edges[0]))-exact_front)),
        diagnostics=dict(quiet_jet_discrepancies=quiet,quiet_tolerances=tolerances.copy(),
            knot_jet_discrepancies=defects,maximum_inverse_bracket_width=float(cap.max_inverse_bracket_width),
            original_field_cuts=len(cuts),original_geometry_knots=len(geometry),retained_knots=len(V)),
        curvature_method='CapLeadingField.profile(second=True): differentiated continuum integral and arrival chain',
        charged_evaluations=budget.kernel_values-before,numerical_error_enclosure=False,
        scope='prepared leading field retained as a computational C2 history with its analytic quartic onset; no finite-amplitude evolution or interpolation error enclosure')
