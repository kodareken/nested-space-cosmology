"""One fixed-coefficient SSPRK3 step and its continuous Hermite residual.

The stage owns source lookup, geometry and physical rates. This adapter owns
only time combinations. A residual at selected points is not an error bound.
Historical Ip/Iq drivers and their authenticated source bytes stay intact.
"""
import copy

import numpy as np

from .nsc_atlas_background_ssprk3 import _times, TIME_SCHEME
from .nsc_atlas_coupled_evolution import _combine, _finite, CoupledEvolutionError
from .nsc_atlas_fixed_cap_stage import fixed_state
from .nsc_atlas_stage_ssprk3 import _shared_budget


def _size(state):
    return sum(np.size(v) for group in state.values() for v in group.values())


def _sum(a,b,wa,wb,budget):
    budget.check(3*_size(a))
    return _combine(a,b,wa,wb)


def advance_fixed_ssprk3(initial,x_left,x_mid,x_right,*,basis,exterior_cells,stage,budget):
    """Three RK stages plus the accepted endpoint RHS, with no retries.

    ``stage(state, *, x, budget)`` returns the common stage mapping including
    ``rates``. The endpoint evaluation supplies the right Hermite slope; the
    predictor's right-time evaluation cannot supply that slope. Arrays and
    full diagnostics are returned in memory for the caller to summarize.
    No cap Courant or stability claim is inferred from exterior speeds.
    """
    phase,calls='validation',0
    state,attempt=copy.deepcopy(initial),copy.deepcopy(initial)
    rows,states,rates=[],[],[]
    try:
        _shared_budget(budget)
        if not callable(stage): raise ValueError('stage must be callable')
        left,mid,right=_times(x_left,x_mid,x_right)
        h=right-left
        state=fixed_state(initial,exterior_cells,basis)

        def evaluate(value,x,label):
            nonlocal phase,calls,attempt
            phase=label
            attempt=fixed_state(value,exterior_cells,basis)
            # A callback cannot change a retained state or an earlier slope.
            snapshot=copy.deepcopy(attempt)
            calls+=1
            row=stage(copy.deepcopy(attempt),x=x,budget=budget)
            rate=fixed_state(row['rates'],exterior_cells,basis)
            budget.check(_size(rate))
            rows.append(row);states.append(snapshot);rates.append(rate)
            return rate

        first=evaluate(state,left,'left_stage')
        phase='predictor'
        predictor=_sum(state,first,1.,h,budget)
        second=evaluate(predictor,right,'right_predictor_stage')
        phase='midpoint_combination'
        midpoint=_sum(state,_sum(predictor,second,1.,h,budget),.75,.25,budget)
        third=evaluate(midpoint,mid,'midpoint_stage')
        phase='accepted_combination'
        accepted=fixed_state(_sum(state,_sum(midpoint,third,1.,h,budget),1/3,2/3,budget),
                             exterior_cells,basis)
        evaluate(accepted,right,'accepted_endpoint_stage')
        budget.check()
    except (ValueError,RuntimeError,FloatingPointError,OverflowError,KeyError,TypeError) as error:
        raise CoupledEvolutionError(dict(phase=phase,stage_calls=calls,
            completed_stages=len(rows),charged_evaluations=getattr(budget,'evaluations',None),
            reason=str(error),accepted_state=copy.deepcopy(state),attempt_state=copy.deepcopy(attempt))) from error
    return dict(initial_state=state,accepted_state=accepted,stage_states=states,stages=rows,
                rates=rates,stage_times=[left,right,mid,right],x_left=left,x_mid=mid,x_right=right,
                stage_calls=calls,time_scheme=TIME_SCHEME,retries=0)


def fixed_hermite(step,x,*,basis,exterior_cells,budget):
    """Evaluate the accepted-endpoint cubic and its analytic derivative.

    The polynomial is built from increments to avoid subtracting large
    absolute states. Exact endpoint queries return the copied endpoint and
    its slope. No extrapolation is admitted.
    """
    _shared_budget(budget)
    left,_,right=_times(step['x_left'],step['x_mid'],step['x_right'])
    x=float(_finite(x,'x',()))
    if not left<=x<=right: raise ValueError('Hermite query leaves the accepted interval')
    y0=fixed_state(step['initial_state'],exterior_cells,basis)
    y1=fixed_state(step['accepted_state'],exterior_cells,basis)
    f0=fixed_state(step['rates'][0],exterior_cells,basis)
    f1=fixed_state(step['rates'][3],exterior_cells,basis)
    if x==left: return dict(state=y0,derivative=f0,x=x)
    if x==right: return dict(state=y1,derivative=f1,x=x)
    h=right-left;t=(x-left)/h
    delta=_sum(y1,y0,1.,-1.,budget)
    a=_sum(delta,_sum(f0,f1,2.,1.,budget),3.,-h,budget)
    b=_sum(delta,_sum(f0,f1,1.,1.,budget),-2.,h,budget)
    value=_sum(_sum(y0,f0,1.,h*t,budget),_sum(a,b,1.,t,budget),1.,t*t,budget)
    derivative=_sum(f0,_sum(a,b,2.,3*t,budget),1.,t/h,budget)
    return dict(state=fixed_state(value,exterior_cells,basis),
                derivative=fixed_state(derivative,exterior_cells,basis),x=x)


def fixed_time_residual(step,x,*,basis,exterior_cells,stage,budget):
    """One independent RHS call on the continuous path, with raw residual.

    Component residuals have their component's derivative units. They are
    not combined into an arbitrary mixed-unit norm or a uniform enclosure.
    Physical junction/constraint diagnostics remain in ``stage``.
    """
    if not callable(stage): raise ValueError('stage must be callable')
    path=fixed_hermite(step,x,basis=basis,exterior_cells=exterior_cells,budget=budget)
    row=stage(copy.deepcopy(path['state']),x=path['x'],budget=budget)
    rate=fixed_state(row['rates'],exterior_cells,basis)
    residual=_sum(path['derivative'],rate,1.,-1.,budget)
    budget.check()
    return dict(x=path['x'],state=path['state'],derivative=path['derivative'],
                residual=residual,stage=row,stage_calls=1,uniform_error_bound=False)
