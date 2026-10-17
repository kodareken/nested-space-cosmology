"""Exact initial-reference enclosure around a retained floating seed.

The inverse residual is integrated with balls at one prescribed seed, then
divided by a proved derivative floor. Clock integrals and the quiet-history
margin follow on that enclosure. This certifies an initial event, not every
floating gauge evaluation or the subsequent sourced evolution.
"""
from fractions import Fraction as Q

from flint import acb, arb, ctx

from .nsc_atlas_background_cap import _budget
from .nsc_atlas_reference_gauge_bound import _ball, _integral, _record, _state
from .nsc_atlas_source_cell_bound import _coordinate


def _from_record(row):
    def endpoint(name):
        value=row[name]
        if type(value['exponent']) is not int:
            raise ValueError('invalid dyadic exponent')
        return Q(int(value['mantissa']))*Q(2)**value['exponent']
    lo,hi=endpoint('lower'),endpoint('upper')
    if lo>hi: raise ValueError('reversed enclosure')
    return _ball(lo).union(_ball(hi))


def inverse_seed_enclosure(seed,residual,derivative_floor):
    """Mean-value radius bound; a common monotone bracket is caller input."""
    if not residual.is_finite() or not derivative_floor>0:
        raise ValueError('finite residual and positive derivative floor required')
    error=(abs(residual)/derivative_floor).upper()
    return seed+(-error).union(error),error


def _clock_integral(seed,kind,budget):
    def integrand(R,analytic):
        budget.check(64)
        bm,_,q=_state(R)
        value=1/q.sqrt(analytic=analytic)
        return value if kind=='s' else bm*value/(1-R*R)
    if kind not in ('s','T'): raise ValueError('unknown reference clock')
    result=acb.integral(integrand,acb(_ball(Q(9,20))),acb(_ball(seed)),
                        rel_tol=arb(2)**-80,abs_tol=arb(2)**-80,
                        deg_limit=64,eval_limit=512,depth_limit=20)
    budget.check()
    if not result.is_finite() or not result.imag.contains(0):
        raise ValueError('clock integral did not enclose a finite real value')
    return result.real


def reference_initial_bound(gauge_certificate,*,seed,retained_wall,history_offset,
                            history_onset,cap_knots,source_geometry,budget):
    """Use the authenticated full-cell bracket to certify its left event."""
    budget=_budget(budget);before=budget.evaluations
    for value in (seed,history_offset,history_onset,*retained_wall.values(),*cap_knots):
        _ball(value)
    seed=Q(seed);inputs=gauge_certificate['inputs']
    left=Q(inputs['left']);rlo=Q(inputs['radius_lower']);rhi=Q(inputs['radius_upper'])
    if not rlo<seed<rhi or set(retained_wall)!={'r','L','T','s','c'}:
        raise ValueError('seed and retained wall must describe the original initial event')
    if len(cap_knots)!=3 or not Q(cap_knots[0])==0<Q(cap_knots[1])<Q(cap_knots[2]):
        raise ValueError('require the original two-cell fixed cap basis')
    budget.check(3072)
    with ctx.workprec(160):
        if not all(_from_record(row)>0 for row in gauge_certificate['bracket_margins']):
            raise ValueError('whole-cell inverse bracket was not proved')
        # These are the same exact reference normalization and real path.
        if gauge_certificate['normalization']!='exact rational corner x0=2/5,R0=9/20,a=4/5':
            raise ValueError('reference normalization differs')
        x,a,Rhat=_ball(left),_ball(Q(4,5)),_ball(seed)
        source_quiet=-1-_coordinate(x,source_geometry)[0]
        if not source_quiet>0: raise ValueError('original source is not strictly before its exact front')
        bracket=_ball(rlo).union(_ball(rhi))
        bm,bp,q=_state(bracket);v=q.sqrt()
        if not q>0 or not a-bracket>0: raise ValueError('reference derivative path not regular')
        derivative=(bracket/(a-bracket)*(1+bp/v)).lower()
        target=_from_record(gauge_certificate['endpoint_targets'][0])
        residual=_integral(seed,budget)-target
        R,error=inverse_seed_enclosure(Rhat,residual,derivative)
        if not R>_ball(rlo) or not R<_ball(rhi):
            raise ValueError('initial correction leaves its certified reference bracket')
        clocks={}
        for kind,rate in (('s',1/v),('T',bm/((1-bracket*bracket)*v))):
            integral=_clock_integral(seed,kind,budget)
            error_clock=(abs(rate).upper()*error).upper()
            clocks[kind]=integral+(-error_clock).union(error_clock)
        Z=R.atanh();C=x/R*((x-R)/a).exp();L=C.log()
        physical=dict(r=R,L=L,s=clocks['s'],T=clocks['T'],c=arb(0))
        margin=_ball(history_onset)-(clocks['T']+_ball(history_offset)+Z)
        cut=(Z-_ball(cap_knots[1]))/(_ball(cap_knots[2])-_ball(cap_knots[1]))
        if not margin>0: raise ValueError('retained history is not proved quiet on the initial cap')
        if not cut>_ball(Q(1,4)) or not cut<_ball(Q(3,4)):
            raise ValueError('corrected wall exits the admitted fixed-basis collar')
        budget.check()
        return dict(precision_bits=160,seed=str(seed),left=str(left),
                    integral_residual=_record(residual),inverse_derivative_floor=_record(derivative),
                    seed_radius_error=_record(error),physical={k:_record(v) for k,v in physical.items()},
                    retained_wall_error={k:_record(abs(v-_ball(retained_wall[k]))) for k,v in physical.items()},
                    quiet_history_margin=_record(margin),fixed_cap_cut=_record(cut),
                    source_front_margin=_record(source_quiet),
                    history_offset=str(history_offset),history_onset=str(history_onset),
                    charged_values=budget.evaluations-before,
                    scope='exact initial reference and retained-wall differences; zero-departure compatibility requires the authenticated quiet profile; later residuals remain external')
