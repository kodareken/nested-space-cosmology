"""Generic retained tail wrapper around the frozen group22 analytic proof."""
from functools import lru_cache
import json
import mpmath as mp

from .nsc_incoming_tail_quadrature_bound import CompactTailKernel, analytic_bloch, endpoint_certificate
from .nsc_incoming_source_quadrature_bound import certify_gauss_brackets, pack, unpack, _conj_coeff
from .nsc_incoming_projector_energy_bound import incoming_riccati_interval_coefficients
from .nsc_incoming_middle_bound import _parameters
from .nsc_incoming_vacuum_tail_bound import _precision, _lo, _hi, _range, _up_float


def endpoint(channel):
    group = channel.get('index')
    if isinstance(group,bool) or group not in range(1,33): raise ValueError('retained non-LLL group required')
    return 320. if group in (10,11,12,31,32) else 160.


def signs(channel): return (1,) if not channel['angular_eigenvalue'] else (1,-1)


@lru_cache(maxsize=2)
def _cached_nodes(serialized_rule):
    rule=json.loads(serialized_rule)
    if (rule['points'],rule['precision'])!=(48,80): raise ValueError('frozen48/80 rule required')
    with _precision(80):
        nodes=certify_gauss_brackets(rule['roots'],48)
        if [pack(w) for _,w in nodes]!=rule['weights']: raise ValueError('saved Gauss weight changed')
        return nodes


class RetainedTailKernel(CompactTailKernel):
    """Real generic initializer; inherited analytic kernel remains unchanged."""
    def __init__(self,channel,sign,positive_coefficients):
        endpoint(channel)
        if sign not in signs(channel) or len(positive_coefficients) != 16:
            raise ValueError('actual sign and order16 coefficients required')
        self.a,self.r,self.mass,angular,self.factor = _parameters(channel)
        self.factor /= len(signs(channel)); self.ell = sign*angular
        self.c = positive_coefficients if sign == 1 else [(-1)**(j+1)*_conj_coeff(v) for j,v in enumerate(positive_coefficients)]
        self.cb = [_conj_coeff(v) for v in self.c]
        self.bloch = analytic_bloch()


def prepare_tail(channel,rule,*,progress=None):
    """One new numerical integral; group22 must reuse its authenticated pilot."""
    lower = endpoint(channel)
    if channel['index'] == 22: raise ValueError('reuse the completed group22 pilot payload')
    if (rule['points'],rule['precision']) != (48,80): raise ValueError('frozen48/80 rule required')
    with _precision(80):
        nodes = _cached_nodes(json.dumps(rule,sort_keys=True))
        coefficients = incoming_riccati_interval_coefficients(channel,order=16)
        h = mp.iv.mpf(1)/(2*lower); radius=2*h; rows=[]
        for sign in signs(channel):
            kernel=RetainedTailKernel(channel,sign,coefficients); margins=kernel.margins(3*h)
            maxima=[mp.mpf(0)]*4
            for arc in range(32):
                theta=2*mp.iv.pi*_range(mp.mpf(arc)/32,mp.mpf(arc+1)/32)
                x=mp.iv.mpc(h+radius*mp.iv.cos(theta),radius*mp.iv.sin(theta))
                for j,v in enumerate(kernel(x)): maxima[j]=max(maxima[j],_hi(abs(v)))
            if not all(mp.isfinite(v) for v in maxima): raise ArithmeticError('nonfinite fixed-resolution tail circle')
            samples=[]
            for node,_ in nodes:
                values=kernel(h+h*node)
                if any(not _lo(v.imag)<=0<=_hi(v.imag) for v in values): raise ArithmeticError('tail sample not real')
                samples.append([pack(v.real) for v in values])
            rows.append({'sign':sign,'physical_factor':pack(kernel.factor),'margins':margins,
                         'circle_absolute_upper':[pack(mp.iv.mpf(v)) for v in maxima],'samples':samples})
            if progress: progress(sign)
        return {'group':channel['index'],'lower_energy':lower,'riccati_order':16,'reference_order':4,
                'precision':80,'points':48,'circle_arcs':32,'halfwidth':pack(h),'radius':pack(radius),
                'rule':rule,'signs':rows,'endpoint_certificate':endpoint_certificate(),
                'analytic_adapter':analytic_bloch().adapter_metadata}


def replay_tail(payload,channel):
    lower=endpoint(channel)
    if (payload['group'],payload['lower_energy'],payload['riccati_order'],payload['reference_order'],payload['points'],payload['precision']) != (channel['index'],lower,16,4,48,80):
        raise ValueError('matching actual retained group/order/endpoint required')
    if payload['endpoint_certificate'] != endpoint_certificate(): raise ValueError('exact endpoint certificate changed')
    with _precision(80):
        nodes=_cached_nodes(json.dumps(payload['rule'],sort_keys=True))
        if [pack(w) for _,w in nodes] != payload['rule']['weights']: raise ValueError('saved Gauss weight changed')
        h,radius=unpack(payload['halfwidth']),unpack(payload['radius'])
        if pack(h)!=pack(mp.iv.mpf(1)/(2*lower)) or pack(radius)!=pack(2*h): raise ValueError('compactified interval changed')
        if [r['sign'] for r in payload['signs']] != list(signs(channel)): raise ValueError('actual sign inventory changed')
        total=[mp.iv.mpf(0) for _ in range(4)]; errors=[mp.iv.mpf(0) for _ in range(4)]
        for row in payload['signs']:
            if len(row['samples'])!=48 or any(_hi(unpack(v))>=1 for v in row['margins'].values()): raise ValueError('sample or analytic margin failure')
            factor=unpack(row['physical_factor'])
            expected_factor=_parameters(channel)[4]/len(signs(channel))
            if pack(factor)!=pack(expected_factor): raise ValueError('physical signed factor changed')
            for j in range(4):
                quadrature=h*sum((w*unpack(sample[j]) for (_,w),sample in zip(nodes,row['samples'])),mp.iv.mpf(0))
                M=unpack(row['circle_absolute_upper'][j]); error=4*h*M*(h/radius)**96/(1-h/radius)
                total[j]+=factor*quadrature; errors[j]+=factor*error
        return {'integral_interval':[pack(v+_range(-_hi(e),_hi(e))) for v,e in zip(total,errors)],
                'quadrature_interval':[pack(v) for v in total],
                'analytic_quadrature_error_upper':[pack(v) for v in errors]}


def thermal_tail(channel,config):
    lower=endpoint(channel)
    with _precision(80):
        a,r,m,ell,factor=_parameters(channel)
        kappa,omega=mp.iv.mpf(config['surface_gravity']),mp.iv.mpf(config['omega'])
        m0,m1=mp.iv.mpf(0),mp.iv.mpf(0)
        for coefficient,alpha in ((2,mp.iv.pi/kappa),(1,2*mp.iv.pi/(omega*kappa))):
            exponential=coefficient*mp.iv.exp(-alpha*lower)
            m0+=exponential/alpha; m1+=exponential*(lower/alpha+1/alpha**2)
        M=mp.iv.sqrt(m*m+(ell/r)**2)
        values=[factor*(m1/a+M*m0),factor*m1/a,factor*m1/a,factor*ell*m0/(2*r)]
        zero=[i for i,v in enumerate(values) if _lo(v)==_hi(v)==0]
        if zero != ([3] if channel['angular_eigenvalue']==0 else []) or any(_lo(v)<=0 for i,v in enumerate(values) if i not in zero):
            raise ArithmeticError('thermal zero/nonzero inventory differs')
        return {'intervals':[pack(v) for v in values], 'upper_display':[mp.nstr(_hi(v),25) for v in values],
                'display_is_not_proof':'lossless interval endpoints own outward rounding',
                'source_law':'2*exp(-pi*E/kappa)+exp(-2*pi*E/(Omega*kappa))', 'exact_zero_components':zero}


def summarize_tail(integral,thermal,archived,physical,channel):
    with _precision(80):
        values=[unpack(v) for v in integral['integral_interval']]
        errors=[mp.iv.mpf(_hi(abs(v-mp.iv.mpf(old)))) for v,old in zip(values,archived)]
        a,r,_,_,_=_parameters(channel)
        vacuum=[mp.iv.mpf(v) for v in physical]
        th=[unpack(v) for v in thermal['intervals']]
        numerical=[4*mp.iv.pi*a*r*r*errors[0],4*mp.iv.pi*a*a*r*r*errors[2]]
        complete=[4*mp.iv.pi*a*r*r*(errors[0]+vacuum[0]+th[0]),4*mp.iv.pi*a*a*r*r*(errors[2]+vacuum[2]+th[2])]
        up=lambda v:0. if _hi(v)==0 else _up_float(v)
        return {'group':channel['index'],'lower_energy':endpoint(channel),'actual_signs':list(signs(channel)),
                'archived_order13_source':archived,'archived_source_unchanged':True,
                'numerical_stress_error_upper_binary':[pack(v) for v in errors],
                'numerical_action_error_upper':[up(v) for v in numerical],
                'combined_action_error_upper':[up(v) for v in complete],
                'combined_action_upper_binary':[pack(mp.iv.mpf(_hi(v))) for v in complete],
                'thermal_tail':thermal,'physical_stress_error_upper':physical,
                'status':'PASS' if max(_hi(v) for v in complete)<=mp.mpf('3e-11') else 'OPEN',
                'stationarity_tolerance':3e-11}
