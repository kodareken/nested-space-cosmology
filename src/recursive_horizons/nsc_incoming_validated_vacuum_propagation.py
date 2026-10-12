"""Fixed-E unitary propagation with directed continuous Taylor residuals.

Degree36, real step<=1/32, analytic disk radius1/8. Interval arithmetic
encloses generator coefficients; midpoint columns carry a separate error
radius, so propagated interval states never wrap around the unit sphere.
"""
from functools import lru_cache
import json
import math
from pathlib import Path

import mpmath as mp

from .nsc_incoming_matched_horizon_initializer import _geometry
from .nsc_incoming_source_quadrature_bound import pack, unpack, float_enclosure
from .nsc_incoming_vacuum_tail_bound import _precision, _lo, _hi, _range
from .nsc_incoming_state_moments import incoming_group_factor


DEGREE, NODES, PRECISION = 36, 128, 80
RADIUS, SAMPLE_RADIUS, MAX_STEP = mp.mpf(1)/8, mp.mpf(1)/16, mp.mpf(1)/32


def point(value):
    return list(mp.mpf(value)._mpf_)


def unpoint(value):
    return mp.make_mpf(tuple(value))


def cpoint(value):
    return [point(value.real), point(value.imag)]


def uncpoint(value):
    return mp.mpc(unpoint(value[0]), unpoint(value[1]))


def _mid(value):
    return (_lo(value)+_hi(value))/2


def _cmid(value):
    return mp.mpc(_mid(value.real), _mid(value.imag))


def _ivc(value):
    if hasattr(value, '_mpci_'):
        return value
    return mp.iv.mpc(value.real, value.imag)


def _norm(vector):
    return mp.iv.sqrt(sum((abs(value)**2 for value in vector), mp.iv.mpf(0)))


def _upper_float(value):
    return float_enclosure(value)[1] if _hi(value) else 0.


def _sinc_minus_one(z):
    """Whole complex rectangle enclosure of sinc(z)-1, including its tail."""
    z2 = z*z
    value = mp.iv.mpc(0)
    for n in range(18, 0, -1):
        value = (value+mp.iv.mpf((-1)**n)/math.factorial(2*n+1))*z2
    a = mp.iv.mpf(_hi(abs(z)))
    ratio = a*a/(40*41)
    if _hi(ratio) >= 1:
        raise ArithmeticError('fixed sinc-series remainder is outside its certified disk')
    remainder = a**38/math.factorial(39)/(1-ratio)
    error = _range(-_hi(remainder), _hi(remainder))
    return value+mp.iv.mpc(error, error)


class ExactProfileGenerator:
    """Analytic y=log(delta) generator of the same fixed canonical vacuum."""

    def __init__(self, channel, config, sign, energy=24):
        if sign not in (1, -1) or energy != 24:
            raise ValueError('fixed E24 and actual angular sign required')
        self.g = _geometry(channel, config)
        self.channel, self.sign, self.energy = dict(channel), sign, energy
        # Explicitly the frozen binary channel labels, not an unannounced
        # substitution of symbolic pi/sqrt values in the propagated operator.
        self.mass = mp.iv.mpf(channel['compact_mass'])
        self.angular = mp.iv.mpf(sign*channel['angular_eigenvalue'])
        self.qh = self.g['qh']
        self.b0 = 2*self.g['kg']
        s, c = mp.iv.sin(2*self.qh), mp.iv.cos(2*self.qh)
        self.linear, self.sinc_coefficient = -3*s-c, -3*c+s

    def values(self, y, *, proof=False):
        delta = mp.iv.exp(y)
        s1, s2 = _sinc_minus_one(delta), _sinc_minus_one(2*delta)
        # Exact W(qh-delta)/delta, using W(qh)=0. Stable at tiny delta.
        B = self.b0+self.linear*delta*(1+s1)**2+self.sinc_coefficient*s2
        sine = mp.iv.sin(self.qh-delta)
        if _lo(B.real) <= 0 or _lo(sine.real) <= 0:
            raise ArithmeticError('analytic geometry disk has no certified nonzero/branch margin')
        radius = 1/sine
        factor = mp.iv.exp(y/2)/(B**mp.iv.mpf('0.5'))
        D = self.energy/B
        plus = factor*(-self.mass*radius+mp.iv.j*self.angular)
        minus = factor*(-self.mass*radius-mp.iv.j*self.angular)
        values = (D, plus, minus)
        if proof:
            norm = mp.iv.sqrt(2*abs(D)**2+abs(plus)**2+abs(minus)**2)
            return values, {'B_real_lower': float_enclosure(B.real)[0],
                            'sin_real_lower': float_enclosure(sine.real)[0],
                            'generator_norm_upper': _upper_float(norm)}
        return values

    def disk(self, center):
        y = mp.iv.mpc(_range(center-RADIUS, center+RADIUS), _range(-RADIUS, RADIUS))
        return self.values(y, proof=True)[1]


@lru_cache(maxsize=2)
def _roots(precision):
    if mp.iv.dps != precision:
        raise ValueError('Fourier-root precision differs from active directed context')
    return tuple(mp.iv.mpc(mp.iv.cos(2*mp.iv.pi*j/NODES), mp.iv.sin(2*mp.iv.pi*j/NODES))
                 for j in range(NODES))


def _fft(values):
    """Directed radix-two DFT, negative Fourier sign."""
    n = len(values)
    if n != NODES:
        raise ValueError('the fixed128-node directed Fourier rule is required')
    bits = n.bit_length()-1
    result = [values[int(f'{j:0{bits}b}'[::-1], 2)] for j in range(n)]
    roots = _roots(mp.iv.dps)
    size = 2
    while size <= n:
        half = size//2
        for start in range(0, n, size):
            for j in range(half):
                root = roots[(n-j*(n//size)) % n]
                before = result[start+j]
                other = root*result[start+j+half]
                result[start+j], result[start+j+half] = before+other, before-other
        size *= 2
    return result


def generator_coefficients(generator, center, disk):
    """Cauchy/DFT coefficient intervals with a rigorous analytic alias tail."""
    values = [generator.values(mp.iv.mpc(center)+SAMPLE_RADIUS*z) for z in _roots(mp.iv.dps)]
    transforms = [_fft([row[j] for row in values]) for j in range(3)]
    M = mp.iv.mpf(disk['generator_norm_upper'])
    ratio = mp.iv.mpf(SAMPLE_RADIUS)/RADIUS
    answer = []
    for j in range(DEGREE+1):
        alias = M/mp.iv.mpf(RADIUS)**j*ratio**NODES/(1-ratio**NODES)
        error = _range(-_hi(alias), _hi(alias))
        answer.append(tuple(f[j]/(NODES*mp.iv.mpf(SAMPLE_RADIUS)**j)+mp.iv.mpc(error,error)
                            for f in transforms))
    return answer


def _apply(g, v):
    D, plus, minus = g
    return [-D*v[0]+minus*v[1], plus*v[0]+D*v[1]]


def _sum_vectors(values):
    result = [mp.iv.mpc(0), mp.iv.mpc(0)]
    for value in values:
        result = [a+b for a,b in zip(result,value)]
    return result


def _step_error(row):
    """Replay the full continuous polynomial residual integral, directed."""
    h = mp.iv.mpf(unpoint(row['step']))
    coefficients = row['residual_coefficient_norm_upper']
    low = sum((mp.iv.mpf(v)*h**(j+1)/(j+1) for j,v in enumerate(coefficients[:DEGREE])),mp.iv.mpf(0))
    high = sum((mp.iv.mpf(v)*h**(j+1)/(j+1) for j,v in enumerate(coefficients[DEGREE:],DEGREE)),mp.iv.mpf(0))
    ratio = h/RADIUS
    tail = mp.iv.mpf(row['geometry_disk']['generator_norm_upper'])*mp.iv.mpf(row['column_polynomial_norm_upper'])*h*ratio**(DEGREE+1)/((DEGREE+2)*(1-ratio))
    centered = mp.iv.mpf(row['centering_error_upper'])
    return low+high+tail+centered, {'coefficient_arithmetic_and_input_integral_upper':_upper_float(low),
                                 'polynomial_residual_tail_integral_upper':_upper_float(high),
                                 'analytic_generator_tail_integral_upper':_upper_float(tail)}


def taylor_step(generator, center, step, column):
    """One midpoint polynomial with continuous residual and endpoint enclosure."""
    disk = generator.disk(center)
    gc = generator_coefficients(generator, center, disk)
    coefficients = [column]
    for n in range(DEGREE):
        value = _sum_vectors(_apply(gc[j],coefficients[n-j]) for j in range(n+1))
        coefficients.append([_cmid(-mp.iv.j*v/(n+1)) for v in value])
    norms = []
    for n in range(2*DEGREE+1):
        value = _sum_vectors(_apply(gc[j],coefficients[n-j])
                             for j in range(max(0,n-DEGREE),min(DEGREE,n)+1))
        value = [mp.iv.j*v for v in value]
        if n < DEGREE:
            value = [v+(n+1)*_ivc(c) for v,c in zip(value,coefficients[n+1])]
        norms.append(_upper_float(_norm(value)))
    h = mp.iv.mpf(step)
    output = [mp.iv.mpc(0),mp.iv.mpc(0)]
    for coefficient in reversed(coefficients):
        output = [v*h+_ivc(c) for v,c in zip(output,coefficient)]
    next_column = [_cmid(v) for v in output]
    centering = _norm([v-_ivc(c) for v,c in zip(output,next_column)])
    poly_norm = sum((_norm([_ivc(v) for v in c])*h**j for j,c in enumerate(coefficients)),mp.iv.mpf(0))
    row = {'center':point(center),'step':point(step),'geometry_disk':disk,
           'residual_coefficient_norm_upper':norms,
           'column_polynomial_norm_upper':_upper_float(poly_norm),
           'centering_error_upper':_upper_float(centering),
           'column_end':[cpoint(v) for v in next_column]}
    error,components = _step_error(row)
    row.update(components);row['local_column_error_upper']=_upper_float(error)
    return next_column,row


def _initial_data(channel,config,sign,initializer_record):
    cert = initializer_record['certificate']
    if cert['group']!=32 or not cert['mathematical_endpoint_budget_pass']:
        raise ValueError('certified same-vacuum initializer required')
    points=[p for p in initializer_record['point_representations'] if p['sign']==sign and p['energy']==24]
    if len(points)!=1:raise ValueError('one certified E24 initializer representation required')
    p=points[0];w=mp.iv.mpc(*p['binary_ratio'])
    denominator=mp.iv.sqrt(1+w.real*w.real+w.imag*w.imag)
    exact=[1/denominator,w/denominator]
    column=[_cmid(v+mp.iv.mpc(0)) for v in exact]
    error=_norm([v-_ivc(c) for v,c in zip(exact,column)])
    mathematical=[r for r in cert['per_sign'] if r['sign']==sign][0]['projector_error_upper']
    projector=unpack(mathematical['binary_interval'])+unpack(p['projector_representation_error_upper']['binary_interval'])
    return column,_upper_float(error),_upper_float(projector)


def propagate_sign(channel,config,sign,initializer_record,*,checkpoint=None,checkpoint_key=None,progress=None):
    """Fixed-resolution real propagation; completed steps resume exactly."""
    with _precision(PRECISION):
        generator=ExactProfileGenerator(channel,config,sign)
        begin=mp.iv.ln(generator.g['delta0'])
        end=mp.iv.ln(generator.qh-3*mp.iv.pi/4)
        start,stop=_mid(begin),_mid(end)
        count=int(mp.ceil((stop-start)/MAX_STEP))
        grid=[start+(stop-start)*j/count for j in range(count+1)];grid[-1]=stop
        column,initial_error,initial_projector_error=_initial_data(channel,config,sign,initializer_record)
        start_margin=_upper_float(mp.iv.mpf(generator.disk(start)['generator_norm_upper'])*abs(begin-mp.iv.mpf(start)))
        endpoint_margin=_upper_float(mp.iv.mpf(generator.disk(stop)['generator_norm_upper'])*abs(end-mp.iv.mpf(stop)))
        binding={'key':checkpoint_key,'group':32,'sign':sign,'energy':24,'degree':DEGREE,
                 'count':count,'start':point(start),'stop':point(stop),'channel':channel,'config':config}
        rows=[]
        if checkpoint and Path(checkpoint).exists():
            saved=json.loads(Path(checkpoint).read_text())
            if saved['binding']!=binding:raise ValueError('propagation checkpoint belongs to changed inputs')
            rows=saved['steps']
            if rows:column=[uncpoint(v) for v in rows[-1]['column_end']]
            for j,row in enumerate(rows):
                if row['center']!=point(grid[j]) or row['step']!=point(grid[j+1]-grid[j]):
                    raise ValueError('checkpoint real-domain coverage changed')
                if _upper_float(_step_error(row)[0])!=row['local_column_error_upper']:
                    raise ValueError('checkpoint residual error replay changed')
        for j in range(len(rows),count):
            column,row=taylor_step(generator,grid[j],grid[j+1]-grid[j],column)
            rows.append(row)
            if (j+1)%16==0 or j+1==count:
                if checkpoint:
                    target=Path(checkpoint);temporary=target.with_suffix('.writing')
                    temporary.write_text(json.dumps({'binding':binding,'steps':rows},sort_keys=True,allow_nan=False))
                    temporary.replace(target)
                if progress:progress({'sign':sign,'steps':j+1,'total_steps':count,
                                      'accumulated_local_error':sum(r['local_column_error_upper'] for r in rows)})
        return {'binding':binding,'steps':rows,'initial_column_error_upper':initial_error,
                'initial_projector_error_upper':initial_projector_error,
                'start_coordinate_error_upper':start_margin,'target_coordinate_error_upper':endpoint_margin,
                'column_end':[cpoint(v) for v in column],
                'coordinate_policy':'exact-profile qh; delta target=qh-3pi/4 at rho1; physical chart interval, not assigned history duration',
                'Hamiltonian_label_policy':'frozen binary compact_mass and signed angular_eigenvalue; geometric horizon/slope exact-profile enclosures'}


def replay_propagation(payload):
    """Sum saved continuous residual envelopes; no generator or field run."""
    with _precision(PRECISION):
        binding=payload['binding'];rows=payload['steps']
        if binding['degree']!=DEGREE or binding['energy']!=24 or len(rows)!=binding['count']:
            raise ValueError('complete fixed pilot propagation required')
        error=mp.iv.mpf(payload['initial_column_error_upper'])+payload['start_coordinate_error_upper']+payload['target_coordinate_error_upper']
        at=unpoint(binding['start'])
        maxima={'coefficient_arithmetic_and_input_integral_upper':0.,'polynomial_residual_tail_integral_upper':0.,'analytic_generator_tail_integral_upper':0.}
        for row in rows:
            h=unpoint(row['step'])
            if unpoint(row['center'])!=at or not 0<h<=MAX_STEP or len(row['residual_coefficient_norm_upper'])!=2*DEGREE+1:
                raise ValueError('fixed-step continuous field domain has a gap or overlap')
            if min(row['geometry_disk']['B_real_lower'],row['geometry_disk']['sin_real_lower'])<=0:
                raise ValueError('geometry analytic disk margin failed')
            step,components=_step_error(row)
            if _upper_float(step)!=row['local_column_error_upper'] or any(row[k]!=v for k,v in components.items()):
                raise ValueError('saved continuous residual certificate changed')
            error+=mp.iv.mpf(row['local_column_error_upper']);at+=h
            for k,v in components.items():maxima[k]=max(maxima[k],v)
        if at!=unpoint(binding['stop']) or rows[-1]['column_end']!=payload['column_end']:
            raise ValueError('actual incoming endpoint changed')
        projector=mp.iv.mpf(payload['initial_projector_error_upper'])+(2+error)*error
        column=[_ivc(uncpoint(v)) for v in payload['column_end']]
        norm=sum((abs(v)**2 for v in column),mp.iv.mpf(0))
        raw_norm=_upper_float(abs(norm-1))
        if _hi(abs(norm-1))>_hi((2+error)*error):
            raise ArithmeticError('raw column norm defect exceeds its certified propagated error')
        channel=binding['channel'];m=mp.iv.mpf(channel['compact_mass']);ell=mp.iv.mpf(channel['angular_eigenvalue'])
        a=mp.iv.sqrt(3*mp.iv.pi/2-4);r=mp.iv.sqrt(2);H=mp.iv.sqrt(m*m+(ell/r)**2+(24/a)**2)
        exact_factor=channel['copy_count']*channel['degeneracy']/(8*mp.iv.pi**2*r*r*a)
        stored=mp.mpf(incoming_group_factor(channel,(1,-1)))
        factor=_range(min(_lo(exact_factor),stored),max(_hi(exact_factor),stored))
        lapse=4*mp.iv.pi*a*r*r*factor*2*H*projector
        return {'sign':binding['sign'],'energy':24,'steps':len(rows),
                'column_error_upper':_upper_float(error),'projector_error_upper':_upper_float(projector),
                'pointwise_weighted_lapse_kernel_error_upper':_upper_float(lapse),
                'raw_column_norm_defect_upper':raw_norm,'residual_component_step_maxima':maxima,
                'column_end':payload['column_end'],
                'units':'lapse action-gradient spectral kernel per positive-energy dE, inherited signed-family factor included',
                'spectral_integral_or_source_quadrature_claimed':False}
