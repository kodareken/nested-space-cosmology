"""Portability of a saved finite nodal diagnostic, never a physical budget.

The arithmetic component propagates binary64 operation error through the
finite Chebyshev profiles, geometry, phase contraction and authenticated
60-array sum. Configured chart integration and root accuracy indicators are
reported separately. Their tolerances are algorithm settings, not proved
quadrature errors; combining them permits diagnostic comparison only.
The unchanged physical gate and its missing errors remain OPEN.
"""
from math import comb

import numpy as np
from scipy.special import expit

from .nsc_local_incoming_family import LocalIncomingFamily

_EPS = np.finfo(float).eps
_SUBNORMAL = np.nextafter(0., 1.)


def _real_array(value):
    raw = np.asarray(value)
    if raw.dtype.kind not in 'fiu' or not np.isfinite(raw).all():
        raise ValueError("finite real numeric diagnostic input required")
    return np.asarray(raw,float)


def _gamma(n):
    n = int(n)
    if n < 0 or n * _EPS >= 1:
        raise ValueError("finite operation count required")
    return n * _EPS / (1 - n * _EPS)


def _up(value):
    # Also charge the short positive arithmetic used to assemble a radius.
    return np.nextafter(np.asarray(value, float)*(1+_gamma(16))+_SUBNORMAL, np.inf)


class _Ball:
    """Forward absolute operation error, with no cancellation-based scale."""
    def __init__(self, value, error=0.):
        self.v, self.e = (x.copy() for x in np.broadcast_arrays(_real_array(value), _real_array(error)))
        if not np.isfinite(self.v).all() or not np.isfinite(self.e).all() or np.any(self.e < 0):
            raise ValueError("finite nonnegative diagnostic error required")

    def __add__(self, other):
        other = _as_ball(other)
        return _Ball(self.v + other.v, _up(self.e + other.e +
                     _EPS * (abs(self.v) + abs(other.v)) + _SUBNORMAL))

    __radd__ = __add__

    def __neg__(self):
        return _Ball(-self.v, self.e)

    def __sub__(self, other):
        return self + -_as_ball(other)

    def __rsub__(self, other):
        return _as_ball(other) + -self

    def __mul__(self, other):
        other = _as_ball(other)
        mass = (abs(self.v) + self.e) * (abs(other.v) + other.e)
        return _Ball(self.v * other.v, _up(abs(self.v)*other.e + abs(other.v)*self.e +
                     self.e*other.e + _EPS*mass + _SUBNORMAL))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = _as_ball(other)
        lower = np.nextafter(abs(other.v) - other.e, -np.inf)
        if np.any(lower <= 0):
            raise ValueError("diagnostic denominator error reaches zero")
        value = self.v / other.v
        return _Ball(value, _up((self.e + abs(value)*other.e)/lower +
                               _EPS*(abs(self.v)+self.e)/lower + _SUBNORMAL))

    def __rtruediv__(self, other):
        return _as_ball(other) / self

    def __pow__(self, exponent):
        if not isinstance(exponent, int) or exponent < 0:
            raise ValueError("nonnegative integer operation power required")
        value = _Ball(np.ones_like(self.v))
        for _ in range(exponent):
            value = value * self
        return value


def _as_ball(value):
    return value if isinstance(value, _Ball) else _Ball(value)


def _logistic(value):
    # expit is charged four ulps as an elementary-function operation. The
    # derivative bound is local to the input interval, retaining suppression
    # near the flat cutoff instead of applying a global 1/4 floor.
    lower, upper = expit(value.v-value.e), expit(value.v+value.e)
    derivative = np.maximum(lower*(1-lower), upper*(1-upper))
    derivative = np.where((value.v-value.e <= 0) & (value.v+value.e >= 0), .25, derivative)
    center = expit(value.v)
    return _Ball(center, _up(derivative*value.e + 4*_EPS*abs(center) + _SUBNORMAL))


def _window(x, inner, outer, order):
    distance = abs(x.v)
    result = [_Ball(np.zeros_like(distance)) for _ in range(order+1)]
    plateau = distance <= inner
    result[0].v[plateau] = 1.
    active = (distance > inner) & (distance < outer)
    # Branch uncertainty is confined to the exponentially flat endpoint.
    uncertain = ((abs(distance-inner) <= x.e) | (abs(distance-outer) <= x.e))
    if np.any(uncertain & (x.e > (outer-inner)/4096)):
        raise ValueError("coordinate indicator crosses a cutoff branch")
    result[0].v[abs(distance-inner)<=x.e] = 1.
    active &= ~uncertain
    # Within 1/4096 of an endpoint, exp(2-2048) times every inverse-u
    # polynomial through order three is below the least binary64 subnormal.
    # Flat continuation thus receives a subnormal derivative allowance.
    for j in range(order+1):
        result[j].e[uncertain] = _SUBNORMAL/(outer-inner)**j
    if not np.any(active):
        return result
    signed = np.sign(x.v[active])
    u = (_Ball(distance[active], x.e[active])-inner)/(outer-inner)
    if np.any((u.v-u.e <= 0) | (u.v+u.e >= 1)):
        # At an ulp-sized endpoint the flat exponential and every represented
        # derivative underflow; no nonzero profile is manufactured there.
        tiny = ((u.v-u.e <= 0) | (u.v+u.e >= 1))
        active_indices = np.flatnonzero(active)
        active[active_indices[tiny]] = False
        if not np.any(active):
            return result
        signed = np.sign(x.v[active])
        u = (_Ball(distance[active], x.e[active])-inner)/(outer-inner)
    q = 1/(1-u)-1/u
    eta, complement = _logistic(-q), _logistic(q)
    product = eta*complement
    q1 = 1/(1-u)**2 + 1/u**2
    q2 = 2*(1/(1-u)**3 - 1/u**3)
    q3 = 6*(1/(1-u)**4 + 1/u**4)
    f1 = -product
    f2 = product*(1-2*eta)
    f3 = -product*(1-6*eta+6*eta**2)
    jets = [eta, f1*q1, f2*q1**2+f1*q2,
            f3*q1**3+3*f2*q1*q2+f1*q3]
    for j in range(order+1):
        scaled = jets[j]*(signed/(outer-inner))**j
        result[j].v[active], result[j].e[active] = scaled.v, scaled.e
    return result


def _cheb_derivative(coefficients, order, scale):
    c = [_Ball(x) for x in coefficients]
    if order >= len(c):
        return [_Ball(0.)]
    for _ in range(order):
        n = len(c)-1
        c = [item*scale for item in c]
        derivative = [_Ball(0.) for _ in range(n)]
        for j in range(n, 2, -1):
            derivative[j-1] = (2*j)*c[j]
            c[j-2] = c[j-2]+j*c[j]/(j-2)
        if n > 1:
            derivative[1] = 4*c[2]
        derivative[0] = c[1]
        c = derivative
    return c


def _cheb_value(coefficients, argument):
    if len(coefficients) == 1:
        return coefficients[0]+_Ball(np.zeros_like(argument.v))
    c0, c1 = coefficients[-2], coefficients[-1]
    for item in reversed(coefficients[:-2]):
        c0, c1 = item-c1, c0+c1*(2*argument)
    return c0+c1*argument


def _profile(coefficients, family, z, order):
    left, right = _Ball(family.center)-family.axial_outer, _Ball(family.center)+family.axial_outer
    scale = 2/(right-left)
    argument = -(left+right)/(right-left)+scale*z
    window = _window(z-family.center, family.axial_inner, family.axial_outer, order)
    result = _Ball(np.zeros_like(z.v))
    for j in range(order+1):
        polynomial = _cheb_value(_cheb_derivative(coefficients, order-j, scale), argument)
        result = result+comb(order,j)*window[j]*polynomial
    return result


def coefficients_midpoint_roundoff(coeff_record, branch_record):
    """Charge the actual endpoint mass used in each recorded midpoint."""
    c, b = coeff_record['coefficients'], branch_record['branch']
    intervals = {**{k: b['intrinsic_intervals'][k] for k in ('a','r','Hr')},
                 **{k: b['imported_unchanged'][k] for k in ('A0','d','e')},
                 'A1': b['intervals']['A1_total'], 'C': c['C_interval']}
    def error(pair):
        x = _real_array(pair)
        if x.shape != (2,) or not np.isfinite(x).all():
            raise ValueError("two finite midpoint endpoints required")
        return float(_up(_gamma(2)*np.sum(abs(x))/2))
    errors = {k: error(v) for k,v in intervals.items()}
    errors['D'] = np.r_[0., [error(c['polynomial_intervals']['D'][str(j)]) for j in range(1,5)]]
    errors['F'] = np.array([error(c['polynomial_intervals']['F'][str(j)]) for j in range(3)])
    return errors


def _polyval(coefficients, argument):
    value = coefficients[-1]
    for item in reversed(coefficients[:-1]):
        value = item+value*argument
    return value


def _geometry(family, z, coefficients, errors):
    # Each original unit-basis profile is followed by a 2n-term projection.
    # The full-polynomial representation below evaluates the same exact
    # expression; charge the basis-projection operation mass independently.
    w = [_profile(family.coefficients[0],family,z,j) for j in range(4)]
    U = [_profile(family.coefficients[1],family,z,j) for j in range(2)]
    for values, coefficients_row in ((w,family.coefficients[0]), (U,family.coefficients[1])):
        for j, value in enumerate(values):
            mass = np.zeros_like(z.v)
            basis_error = np.zeros_like(z.v)
            for k, amplitude in enumerate(coefficients_row):
                unit = np.zeros(len(coefficients_row)); unit[k]=1.
                basis = _profile(unit,family,z,j)
                mass += abs(amplitude)*(abs(basis.v)+basis.e)
                basis_error += abs(amplitude)*basis.e
            value.e = _up(value.e+basis_error+_gamma(4*len(family.coefficients[0]))*mass)
    c = {k: _Ball(v,errors[k]) for k,v in coefficients.items() if k not in ('D','F')}
    D = [_Ball(v,e) for v,e in zip(coefficients['D'],errors['D'],strict=True)]
    F = [_Ball(v,e) for v,e in zip(coefficients['F'],errors['F'],strict=True)]
    A = c['A0']+c['A1']*w[0]
    B = -(2*A+(c['Hr']+w[0]/c['r'])*c['d'])/c['a']**2
    lapse = A*U[0]+B*w[2]+c['C']*w[1]**2+_polyval(D,w[0])
    shift = c['d']*U[1]+c['e']*w[3]+_polyval(F,w[0])*w[1]
    return np.column_stack((lapse.v,shift.v)), np.column_stack((lapse.e,shift.e))


def _phase(family, z, weights, axial, axial_error, rho_up, count, *, indicators):
    from .nsc_dirac_source_phase import _panels, _gauss_nodes, _normal_coordinate, _characteristic_distance
    from .nsc_ks_spacetime_variation import chart_coordinates
    from .nsc_pg_ks_metric_pullback import reference_chart
    from .nsc_ks_value_evaluator import one_direction_metric
    metric = one_direction_metric(family)
    panels = _panels(metric,rho_up)
    rho, quadrature_weights = _gauss_nodes(panels,count)
    # This finite recipe includes the Legendre eigensolve, Newton correction
    # and polynomial weights. A conservative n^2 operation allowance is
    # charged before applying the sensitive profile/reciprocal arithmetic.
    legendre = _gamma(64*count*count)
    reference_T = chart_coordinates(1.)[0]
    flux = [_Ball(np.zeros_like(z.v)), _Ball(np.zeros_like(z.v))]
    for point, weight in zip(rho,quadrature_weights,strict=True):
        a = reference_chart(float(point))[2]
        root_indicator = (2e-14+1e-12*abs(point)) if indicators else 0.
        rho_error = legendre*(rho_up-1)+root_indicator
        normal = _normal_coordinate(float(point))
        normal_error = _gamma(8)*(abs(chart_coordinates(float(point))[0])+abs(reference_T))
        distance = _characteristic_distance(float(point))
        distance_error = _gamma(8)*abs(distance)
        if indicators:
            # Existing integration accuracy indicators, kept apart from the
            # arithmetic component and never asserted as rigorous errors.
            normal_error += 2e-13*(max(1.,abs(reference_T))+max(1.,abs(reference_T+normal)))
            distance_error += 2e-13*max(1.,abs(distance))
        normal_error += rho_error/a
        distance_error += rho_error/(a*a)
        s = _Ball(np.full_like(z.v,normal),np.full_like(z.v,normal_error))
        window = _window(s,family.normal_inner,family.normal_outer,0)[0]
        r_ref = np.sqrt(1+point*point)
        rref = _Ball(r_ref,_gamma(4)*r_ref+rho_error*abs(point)/r_ref)
        # Panel-length variation and weight arithmetic are distinct from a
        # quadrature truncation error, which remains unconstructed.
        wgt = _Ball(weight,legendre*abs(weight)+
                    (2*root_indicator*abs(weight)/(rho_up-1) if indicators else 0.))
        for index, sign in enumerate((1,-1)):
            space = z-sign*_Ball(distance,distance_error)
            w0,w1 = (_profile(family.coefficients[0],family,space,j) for j in (0,1))
            u0,u1 = (_profile(family.coefficients[1],family,space,j) for j in (0,1))
            radius = rref+window*(s*w0+s**3*u0/6)
            radius_z = window*(s*w1+s**3*u1/6)
            flux[index] = flux[index]-wgt*radius_z/radius**3
    plus,minus = (_Ball(weight,_gamma(63)*abs(weight))*value
                  for weight,value in zip(weights,flux,strict=True))
    pi = _Ball(np.pi,_EPS*np.pi)
    lapse = (plus-minus)/_Ball(axial,axial_error)/(2*pi)
    shift = -(plus+minus)/(2*pi)
    return np.column_stack((lapse.v,shift.v)),np.column_stack((lapse.e,shift.e))


def saved_value_replay_envelope(family, nodes, coefficients, baseline, geometry, edge,
                               family_matter_arrays, signed_weights, rho_up, phase_nodes,
                               *, coefficient_errors):
    """Compare one unchanged finite recipe; no physical criterion is read."""
    if not isinstance(family,LocalIncomingFamily):
        raise TypeError("original LocalIncomingFamily required")
    nodes = _real_array(nodes)
    shape = (len(nodes),2)
    arrays = [_real_array(a) for a in family_matter_arrays]
    weights = _real_array(signed_weights)
    if (len(arrays)!=60 or any(a.shape!=shape or not np.isfinite(a).all() for a in arrays)
            or weights.shape!=(2,) or not np.isfinite(weights).all() or np.any(weights<0)
            or isinstance(phase_nodes,bool) or not isinstance(phase_nodes,(int,np.integer)) or phase_nodes<1):
        raise ValueError("60 authenticated arrays and original finite signed phase recipe required")
    base = np.broadcast_to(_real_array(baseline),shape)
    geometry,edge = _real_array(geometry),_real_array(edge)
    if any(a.shape!=shape or not np.isfinite(a).all() for a in (base,geometry,edge)):
        raise ValueError("finite original nodal parts required")
    z = _Ball(nodes)
    predicted_geometry, geometry_error = _geometry(family,z,coefficients,coefficient_errors)
    predicted_edge, edge_error = _phase(family,z,weights,coefficients['a'],coefficient_errors['a'],
                                      float(rho_up),int(phase_nodes),indicators=False)
    indicated_edge, indicated_error = _phase(family,z,weights,coefficients['a'],coefficient_errors['a'],
                                            float(rho_up),int(phase_nodes),indicators=True)
    if np.any(abs(geometry-predicted_geometry)>geometry_error):
        raise ValueError("current geometry changed beyond finite operation model")
    if np.any(abs(edge-predicted_edge)>edge_error):
        raise ValueError("current phase edge changed beyond finite operation model")
    # Both executions are charged; current parts must first lie within the
    # independently evaluated model. A mutation cannot enlarge its own
    # allowance by becoming a larger discrepancy from that model.
    matter_mass = np.sum([abs(a) for a in arrays],axis=0)
    summation = _gamma(60)*matter_mass
    assembly = _gamma(3)*(abs(base)+abs(geometry)+matter_mass+abs(edge))
    arithmetic = _up(2*(geometry_error+edge_error+summation+assembly))
    indicator = _up(2*np.maximum(0,indicated_error-edge_error)+abs(indicated_edge-predicted_edge))
    return {'arithmetic':arithmetic,'coordinate_accuracy_indicator':indicator,
            'summary_envelope':_up(arithmetic+indicator),
            'coordinate_settings':{'quad_epsabs':2e-13,'quad_epsrel':2e-13,
                                   'root_xtol':2e-14,'root_rtol':1e-12},
            'scope':'saved finite nodal diagnostic; configured indicators are not proved errors',
            'physical_error_bound':None,'phase_quadrature_error_bound':None,
            'physical_constraint_status':'OPEN'}


def require_saved_value_summaries(record, residual, envelope):
    """Only the two maxima and their joint sup receive diagnostic allowances."""
    residual,envelope = _real_array(residual),_real_array(envelope)
    if (residual.ndim!=2 or residual.shape[1]!=2 or residual.shape!=envelope.shape
            or not np.isfinite(residual).all() or not np.isfinite(envelope).all()
            or np.any(envelope<0)):
        raise ValueError("finite nodal residual and nonnegative envelope required")
    actual = np.max(abs(residual),axis=0)
    maxima = record['constraint_maxima']
    merit = record['merit']
    if (not isinstance(maxima,list) or len(maxima)!=2
            or any(type(x) is not float or not np.isfinite(x) or x<0 for x in maxima)
            or type(merit) is not float or not np.isfinite(merit) or merit<0
            or merit!=max(maxima)):
        raise ValueError("saved finite diagnostic summary structure changed")
    bounds = np.max(envelope,axis=0)
    if np.any(abs(actual-np.asarray(maxima))>bounds) or abs(float(max(actual))-merit)>max(bounds):
        raise ValueError("saved finite diagnostic changed beyond replay envelope")
    return {'constraint_maxima_bounds':bounds.tolist(),'merit_bound':float(max(bounds))}


def require_saved_value_receipt(actual, saved, envelope):
    """Only declared diagnostic scalars may differ in two complete receipts."""
    def identical(left,right):
        if type(left) is not type(right):
            return False
        if isinstance(left,dict):
            return left.keys()==right.keys() and all(identical(left[k],right[k]) for k in left)
        if isinstance(left,list):
            return len(left)==len(right) and all(identical(a,b) for a,b in zip(left,right))
        return left==right
    if not isinstance(actual,dict) or not isinstance(saved,dict) or actual.keys()!=saved.keys():
        raise ValueError("saved finite receipt fields changed")
    allowed = {'merit','constraint_maxima'}
    if any(not identical(actual[k],saved[k]) for k in actual.keys()-allowed):
        raise ValueError("non-diagnostic saved finite receipt changed")
    # A maxima row represents exactly the diagnostic summary being compared.
    residual = np.asarray([actual['constraint_maxima']],float)
    require_saved_value_summaries(actual,residual,np.zeros_like(residual))
    bound = np.max(np.asarray(envelope,float),axis=0,keepdims=True)
    return require_saved_value_summaries(saved,residual,bound)
