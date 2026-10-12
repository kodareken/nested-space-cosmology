"""Exact rational domain-of-dependence margin for the owned KS preparation.

This certifies a continuum geometric embedding. It does not make a finite
Fourier discretization causal or bound its numerical/source errors.
"""
from fractions import Fraction


def _q(value,name):
    if value is None or isinstance(value,bool):raise ValueError('explicit finite '+name+' required')
    try:return Fraction(value)
    except (TypeError,ValueError,OverflowError,ZeroDivisionError) as error:
        raise ValueError('finite '+name+' required') from error


def continuum_local_embedding(*,rho_up,rho_sigma,axial_lower,period_left,period_length,
                               physical_interval,axial_support,radius_positive):
    """Distance <= (rho_up-rho_sigma)/a_min^2 in each characteristic direction.

    The same initial covariance must be supplied on the upstream slice. The
    periodic and original coefficients agree in the central cell, and the
    compact support is strictly separated from its boundary by >2D.
    """
    up,sigma,a,left,length=(_q(v,n) for v,n in (
        (rho_up,'rho_up'),(rho_sigma,'rho_sigma'),(axial_lower,'axial lower'),
        (period_left,'period origin'),(period_length,'period length')))
    if sigma!=1 or not 1<up<=Fraction(33,32) or a<=0 or a>Fraction(4,5) or length<=0:
        raise ValueError('owned preparation slab and conservative a_min<=4/5 required')
    if radius_positive is not True:raise ValueError('radius positivity must be established for the history')
    if len(physical_interval)!=2 or len(axial_support)!=2:raise ValueError('two interval endpoints required')
    il,ir=(_q(v,'incoming interval') for v in physical_interval)
    sl,sr=(_q(v,'axial support') for v in axial_support)
    right=left+length
    if not left<sl<=il<ir<=sr<right:raise ValueError('positive incoming I inside compact support inside numerical cell required')
    distance=(up-sigma)/(a*a)
    padding=(sl-left,right-sr)
    cone=(il-distance,ir+distance)
    margins=(cone[0]-left,right-cone[1])
    if min(padding)<=2*distance:raise ValueError('owned causal buffer requires both paddings >2D')
    if min(margins)<=0:raise ValueError('incoming backward cone reaches another numerical cell')
    return {'schema':'NSC-KS-CONTINUUM-LOCAL-EMBEDDING-v1',
        'coordinate_interval':[str(il),str(ir)],'period_left':str(left),'period_length':str(length),
        'rho_up':str(up),'rho_sigma':str(sigma),'axial_lower':str(a),
        'characteristic_distance_upper':str(distance),'support_padding':[str(v) for v in padding],
        'incoming_backward_cone':[str(v) for v in cone],'cone_cell_margin':[str(v) for v in margins],
        'two_distance_buffer_margin':[str(v-2*distance) for v in padding],
        'continuum_local_embedding':True,'same_upstream_covariance_required':True,
        'physical_periodic_boundary_condition':False,'global_matching_required':False,
        'finite_Fourier_causality_certified':False,'numerical_field_error_bound':None,
        'source_accuracy_bound':None,'physical_local_gate':'OPEN','metric_timestep':False}
