"""Conditional large-transfer bound for the linear source-fixed response.

Both real Fourier energies are integrated. This bounds only the linearized
raw matter insertion at the reference, not a finite-history error or a full
constraint. The full coherent covariance has two response terms.
"""
from fractions import Fraction


def _rational(value,name):
    if isinstance(value,bool):raise ValueError(name+' must be a finite rational bound')
    try:result=Fraction(value)
    except (ValueError,TypeError,OverflowError) as error:raise ValueError(name+' must be a finite rational bound') from error
    return result


def two_energy_tail_coefficients(order):
    """Exact upper coefficients after both E signs and both omega tails.

    Integrate |omega|^-p and (2|E|+|omega|)|omega|^-p over
    |E|>=E0, |omega|>=|E|/2. Density scales E0^(2-p), current E0^(3-p).
    """
    if isinstance(order,bool) or not isinstance(order,int) or order<=3:
        raise ValueError('integer derivative order >3 required for the weighted double tail')
    p=order
    density=Fraction(2**(p+1),(p-1)*(p-2))
    current=Fraction(2**p,p-3)*(Fraction(4,p-1)+Fraction(1,p-2))
    return {'density':density,'current':current,'density_power':2-p,'current_power':3-p}


def normal_insertion_bounds(normal_outer=Fraction(3,100),reference_radius_squared_lower=2):
    s=_rational(normal_outer,'normal outer radius')
    r2=_rational(reference_radius_squared_lower,'reference radius squared lower bound')
    if min(s,r2)<=0:raise ValueError('positive normal support and reference radius bounds required')
    return {'w':s*s/(2*r2),'U':s**4/(24*r2)}


def linear_large_transfer_bound(*,order,energy_cut,w_derivative_l1,U_derivative_l1,
                               mass_upper,angular_abs_upper,multiplicity,
                               axial_lower,radius_lower,normal_outer=Fraction(3,100),
                               reference_radius_squared_lower=2):
    """Exact-rational sufficient upper bound, conditional on supplied norms.

    w_derivative_l1/U_derivative_l1 bound the pth axial derivatives in L1(R).
    Inputs are bounds, not fitted action/source parameters. A conservative
    pi>3 gives 1/(2pi)^2 <1/36. Missing profile norms are never set to zero.
    """
    p=two_energy_tail_coefficients(order)
    values={k:_rational(v,k) for k,v in {
        'energy_cut':energy_cut,'w_derivative_l1':w_derivative_l1,
        'U_derivative_l1':U_derivative_l1,'mass_upper':mass_upper,
        'angular_abs_upper':angular_abs_upper,'multiplicity':multiplicity,
        'axial_lower':axial_lower,'radius_lower':radius_lower}.items()}
    E,W,V,m,ell,mu,a,r=(values[k] for k in (
        'energy_cut','w_derivative_l1','U_derivative_l1','mass_upper',
        'angular_abs_upper','multiplicity','axial_lower','radius_lower'))
    if min(E,mu,a,r)<=0 or min(W,V,m,ell)<0:raise ValueError('positive scales and nonnegative norm/mass bounds required')
    collar=normal_insertion_bounds(normal_outer,reference_radius_squared_lower)
    profile=collar['w']*W+collar['U']*V
    density=p['density']*E**p['density_power']
    current=p['current']*E**p['current_power']
    factor=mu*ell*profile/Fraction(36)
    return {
        'N_upper':factor*(4*(m+ell/r)*density+2*current/a),
        'beta_upper':factor*2*current,
        'collar_bounds':collar,'density_tail_factor':density,'current_tail_factor':current,
        'scope':'linear raw matter large-transfer region at the reference only',
        'region':'|E_i|>=E0, |E_o-E_i|>=|E_i|/2',
        'physical_local_gate':'OPEN','finite_history_remainder_bound':None,
        'nonresonant_remainder_bound':None,'profile_norms_certified_here':False,
        'two_source_energy_signs_included':True,'angular_multiplicity_applied_once':True,
        'coherence_retained':True,'source_law_changed':False}
