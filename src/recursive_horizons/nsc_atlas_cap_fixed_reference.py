"""One constrained mass/gradient reference in the fixed physical cap space.

Reuses the radial forward/adjoint engine. Total fields are fixed physical
polynomials; the wall only clips their domain. No time trajectory is solved.
"""
from fractions import Fraction as Q
from collections.abc import Mapping
from types import MappingProxyType

import numpy as np
from numpy.polynomial.polynomial import polyval

from .nsc_atlas_background_cap import _budget
from .nsc_atlas_cap_fixed_basis import FixedCapBasis,exact_number
from .nsc_atlas_cap_fixed_bound import fixed_metric_bound
from .nsc_atlas_cap_geometric_reference import _positive_gap
from .nsc_atlas_cap_radial_reference import reference_potential
from .nsc_atlas_cap_variational_reference import _radial_variations,GRADIENT_ABSOLUTE_TOLERANCE
from .nsc_atlas_polar_cap import _scalar,_vector


def _delta(power,left,width):
    """Exact local powers of z*dP/dz-P, with z=left+width*xi."""
    return tuple((p-1)*power[p]+(p+1)*left/width*power[p+1] if p<5 else 4*power[5]
                 for p in range(6))


class _Fields:
    def __init__(self,basis,chi,density,budget):
        self.basis=basis;self.powers=[];self.deltas=[];self.shape_deltas=[];self.shapes=[]
        self.left=[];self.width=[]
        exact=(tuple(Q(float(v)) for v in chi),tuple(Q(float(v)) for v in density))
        for cell in range(basis.cells):
            left,right=basis.exact_knots[cell:cell+2];width=right-left
            shapes=basis.exact_shapes(cell);idx=basis.indices(cell)
            powers=tuple(tuple(sum(row[j]*coeff[int(i)] for j,i in enumerate(idx)) for row in shapes)
                         for coeff in exact)
            columns=tuple(zip(*shapes))
            self.powers.append(np.array(powers,dtype=float))
            self.deltas.append(np.array(_delta(powers[0],left,width),dtype=float))
            self.shape_deltas.append(np.array([_delta(column,left,width) for column in columns],dtype=float).T)
            self.shapes.append(np.array(shapes,dtype=float))
            self.left.append(float(left));self.width.append(float(width))
            budget.check(768)
        if not all(np.all(np.isfinite(v)) for arrays in
                   (self.powers,self.deltas,self.shape_deltas,self.shapes) for v in arrays):
            raise ValueError('fixed field polynomial conversion is nonfinite')

    def evaluate(self,cell,xi,budget):
        row=self.basis.values(cell,xi,budget=budget)
        budget.check(160)
        z=self.left[cell]+self.width[cell]*xi
        r,A0=float(np.tanh(z)),float(1/np.cosh(z)**2)
        powers=self.powers[cell];shapes=self.shapes[cell]
        chi,k=(float(polyval(xi,p)) for p in powers)
        chi_z=float(polyval(xi,np.arange(1,6)*powers[0,1:])/self.width[cell])
        pb=row['value']/A0
        if cell==0:
            q=shapes[1]+xi*xi*(shapes[3]+xi*xi*shapes[5])
            qfield=powers[0,1]+xi*xi*(powers[0,3]+xi*xi*powers[0,5])
            ratio=z/r if z else 1.
            fb=q*ratio/self.width[cell];f=qfield*ratio/self.width[cell]
            gap=_positive_gap(z) if z else 0.
            hb=((2*shapes[3]*xi*xi+4*shapes[5]*xi**4)/self.width[cell]+gap*row['value'])/A0
            h=((2*powers[0,3]*xi*xi+4*powers[0,5]*xi**4)/self.width[cell]+gap*chi)/A0
        else:
            gap=_positive_gap(z)
            fb=row['value']/r;f=chi/r
            hb=(polyval(xi,self.shape_deltas[cell])/z+gap*row['value'])/A0
            h=(float(polyval(xi,self.deltas[cell]))/z+gap*chi)/A0
        return dict(r=r,A0=A0,z=z,indices=row['indices'],value=row['value'],first=row['first'],
                    fb=fb,pb=pb,hb=hb,f=float(f),p=k/A0,h=float(h),chi=chi,k=k,chi_z=chi_z)


def fixed_cap_reference(basis,chi,canonical_density,Z,*,alpha,T_x=None,radius_rate=None,
                        incoming_w_plus=None,budget,rtol=2e-13,atol=2e-15,wall_coupler=None):
    """One instantaneous radial solve and weak coefficient RHS.

    The metric certificate covers exact rational shapes and binary inputs.
    The floating operator, forward/adjoint agreement and Gram residual are
    numerical consistency checks, not uniform truncation/error enclosures.
    A wall_coupler receives one immutable scalar metric/trace packet after
    the radial solve. It supplies T_x, radius_rate and incoming_w_plus from
    the same junction, instead of independently supplied diagnostic rates.
    """
    budget=_budget(budget);before=budget.evaluations
    if not isinstance(basis,FixedCapBasis): raise ValueError('require a fixed local basis')
    c=_vector(chi,'chi').copy();d=_vector(canonical_density,'canonical density').copy()
    if c.shape!=(basis.dimension,) or d.shape!=c.shape:
        raise ValueError('fixed coefficient dimensions differ')
    Z=exact_number(Z,'wall coordinate')
    alpha=_scalar(alpha,'alpha',nonnegative=True)
    if wall_coupler is None:
        tx=_scalar(T_x,'T_x',positive=True)
        rx=_scalar(radius_rate,'radius rate');incoming=_scalar(incoming_w_plus,'incoming wall field')
    elif (not callable(wall_coupler) or
          any(v is not None for v in (T_x,radius_rate,incoming_w_plus))):
        raise ValueError('wall coupler and prescribed rates are mutually exclusive')
    rtol=_scalar(rtol,'rtol',positive=True);atol=_scalar(atol,'atol',positive=True)
    if not 64*np.finfo(float).eps<=rtol<=1e-6 or not 1e-15<=atol<=1e-7:
        raise ValueError('fixed reference radial tolerances outside the numerical domain')
    if alpha>0 and np.sqrt(alpha/6)==0:
        raise ValueError('fixed reference alpha scale underflows')
    certificate=fixed_metric_bound(basis,c,d,Z,alpha=alpha,budget=budget)
    gram=basis.gram(Z,budget=budget)
    count=basis.dimension
    cut=float((Z-basis.exact_knots[-2])/(basis.exact_knots[-1]-basis.exact_knots[-2]))
    # One coordinate per element avoids subtracting almost equal physical
    # knots. The physical field space is unchanged; only quadrature is mapped.
    breaks=np.arange(basis.cells+1,dtype=float)
    breaks[-1]=basis.cells-1+cut
    fields=_Fields(basis,c,d,budget)

    def coefficients(cell,point,mass):
        if not 0<=cell<basis.cells or not breaks[cell]<=point<=breaks[cell+1]:
            raise ValueError('fixed radial evaluation escaped its physical cell')
        xi=point-cell
        if point==breaks[cell+1]: xi=cut if cell==basis.cells-1 else 1.
        row=fields.evaluate(cell,xi,budget)
        r,A0=row['r'],row['A0'];jacobian=fields.width[cell]*A0
        partials=np.zeros(2*count);budget.check(2*count+48)
        if r==0: return r,jacobian,1.,0.,0.,partials
        f,p,h=row['f'],row['p'],row['h']
        potential=reference_potential(np.array([f]),alpha,budget=budget)
        P,U=float(potential['Pdiff'][0]),float(potential['Pdiff_f'][0])
        A=A0-2*alpha*mass/r
        if not np.isfinite(A) or A<=0:
            raise ValueError('numerical fixed cap left the positive metric domain')
        square=p*p+h*h;a=alpha*square/(2*r)
        rho=r*r*P/2+A*square/4
        idx=row['indices']
        partials[idx]=r*r*U*row['fb']/2+A*h*row['hb']/2
        partials[count+idx]=A*p*row['pb']/2
        return r,jacobian,A,a,rho,partials

    try:
        with np.errstate(over='raise',invalid='raise',divide='raise'):
            centre_f=float(c[0]/fields.width[0]);centre_pi=float(d[0]/fields.width[0])
            q=np.sqrt(alpha/6)*centre_f
            centre_lapse=-float(np.log1p(2*np.sinh(q/2)**2))
            reference=_radial_variations(breaks,coefficients,centre_lapse=centre_lapse,
                derivative_count=2*count,nonzero_field=bool(np.any(c) or np.any(d)),
                budget=budget,rtol=rtol,atol=atol)
            mass,lapse=reference['mass'],reference['lapse']
            gc,gd=reference['gradient'][:count],reference['gradient'][count:]
            wall=fields.evaluate(basis.cells-1,cut,budget)
            R,A0=wall['r'],wall['A0']
            _,_,A,_,rho,_=coefficients(basis.cells-1,breaks[-1],mass)
            mu=float(np.exp(lapse));B=A*mu
            if wall_coupler is not None:
                packet=MappingProxyType(dict(radius=R,f=wall['f'],p=wall['p'],h=wall['h'],
                    A=A,B=B,delta=lapse,mass=mass,rho=rho,alpha=alpha))
                supplied=wall_coupler(packet)
                if not isinstance(supplied,Mapping) or set(supplied)!={'T_x','radius_rate','incoming_w_plus'}:
                    raise ValueError('wall coupler must return exactly the three junction inputs')
                tx=_scalar(supplied['T_x'],'coupled T_x',positive=True)
                rx=_scalar(supplied['radius_rate'],'coupled radius rate')
                incoming=_scalar(supplied['incoming_w_plus'],'coupled incoming wall field')
            if not abs(rx)<tx*B: raise ValueError('fixed cap wall speed is not timelike')
            V=2*mu*gram.solve(gd,budget=budget)
            idx=gram.wall['indices'];e=gram.wall['value']
            vw=float(e@V[idx]);drive=B*R*incoming
            force=-2*mu*gc.copy();force[idx]+=e*(drive-vw)
            c_rate=tx*V;d_rate=tx*gram.solve(force,budget=budget)
            # Propagate the independent gradient comparison through the
            # actual Gram solves; a small unscaled gradient gap can amplify.
            gradient_gap=reference['gradient']-reference['adjoint_gradient']
            velocity_gap=2*mu*gram.solve(gradient_gap[count:],budget=budget)
            wall_velocity_gap=float(e@velocity_gap[idx])
            force_gap=-2*mu*gradient_gap[:count].copy()
            force_gap[idx]-=e*wall_velocity_gap
            density_rate_gap=tx*gram.solve(force_gap,budget=budget)
            weak_flux=vw*(drive-vw)/(2*mu)
            mass_rate=float(gc@c_rate+gd@d_rate+rx*rho)
            expected=tx*weak_flux+rx*rho
            pweak=vw/B;hweak=R*incoming-pweak
            budget.check(80*count+128)
            result=dict(mass=mass,lapse=lapse,alpha=alpha,radius=R,wall_coordinate=float(Z),
                T_x=tx,radius_rate=rx,dimension=count,cells=basis.cells,
                gradient=dict(chi=gc.copy(),canonical_density=gd.copy()),
                adjoint_gradient=dict(chi=reference['adjoint_gradient'][:count],
                                      canonical_density=reference['adjoint_gradient'][count:]),
                gradient_discrepancy=reference['gradient_discrepancy'],
                gradient_absolute_tolerance=GRADIENT_ABSOLUTE_TOLERANCE,
                projected_gradient_discrepancy=dict(
                    chi_rate_infinity=float(np.max(abs(tx*velocity_gap))),
                    canonical_density_rate_infinity=float(np.max(abs(density_rate_gap))),
                    chi_rate_L2=float(np.sqrt(gram.squared_norm(tx*velocity_gap,budget=budget))),
                    canonical_density_rate_L2=float(np.sqrt(gram.squared_norm(density_rate_gap,budget=budget))),
                    wall_velocity=wall_velocity_gap,
                    scope='forward minus adjoint rates on the same numerical metric; not an error enclosure'),
                shape_derivative_Z=A0*rho,rates=dict(chi=c_rate,canonical_density=d_rate),
                mass_rate=mass_rate,physical_weak_flux=float(weak_flux),moving_density=float(rx*rho),
                weak_balance_defect=float(mass_rate-expected),incoming_flux_bound=float(tx*drive*drive/(8*mu)),
                wall=dict(A=A,B=B,f=wall['f'],p_trace=wall['p'],h_trace=wall['h'],
                    q_trace=wall['chi_z']/A0,p_weak=float(pweak),h_weak=float(hweak),velocity=vw,
                    velocity_discrepancy=float(vw-B*wall['p']),momentum_discrepancy=float(wall['p']-pweak),
                    normal_discrepancy=float(wall['h']-hweak),
                    incoming_discrepancy=float((wall['p']+wall['h'])/R-incoming),
                    polynomial_value_discrepancy=float(wall['chi']-e@c[idx]),
                    polynomial_density_discrepancy=float(wall['k']-e@d[idx]),
                    polynomial_slope_discrepancy=float(wall['chi_z']-gram.wall['first']@c[idx])),
                wall_coordinate_roundoff=float(Q(wall['z'])-Z),
                centre=dict(f=centre_f,Pi=centre_pi,Phi=0.,delta=centre_lapse),
                metric_certificate=certificate,scaled_gram_condition=gram.condition,
                sampled_minimum_A=reference['sampled_minimum_A'],field_sites=reference['field_sites'],
                requested_tolerances=dict(rtol=rtol,atol=atol),adjoint_order=reference['adjoint_order'],
                coefficient_evaluations=budget.evaluations-before,numerical_error_enclosure=False,
                domain='one fixed-physical-space constrained mass and weak RHS; no time evolution')
    except (FloatingPointError,np.linalg.LinAlgError) as error:
        raise ValueError('unresolved fixed cap arithmetic') from error
    def finite(value):
        if isinstance(value,dict): return all(finite(v) for v in value.values())
        return isinstance(value,str) or bool(np.all(np.isfinite(value)))
    if not finite(result): raise ValueError('nonfinite fixed cap reference')
    budget.check()
    return result
