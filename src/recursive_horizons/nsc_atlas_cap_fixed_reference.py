"""One constrained mass/gradient reference in the fixed physical cap space.

Reuses the radial forward/adjoint engine. An optional retained incoming
history carries the background while fixed polynomials describe departures.
The wall clips their domain. No time trajectory is solved.
"""
from fractions import Fraction as Q
from collections.abc import Mapping
from types import MappingProxyType

import numpy as np
from numpy.polynomial.polynomial import polyval
from numpy.polynomial.legendre import leggauss

from .nsc_atlas_background_cap import _budget
from .nsc_atlas_cap_fixed_basis import FixedCapBasis,exact_number
from .nsc_atlas_cap_fixed_bound import fixed_metric_bound,fixed_background_metric_bound
from .nsc_atlas_cap_background_jets import RetainedCapBackground,SMALL_Z
from .nsc_atlas_cap_panel_jets import history_piece_indices
from .nsc_atlas_cap_geometric_reference import _positive_gap
from .nsc_atlas_cap_potential_pair import reference_potential_pair
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


def _background_partition(basis,Z,T,background):
    """History polynomial breaks and formula switch in fixed coordinates.

    The onset constructor proves all preceding pieces identically zero.
    Their interior knots need not split numerical integration. Preserve the
    onset; a certified zero extension also removes the redundant quiet
    domain endpoint. The metric certificate retains every cut.
    """
    labels=background.inputs['labels'];onset=background.onset
    if onset is not None and onset['index']>=1:
        # The certified quiet extension and every pre-onset piece are the
        # same zero polynomial. Their domain endpoint is not a formula
        # break. Keeping it can create a one-ulp panel at the initial wall
        # that collapses when represented as cell_index + local_coordinate.
        # Actual history-domain checks remain in the field and metric owners.
        labels=(labels[onset['index']:] if background.profile.quiet_extension
                else np.r_[labels[:1],labels[onset['index']:]])
    knots=tuple(Q(float(v)) for v in labels);time=Q(T)
    events={Q(SMALL_Z),*(v-time for v in knots),*(time-v for v in knots)}
    breaks=[0.];owners=[]
    for cell in range(basis.cells):
        lo,top=basis.exact_knots[cell:cell+2];width=top-lo;hi=min(top,Z)
        ends=sorted(v for v in events if lo<v<hi)+[hi]
        for end in ends:
            value=cell+float((end-lo)/width)
            if value<=breaks[-1]:
                raise ValueError('background integration panel is unresolved in fixed coordinates')
            breaks.append(value);owners.append(cell)
    return np.array(breaks),np.array(owners,dtype=int)


def _background_time_projection(basis,fields,gram,breaks,owners,background,T,budget):
    # Each knot-split product has degree <=11 with the degree-seven onset
    # (<=9 for quintic history): six Gauss points suffice in
    # real arithmetic. Evaluation and Gram rounding remain numerical errors.
    nodes,weights=leggauss(6);loads=np.zeros((2,basis.dimension))
    for panel,cell in enumerate(owners):
        left,right=breaks[panel:panel+2];half=(right-left)/2
        for point,weight in zip(left+half*(nodes+1),weights):
            xi=float(point-cell);row=basis.values(int(cell),xi,budget=budget)
            z=fields.left[cell]+fields.width[cell]*xi
            jet=background.reduced_jets(T,z,budget=budget)
            loads[:,row['indices']]+=half*fields.width[cell]*weight*np.array([jet['k'],jet['k_T']])[:,None]*row['value']
            budget.check(8*len(row['indices']))
    return tuple(gram.solve(load,budget=budget) for load in loads)


def fixed_cap_reference(basis,chi,canonical_density,Z,*,alpha,T_x=None,radius_rate=None,
                        incoming_w_plus=None,budget,rtol=2e-13,atol=2e-15,wall_coupler=None,
                        background=None,background_time=None,radial_step_fraction=None,
                        background_panel_sides=False):
    """One instantaneous radial solve and weak coefficient RHS.

    The metric certificate covers exact rational shapes and binary inputs.
    The floating operator, forward/adjoint agreement and Gram residual are
    numerical consistency checks, not uniform truncation/error enclosures.
    A wall_coupler receives one immutable scalar metric/trace packet after
    the radial solve. It supplies T_x, radius_rate and incoming_w_plus from
    the same junction, instead of independently supplied diagnostic rates.
    With a retained background, coefficients are departures. Its explicit
    time gradient and the actual affine mass remainder remain in the result.
    """
    budget=_budget(budget);before=budget.evaluations
    if type(background_panel_sides) is not bool or (background_panel_sides and background is None):
        raise ValueError('background panel sides require a retained background and a boolean option')
    if not isinstance(basis,FixedCapBasis): raise ValueError('require a fixed local basis')
    c=_vector(chi,'chi').copy();d=_vector(canonical_density,'canonical density').copy()
    if c.shape!=(basis.dimension,) or d.shape!=c.shape:
        raise ValueError('fixed coefficient dimensions differ')
    Z=exact_number(Z,'wall coordinate')
    alpha=_scalar(alpha,'alpha',nonnegative=True)
    if background is not None:
        if not isinstance(background,RetainedCapBackground):
            raise ValueError('fixed cap requires original retained background data')
        time=_scalar(background_time,'background time')
    elif background_time is not None:
        raise ValueError('background time requires a retained background')
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
    certificate=(fixed_metric_bound(basis,c,d,Z,alpha=alpha,budget=budget) if background is None else
        fixed_background_metric_bound(basis,c,d,Z,background=background,T=time,alpha=alpha,budget=budget))
    gram=basis.gram(Z,budget=budget)
    count=basis.dimension
    cut=float((Z-basis.exact_knots[-2])/(basis.exact_knots[-1]-basis.exact_knots[-2]))
    # One coordinate per element avoids subtracting almost equal physical
    # knots. The physical field space is unchanged; only quadrature is mapped.
    breaks=np.arange(basis.cells+1,dtype=float)
    breaks[-1]=basis.cells-1+cut
    fields=_Fields(basis,c,d,budget)
    owners=np.arange(basis.cells)
    if background is not None:
        breaks,owners=_background_partition(basis,Z,time,background)
        a_h,b_h=_background_time_projection(basis,fields,gram,breaks,owners,background,time,budget)
    panel_indices=None
    if background_panel_sides:
        panel_indices=[]
        for panel,cell in enumerate(owners):
            midpoint=breaks[panel]+(breaks[panel+1]-breaks[panel])/2
            z=fields.left[cell]+fields.width[cell]*(midpoint-cell)
            panel_indices.append(history_piece_indices(background.profile,[time-z,time+z],budget=budget))
    nonzero_total_seen=False
    # A full-panel startup is exact for certified zero fields. For nonzero
    # data retain the existing solver heuristic: changing its initial step
    # amplified roundoff through the Gram solve in the affine control.
    exactly_quiet=not np.any(c) and not np.any(d)
    if exactly_quiet and background is not None:
        onset=background.onset
        rounding_margin=Q(16*np.finfo(float).eps*max(1.,abs(time),abs(float(Z))))
        exactly_quiet=bool(onset is not None and onset['index']>=1 and
            background.profile.quiet_extension and
            Q(time)+Z+rounding_margin<Q(float(background.inputs['labels'][onset['index']])))

    def total_fields(cell,xi,polynomial_indices=None):
        nonlocal nonzero_total_seen
        row=fields.evaluate(cell,xi,budget)
        if background is not None:
            bg=background.reduced_jets(time,row['z'],budget=budget,
                                       polynomial_indices=polynomial_indices)
            for key in ('f','p','h','chi','k','chi_z'): row[key]+=bg[key]
            row['background']=bg
        nonzero_total_seen |= any(row[key]!=0 for key in ('f','p','h'))
        return row

    def coefficients(panel,point,mass):
        if not 0<=panel<len(owners):
            raise ValueError('fixed radial evaluation escaped its physical cell')
        left,right=breaks[panel:panel+2]
        # A Runge--Kutta endpoint formed as t+c*h can round one ulp past
        # its declared panel. Restore that endpoint, without extrapolation.
        if not np.nextafter(left,-np.inf)<=point<=np.nextafter(right,np.inf):
            raise ValueError('fixed radial evaluation escaped its physical cell')
        point=min(right,max(left,point))
        cell=int(owners[panel])
        xi=point-cell
        if point==breaks[-1]: xi=cut
        row=total_fields(cell,xi,None if panel_indices is None else panel_indices[panel])
        r,A0=row['r'],row['A0'];jacobian=fields.width[cell]*A0
        partials=np.zeros(2*count+(background is not None));budget.check(len(partials)+48)
        if r==0: return r,jacobian,1.,0.,0.,partials
        f,p,h=row['f'],row['p'],row['h']
        P,U=reference_potential_pair(f,alpha,budget=budget)
        A=A0-2*alpha*mass/r
        if not np.isfinite(A) or A<=0:
            raise ValueError('numerical fixed cap left the positive metric domain')
        square=p*p+h*h;a=alpha*square/(2*r)
        rho=r*r*P/2+A*square/4
        idx=row['indices']
        partials[idx]=r*r*U*row['fb']/2+A*h*row['hb']/2
        partials[count+idx]=A*p*row['pb']/2
        if background is not None:
            bg=row['background']
            partials[-1]=r*r*U*bg['f_T']/2+A*(p*bg['p_T']+h*bg['h_T'])/2
        return r,jacobian,A,a,rho,partials

    try:
        with np.errstate(over='raise',invalid='raise',divide='raise'):
            centre_f=float(c[0]/fields.width[0]);centre_pi=float(d[0]/fields.width[0])
            if background is not None:
                axis=background.reduced_jets(time,0.,budget=budget)
                centre_f+=axis['f'];centre_pi+=axis['f_T']
            q=np.sqrt(alpha/6)*centre_f
            centre_lapse=-float(np.log1p(2*np.sinh(q/2)**2))
            reference=_radial_variations(breaks,coefficients,centre_lapse=centre_lapse,
                derivative_count=2*count+(background is not None),
                nonzero_field=bool(np.any(c) or np.any(d)) if background is None else False,
                budget=budget,rtol=rtol,atol=atol,
                first_step_fraction=1.0 if exactly_quiet else None,
                max_step_fraction=None if exactly_quiet else radial_step_fraction)
            mass,lapse=reference['mass'],reference['lapse']
            if mass==0 and nonzero_total_seen: raise ValueError('positive total background mass is unresolved')
            gc,gd=reference['gradient'][:count],reference['gradient'][count:2*count]
            explicit_time=float(reference['gradient'][-1]) if background is not None else 0.
            wall=total_fields(basis.cells-1,cut)
            R,A0=wall['r'],wall['A0']
            _,_,A,_,rho,_=coefficients(len(owners)-1,breaks[-1],mass)
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
            if background is not None: V-=a_h
            vw=float(e@V[idx])
            if background is not None: vw+=wall['background']['k']
            drive=B*R*incoming
            force=-2*mu*gc.copy();force[idx]+=e*(drive-vw)
            density_rate=gram.solve(force,budget=budget)
            if background is not None: density_rate-=b_h
            c_rate=tx*V;d_rate=tx*density_rate
            # Propagate the independent gradient comparison through the
            # actual Gram solves; a small unscaled gradient gap can amplify.
            gradient_gap=reference['gradient']-reference['adjoint_gradient']
            velocity_gap=2*mu*gram.solve(gradient_gap[count:2*count],budget=budget)
            wall_velocity_gap=float(e@velocity_gap[idx])
            force_gap=-2*mu*gradient_gap[:count].copy()
            force_gap[idx]-=e*wall_velocity_gap
            density_rate_gap=tx*gram.solve(force_gap,budget=budget)
            weak_flux=vw*(drive-vw)/(2*mu)
            mass_rate=float(tx*explicit_time+gc@c_rate+gd@d_rate+rx*rho)
            expected=tx*weak_flux+rx*rho
            pweak=vw/B;hweak=R*incoming-pweak
            budget.check(80*count+128)
            result=dict(mass=mass,lapse=lapse,alpha=alpha,radius=R,wall_coordinate=float(Z),
                T_x=tx,radius_rate=rx,dimension=count,cells=basis.cells,
                gradient=dict(chi=gc.copy(),canonical_density=gd.copy()),
                adjoint_gradient=dict(chi=reference['adjoint_gradient'][:count],
                                      canonical_density=reference['adjoint_gradient'][count:2*count]),
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
                    polynomial_value_discrepancy=float(wall['chi']-(wall['background']['chi'] if background is not None else 0.)-e@c[idx]),
                    polynomial_density_discrepancy=float(wall['k']-(wall['background']['k'] if background is not None else 0.)-e@d[idx]),
                    polynomial_slope_discrepancy=float(wall['chi_z']-(wall['background']['chi_z'] if background is not None else 0.)-gram.wall['first']@c[idx])),
                wall_coordinate_roundoff=float(Q(wall['z'])-Z),
                centre=dict(f=centre_f,Pi=centre_pi,Phi=0.,delta=centre_lapse),
                metric_certificate=certificate,scaled_gram_condition=gram.condition,
                sampled_minimum_A=reference['sampled_minimum_A'],field_sites=reference['field_sites'],
                requested_tolerances=dict(rtol=rtol,atol=atol),adjoint_order=reference['adjoint_order'],
                coefficient_evaluations=budget.evaluations-before,numerical_error_enclosure=False,
                domain='one fixed-physical-space constrained mass and weak RHS; no time evolution')
            if background is not None:
                wall_gap=float(wall['background']['k']-e@a_h[idx])
                remainder=float(explicit_time-gc@a_h-gd@b_h-wall_gap*(drive-vw)/(2*mu))
                result['gradient']['T']=explicit_time
                result['adjoint_gradient']['T']=float(reference['adjoint_gradient'][-1])
                result['background']=dict(time=time,radial_intervals=len(owners),
                    projected_time_chi=a_h,projected_time_density=b_h,wall_time_projection_gap=wall_gap)
                result['explicit_time_remainder']=remainder
                result['remainder_identity_defect']=float(result['weak_balance_defect']-tx*remainder)
                result['domain']='one retained-history plus fixed-departure constrained mass and weak RHS; no time evolution'
    except (FloatingPointError,np.linalg.LinAlgError) as error:
        raise ValueError('unresolved fixed cap arithmetic') from error
    def finite(value):
        if isinstance(value,dict): return all(finite(v) for v in value.values())
        return isinstance(value,str) or bool(np.all(np.isfinite(value)))
    if not finite(result): raise ValueError('nonfinite fixed cap reference')
    budget.check()
    return result
