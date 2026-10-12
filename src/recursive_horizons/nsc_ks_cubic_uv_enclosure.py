"""Directed Volterra enclosure for the reduced massless cubic UV current.

The original profile, chart and source labels are retained. This encloses a
formal coefficient; it supplies neither the higher-energy remainder nor the
physical incoming gate. No refinement difference is used as a bound.
"""
from dataclasses import dataclass
from math import factorial

from flint import acb, arb, arb_series, ctx

from . import nsc_ks_ball_geometry as G
from .nsc_ks_ball_operator import slab_ball
from .nsc_ks_chebyshev_jet_bound import (
    chebyshev_endpoint_derivative_bounds, local_axial_interval_series,
)
from .nsc_local_incoming_family import LocalIncomingFamily


def _distance_integrand(z, analytic):
    if analytic and not z.real > 0:
        return acb('nan')
    a2 = 3*((1+z*z)*(acb.pi()/2-z.atan())-z)-1
    return 1/a2


def characteristic_distance(rho):
    """Enclose int_1^rho a^-2, including the whole supplied rho ball."""
    rho = G._rho_ball(rho)
    value = acb.integral(_distance_integrand,acb(1),acb(rho.mid()))
    if not value.is_finite() or not value.imag.contains(0):
        raise ArithmeticError('characteristic distance integral failed')
    return value.real+arb(25)/16*rho.rad()*arb(0,1)


@dataclass(frozen=True)
class CubicForcing:
    qz: arb
    qzz: arb
    current: arb
    coupling: arb


class CubicGeometry:
    """Exact-input ball jets, using the existing cutoff and Chebyshev owners."""
    def __init__(self,family,angular,*,bits=160):
        if not isinstance(family,LocalIncomingFamily):
            raise TypeError('owned local history required')
        if type(bits) is not int or bits < 80:
            raise ValueError('at least 80 precision bits required')
        self.family, self.bits = family,bits
        with ctx.workprec(bits):
            self.angular = arb(angular)
            if not self.angular.is_finite() or self.angular.contains(0):
                raise ValueError('finite nonzero angular label required')
        self.profiles = family.functions
        self.profile_bounds = tuple(chebyshev_endpoint_derivative_bounds(p,
            max(16,len(p.coefficients)-1)) for p in self.profiles)

    def forcing_series(self,rho,z_in,sign,order):
        """Taylor jets along the characteristic, uniform on the rho ball.

        First construct coordinate partials of b on a two-variable jet.
        Only then compose z(rho)=z_in-s*D(rho). Mixing these two operations
        would confuse b_rho with T_s b in the proved current equation.
        """
        if type(sign) is not int or sign not in (-1,1):
            raise ValueError('explicit characteristic sign +/-1 required')
        if type(order) is not int or not 0 <= order <= G.MAX_JET_ORDER-3:
            raise ValueError('forcing degree must leave three mixed derivative orders')
        degree=order+3
        with ctx.workprec(self.bits),G._JetWork(degree):
            rho=G._rho_ball(rho)
            bg=G.background_series(rho,degree)
            D=characteristic_distance(rho)
            z=arb(z_in)-sign*D
            shift=G.shift_series(rho,degree)
            window=G.plateau_series(-shift,self.family.normal_inner,
                                    self.family.normal_outer,degree)
            first,third=window*shift,window*shift**3/6
            w,U=(local_axial_interval_series(p,z,degree,bounds)
                 for p,bounds in zip(self.profiles,self.profile_bounds))
            pairs=[(i,n-i) for n in range(degree+1) for i in range(n+1)]
            radius={(i,j):(bg.r[i] if j==0 else arb(0))+first[i]*w[j]+third[i]*U[j]
                    for i,j in pairs}
            if not radius[(0,0)] > 0:
                raise G.SubdivisionNeeded('positive bivariate radius not enclosed')
            inverse={(0,0):1/radius[(0,0)]}
            for i,j in pairs[1:]:
                inverse[(i,j)]=-inverse[(0,0)]*sum(
                    (radius[(p,q)]*inverse[(i-p,j-q)]
                     for p in range(i+1) for q in range(j+1) if p or q),arb(0))
            coefficients={(i,j):self.angular/2*sum(
                (bg.a[p]*inverse[(i-p,j)] for p in range(i+1)),arb(0)) for i,j in pairs}
            # Both increments have identically zero constant, including at
            # interval basepoints. Their coefficients enclose total derivatives.
            x=arb_series([0,1],prec=order+1)
            dz=arb_series([0]+[-sign*bg.inv_a2[n-1]/n for n in range(1,order+1)],prec=order+1)
            xp=[arb_series([1],prec=order+1)];zp=[xp[0]]
            for _ in range(order):
                xp.append(xp[-1]*x);zp.append(zp[-1]*dz)
            def derivative(p,q):
                return sum((coefficients[(i+p,j+q)]
                    *factorial(i+p)/factorial(i)*factorial(j+q)/factorial(j)
                    *xp[i]*zp[j]
                    for n in range(order+1) for i in range(n+1) for j in (n-i,)),
                    arb_series([0],prec=order+1))
            b,bz,bp,bzz,bpz,bppz,bpzz=(derivative(p,q) for p,q in
                ((0,0),(0,1),(1,0),(0,2),(1,1),(2,1),(1,2)))
            a=arb_series(G._coeffs(bg.a,order),prec=order+1)
            ap=arb_series([(n+1)*bg.a[n+1] for n in range(order+1)],prec=order+1)
            current=(-a*a*b*bppz+a*a*bp*bpz-2*a*ap*b*bpz
                     -sign*b*bpzz+sign*bz*bpz-16*b**3*bz/a**2)
            result={'qz':-4*b*bz/a**2,'qzz':-4*(bz*bz+b*bzz)/a**2,
                    'current':current,'coupling':4*sign*b*b/a**2}
            if any(not v.is_finite() for series in result.values() for v in G._coeffs(series,order)):
                raise ArithmeticError('nonfinite characteristic forcing jet')
            return result

    def jets(self,rho,z):
        """Return b's coordinate partial derivatives on a rho/z rectangle.

        The rho and z variables are independent in these jets. The caller
        inserts the characteristic position AFTER choosing the rho rectangle.
        The bivariate reciprocal uses Taylor coefficients with factorials.
        """
        with ctx.workprec(self.bits),G._JetWork(3):
            rho,z = G._rho_ball(rho),G._binary_arb(z,'z')
            bg = G.background_series(rho,3)
            shift = G.shift_series(rho,3)
            window = G.plateau_series(-shift,self.family.normal_inner,
                                      self.family.normal_outer,3)
            first,third = shift*window,shift**3*window/6
            w,U = (local_axial_interval_series(p,z,3,bounds)
                   for p,bounds in zip(self.profiles,self.profile_bounds))
            radius = [first*w[j]+third*U[j] for j in range(4)]
            radius[0] += bg.r
            if not radius[0][0] > 0:
                raise G.SubdivisionNeeded('positive radius not enclosed')
            inverse = [radius[0].inv()]
            for j in range(1,4):
                inverse.append(-inverse[0]*sum((radius[k]*inverse[j-k]
                    for k in range(1,j+1)),arb_series([0],prec=4)))
            b = [bg.a*self.angular/2*v for v in inverse]
            derivative = lambda p,q:b[q][p]*factorial(p)*factorial(q)
            return {'a':bg.a[0],'a_rho':bg.a[1],
                    'b':derivative(0,0),'b_z':derivative(0,1),
                    'b_rho':derivative(1,0),'b_zz':derivative(0,2),
                    'b_rhoz':derivative(1,1),'b_rhorhoz':derivative(2,1),
                    'b_rhozz':derivative(1,2),
                    'b_reference_rho':(bg.a/bg.r*self.angular/2)[1]}

    def forcing(self,rho,z_in,sign):
        if type(sign) is not int or sign not in (-1,1):
            raise ValueError('explicit characteristic sign +/-1 required')
        with ctx.workprec(self.bits):
            z = arb(z_in)-sign*characteristic_distance(rho)
            j = self.jets(rho,z)
            a,ap,b,bz,bp,bzz,bpz,bppz,bpzz = (j[name] for name in
                ('a','a_rho','b','b_z','b_rho','b_zz','b_rhoz','b_rhorhoz','b_rhozz'))
            current = (-a*a*b*bppz+a*a*bp*bpz-2*a*ap*b*bpz
                       -sign*b*bpzz+sign*bz*bpz-16*b**3*bz/a**2)
            values = CubicForcing(-4*b*bz/a**2,
                -4*(bz*bz+b*bzz)/a**2,current,4*sign*b*b/a**2)
            if not all(v.is_finite() for v in vars(values).values()):
                raise ArithmeticError('nonfinite reduced forcing enclosure')
            return values


def volterra_cell(state,forcing,start,end):
    """Enclose a triangular affine system over one decreasing-rho cell.

    Every forcing entry is uniform on the cell. The intermediate qzz tube
    covers every partial cell integral. Thus current gets the nested integral
    as well as the direct forcing; its endpoint interval is not a sampled sum.
    """
    if not isinstance(forcing,CubicForcing) or len(state)!=3:
        raise TypeError('three state intervals and a uniform forcing required')
    start,end = arb(start),arb(end)
    if not start > end:
        raise ValueError('strictly decreasing cell required')
    qz,qzz,current = map(arb,state)
    h = end-start
    if not all(v.is_finite() for v in (qz,qzz,current,*vars(forcing).values())):
        raise ValueError('finite state and forcing required')
    qzz_tube = qzz+h*arb('0.5','0.5')*forcing.qzz
    return (qz+h*forcing.qz,qzz+h*forcing.qzz,
            current+h*(forcing.current+forcing.coupling*qzz_tube))


def taylor_volterra_cell(state,point,ranges,start,end,order):
    """Integral Taylor polynomial plus a directed Lagrange remainder.

    point/ranges contain derivatives divided by factorials at the starting
    point / throughout the cell. The qzz value needed in the current's highest
    derivative is enclosed on the entire cell using its uniform forcing.
    """
    start,end=arb(start),arb(end)
    if not start>end or type(order) is not int or order<1:
        raise ValueError('decreasing cell and positive Taylor order required')
    h=end-start
    qz,qzz,current=map(arb,state)
    qzz_tube=qzz+h*arb('0.5','0.5')*ranges['qzz'][0]
    initial_q=[qzz]+[point['qzz'][j-1]/j for j in range(1,order+1)]
    range_q=[qzz_tube]+[ranges['qzz'][j-1]/j for j in range(1,order+1)]
    def current_coefficient(data,qcoeff,j):
        return data['current'][j]+sum((data['coupling'][k]*qcoeff[j-k]
                                      for k in range(j+1)),arb(0))
    for j in range(order):
        weight=h**(j+1)/(j+1)
        qz+=point['qz'][j]*weight
        qzz+=point['qzz'][j]*weight
        current+=current_coefficient(point,initial_q,j)*weight
    weight=h**(order+1)/(order+1)
    return (qz+ranges['qz'][order]*weight,
            qzz+ranges['qzz'][order]*weight,
            current+current_coefficient(ranges,range_q,order)*weight)


def enclose_characteristic(model,z_in,sign,rho_up,*,cells=128,max_depth=16,order=0):
    """Finite full-path enclosure; unresolved walls fail rather than pass.

    A supplied z interval is covered continuously when all its boxes succeed.
    Wide axial boxes can need a separate z subdivision by the caller. This
    owner never labels point coverage as the complete incoming interval.
    """
    if not isinstance(model,CubicGeometry):
        raise TypeError('validated cubic geometry model required')
    if type(cells) is not int or cells < 1 or type(max_depth) is not int or max_depth<0:
        raise ValueError('positive cell count and nonnegative subdivision depth required')
    if type(order) is not int or not 0<=order<=G.MAX_JET_ORDER-3:
        raise ValueError('admissible forcing Taylor order required')
    with ctx.workprec(model.bits):
        upper = G._rho_ball(rho_up)
        if not upper > 1 or not upper.is_exact():
            raise ValueError('exact upstream point in the declared slab required')
        if not -G.shift_series(upper,0)[0] >= arb(model.family.normal_outer):
            raise ValueError('history must be flat in the unchanged upstream preparation')
        z_in = G._binary_arb(z_in,'incoming z')
        state = (arb(0),arb(0),arb(0))
        visited = 0
        def advance(state,start,end,depth):
            nonlocal visited
            try:
                if order:
                    forcing=model.forcing_series(slab_ball(end,start),z_in,sign,order)
                    point=model.forcing_series(start,z_in,sign,order)
                else:
                    forcing = model.forcing(slab_ball(end,start),z_in,sign)
            except G.SubdivisionNeeded:
                if depth >= max_depth:
                    raise G.SubdivisionNeeded('maximum cell subdivision reached; no enclosure emitted')
                mid = (start+end)/2
                state = advance(state,start,mid,depth+1)
                return advance(state,mid,end,depth+1)
            visited += 1
            if order:
                return taylor_volterra_cell(state,point,forcing,start,end,order)
            return volterra_cell(state,forcing,start,end)
        for j in range(cells):
            start = upper-(upper-1)*j/cells
            end = upper-(upper-1)*(j+1)/cells
            state = advance(state,start,end,0)
        g = model.jets(arb(1),z_in)
        a,b = g['a'],g['b']
        N = (sign*state[2]/a+a**3/2*(g['b_rho']**2-g['b_reference_rho']**2)
             -4*sign*b*b*state[0]/a+sign*a*b*g['b_rhoz'])
        return {'q_z':state[0],'q_zz':state[1],'J3':state[2],
                'N_bracket3':N,'accepted_cells':visited,
                'z_domain':z_in,'rho_up':upper,'sign':sign,
                'coefficient_scope':'history-minus-reference massless cubic; common homogeneous upstream differences zero',
                'physical_local_gate':'OPEN','higher_uv_remainder_bound':None}
