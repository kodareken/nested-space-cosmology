"""Exact mixed-normal shift response of the existing spatial/reference action."""
from functools import lru_cache

import mpmath as mp
import numpy as np
import sympy as sp

from .nsc_incoming_cauchy_jets import incoming_cauchy_jets
from .nsc_incoming_lapse_coefficient import _hull, _interval
from .nsc_incoming_vacuum_tail_bound import _precision
from .nsc_magnetic_light_reference import magnetic_light_spectrum


@lru_cache(maxsize=1)
def weyl_trace_identity():
    """Owned first star-product sign; ordinary rank-one blocks fix Tr R."""
    I = sp.eye(2); sigma = (sp.Matrix([[0,1],[1,0]]), sp.Matrix([[0,-sp.I],[sp.I,0]]), sp.diag(1,-1))
    b = sp.Matrix(sp.symbols('b0:3', real=True)); z = sp.Matrix(sp.symbols('z0:3', real=True)); k = sp.Matrix(sp.symbols('k0:3', real=True))
    x = sp.Matrix(sp.symbols('x0:3', real=True)); t = sp.Symbol('t', real=True)
    pauli = lambda v: sum((v[i]*sigma[i] for i in range(3)), sp.zeros(2))
    Q = (I+pauli(b))/2; R = (t*I+pauli(x))/2
    star1 = sp.I/2*((pauli(z)/2)*(pauli(k)/2)-(pauli(k)/2)*(pauli(z)/2))
    product = sp.simplify(star1+pauli(z.cross(k))/4)
    blocks = sp.simplify(Q*R+R*Q-R-(b.dot(x)*I+t*pauli(b))/2)
    if product != sp.zeros(2) or blocks != sp.zeros(2): raise ArithmeticError('owned Weyl/Pauli sign identity failed')
    return {'first_star_order_residual':'0','rank_one_block_residual':'0',
            'trace_first_spatial_correction':'epsilon*b.dot(b_z.cross(b_k))/2',
            'conditions':['formal b.dot(b)=1','first order in the spatial gradient'],
            'higher_spatial_orders_at_linear_v':'r(T,z)=r_background(T)+v*T*z; each extra spatial derivative introduces another v'}


@lru_cache(maxsize=1)
def reference_shift_expressions():
    """Extract the linear-v current without momentum sampling or quadrature."""
    m,w,r = sp.symbols('m w r', positive=True)
    L = sp.symbols('L0:4', real=True); p = sp.symbols('p0:4', real=True)
    Ha,Hr,A2,R2,A3,R3 = sp.symbols('H_a H_r A2 R2 A3 R3', real=True)
    D = m*m+L[0]*L[0]+p[0]*p[0]
    def reduce(value):
        num,den = sp.fraction(sp.cancel(value))
        return sp.cancel(sp.rem(num,w*w-D,w)/den)
    def dt(value):
        return reduce(sum(sp.diff(value,L[j])*L[j+1]+sp.diff(value,p[j])*p[j+1] for j in range(3))
                      +sp.diff(value,w)*(L[0]*L[1]+p[0]*p[1])/w)
    h = sp.Matrix([-m,L[0],p[0]]); b = [-h/w]
    for n in range(1,4):
        derivative = b[-1].applyfunc(dt)
        normalization = sum((b[j].dot(b[n-j]) for j in range(1,n)),sp.S.Zero)
        b.append((-h.cross(derivative)/(2*w*w)-normalization*b[0]/2).applyfunc(reduce))
    normalization_residuals=[reduce(sum((b[i].dot(b[n-i]) for i in range(n+1)),sp.S.Zero)-(1 if n==0 else 0)) for n in range(4)]
    if normalization_residuals!=[0]*4: raise ArithmeticError('homogeneous formal Bloch normalization failed')
    replacement = {L[1]:-L[0]*Hr,L[2]:L[0]*(2*Hr**2-R2),L[3]:L[0]*(-6*Hr**3+6*Hr*R2-R3),
                   p[1]:-p[0]*Ha,p[2]:p[0]*(2*Ha**2-A2),p[3]:p[0]*(-6*Ha**3+6*Ha*A2-A3)}
    Lz = [0,-L[0]/r,4*L[0]*Hr/r,L[0]*(6*R2-18*Hr**2)/r]
    bz = [v.applyfunc(lambda f: reduce(sum(sp.diff(f,L[j])*Lz[j] for j in range(4))).subs(replacement, simultaneous=True).cancel()) for v in b]
    base = [v.applyfunc(lambda f: reduce(f.subs(replacement, simultaneous=True))) for v in b]
    bp = [v.applyfunc(lambda f: reduce(sp.diff(f,p[0])+sp.diff(f,w)*p[0]/w)) for v in base]
    kernels = {}
    for order in (1,2,3,4):
        value = sum((base[i].dot(bz[j].cross(bp[k])) for i in range(4) for j in range(4) for k in range(4)
                     if i+j+k==order-1),sp.S.Zero)
        kernels[order] = sp.factor(reduce(-p[0]*value/2))
    paired = {n:sp.factor((v+v.subs(L[0],-L[0]))/2) for n,v in kernels.items()}
    if paired[1] != 0 or paired[3] != 0: raise ArithmeticError('odd formal current orders did not vanish in the actual angular pair')
    if any(v.has(A3,R3) for v in kernels.values()): raise ArithmeticError('unexpected higher background derivative in linear-v response')
    return {'symbols':(m,L[0],p[0],w,r,Ha,Hr,A2,R2),'kernels':kernels,'paired_kernels':paired,
            'bloch':base,'bloch_spatial_response':bz,'normalization_residuals':[str(v) for v in normalization_residuals]}


@lru_cache(maxsize=1)
def local_shift_expressions():
    """Full2D beta Euler operator, then the existing mixed-normal plane.

    The displayed scalar contractions are the existing spherical_invariants
    contractions in the gauge N=1,a=a(T). Beta is varied BEFORE it is set0.
    Radius mixed jets are retained until after both Euler derivatives.
    """
    a = sp.symbols('a0:5', real=True); r = {(t,z):sp.Symbol(f'r{t}{z}',real=True)
        for t in range(5) for z in range(5-t)}
    B,Bt,Bz,Btz = sp.symbols('B B_T B_z B_Tz',real=True)
    betas = (B,Bt,Bz,Btz); zerob = dict.fromkeys(betas,0)
    A,CW,CE,Cg,Cq,hq,q = sp.symbols('A C_W C_E C_g C_q h_q q',real=True)
    aa,at,att = a[:3]; rr,rt,rz,rtt,rtz,rzz = [r[x] for x in ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))]
    Ha = at/aa
    baseR = -2*att/aa+2*Btz+4*Ha*Bz
    inverse = sp.Matrix([[1,-B],[-B,-1/aa**2]])
    hess = sp.Matrix([[rtt-(Bt+2*Ha*B)*rz, rtz-Ha*rz-aa*at*B*rt],
                      [rtz-Ha*rz-aa*at*B*rt, rzz-aa*at*rt+aa**2*Bz*rt+aa*at*B*rz]])
    gradient = sp.Matrix([rt,rz]); raised = inverse*gradient
    box = rtt+Ha*rt-rzz/aa**2-2*B*rtz-(Bt+Ha*B)*rz-Bz*rt
    grad2 = rt**2-rz**2/aa**2-2*B*rt*rz
    hess2 = sp.trace((inverse*hess)**2)
    hrr = (raised.T*hess*raised)[0]
    scalar = baseR-4*box/rr-2*(1+grad2)/rr**2
    weyl = (baseR+2*box/rr-2*(1+grad2)/rr**2)**2/3
    euler = 8*(box**2-hess2)/rr**2-4*baseR*(1+grad2)/rr**2
    barR2 = rr**2*baseR+2*rr*box-2*grad2
    s = sp.Symbol('s',real=True)
    weighted_euler = (8*s*s*rr*rr*(box*box-hess2)+16*s*s*(s-1)*rr*hrr
                      -4*(rr*rr*baseR-2*(s-1)*rr*box+2*(s-1)*grad2)*(1+s*s*grad2))
    # This polynomial is the contraction of the existing conformal path
    # g_s=r^(2s-2)g, r_s=r^s; integrate the formal parameter exactly.
    euler_average = sp.integrate(weighted_euler,(s,0,1))
    Hbar = rr*Ha-rt+B*rz-rr*Bz
    densities = {
        'compact/einstein_bulk':-4*sp.pi*A*aa*rr**2*scalar,
        'compact/weyl_bulk':-4*sp.pi*CW*aa*rr**2*weyl,
        'compact/euler_bulk_diagnostic':-4*sp.pi*CE*aa*rr**2*euler,
        'compact/maxwell_bulk':-2*sp.pi*Cg*q*q*aa/rr**2,
        'light/cylinder':-4*sp.pi*Cq*aa/rr**2,
        'light/bar_radial_R2_squared':hq*aa*barR2**2/(240*sp.pi*rr**2),
        'light/LLL_geometry_bar':q/(24*sp.pi)*(-aa*Hbar**2/rr**2+rz**2/(aa*rr**2)),
        'light/WZ_Weyl':aa*rr**2*sp.log(rr)*weyl/(80*sp.pi),
        'light/WZ_Euler':-11*aa*sp.log(rr)*euler_average/(1440*sp.pi*rr**2),
        'light/WZ_boxR':aa*((barR2-2)**2-rr**4*scalar**2)/(1440*sp.pi*rr**2),
        'light/WZ_gauge':-aa*q*q*sp.log(rr)/(12*sp.pi*rr**2),
    }
    def dt(f):
        return sum(sp.diff(f,a[i])*a[i+1] for i in range(4))+sum(sp.diff(f,x)*r[t+1,z] for (t,z),x in r.items() if t+z<4)
    def dz(f): return sum(sp.diff(f,x)*r[t,z+1] for (t,z),x in r.items() if t+z<4)
    v = sp.Symbol('v',real=True)
    plane = {x:(v if (t,z)==(1,1) else 0) for (t,z),x in r.items() if z}
    coefficients = {}; density_linear = {}; plane_residual = {}
    for name,density in densities.items():
        terms = [sp.factor(sp.diff(density,b).subs(zerob)) for b in betas]
        E = terms[0]-dt(terms[1])-dz(terms[2])+dt(dz(terms[3]))
        E = sp.factor(E.subs(plane,simultaneous=True))
        coefficient = sp.factor(sp.diff(E,v).subs(v,0))
        remainder = sp.simplify(E-coefficient*v)
        if remainder != 0: raise ArithmeticError('local beta response is not the established linear plane: '+name)
        coefficients[name] = coefficient; density_linear[name] = terms; plane_residual[name] = str(remainder)
    if coefficients['compact/euler_bulk_diagnostic'] != 0: raise ArithmeticError('constant Euler bulk identity failed')
    return {'axial_jets':a,'radius_jets':r,'parameters':(A,CW,CE,Cg,Cq,hq,q),
            'beta_jets':betas,'density_linear':density_linear,'coefficients':coefficients,
            'plane_residuals':plane_residual,'densities':densities}


@lru_cache(maxsize=1)
def integrated_reference_expressions():
    data = reference_shift_expressions(); m,L,p,w,r,Ha,Hr,A2,R2 = data['symbols']
    M = sp.Symbol('M',positive=True); a,d = sp.symbols('a d',positive=True)
    D = m*m+L*L+p*p; expressions = {}; moment_identities = {}; canonical = {}
    for order,power in ((2,5),(4,11)):
        raw = data['paired_kernels'][order]
        numerator = sp.factor(sp.cancel((raw*w**power).subs(w,sp.sqrt(D))))
        if numerator.has(w) or sp.simplify(raw.subs(w,sp.sqrt(D))-numerator/D**sp.Rational(power,2)) != 0:
            raise ArithmeticError('reference current gap-power reduction failed')
        poly = sp.Poly(sp.expand(numerator),p); integral = 0
        for (degree,),coefficient in poly.terms():
            if degree%2: raise ArithmeticError('even signed reference coefficient expected')
            if power<=degree+1: raise ArithmeticError('reference coefficient full-line moment is not integrable')
            moment = sp.simplify((M**(degree+1-power)*sp.beta(sp.Rational(degree+1,2),sp.Rational(power-degree-1,2))).rewrite(sp.gamma))
            moment_identities[f'p{degree}_gap{power}'] = str(moment)
            integral += coefficient*moment
        expressions[order] = sp.factor(a*d/(2*sp.pi)*integral.subs(M,sp.sqrt(m*m+L*L)))
        canonical[order] = numerator/D**sp.Rational(power,2)
    odd = sp.simplify(data['kernels'][3]+data['kernels'][3].subs(p,-p))
    if odd != 0: raise ArithmeticError('third-order current is not odd in momentum')
    M2=m*m+L*L; S=A2+R2+2*Hr**2-5*Ha*Hr+2*Ha**2
    compact={2:d*a*L*L/(12*sp.pi*r*M2),
             4:d*a*L*L/(120*sp.pi*r)*(S/M2**2+2*m*m*(Hr**2+4*R2)/M2**3-20*m**4*Hr**2/M2**4)}
    if any(sp.factor(expressions[n]-compact[n])!=0 for n in (2,4)):
        raise ArithmeticError('compact closed reference coefficient differs from full-line moments')
    return {'symbols':(m,L,r,a,d,Ha,Hr,A2,R2),'coefficients':expressions,
            'compact_coefficients':compact,'compact_reduction_residuals':['0','0'],
            'point_kernels':canonical,'moments':moment_identities,
            'odd_order_full_line_identity':str(odd)}


def closed_shift_coefficient(channels,ledger,*,precision=50):
    """Directed evaluation of the exact extracted coefficient, with input hulls."""
    if ledger['magnetic_flux']!=4 or [c['index'] for c in channels]!=list(range(33)):
        raise ValueError('unchanged q4 retained33 inventory required')
    domain=incoming_cauchy_jets(); normal=domain.normal_geometry()
    fields=domain.fields; scalar=lambda x:float(x.value[0,0,0].real)
    storedA2=scalar(fields[2].derivative(t=2))/normal['intrinsic_a']
    storedR2=scalar(fields[3].derivative(t=2))/normal['intrinsic_r']
    with _precision(precision):
        pi=_hull(mp.iv.pi,np.pi)
        a=_hull(mp.iv.sqrt(3*mp.iv.pi/2-4),normal['intrinsic_a']); r=_hull(mp.iv.sqrt(2),normal['intrinsic_r'])
        Ha=_hull((6-3*mp.iv.pi/2)/(2*a),normal['k_parallel']); Hr=_hull(-a/2,normal['k_perp'])
        A2=_hull((3*mp.iv.pi/2-3)/2,storedA2); R2=_hull((3*mp.iv.pi-10)/4,storedR2)
        hq=_hull(mp.iv.euler/2-mp.iv.mpf(25)/48,float(magnetic_light_spectrum(4).fourth_order_harmonic))
        CW=mp.iv.mpf(float(ledger['C_Weyl'])); A=mp.iv.mpf(float(ledger['A'])); q=mp.iv.mpf(ledger['magnetic_flux'])
        S=R2+A2+2*Hr**2-5*Ha*Hr+2*Ha**2
        totals={2:mp.iv.mpf(0),4:mp.iv.mpf(0)}; groups=[]
        for ch in channels[1:]:
            m=_hull(ch['compact_level']*mp.iv.pi/2,ch['compact_mass'])
            ell=_hull(mp.iv.sqrt(ch['angular_level']*(ch['angular_level']+4)),ch['angular_eigenvalue'])
            d=ch['copy_count']*ch['degeneracy']; L=ell/r; M2=m*m+L*L
            values={2:d*a*L*L/(12*pi*r*M2),
                    4:d*a*L*L/(120*pi*r)*(S/M2**2+2*m*m*(Hr**2+4*R2)/M2**3-20*m**4*Hr**2/M2**4)}
            for n,value in values.items():totals[n]+=value
            groups.append({'group':ch['index'],'reference_order2_interval':_interval(values[2]),'reference_order4_interval':_interval(values[4])})
        local={
            'compact/einstein_bulk':-16*pi*A*a*r,
            'compact/weyl_bulk':32*pi*CW*a*r*S/3,
            'light/bar_radial_R2_squared':-hq*a*r*S/(30*pi),
            'light/LLL_geometry_bar':a*q/(12*pi*r),
            'light/WZ_Weyl':-a*r*(S*mp.iv.ln(r)-R2+3*Hr**2+1/r**2-2*Ha*Hr+A2)/(30*pi),
            'light/WZ_Euler':11*a*r*(Hr**2+1/r**2)/(180*pi),
            'light/WZ_boxR':a*r*(R2-2*Hr**2+Ha*Hr+A2-2*Ha**2)/(60*pi),
        }
        ref=sum(totals.values(),mp.iv.mpf(0)); loc=sum(local.values(),mp.iv.mpf(0)); total=ref+loc
        return {'reference_order_intervals':{str(n):_interval(v) for n,v in totals.items()},
            'reference_total_interval':_interval(ref),'per_group_reference':groups,
            'local_channel_intervals':{n:_interval(v) for n,v in local.items()},'local_total_interval':_interval(loc),
            'total_coefficient_interval':_interval(total),'strictly_positive':total.a>0,
            'zero_local_channels':['compact/euler_bulk_diagnostic','compact/maxwell_bulk','light/cylinder','light/WZ_gauge'],
            'input_interpretation':{'geometry':'hulls contain exact owner expressions and stored incoming jet values',
                'ledger':'unchanged declared binary A/C_Weyl coefficients','harmonic':'exact EulerGamma/2-H4/4 plus stored value',
                'mass_angular':'existing definitions and stored labels; zero-angular reference coefficient is zero'},
            'precision':precision}


@lru_cache(maxsize=1)
def local_compact_identities():
    data=local_shift_expressions(); aa=data['axial_jets']; rr=data['radius_jets']
    A,CW,CE,Cg,Cq,hq,q=data['parameters']
    a,r,Ha,Hr,A2,R2=sp.symbols('a r H_a H_r A2 R2',real=True,nonzero=True)
    substitutions={aa[0]:a,aa[1]:a*Ha,aa[2]:a*A2,rr[0,0]:r,rr[1,0]:r*Hr,rr[2,0]:r*R2}
    S=A2+R2+2*Hr**2-5*Ha*Hr+2*Ha**2
    compact={
        'compact/einstein_bulk':-16*sp.pi*A*a*r,
        'compact/weyl_bulk':32*sp.pi*CW*a*r*S/3,
        'light/bar_radial_R2_squared':-hq*a*r*S/(30*sp.pi),
        'light/LLL_geometry_bar':a*q/(12*sp.pi*r),
        'light/WZ_Weyl':-a*r*(S*sp.log(r)-R2+3*Hr**2+1/r**2-2*Ha*Hr+A2)/(30*sp.pi),
        'light/WZ_Euler':11*a*r*(Hr**2+1/r**2)/(180*sp.pi),
        'light/WZ_boxR':a*r*(R2-2*Hr**2+Ha*Hr+A2-2*Ha**2)/(60*sp.pi),
    }
    residuals={name:str(sp.simplify(expr.subs(substitutions)-compact.get(name,0))) for name,expr in data['coefficients'].items()}
    if set(residuals.values())!={'0'}: raise ArithmeticError('compact local formula differs from full2D Euler coefficient')
    return {'closed_form_residuals':residuals,'coefficients':{k:str(v) for k,v in compact.items()},
            'plane_residuals':data['plane_residuals']}


def local_density_controls(ledger):
    """New algebra-only bridge to original densities; no Euler/probe generator."""
    from .nsc_incoming_local_constraints import local_action_densities
    from .nsc_light_restoration_action import LightRestorationAction
    from .nsc_spherical_local_history import LockedSphericalLocalAction
    spectrum=magnetic_light_spectrum(4); light=LightRestorationAction(spectrum); compact=LockedSphericalLocalAction(ledger)
    data=local_shift_expressions(); a=data['axial_jets']; r=data['radius_jets']
    raw=np.array([[1,0,0,0,0,0],[0,0,0,0,0,0],[.9,.2,0,.3,0,0],[1.4,-.4,.15,.18,-.11,.13]],float)
    values={a[0]:.9,a[1]:.2,a[2]:.3,r[0,0]:1.4,r[1,0]:-.4,r[0,1]:.15,r[2,0]:.18,r[1,1]:-.11,r[0,2]:.13}
    A,CW,CE,Cg,Cq,hq,q=data['parameters']
    values.update({A:ledger['A'],CW:ledger['C_Weyl'],CE:ledger['C_Euler'],Cg:ledger['C_gauge'],
                   Cq:float(spectrum.cylinder_density()),hq:float(spectrum.fourth_order_harmonic),q:4})
    errors={}; step=1e-20
    for j,slot in enumerate((0,1,2,4)):
        perturbed=raw.astype(complex);perturbed[1,slot]+=1j*step
        actual=local_action_densities(perturbed,light,compact)
        errors[str(slot)]={name:abs(float(value.imag/step)-float(data['density_linear'][name][j].subs(values))) for name,value in actual.items()}
    return {'beta_two_jet_derivative_errors':errors,'maximum':max(v for row in errors.values() for v in row.values()),
            'control_kind':'local algebraic density derivatives; no coordinate Euler contour or physical state/history',
            'complex_step':step}
