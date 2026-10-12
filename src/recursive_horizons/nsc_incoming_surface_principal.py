"""New principal entries of compatible normal functions; prior cu/cv reused."""
from functools import lru_cache

import numpy as np
import mpmath as mp
import sympy as sp

from .nsc_incoming_lapse_coefficient import _hull,_interval
from .nsc_incoming_vacuum_tail_bound import _precision,_lo,_hi
from .nsc_incoming_cauchy_jets import incoming_cauchy_jets
from .nsc_magnetic_light_reference import magnetic_light_spectrum
from .nsc_incoming_surface_reference_lapse import reference_lapse_group


@lru_cache(maxsize=1)
def local_principal_expressions():
    """Extract mixed density derivatives dictated by the FULL2D Euler sum.

    B=-L_(N_T,r_zz)-L_(N_z,r_Tz)+L_(N_zz,r_T).
    d=L_(beta_Tz,r_TT), e=L_(beta_Tz,r_zz).
    The beta_Tz dependence is2 times the base-curvature dependence.
    """
    a,r=sp.symbols('a r',positive=True)
    at,rt,rz,rtt,rtz,rzz,K=sp.symbols('a_T r_T r_z r_TT r_Tz r_zz K',real=True)
    nt,nz,nzz=sp.symbols('N_T N_z N_zz',real=True)
    A,CW,CE,Cg,Cq,hq,q=sp.symbols('A C_W C_E C_g C_q h_q q',real=True)
    R2=K+2*at/a*nt+2*nzz/a**2
    hess=sp.Matrix([[rtt-nt*rt-nz*rz/a**2,rtz-nz*rt-at/a*rz],
                    [rtz-nz*rt-at/a*rz,rzz-a*at*rt]])
    inverse=sp.diag(1,-1/a**2)
    box=rtt+(at/a-nt)*rt-rzz/a**2-nz*rz/a**2
    grad=rt**2-rz**2/a**2; raised=sp.Matrix([rt,-rz/a**2])
    hess2=sp.trace((inverse*hess)**2);hrr=(raised.T*hess*raised)[0]
    scalar=R2-4*box/r-2*(1+grad)/r**2
    weyl=(R2+2*box/r-2*(1+grad)/r**2)**2/3
    euler=8*(box**2-hess2)/r**2-4*R2*(1+grad)/r**2
    barR2=r*r*R2+2*r*box-2*grad
    s=sp.Symbol('s',real=True)
    weighted_euler=8*s*s*r*r*(box**2-hess2)+16*s*s*(s-1)*r*hrr-4*(r*r*R2-2*(s-1)*r*box+2*(s-1)*grad)*(1+s*s*grad)
    average=sp.integrate(weighted_euler,(s,0,1))
    densities={
        'compact/einstein_bulk':-4*sp.pi*A*a*r*r*scalar,
        'compact/weyl_bulk':-4*sp.pi*CW*a*r*r*weyl,
        'compact/euler_bulk_diagnostic':-4*sp.pi*CE*a*r*r*euler,
        'compact/maxwell_bulk':-2*sp.pi*Cg*q*q*a/r**2,
        'light/cylinder':-4*sp.pi*Cq*a/r**2,
        'light/bar_radial_R2_squared':hq*a*barR2**2/(240*sp.pi*r*r),
        'light/LLL_geometry_bar':q/(24*sp.pi)*(-a*(r*at/a-rt)**2/r**2-2*nz*rz/(a*r)+rz**2/(a*r*r)),
        'light/WZ_Weyl':a*r*r*sp.log(r)*weyl/(80*sp.pi),
        'light/WZ_Euler':-11*a*sp.log(r)*average/(1440*sp.pi*r*r),
        'light/WZ_boxR':a*((barR2-2)**2-r**4*scalar**2)/(1440*sp.pi*r*r),
        'light/WZ_gauge':-a*q*q*sp.log(r)/(12*sp.pi*r*r),
    }
    base={nt:0,nz:0,nzz:0,rz:0,rtz:0,rzz:0}; result={};ratio={}
    for name,L in densities.items():
        B=sp.factor((-sp.diff(L,nt,rzz)-sp.diff(L,nz,rtz)+sp.diff(L,nzz,rt)).subs(base))
        d=sp.factor((2*sp.diff(L,K,rtt)).subs(base));e=sp.factor((2*sp.diff(L,K,rzz)).subs(base))
        result[name]={'B':B,'d':d,'e':e};ratio[name]=str(sp.simplify(e+d/a**2))
    if set(ratio.values())!={'0'}:raise ArithmeticError('local highest beta derivative ratio failed')
    return {'symbols':(a,r,at,rt,rtt,K,A,CW,CE,Cg,Cq,hq,q),
            'entries':result,'beta_wave_ratio_residuals':ratio,'densities':densities,
            'raw_derivative_symbols':(nt,nz,nzz,rt,rtt,rtz,rzz,rz)}


@lru_cache(maxsize=1)
def reference_shift_principal():
    """All spatial orders from the exact linear Weyl-shift response.

    Formal exp(-i omega T+i q z) gives H star dP=H(k+q/2)dP
    and dP star H=dP H(k-q/2); no physical frequency is selected.
    """
    m,L,p,a,r=sp.symbols('m L p a r',positive=True)
    W=sp.sqrt(m*m+L*L+p*p);s=sp.Symbol('s',real=True)
    Ep=sp.sqrt(m*m+L*L+(p+s)**2);Em=sp.sqrt(m*m+L*L+(p-s)**2)
    # The scalar-current term surviving the full momentum/angular trace is
    # a*p*L^2/r * (1/Em-1/Ep)*omega / ((Ep+Em)^2-omega^2).
    first=sp.simplify(sp.diff(1/Em-1/Ep,s).subs(s,0))
    third=sp.simplify(sp.diff(1/Em-1/Ep,s,3).subs(s,0)/6)
    sum_second=sp.simplify(sp.diff(Ep+Em,s,2).subs(s,0)/2)
    d_kernel=sp.factor(-a*p*L*L/r*first/(16*W**4)/(2*a))
    e_kernel=sp.factor(-a*p*L*L/r*(third/(4*W**2)-first*sum_second/(4*W**3))/(2*a)**3)
    expected_d=-L*L*p*p/(16*r*W**7)
    expected_e=L*L*p*p*(5*(m*m+L*L)-2*p*p)/(32*r*a*a*W**9)
    residuals=[sp.simplify(d_kernel-expected_d),sp.simplify(e_kernel-expected_e)]
    if residuals!=[0,0]:raise ArithmeticError('third spatial Weyl/normal principal extraction failed')
    M=sp.Symbol('M',positive=True);D=sp.Symbol('D',positive=True)
    moment=lambda degree,power:sp.simplify((M**(degree+1-power)*sp.beta(sp.Rational(degree+1,2),sp.Rational(power-degree-1,2))).rewrite(sp.gamma))
    d_int=-D*a*L*L/(32*sp.pi*r)*moment(2,7)
    e_int=D*L*L/(64*sp.pi*r*a)*(5*M*M*moment(2,9)-2*moment(4,9))
    d_closed=-D*a*L*L/(120*sp.pi*r*M**4);e_closed=-d_closed/a**2
    if [sp.simplify(d_int-d_closed),sp.simplify(e_int-e_closed)]!=[0,0]:raise ArithmeticError('reference principal full-line moments failed')
    return {'symbols':(m,L,p,a,r,D),'kernels':{'d':d_kernel,'e':e_kernel},
            'closed':{'d':d_closed.subs(M,sp.sqrt(m*m+L*L)),'e':e_closed.subs(M,sp.sqrt(m*m+L*L))},
            'kernel_identity_residuals':['0','0'],'moment_identity_residuals':['0','0'],
            'beta_wave_ratio_residual':'0','frequency_labels_are_formal_only':True}


@lru_cache(maxsize=1)
def linear_weyl_response_identity():
    """Exact off-band linear transport/idempotence including all spatial shifts."""
    m,L,p,s,a,r,Eplus,Eminus,omega=sp.symbols('m L p s a r E_plus E_minus omega',real=True,nonzero=True)
    I=sp.eye(2);sx=sp.Matrix([[0,1],[1,0]]);sy=sp.Matrix([[0,-sp.I],[sp.I,0]]);sz=sp.diag(1,-1)
    Hplus=-m*sx+L*sy+(p+s)*sz;Hminus=-m*sx+L*sy+(p-s)*sz
    Pplus=(I-Hplus/Eplus)/2;Pminus=(I-Hminus/Eminus)/2
    Qplus=I-Pplus;Qminus=I-Pminus;V=-L*sy/r;S=Eplus+Eminus
    response=-Pplus*V*Qminus/(S+omega)-Qplus*V*Pminus/(S-omega)
    def reduce(value):
        num,den=sp.fraction(sp.cancel(value))
        num=sp.rem(num,Eplus**2-(m*m+L*L+(p+s)**2),Eplus)
        num=sp.rem(num,Eminus**2-(m*m+L*L+(p-s)**2),Eminus)
        return sp.factor(num/den)
    transport=(omega*response-Hplus*response+response*Hminus-(V*Pminus-Pplus*V)).applyfunc(reduce)
    projector=(Pplus*response+response*Pminus-response).applyfunc(reduce)
    current=sp.factor(-a*p*sp.trace(response))
    even=a*p*L*L/r*(1/Eminus-1/Eplus)*omega/(S*S-omega*omega)
    odd=2*sp.I*a*p*L*m*s*S/(r*Eplus*Eminus*(S*S-omega*omega))
    trace=sp.factor(current-even-odd)
    odd_pair=sp.factor(odd+odd.subs({p:-p,Eplus:Eminus,Eminus:Eplus},simultaneous=True))
    if transport!=sp.zeros(2) or projector!=sp.zeros(2) or trace!=0 or odd_pair!=0:
        raise ArithmeticError('full linear Weyl reference response failed')
    # The owned star_order coefficient gives f(p) star exp(i q z)
    # = f(p+q/(2a)) exp(i q z); check a degree4 test polynomial exactly.
    q=sp.Symbol('q_formal');f=sp.Function('f')
    polynomial=sum(sp.Symbol(f'c{i}')*p**i for i in range(5))
    star_shift=sum((sp.I/2)**n*(-1)**n*(sp.I*q)**n*sp.diff(polynomial,p,n)/(sp.factorial(n)*a**n) for n in range(5))
    shift=sp.expand(star_shift-polynomial.subs(p,p+q/(2*a)))
    if shift!=0:raise ArithmeticError('owned Weyl shift orientation failed')
    return {'linear_transport_residual':'0','linear_idempotency_residual':'0','scalar_trace_residual':'0',
            'odd_full_momentum_residual':'0','owned_Weyl_shift_residual':'0',
            'formal_current_even':str(even),'formal_current_odd':str(odd),
            'formal_mode':'exp(-i*omega*T+i*q_formal*z); coefficient labels only, no physical frequency/domain',
            'spectral_gap_conditions':'E_plus^2=m^2+L^2+(p+s)^2; E_minus^2=m^2+L^2+(p-s)^2; s=q_formal/(2a)'}


def local_principal_controls(ledger,step):
    """Exact-degree density stencils for NEW entries, never old cu/cv probes."""
    from .nsc_incoming_cauchy_jets import incoming_cauchy_jets
    from .nsc_incoming_local_constraints import scalar_taylor_coefficients,taylor_two_jets,local_action_densities
    from .nsc_light_restoration_action import LightRestorationAction
    from .nsc_magnetic_light_reference import magnetic_light_spectrum
    from .nsc_spherical_local_history import LockedSphericalLocalAction
    if step not in (.5,1.):raise ValueError('bounded exact polynomial stencil controls only')
    domain=incoming_cauchy_jets();raw=taylor_two_jets(scalar_taylor_coefficients(domain.fields),0.,0.)[0]
    light=LightRestorationAction(magnetic_light_spectrum(4));compact=LockedSphericalLocalAction(ledger)
    def mixed(field,slot,rslot):
        values={}
        for x in (-1,1):
            for y in (-1,1):
                data=raw.copy();data[field,slot]+=x*step;data[3,rslot]+=y*step
                for name,v in local_action_densities(data,light,compact).items():
                    values[name]=values.get(name,0)+x*y*v/(4*step*step)
        return values
    ntrzz=mixed(0,1,5);nzrtz=mixed(0,2,4);nzzrt=mixed(0,5,1)
    bttrtt=mixed(1,4,3);bttrzz=mixed(1,4,5)
    return {name:{'B':float((-ntrzz[name]-nzrtz[name]+nzzrt[name]).real),
                  'd':float(bttrtt[name].real),'e':float(bttrzz[name].real)} for name in ntrzz}


def directed_d_components(channels,ledger,*,precision=50):
    """New d/e entries only. No old cu or cv coefficient producer is called."""
    if ledger['magnetic_flux']!=4 or [c['index'] for c in channels]!=list(range(33)):
        raise ValueError('same q4 retained33 inventory required')
    normal=incoming_cauchy_jets().normal_geometry()
    with _precision(precision):
        pi=_hull(mp.iv.pi,np.pi);a=_hull(mp.iv.sqrt(3*mp.iv.pi/2-4),normal['intrinsic_a'])
        r=_hull(mp.iv.sqrt(2),normal['intrinsic_r']);Ha=_hull((6-3*mp.iv.pi/2)/(2*a),normal['k_parallel'])
        Hr=_hull(-a/2,normal['k_perp']);hq=_hull(mp.iv.euler/2-mp.iv.mpf(25)/48,float(magnetic_light_spectrum(4).fourth_order_harmonic))
        CW=mp.iv.mpf(float(ledger['C_Weyl']));ref=mp.iv.mpf(0);refB=mp.iv.mpf(0);groups=[]
        for ch in channels[1:]:
            m=_hull(ch['compact_level']*mp.iv.pi/2,ch['compact_mass'])
            ell=_hull(mp.iv.sqrt(ch['angular_level']*(ch['angular_level']+4)),ch['angular_eigenvalue'])
            L=ell/r;M2=m*m+L*L;D=ch['copy_count']*ch['degeneracy']
            value=-D*a*L*L/(120*pi*r*M2**2);ref+=value
            B=reference_lapse_group(m,L,a,r,Ha,Hr,D,pi);refB+=B
            groups.append({'group':ch['index'],'d_interval':_interval(value),'B_interval':_interval(B)})
        factors={'compact/weyl_bulk':-32*pi*CW/3,'light/bar_radial_R2_squared':hq/(30*pi),
                 'light/WZ_Weyl':mp.iv.ln(r)/(30*pi),'light/WZ_boxR':1/(60*pi)}
        local={};loc=mp.iv.mpf(0);locB=mp.iv.mpf(0)
        for name,K in factors.items():
            d=a*r*K;loc+=d
            B=r/a*K*((2*Ha+Hr) if name=='light/WZ_boxR' else (2*Ha-3*Hr))
            locB+=B
            local[name]={'B_interval':_interval(B),'d_interval':_interval(d),'e_interval':_interval(-d/a**2)}
        return {'reference_d_interval':_interval(ref),'local_d_interval':_interval(loc),
                'total_d_interval':_interval(ref+loc),'total_e_interval':_interval(-(ref+loc)/a**2),
                'reference_B_interval':_interval(refB),'local_B_interval':_interval(locB),
                'direct_total_B_interval':_interval(refB+locB),
                'local_entries':local,'per_group_reference_d':groups,
                'zero_local_channels':['compact/einstein_bulk','compact/euler_bulk_diagnostic','compact/maxwell_bulk',
                    'light/cylinder','light/LLL_geometry_bar','light/WZ_Euler','light/WZ_gauge'],
                'a_interval':_interval(a),'H_r_interval':_interval(Hr),'precision':precision}


@lru_cache(maxsize=1)
def local_euler_principal_identity():
    """Independent formal total-derivative extraction of the same new slots."""
    data=local_principal_expressions(); a,r,at,rt,rtt,K,*_=data['symbols']
    nt,nz,nzz,_,_,rtz,rzz,rz=data['raw_derivative_symbols']
    jets={(t,z):sp.Symbol(f'j{t}{z}',real=True) for t in range(5) for z in range(5-t)}
    aj=sp.symbols('aj0:5',real=True)
    substitution={a:aj[0],at:aj[1],K:-2*aj[2]/aj[0],r:jets[0,0],rt:jets[1,0],rz:jets[0,1],
                  rtt:jets[2,0],rtz:jets[1,1],rzz:jets[0,2]}
    bg={x:0 for (t,z),x in jets.items() if z}
    def dt(f):return sum(sp.diff(f,aj[j])*aj[j+1] for j in range(4))+sum(sp.diff(f,x)*jets[t+1,z] for (t,z),x in jets.items() if t+z<4)
    def dz(f):return sum(sp.diff(f,x)*jets[t,z+1] for (t,z),x in jets.items() if t+z<4)
    residuals={}
    for name,L in data['densities'].items():
        term=lambda x:sp.diff(L,x).subs({nt:0,nz:0,nzz:0}).subs(substitution,simultaneous=True)
        lapse=-dt(term(nt))-dz(term(nz))+dz(dz(term(nzz)))
        shift=dt(dz((2*sp.diff(L,K)).subs({nt:0,nz:0,nzz:0}).subs(substitution,simultaneous=True)))
        extracted={'B':sp.diff(lapse,jets[1,2]).subs(bg),
                   'd':sp.diff(shift,jets[3,1]).subs(bg),'e':sp.diff(shift,jets[1,3]).subs(bg)}
        residuals[name]={key:str(sp.simplify(value-data['entries'][name][key].subs(substitution,simultaneous=True).subs(bg)))
                         for key,value in extracted.items()}
    if any(v!='0' for row in residuals.values() for v in row.values()):raise ArithmeticError('full local Euler principal extraction failed')
    return residuals


@lru_cache(maxsize=1)
def combined_principal_identities():
    """Use the already-certified cu equations; do not derive them again."""
    a,r,Ha,Hr,m,L,D,CW,hq=sp.symbols('a r H_a H_r m L D C_W h_q',real=True,nonzero=True)
    M2=m*m+L*L
    # Imported equation from nsc-incoming-lapse-coefficient.json/identity.
    old_Aref=-D*a*L*L/(120*sp.pi*r)*(4*m*m*Hr/M2**3+(Hr-Ha)/M2**2)
    dref=-D*a*L*L/(120*sp.pi*r*M2**2)
    Bref=reference_lapse_group(m,L,a,r,Ha,Hr,D,sp.pi)
    ref=sp.factor(Bref+(2*old_Aref+Hr*dref)/a**2)
    old_local={
        'compact/weyl_bulk':32*sp.pi*CW*a*r*(Ha-Hr)/3,
        'light/bar_radial_R2_squared':-hq*a*r*(Ha-Hr)/(30*sp.pi),
        'light/WZ_Weyl':-sp.log(r)*a*r*(Ha-Hr)/(30*sp.pi),
        'light/WZ_boxR':-a*r*(Ha+Hr)/(60*sp.pi),
    }
    local=local_principal_expressions();aa,rr,at,rt,rtt,K,A_,CW_,CE_,Cg_,Cq_,hq_,q_=local['symbols']
    replacements={aa:a,rr:r,at:a*Ha,rt:r*Hr,CW_:CW,hq_:hq}
    residues={}
    for name,row in local['entries'].items():
        B=row['B'].subs(replacements,simultaneous=True);d=row['d'].subs(replacements,simultaneous=True)
        residues[name]=str(sp.simplify(B+(2*old_local.get(name,0)+Hr*d)/a**2))
    if ref!=0 or set(residues.values())!={'0'}:raise ArithmeticError('new B does not satisfy the reused-cu principal identity')
    A,d=sp.symbols('A d',real=True,nonzero=True);B=-(2*A+Hr*d)/a**2;e=-d/a**2
    determinant=sp.factor(A*e-B*d);eliminated=sp.factor(e-d*B/A)
    expected=d*(A+Hr*d)/a**2
    if sp.simplify(determinant-expected)!=0 or sp.simplify(eliminated-expected/A)!=0:
        raise ArithmeticError('principal determinant/elimination identity failed')
    return {'reference_reused_cu_relation':str(ref),'local_reused_cu_relations':residues,
            'B_relation':'B=-(2*A+H_r*d)/a^2','e_relation':'e=-d/a^2',
            'determinant':'d*(A+H_r*d)/a^2','eliminated_w3_coefficient':'d*(A+H_r*d)/(a^2*A)',
            'determinant_identity_residual':'0','elimination_identity_residual':'0',
            'A_is_reused_certified_c_u':True,'c_u_or_c_v_producer_called':False}


def directed_principal_matrix(channels,ledger,cu_interval,*,precision=50):
    """New entries and determinant, with A imported from its prior certificate."""
    components=directed_d_components(channels,ledger,precision=precision)
    with _precision(precision):
        A=mp.iv.mpf(cu_interval);a=mp.iv.mpf(components['a_interval']);Hr=mp.iv.mpf(components['H_r_interval'])
        d=mp.iv.mpf(components['total_d_interval']);B=-(2*A+Hr*d)/a**2;e=-d/a**2
        directB=mp.iv.mpf(components['direct_total_B_interval'])
        lo=max(_lo(B),_lo(directB));hi=min(_hi(B),_hi(directB))
        if lo>hi:raise ArithmeticError('independent direct B and reused-cu relation do not overlap')
        B=mp.iv.mpf([lo,hi]);det=d*(A+Hr*d)/a**2;direct_det=A*e-B*d
        if _lo(A)<=0:raise ArithmeticError('reused positive cu certificate lost its sign')
        schur=det/A
        return {'components':components,'matrix_intervals':[[_interval(A),_interval(B)],[_interval(d),_interval(e)]],
            'row_order':['E_N','E_beta'],'column_order':['U or U_z','w_zz or w_zzz'],
            'reused_cu_interval':list(cu_interval),'determinant_interval':_interval(det),
            'direct_matrix_determinant_interval':_interval(direct_det),
            'eliminated_w3_coefficient_interval':_interval(schur),
            'strictly_nonzero':_hi(det)<0 or _lo(det)>0,'strictly_negative':_hi(det)<0,
            'B_enclosure_method':'intersection of direct coefficient evaluation and exact relation to saved cu',
            'input_interpretation':'same geometric/harmonic definition and stored-value hulls; unchanged binary local ledger',
            'precision':precision}
