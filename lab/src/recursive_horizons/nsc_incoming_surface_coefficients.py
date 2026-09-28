"""Fixed polynomial coefficients of the compatible included surface operators."""
from functools import lru_cache
import sympy as sp


@lru_cache(maxsize=1)
def homogeneous_reference_energy():
    """Only the new homogeneous energy polynomial; no old Cv/probe producer."""
    m,W=sp.symbols('m W',positive=True);L=sp.symbols('L0:4',real=True);p=sp.symbols('p0:4',real=True)
    Ha,Hr,A2,R2,A3,R3=sp.symbols('H_a H_r A2 R2 A3 R3',real=True)
    gap=m*m+L[0]*L[0]+p[0]*p[0]
    def reduce(f):
        num,den=sp.fraction(sp.cancel(f));return sp.cancel(sp.rem(num,W*W-gap,W)/den)
    def dt(f):return reduce(sum(sp.diff(f,L[j])*L[j+1]+sp.diff(f,p[j])*p[j+1] for j in range(3))
                            +sp.diff(f,W)*(L[0]*L[1]+p[0]*p[1])/W)
    h=sp.Matrix([-m,L[0],p[0]]);h1=sp.Matrix([0,L[1],p[1]])
    b0=-h/W;b1=h.cross(h1)/(2*W**3)
    b2=(-h.cross(b1.applyfunc(dt))/(2*W**2)-b1.dot(b1)*b0/2).applyfunc(reduce)
    E2=reduce(W*b1.dot(b1)/2);E3=reduce(W*b1.dot(b2))
    E4=reduce(h.cross(b1).dot(b2.applyfunc(dt))/(2*W)+W*b2.dot(b2)/2)
    replacement={L[1]:-L[0]*Hr,L[2]:L[0]*(2*Hr**2-R2),L[3]:L[0]*(-6*Hr**3+6*Hr*R2-R3),
        p[1]:-p[0]*Ha,p[2]:p[0]*(2*Ha**2-A2),p[3]:p[0]*(-6*Ha**3+6*Ha*A2-A3)}
    energies={n:sp.factor(expr.subs(replacement,simultaneous=True)) for n,expr in ((2,E2),(3,E3),(4,E4))}
    odd=reduce(energies[3]+energies[3].subs(p[0],-p[0]))
    if odd!=0:raise ArithmeticError('homogeneous third energy order did not cancel on the full momentum line')
    a,r,D,w=sp.symbols('a r D w',real=True,nonzero=True);M=sp.Symbol('M',positive=True)
    integrals={};moments={}
    for order,power in ((2,5),(4,13)):
        numerator=sp.factor(sp.cancel((energies[order]*W**power).subs(W,sp.sqrt(gap))))
        poly=sp.Poly(sp.expand(numerator),p[0]);integral=0
        for (degree,),coefficient in poly.terms():
            if degree%2 or power<=degree+1:raise ArithmeticError('nonintegrable/odd retained even-order energy moment')
            moment=sp.simplify((M**(degree+1-power)*sp.beta(sp.Rational(degree+1,2),sp.Rational(power-degree-1,2))).rewrite(sp.gamma))
            moments[f'p{degree}_gap{power}']=str(moment);integral+=coefficient*moment
        integrals[order]=sp.factor(D*a/(2*sp.pi)*integral.subs(M,sp.sqrt(m*m+L[0]*L[0])))
    change=sp.Poly(sp.expand(sum(integrals.values()).subs(Hr,Hr+w/r)-sum(integrals.values())),w)
    if change.degree()>4 or change.coeff_monomial(1)!=0:raise ArithmeticError('unexpected homogeneous energy polynomial degree/constant')
    # Reuse the already-certified Cu expression solely as a sign/measure check.
    old_cu=-D*a*L[0]**2/(120*sp.pi*r)*(4*m*m*Hr/(m*m+L[0]**2)**3+(Hr-Ha)/(m*m+L[0]**2)**2)
    cu_relation=sp.factor(sp.diff(integrals[4],R3)/r-old_cu)
    if cu_relation!=0:raise ArithmeticError('new energy polynomial disagrees with reused Cu equation')
    return {'symbols':(m,L[0],p[0],W,a,r,D,Ha,Hr,A2,R2,A3,R3,w),
        'point_energy_orders':energies,'integrated_energy_orders':integrals,
        'D_reference_coefficients':{j:sp.factor(change.coeff_monomial(w**j)) for j in range(1,5)},
        'standard_moments':moments,'odd_order_full_line_residual':str(odd),'reused_Cu_equation_residual':str(cu_relation)}


@lru_cache(maxsize=1)
def local_density_expressions():
    """The original covariant densities with full N,N_T,N_z,N_zz dependence."""
    N,a,r=sp.symbols('N a r',positive=True)
    nt,nz,nzz,at,att,rt,rz,rtt,rtz,rzz=sp.symbols('N_T N_z N_zz a_T a_TT r_T r_z r_TT r_Tz r_zz',real=True)
    A,CW,CE,Cg,Cq,hq,q=sp.symbols('A C_W C_E C_g C_q h_q q',real=True)
    baseR=-2*att/(N*N*a)+2*at*nt/(N**3*a)+2*nzz/(a*a*N)
    hess=sp.Matrix([[rtt-nt*rt/N-N*nz*rz/a**2,rtz-nz*rt/N-at*rz/a],
                    [rtz-nz*rt/N-at*rz/a,rzz-a*at*rt/N**2]])
    inverse=sp.diag(1/N**2,-1/a**2);raised=sp.Matrix([rt/N**2,-rz/a**2])
    box=sp.trace(inverse*hess);grad=rt*rt/N**2-rz*rz/a**2
    hess2=sp.trace((inverse*hess)**2);hrr=(raised.T*hess*raised)[0]
    scalar=baseR-4*box/r-2*(1+grad)/r**2
    weyl=(baseR+2*box/r-2*(1+grad)/r**2)**2/3
    euler=8*(box*box-hess2)/r**2-4*baseR*(1+grad)/r**2
    barR2=r*r*baseR+2*r*box-2*grad;s=sp.Symbol('s',real=True)
    weighted=8*s*s*r*r*(box*box-hess2)+16*s*s*(s-1)*r*hrr-4*(r*r*baseR-2*(s-1)*r*box+2*(s-1)*grad)*(1+s*s*grad)
    average=sp.integrate(weighted,(s,0,1))
    densities={
        'compact/einstein_bulk':-4*sp.pi*A*N*a*r*r*scalar,
        'compact/weyl_bulk':-4*sp.pi*CW*N*a*r*r*weyl,
        'compact/euler_bulk_diagnostic':-4*sp.pi*CE*N*a*r*r*euler,
        'compact/maxwell_bulk':-2*sp.pi*Cg*q*q*N*a/r**2,
        'light/cylinder':-4*sp.pi*Cq*N*a/r**2,
        'light/bar_radial_R2_squared':hq*N*a*barR2**2/(240*sp.pi*r*r),
        'light/LLL_geometry_bar':q/(24*sp.pi)*(-a*(r*at/a-rt)**2/(N*r*r)-2*nz*rz/(a*r)+N*rz*rz/(a*r*r)),
        'light/WZ_Weyl':N*a*r*r*sp.log(r)*weyl/(80*sp.pi),
        'light/WZ_Euler':-11*N*a*sp.log(r)*average/(1440*sp.pi*r*r),
        'light/WZ_boxR':N*a*((barR2-2)**2-r**4*scalar**2)/(1440*sp.pi*r*r),
        'light/WZ_gauge':-N*a*q*q*sp.log(r)/(12*sp.pi*r*r),
    }
    return {'symbols':(N,a,r,nt,nz,nzz,at,att,rt,rz,rtt,rtz,rzz,A,CW,CE,Cg,Cq,hq,q),'densities':densities}


@lru_cache(maxsize=1)
def local_surface_polynomials():
    data=local_density_expressions();N,a,r,nt,nz,nzz,at,att,rt,rz,rtt,rtz,rzz,A,CW,CE,Cg,Cq,hq,q=data['symbols']
    aj=sp.symbols('aj0:5',real=True);rj={(t,z):sp.Symbol(f'rj{t}{z}',real=True) for t in range(5) for z in range(5-t)}
    replace={a:aj[0],at:aj[1],att:aj[2],r:rj[0,0],rt:rj[1,0],rz:rj[0,1],rtt:rj[2,0],rtz:rj[1,1],rzz:rj[0,2]}
    def dt(f):return sum(sp.diff(f,aj[j])*aj[j+1] for j in range(4))+sum(sp.diff(f,x)*rj[t+1,z] for (t,z),x in rj.items() if t+z<4)
    def dz(f):return sum(sp.diff(f,x)*rj[t,z+1] for (t,z),x in rj.items() if t+z<4)
    gauge={N:1,nt:0,nz:0,nzz:0};spatial_zero={x:0 for (t,z),x in rj.items() if z}
    Ha,Hr,A2,R2,A3,R3,w,v=sp.symbols('H_a H_r A2 R2 A3 R3 w v',real=True)
    homogeneous={aj[0]:a,aj[1]:a*Ha,aj[2]:a*A2,aj[3]:a*A3,
                 rj[0,0]:r,rj[1,0]:r*Hr,rj[2,0]:r*R2,rj[3,0]:r*R3}
    zero_derivatives={**{x:0 for x in aj[1:]},**{x:0 for (t,z),x in rj.items() if t+z}}
    zero_derivatives[rj[1,1]]=v;zero_derivatives[aj[0]]=a;zero_derivatives[rj[0,0]]=r
    Dlocal={};Clocal={};homogeneous_E={};C_residuals={}
    for name,L in data['densities'].items():
        derivative=lambda variable:sp.diff(L,variable).subs(gauge).subs(replace,simultaneous=True)
        EN=derivative(N)-dt(derivative(nt))-dz(derivative(nz))+dz(dz(derivative(nzz)))
        EH=sp.factor(EN.subs(spatial_zero).subs(homogeneous,simultaneous=True));homogeneous_E[name]=EH
        difference=sp.Poly(sp.expand(EH.subs(Hr,Hr+w/r)-EH),w)
        if difference.degree()>4 or difference.coeff_monomial(1)!=0:raise ArithmeticError('local homogeneous source polynomial changed its degree/constant')
        Dlocal[name]={j:sp.factor(difference.coeff_monomial(w**j)) for j in range(1,5)}
        pure=sp.factor(EN.subs(zero_derivatives,simultaneous=True));C=sp.factor(sp.diff(pure,v,2)/2)
        if C.has(v):raise ArithmeticError('local v² coefficient is not grade4')
        Clocal[name]=C
    expected={'compact/weyl_bulk':64*sp.pi*CW/(3*a),
              'light/bar_radial_R2_squared':-hq/(15*sp.pi*a),'light/WZ_Weyl':-sp.log(r)/(15*sp.pi*a)}
    C_residuals={name:str(sp.factor(value-expected.get(name,0))) for name,value in Clocal.items()}
    if set(C_residuals.values())!={'0'}:raise ArithmeticError('full2D local C extraction mismatch')
    return {'symbols':(a,r,Ha,Hr,A2,R2,A3,R3,w,A,CW,CE,Cg,Cq,hq,q),
            'D_local_coefficients':Dlocal,'C_local_channels':Clocal,'C_local_residuals':C_residuals,
            'homogeneous_lapse_channels':homogeneous_E}


def polynomial_proof(saved_shift):
    """Reuse authenticated general Cv formulas; compute only their w dependence."""
    ref=homogeneous_reference_energy();local=local_surface_polynomials()
    names='m L0 p0 W a r D d H_a H_r A2 R2 A3 R3 w A C_W C_E C_g C_q h_q q'
    symbols={name:sp.Symbol(name) for name in names.split()}
    canon=lambda expr:expr.xreplace({x:symbols[str(x)] for x in expr.free_symbols})
    parse=lambda text:sp.sympify(text,locals=symbols)
    Hr,w,r=symbols['H_r'],symbols['w'],symbols['r']
    shiftref=sum(parse(x) for x in saved_shift['proof']['reference_coefficients'].values())
    shiftlocal={n:parse(x) for n,x in saved_shift['proof']['local']['coefficients'].items()}
    coeff=lambda expr,n:{str(j):str(sp.factor(sp.expand(expr.subs(Hr,Hr+w/r)-expr).coeff(w,j))) for j in range(1,n+1)}
    from .nsc_incoming_surface_quadratic_lapse import reference_quadratic_group
    m,L,a,D=[symbols[x] for x in ('m','L0','a','D')];M2=m*m+L*L
    Cref=reference_quadratic_group(m,L,a,r,D,sp.pi)
    Aref=-D*a*L*L/(120*sp.pi*r)*(4*m*m*Hr/M2**3+(Hr-symbols['H_a'])/M2**2)
    dref=-D*a*L*L/(120*sp.pi*r*M2**2)
    K=-32*sp.pi*symbols['C_W']/3+(symbols['h_q']+sp.log(r))/(30*sp.pi)
    Clocal=sum(canon(x) for x in local['C_local_channels'].values())
    identities={'reference_C_and_affine_A':str(sp.factor(Cref+(sp.diff(Aref,Hr)/r+dref/r)/a**2)),
                'local_C_and_affine_A':str(sp.factor(Clocal+(a*(K-1/(60*sp.pi))+a*(K+1/(60*sp.pi)))/a**2)),
                'homogeneous_reference_reused_Cu':ref['reused_Cu_equation_residual'],
                'odd_reference_energy_full_line':ref['odd_order_full_line_residual'],
                **{'full2D_C/'+k:v for k,v in local['C_local_residuals'].items()}}
    alpha0,alpha1,cc,dd=sp.symbols('alpha0 alpha1 C_total d_total')
    affineA=alpha0+alpha1*w;affineB=-(2*affineA+(Hr+w/r)*dd)/a**2
    delta=affineA*(-dd/a**2)-affineB*dd
    identities['Delta_w_plus_dC']=str(sp.factor((sp.diff(delta,w)+dd*cc).subs(alpha1,-a*a*cc-dd/r)))
    if set(identities.values())!={'0'}:raise ArithmeticError('surface coefficient exact identity failed')
    return {'D_reference':{str(k):str(v) for k,v in ref['D_reference_coefficients'].items()},
        'D_local':{n:{str(k):str(v) for k,v in row.items()} for n,row in local['D_local_coefficients'].items()},
        'F_reference':coeff(shiftref,2),'F_local':{n:coeff(v,2) for n,v in shiftlocal.items()},
        'C_reference':str(Cref),'C_local':{k:str(v) for k,v in local['C_local_channels'].items()},
        'standard_moments':ref['standard_moments'],'identities':identities,
        'point_energy_orders':{str(k):str(v) for k,v in ref['point_energy_orders'].items()},
        'homogeneous_local_lapse':{k:str(v) for k,v in local['homogeneous_lapse_channels'].items()},
        'definition':'D(w)=S_N+sum(D_j*w^j,j=1..4); F(w)=saved_cv+F1*w+F2*w²; R2,R3 fixed'}


def directed_coefficients(proof,channels,ledger,saved_cv,*,precision=50):
    """Directed evaluation of stored exact expressions; no symbolic derivation."""
    import mpmath as mp
    import numpy as np
    from .nsc_incoming_cauchy_jets import incoming_cauchy_jets
    from .nsc_incoming_lapse_coefficient import _hull,_interval
    from .nsc_incoming_vacuum_tail_bound import _precision
    from .nsc_magnetic_light_reference import magnetic_light_spectrum
    if ledger['magnetic_flux']!=4 or [x['index'] for x in channels]!=list(range(33)):
        raise ValueError('unchanged retained33 q4 inventory required')
    domain=incoming_cauchy_jets();normal=domain.normal_geometry()
    raw=lambda field,n:float(domain.fields[field].derivative(t=n).value[0,0,0].real)
    def evaluate(text,env):
        expr=sp.sympify(text,locals={x:sp.Symbol(x) for x in 'm L0 a r D d H_a H_r A2 R2 A3 R3 A C_W C_E C_g C_q h_q q'.split()})
        def walk(x):
            if x==sp.pi:return env['pi']
            if x.is_Symbol:return env[str(x)]
            if x.is_Rational:return mp.iv.mpf(int(x.p))/int(x.q)
            if x.is_Add:return sum((walk(y) for y in x.args),mp.iv.mpf(0))
            if x.is_Mul:
                v=mp.iv.mpf(1)
                for y in x.args:v*=walk(y)
                return v
            if x.is_Pow and x.exp.is_Integer:return walk(x.base)**int(x.exp)
            if x.func==sp.log:return mp.iv.ln(walk(x.args[0]))
            raise ValueError('unsupported interval expression: '+str(x))
        return walk(expr)
    with _precision(precision):
        iv=mp.iv;pi=_hull(iv.pi,np.pi);a=_hull(iv.sqrt(3*iv.pi/2-4),normal['intrinsic_a']);r=_hull(iv.sqrt(2),normal['intrinsic_r'])
        Ha=_hull((6-3*iv.pi/2)/(2*a),normal['k_parallel']);Hr=_hull(-a/2,normal['k_perp'])
        P1=6-3*iv.pi/2;P2=3-3*iv.pi/2
        A2=_hull((3*iv.pi/2-3)/2,raw(2,2)/raw(2,0));R2=_hull((3*iv.pi-10)/4,raw(3,2)/raw(3,0))
        A3=_hull((-P1*P2/4+3*a*a/2)/a,raw(2,3)/raw(2,0))
        R3=_hull(a*P2/4+3*a*P1/8+3*a**3/8,raw(3,3)/raw(3,0))
        hq=_hull(iv.euler/2-iv.mpf(25)/48,float(magnetic_light_spectrum(4).fourth_order_harmonic))
        env={'pi':pi,'a':a,'r':r,'H_a':Ha,'H_r':Hr,'A2':A2,'R2':R2,'A3':A3,'R3':R3,'h_q':hq,'q':iv.mpf(4),
             'A':iv.mpf(float(ledger['A'])),'C_W':iv.mpf(float(ledger['C_Weyl'])),'C_E':iv.mpf(float(ledger['C_Euler'])),
             'C_g':iv.mpf(float(ledger['C_gauge'])),'C_q':iv.mpf(float(magnetic_light_spectrum(4).cylinder_density()))}
        localC={n:evaluate(v,env) for n,v in proof['C_local'].items()};refC=iv.mpf(0)
        local={family:{str(j):sum((evaluate(row[str(j)],env) for row in proof[family+'_local'].values()),iv.mpf(0)) for j in range(1,degree+1)} for family,degree in [('D',4),('F',2)]}
        refs={family:{str(j):iv.mpf(0) for j in range(1,degree+1)} for family,degree in [('D',4),('F',2)]};groups=[]
        for ch in channels[1:]:
            m=_hull(ch['compact_level']*iv.pi/2,ch['compact_mass']);ell=_hull(iv.sqrt(ch['angular_level']*(ch['angular_level']+4)),ch['angular_eigenvalue'])
            groupenv={**env,'m':m,'L0':ell/r,'D':iv.mpf(ch['copy_count']*ch['degeneracy']),'d':iv.mpf(ch['copy_count']*ch['degeneracy'])}
            C=evaluate(proof['C_reference'],groupenv);refC+=C;row={'group':ch['index'],'C_interval':_interval(C)}
            for family,degree in [('D',4),('F',2)]:
                row[family]={}
                for j in range(1,degree+1):
                    value=evaluate(proof[family+'_reference'][str(j)],groupenv);refs[family][str(j)]+=value;row[family][str(j)]=_interval(value)
            groups.append(row)
        totals={family:{j:_interval(local[family][j]+refs[family][j]) for j in local[family]} for family in ('D','F')}
        totals['F']['0']=list(saved_cv);totals['D']['0']='S_N (exact formal finite source constant)'
        C_local=sum(localC.values(),iv.mpf(0))
        return {'polynomial_intervals':totals,'C_interval':_interval(C_local+refC),'C_local_interval':_interval(C_local),'C_reference_interval':_interval(refC),
                'local_C_channels':{n:_interval(v) for n,v in localC.items()},
                'local_polynomials':{f:{j:_interval(v) for j,v in row.items()} for f,row in local.items()},
                'reference_polynomials':{f:{j:_interval(v) for j,v in row.items()} for f,row in refs.items()},
                'reference_groups':groups,'geometry_hulls':{k:_interval(env[k]) for k in ('a','r','H_a','H_r','A2','R2','A3','R3')},
                'precision':precision,'input_interpretation':'same geometric definition/stored-value hulls and binary local ledger; angular signs included once in D'}


def new_algebra_controls(ledger):
    """One w-direction point control, plus raw local-density derivatives.

    No momentum integration or old constraint probe is performed.
    """
    import numpy as np
    from .nsc_incoming_cauchy_jets import incoming_cauchy_jets,IncomingNormalJetChange
    from .nsc_spatial_reference_symbol import reference_projector,SIGMA
    from .nsc_incoming_local_constraints import local_action_densities
    from .nsc_light_restoration_action import LightRestorationAction
    from .nsc_spherical_local_history import LockedSphericalLocalAction
    from .nsc_magnetic_light_reference import magnetic_light_spectrum
    data=homogeneous_reference_energy();m,L,p,W,a,r,D,Ha,Hr,A2,R2,A3,R3,w=data['symbols']
    domain=incoming_cauchy_jets((IncomingNormalJetChange('r',1,0,.02),))
    value=lambda field,n=0:float(domain.fields[field].derivative(t=n).value[0,0,0].real)
    aa,rr=value(2),value(3)
    ks=np.array([-1.3,-.4,.4,1.3]);ms=np.full(4,np.pi/2);ells=np.array([-np.sqrt(5),np.sqrt(5),-np.sqrt(5),np.sqrt(5)])
    exact=reference_projector(*domain.fields,ks,ms,ells)
    H=sum(x[:,None,None]*s for x,s in zip((-ms,ells/rr,ks/aa),SIGMA))
    common={a:aa,r:rr,Ha:value(2,1)/aa,Hr:value(3,1)/rr,A2:value(2,2)/aa,R2:value(3,2)/rr,A3:value(2,3)/aa,R3:value(3,3)/rr}
    errors={}
    for n,expr in data['point_energy_orders'].items():
        actual=np.trace(H@exact['orders'][n],axis1=-2,axis2=-1).real
        expected=np.array([float(expr.subs({**common,m:mm,L:ll/rr,p:kk/aa,W:np.sqrt(mm*mm+(ll/rr)**2+(kk/aa)**2)})) for kk,mm,ll in zip(ks,ms,ells)])
        errors[str(n)]=float(np.max(abs(actual-expected)))
    spectrum=magnetic_light_spectrum(4);light=LightRestorationAction(spectrum);compact=LockedSphericalLocalAction(ledger)
    local=local_density_expressions();syms=local['symbols'];N,aaS,rrS,nt,nz,nzz,at,att,rt,rz,rtt,rtz,rzz,A,CW,CE,Cg,Cq,hq,q=syms
    raw=np.array([[1.,0,0,0,0,0],[0,0,0,0,0,0],[aa,value(2,1),0,value(2,2),0,0],[rr,value(3,1),.15,value(3,2),-.11,.13]])
    env={N:1,aaS:aa,rrS:rr,nt:0,nz:0,nzz:0,at:raw[2,1],att:raw[2,3],rt:raw[3,1],rz:.15,rtt:raw[3,3],rtz:-.11,rzz:.13,
         A:ledger['A'],CW:ledger['C_Weyl'],CE:ledger['C_Euler'],Cg:ledger['C_gauge'],Cq:float(spectrum.cylinder_density()),hq:float(spectrum.fourth_order_harmonic),q:4}
    localerrors={};step=1e-20
    for symbol,slot in ((N,0),(nt,1),(nz,2),(nzz,5)):
        perturbed=raw.astype(complex);perturbed[0,slot]+=1j*step
        actual=local_action_densities(perturbed,light,compact)
        localerrors[str(slot)]={name:abs(float(value.imag/step)-float(sp.diff(local['densities'][name],symbol).subs(env))) for name,value in actual.items()}
    return {'new_w_point_energy_errors':errors,'new_w_point_energy_maximum':max(errors.values()),
            'local_density_derivative_errors':localerrors,'local_density_maximum':max(v for row in localerrors.values() for v in row.values()),
            'new_reference_point_count':4,'new_reference_momentum_integrations':0,'w_control':.02,
            'kind':'new homogeneous w point-kernel bridge and original local density partial derivatives; no constraint/probe or full quadrature'}
