"""Directed subgap source covariance transport in the exact metric.

The three real Bloch components represent one complete 2x2 source covariance
at a fixed energy. This is not a closure of the evolved field's energy or
spatial dependence. The original source columns are never replaced.
"""
from math import comb, factorial

import numpy as np
from flint import acb,arb,arb_series,ctx
from scipy.integrate import solve_ivp

from .nsc_massive_jost_phase_bound import _series_context


def _coefficients(q,terms):
    c=arb(3)/2*(2*q).sin()+arb(1)/2*(2*q).cos()
    s=arb(3)/2*(2*q).cos()-arb(1)/2*(2*q).sin()
    return [3-2*s]+[-arb(2)**(n+1)/factorial(n+1)*(
        c*(0,-1,0,1)[n%4]+s*(1,0,-1,0)[n%4]) for n in range(1,terms)]


def _polynomial(coefficients,x,order):
    value=arb_series([0],prec=order+1)
    for coefficient in reversed(coefficients):
        value=value*x+coefficient
    return value


def metric_H(q,u,order,*,terms=48):
    """H(u)=-W(q_h+u)/u, with the entire omitted sine/cosine series.

    Coefficient |H_n|<=4*2^n/(n+1)! for n>=1. Derivative tail majorants are
    composed into the supplied u-series before being added to its jets.
    """
    if type(terms) is not int or terms <= order or terms<8:
        raise ValueError('metric tail order must exceed jet order')
    U=u[0].abs_upper()
    if not U<arb('0.4'):
        raise ValueError('source horizon metric domain requires |u|<0.4')
    result=_polynomial(_coefficients(q,terms),u,order)
    increment=arb_series([0]+[u[j].abs_upper() for j in range(1,order+1)],prec=order+1)
    power=arb_series([1],prec=order+1)
    tail=arb_series([0],prec=order+1)
    x=2*U
    for j in range(order+1):
        k=terms-j
        exp_tail=x**k/factorial(k)/(1-x/(k+1))
        upper=4*arb(2)**j/factorial(j)/(terms+1)*exp_tail
        tail+=upper*power
        power*=increment
    return arb_series([result[j]+arb(0,tail[j].abs_upper())
                       for j in range(order+1)],prec=order+1)


class BlochSource:
    def __init__(self,q_h,energy,mass,angular,*,bits=192,metric_terms=48):
        if type(bits) is not int or bits<80:
            raise ValueError('directed precision >=80 required')
        self.bits,self.metric_terms=bits,metric_terms
        with ctx.workprec(bits):
            self.q,self.energy,self.mass,self.angular=map(arb,(q_h,energy,mass,angular))
            if not all(v.is_finite() for v in (self.q,self.energy,self.mass,self.angular)):
                raise ValueError('finite source channel required')
            if not self.energy>0 or not self.mass>self.energy:
                raise ValueError('this source owner requires the positive subgap channel')
        self.numeric_coefficients=np.array([float(x) for x in _coefficients(self.q,metric_terms)])

    def hamiltonian_series(self,y,order):
        delta=y.exp()
        H=metric_H(self.q,-delta,order,terms=self.metric_terms)
        if not H[0]>0:
            raise ArithmeticError('positive interior metric not enclosed')
        radius=1/(self.q-delta).sin()
        scale=(delta/H).sqrt()
        return (-self.mass*radius*scale,self.angular*scale,-self.energy/H)

    def rhs_numeric(self,y,n):
        delta=np.exp(y)
        H=np.polynomial.polynomial.polyval(-delta,self.numeric_coefficients)
        radius=1/np.sin(float(self.q)-delta)
        scale=np.sqrt(delta/H)
        h=np.array([-float(self.mass)*radius*scale,float(self.angular)*scale,-float(self.energy)/H])
        return 2*np.cross(h,n)


def _cross(left,right):
    return [left[1]*right[2]-left[2]*right[1],
            left[2]*right[0]-left[0]*right[2],
            left[0]*right[1]-left[1]*right[0]]


def initial_bloch(frame,reflection,kappa_source,y_start):
    """Full coherent affine source in the direct compact-coordinate frame.

    The extra legacy interior phase is diag(exp(-i*pi/4),exp(i*pi/4)). It
    turns the sewn covariance's i*s*conj(R) into s*conj(R). The common
    exterior/interior phase cancels exactly with the transformed reflection.
    """
    with ctx.workprec(frame.bits):
        kappa_source=arb(kappa_source)
        if not kappa_source>0:
            raise ValueError('fixed positive source temperature required')
        F,_=frame.evaluate(arb(y_start).exp(),interior=True)
        x=arb.pi()*frame.energy/kappa_source
        f=1/(1+(2*x).exp())
        coherence=1/(2*x.cosh())
        R=acb(reflection)
        if not abs(R).contains(1):
            raise ValueError('certified subgap reflection must enclose its unit-circle value')
        # Zero exterior current proves the TRUE subgap reflection is unit.
        # Enclose its phase; do not normalize any original numerical columns.
        phase=R.arg()
        if not phase.is_finite():
            raise ArithmeticError('subgap reflection phase is unresolved')
        R=acb(0,phase).exp()
        C=[[acb(1-f),coherence*R.conjugate()],
           [coherence*R,acb(f)]]
        Q=[[sum((F[i][a]*C[a][b]*F[j][b].conjugate()
                  for a in range(2) for b in range(2)),acb(0))
             for j in range(2)] for i in range(2)]
        return [2*Q[0][1].real,-2*Q[0][1].imag,(Q[0][0]-Q[1][1]).real]


def capture_bloch(model,initial,y_start,y_end,*,max_step=.05):
    numeric=np.array([float(v.mid()) for v in initial])
    run=solve_ivp(model.rhs_numeric,(float(y_start),float(y_end)),numeric,
                  method='DOP853',rtol=2e-13,atol=2e-15,max_step=max_step,dense_output=True)
    if not run.success:
        raise ArithmeticError(run.message)
    # Preserve anchored binary endpoints and the six Hairer correction rows.
    rows=[]
    for j,polynomial in enumerate(run.sol.interpolants):
        rows.append(np.r_[run.t[j],run.t[j+1],run.y[:,j],run.y[:,j+1],polynomial.F[1:].ravel()])
    return np.asarray(rows),int(run.nfev)


def _anchored_coefficients(row):
    if len(row)!=26 or not np.isfinite(row).all() or not row[1]>row[0]:
        raise ValueError('complete increasing anchored three-component cell required')
    out=[]
    corrections=np.asarray(row[8:]).reshape(6,3)
    for component in range(3):
        values=[arb(0) for _ in range(8)]
        values[0]=arb(float(row[2+component]))
        values[1]=arb(float(row[5+component]))-values[0]
        for j in range(1,7):
            p,q=(j+2)//2,(j+1)//2
            for k in range(q+1):
                values[p+k]+=arb(float(corrections[j-1,component]))*(-1)**k*comb(q,k)
        out.append(values)
    return out


def cell_defect(model,row,*,degree=12):
    """Enclose integral ||p_x-h A(y)p|| over x in [0,1]. No extra h."""
    if type(degree) is not int or degree<7 or degree>=model.metric_terms:
        raise ValueError('resolved degree-7 interpolant and metric tail required')
    with _series_context(model.bits,degree):
        coeffs=_anchored_coefficients(row)
        h=arb(float(row[1]))-arb(float(row[0]))
        def residual(center):
            x=arb_series([center,1],prec=degree+1)
            y=arb(float(row[0]))+h*x
            H=model.hamiltonian_series(y,degree)
            p=[_polynomial(c,x,degree) for c in coeffs]
            px=[_polynomial([j*c[j] for j in range(1,8)],x,degree) for c in coeffs]
            motion=_cross(H,p)
            return [px[j]-2*h*motion[j] for j in range(3)]
        point=residual(arb(1)/2)
        whole=residual(arb('0.5','0.5'))
        norms=[sum((p[j].abs_upper()/arb(2)**j for j in range(degree)),arb(0))
               +w[degree].abs_upper()/arb(2)**degree for p,w in zip(point,whole)]
        bound=sum((v*v for v in norms),arb(0)).sqrt().upper()
        if not bound.is_finite():
            raise ArithmeticError('nonfinite Bloch defect bound')
        return bound


def validate_bloch(model,initial,trace,y_target,*,degree=12):
    with ctx.workprec(model.bits):
        trace=np.asarray(trace,float)
        if trace.ndim!=2 or trace.shape[1]!=26 or not len(trace):
            raise ValueError('complete Bloch trajectory required')
        if not np.array_equal(trace[:-1,1],trace[1:,0]) or not np.array_equal(trace[:-1,5:8],trace[1:,2:5]):
            raise ValueError('Bloch cells do not join at identical endpoints')
        error=sum(((v-arb(float(n))).abs_upper()**2 for v,n in zip(initial,trace[0,2:5])),arb(0)).sqrt()
        initial_error=error.upper()
        defect=arb(0)
        for row in trace:
            defect+=cell_defect(model,row,degree=degree)
        # Exact endpoint may differ from the numerical log/exp round trip.
        end=arb(float(trace[-1,1]));target=arb(y_target)
        with _series_context(model.bits,0):
            H=model.hamiltonian_series(arb_series([end.union(target)],prec=1),0)
            generator=2*sum((v[0].abs_upper()**2 for v in H),arb(0)).sqrt()
        initial_norm=sum((v.abs_upper()**2 for v in initial),arb(0)).sqrt()
        bridge=(generator*initial_norm*(target-end).abs_upper()).upper()
        total=(error+defect+bridge).upper()
        return {'initial_error':initial_error,'defect_integral':defect.upper(),
                'endpoint_bridge':bridge,'bloch_error':total,
                'endpoint':tuple(arb(float(v)) for v in trace[-1,5:8]),
                'cells':len(trace),'bits':model.bits,'physical_local_gate':'OPEN'}


def original_covariance_error(prepared_columns,source_covariance,validated):
    """Operator-norm error against the unchanged numerical 2x3 source data."""
    return original_covariance_distance_bounds(prepared_columns,source_covariance,validated)['upper']


def original_covariance_distance_bounds(prepared_columns,source_covariance,validated):
    """Two-sided distance control, including a lower bound for mutations."""
    with ctx.workprec(validated['bits']):
        return _original_covariance_error(prepared_columns,source_covariance,validated)


def _original_covariance_error(prepared_columns,source_covariance,validated):
    F=np.asarray(prepared_columns,complex);C=np.asarray(source_covariance,complex)
    if F.shape!=(2,3) or C.shape!=(3,3) or not np.isfinite(F).all() or not np.isfinite(C).all():
        raise ValueError('original finite 2x3 columns and 3x3 covariance required')
    if not np.array_equal(C,C.conj().T):
        raise ValueError('the retained source covariance must be Hermitian')
    if not validated['bloch_error']>=0:
        raise ValueError('a nonnegative directed Bloch error is required')
    f=[[acb(v) for v in row] for row in F]
    c=[[acb(v) for v in row] for row in C]
    Q=[[sum((f[i][a]*c[a][b]*f[j][b].conjugate() for a in range(3) for b in range(3)),acb(0))
        for j in range(2)] for i in range(2)]
    trace=(Q[0][0]+Q[1][1]).real
    n=[2*Q[0][1].real,-2*Q[0][1].imag,(Q[0][0]-Q[1][1]).real]
    differences=[v-q for v,q in zip(validated['endpoint'],n)]
    upper_norm=sum((v.abs_upper()**2 for v in differences),arb(0)).sqrt()
    lower_norm=sum((v.abs_lower()**2 for v in differences),arb(0)).sqrt()
    upper=((abs(1-trace)+upper_norm+validated['bloch_error'])/2).upper()
    lower=max(arb(0),((abs(1-trace).lower()+lower_norm-validated['bloch_error'])/2).lower())
    return {'lower':lower,'upper':upper}
