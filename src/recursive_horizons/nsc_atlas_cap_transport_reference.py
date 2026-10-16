"""Independent mass/lapse/residual integration on one smooth cap piece.

The caller owns the field jets, endpoint traces and interface inventory.
This component never reads the native mass tangent or its radial map error.
"""
import numpy as np
from scipy.integrate import solve_ivp

from .nsc_atlas_background_cap import _budget
from .nsc_atlas_polar_cap import _scalar


STATE_KEYS = ('n','delta','Z_f','Z_Pi','Z_Phi','Z_jump')
FIELD_KEYS = frozenset(('f_x','Pi','Phi','Pi_x','Phi_x','Pi_r','Phi_r','Pdiff','Pdiff_f'))


def _density(r, n, delta, fields, alpha, tx):
    """Six regular RHS components; r>0, with weighted Pi residual."""
    pi,phi=fields['Pi'],fields['Phi']
    exp_delta=np.exp(delta)
    if not np.isfinite(exp_delta) or exp_delta<=0:
        raise ValueError('reference lapse exponential must be finite and positive')
    energy=pi*pi+phi*phi
    a=alpha*r*energy/2
    b=r*r*(fields['Pdiff']/2+(1-r*r)*energy/4)
    rho=b-a*n
    A=1-r*r-2*alpha*n/r
    Ar=-2*r-2*alpha*(rho/r-n/(r*r))
    B=A*exp_delta
    Br=exp_delta*(Ar+A*a)
    rf=fields['f_x']-tx*B*pi
    # Multiply the Pi equation by r^2 before evaluating; no singular
    # 2/r or 1/r^2 field residual is needed by the mass identity.
    r2rp=r*r*fields['Pi_x']-tx*(2*r*B*phi+r*r*(Br*phi+B*fields['Phi_r']))
    r2rp+=tx*r*r*exp_delta*fields['Pdiff_f']
    rq=fields['Phi_x']-tx*(Br*pi+B*fields['Pi_r'])
    sources=np.array([r*r*fields['Pdiff_f']*rf/2,A*pi*r2rp/2,r*r*A*phi*rq/2])
    if not np.all(np.isfinite([rho,a,A,*sources])):
        raise ValueError('nonfinite reference residual density')
    return rho,a,A,sources


def transport_radial_piece(fields, *, left, right, initial, alpha, clock_rate,
                           budget, rtol=2e-11, atol=2e-13):
    """Integrate n,delta and four attenuated channels on a declared piece.

    fields(r) supplies finite scalar fixed-radius x/r jets and the smooth
    potential/derivative. The caller charges that evaluation. A regular
    centre starts with n=Z_i=0 and supplied delta; its RHS is exactly zero,
    without querying axis field jets. A zero-width piece is the identity;
    its neighboring interface jumps are the caller's responsibility.

    No division by A occurs. The minimum A is over RHS samples only; neither
    a positive sampled margin nor requested ODE tolerances enclose a domain.
    Each RHS charges 24 algebraic components plus the caller's field work.
    """
    budget=_budget(budget)
    left,right=(_scalar(v,k,nonnegative=True) for k,v in (('left',left),('right',right)))
    alpha=_scalar(alpha,'alpha',nonnegative=True); tx=_scalar(clock_rate,'clock rate')
    rtol=_scalar(rtol,'rtol',positive=True); atol=_scalar(atol,'atol',positive=True)
    if right<left or right>=1:
        raise ValueError('transport piece requires 0 <= left <= right < 1')
    if not 1e-13<=rtol<=1e-5 or not 1e-15<=atol<=1e-7:
        raise ValueError('transport solver tolerances outside the declared range')
    if not callable(fields) or not isinstance(initial,dict) or set(initial)!=set(STATE_KEYS):
        raise ValueError('field callback and exactly six initial components required')
    values=np.array([_scalar(initial[k],k) for k in STATE_KEYS])
    if left==0 and np.any(values[[0,2,3,4,5]]!=0):
        raise ValueError('regular centre requires zero mass and accumulated channels')
    before=budget.evaluations; calls=rhs_calls=0
    budget.check(8)
    minimum_A=1. if left==0 else 1-left*left-2*alpha*values[0]/left
    width=right-left

    def rhs(u,state):
        nonlocal calls,rhs_calls,minimum_A
        budget.check(24); rhs_calls+=1
        r=left if u==0 else right if u==1 else left+u*width
        if r==0:
            return np.zeros(6)
        row=fields(float(r)); calls+=1
        if not isinstance(row,dict) or not FIELD_KEYS<=set(row):
            raise ValueError('reference field callback lacks required jets')
        row={k:_scalar(row[k],k) for k in FIELD_KEYS}
        with np.errstate(over='raise',divide='raise',invalid='raise'):
            rho,a,A,sources=_density(r,state[0],state[1],row,alpha,tx)
            out=width*np.r_[rho,a,sources-a*state[2:5],-a*state[5]]
        minimum_A=min(minimum_A,A)
        if not np.all(np.isfinite(out)):
            raise ValueError('nonfinite reference transport RHS')
        return out

    if width:
        solution=solve_ivp(rhs,(0.,1.),values,method='DOP853',rtol=rtol,atol=atol,max_step=.25)
        if not solution.success or solution.t[-1]!=1.:
            raise RuntimeError('reference transport integration failed: '+solution.message)
        values=solution.y[:,-1]
    if not np.all(np.isfinite(values)) or not np.isfinite(minimum_A):
        raise ValueError('nonfinite reference transport output')
    budget.check()
    return dict(zip(STATE_KEYS,map(float,values)),coefficient_calls=calls,rhs_calls=rhs_calls,
        evaluations=budget.evaluations-before,minimum_sampled_A=float(minimum_A),
        sampled_A_positive=bool(minimum_A>0),born=width==0,
        requested_rtol=rtol,requested_atol=atol,uniform_error_enclosure=False)
