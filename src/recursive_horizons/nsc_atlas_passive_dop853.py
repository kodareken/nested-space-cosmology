"""Passive DOP853 quadrature from the solver's actual RHS stage calls.

The work value is recorded as a side channel of the existing RHS. It is not
an adaptive ODE component and cannot influence step-size/error control.
The accepted trial uses its eleven internal evaluations and endpoint, plus
the cached old-endpoint work. Rejected trials and old dense-output stages
are ignored by taking the final twelve records at each accepted step.

Stage order and dense coefficients follow the project's installed SciPy
``integrate/_ivp/rk.py`` DOP853 implementation. Accuracy of an arbitrary
passive work function still needs an independent resolution check: the
main ODE controller does not control its quadrature error.
"""
from collections import deque
import math

import numpy as np
from scipy.integrate import DOP853


class PassiveDOP853Quadrature:
    """A bounded scalar work ledger for one live DOP853 solver.

    Call record(work) exactly once inside every actual RHS call. After
    constructing the solver call initialize(solver). Call accepted_step
    immediately after every successful nonzero solver.step(), before dense
    output. For dense values first call solver.dense_output(), then this
    object's dense_value(solver,t). This object never calls the RHS itself.
    Only the last accepted interval has a dense quadrature cache.
    """

    def __init__(self,initial=0.):
        initial=float(initial)
        if not math.isfinite(initial):raise ValueError('finite initial quadrature required')
        self.value=initial
        self._records=deque(maxlen=16)
        self._total=0
        self._first_work=None
        self._solver_id=None
        self._f_old=None
        self._next_t=None
        self._mark=0
        self._nfev=0
        self._interval=None
        self._K=None
        self._I_old=None
        self._delta=None
        self._dense_F=None
        self._dense_mark=None

    def record(self,work):
        """Record one real scalar at each existing solver RHS invocation."""
        work=float(work)
        if not math.isfinite(work):raise ValueError('finite passive work value required')
        if self._first_work is None:self._first_work=work
        self._records.append(work)
        self._total+=1

    def initialize(self,solver):
        """Cache the true initial work, not a step-selection probe value."""
        if self._solver_id is not None:raise RuntimeError('quadrature is already initialized')
        if not isinstance(solver,DOP853) or solver.n_stages!=12 or len(solver.B)!=12:
            raise TypeError('this stage ledger requires DOP853 with twelve stages')
        if len(solver.A_EXTRA)!=3 or np.shape(solver.D)!=(4,16):
            raise TypeError('the DOP853 dense-stage layout is unsupported')
        if solver.t_old is not None:raise RuntimeError('initialize before the first solver step')
        if self._first_work is None or self._total!=solver.nfev:
            raise RuntimeError('record work once for every constructor RHS call')
        self._solver_id=id(solver)
        self._f_old=self._first_work
        self._next_t=float(solver.t)
        self._mark=self._total;self._nfev=solver.nfev
        self._records.clear()
        return self.value

    def _check_solver(self,solver):
        if self._solver_id!=id(solver):raise RuntimeError('quadrature belongs to its initialized solver')

    def accepted_step(self,solver):
        """Accumulate the final accepted trial; no solver field is modified."""
        self._check_solver(solver)
        if solver.status=='failed' or solver.t_old is None or solver.t==solver.t_old:
            raise RuntimeError('a successful nonzero DOP853 step is required')
        interval=(float(solver.t_old),float(solver.t))
        if interval==self._interval:raise RuntimeError('this accepted step was already recorded')
        if interval[0]!=self._next_t:raise RuntimeError('an accepted solver step was skipped')
        added=self._total-self._mark
        if added!=solver.nfev-self._nfev or len(self._records)<12:
            raise RuntimeError('record work exactly once at every actual RHS call')
        if self._interval is None and added%12:
            raise RuntimeError('accepted_step must precede dense output for that step')
        new=np.asarray(tuple(self._records)[-12:],dtype=float)
        K=np.r_[self._f_old,new]
        h=interval[1]-interval[0]
        delta=float(h*np.dot(solver.B,K[:12]))
        self._I_old=self.value
        self.value+=delta
        self._delta=delta;self._K=K
        self._f_old=float(K[-1])
        self._next_t=interval[1];self._interval=interval
        self._mark=self._total;self._nfev=solver.nfev
        self._dense_F=None;self._dense_mark=None
        self._records.clear()
        return self.value

    def dense_value(self,solver,t):
        """Evaluate scalar dense output after the solver built its own output.

        The solver's three extra RHS calls provide the remaining work stages.
        Repeated evaluations reuse the same coefficients and add no RHS calls.
        A second solver.dense_output() call is supported by rebuilding from its
        latest three work records. Earlier intervals are not extrapolated.
        """
        self._check_solver(solver)
        if self._interval is None or (float(solver.t_old),float(solver.t))!=self._interval:
            raise RuntimeError('record the current accepted step before dense evaluation')
        extras=self._total-self._mark
        if extras<3 or extras%3 or extras!=solver.nfev-self._nfev:
            raise RuntimeError('call solver.dense_output() before passive dense evaluation')
        if self._dense_mark!=self._total:
            if len(self._records)<3:raise RuntimeError('dense work stages are missing')
            extended=np.r_[self._K,np.asarray(tuple(self._records)[-3:],dtype=float)]
            h=self._interval[1]-self._interval[0]
            F=np.empty(7)
            F[0]=self._delta
            F[1]=h*self._K[0]-self._delta
            F[2]=2*self._delta-h*(self._K[-1]+self._K[0])
            F[3:]=h*np.dot(solver.D,extended)
            self._dense_F=F;self._dense_mark=self._total
        times=np.asarray(t,dtype=float)
        if not np.all(np.isfinite(times)):raise ValueError('finite interpolation time required')
        lower,upper=sorted(self._interval)
        tolerance=16*np.finfo(float).eps*max(1.,abs(lower),abs(upper))
        if np.any(times<lower-tolerance) or np.any(times>upper+tolerance):
            raise ValueError('passive dense value is only defined on the last accepted interval')
        x=(times-self._interval[0])/(self._interval[1]-self._interval[0])
        value=np.zeros_like(times)
        for i,coefficient in enumerate(self._dense_F[::-1]):
            value+=coefficient
            value*=x if i%2==0 else 1-x
        value+=self._I_old
        return float(value) if value.ndim==0 else value
