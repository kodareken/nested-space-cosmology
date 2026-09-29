"""Continuous proof of a correction to the homogeneous vacuum expansion.

The original preparation is an input, not overwritten by this proof witness.
The correction's generator is skew; the finite expansion defect is retained
as an explicit forcing, avoiding subtraction of large oscillatory columns.
"""
import time

import numpy as np
from flint import arb,arb_series,ctx
from scipy.integrate import solve_ivp

from .nsc_massive_jost_phase_bound import _series_context
from .nsc_subgap_source_covariance import (
    BlochSource,_anchored_coefficients,_coefficients,_cross,_polynomial,
)
from .nsc_vacuum_source_remainder import VacuumSourceExpansion,_real


class VacuumCorrection:
    def __init__(self,expansion,energy):
        if not isinstance(expansion,VacuumSourceExpansion):
            raise TypeError('owned homogeneous vacuum expansion required')
        self.expansion=expansion
        self.bits=expansion.bits;self.metric_terms=expansion.metric_terms
        self.q=expansion.q;self.mass=expansion.mass;self.angular=expansion.angular
        with ctx.workprec(self.bits):
            self.energy=_real(energy,'energy')
            if not self.energy.is_finite() or not self.energy>0:
                raise ValueError('positive finite energy required')
            self.numeric_coefficients=np.array([float(v) for v in _coefficients(self.q,self.metric_terms)])

    def hamiltonian_series(self,y,order):
        # The existing local generator is valid at every positive E. Only
        # BlochSource's separate subgap initializer requires E<m.
        return BlochSource.hamiltonian_series(self,y,order)

    def defect_series(self,y,order):
        """Exact expansion defect (-2 i t P_(M+1)/(H E^M),0), in y jets."""
        with _series_context(self.bits,order):
            delta=y.exp()
            H,ps,_=self.expansion.jets(delta[0],remainder_jet_order=max(1,order))
            increment=arb_series([0]+[delta[j] for j in range(1,order+1)],prec=order+1)
            def compose(series):
                result=arb_series([0],prec=order+1)
                for j in range(order,-1,-1):result=result*increment+series[j]
                return result
            real,imag=(compose(p) for p in ps[-1])
            factor=2*delta.sqrt()/(compose(H)*self.energy**self.expansion.order)
            return (factor*imag,-factor*real,arb_series([0],prec=order+1))

    def rhs_numeric(self,y,error):
        rotation=BlochSource.rhs_numeric(self,y,error)
        with _series_context(self.bits,0):
            defect=self.defect_series(arb_series([arb(float(y))],prec=1),0)
            forcing=np.array([float(v[0].mid()) for v in defect])
        return rotation-forcing

    def expansion_at_log_distance(self,y):
        with ctx.workprec(self.bits):
            delta=arb(y).exp()
            _,ps,zs=self.expansion.jets(delta)
            E=self.energy;M=self.expansion.order
            transverse=[delta.sqrt()*sum((ps[j][k][0]/E**j for j in range(1,M+1)),arb(0))
                        for k in (0,1)]
            return (*transverse,sum((zs[j][0]/E**j for j in range(M+1)),arb(0)))

    def initial_error(self,y_start):
        """Nonzero bound inherited from the exact horizon limit, never set to zero."""
        with ctx.workprec(self.bits):
            delta=arb(y_start).exp()
            rho=(self.q-delta-arb.pi()/2).tan()
            constant=self.expansion.remainder_constant(rho,cells=8)
            return (constant/self.energy**self.expansion.order).upper()


def capture_correction(model,y_start,y_end,*,rtol=1e-9,atol=1e-21,max_step=.1,cpu_limit=30.):
    started=time.process_time()
    if (not np.isfinite([y_start,y_end,rtol,atol,max_step,cpu_limit]).all()
            or y_end<=y_start or min(rtol,atol,max_step,cpu_limit)<=0):
        raise ValueError('increasing finite times and explicit positive controls required')
    def rhs(y,error):
        if time.process_time()-started>cpu_limit:
            raise TimeoutError('vacuum-correction proof-witness CPU budget exhausted')
        return model.rhs_numeric(y,error)
    run=solve_ivp(rhs,(float(y_start),float(y_end)),np.zeros(3),method='DOP853',
                  rtol=rtol,atol=atol,max_step=max_step,dense_output=True)
    if not run.success:raise ArithmeticError(run.message)
    rows=[np.r_[run.t[j],run.t[j+1],run.y[:,j],run.y[:,j+1],poly.F[1:].ravel()]
          for j,poly in enumerate(run.sol.interpolants)]
    return np.asarray(rows),int(run.nfev)


def correction_cell_defect(model,row,*,degree=12):
    """Integral defect on normalized x: p_x-2 h hvec cross p+h r_M."""
    if type(degree) is not int or degree<7:
        raise ValueError('at least degree seven for the anchored interpolant required')
    with _series_context(model.bits,degree):
        coefficients=_anchored_coefficients(row)
        width=arb(float(row[1]))-arb(float(row[0]))
        def residual(center):
            x=arb_series([center,1],prec=degree+1)
            y=arb(float(row[0]))+width*x
            h=model.hamiltonian_series(y,degree)
            forcing=model.defect_series(y,degree)
            p=[_polynomial(c,x,degree) for c in coefficients]
            px=[_polynomial([j*c[j] for j in range(1,8)],x,degree) for c in coefficients]
            rotation=_cross(h,p)
            return [px[j]-2*width*rotation[j]+width*forcing[j] for j in range(3)]
        middle=residual(arb(1)/2);whole=residual(arb('0.5','0.5'))
        errors=[sum((p[j].abs_upper()/arb(2)**j for j in range(degree)),arb(0))
                +q[degree].abs_upper()/arb(2)**degree for p,q in zip(middle,whole)]
        result=sum((v*v for v in errors),arb(0)).sqrt().upper()
        if not result.is_finite():raise ArithmeticError('nonfinite correction defect')
        return result


def validate_correction(model,trace,y_start,rho_target,*,degree=12):
    trace=np.asarray(trace,float)
    if (trace.ndim!=2 or trace.shape[1]!=26 or not len(trace)
            or not np.isfinite(trace).all() or trace[0,0]!=float(y_start)):
        raise ValueError('complete correction trajectory at its declared start required')
    if (not np.array_equal(trace[:-1,1],trace[1:,0])
            or not np.array_equal(trace[:-1,5:8],trace[1:,2:5])):
        raise ValueError('correction cells must join at identical binary endpoints')
    with ctx.workprec(model.bits):
        seed=model.initial_error(float(y_start))
        initial=(seed+sum((arb(float(v))**2 for v in trace[0,2:5]),arb(0)).sqrt()).upper()
        defect=sum((correction_cell_defect(model,row,degree=degree) for row in trace),arb(0)).upper()
        end=arb(float(trace[-1,1]));target=model.expansion.target_distance(rho_target).log()
        with _series_context(model.bits,0):
            h=model.hamiltonian_series(arb_series([end.union(target)],prec=1),0)
            # The exact VACUUM Bloch vector has norm one. No constant norm
            # is assumed for the driven correction or finite-occupation state.
            bridge=(2*sum((v[0].abs_upper()**2 for v in h),arb(0)).sqrt()
                    *(target-end).abs_upper()).upper()
        approximate=model.expansion_at_log_distance(end)
        endpoint=tuple(n+arb(float(v)) for n,v in zip(approximate,trace[-1,5:8]))
        return {'bits':model.bits,'endpoint':endpoint,'bloch_error':(initial+defect+bridge).upper(),
                'initial_error':initial,'normalized_defect_integral':defect,
                'endpoint_bridge':bridge,'cells':len(trace),'physical_local_gate':'OPEN'}
