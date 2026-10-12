"""CF4 evolution of the same prepared Dirac source and its retarded tangent.

Blanes--Moan (2006), equation43: two Gauss nodes and two exponentials.
No midpoint-owner bytes or physics are changed. Source phases stay inside
the augmented matrix and axial derivatives come from the actual node PDE.
"""
from dataclasses import dataclass

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply

from .nsc_exact_phase_prepared_state import (
    ExactPhaseIncoming, _finite_array, _finite_real, _finite_complex,
    _digest_arrays, _fixed_preparation_digest, _require_fixed_source,
    _require_reference_support, _inflow_amplitude, _sample_history_profiles,
    _augmented_generator, _node_pde, _restrict,
)
from .nsc_evolved_incoming_state import (
    AmplitudeOnlyMetric, CachedCompatibleIncomingMetric, EvolvedIncomingState,
    HistoryBinding, _as_compatible_family, _assert_fixed_intrinsic, _coverage,
    _prepare_provider, _reject_legacy_kwargs, rho1_node_index,
)
from .nsc_ks_spacetime_variation import chart_coordinates
from .nsc_transmitting_history_jets import generator_direction

CF4_SOURCE='https://personales.upv.es/~serblaza/2006APNUM.pdf'
A1=(3-2*np.sqrt(3.))/12
A2=(3+2*np.sqrt(3.))/12
C1=.5-np.sqrt(3.)/6
C2=.5+np.sqrt(3.)/6


def _gauss_profiles(provider,family,times,x):
    full=provider.inner if isinstance(provider,AmplitudeOnlyMetric) else provider
    if not isinstance(full,CachedCompatibleIncomingMetric):full=CachedCompatibleIncomingMetric(x,family)
    times=np.asarray(times)
    stages=(times[:-1,None]+np.diff(times)[:,None]*np.array([C1,C2])).ravel()
    metrics=np.array([full.values(t,x) for t in stages])
    directions=np.array([full.log_directions(t,x,m) for t,m in zip(stages,metrics)])
    return stages,metrics,directions,_digest_arrays(stages,x,metrics,directions)


@dataclass(frozen=True)
class CF4Incoming(ExactPhaseIncoming):
    gauss_times: object
    gauss_metrics: object
    gauss_directions: object
    gauss_fingerprint: str

    def __post_init__(self):
        for name in ('gauss_times','gauss_metrics','gauss_directions'):
            object.__setattr__(self,name,_finite_real(getattr(self,name),name))
        super().__post_init__()

    def validate(self):
        """Same source/history invariants with this explicitly different scheme.

        The immutable midpoint owner's validator deliberately rejects CF4.
        This subclass validates the common bindings and both actual Gauss
        stages instead of falsely reporting that it used midpoint evolution.
        """
        self.state.validate();h=self.state.history
        expected=_digest_arrays(h.times,.5*(h.times[:-1]+h.times[1:]),h.spatial_grid,
                                self.sampled_node_metrics,self.sampled_node_directions,
                                self.sampled_midpoint_metrics,self.sampled_midpoint_directions)
        if expected!=self.sampled_history_fingerprint:raise ValueError('node/midpoint history binding changed')
        if h.source_digest!=_digest_arrays(self.state.source_covariance,self.state.column_weights,self.state.source_energies):
            raise ValueError('fixed source binding changed')
        if self.fixed_preparation_digest!=self.diagnostics.get('fixed_preparation_digest'):
            raise ValueError('fixed preparation binding changed')
        nstage=2*(len(h.times)-1);nfam=len(h.amplitudes);n=len(h.spatial_grid)
        if (self.gauss_times.shape!=(nstage,) or self.gauss_metrics.shape!=(nstage,4,n)
                or self.gauss_directions.shape!=(nstage,nfam,4,n)):
            raise ValueError('both actual Gauss stages required for every step')
        expected_times=(h.times[:-1,None]+np.diff(h.times)[:,None]*np.array([C1,C2])).ravel()
        if not np.array_equal(expected_times,self.gauss_times):raise ValueError('CF4 stage times changed')
        if _digest_arrays(self.gauss_times,h.spatial_grid,self.gauss_metrics,self.gauss_directions)!=self.gauss_fingerprint:
            raise ValueError('CF4 Gauss profile binding changed')
        if self.diagnostics.get('CF4_time_integrator') is not True:raise ValueError('CF4 scheme declaration required')
        if self.diagnostics.get('physical_history_status')!='OPEN':raise ValueError('physical history remains OPEN')
        if self.diagnostics.get('continuum_error_bound') is not None:raise ValueError('no continuum bound is supplied here')
        for flag in ('metric_evolution','midpoint_held_forcing','outgoing_momentum_k_equals_minus_E','callback_identity_claimed'):
            if self.diagnostics.get(flag):raise ValueError('unsupported '+flag)
        return True

    def require_history(self,provider,times,x):
        super().require_history(provider,times,x)
        family=_as_compatible_family(provider)
        profiles=_gauss_profiles(provider,family,self.state.history.times,self.state.history.spatial_grid)
        if profiles[3]!=self.gauss_fingerprint:raise ValueError('actual Gauss stage history changed')
        return True


def prepare_cf4_incoming(owner,times,provider,initial_fields,source,**kwargs):
    _reject_legacy_kwargs(kwargs);source=_require_fixed_source(source)
    times=_finite_real(times,'time nodes');x=_finite_real(owner.x,'spatial grid')
    if times.ndim!=1 or len(times)<2 or np.any(np.diff(times)<=0):raise ValueError('increasing time nodes required')
    phi=_finite_complex(initial_fields,'initial columns');rows=2*len(x);ns=len(source.energies)
    if phi.shape!=(rows,ns):raise ValueError('initial/source dimensions differ')
    provider=_prepare_provider(owner,provider);family=_as_compatible_family(provider)
    if family is None:raise ValueError('owned compatible family required')
    rho1=rho1_node_index(x);R,metric1=_assert_fixed_intrinsic(provider,times,x,rho1)
    _,initial_dirs=_require_reference_support(owner,provider,times,x);nd=len(initial_dirs)
    B_dense=_inflow_amplitude(owner,phi);B=csr_matrix(B_dense)
    take=np.array([rho1,len(x)+rho1]);Q0=np.eye(ns,dtype=complex)
    state=np.vstack((phi,np.zeros((nd*rows,ns),complex),Q0))
    F=np.empty((len(times),2,ns),complex);Fz=np.empty_like(F)
    dF=np.empty((nd,len(times),2,ns),complex);dFz=np.empty_like(dF)
    node_fields=np.empty((len(times),rows,ns),complex)
    node_tangents=np.empty((nd,len(times),rows,ns),complex)
    phases=np.empty((len(times),ns,ns),complex);max_flux=max_dflux=0.

    def node(index):
        nonlocal max_flux,max_dflux
        field=state[:rows];jets=state[rows:(nd+1)*rows].reshape(nd,rows,ns);phase=state[(nd+1)*rows:]
        metric=provider.values(times[index],x);dirs=provider.log_directions(times[index],x,metric)
        _,rhs,drhs,checks,dflux=_node_pde(owner,metric,dirs,field,jets,phase,B_dense)
        F[index],Fz[index],dF[:,index],dFz[:,index]=_restrict(R,take,field,rhs,jets,drhs)
        node_fields[index]=field;node_tangents[:,index]=jets;phases[index]=phase
        max_flux=max(max_flux,checks['norm_flux_algebra']);max_dflux=max(max_dflux,dflux)

    node(0)
    for step,(left,right) in enumerate(zip(times[:-1],times[1:])):
        matrices=[]
        for c in (C1,C2):
            t=left+c*(right-left);metric=provider.values(t,x);L,checks=owner.generator(metric)
            dirs=provider.log_directions(t,x,metric);dLs=[]
            for d in dirs:
                dL,error=generator_direction(owner,metric,d);dLs.append(dL);max_dflux=max(max_dflux,error)
            matrices.append(_augmented_generator(L,dLs,B,source.energies))
            max_flux=max(max_flux,checks['norm_flux_algebra'])
        early=A2*matrices[0]+A1*matrices[1]
        late=A1*matrices[0]+A2*matrices[1]
        for matrix in (early,late):
            scaled=(right-left)*matrix
            state=expm_multiply(scaled,state,traceA=scaled.diagonal().sum())
        if not np.isfinite(state).all():raise ArithmeticError('nonfinite CF4 state')
        node(step+1)

    profiles=_sample_history_profiles(provider,family,times,x)
    gt,gm,gd,gh=_gauss_profiles(provider,family,times,x)
    S1=float(chart_coordinates(1.)[1])
    history=HistoryBinding(tuple(family.amplitudes),tuple(d.inner_radius for d in family.directions),
                           tuple(d.outer_radius for d in family.directions),times,x,rho1,metric1,R,0.,S1,source.digest)
    preparation=_fixed_preparation_digest(source,phi,x,owner.mass,owner.angular,times[0])
    expected=np.eye(ns)[None]*np.exp(-1j*(times-times[0])[:,None,None]*source.energies)
    common={'norm_flux_residual':max_flux,'metric_flux_tangent_residual':max_dflux,
            'fixed_preparation_digest':preparation,'physical_history_status':'OPEN',
            'physical_preparation_status':'OPEN','source_covariance_tangent_explicit':0.,
            'restriction_parameter_derivative':0.,'exact_harmonic_phase':True,
            'CF4_time_integrator':True,'midpoint_held_forcing':False,
            'outgoing_momentum_k_equals_minus_E':False,'callback_identity_claimed':False,
            'metric_evolution':False,'continuum_error_bound':None,'full_source_error_bound':None}
    evolved=EvolvedIncomingState(times+S1,F,dF,source.covariance,source.column_weights,source.energies,
                                 history,_coverage(ns,len(np.unique(source.energies))),common)
    diagnostics={**common,'phase':_finite_complex(phases,'phase'),
                 'expected_phase':_finite_complex(expected,'expected phase'),
                 'exact_harmonic_phase_residual':float(np.max(abs(phases-expected))),
                 'node_fields':_finite_complex(node_fields,'node fields'),
                 'node_tangents':_finite_complex(node_tangents,'node tangents'),
                 'CF4_coefficients':[float(A1),float(A2)],'CF4_Gauss_nodes':[float(C1),float(C2)],
                 'CF4_source':CF4_SOURCE,'steps':len(times)-1,
                 'temporal_prescription':'Blanes-Moan equation43, two actual Gauss stages and two augmented exponentials',
                 'initial_and_incident_tangents':'fixed supplied preparation; earlier retarded support is a caller contract'}
    return CF4Incoming(evolved,Fz,dFz,profiles['node_metrics'],profiles['midpoint_metrics'],
                       profiles['node_directions'],profiles['midpoint_directions'],profiles['fingerprint'],
                       preparation,diagnostics,gt,gm,gd,gh)
