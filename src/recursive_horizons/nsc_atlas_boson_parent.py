"""Finite canonical matter probes on the resolved inherited GH history.

The adapter reuses the dense background polynomials and their recorded
quadrature.  It never solves new initial constraints, resets a field state
or interpolates the sparse quarter-interval snapshots.  Its retained
fluctuation sector is S(++), chi(--), with a declared harmonic cutoff.
"""
import math

import numpy as np
from scipy.optimize import brentq

from .nsc_atlas_boson_hamiltonian import quadratic_hamiltonian
from .nsc_atlas_gh_history import DenseGHHistory
from .nsc_atlas_parent_gh import PAIRS,ParentGH
from .nsc_atlas_s3_frame import point_values
from .nsc_atlas_s3_hopf import build_hopf_frame_basis


class ParentBosonProbe:
    """Build G(t) for canonical fluctuations on one recorded parent.

    ``hamiltonian(t)`` returns the q-then-p quadratic matrix for transport.
    ``hamiltonian_details(t)`` also returns lapse, shift and mass diagnostics.
    Canonical modes use physical harmonics r^(3/2)*Y_unit; their momentum
    coefficients represent density relative to the physical reference volume.
    No incoming covariance or quantum state is chosen by this adapter.
    """

    def __init__(self,history,probe_degree,*,decoupled_control=False):
        if not isinstance(history,DenseGHHistory) or not history.bounds:
            raise ValueError('a nonempty DenseGHHistory is required')
        if isinstance(probe_degree,(bool,np.bool_)) or not isinstance(probe_degree,(int,np.integer)) or probe_degree<2:
            raise ValueError('probe_degree must be an integer >=2 for both inherited source modes')
        if not isinstance(decoupled_control,(bool,np.bool_)):
            raise ValueError('decoupled_control must be boolean')
        self.history=history;self.metadata=dict(history.metadata)
        required=('degree','epsilon','gauge','backend','physical_coordinate_endpoint','proper_endpoint')
        if any(key not in self.metadata for key in required):
            raise ValueError('background metadata must declare model, gauge and physical endpoints')
        if self.metadata['backend']!='hopf':
            raise ValueError('this adapter requires a Hopf background to preserve its exact quadrature')
        self.system=ParentGH(self.metadata['degree'],self.metadata['epsilon'],
            gauge=self.metadata['gauge'],backend=self.metadata['backend'],
            quadrature_order=self.metadata.get('quadrature_order'),
            kick_amplitude=self.metadata.get('kick_amplitude',0.),kick_source=self.metadata.get('kick_source',1))
        self.initial_time=float(history.bounds[0][0])
        self.coordinate_endpoint=float(self.metadata['physical_coordinate_endpoint'])
        if not math.isfinite(self.coordinate_endpoint) or not self.initial_time<self.coordinate_endpoint<=history.bounds[-1][1]:
            raise ValueError('physical endpoint must be inside the recorded history')
        state=self._state(self.initial_time)
        if len(state)!=self.system.dof+10:
            raise ValueError('background state size does not match its recorded model')
        self.probe_degree=int(probe_degree);self.decoupled_control=bool(decoupled_control)
        shape=self.system.basis.quadrature_shape
        self.basis=build_hopf_frame_basis(self.probe_degree,(shape[0],shape[1]))
        if not np.array_equal(self.basis.points,self.system.basis.points) or not np.array_equal(self.basis.weights,self.system.basis.weights):
            raise ValueError('probe and background quadratures differ')
        selectedS=np.flatnonzero(np.all(self.basis.characters==[1,1],axis=1))
        selectedChi=np.flatnonzero(np.all(self.basis.characters==[-1,-1],axis=1))
        self.s_indices=np.arange(len(selectedS))
        self.chi_indices=len(selectedS)+np.arange(len(selectedChi))
        self.selected_spatial_indices=np.r_[selectedS,selectedChi]
        self.fields=np.r_[np.zeros(len(selectedS),dtype=int),np.ones(len(selectedChi),dtype=int)]
        self.labels=tuple(self.basis.labels[index] for index in self.selected_spatial_indices)
        self.counts=dict(S=len(selectedS),chi=len(selectedChi),total=len(self.fields))
        self.n=len(self.fields);self.phase_dimension=2*self.n
        columns=np.zeros((len(self.basis.labels),self.n))
        columns[self.selected_spatial_indices,np.arange(self.n)]=1.
        self.Y=self.basis.evaluate(columns)
        self.DY=np.stack([self.system.r*self.basis.evaluate(generator@columns)
                         for generator in self.basis.generators],axis=1)
        self._source_coefficients={};self.source_projection_residuals={}
        points=self.basis.points
        profiles={1:math.sqrt(2)/math.pi*points[:,0],
                  2:math.sqrt(12)/math.pi*points[:,1]*points[:,3]}
        for number,profile in profiles.items():
            coefficients=np.zeros(self.n)
            coefficients[self.chi_indices]=self.Y[:,self.chi_indices].T@(self.basis.weights*profile)
            norm=float(np.linalg.norm(coefficients))
            if not math.isfinite(norm) or norm<1e-12:
                raise ValueError('the inherited source mode is absent from the retained probe sector')
            coefficients/=norm
            residual=self.Y@coefficients-profile
            self.source_projection_residuals[number]=float(math.sqrt(self.basis.weights@(residual*residual)))
            self._source_coefficients[number]=coefficients
        self.proper_start=float(state[self.system.dof+8])
        self.proper_endpoint=float(self.metadata['proper_endpoint'])
        actual_end=float(self._state(self.coordinate_endpoint)[self.system.dof+8])
        if not math.isfinite(self.proper_endpoint) or self.proper_endpoint<=self.proper_start or not math.isclose(
                actual_end,self.proper_endpoint,rel_tol=0,abs_tol=1e-9):
            raise ValueError('proper endpoint does not match the dense physical endpoint')
        self.proper_checkpoints=[]
        for index,target in enumerate(np.linspace(self.proper_start,self.proper_endpoint,5)):
            if index==0:t=self.initial_time
            elif index==4:t=self.coordinate_endpoint
            else:
                t=brentq(lambda coordinate:self._state(coordinate)[self.system.dof+8]-target,
                    self.initial_time,self.coordinate_endpoint,xtol=2e-12)
            self.proper_checkpoints.append(dict(coordinate_t=float(t),
                observer_proper_duration=float(self._state(t)[self.system.dof+8])))
        self.checkpoint_times=np.array([point['coordinate_t'] for point in self.proper_checkpoints])

    def _state(self,t):
        t=float(t)
        if not math.isfinite(t) or t<self.initial_time or t>self.coordinate_endpoint:
            raise ValueError('probe time is outside the physical background interval')
        return self.history.evaluate(t)

    def background_nodes(self,t):
        """Evaluate only metric and background fields; no GH jets or RHS."""
        state=self._state(t);coefficients=self.system.unpack(state)
        selected=np.r_[np.arange(10),[20,21]]
        values=self.system.basis.evaluate(coefficients[selected].T)
        g=np.empty((len(values),4,4))
        for index,(a,b) in enumerate(PAIRS):g[:,a,b]=g[:,b,a]=values[:,index]
        return g,values[:,10:12]

    def hamiltonian_details(self,t):
        g,background=self.background_nodes(t)
        result=quadratic_hamiltonian(g,background,self.Y,self.DY,self.basis.weights,self.fields,
            sstar=self.system.sstar)
        if self.decoupled_control:
            K=result['K'].copy();K[self.fields[:,None]!=self.fields[None,:]]=0.
            G=result['G'].copy();G[:self.n,:self.n]=K
            result.update(K=K,G=G)
        result['decoupled_control']=self.decoupled_control
        return result

    def hamiltonian(self,t):
        """Return the full symmetric G(t); state covariance stays with caller."""
        return self.hamiltonian_details(t)['G']

    def source_kick(self,number):
        """Unit canonical momentum kick with fixed normalized Y1 or Y2 profile."""
        if isinstance(number,(bool,np.bool_)) or not isinstance(number,(int,np.integer)) or number not in (1,2):
            raise ValueError('source mode must be integer1 or2')
        return np.r_[np.zeros(self.n),self._source_coefficients[int(number)]]

    def observer_rows(self,t):
        """Rows on q giving physical S, chi and delta log(g_KK) at the observer."""
        state=self._state(t);point=state[self.system.dof:self.system.dof+4]
        values=self.system.r**1.5*point_values(point,self.labels)
        S=values*(self.fields==0);chi=values*(self.fields==1)
        return dict(S=S,chi=chi,log_g_KK=-7/math.sqrt(self.system.K)*S)
