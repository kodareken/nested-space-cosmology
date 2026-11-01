"""One declared counterrotating preparation and its coupled initial constraints.

chi=A X1, Pi_chi=lambda*f/psi6, f=H^1.5 sqrt(12)/pi X2 X4;
K_reference=lambda*psi^-2*Abar, trK=0, freely supplied conformal TT=0.
The amplitude path is preparation data, not time evolution or selection by
inheritance. Old mixed-source initializers are not modified or invoked.
"""
from dataclasses import dataclass,field
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.linalg import lu_factor,lu_solve
from scipy.sparse.linalg import LinearOperator,gmres
from scipy.special import eval_gegenbauer
from threadpoolctl import threadpool_limits

from .nsc_atlas_parent_circulation import Parent
from .nsc_atlas_valve_neck import solve_necks
from .nsc_atlas_finite_parent_data import FiniteParentData,positivity_lower_bound
from .nsc_atlas_kick_york import solve_longitudinal,inherited_polynomial_inverse,YorkKick,tensor_divergence
from .nsc_atlas_s3_frame import build_frame_basis,arbitrary_points_values
from .nsc_atlas_s3_spectral import SurfaceBasis,harmonic_values,_nodes,_quadrature_from_nodes,_factors,_radial_normalization


def _points(points):
    x=np.asarray(points,dtype=float)
    if x.ndim!=2 or x.shape[1]!=4 or np.any(~np.isfinite(x)) or np.any(abs(np.sum(x*x,axis=1)-1)>2e-10):
        raise ValueError('finite unit-sphere points with shape(N,4) required')
    return x


def rotated_points(points):
    x=_points(points)
    return np.column_stack([-x[:,3],x[:,2],-x[:,1],x[:,0]])


@dataclass
class CounterrotatingData:
    parent: object
    coefficients: np.ndarray
    labels: tuple
    A: float
    kick_amplitude: float
    tensor_coefficients: np.ndarray
    tensor_labels: tuple
    metadata: dict=field(default_factory=dict)
    provenance: dict=field(default_factory=dict)

    def __post_init__(self):
        self.labels=tuple(tuple(int(v) for v in label) for label in self.labels)
        self.tensor_labels=tuple((int(l),int(r),int(s),str(kind)) for l,r,s,kind in self.tensor_labels)
        self.coefficients=np.asarray(self.coefficients,dtype=float).copy()
        self.tensor_coefficients=np.asarray(self.tensor_coefficients,dtype=float).copy()
        if (self.coefficients.shape!=(len(self.labels),) or self.tensor_coefficients.shape!=(3,3,len(self.tensor_labels)) or
                np.any(~np.isfinite(self.coefficients)) or np.any(~np.isfinite(self.tensor_coefficients)) or
                not math.isfinite(self.A) or self.A<=0 or not math.isfinite(self.kick_amplitude) or not 0<=self.kick_amplitude<=1):
            raise ValueError('finite coefficients and declared positive A,lambda in[0,1] required')
        self.coefficients.setflags(write=False);self.tensor_coefficients.setflags(write=False)

    @property
    def degree(self):return max(label[0] for label in self.labels)

    @property
    def epsilon(self):return 1.

    @property
    def york_basis(self):
        basis=build_frame_basis(4,7)
        if basis.labels!=self.tensor_labels:raise ValueError('York labels do not match the declared polynomial frame')
        return basis

    def values(self,points,*,chunk=8192):
        x=rotated_points(points);result=np.empty(len(x))
        for lo in range(0,len(x),chunk):result[lo:lo+chunk]=harmonic_values(x[lo:lo+chunk],self.labels)@self.coefficients
        return result

    def tensor(self,points,*,chunk=8192):
        """Unit-lambda Abar in the ORIGINAL right reference frame."""
        x=_points(points);result=np.empty((len(x),3,3));flat=self.tensor_coefficients.reshape(9,-1).T
        for lo in range(0,len(x),chunk):
            result[lo:lo+chunk]=(arbitrary_points_values(x[lo:lo+chunk],self.tensor_labels)@flat).reshape(-1,3,3)
        return result

    def f(self,points):
        x=_points(points)
        return self.parent.H**1.5*math.sqrt(12)/math.pi*x[:,1]*x[:,3]

    def chi(self,points):return self.A*_points(points)[:,0]

    def momentum_profile(self,points):return self.f(points)

    def squared_norm(self,points):return np.sum(self.tensor(points)**2,axis=(1,2))

    @classmethod
    def from_record(cls,recordpath):
        """SHA-authenticated state load, with no solver or source substitution.

        Current captured source pins are checked. A changed implementation
        requires its declared replay/regeneration rather than silent reuse.
        """
        path=Path(recordpath).resolve();record=json.loads(path.read_text(encoding='utf-8'))
        if record.get('schema')!='nsc-atlas-valve-counterrotation-v1':raise ValueError('unsupported counterrotation record schema')
        root=(path.parent/record['repository_root_relative_to_record']).resolve()
        def member(name):
            p=(root/name).resolve()
            if not p.is_relative_to(root):raise ValueError('record member leaves the declared repository')
            return p
        pins=record['provenance']['source_sha256_after']
        if pins!=record['provenance']['source_sha256_before'] or not record['provenance']['source_unchanged']:
            raise ValueError('mixed-source counterrotation record')
        for name,digest in pins.items():
            if hashlib.sha256(member(name).read_bytes()).hexdigest()!=digest:raise ValueError('counterrotation source SHA mismatch: '+name)
        raw=record['raw_cache'];cache=member(raw['path'])
        if hashlib.sha256(cache.read_bytes()).hexdigest()!=raw['sha256']:raise ValueError('counterrotation raw-cache SHA mismatch')
        with np.load(cache,allow_pickle=False) as arrays:
            meta=json.loads(str(arrays['metadata']))
            if meta['source_sha256']!=pins or meta['preparation_id']!='one-mode-counterrotation-v1':
                raise ValueError('counterrotation cache metadata mismatch')
            if meta['data_metadata']!=record['metadata']:
                raise ValueError('counterrotation record/cache diagnostics mismatch')
            p=Parent(float(arrays['vacuum']),float(arrays['mass_squared']),float(arrays['offset']),float(arrays['kappa']))
            labels=tuple(map(tuple,arrays['labels'].tolist()))
            tensor_labels=tuple((*map(int,lrs),str(kind)) for lrs,kind in zip(arrays['tensor_lrs'],arrays['tensor_kind']))
            result=cls(p,arrays['coefficients'],labels,float(arrays['A']),float(arrays['kick_amplitude']),
                       arrays['tensor_coefficients'],tensor_labels,record['metadata'],record['provenance'])
        prep=record['preparation']
        if result.A!=prep['A'] or result.kick_amplitude!=prep['kick_amplitude']:
            raise ValueError('record/cache preparation parameter mismatch')
        return result


def original_points(z):
    # Inverse of the frame-preserving left map Z=(-X4,X3,-X2,X1).
    return np.column_stack([z[:,3],-z[:,2],z[:,1],-z[:,0]])

def _make_basis(degree=40,max_j=4,orders=(141,13,25)):
    labels=tuple((l,j,m) for l in range(0,degree+1,2)
                 for j in range(0,min(l,max_j)+1,2) for m in range(0,j+1,2))
    u,wu,v,wv,phi=_nodes(*orders);points,weights=_quadrature_from_nodes(u,wu,v,wv,phi)
    radial,angular,azimuth=_factors(u,v,phi,labels)
    values=np.empty((len(u),len(v),len(phi),len(labels)))
    for iu in range(len(u)):
        for col,(l,j,m) in enumerate(labels):
            values[iu,:,:,col]=radial[l,j][iu]*angular[j,m][:,None]*azimuth[m][None,:]
    eigenvalues=np.array([l*(l+2) for l,_,_ in labels],float)
    return SurfaceBasis(degree,labels,eigenvalues,points,weights,values.reshape(len(points),len(labels)))

def _build_kick(parent,A):
    basis=build_frame_basis(4,7);x=basis.points
    cc=basis.project(A*x[:,0])
    gradient=np.stack([basis.evaluate(parent.H*(D@cc)) for D in basis.generators],axis=1)
    f=parent.H**1.5*math.sqrt(12)/math.pi*x[:,1]*x[:,3]
    source=-parent.kappa*f[:,None]*gradient;coeff=basis.project(source).T
    solved=solve_longitudinal(basis,coeff,r=parent.H)
    closed,annihilator=inherited_polynomial_inverse(basis,coeff,r=parent.H)
    solved['diagnostics'].update(source_projection_L2_residual=float(np.linalg.norm(basis.evaluate(coeff.T)-source)),
        closed_inverse_relative_difference=float(np.linalg.norm(closed-solved['vector_coefficients'])/np.linalg.norm(closed)),
        source_annihilator_roundoff=annihilator,source='chi=A X1,f=H^1.5 sqrt12/pi X2 X4',
        maximum_source_polynomial_degree=3,maximum_Q_polynomial_degree=6)
    return YorkKick(parent,1.,2,basis,solved['vector_coefficients'],solved['tensor_coefficients'],coeff,solved['diagnostics'])

def _sources(kick,x):
    f=kick.momentum_profile(x);A2=np.empty(len(x))
    for lo in range(0,len(x),8192):A2[lo:lo+8192]=kick.squared_norm(x[lo:lo+8192])
    return f,A2,(A2+kick.parent.kappa*f*f)/kick.parent.H**2

def _radial_preconditioner(basis,neck):
    p=neck.parent;u,w,*_=_nodes(141,1,1);a=24/(p.mass_squared/p.H**2+3)
    psi=neck.values(u);mult=a*(1-u*u)+5*(6+a*p.mass_squared/p.H**2*u*u)*psi**4
    groups={}
    for index,(l,j,m) in enumerate(basis.labels):groups.setdefault((j,m),[]).append(index)
    blocks=[];minimum=math.inf
    for (j,m),indices in groups.items():
        degrees=[basis.labels[i][0] for i in indices]
        R=np.column_stack([_radial_normalization(l,j)*(1-u*u)**(j/2)*eval_gegenbauer(l-j,j+1,u) for l in degrees])
        J=np.diag(8*np.array([l*(l+2) for l in degrees])+6)-R.T@((w*mult)[:,None]*R)
        minimum=min(minimum,float(min(abs(np.linalg.eigvalsh(J)))))
        blocks.append((np.asarray(indices),lu_factor(J)))
    def apply(v):
        out=np.empty_like(v)
        for indices,factor in blocks:out[indices]=lu_solve(factor,v[indices])
        return out
    return LinearOperator((len(basis.labels),)*2,matvec=apply,dtype=float),minimum

def _accounting(p,a,lam,psi,u,w,f,A2):
    volume=float(w@psi**6)/p.H**3
    normal=.5*lam*lam*float(w@(f*f*psi**-6))/p.H**3
    gradient=.5*a/p.kappa*float(w@((1-u*u)*psi**2))/p.H
    scalar_mass=.5*p.mass_squared*a/p.kappa*float(w@(u*u*psi**6))/p.H**3
    vacuum=p.vacuum*volume
    shear_integral=lam*lam*float(w@(A2*psi**-6))/p.H**3
    Pi=lam*f*psi**-6;K2=lam*lam*A2*psi**-12
    return dict(proper_volume=volume,total_proper_matter_energy=normal+gradient+scalar_mass+vacuum,
                normal_Pi_energy=normal,spatial_gradient_energy=gradient,scalar_mass_energy=scalar_mass,vacuum_energy=vacuum,
                normal_Pi_sample_max=float(max(abs(Pi))),normal_Pi_volume_RMS=math.sqrt(2*normal/volume),
                integrated_KijKij=shear_integral,shear_constraint_energy=shear_integral/(2*p.kappa),
                shear_volume_RMS=math.sqrt(shear_integral/volume),shear_sample_max=math.sqrt(float(max(K2))))

def _independent_validation(basis,c,kick,a,lam,orders=(221,43,89)):
    p=kick.parent;u,wu,v,wv,phi=_nodes(*orders)
    york_basis=kick.york_basis
    radial,angular,azimuth=_factors(u,v,phi,basis.labels)
    ang=np.column_stack([(angular[j,m][:,None]*azimuth[m][None,:]).ravel() for l,j,m in basis.labels])
    eig=basis.eigenvalues;divcoeff=tensor_divergence(york_basis,kick.tensor_coefficients,r=p.H)
    sums=dict(R2=0.,physicalH2=0.,V6=0.,physicalM2=0.)
    minima=math.inf;maximum=0.;maxR=0.
    totals=None
    for iu,coord in enumerate(u):
        Y=ang*np.array([radial[l,j][iu] for l,j,m in basis.labels])[None,:]
        psi=Y@c;lap=Y@(-eig*c)
        Z,w=_quadrature_from_nodes(np.array([coord]),np.array([wu[iu]]),v,wv,phi)
        x=original_points(Z);f,A2,Q=_sources(kick,x)
        R=-8*lap+6*psi-a*(1-coord*coord)*psi-(6+a*p.mass_squared/p.H**2*coord*coord)*psi**5-lam*lam*Q*psi**-7
        A=math.sqrt(a/p.kappa)
        gradient=p.H*A*np.column_stack([x[:,3],-x[:,2],x[:,1]])
        div=arbitrary_points_values(x,york_basis.labels)@divcoeff.T
        momentum=div+p.kappa*f[:,None]*gradient
        sums['R2']+=float(w@(R*R));sums['physicalH2']+=float(w@(R*R*psi**-4));sums['V6']+=float(w@psi**6)
        sums['physicalM2']+=float(w@(np.sum(momentum*momentum,axis=1)*psi**-10))
        minima=min(minima,float(min(psi)));maximum=max(maximum,float(max(psi)));maxR=max(maxR,float(max(abs(R))))
        acc=_accounting(p,a,lam,psi,np.full(len(x),coord),w,f,A2)
        if totals is None:totals={key:0. for key in ('proper_volume','total_proper_matter_energy','normal_Pi_energy',
                 'spatial_gradient_energy','scalar_mass_energy','vacuum_energy','integrated_KijKij','shear_constraint_energy')}
        for key in totals:totals[key]+=acc[key]
    return dict(orders=list(orders),quadrature_points=len(u)*len(v)*len(phi),
                full_conformal_residual_RMS=math.sqrt(sums['R2']/(2*math.pi**2)),
                full_conformal_residual_sample_max=maxR,
                physical_Hamiltonian_RMS_over_H_squared=math.sqrt(sums['physicalH2']/sums['V6']),
                physical_momentum_RMS_over_H_squared=abs(lam)*math.sqrt(sums['physicalM2']/sums['V6'])/p.H**2,
                sampled_psi_min=minima,sampled_psi_max=maximum,accounting=totals,
                momentum_source='Independently evaluated divbarA+kappa*f*H*A*(X4,-X3,X2), not subtraction of projected source coefficients.',
                scope='Independent full rational constraint residuals, including omitted angular harmonics; numerical estimates, not interval-certified.')


def solve_counterrotating_initial_data(*,max_j=8,lambda_path=(0.,.25,.5,.75,1.),stage_budget_seconds=120.):
    """One positive-branch amplitude continuation from the existing neck40.

    Fixed defaults are the already checked finite-amplitude experiment.
    Failure raises with its stage; no grid refinement or change of source
    is attempted automatically. The independent full residual is retained.
    """
    if isinstance(max_j,bool) or not isinstance(max_j,int) or max_j<4 or max_j>12 or max_j%2:
        raise ValueError('even max_j in[4,12] required; no automatic refinement')
    path=tuple(float(v) for v in lambda_path);budget=float(stage_budget_seconds)
    if not path or path[0]!=0 or any(not math.isfinite(v) or not 0<=v<=1 for v in path) or any(b<=a for a,b in zip(path,path[1:])):
        raise ValueError('one strictly increasing lambda path starting0 in[0,1] required')
    if not math.isfinite(budget) or budget<=0:raise ValueError('positive finite per-stage budget required')
    started=time.perf_counter();stages=[]
    with threadpool_limits(limits=6):
        neck=solve_necks(40,(1.,))[0];p=neck.parent;z=p.mass_squared/p.H**2;a=24/(z+3);A=math.sqrt(a/p.kappa)
        kick=_build_kick(p,A);orders=(141,3*max_j+1,6*max_j+1)
        basis=_make_basis(40,max_j,orders);x=original_points(basis.points);u=x[:,0]
        f,A2,Q=_sources(kick,x);c=np.zeros(len(basis.labels))
        for l,value in zip(range(0,41,2),neck.coefficients):c[basis.labels.index((l,0,0))]=value
        initial_coefficients=c.copy();linear=8*basis.eigenvalues+6;alpha=a*(1-u*u);B=6+a*z*u*u
        preconditioner,minimum_J=_radial_preconditioner(basis,neck)
        def residual(coef,lam):
            psi=basis.evaluate(coef)
            return linear*coef-basis.project(alpha*psi+B*psi**5+lam*lam*Q*psi**-7)
        for lam in path:
            stage_start=time.perf_counter();calls=0
            for iteration in range(18):
                if time.perf_counter()-stage_start>budget:raise TimeoutError(f'counterrotation lambda{lam}: stage budget exceeded')
                psi=basis.evaluate(c);R=residual(c,lam);norm=float(np.linalg.norm(R))
                if norm<3e-11:break
                mult=alpha+5*B*psi**4-7*lam*lam*Q*psi**-8
                operator=LinearOperator((len(c),)*2,matvec=lambda v:linear*v-basis.project(mult*basis.evaluate(v)),dtype=float)
                def counted(_):
                    nonlocal calls
                    calls+=1
                step,info=gmres(operator,-R,M=preconditioner,rtol=2e-10,atol=2e-14,
                                restart=45,maxiter=10,callback=counted,callback_type='pr_norm')
                if info:raise RuntimeError(f'counterrotation lambda{lam}: Jacobian GMRES failed {info}')
                for halve in range(16):
                    trial=c+step*2.**-halve
                    if min(basis.evaluate(trial))>0 and np.linalg.norm(residual(trial,lam))<norm*(1-1e-4*2.**-halve):
                        c=trial;break
                else:raise RuntimeError(f'counterrotation lambda{lam}: positive branch line search failed')
            else:raise RuntimeError(f'counterrotation lambda{lam}: Newton18-step budget exhausted')
            psi=basis.evaluate(c);legacy=FiniteParentData(p,basis,1.,c,stages)
            stages.append(dict(lambda_value=lam,Newton_steps=iteration,Krylov_iterations=calls,
                projected_residual_norm=float(np.linalg.norm(residual(c,lam))),
                analytic_finite_polynomial_psi_lower_bound=positivity_lower_bound(legacy),
                sampled_psi_min=float(min(psi)),sampled_psi_max=float(max(psi)),
                accounting=_accounting(p,a,lam,psi,u,basis.weights,f,A2),
                elapsed_seconds=time.perf_counter()-stage_start))
        metadata=dict(preparation_id='one-mode-counterrotation-v1',radial_degree=40,angular_max_j=max_j,
            mode_count=len(c),quadrature_orders=list(orders),solve_quadrature_points=len(basis.points),
            lambda_path=list(path),stages=stages,branch_completed=True,
            initial_radial_block_Jacobian_minimum_absolute_eigenvalue=minimum_J,
            York=kick.diagnostics,coefficient_change_from_neck=float(np.linalg.norm(c-initial_coefficients)),
            source_symmetry='psi even under all four coordinate sign flips; l,j,m all even and cosine modes.',
            exact_CKV_compatibility='All six integrals X2X4*(Xa delta_b1-Xb delta_a1) and four integrals X2X4*(delta_a1-XaX1) vanish by odd sphere moments.',
            constraint='-8Delta_gamma psi+[6-a(1-X1²)]psi-[6+a*z*X1²]psi5-lambda²Q psi^-7=0; Q=(|Abar|²+kappa f²)/H².',
            scope='Declared finite initial preparation, not evolution or inherited state selection; residuals are numerical estimates, not a continuum existence certificate.')
        data=CounterrotatingData(p,c,basis.labels,A,path[-1],kick.tensor_coefficients,kick.basis.labels,metadata,{})
        metadata['independent_validation']=_independent_validation(basis,c,data,a,path[-1],(221,5*max_j+3,10*max_j+9))
        metadata['elapsed_seconds']=time.perf_counter()-started
        return data
