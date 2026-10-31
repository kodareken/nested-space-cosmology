"""Full parent geometry at arbitrary S3 points for closed-surface diagnostics.

The provider contract is h, dh_kij=E_k h_ij and K=-Lie_n(h)/2, with
E_i=r*A_i*X and [E_i,E_j]=2*r*epsilon_ijk E_k. No spherical average,
waveform differentiation or new GH evolution is performed.
"""
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess

import numpy as np

from .nsc_atlas_finite_parent_data import source_profile
from .nsc_atlas_gh_constraints import adm_constraints
from .nsc_atlas_gh_history import DenseGHHistory
from .nsc_atlas_parent_gh import PAIRS,ParentGH
from .nsc_atlas_s3_frame import arbitrary_points_values
from .nsc_atlas_s3_hopf import build_hopf_frame_basis

ROOT=Path(__file__).resolve().parents[2]
EPS=np.zeros((3,3,3))
for _i,_j,_k in ((0,1,2),(1,2,0),(2,0,1)):
    EPS[_i,_j,_k]=1.;EPS[_j,_i,_k]=-1.


def _points(value):
    X=np.asarray(value)
    if X.dtype.kind not in 'iuf' or X.ndim!=2 or X.shape[1]!=4 or not np.all(np.isfinite(X)):
        raise ValueError('surface points require finite real shape(N,4)')
    X=np.ascontiguousarray(X,dtype=float)
    if np.any(abs(np.sum(X*X,axis=1)-1)>2e-10):
        raise ValueError('surface points must lie on unit S3')
    return X


def _chunk(value):
    if isinstance(value,(bool,np.bool_)) or not isinstance(value,(int,np.integer)) or value<1:
        raise ValueError('positive integer chunk size required')
    return int(value)


def _sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _record(value):
    if isinstance(value,dict):return dict(value)
    return json.loads(Path(value).read_bytes())


def _verify_record(record):
    """Authenticate producer objects and require strict current-source identity."""
    if record.get('schema') not in ('nsc-atlas-parent-gh-v2','nsc-atlas-parent-gh-v3'):
        raise ValueError('recognized GH state schema required')
    commit=record.get('producer_commit','')
    if not re.fullmatch('[0-9a-f]{40}',commit) or not record.get('source_bound_to_commit') or not record.get('producer_source_unchanged_during_run'):
        raise ValueError('source-bound unchanged producer record required')
    sources=record.get('source_sha256')
    required={'src/recursive_horizons/nsc_atlas_parent_gh.py',
        'src/recursive_horizons/nsc_atlas_s3_frame.py',
        'src/recursive_horizons/nsc_atlas_parent_circulation.py',
        'src/recursive_horizons/nsc_atlas_parent_scale_selector.py'}
    if record.get('spatial_transform')=='hopf':required.add('src/recursive_horizons/nsc_atlas_s3_hopf.py')
    if not isinstance(sources,dict) or not required.issubset(sources):
        raise ValueError('producer record lacks required decoder source pins')
    for name,expected in sources.items():
        relative=Path(name)
        if relative.is_absolute() or '..' in relative.parts or not re.fullmatch('[0-9a-f]{64}',str(expected)):
            raise ValueError('invalid producer source pin')
        stored=subprocess.run(['git','show',commit+':'+relative.as_posix()],cwd=ROOT,capture_output=True)
        if stored.returncode or hashlib.sha256(stored.stdout).hexdigest()!=expected:
            raise ValueError('producer Git source pin failed: '+name)
        if not (ROOT/relative).is_file() or _sha(ROOT/relative)!=expected:
            raise ValueError('current decoder source differs from producer: '+name)
    return commit


def _cache_system(metadata,record):
    if metadata.get('source_sha256')!=record['source_sha256']:
        raise ValueError('cache source pins differ from producer record')
    for key,expected in (('degree',record.get('degree')),('epsilon',record.get('epsilon')),
            ('gauge',record.get('prescribed_GH_gauge')),('backend',record.get('spatial_transform'))):
        if metadata.get(key)!=expected:raise ValueError('cache model metadata differs: '+key)
    system=ParentGH(metadata['degree'],metadata['epsilon'],gauge=metadata['gauge'],
        backend=metadata['backend'],quadrature_order=metadata.get('quadrature_order'),
        kick_amplitude=metadata.get('kick_amplitude',0.),kick_source=metadata.get('kick_source',1))
    model=record.get('model_inputs',{})
    for name,actual in (('H_reference',system.r),('vacuum',system.parent.vacuum),
            ('scalar_mass_squared',system.parent.mass_squared),('volume_kinetic',system.K),('sstar',system.sstar)):
        if name not in model or not math.isclose(float(model[name]),actual,rel_tol=2e-13,abs_tol=2e-14):
            raise ValueError('current model constants differ from producer: '+name)
    for key,default in (('kick_amplitude',0.),('kick_source',1)):
        if metadata.get(key,default)!=model.get(key,default):
            raise ValueError('cache kick metadata differs: '+key)
    return system


def _first_geometry(g,v,dg,r):
    """ADM normal and K from the spacetime connection difference, first jets only."""
    h=g[:,1:,1:]
    if np.any(np.linalg.eigvalsh(h)<=0):raise ValueError('spatial metric is not positive definite')
    inverse=np.linalg.inv(g)
    if not np.all(np.isfinite(inverse)) or np.any(inverse[:,0,0]>=0):
        raise ValueError('slice does not have positive real lapse')
    lapse=1/np.sqrt(-inverse[:,0,0]);shift=np.einsum('nij,nj->ni',np.linalg.inv(h),g[:,0,1:])
    covariant=np.empty((len(g),4,4,4));covariant[:,0]=v
    for i in range(3):
        connection=np.zeros((4,4))
        for a in range(3):
            for b in range(3):connection[a+1,b+1]=r*EPS[i,b,a]
        covariant[:,i+1]=dg[:,i]-connection.T@g-g@connection
    K=np.empty((len(g),3,3))
    for i in range(3):
        for j in range(3):
            lower=covariant[:,i+1,:,j+1]+covariant[:,j+1,:,i+1]-covariant[:,:,i+1,j+1]
            K[:,i,j]=-.5*lapse*np.einsum('na,na->n',inverse[:,0,:],lower)
    if not np.all(np.isfinite(K)):raise ValueError('extrinsic curvature is not finite')
    L=np.linalg.cholesky(h);left=np.linalg.solve(L,K)
    orthonormal=np.linalg.solve(L,left.swapaxes(1,2)).swapaxes(1,2)
    eigenvalues=np.linalg.eigvalsh(orthonormal)
    return dict(h=h,dh=dg[:,:,1:,1:],K=K,lapse=lapse,shift=shift,
        K_eigenvalues=eigenvalues,future_tangential_capacity=eigenvalues[:,1]+eigenvalues[:,2])


class GHSurfaceGeometry:
    """Callable exact spectral first jets of one accepted GH state.

    A direct constructor accepts caller-owned state data explicitly. Use
    the strict cache factories to authenticate archived inputs. The state
    is copied and never changed by surface iteration.
    """
    def __init__(self,system,state,*,coordinate_t=None,chunk=1024,provenance=None):
        if not isinstance(system,ParentGH):raise ValueError('existing ParentGH system required')
        y=np.asarray(state)
        if y.dtype.kind not in 'iuf' or y.shape!=(system.dof+10,) or not np.all(np.isfinite(y)):
            raise ValueError('finite full GH state of the declared system required')
        self.system=system;self.r=system.r;self.chunk=_chunk(chunk)
        self.state=np.array(y,dtype=float,copy=True);self.state.setflags(write=False)
        self.coefficients=system.unpack(self.state);self.coefficients.setflags(write=False)
        self.degree=system.basis.degree;self.coordinate_t=coordinate_t
        self.proper_duration=float(self.state[system.dof+8])
        self.provenance=dict(provenance or {'source_binding':'explicit_caller_owned_state'})
        self.provenance.update(adapter_sha256=_sha(__file__),state_sha256=hashlib.sha256(self.state.tobytes()).hexdigest())
        self._first=np.array([self.r*(D@self.coefficients.T).T for D in system.basis.generators])

    @classmethod
    def from_state_cache(cls,path,producer_record,*,expected_sha256=None,snapshot=-1,chunk=1024):
        record=_record(producer_record);commit=_verify_record(record);path=Path(path)
        expected=record.get('state_cache_sha256',expected_sha256)
        if expected is None:
            raise ValueError('sparse state cache lacks producer array hash; an explicitly admitted expected_sha256 is required')
        if not re.fullmatch('[0-9a-f]{64}',str(expected)) or _sha(path)!=expected:
            raise ValueError('state cache SHA256 pin failed')
        with np.load(path,allow_pickle=False) as data:
            times=data['times'].copy();states=data['states'].copy();metadata=json.loads(str(data['metadata']))
        system=_cache_system(metadata,record)
        points=record.get('snapshots',[])
        if times.ndim!=1 or states.shape!=(len(times),system.dof+10) or len(points)!=len(times) or not len(times):
            raise ValueError('state-cache snapshots do not match producer shape')
        if not np.all(np.isfinite(times)) or not np.all(np.isfinite(states)) or np.any(np.diff(times)<0):
            raise ValueError('state cache has nonfinite or unordered samples')
        for time,state,point in zip(times,states,points):
            if not math.isclose(float(time),point['coordinate_t'],rel_tol=0,abs_tol=2e-12) or not math.isclose(
                    float(state[system.dof+8]),point['observer_proper_duration'],rel_tol=0,abs_tol=2e-10):
                raise ValueError('state-cache event does not match producer record')
        if isinstance(snapshot,(bool,np.bool_)) or not isinstance(snapshot,(int,np.integer)) or not -len(times)<=snapshot<len(times):
            raise ValueError('snapshot index is outside the admitted cache')
        return cls(system,states[snapshot],coordinate_t=float(times[snapshot]),chunk=chunk,
            provenance=dict(source_binding='strict_producer_and_current_source_identity',producer_commit=commit,
                producer_status=record.get('status'),cache_sha256=expected,
                cache_pin_origin='producer_record' if 'state_cache_sha256' in record else 'explicit_caller_admission',snapshot=int(snapshot)))

    @classmethod
    def from_history(cls,path,producer_record,*,coordinate_t,chunk=1024,retain_history=False):
        record=_record(producer_record);commit=_verify_record(record);path=Path(path)
        expected=record.get('dense_history',{}).get('sha256')
        if expected is None or _sha(path)!=expected:raise ValueError('dense history SHA256 pin failed')
        history=DenseGHHistory.load(path);system=_cache_system(history.metadata,record)
        if history.metadata.get('producer_commit')!=commit:raise ValueError('dense history producer mismatch')
        t=float(coordinate_t);endpoint=float(history.metadata['physical_coordinate_endpoint'])
        last=record.get('snapshots',[{}])[-1]
        if not math.isclose(endpoint,float(last.get('coordinate_t',math.nan)),rel_tol=0,abs_tol=2e-12) or not math.isclose(
                float(history.metadata.get('proper_endpoint',math.nan)),float(last.get('observer_proper_duration',math.nan)),rel_tol=0,abs_tol=2e-10):
            raise ValueError('dense history physical endpoint differs from producer record')
        if not math.isfinite(t) or not history.bounds[0][0]<=t<=endpoint:
            raise ValueError('surface event lies outside the declared physical history')
        provider=cls(system,history.evaluate(t),coordinate_t=t,chunk=chunk,
            provenance=dict(source_binding='strict_producer_and_current_source_identity',producer_commit=commit,
                producer_status=record.get('status'),cache_sha256=expected,cache_pin_origin='producer_record'))
        if retain_history:
            # Passive event searches own this same authenticated archive;
            # no loader patch, repeated decompression or new evolution.
            provider.retained_history=history
        return provider

    def _jets(self,points,*,second=False):
        Y=arbitrary_points_values(points,self.system.basis.labels)
        values=Y@self.coefficients.T
        selected=np.arange(24) if second else np.r_[np.arange(10),[20,21]]
        first=np.array([Y@coeff[selected].T for coeff in self._first]).transpose(1,0,2)
        count=len(points);g=np.empty((count,4,4));v=np.empty_like(g)
        dg=np.empty((count,3,4,4));dv=np.empty_like(dg) if second else None
        for k,(a,b) in enumerate(PAIRS):
            g[:,a,b]=g[:,b,a]=values[:,k];v[:,a,b]=v[:,b,a]=values[:,k+10]
            dg[:,:,a,b]=dg[:,:,b,a]=first[:,:,k]
            if second:dv[:,:,a,b]=dv[:,:,b,a]=first[:,:,k+10]
        result=dict(g=g,v=v,dg=dg,q=values[:,20:22],p=values[:,22:24],
            dq=first[:,:,20:22] if second else first[:,:,10:12])
        if second:
            result['dv']=dv
            selected=np.r_[np.arange(10),[20,21]];C=self.coefficients[selected]
            d2=np.empty((count,3,3,12))
            for i,Di in enumerate(self.system.basis.generators):
                for j,Dj in enumerate(self.system.basis.generators):
                    d2[:,i,j]=Y@(self.r**2*(Di@(Dj@C.T)))
            d2g=np.empty((count,3,3,4,4))
            for k,(a,b) in enumerate(PAIRS):d2g[:,:,:,a,b]=d2g[:,:,:,b,a]=d2[:,:,:,k]
            result['d2g']=d2g
        return result

    def __call__(self,points):
        X=_points(points);out={name:np.empty(shape) for name,shape in
            (('h',(len(X),3,3)),('dh',(len(X),3,3,3)),('K',(len(X),3,3)),
             ('lapse',(len(X),)),('shift',(len(X),3)),('q',(len(X),2)),
             ('Pi',(len(X),2)),('dq',(len(X),3,2)),('W',(len(X),)),
             ('K_eigenvalues',(len(X),3)),('future_tangential_capacity',(len(X),)))}
        for lo in range(0,len(X),self.chunk):
            hi=min(lo+self.chunk,len(X));jet=self._jets(X[lo:hi]);geometry=_first_geometry(jet['g'],jet['v'],jet['dg'],self.r)
            geometry.update(q=jet['q'],dq=jet['dq'],W=self.system.potential(jet['q']),
                Pi=(jet['p']-np.einsum('ni,nia->na',geometry['shift'],jet['dq']))/geometry['lapse'][:,None])
            for name in out:out[name][lo:hi]=geometry[name]
        return out

    def adm_sample(self,points):
        """Independent ADM constraints at the supplied surface-neighborhood nodes."""
        X=_points(points)
        if not len(X):raise ValueError('nonempty ADM sample required')
        outputs=[]
        for lo in range(0,len(X),self.chunk):
            jet=self._jets(X[lo:lo+self.chunk],second=True)
            outputs.append(adm_constraints(**jet,r=self.r,kappa=self.system.parent.kappa,
                potential=self.system.potential(jet['q'])))
        return {name:np.concatenate([part[name] for part in outputs],axis=0) for name in outputs[0]}


class TimeSymmetricParentGeometry:
    """Exact K=0 provider from the existing finite conformal-data interface."""
    def __init__(self,data,*,chunk=1024):
        self.data=data;self.r=data.parent.H;self.degree=data.basis.degree;self.chunk=_chunk(chunk)
        self.basis=build_hopf_frame_basis(self.degree)
        samples=np.empty(len(self.basis.points))
        for lo in range(0,len(samples),self.chunk):samples[lo:lo+self.chunk]=data.values(self.basis.points[lo:lo+self.chunk])
        self.coefficients=self.basis.project(samples)
        residual=self.basis.evaluate(self.coefficients)-samples
        relative=np.linalg.norm(residual)/max(np.linalg.norm(samples),np.finfo(float).tiny)
        if relative>2e-11:raise ValueError('conformal field is not represented by its declared finite degree')
        self._gradient=np.array([self.r*(D@self.coefficients) for D in self.basis.generators])
        self._laplacian=-self.r**2*np.array([label[0]*(label[0]+2) for label in self.basis.labels])*self.coefficients
        f,_=source_profile(self.basis.points,data.parent)
        self._chi_coefficients=self.basis.project(data.epsilon*f)
        self._chi_gradient=np.array([self.r*(D@self._chi_coefficients) for D in self.basis.generators])
        self.proper_duration=0.;self.coordinate_t=0.
        self.provenance=dict(source_binding='explicit_current_conformal_data',
            conversion_relative_residual=float(relative),adapter_sha256=_sha(__file__))

    def __call__(self,points):
        X=_points(points);h=np.zeros((len(X),3,3));dh=np.zeros((len(X),3,3,3))
        q=np.zeros((len(X),2));dq=np.zeros((len(X),3,2))
        for lo in range(0,len(X),self.chunk):
            hi=min(lo+self.chunk,len(X));Y=arbitrary_points_values(X[lo:hi],self.basis.labels)
            psi=Y@self.coefficients;gradient=Y@self._gradient.T
            if np.any(psi<=0):raise ValueError('conformal metric is not positive')
            for i in range(3):h[lo:hi,i,i]=psi**4;dh[lo:hi,:,i,i]=4*psi[:,None]**3*gradient
            q[lo:hi,1]=Y@self._chi_coefficients;dq[lo:hi,:,1]=Y@self._chi_gradient.T
        return dict(h=h,dh=dh,K=np.zeros_like(h),lapse=np.ones(len(X)),shift=np.zeros((len(X),3)),
            q=q,dq=dq,Pi=np.zeros_like(q),W=self.data.parent.vacuum+self.data.parent.mass_squared*q[:,1]**2/2,
            K_eigenvalues=np.zeros((len(X),3)),future_tangential_capacity=np.zeros(len(X)))

    def adm_sample(self,points):
        X=_points(points)
        if not len(X):raise ValueError('nonempty ADM sample required')
        R=np.empty(len(X));rho=np.empty(len(X));eps=self.data.epsilon
        for lo in range(0,len(X),self.chunk):
            hi=min(lo+self.chunk,len(X));Y=arbitrary_points_values(X[lo:hi],self.basis.labels);psi=Y@self.coefficients
            if np.any(psi<=0):raise ValueError('conformal metric is not positive')
            R[lo:hi]=psi**-5*(-8*(Y@self._laplacian)+6*self.r**2*psi)
            f,D=source_profile(X[lo:hi],self.data.parent)
            rho[lo:hi]=self.data.parent.vacuum+eps**2*(self.r**2*psi**-4*D+self.data.parent.mass_squared*f*f)/2
        return dict(R3=R,rho=rho,H=R-2*self.data.parent.kappa*rho,M=np.zeros((len(X),3)),Kij=np.zeros((len(X),3,3)))


def compare_surface_neighborhood(providers,points,*,weights=None):
    """Physical sampled tensor gaps on common points; no horizon error bound."""
    if len(providers)<2:raise ValueError('at least two explicit geometry providers required')
    X=_points(points)
    if not len(X):raise ValueError('nonempty neighborhood required')
    w=np.ones(len(X)) if weights is None else np.asarray(weights,dtype=float)
    if w.shape!=(len(X),) or not np.all(np.isfinite(w)) or np.any(w<0) or not w.sum()>0:
        raise ValueError('nonnegative finite sampling weights with positive sum required')
    w=w/w.sum();geometry=[provider(X) for provider in providers];rows=[]
    for index,(old,new) in enumerate(zip(geometry,geometry[1:])):
        r=float(providers[index+1].r)
        if not math.isclose(float(providers[index].r),r,rel_tol=2e-13,abs_tol=0):
            raise ValueError('geometry comparisons require a common reference ruler')
        inverse=np.linalg.inv(new['h']);delta_h=new['h']-old['h'];delta_K=new['K']-old['K']
        strain=np.einsum('nik,njl,nij,nkl->n',inverse,inverse,delta_h,delta_h)
        curvature=np.einsum('nik,njl,nij,nkl->n',inverse,inverse,delta_K,delta_K)
        old_norm=np.einsum('nik,njl,nij,nkl->n',inverse,inverse,new['K'],new['K'])
        if min(float(np.min(strain)),float(np.min(curvature)))<-1e-12:
            raise ValueError('sampled physical tensor norm became negative')
        metric_rms=math.sqrt(float(w@np.maximum(strain,0)));K_rms=math.sqrt(float(w@np.maximum(curvature,0)))
        Kscale=math.sqrt(float(w@np.maximum(old_norm,0)))
        gradient=new['dh']-old['dh']
        for k in range(3):
            connection=r*EPS[k].T
            gradient[:,k]-=connection.T@delta_h+delta_h@connection
        gradnorm=np.einsum('nka,nib,njc,nkij,nabc->n',inverse,inverse,inverse,gradient,gradient,optimize=True)
        if float(np.min(gradnorm))<-1e-12:raise ValueError('sampled gradient norm became negative')
        rows.append(dict(from_provider=index,to_provider=index+1,
            metric_relative_strain_RMS=metric_rms,metric_relative_strain_max=math.sqrt(max(0.,float(max(strain)))),
            metric_reference_covariant_gradient_RMS_over_reference_scale=math.sqrt(float(w@np.maximum(gradnorm,0)))/r,
            extrinsic_curvature_physical_RMS=K_rms,
            extrinsic_curvature_RMS_over_reference_scale=K_rms/r,
            extrinsic_curvature_relative_RMS=K_rms/Kscale if Kscale else None,
            future_tangential_capacity_max_abs_gap=float(np.max(abs(new['future_tangential_capacity']-old['future_tangential_capacity'])))))
    return dict(points=len(X),comparisons=rows,
        scope='Weighted finite samples in a common declared surface neighborhood. Same-event/protocol identification is owned by the caller; no certified global bound or horizon classification follows.')
