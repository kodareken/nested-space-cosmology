"""Explicitly versioned historical and matched compact seed generations.

Both paths start from the owned horizon/incoming state. The same inherited
6000-step commutator-free SU(2) evolution is shared between initial frames.
No seed covariance is used to construct either result.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from math import exp,log,pi,sqrt

import numpy as np

from .nsc_compact_ctp_neck import _su2_step
from .nsc_paired_horizon_preparation import PairedHorizonSeedMap,source_covariance
from .nsc_transmitting_dirac_domain import S3


HISTORICAL='historical_radial_r_h'
MATCHED='matched_delta_q'
VERSIONS={HISTORICAL:'compact-seed-v1-historical',MATCHED:'compact-seed-v2-matched-delta-q'}


def require_convention(value):
    if value not in VERSIONS:raise ValueError('explicit historical_radial_r_h or matched_delta_q convention required')
    return value


def inherited_propagator(labels,config,*,progress=None):
    mass=labels[:,3,None];angular=labels[:,4,None];frequency=labels[:,5,None]
    count=len(labels)
    owner=PairedHorizonSeedMap(config['horizon_rho'],config['surface_gravity'],config['omega'],
        config['horizon_offset'],config['scattering_tolerance'],config['outer_floor'])
    b=owner.background;begin=log(config['horizon_offset']);end=log(b.horizon_q-pi/2)
    steps=int(config['evolution_steps']);step=(end-begin)/steps
    aa=(3-2*sqrt(3))/12;bb=(3+2*sqrt(3))/12
    c1=.5-sqrt(3)/6;c2=.5+sqrt(3)/6
    upper=np.broadcast_to(np.array([1.,0.],complex),(count,2)).copy()
    lower=np.broadcast_to(np.array([0.,1.],complex),(count,2)).copy()
    def h(value):
        delta=exp(value);W=b.interior_W(delta);r=1/np.sin(b.horizon_q-delta)
        return -mass*r*delta/sqrt(W),angular*delta/sqrt(W),-frequency*delta/W
    for index in range(steps):
        y=begin+index*step
        x1,y1,z1=h(y+c1*step);x2,y2,z2=h(y+c2*step)
        for hx,hy,hz in ((bb*x1+aa*x2,bb*y1+aa*y2,bb*z1+aa*z2),
                         (aa*x1+bb*x2,aa*y1+bb*y2,aa*z1+bb*z2)):
            upper,lower=_su2_step(upper,lower,hx,hy,hz,step)
        if progress is not None and (index+1)%1000==0:progress(index+1,steps)
    return np.stack((upper,lower),axis=1)


def assemble_versions(scattering,config,*,progress=None):
    labels=scattering['labels'];count=len(labels)
    u=inherited_propagator(labels,config,progress=progress)
    owner=PairedHorizonSeedMap(config['horizon_rho'],config['surface_gravity'],config['omega'],
        config['horizon_offset'],config['scattering_tolerance'],config['outer_floor'])
    source=np.array([source_covariance(row[5],config['surface_gravity'],config['omega'],row[3]) for row in labels])
    negative_source=np.array([source_covariance(-row[5],config['surface_gravity'],config['omega'],row[3]) for row in labels])
    sewing=np.zeros((count,2,3),complex)
    sewing[:,0,1]=1.;sewing[:,1,0]=scattering['reflection'];sewing[:,1,2]=np.sqrt(scattering['transmission'])
    lookup={(int(row[0]),int(row[1]),int(row[2])):i for i,row in enumerate(labels)}
    partner=np.array([lookup[(int(row[0]),-int(row[1]) if row[4]!=0 else 1,int(row[2]))] for row in labels],dtype=np.int64)
    arrays={'labels':labels,'unitary':u,'source_covariance_positive':source,'source_covariance_negative':negative_source,
            'signed_partner_index':partner}
    for convention in (HISTORICAL,MATCHED):
        frames=np.array([owner.working_frame(log(config['horizon_offset']),row[5],row[3],row[4],
                         legacy_compact_distance=convention==HISTORICAL) for row in labels])
        modes=u@frames@sewing
        covariance=modes@source@modes.swapaxes(-1,-2).conj()
        negative_modes=S3@modes[partner].conj()
        negative=negative_modes@negative_source@negative_modes.swapaxes(-1,-2).conj()
        arrays['mode_map_'+convention]=modes
        arrays['covariance_'+convention]=covariance
        arrays['negative_mode_map_'+convention]=negative_modes
        arrays['negative_covariance_'+convention]=negative
    return arrays


@dataclass(frozen=True)
class CompactSeedGeneration:
    """A seed-frequency checkpoint; never implicitly a spatial PG state."""
    start_convention: str
    labels: np.ndarray
    covariance: np.ndarray
    negative_covariance: np.ndarray
    mode_map: np.ndarray

    @classmethod
    def select(cls,arrays,*,start_convention):
        name=require_convention(start_convention)
        return cls(name,arrays['labels'],arrays['covariance_'+name],arrays['negative_covariance_'+name],arrays['mode_map_'+name])

    @classmethod
    def load(cls,path,*,start_convention,expected_sha256):
        require_convention(start_convention)
        if hashlib.sha256(path.read_bytes()).hexdigest()!=expected_sha256:
            raise ValueError('versioned seed artifact digest differs')
        with np.load(path,allow_pickle=False) as archive:
            metadata=json.loads(archive['metadata_json'].tobytes())
            if metadata['schema']!='NSC-COMPACT-SEED-GENERATIONS-v2' or metadata['start_conventions']!=VERSIONS:
                raise ValueError('explicit convention/version declaration differs')
            arrays={k:archive[k].copy() for k in archive.files}
        return cls.select(arrays,start_convention=start_convention)

    @property
    def version(self):return VERSIONS[self.start_convention]

    def as_spatial_PG_covariance(self):
        raise ValueError('a versioned seed-frequency checkpoint still requires the global PG Cauchy mode/current normalization')
