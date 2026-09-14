"""One-time physical scattering input for the versioned compact restart.

Only labels/configuration are accepted. Seed covariances are not inputs.
The unchanged exterior massive Dirac owner generates each signed-angular
reflection coefficient; its output is persisted so modal checks do not rerun it.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path

import numpy as np

from .nsc_compact_ctp_neck import _frequency_grid, _massive_reflection
from .nsc_unruh_state import ParentDirac
from .nsc_mode_resolved_cauchy_state import deterministic_npz_bytes


SOURCE_PATHS=(
    'src/recursive_horizons/nsc_compact_restart_scattering.py',
    'src/recursive_horizons/nsc_compact_ctp_neck.py',
    'src/recursive_horizons/nsc_unruh_state.py',
    'src/recursive_horizons/nsc_gauge_source.py',
    'src/recursive_horizons/nsc_mode_resolved_cauchy_state.py',
)


def labels_from_declaration(channels,config):
    rows=[]
    for channel in channels:
        if channel['compact_level']==0:continue
        mass=float(channel['compact_mass']);lam=float(channel['angular_eigenvalue'])
        frequencies,weights,_=_frequency_grid(mass,config['frequency_max'],config['frequency_points'])
        if len(frequencies)!=channel['sample_count']:raise ValueError('declared frequency count changed')
        for sign in ((1,-1) if lam else (1,)):
            for j,(frequency,weight) in enumerate(zip(frequencies,weights)):
                rows.append((channel['index'],sign,channel['sample_offset']+j,mass,sign*lam,float(frequency),float(weight)))
    return np.asarray(rows,float)


def _one(job):
    row,config=job
    b=ParentDirac(config['horizon_rho'],config['surface_gravity'])
    result=_massive_reflection(b,row[5],row[4],row[3],inner_offset=config['horizon_offset'],
        tolerance=config['scattering_tolerance'],outer_floor=config['outer_floor'])
    return result['reflection'],result['transmission'],result['current_defect'],result['function_evaluations']


def generate(root:Path,channels,config,*,workers=4,progress=None):
    labels=labels_from_declaration(channels,config)
    values=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i,value in enumerate(pool.map(_one,((row,config) for row in labels),chunksize=8),1):
            values.append(value)
            if progress is not None and (i%64==0 or i==len(labels)):progress(i,len(labels))
    reflection=np.array([v[0] for v in values],complex)
    transmission=np.array([v[1] for v in values],float)
    metadata={
        'schema':'NSC-COMPACT-RESTART-SCATTERING-v1',
        'source_hashes':{p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'config':config,
        'label_columns':['channel_index','angular_sign','seed_control_index','compact_mass','signed_angular_eigenvalue','frequency_abs','quadrature_weight'],
        'source':'unchanged _massive_reflection on both required angular signs',
        'exterior_start_convention':'radial rho offset; exterior inverse-r_h factor retained',
        'seed_covariance_used':False,
    }
    arrays={
        'labels':labels,'reflection':reflection,'transmission':transmission,
        'primitive_current_defect':np.array([v[2] for v in values],float),
        'function_evaluations':np.array([v[3] for v in values],dtype=np.int64),
        'metadata_json':np.frombuffer(json.dumps(metadata,sort_keys=True,separators=(',',':')).encode(),dtype=np.uint8),
    }
    data=deterministic_npz_bytes(arrays);digest=hashlib.sha256(data).hexdigest()
    path=root/f'results/development/artifacts/nsc-compact-restart-scattering.{digest}.npz'
    if path.exists():
        if path.read_bytes()!=data:raise ValueError('existing scattering artifact differs')
    else:path.write_bytes(data)
    return path


def load_authenticated(root:Path,path:Path,channels,config):
    data=path.read_bytes();digest=hashlib.sha256(data).hexdigest()
    if path.name!=f'nsc-compact-restart-scattering.{digest}.npz':raise ValueError('scattering artifact name/hash mismatch')
    with np.load(path,allow_pickle=False) as archive:arrays={k:archive[k].copy() for k in archive.files}
    metadata=json.loads(arrays['metadata_json'].tobytes())
    if metadata['schema']!='NSC-COMPACT-RESTART-SCATTERING-v1' or metadata['config']!=config:raise ValueError('scattering context differs')
    for name,expected in metadata['source_hashes'].items():
        if hashlib.sha256((root/name).read_bytes()).hexdigest()!=expected:raise ValueError('scattering source changed: '+name)
    if not np.array_equal(arrays['labels'],labels_from_declaration(channels,config)):raise ValueError('scattering labels differ')
    defect=np.max(abs(abs(arrays['reflection'])**2+arrays['transmission']-1))
    if defect>3e-11 or not np.isfinite(arrays['reflection']).all():raise ValueError('invalid physical scattering coefficients')
    return arrays,metadata,{'path':str(path.relative_to(root)),'sha256':digest,'bytes':len(data),'current_normalization_residual':float(defect)}
