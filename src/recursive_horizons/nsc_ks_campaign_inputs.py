"""Explicit raw campaign availability and finite linear-algebra controls.

Descriptors are metadata, never substitutes for absent numerical arrays. These
checks neither fetch nor generate data and never authenticate numeric data from
Git. Finite linear-algebra residual allowances describe rounding of the stated
operands; they are not physical response errors or local-gate certificates.
"""
from hashlib import sha256
import json
from pathlib import Path
import re

import numpy as np

from .nsc_ks_saved_value_arithmetic import _gamma, _up
from .nsc_frozen_sources import authenticate_source, is_implementation_source

_STAGE0_SOURCE_CONTEXT='329ef809e466a86ed4a290ee32c19485f85add9d'
_STAGE0_ORIGINAL_PINS={
    'src/recursive_horizons/nsc_ks_difference_envelope.py':(
        'e62d28f3230ae6ddb2712df71eb7c61a8ed65ad6',
        '77c89c8e720ca97336ccdee441a645b5c73472cb6f15b09e8120f68ab53db196'),
    'src/recursive_horizons/nsc_ks_energy_propagator.py':(
        'c7dc0f98cacb28afd030b089bdb9e1369288120c',
        '5f75023efc8c2b637fc4983b414cab50f15f90f8abc1a24dbc1d23ecef9ecf94'),
    'src/recursive_horizons/nsc_local_incoming_family.py':(
        '09f3faa364ca88b8ff482c413022ce38f40bef0d',
        '9d2b9fc57e2057d6c34d8943384428341b314370021bcd4156af88320c3de04b'),
}
_STAGE0_UNAVAILABLE_SOURCE=(
    'src/recursive_horizons/nsc_ks_history_evaluator.py',
    '86ed312db133f201a3cce11cacc96b62dc89d6c3fe434bdac0f77ca54d023cd6')


def stage0_capture_source_provenance(root, hashes):
    """Source-only report for one declared original metadata capture.

    Three inherited byte identities have explicit existing Merkle pins. One
    declared inherited source remains unavailable; the report keeps that exact
    SHA and does not certify the captured old mathematics or numerical replay.
    Strict numerical replay uses require_original_campaign_sources instead.
    """
    producer='scripts/derive_nsc_ks_coupled_newton_n32_truncated_iterate.py'
    missing_name,missing_sha=_STAGE0_UNAVAILABLE_SOURCE
    if not isinstance(hashes,dict) or producer not in hashes or missing_name not in hashes:
        raise ValueError('captured proposal source scope changed')
    bindings=[]
    for name,expected in hashes.items():
        _sha(expected,name)
        if name==missing_name:
            if expected!=missing_sha:
                raise ValueError('declared unavailable inherited source SHA256 changed')
            bindings.append({'path':name,'expected_sha256':expected,'authenticated':False,
                'reason':'declared inherited history-evaluator bytes unavailable in bounded preserved source history'})
            continue
        context=_STAGE0_SOURCE_CONTEXT
        if name in _STAGE0_ORIGINAL_PINS:
            context,pinned_sha=_STAGE0_ORIGINAL_PINS[name]
            if expected!=pinned_sha:raise ValueError('declared original source pin changed: '+name)
        authenticate_source(root,name,expected,context=context)
        bindings.append({'path':name,'expected_sha256':expected,'authenticated':True,'context':context})
    return {'scope':'source byte provenance of captured metadata only; no old numerical replay',
            'creating_CLI_authenticated':True,'inherited_math_authenticated':False,
            'bindings':bindings,'physical_local_gate':'OPEN'}


def require_original_campaign_sources(root, hashes, *, source_context=None):
    """Strict source/data binding check; no unavailable source is accepted."""
    missing_name,missing_sha=_STAGE0_UNAVAILABLE_SOURCE
    for name,expected in hashes.items():
        if name==missing_name and expected==missing_sha:
            raise ValueError('historical numerical replay original inherited source unavailable: '+name)
        if is_implementation_source(name):
            context=source_context
            if source_context is not None and name in _STAGE0_ORIGINAL_PINS:
                pin,pinned_sha=_STAGE0_ORIGINAL_PINS[name]
                if expected!=pinned_sha:raise ValueError('declared original source pin changed: '+name)
                context=pin
            authenticate_source(root,name,expected,context=context)
        else:
            _sha(expected,name)
            path=Path(root)/name
            if not path.is_file() or sha256(path.read_bytes()).hexdigest()!=expected:
                raise ValueError('current local campaign numeric binding changed: '+name)


class ExternalRawInputsUnavailable(FileNotFoundError):
    """An explicitly declared campaign's strict local raw input is absent."""


def _sha(value, label):
    if not isinstance(value,str) or re.fullmatch(r'[0-9a-f]{64}',value) is None:
        raise ValueError('invalid campaign descriptor SHA256: '+label)
    return value


def descriptor_inventory(root, directory, *, expected_count=60, expected_digest=None):
    """Validate a local descriptor census, without claiming raw array integrity."""
    root=Path(root).resolve();directory=Path(directory).resolve()
    if root not in directory.parents:
        raise ValueError('campaign directory escapes repository')
    if type(expected_count) is not int or expected_count <= 0:
        raise ValueError('positive exact descriptor count required')
    files=sorted(directory.glob('family-*.json'))
    if len(files)!=expected_count:
        raise ValueError('campaign descriptor coverage changed: '+directory.name)
    census=sha256();rows=[];keys=set()
    for file in files:
        raw=file.read_bytes()
        census.update(file.name.encode('utf8')+b'\0'+sha256(raw).digest())
        record=json.loads(raw)
        if not isinstance(record,dict):raise ValueError('campaign descriptor must be a mapping')
        key=record.get('positive_family')
        if (not isinstance(key,list) or len(key)!=2 or any(type(x) is not int for x in key)
                or not 0<=key[0]<=99 or key[1] not in (-1,1) or tuple(key) in keys):
            raise ValueError('campaign positive-family coverage changed')
        keys.add(tuple(key))
        if file.name!=f'family-{key[0]:02d}-{key[1]:+d}.json':
            raise ValueError('campaign family filename identity changed')
        _sha(record.get('profile_identity'),'profile identity')
        hashes=record.get('source_hashes')
        if not isinstance(hashes,dict) or not hashes:
            raise ValueError('campaign descriptor source metadata missing')
        for name,digest in hashes.items():
            if not isinstance(name,str):raise ValueError('invalid descriptor source path')
            _sha(digest,name)
        payload=record.get('payload')
        if not isinstance(payload,dict) or set(payload)!= {'path','bytes','sha256'}:
            raise ValueError('campaign descriptor payload keys changed')
        expected=file.with_suffix('.npz').relative_to(root).as_posix()
        if payload['path']!=expected:
            raise ValueError('campaign descriptor payload path changed')
        if type(payload['bytes']) is not int or payload['bytes']<=0:
            raise ValueError('campaign descriptor payload byte count changed')
        _sha(payload['sha256'],'payload')
        rows.append({'descriptor_path':file.relative_to(root).as_posix(),
                     'record':record,'payload':payload})
    if expected_digest is not None:
        if census.hexdigest()!=_sha(expected_digest,'descriptor census'):
            raise ValueError('campaign descriptor census differs from original capture')
    return {'rows':rows,'descriptor_sha256':census.hexdigest(),
            'scope':'descriptor byte integrity and coverage only; raw arrays not authenticated',
            'raw_array_integrity_claimed':False}


def raw_input_availability(root, directory, *, expected_count=60):
    """Check present bytes strictly; report missing original arrays explicitly."""
    root=Path(root).resolve()
    inventory=descriptor_inventory(root,directory,expected_count=expected_count)
    missing=[];present=0
    for row in inventory['rows']:
        payload=row['payload'];path=root/payload['path']
        if not path.is_file():
            missing.append(payload['path']);continue
        if path.stat().st_size!=payload['bytes']:
            raise ValueError('local campaign raw payload byte count changed: '+payload['path'])
        if sha256(path.read_bytes()).hexdigest()!=payload['sha256']:
            raise ValueError('local campaign raw payload SHA256 changed: '+payload['path'])
        present+=1
    return {'missing':missing,'present_authenticated':present,
            'required':len(inventory['rows']),'all_available':not missing,
            'scope':'strict local raw bytes only; no numerical replay performed'}


def require_campaign_raw_inputs(root, directory, *, expected_count=60):
    result=raw_input_availability(root,directory,expected_count=expected_count)
    if result['missing']:
        raise ExternalRawInputsUnavailable(
            'external raw campaign inputs unavailable: '+Path(directory).name+
            f" ({len(result['missing'])}/{result['required']} NPZs absent; first: "+
            result['missing'][0]+'); supply original bytes with recorded SHA256; no fetch or regeneration')
    return result


def validate_truncated_linear_step(gradient, jacobian, delta, cutoff, kept_modes):
    """Check a current truncated-SVD step's finite linear system independently.

    The SVD's matrix reconstruction and orthogonality are checked as numerical
    backward residuals. The resulting step must solve the retained left-space
    projection and lie in the retained right-space. Allowances charge actual
    uncancelled matrix products, the observed SVD backward residual, and the
    original division/multiplication recipe. No old bit pattern is required.
    """
    arrays=[np.asarray(x) for x in (gradient,jacobian,delta)]
    raw_cutoff=np.asarray(cutoff)
    if (any(a.dtype.kind not in 'fiu' or not np.isfinite(a).all() for a in arrays)
            or raw_cutoff.shape!=() or raw_cutoff.dtype.kind not in 'fiu'
            or not np.isfinite(raw_cutoff) or raw_cutoff<=1):
        raise ValueError('invalid truncated-step real numeric input')
    g,J,d=(a.astype(float,copy=False) for a in arrays)
    if (g.ndim!=2 or g.shape[0]<=0 or g.shape[1]!=2 or J.ndim!=3 or J.shape[1:]!=g.shape
            or J.shape[0]<=0 or d.shape!=(J.shape[0],)
            or type(kept_modes) is not int or not 0<kept_modes<=J.shape[0]):
        raise ValueError('invalid truncated-step finite input layout')
    residual=np.concatenate([g[:,0],g[:,1]])
    matrix=np.vstack([J[:,:,0].T,J[:,:,1].T])
    row=np.linalg.norm(matrix,axis=1);col=np.linalg.norm(matrix,axis=0)
    row=np.where(row>0,row,1.);col=np.where(col>0,col,1.)
    A=matrix/row[:,None]/col[None,:];rhs=-residual/row
    U,s,V=np.linalg.svd(A,full_matrices=False)
    m,n=A.shape;rank=len(s)
    if len(s)==0 or s[0]<=0:raise ValueError('truncated step has no nonzero spectrum')
    keep=s>s[0]/float(cutoff)
    if int(np.count_nonzero(keep))!=kept_modes or kept_modes<=0:
        raise ValueError('truncated step cutoff coverage changed')
    reconstruction=(U*s)@V
    mass=(abs(U)*s)@abs(V)+abs(A)
    # Each reconstruction dot has rank products/additions; charge the finite
    # orthogonal-transform recipe through its m+n rows as a residual policy.
    svd_cap=_up(_gamma(8*rank+16*(m+n))*mass)
    backward=abs(A-reconstruction)
    # Orthogonal SVD transformations have a normwise backward residual;
    # structural zero entries do not receive a componentwise relative bound.
    # Compare global Frobenius masses with the same finite operation count.
    if np.linalg.norm(backward)>np.linalg.norm(svd_cap):
        raise ValueError('truncated SVD matrix backward residual failed')
    for left,right in ((U.T,U),(V,V.T)):
        identity=left@right
        cap=_up(_gamma(8*left.shape[1]+8)*(abs(left)@abs(right)+np.eye(rank)))
        if np.linalg.norm(identity-np.eye(rank))>np.linalg.norm(cap):
            raise ValueError('truncated SVD orthogonality residual failed')
    Uk=U[:,keep];Vk=V[keep].T;sk=s[keep]
    projection=Uk.T@rhs
    q=projection/sk
    x=d*col
    recipe_x=Vk@q
    q_error=_up(_gamma(8*m+8)*(abs(Uk).T@abs(rhs))/sk)
    x_error=_up(abs(Vk)@q_error+_gamma(8*kept_modes+8)*(abs(Vk)@abs(q))+
                _gamma(4)*abs(recipe_x))
    if np.any(abs(x-recipe_x)>x_error):
        raise ValueError('truncated step retained coefficient equation failed')
    projected=Uk@projection
    projected_error=_up(abs(Uk)@(_gamma(8*m+8)*(abs(Uk).T@abs(rhs)))+
        _gamma(8*kept_modes+8)*(abs(Uk)@abs(projection)))
    linear_residual=abs(A@x-projected)
    linear_cap=_up(abs(A)@x_error+backward@abs(recipe_x)+projected_error+
        _gamma(8*n+8)*(abs(A)@abs(x)+abs(projected)))
    if np.any(linear_residual>linear_cap):
        raise ValueError('truncated step projected linear residual failed')
    if np.any(~keep):
        outside=abs(V[~keep]@x)
        outside_cap=_up(abs(V[~keep])@x_error+
            _gamma(8*n+8)*(abs(V[~keep])@abs(x)))
        if np.any(outside>outside_cap):
            raise ValueError('truncated step lies outside retained right-space')
    return {'kept_modes':kept_modes,'matrix_backward_residual_max':float(np.max(backward)),
            'matrix_backward_residual_Frobenius':float(np.linalg.norm(backward)),
            'matrix_backward_operation_allowance_Frobenius':float(np.linalg.norm(svd_cap)),
            'projected_linear_residual_max':float(np.max(linear_residual)),
            'projected_linear_operation_allowance_max':float(np.max(linear_cap)),
            'scope':'finite numerical linear-algebra residuals only',
            'physical_error_bound':None,'physical_local_gate':'OPEN'}
