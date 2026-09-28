"""Exact analytic profile identity, independent of numerical mesh and source."""
from hashlib import sha256
import json

from .nsc_ks_source_envelope import _as_metric
from .nsc_local_incoming_family import LocalAxialFunction


def _profile(profile):
    if not isinstance(profile,LocalAxialFunction):
        raise TypeError('owned analytic LocalAxialFunction required')
    off,scale=profile._derivatives[0].mapparms()
    return {'coefficients':[float(v).hex() for v in profile.coefficients],
            'center':float(profile.center).hex(),'inner':float(profile.inner).hex(),
            'outer':float(profile.outer).hex(),'map_offset':float(off).hex(),'map_scale':float(scale).hex()}


def profile_description(family, *, include_normal_window=True):
    _,metric=_as_metric(family)
    rows=[]
    for amplitude,direction in zip(metric.amplitudes,metric.directions):
        row={'amplitude':float(amplitude).hex(),'w':_profile(direction.w),'U':_profile(direction.U)}
        if include_normal_window:
            row['normal_inner']=float(direction.inner_radius).hex()
            row['normal_outer']=float(direction.outer_radius).hex()
        rows.append(row)
    return {'schema':'NSC-ANALYTIC-RADIUS-PROFILE-v1','directions':rows,
            'normal_windows_included':bool(include_normal_window)}


def profile_identity(family, *, include_normal_window=True):
    description=profile_description(family,include_normal_window=include_normal_window)
    return sha256(json.dumps(description,sort_keys=True,separators=(',',':')).encode()).hexdigest()
