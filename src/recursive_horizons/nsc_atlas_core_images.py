"""Fixed-worldline images of the existing evolving scalar concentration.

The field and ray steppers are imported unchanged. The new physical input
is only a receiving worldline; the numerical angular patches invert the
previously calculated landing map, retaining its full arrival events.
"""
from dataclasses import asdict
import hashlib
import math
from pathlib import Path
import time

import numpy as np
from scipy.interpolate import CubicHermiteSpline, CubicSpline
from scipy.optimize import brentq

from .nsc_atlas_core_optics import OpticalRun, stream_rays
from .nsc_atlas_core_evolution import initial_state, moving_mesh, evolve_arrays

SCREEN = .6
OFFSET = .0002
EPOCH = 2.8
# Centres are estimates from the old, frozen fan; they are numerical
# search coordinates, not physical parameters or final image roots.
CENTRES = {1: .006517779515298096, -1: .005404164054473755}
AGES = EPOCH + np.array([-.0002, -.0001, 0., .0001, .0002])


def patches():
    angles = []; epochs = []; branches = []
    for age in AGES:
        for branch, centre in CENTRES.items():
            for angle in centre + np.linspace(-.00008, .00008, 5):
                angles.append(angle); epochs.append(age); branches.append(branch)
    return np.array(angles), np.array(epochs), np.array(branches)


def run_images(config=OpticalRun(), *, cache_dir=None):
    """One batched source replay for all inverse patches and source ages."""
    c = config; cpu = time.process_time(); wall = time.perf_counter()
    if c.kind != 'diffuse' or c.screen != SCREEN or c.prefix != 2.5 or c.end != 5.:
        raise ValueError('outside the declared image experiment')
    r = np.arange(round(c.outer/c.spacing)+1)*c.spacing
    delta, Pi, bg, model, note = initial_state(c.kind, r)
    a0 = bg[2]; mesh = moving_mesh(r)
    prefix = evolve_arrays(delta, Pi, bg, c.spacing, c.step, round(c.prefix/c.step),
        2000., 1., 9/32, c.screen, np.empty(0), round(c.prefix/c.step), mesh, a0, True, 0.)
    angles, epochs, branches = patches()
    out, status, delta, Pi, bg = stream_rays(prefix[0], prefix[1], prefix[2], c.spacing,
        c.step, round((c.end-c.prefix)/c.step), c.block_steps, mesh, a0, c.screen,
        epochs, angles, prefix[3][1], c.prefix)
    if not np.all(status == 2):
        raise RuntimeError('image rays failed to return: '+str(np.where(status != 2)[0]))
    radius = r+(bg[2]/a0-1)*mesh[0]
    digest = hashlib.sha256(np.asarray(np.column_stack((radius, delta, Pi)), dtype='<f8').tobytes()).hexdigest()
    names = ('arrival_time', 'arrival_proper_age', 'landing_angle', 'received_frequency',
        'received_radial_cosine', 'jacobi_perpendicular', 'jacobi_parallel', 'signed_beam_area',
        'rotation_jacobi_gap', 'frequency_transport_gap', 'minimum_radius', 'affine_distance',
        'rotation_jacobi', 'emission_time', 'emission_normal_energy', 'sky_angle')
    rows = [dict(zip(names, map(float, values)), emission_proper_age=float(age), branch=int(branch))
        for values, age, branch in zip(out, epochs, branches)]
    record = dict(config=asdict(c), rays=rows, final_field_sha256=digest,
        cpu_seconds=time.process_time()-cpu, wall_seconds=time.perf_counter()-wall)
    if cache_dir:
        import json
        path = Path(cache_dir); path.mkdir(parents=True, exist_ok=True)
        (path/'completed.json').write_text(json.dumps(record, indent=2, allow_nan=False)+'\n',
            encoding='utf-8', newline='\n')
    return record


def invert_patch(rows, *, offset=OFFSET, screen=SCREEN, stride=1):
    """Invert one monotone patch, checking curvature against the landing map.

    Plain cubic interpolation and curvature-slope Hermite interpolation
    are independent numerical reconstructions. Nested angular sampling
    supplies a further finite interpolation check.
    """
    rows = sorted(rows, key=lambda r: r['sky_angle'])[::stride]
    branch = rows[0]['branch']; target = math.pi-branch*offset
    x = np.array([r['sky_angle'] for r in rows])
    phi = np.array([r['landing_angle'] for r in rows])
    slope = np.array([-r['jacobi_parallel']/(screen*r['received_radial_cosine']) for r in rows])
    if not np.all(np.diff(phi) < 0):
        raise ValueError('image patch is not monotone')
    plain = CubicSpline(x, phi); curvature = CubicHermiteSpline(x, phi, slope)
    root = brentq(lambda a: plain(a)-target, x[0], x[-1], xtol=5e-16)
    independent = brentq(lambda a: curvature(a)-target, x[0], x[-1], xtol=5e-16)
    def at(key): return float(CubicSpline(x, [r[key] for r in rows])(root))
    frequency = at('received_frequency'); area = at('signed_beam_area')
    # Positive sky orientation is opposite increasing coordinate azimuth:
    # looking back reverses the photon's tangential momentum.
    image_angle = branch*math.asin(math.sin(root)/frequency)
    arrival_spline = CubicSpline(x, [r['arrival_proper_age'] for r in rows])
    endpoint_slope = screen*math.sin(root)/frequency*float(plain(root, 1))
    rotation = screen*math.sin(target)/math.sin(root)
    parallel = -screen*at('received_radial_cosine')*float(plain(root, 1))
    result = dict(branch=branch, emission_proper_age=rows[0]['emission_proper_age'],
        target_unwrapped_angle=target, launch_angle=root, received_sky_angle=image_angle,
        arrival_proper_age=at('arrival_proper_age'), received_frequency=frequency,
        jacobi_perpendicular=at('jacobi_perpendicular'), jacobi_parallel=at('jacobi_parallel'),
        signed_beam_area=area, parity=1 if area>0 else -1,
        fluence_weight=frequency/abs(area), luminosity_flux_weight=frequency**2/abs(area),
        angular_bracket=[float(x[0]), float(x[-1])],
        independent_root_gap=root-independent, landing_residual=float(plain(root)-target),
        endpoint_variation_gap=float(arrival_spline(root, 1))-endpoint_slope,
        rotation_jacobi_gap=at('jacobi_perpendicular')-rotation,
        parallel_jacobi_gap=at('jacobi_parallel')-parallel,
        frequency_transport_gap=at('frequency_transport_gap'))
    return result


def paired_observables(positive, negative):
    """Same emitted pulse; two arrivals on one receiver's proper clock."""
    return dict(image_separation=positive['received_sky_angle']-negative['received_sky_angle'],
        delay_negative_minus_positive=negative['arrival_proper_age']-positive['arrival_proper_age'],
        frequency_ratio_negative_over_positive=negative['received_frequency']/positive['received_frequency'],
        fluence_ratio_negative_over_positive=negative['fluence_weight']/positive['fluence_weight'],
        luminosity_flux_ratio_negative_over_positive=negative['luminosity_flux_weight']/positive['luminosity_flux_weight'],
        predicted_delay_drift=1/negative['received_frequency']-1/positive['received_frequency'])


def reduce_run(record):
    pairs = []; images = []
    for age in AGES:
        pair = []
        for branch in (1, -1):
            rows = [r for r in record['rays'] if r['emission_proper_age']==age and r['branch']==branch]
            image = invert_patch(rows); sparse = invert_patch(rows, stride=2)
            image['nested_angular_gaps'] = {key: image[key]-sparse[key] for key in
                ('launch_angle', 'arrival_proper_age', 'received_frequency', 'fluence_weight', 'received_sky_angle')}
            pair.append(image); images.append(image)
        pairs.append(dict(emission_proper_age=float(age), **paired_observables(*pair)))
    delays = np.array([p['delay_negative_minus_positive'] for p in pairs])
    d1 = (delays[4]-delays[0])/.0004; d2 = (delays[3]-delays[1])/.0002
    derivative = (4*d2-d1)/3
    check = dict(delay_derivative=float(derivative), predicted=pairs[2]['predicted_delay_drift'],
        identity_gap=float(derivative-pairs[2]['predicted_delay_drift']),
        derivative_refinement_gap=float(d2-d1))
    return dict(images=images, pairs=pairs, delay_drift_check=check)
