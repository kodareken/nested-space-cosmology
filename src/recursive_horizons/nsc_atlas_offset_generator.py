"""One spatial interval and one accumulated offset: explicit metric completions.

Retained history, the distance kernel and metric completion are declared
construction rules. This is not a physical dimension-selection or time-origin
theorem, nor an assertion that a finite drawing fills a volume.
"""
import math
import mpmath as mp
import numpy as np
from scipy.spatial import SphericalVoronoi,cKDTree


def opposite_winding_path(position, contraction, turns=1, radius=1.):
    """A closed forward/reverse loop and its null homotopy, not a time law.

    position is the retained drawing coordinate in [0,2]. At contraction=0
    the first half has degree turns and the second degree -turns. During
    contraction only the combined loop stays closed. The two halves need
    not remain individually closed: their contact is allowed to reconnect.
    """
    if type(turns) is not int or turns == 0:
        raise ValueError('nonzero integer winding required')
    if not math.isfinite(radius) or radius <= 0:
        raise ValueError('positive finite radius required')
    if not math.isfinite(contraction) or not 0 <= contraction <= 1:
        raise ValueError('contraction must lie in [0,1]')
    x = np.asarray(position, dtype=float)
    if not np.all(np.isfinite(x)) or np.any((x < 0) | (x > 2)):
        raise ValueError('drawing coordinate must lie in [0,2]')
    lift = np.minimum(x, 2 - x)
    return radius * np.exp(2j * np.pi * turns * (1 - contraction) * lift)


def winding_variation(contraction, turns=1, radius=1.):
    """Half the integral of |dz/dx|^2 over [0,2], in drawing units.

    This geometric functional has no supplied conversion to physical energy.
    It is the minimum for a fixed circle radius, endpoints and midpoint lift.
    """
    opposite_winding_path(0., contraction, turns, radius)
    return 4 * math.pi**2 * radius**2 * turns**2 * (1 - contraction)**2


def balanced_output_directions(count):
    """Competing balanced channel geometries; count is supplied, not selected."""
    if type(count) is not int or count < 2:
        raise ValueError('at least two output channels required')
    angles = 2 * np.pi * np.arange(count) / count
    return np.c_[np.cos(angles), np.sin(angles), np.zeros(count)]


def winding_calculation():
    """Geometry-only reference calculation: no parent action or observed scale."""
    cases = []
    for count in (2, 3, 4):
        directions = balanced_output_directions(count)
        cases.append(dict(channels=count, fraction_each=1/count,
                          normalized_vector_sum=(directions.sum(axis=0)/count).tolist(),
                          adjacent_angle_degrees=360/count))
    return dict(
        schema='nsc-offset-winding-geometry-v1',
        premises=['circle target with unit radius for this dimensionless example',
                  'integer winding of a continuous closed loop',
                  'concatenation joins opposite loops at one contact',
                  'squared-derivative functional in the stated drawing coordinate'],
        pair=dict(degrees=[1, -1], combined_degree=0,
                  variation_each=2*math.pi**2, variation_pair=4*math.pi**2,
                  conjugation_preserves_positive_variation=True),
        contraction=[dict(parameter=h, variation=winding_variation(h))
                     for h in (0., .25, .5, .75, 1.)],
        two_channel_consequence=dict(
            assumptions=['exactly two unit-direction outputs',
                         'positive scalar budgets with total one',
                         'zero weighted direction sum, with no other receiver'],
            fractions=[.5, .5], ratio=1., separation_degrees=180.),
        competing_balanced_channels=cases,
        missing_connection='contact composition/evolution law selecting the channels and transfer',
        physical_energy_unit_selected=False,
        physical_duration_selected=False,
        electron_or_positron_identification_derived=False,
        no_observational_input=True)


def offset_value(kind):
    if kind=='golden':return (3-mp.sqrt(5))/2
    if kind=='action':
        x=mp.findroot(lambda x:25*x**6-140*x+112,(mp.mpf('1.07'),mp.mpf('1.08')))
        return mp.sqrt((475*x-468)/(12800+75*x-84))
    raise ValueError('golden or retained action offset required')


def phase_history(count,kind='golden'):
    if type(count) is not int or count<1:raise ValueError('positive integer count required')
    with mp.workdps(80):
        alpha=offset_value(kind)
        return np.array([[float(mp.frac(n*alpha)),float(mp.frac(n*(n-1)*alpha/2))] for n in range(count)])


def directions(phases):
    theta=2*np.pi*phases[:,0];phi=2*np.pi*phases[:,1]
    return np.c_[np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)]


def euclidean_kernel_distance(r,s,theta,phi,theta2,phi2):
    C=math.cos(theta)*math.cos(theta2)+math.sin(theta)*math.sin(theta2)*math.cos(phi-phi2)
    return math.sqrt(max(0.,r*r+s*s-2*r*s*C))


def closed_parent_point(radial,direction,radius=1.):
    if not 0<=radial<=math.pi*radius:raise ValueError('north-to-south seed interval required')
    u=radial/radius
    return radius*np.r_[math.cos(u),math.sin(u)*np.asarray(direction)]


def sphere_cover(points):
    points=np.asarray(points)
    diagram=SphericalVoronoi(points,radius=1,center=np.zeros(3),threshold=1e-8)
    chord=cKDTree(points).query(diagram.vertices,workers=8)[0]
    return float(np.max(2*np.arcsin(np.clip(chord/2,0,1))))


def finite_record(count,kind):
    phases=phase_history(count,kind);v=directions(phases)
    line_cover=sphere_cover(np.r_[v,-v]);oriented_cover=sphere_cover(v)
    weights=np.abs(np.sin(2*np.pi*phases[:,0]))
    raw=(v.T@v)/count;weighted=(v.T@(v*weights[:,None]))/weights.sum()
    return dict(passes=count,
        line_direction_cover_radians=line_cover,
        ball_Hausdorff_distance_over_radius=math.sin(line_cover),
        oriented_direction_cover_radians=oriented_cover,
        S3_meridian_Hausdorff_distance_over_radius=min(oriented_cover,math.pi/2),
        raw_direction_second_moment=np.diag(raw).tolist(),
        area_weighted_direction_second_moment=np.diag(weighted).tolist(),
        finite_phase_Fourier_moments={f'{k},{l}':float(abs(np.mean(np.exp(2j*np.pi*(k*phases[:,0]+l*phases[:,1])))))
                                     for k,l in ((1,0),(0,1),(1,1),(2,-1))})


def calculation():
    # This is the retained action-based comparison, not a foundation input.
    from .nsc_atlas_parent_scale_selector import VolumeParent
    cases=[]
    for kind in ('action','golden'):
        with mp.workdps(80):alpha=mp.nstr(offset_value(kind),50)
        rows=[finite_record(n,kind) for n in (100,200,1000,4000)]
        cases.append(dict(offset=kind,alpha=alpha,angle=2*math.pi*float(alpha),coverage=rows))
    minimum=VolumeParent().vacuum_roots()[1];Lambda=minimum['vacuum'];R=math.sqrt(3/Lambda)
    return dict(schema='nsc-atlas-offset-generator-v1',
        generator=dict(seed='One continuous spatial interval; n counts retained offset updates, not a supplied proper-time coordinate.',
            rule='theta[n+1]=(theta[n]+alpha) mod1; phi[n+1]=(phi[n]+theta[n]) mod1; theta0=phi0=0',
            closed_form='theta_n=n alpha mod1; phi_n=n(n-1)alpha/2 mod1',
            repeated_operation='Delta theta=alpha; Delta phi=theta: the same discrete accumulation used twice.',
            inverse='theta_old=(theta_new-alpha) mod1; phi_old=(phi_new-theta_old) mod1',
            proof='Every nonzero integer character is a nonconstant polynomial in n with an irrational nonconstant coefficient; Weyl criterion gives joint equidistribution.',
            reversible=True,physical_duration_selected=False),
        metric_completion=dict(
            phase_angles='Theta=2pi theta; Phi=2pi phi',
            cosine_kernel='C_nm=cosTheta_n cosTheta_m+sinTheta_n sinTheta_m cos(Phi_n-Phi_m)',
            disk='Retain the seed interval rotated by Theta_n in one plane; its completion is the disk.',
            ball='d²=r²+s²-2rs C_nm; quotient coincident labels and complete: the Euclidean 3-ball.',
            closed_parent='d=R arccos[cos(r/R)cos(s/R)+sin(r/R)sin(s/R)C_nm], 0<=r,s<=piR; completion is round S3.',
            distinction='S3 has three intrinsic spatial dimensions despite a four-component unit-vector embedding. It is not a five-dimensional spacetime.',
            assumptions=['retain earlier drawing labels','use the stated cross-pass distance kernel','take the metric completion'],
            finite_history_has_nonzero_volume=False,
            travel_only_along_spokes_gives_same_completion=False,
            unique_three_dimension_selection=False),
        measure=dict(raw_phase_limit_second_moment=[.25,.25,.5],
            sphere_area_weight='abs(sinTheta)',area_weighted_limit_second_moment=[1/3]*3,
            ball_volume='r² abs(sinTheta) dr dTheta dPhi',
            closed_volume='R² sin²(r/R) abs(sinTheta) dr dTheta dPhi',
            physical_energy_measure_selected=False),
        cases=cases,
        inherited_parent=dict(Lambda=Lambda,radius=R,scalar_curvature=6/R**2,
            Hamiltonian_residual=6/R**2-2*Lambda,extrinsic_curvature=0,
            first_scalar_Laplacian_eigenvalues=[3/R**2,8/R**2],
            normalized_eigenvalues=[3,8],
            law='Under radius rescaling, R²(-Delta_S3_R) has the same dimensionless spectrum.',
            scope='The generated S3 can carry the existing constant-field, time-symmetric parent data. This is compatible initial geometry, not a selected formation history.'),
        physical_open=['Why the offset dynamics and metric kernel are selected',
                       'Why two accumulation levels followed by nesting are physically selected',
                       'Regional scale/energy transfer, a duration readout and an arrow',
                       'Birth of a daughter and observational abundances'],
        no_observational_fit=True)
