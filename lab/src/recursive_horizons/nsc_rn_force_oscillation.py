"""Read-only zero-crossing diagnostics of recorded RN pressures and spinor parts.

No reconstruction of the producer, no smoothing, and no power-law fit.
Crossing positions use linear interpolation; exact-zero plateaux use their
midpoint and only count when the signs on the two sides differ.
"""
from pathlib import Path
import hashlib
import re
import numpy as np
from numba import njit, prange


@njit(parallel=True, fastmath=True, cache=True)
def crossing_kernel(values, relative_floor):
    """Rows are real fields. Return padded cell positions and row counts.

    Nonzero anchors must exceed relative_floor * row max(abs(values)).
    Subthreshold runs are treated as zero plateaux, never as extra crossings.
    """
    rows, n = values.shape
    positions = np.full((rows, n), np.nan)
    counts = np.zeros(rows, dtype=np.int64)
    for row in prange(rows):
        scale = 0.0
        for j in range(n):
            scale = max(scale, abs(values[row, j]))
        floor = relative_floor * scale
        previous = -1
        count = 0
        for j in range(n):
            v = values[row, j]
            if abs(v) <= floor:
                continue
            if previous >= 0:
                a = values[row, previous]
                if (a < 0.0) != (v < 0.0):
                    if j == previous + 1:
                        x = previous + abs(a) / (abs(a) + abs(v))
                    else:
                        x = (previous + j) / 2.0
                    positions[row, count] = x
                    count += 1
            previous = j
        counts[row] = count
    return positions, counts


@njit(parallel=True, fastmath=True, cache=True)
def association_kernel(masks):
    """Pressure vs four spinor-part edge masks: phi coefficient and proximity.

    Metrics: same-edge overlap, binary Pearson (NaN if constant), fraction of
    force crossings within two edges of a spinor crossing. Null expectation
    of the last metric uses the spinor's dilated edge-mask coverage.
    """
    result = np.empty((2, 4, 5))
    n = masks.shape[1]
    for pair in prange(8):
        f, s = pair // 4, pair % 4 + 2
        nf, ns, overlap, near, coverage = 0.0, 0.0, 0.0, 0.0, 0.0
        for j in range(n):
            a, b = masks[f, j], masks[s, j]
            nf += a
            ns += b
            overlap += a * b
            hit = False
            for k in range(max(0, j - 2), min(n, j + 3)):
                if masks[s, k]:
                    hit = True
            if hit:
                coverage += 1
                near += a
        denom = np.sqrt(nf * (n - nf) * ns * (n - ns))
        result[f, s - 2, 0] = overlap
        result[f, s - 2, 1] = (n * overlap - nf * ns) / denom if denom > 0 else np.nan
        result[f, s - 2, 2] = near / nf if nf > 0 else np.nan
        result[f, s - 2, 3] = coverage / n
        result[f, s - 2, 4] = nf * ns / n
    return result


NAMES = ('radial_pressure', 'angular_pressure', 'phi0_real', 'phi0_imag',
         'phi1_real', 'phi1_imag')
FLOORS = (0.0, 1e-12, 1e-8, 1e-4)


def finite_or_none(x):
    return float(x) if np.isfinite(x) else None


def spacing_summary(spacing, dr):
    if not len(spacing):
        return {'intervals': 0, 'mean_cells': None, 'mean_r_m': None,
                'median_cells': None, 'min_cells': None, 'max_cells': None,
                'le3_cells': 0, 'gt3_le8_cells': 0, 'gt8_cells': 0}
    return {'intervals': len(spacing), 'mean_cells': float(np.mean(spacing)),
            'mean_r_m': float(np.mean(spacing) * dr),
            'median_cells': float(np.median(spacing)),
            'min_cells': float(np.min(spacing)), 'max_cells': float(np.max(spacing)),
            'le3_cells': int(np.sum(spacing <= 3)),
            'gt3_le8_cells': int(np.sum((spacing > 3) & (spacing <= 8))),
            'gt8_cells': int(np.sum(spacing > 8))}


def profile_summary(crossings, radius):
    dr = (radius[-1] - radius[0]) / (len(radius) - 1)
    radii = radius[0] + crossings * dr
    spacing = np.diff(crossings)
    midpoints = (radii[1:] + radii[:-1]) / 2
    edges = np.array([radius[0], radius[0] + 1, 4, 8, 16, radius[-1]])
    regions = []
    for j in range(len(edges) - 1):
        lo, hi = edges[j:j+2]
        include = (radii >= lo) & ((radii <= hi) if j == len(edges)-2 else (radii < hi))
        intervals = (midpoints >= lo) & ((midpoints <= hi) if j == len(edges)-2 else (midpoints < hi))
        regions.append({'r_m': [float(lo), float(hi)], 'crossings': int(include.sum()),
                        **spacing_summary(spacing[intervals], dr)})
    return {'sign_changes': len(crossings), 'crossing_cells': crossings.tolist(),
            'crossing_r_m': radii.tolist(), **spacing_summary(spacing, dr), 'regions': regions}


def analyze_station(path):
    path = Path(path)
    with np.load(path, allow_pickle=False) as d:
        radius = np.asarray(d['radius'], dtype=np.float64)
        phi = d['coupled_phi']
        if phi.shape != (2, len(radius), 1):
            raise ValueError(f'Unexpected spinor shape: {phi.shape}')
        values = np.ascontiguousarray(np.vstack([
            d['radial_pressure'], d['angular_pressure'], phi[0, :, 0].real,
            phi[0, :, 0].imag, phi[1, :, 0].real, phi[1, :, 0].imag]), dtype=np.float64)
        actual_time = float(d['coupled_t'].item())
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(radius)):
        raise ValueError('Nonfinite station')
    dr = (radius[-1] - radius[0]) / (len(radius) - 1)
    grid_error = float(np.max(np.abs(radius - (radius[0] + np.arange(len(radius)) * dr))))
    if dr <= 0 or grid_error > 1e-10 or abs(radius[-1] - 32) > 1e-10:
        raise ValueError('Expected uniform radial domain ending at 32 r_m')
    positions, counts = crossing_kernel(values, 0.0)
    profiles = {name: profile_summary(positions[i, :counts[i]], radius)
                for i, name in enumerate(NAMES)}
    masks = np.zeros((6, len(radius)-1), dtype=np.int64)
    for i in range(6):
        index = np.minimum(np.floor(positions[i, :counts[i]]).astype(int), len(radius)-2)
        masks[i, index] = 1
    association = association_kernel(masks)
    correlations = {}
    for f in range(2):
        correlations[NAMES[f]] = {}
        for s in range(4):
            a = association[f, s]
            correlations[NAMES[f]][NAMES[s+2]] = {
                'same_edge_overlap': int(a[0]), 'binary_pearson': finite_or_none(a[1]),
                'fraction_within_2_edges': finite_or_none(a[2]),
                'uniform_edge_null_near_fraction': float(a[3]),
                'uniform_edge_null_same_edge_count': float(a[4])}
    sensitivity = {}
    for floor in FLOORS[1:]:
        p, c = crossing_kernel(values, floor)
        sensitivity[str(floor)] = {name: {'sign_changes': int(c[i]),
            **spacing_summary(np.diff(p[i, :c[i]]), dr)} for i, name in enumerate(NAMES)}
    return {'file': path.name, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'station_label': float(re.fullmatch(r'station-([0-9.]+)\.npz', path.name)[1]),
            'actual_t_over_r_m': actual_time, 'points': len(radius),
            'domain_r_m': [float(radius[0]), float(radius[-1])], 'dr_r_m': float(dr),
            'uniform_grid_max_error_r_m': grid_error, 'profiles': profiles,
            'spinor_crossing_association': correlations,
            'relative_amplitude_floor_sensitivity': sensitivity}


def analyze_family(root):
    cases = []
    for directory in sorted(Path(root).glob('mu-*-eps-*')):
        files = [p for p in directory.iterdir() if re.fullmatch(r'station-[0-9.]+\.npz', p.name)]
        files.sort(key=lambda p: float(p.stem.split('-')[1]))
        if not files:
            raise ValueError(f'No stations in {directory}')
        match = re.fullmatch(r'mu-([0-9.]+)-eps-([0-9.]+)', directory.name)
        cases.append({'case': directory.name, 'mu': float(match[1]), 'eps': float(match[2]),
                      'stations': [analyze_station(p) for p in files]})
    if len(cases) != 6:
        raise ValueError(f'Expected six cases, found {len(cases)}')
    return {'schema': 'nsc-rn-force-oscillation-v1',
            'source': 'lab/results/development/nsc-rn-neutral-family-v1',
            'method': 'Unsmooth signed crossings; linear adjacent interpolation; zero-run midpoint; interval midpoint radial bins.',
            'units': 'radius and wavelength in r_m; time in r_m; cell positions on recorded uniform grid',
            'resolution_rule': '<=3 cells grid-scale; >3 to 8 marginal; >8 spatially sampled, not a convergence proof',
            'caveats': ['Pressures are analyzed as stored, not differentiated to construct a new force.',
                        'Complex spinor has no invariant sign. Real/imag crossings are gauge/phase dependent.',
                        'Crossing correlation is descriptive, not a causal or physicality test.',
                        'No finer-grid replay or continuum error estimate is available.'],
            'numba_dispatchers': {fn.__name__: type(fn).__name__ for fn in (crossing_kernel, association_kernel)},
            'cases': cases}
