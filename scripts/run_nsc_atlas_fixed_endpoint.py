#!/usr/bin/env python3
"""One fixed-grid march from the quiet anchor to the original source endpoint.

Summaries, hashes and nine shared source-row cuts only. No raw history and
no final coefficient arrays. The step sizes are a numerical design, not a
stability proof. An instability stops the run; nothing is retried or remeshed.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import io
import json
import math
import sys
from pathlib import Path
import subprocess
import time

import numpy as np

from diagnose_nsc_atlas_boundary_rows import authenticate, constraint_terms
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_characteristic_trace import SourceSnapshot
from recursive_horizons.nsc_atlas_coupled_evolution import CoupledEvolutionError
from recursive_horizons.nsc_atlas_endpoint_comparison import jordan_radius, ssprk3_diagnostic_increment
from recursive_horizons.nsc_atlas_fixed_cap_stage import quiet_fixed_state, coupled_fixed_stage
from recursive_horizons.nsc_atlas_fixed_cap_time import (
    advance_fixed_ssprk3, fixed_time_residual, fixed_derived_residual,
)
from recursive_horizons.nsc_atlas_reference_wall_gauge import ReferenceWallGauge
from run_nsc_atlas_fixed_onset import (
    ROOT, ENRICHED_REFERENCE, frozen_sources, prepare_history, plain, digest_state, producer,
)

ORIGINAL_ALPHA = 0.00012983127142551443
X_ANCHOR = .4
X_LAST = 0.4815042082576364
REFERENCE_RADIUS_ANCHOR = .45
REFERENCE_RADIUS_ENDPOINT = .55
SOURCE_ROWS = 257
CUT_STRIDE = 32
CHARGED_EVALUATIONS = 500_000_000_000
CPU_SECONDS = 10000.
FRAME_CHECK_SCALE = 1e-12
LAST_WIDTHS = (Q('0.8'), Q('0.73'), Q('0.66'))
INITIAL_WALL_CUT = Q('0.54')
LEVELS = (
    dict(cap_cells=2, exterior_cells=32, last_width=LAST_WIDTHS[0], max_dx=1e-4, step_budget=2000),
    dict(cap_cells=4, exterior_cells=64, last_width=LAST_WIDTHS[1], max_dx=1e-4/2, step_budget=4000),
    dict(cap_cells=8, exterior_cells=128, last_width=LAST_WIDTHS[2], max_dx=1e-4/4, step_budget=8000),
)
# Stations that exist on every admitted exterior grid. Indices scale with the cell count.
COMMON_ETA_OVER_32 = (12, 14, 16, 18, 26, 28, 30)
TUBE = dict(r=[.2, .8], L=[-2., 2.], c=[-4., 4.])


def level_spec(level):
    if type(level) is not int or not 0 <= level < len(LEVELS):
        raise ValueError('level must be 0, 1 or 2')
    return LEVELS[level]


def _ulp_gap(left, right):
    return 64*math.ulp(max(1., abs(left), abs(right)))


def _prefer(value, anchors, endpoints):
    for item in endpoints:
        if abs(value-item) <= _ulp_gap(value, item):
            return float(item)
    for item in anchors:
        if abs(value-item) <= _ulp_gap(value, item):
            return float(item)
    return float(value)


def _is_endpoint(value, endpoints):
    return any(abs(value-item) <= _ulp_gap(value, item) for item in endpoints)


def _is_anchor(value, anchors):
    return any(abs(value-item) <= _ulp_gap(value, item) for item in anchors)


def merge_boundaries(points, *, anchors, endpoints):
    """Drop only roundoff-coincident boundaries. Endpoints outrank source rows."""
    anchors = tuple(float(v) for v in anchors)
    endpoints = tuple(float(v) for v in endpoints)
    ordered = sorted(set(float(v) for v in points))
    if not ordered:
        raise ValueError('empty boundary set')
    merged = []
    for value in ordered:
        chosen = _prefer(value, anchors, endpoints)
        if not merged:
            merged.append(chosen)
            continue
        if chosen <= merged[-1] or chosen-merged[-1] <= _ulp_gap(merged[-1], chosen):
            previous = merged[-1]
            if _is_endpoint(chosen, endpoints) or _is_endpoint(previous, endpoints):
                merged[-1] = chosen if _is_endpoint(chosen, endpoints) else previous
            elif _is_anchor(chosen, anchors) and not _is_anchor(previous, anchors):
                merged[-1] = chosen
            elif _is_anchor(previous, anchors) and not _is_anchor(chosen, anchors):
                merged[-1] = previous
            continue
        merged.append(chosen)
    if any(b <= a for a, b in zip(merged, merged[1:])):
        raise ValueError('merged boundaries are not strictly increasing')
    return merged


def time_partition(source_x, front, back, *, x_left, x_right, max_dx, step_budget):
    """Split at source rows, the source front and the source back.

    Near-coincident values collapse so a one-ulp duplicate cannot become a step.
    ``max_dx`` is an admitted numerical width, not a proved stability restriction.
    """
    rows = [float(v) for v in source_x]
    if len(rows) < 2 or any(b <= a for a, b in zip(rows, rows[1:])):
        raise ValueError('source rows must increase')
    x_left, x_right = float(x_left), float(x_right)
    max_dx = float(max_dx)
    if type(step_budget) is not int or step_budget < 1:
        raise ValueError('step budget must be a positive integer')
    if not x_left < x_right or not 0 < max_dx:
        raise ValueError('require a positive endpoint interval and step')
    raw = [float(v) for v in (*rows, float(front), float(back), x_left, x_right)
           if x_left-_ulp_gap(x_left, v) <= float(v) <= x_right+_ulp_gap(x_right, v)]
    bounds = merge_boundaries(raw, anchors=rows, endpoints=(x_left, x_right))
    if abs(bounds[0]-x_left) > _ulp_gap(bounds[0], x_left) or abs(bounds[-1]-x_right) > _ulp_gap(bounds[-1], x_right):
        raise ValueError('partition misses the endpoint interval')
    bounds[0], bounds[-1] = x_left, x_right
    pieces = []
    total = 0
    for left, right in zip(bounds, bounds[1:]):
        width = right-left
        if not width > 0:
            raise ValueError('nonincreasing partition')
        count = max(1, math.ceil(width/max_dx-1e-15))
        while count <= step_budget and width/count > max_dx:
            count += 1
        midpoint = left+width/2
        if not left < midpoint < right:
            raise ValueError('partition piece has no interior midpoint')
        pieces.append((left, right, count))
        total += count
    if total > step_budget:
        raise ValueError('derived time schedule exceeds admitted steps')
    return pieces, total


def step_nodes(left, right, count, max_dx):
    nodes = np.linspace(float(left), float(right), int(count)+1)
    nodes[0], nodes[-1] = float(left), float(right)
    gaps = np.diff(nodes)
    if np.any(gaps <= 0) or float(np.max(gaps)) > float(max_dx)+4*math.ulp(float(max_dx)):
        raise ValueError('time nodes exceed the admitted step')
    return nodes


def common_source_cuts(source_x):
    """Nine cuts shared by every level: stored rows 0,32,...,256."""
    rows = np.asarray(source_x, dtype=float)
    if rows.ndim != 1 or len(rows) != SOURCE_ROWS or np.any(np.diff(rows) <= 0):
        raise ValueError('common cuts require the original 257 increasing source rows')
    cuts = rows[::CUT_STRIDE]
    if len(cuts) != 9 or cuts[0] != rows[0] or cuts[-1] != rows[-1]:
        raise ValueError('stride 32 must select nine source rows including both endpoints')
    return cuts.copy()


def reference_midpoint():
    return Q(float((math.atanh(REFERENCE_RADIUS_ANCHOR)+math.atanh(REFERENCE_RADIUS_ENDPOINT))/2))


def cap_knots(level):
    """Fixed cells with an initial wall cut of 0.54 and equal inner cells.

    Widths .8, .73, .66 refine both the interior and the wall cell while
    keeping the clipped Gram matrix inside its existing accuracy domain.
    This is not an asymptotic uniform-h family.
    """
    spec = level_spec(level)
    width = spec['last_width']
    last_left = Q(math.atanh(REFERENCE_RADIUS_ANCHOR))-INITIAL_WALL_CUT*width
    last_right = last_left+width
    interior = spec['cap_cells']-1
    if not 0 < last_left < last_right or interior < 1:
        raise ValueError('last cell does not leave a positive interior')
    step = last_left/interior
    knots = [Q(0)]+[step*i for i in range(1, interior)]+[last_left, last_right]
    if any(b <= a for a, b in zip(knots, knots[1:])):
        raise ValueError('cap knots do not increase')
    return tuple(knots)


def cell_widths(knots):
    return tuple(b-a for a, b in zip(knots, knots[1:]))


def last_cell_cut(knots, radius):
    radius = float(radius)
    if not 0 < radius < 1:
        raise ValueError('cut radius must lie in (0,1)')
    Z = Q(float(math.atanh(radius)))
    left, right = knots[-2], knots[-1]
    return (Z-left)/(right-left)


def grid_report(level):
    knots = cap_knots(level)
    widths = cell_widths(knots)
    cuts = tuple(last_cell_cut(knots, radius) for radius in
                 (REFERENCE_RADIUS_ANCHOR, REFERENCE_RADIUS_ENDPOINT))
    guard = min(min(cuts)-Q(1, 4), Q(3, 4)-max(cuts))
    if guard <= 0:
        raise ValueError('reference wall does not sit inside the last cell with a finite guard')
    return dict(knots=knots, widths=widths, width_max=max(widths), width_min=min(widths),
                reference_cuts=cuts, cut_guard=guard, uniform_h_family=max(widths) == min(widths),
                remesh=False, adapt=False, activation=False)


def eta_patch(exterior_cells):
    if type(exterior_cells) is not int or exterior_cells % 32:
        raise ValueError('exterior cells must be a positive multiple of 32')
    scale = exterior_cells//32
    return np.array([k*scale for k in COMMON_ETA_OVER_32], dtype=int)


def history_label(central_time, radius, clock_offset):
    return float(central_time)+float(clock_offset)+float(math.atanh(radius))


def require_label_inside(last_label, first_label, upper, lower, *, quiet_extension):
    margin = float(last_label)-float(upper)
    if margin < 0:
        raise ValueError('history label would extrapolate past the prepared domain')
    if not quiet_extension and float(lower) < float(first_label):
        raise ValueError('history label would extrapolate below the prepared domain')
    return margin


def physical_gaps(residual):
    """Physical junction gaps. Weak algebraic balance is not one of them."""
    return dict(mass_actual=float(residual['physical_temporal_mass_actual']),
                mass_matched=float(residual['physical_temporal_mass_matched']),
                flux=float(residual['conformal_wall_flux']))


def matched_mass_rate(*, radius, radius_rate, alpha, moving_density, clock_rate, metric_A, lapse, Pi, Phi):
    """Matched transport dM_E/dx. ``moving_density`` is the stage's rx*rho term."""
    transport = float(clock_rate)*(float(metric_A)**2*math.exp(float(lapse))/2)*float(radius)**2*float(Pi)*float(Phi)
    return 1.5*float(radius)**2*float(radius_rate)+float(alpha)*(float(moving_density)+transport)


def physical_mass(radius, normalized_mass, alpha):
    return float(radius)**3/2+float(alpha)*float(normalized_mass)


def ssprk3_increment(width, left, mid, right):
    """Simpson weights realized by SSPRK3 on a rate sampled along the step."""
    return float(width)*(float(left)+4*float(mid)+float(right))/6


def frame_observables(*, radius, scalar_c, canonical_f, alpha, normalized_mass, central_T, einstein_s):
    """Jordan radius from the invariant scalar, cross-checked against canonical cosh.

    ``einstein_s`` is the stored Einstein wall proper time, not a Jordan clock.
    """
    radius, scalar_c = float(radius), float(scalar_c)
    canonical_f, alpha = float(canonical_f), float(alpha)
    jordan = jordan_radius(radius, scalar_c, alpha)
    q = math.sqrt(alpha/6)*canonical_f
    cosh_radius = radius*math.cosh(q)
    sinhc = math.sinh(q)/q if q else 1.
    scalar_match_defect = scalar_c-radius*canonical_f*sinhc
    gap = jordan-cosh_radius
    scale = FRAME_CHECK_SCALE*max(1., abs(radius), abs(cosh_radius))
    return dict(
        einstein_radius=radius, jordan_radius=jordan, canonical_cosh_radius=cosh_radius,
        frame_formula_gap=gap, scalar_match_defect=scalar_match_defect,
        frame_formula_needs_check=abs(gap) > scale or abs(scalar_match_defect) > scale,
        frame_formula='primary R_J=hypot(R_E, sqrt(alpha/6)*c); cross-check R_E*cosh(sqrt(alpha/6)*f)',
        central_proper_clock_T=float(central_T), einstein_wall_proper_time_s=float(einstein_s),
        wall_proper_time_frame='Einstein', c=scalar_c, f=canonical_f,
        normalized_mass=float(normalized_mass), vacuum_mass=radius**3/2,
        physical_mass=physical_mass(radius, normalized_mass, alpha),
        vacuum_subtraction='M_E=R_E**3/2+alpha*normalized_mass', units='dimensionless')


def stage_observables(row, alpha):
    wall, cap, hybrid = row['physical_wall'], row['connection']['cap'], row['connection']['hybrid']
    rates = row['wall_rates']
    observables = frame_observables(radius=wall['r'], scalar_c=wall['c'], canonical_f=cap['wall']['f'],
                                    alpha=alpha, normalized_mass=cap['mass'], central_T=wall['T'],
                                    einstein_s=wall['s'])
    observables.update(
        matched_mass_rate=matched_mass_rate(radius=wall['r'], radius_rate=rates['r'], alpha=alpha,
                                            moving_density=cap['moving_density'], clock_rate=rates['T'],
                                            metric_A=cap['wall']['A'], lapse=cap['lapse'],
                                            Pi=hybrid['Pi'], Phi=hybrid['Phi']),
        hybrid_c_rate=float(rates['c']), shared_c_rate=float(row['shared_wall_rates']['c']),
        temporal_residual=physical_gaps(row['connection']['residual']),
        boundary_residuals=dict(conformal_wall_flux=float(row['connection']['residual']['conformal_wall_flux']),
                              kinematic={k: float(v) for k, v in row['wall_kinematic_defect'].items()},
                              incoming_jet={k: float(v) for k, v in row['wall_incoming_jet_defect'].items()}),
        weak_balance_defect=float(cap['weak_balance_defect']),
        metric_certificate_floor=float(cap['metric_certificate']['A_lower_float_down']))
    return observables


class RunningTrapezoid:
    """Composite trapezoid on accepted endpoints. Diagnostic only."""

    def __init__(self):
        self._x = None
        self._values = None
        self.signed, self.absolute = {}, {}

    def add(self, x, values):
        x = float(x)
        values = {key: float(value) for key, value in values.items()}
        if self._x is not None:
            dx = x-self._x
            for key, value in values.items():
                self.signed[key] = self.signed.get(key, 0.)+.5*(self._values[key]+value)*dx
                self.absolute[key] = self.absolute.get(key, 0.)+.5*(abs(self._values[key])+abs(value))*dx
        self._x, self._values = x, values


class TransportComparison:
    """Compare the accepted state with an SSPRK3 integral. Never corrects it."""

    def __init__(self):
        self.mass_predicted = None
        self.scalar_predicted = None
        self.max_abs_mass_gap = 0.
        self.max_abs_scalar_gap = 0.
        self.last = None

    def start(self, mass, scalar):
        self.mass_predicted = float(mass)
        self.scalar_predicted = float(scalar)

    def advance(self, mass_increment, scalar_increment, actual_mass, actual_scalar):
        if self.mass_predicted is None:
            raise ValueError('transport comparison has no initial value')
        self.mass_predicted += float(mass_increment)
        self.scalar_predicted += float(scalar_increment)
        mass_gap = float(actual_mass)-self.mass_predicted
        scalar_gap = float(actual_scalar)-self.scalar_predicted
        self.max_abs_mass_gap = max(self.max_abs_mass_gap, abs(mass_gap))
        self.max_abs_scalar_gap = max(self.max_abs_scalar_gap, abs(scalar_gap))
        self.last = dict(mass_predicted=self.mass_predicted, mass_actual=float(actual_mass),
                         mass_signed_gap=mass_gap, mass_absolute_gap=abs(mass_gap),
                         scalar_predicted=self.scalar_predicted, scalar_actual=float(actual_scalar),
                         scalar_signed_gap=scalar_gap, scalar_absolute_gap=abs(scalar_gap),
                         state_corrected=False)


def assert_summary(value, path='report'):
    if isinstance(value, np.ndarray):
        raise ValueError('raw array in '+path)
    if isinstance(value, dict):
        for key, item in value.items():
            assert_summary(item, path+'.'+str(key))
        return
    if isinstance(value, (list, tuple)):
        if len(value) > 16 and value and isinstance(value[0], (int, float)):
            raise ValueError('raw numeric series in '+path)
        for index, item in enumerate(value):
            assert_summary(item, path+'['+str(index)+']')


def _finite_scalar(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('nonfinite summary scalar')
    return value


def _residual_scalars(dense, *, basis, budget, background):
    residual, physical = dense['residual'], dense['stage']
    Z = float(np.arctanh(physical['physical_wall']['r']))
    gram = basis.gram(Z, budget=budget)
    out = {}
    for key in ('chi', 'canonical_density'):
        out['cap_L2.'+key] = math.sqrt(gram.squared_norm(residual['cap'][key], budget=budget))
    for group in ('departure', 'wall'):
        for key, value in residual[group].items():
            out[group+'.'+key] = float(np.max(np.abs(value)))
    derived = fixed_derived_residual(residual, physical, basis=basis, budget=budget, background=background)
    for key, value in derived.items():
        out['derived.'+key] = _finite_scalar(value)
    return out


def initial_state(basis, exterior_cells):
    """Vacuum departure: zero cap coefficients and zero exterior/wall departures."""
    return quiet_fixed_state(basis, exterior_cells=exterior_cells)


def save_failed_stage(state, *, x, level, source_commit, output):
    """Keep one small reproduction seed locally; commit only its identity.

    No trajectory is stored. A single failed stage lets an independent
    radial check inspect the actual discrepancy without rerunning its past.
    """
    payload=dict(x=np.array(float(x)),level=np.array(int(level)),
                 source_commit=np.array(source_commit))
    count=0
    for group,values in state.items():
        for key,value in values.items():
            array=np.asarray(value,dtype=float)
            if not np.all(np.isfinite(array)):
                raise ValueError('failed-stage seed contains nonfinite state')
            payload[group+'.'+key]=array
            count+=array.size
    if count>4096:
        raise ValueError('failed-stage seed exceeds its compact-state budget')
    stream=io.BytesIO()
    np.savez_compressed(stream,**payload)
    raw=stream.getvalue()
    if len(raw)>65536:
        raise ValueError('failed-stage seed exceeds its byte budget')
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('xb') as destination:
        destination.write(raw)
    return dict(path=output.relative_to(ROOT).as_posix(),bytes=len(raw),
                sha256=hashlib.sha256(raw).hexdigest(),numeric_state_values=count,
                raw_trajectory_saved=False,git_tracked=False,source_x=float(x))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--level', type=int, required=True, choices=(0, 1, 2))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--diagnostic-step-limit',type=int,
                        help='admitted replay limit; cannot produce a completed endpoint')
    args = parser.parse_args()
    if args.diagnostic_step_limit is not None and args.diagnostic_step_limit<=0:
        raise ValueError('diagnostic step limit must be positive')
    if args.output.exists():
        raise ValueError('preserve the existing result; choose a new admitted output')
    if not args.output.parent.is_dir():
        raise ValueError('output parent must already exist')
    spec = level_spec(args.level)
    commit, bindings = frozen_sources()
    budget = Budget(maximum_evaluations=CHARGED_EVALUATIONS, cpu_seconds=CPU_SECONDS)
    report = dict(schema='nsc-atlas-fixed-endpoint-v1', status='incomplete', level=args.level,
                  runtime=producer.runtime_versions(), source_commit=commit, source_bindings=bindings,
                  scope=('one enriched-history fixed-cap march on the original source interval; '
                         'sampled defects and diagnostic quadratures; no continuum enclosure, '
                         'stability proof, remesh, or observational identification'),
                  limits=dict(charged_values=CHARGED_EVALUATIONS, cpu_seconds=CPU_SECONDS,
                              accepted_steps=spec['step_budget'], per_level=True, budget_reset=False))
    phase = 'input authentication'
    if args.diagnostic_step_limit is not None:
        report['diagnostic_step_limit']=args.diagnostic_step_limit
    calls = accepted = 0
    maxima = {}
    minimum_A = None
    minimum_label_margin = None
    last_x = None
    history = hashlib.sha256()
    samples = [None]*9
    trapezoid = RunningTrapezoid()
    transport = TransportComparison()

    def maximum(key, value):
        maxima[key] = max(maxima.get(key, 0.), float(np.max(np.abs(value))))

    try:
        _, _, _, identity = authenticate(ROOT, budget)
        raw = (ROOT/producer.BASE_INPUT).read_bytes()
        if hashlib.sha256(raw).hexdigest() != identity['source_bindings'][producer.BASE_INPUT]:
            raise ValueError('source input changed')
        base = json.loads(raw)
        geometry, data = base['geometry'], base['initial_data']['finite']
        if data['alpha'] != ORIGINAL_ALPHA:
            raise ValueError('finite source alpha is not the original amplitude')
        x_right = float(geometry['x_last'])
        if geometry['x_corner'] != X_ANCHOR or abs(x_right-X_LAST) > _ulp_gap(x_right, X_LAST):
            raise ValueError('source anchor or endpoint changed')
        rows = {key: np.asarray(value, dtype=float) for key, value in data['exact_rows'].items()}
        source_x = rows['x']
        if (len(source_x) != SOURCE_ROWS or abs(float(source_x[0])-X_ANCHOR) > _ulp_gap(float(source_x[0]), X_ANCHOR)
                or abs(float(source_x[-1])-x_right) > _ulp_gap(float(source_x[-1]), x_right)):
            raise ValueError('endpoint march requires all 257 original source rows')
        cuts = common_source_cuts(source_x)
        front, back = (float(v) for v in geometry['support_x'])
        pieces, total_steps = time_partition(source_x, front, back, x_left=X_ANCHOR, x_right=x_right,
                                             max_dx=spec['max_dx'], step_budget=spec['step_budget'])
        if args.diagnostic_step_limit is not None and args.diagnostic_step_limit>=total_steps:
            raise ValueError('diagnostic replay must stop before the full endpoint')
        widths = [right-left for left, right, _ in pieces]
        report['schedule'] = dict(piece_count=len(pieces), accepted_step_budget=total_steps,
                                  max_dx=spec['max_dx'], min_piece_width=min(widths), max_piece_width=max(widths),
                                  merged_roundoff_boundaries=len({float(v) for v in (*source_x, front, back, X_ANCHOR, x_right)
                                                                   if X_ANCHOR-_ulp_gap(X_ANCHOR, v) <= float(v) <= x_right+_ulp_gap(x_right, v)})-len(pieces)-1,
                                  split='source rows, source front and source back; ulp-coincident boundaries merged',
                                  stability_proof=False, retries=0)
        phase = 'fixed physical cap'
        grid = grid_report(args.level)
        basis = FixedCapBasis(grid['knots'], budget=budget)
        exterior_cells = spec['exterior_cells']
        eta = np.linspace(0., 1., exterior_cells+1)
        patch = eta_patch(exterior_cells)
        report['geometry'] = dict(x_left=X_ANCHOR, x_right=x_right, reference_radius_endpoint=REFERENCE_RADIUS_ENDPOINT,
                                  cap_cells=spec['cap_cells'], cap_dimension=basis.dimension, exterior_cells=exterior_cells,
                                  exact_cap_knots=[str(v) for v in grid['knots']],
                                  width_min=str(grid['width_min']), width_max=str(grid['width_max']),
                                  reference_last_cell_cuts=[str(v) for v in grid['reference_cuts']],
                                  last_cell_cut_guard=str(grid['cut_guard']), uniform_h_family=False,
                                  remesh=False, adapt=False, activation=False,
                                  common_source_cuts=cuts.tolist(),
                                  common_cut_rule='stored source rows [::32], shared by all three levels',
                                  selected_exterior_eta=(patch/exterior_cells).tolist())
        report['input'] = dict(path=producer.BASE_INPUT,
                               sha256=hashlib.sha256(raw).hexdigest(),
                               operator_sha256=identity['sha256'], manifest_sha256=identity['manifest_sha256'],
                               authenticated_input_files=104, source_rows=SOURCE_ROWS, alpha=ORIGINAL_ALPHA)
        phase = 'one enriched history'
        report['preparation'] = {}
        background, clock_offset = prepare_history(budget, report['preparation'], enrich_ordinary_pieces=True)
        reference_raw = (ROOT/ENRICHED_REFERENCE).read_bytes()
        pinned = subprocess.check_output(['git', 'show', commit+':'+ENRICHED_REFERENCE], cwd=ROOT)
        if reference_raw != pinned:
            raise ValueError('enriched reference differs from committed receipt')
        checked = json.loads(reference_raw)
        if (checked['schema'] != 'nsc-atlas-prepared-endpoint-enriched-v1'
                or checked['status'] != 'completed_reference_comparison'
                or set(checked['acceptance']) != {'history_material_rate', 'conformal_wall_gap',
                                                 'prepared_prefix_gap', 'velocity_channel', 'density_force_L2'}
                or not all(item['passed'] for item in checked['acceptance'].values())
                or any(report['preparation'][key] != checked['preparation'][key]
                       for key in ('input_sha256', 'original_input_sha256'))):
            raise ValueError('enriched history differs from its checked preparation')
        report['enriched_reference'] = dict(path=ENRICHED_REFERENCE,
                                            sha256=hashlib.sha256(reference_raw).hexdigest(),
                                            field_nodes=257, geometry_knots=257)
        phase = 'strict history domain'
        labels = background.inputs['labels']
        first_label, last_label = float(labels[0]), float(labels[-1])
        quiet_extension = bool(background.profile.quiet_extension)
        gauge = ReferenceWallGauge(geometry, order=64, budget=budget)
        snapshot = SourceSnapshot(geometry, alpha=ORIGINAL_ALPHA, x_nodes=source_x, source_rows=rows, budget=budget)
        reference = gauge.sample(source_x)
        upper = reference['T']+clock_offset+np.arctanh(reference['r'])
        lower = reference['T']+clock_offset-np.arctanh(reference['r'])
        # Reference quadrature and history labels have independent rounding.
        # Actual finite-state queries retain the strict no-extrapolation guard.
        reference_margin = float(np.min(last_label-upper))
        onset_label = float(labels[background.onset['index']])
        initial_label = history_label(reference['T'][0], reference['r'][0], clock_offset)
        if not initial_label < onset_label:
            raise ValueError('initial cap history is not entirely quiet')
        report['label_margin'] = dict(history_first_label=first_label, history_last_label=last_label,
                                      quiet_extension=quiet_extension, extrapolation=False,
                                      reference_minimum_upper=reference_margin,
                                      reference_endpoint_upper=float(last_label-upper[-1]),
                                      initial_history_quiet_margin=onset_label-initial_label)
        state = initial_state(basis, exterior_cells)
        last_x = X_ANCHOR
        history.update(digest_state(dict(x=X_ANCHOR, state=state)).encode())
        print(json.dumps(dict(event='prepared', level=args.level, steps=total_steps,
                              cap_cells=spec['cap_cells'], exterior_cells=exterior_cells)), flush=True)

        def stage(value, *, x, budget):
            nonlocal calls, minimum_A, minimum_label_margin
            calls += 1
            source = snapshot.sample(x)
            ref = gauge.sample(x)
            try:
                row = coupled_fixed_stage(value, basis=basis, x=x, eta=eta, source=source, gauge=ref,
                                          alpha=ORIGINAL_ALPHA, tube=TUBE, budget=budget, background=background,
                                          background_clock_offset=clock_offset)
            except (ValueError,RuntimeError,FloatingPointError,OverflowError,KeyError,TypeError):
                try:
                    report['failed_stage_seed']=save_failed_stage(
                        value,x=x,level=args.level,source_commit=commit,
                        output=ROOT/'.build/results'/(args.output.stem+'-failed-stage.npz'))
                except (OSError,ValueError) as capture_error:
                    report['failed_stage_seed_capture_error']=str(capture_error)
                raise
            connection = row['connection']
            cap = connection['cap']
            floor = cap['metric_certificate']['A_lower_float_down']
            minimum_A = floor if minimum_A is None else min(minimum_A, floor)
            for key, item in connection['residual'].items():
                maximum('connection.'+key, item)
            gaps = physical_gaps(connection['residual'])
            maximum('temporal.physical_temporal_mass_actual', gaps['mass_actual'])
            maximum('temporal.physical_temporal_mass_matched', gaps['mass_matched'])
            maximum('boundary.conformal_wall_flux', gaps['flux'])
            for group in ('wall_kinematic_defect', 'shared_wall_update_defect', 'wall_incoming_jet_defect'):
                for key, item in row[group].items():
                    maximum(group+'.'+key, item)
                    if group != 'shared_wall_update_defect':
                        maximum('boundary.'+group+'.'+key, item)
            for key in ('source_constraint', 'source_h_xx_defect', 'source_constraint_jordan',
                        'source_constraint_from_curvature_defect'):
                maximum('source.'+key, source[key])
            maximum('source.stage_minus_snapshot', row['source_constraint']-source['source_constraint'])
            maximum('cap.weak_balance_defect', cap['weak_balance_defect'])
            queried_upper = history_label(row['physical_wall']['T'], row['physical_wall']['r'], clock_offset)
            queried_lower = float(row['physical_wall']['T'])+clock_offset-float(np.arctanh(row['physical_wall']['r']))
            margin = require_label_inside(last_label, first_label, queried_upper, queried_lower,
                                          quiet_extension=quiet_extension)
            minimum_label_margin = margin if minimum_label_margin is None else min(minimum_label_margin, margin)
            if x <= front:
                maximum('quiet.cap_wall_scalar', row['physical_wall']['c'])
                maximum('quiet.exterior_scalar_departure', value['departure']['c'])
            vac = row['exterior_prefix']['vacuum']
            offset = float(source['r'])-float(vac['r'][0])
            delta_r = np.r_[0., value['departure']['r'], row['physical_wall']['r']-vac['r'][-1]-offset]
            second = (delta_r[2:]-2*delta_r[1:-1]+delta_r[:-2])*exterior_cells**2
            pure = constraint_terms({k: row['exterior'][k][patch] for k in ('r', 'c')},
                                    {k: vac[k][patch] for k in ('r_Ubar', 'r_Ubar_Ubar', 'log_C_Ubar')},
                                    {k: row['departure_gradients'][k][patch] for k in ('r', 'L', 'c')},
                                    second[patch-1], float(ref['W']), ORIGINAL_ALPHA, budget)
            maximum('exterior.selected_pure_U_constraint', pure['E_Ubar'])
            maximum('identity.selected_pure_U_decomposition', pure['algebraic_identity_defect'])
            for key in ('r', 'L', 'c'):
                defect = row['wall_kinematic_defect'][key]+float(ref['W_x'])*row['wall_incoming_jet_defect'][key]
                maximum('identity.wall.'+key, defect)
                if abs(defect) > 2e-11:
                    raise ValueError('wall jet/kinematic identity failed')
            radius = row['physical_wall']['r']
            rootD = math.sqrt(float(connection['hybrid']['D']))
            defect = connection['residual']['conformal_wall_flux']-radius*rootD*connection['residual']['scalar_tangent']
            maximum('identity.conformal_scalar', defect)
            if abs(defect) > 2e-11:
                raise ValueError('conformal scalar tangent identity failed')
            return row

        def remember(cut_index, row, *, kind, residual, midpoint_residual=None):
            current = samples[cut_index]
            if current is None:
                observables = stage_observables(row, ORIGINAL_ALPHA)
                observables.update(normalized_source_row=cut_index/8, x=float(cuts[cut_index]),
                                   residual_kind=kind, time_residual=residual,
                                   arriving_midpoint_residual=midpoint_residual,
                                   label_margin=_finite_scalar(last_label-history_label(
                                       row['physical_wall']['T'], row['physical_wall']['r'], clock_offset)))
                samples[cut_index] = observables
                return
            if residual is not None and current.get('time_residual') is None:
                current['time_residual'] = residual
                current['residual_kind'] = kind
            if midpoint_residual is not None and current.get('arriving_midpoint_residual') is None:
                current['arriving_midpoint_residual'] = midpoint_residual

        final_row = None
        for left, right, count in pieces:
            nodes = step_nodes(left, right, count, spec['max_dx'])
            for index in range(count):
                if args.diagnostic_step_limit is not None and accepted>=args.diagnostic_step_limit:
                    phase='declared diagnostic step limit'
                    raise RuntimeError('diagnostic replay reached its admitted accepted-step limit')
                xa, xb = float(nodes[index]), float(nodes[index+1])
                xm = xa+(xb-xa)/2
                phase = 'common step'
                step = advance_fixed_ssprk3(state, xa, xm, xb, basis=basis, exterior_cells=exterior_cells,
                                            stage=stage, budget=budget)
                if final_row is None:
                    first = step['stages'][0]
                    report['initial'] = dict(maximum_absolute_rate=max(float(np.max(abs(v)))
                                                                       for group in first['rates'].values()
                                                                       for v in group.values()),
                                             source_constraint=float(first['source_constraint']),
                                             observables=stage_observables(first, ORIGINAL_ALPHA))
                    if report['initial']['maximum_absolute_rate'] > 1e-10:
                        raise ValueError('initial quiet value/first-jet recovery failed')
                    transport.start(report['initial']['observables']['physical_mass'],
                                    report['initial']['observables']['c'])
                    trapezoid.add(xa, dict(physical_temporal_mass_actual=report['initial']['observables']['temporal_residual']['mass_actual'],
                                           physical_temporal_mass_matched=report['initial']['observables']['temporal_residual']['mass_matched'],
                                           conformal_wall_flux=report['initial']['observables']['temporal_residual']['flux']))
                    remember(0, first, kind='initial_anchor', residual=None)
                arriving = [i for i, cut in enumerate(cuts) if abs(float(cut)-xb) <= _ulp_gap(float(cut), xb)]
                queries = []
                if (accepted+1) % 8 == 0 or arriving:
                    queries.append((xm, 'interior_midpoint'))
                for cut_index, cut in enumerate(cuts):
                    cut = float(cut)
                    if xa < cut < xb:
                        queries.append((cut, 'interior_sample'))
                    elif abs(cut-xb) <= _ulp_gap(cut, xb):
                        queries.append((xb, 'sample_endpoint'))
                    elif accepted == 0 and abs(cut-xa) <= _ulp_gap(cut, xa):
                        queries.append((xa, 'sample_endpoint'))
                seen = set()
                midpoint_residual = None
                phase = 'continuous residual'
                for query, kind in queries:
                    key = float(query)
                    if key in seen:
                        continue
                    seen.add(key)
                    dense = fixed_time_residual(step, key, basis=basis, exterior_cells=exterior_cells,
                                               stage=stage, budget=budget)
                    residual = _residual_scalars(dense, basis=basis, budget=budget, background=background)
                    for name, value in residual.items():
                        maximum('time.'+name, value)
                    if kind == 'interior_midpoint':
                        midpoint_residual = residual
                    for cut_index, cut in enumerate(cuts):
                        if abs(float(cut)-key) <= _ulp_gap(float(cut), key):
                            remember(cut_index, dense['stage'], kind=kind, residual=residual,
                                     midpoint_residual=midpoint_residual if kind == 'interior_midpoint' else None)
                rates = [stage_observables(item, ORIGINAL_ALPHA) for item in step['stages']]
                # Diagnostic quadrature only. The accepted state is not replaced.
                transport.advance(ssprk3_diagnostic_increment(xb-xa, [item['matched_mass_rate'] for item in rates]),
                                  ssprk3_diagnostic_increment(xb-xa, [item['hybrid_c_rate'] for item in rates]),
                                  rates[3]['physical_mass'], rates[3]['c'])
                for item in step['stages']:
                    maximum('exterior_characteristic_courant', (xb-xa)*item['exterior_characteristic_courant_rate'])
                state = step['accepted_state']
                final_row = step['stages'][3]
                last_x = xb
                accepted += 1
                trapezoid.add(xb, dict(physical_temporal_mass_actual=rates[3]['temporal_residual']['mass_actual'],
                                       physical_temporal_mass_matched=rates[3]['temporal_residual']['mass_matched'],
                                       conformal_wall_flux=rates[3]['temporal_residual']['flux']))
                for cut_index in arriving:
                    remember(cut_index, final_row, kind='accepted_endpoint', residual=None,
                             midpoint_residual=midpoint_residual)
                history.update(digest_state(dict(x=xb, state=state)).encode())
                if accepted % 16 == 0:
                    print(json.dumps(dict(event='progress', accepted=accepted, x=xb,
                                          cpu_seconds=time.process_time()-budget.started,
                                          wall_radius=rates[3]['einstein_radius'],
                                          max_mass_gap=max(maxima.get('temporal.physical_temporal_mass_actual', 0.),
                                                           maxima.get('temporal.physical_temporal_mass_matched', 0.)),
                                          max_flux_gap=maxima.get('boundary.conformal_wall_flux', 0.))), flush=True)
                budget.check()
        if last_x != x_right:
            raise ValueError('march did not end at the source endpoint')
        if any(item is None for item in samples):
            raise ValueError('a common source cut was not recorded')
        phase = 'endpoint summary'
        final = stage_observables(final_row, ORIGINAL_ALPHA)
        report['endpoint'] = final
        report['samples'] = samples
        report['label_margin']['evolved_minimum_upper'] = minimum_label_margin
        report['integrated_defect'] = dict(
            quadrature='composite trapezoid on accepted step endpoints',
            fields=['physical_temporal_mass_actual', 'physical_temporal_mass_matched', 'conformal_wall_flux'],
            signed=trapezoid.signed, absolute=trapezoid.absolute,
            weak_balance_included=False,
            note='sampled defect integral, not an error enclosure')
        report['transport_comparison'] = dict(
            quadrature='SSPRK3 weights 1/6, 1/6, 2/3 on left, right-predictor and midpoint stages; accepted endpoint has zero weight',
            mass_definition='M_E=R_E**3/2+alpha*normalized_mass',
            matched_rate='1.5*R^2*r_x+alpha*(moving_density+T_x*A^2*exp(lapse)/2*R^2*Pi*Phi)',
            scalar='physical shared c compared with the integral of the hybrid wall c rate',
            state_corrected=False, **transport.last,
            max_abs_running_mass_gap=transport.max_abs_mass_gap,
            max_abs_running_scalar_gap=transport.max_abs_scalar_gap)
        report['endpoint_state_sha256'] = digest_state(state)
        budget.check()
        report['status'] = 'completed_numerical_endpoint'
    except (CoupledEvolutionError, ValueError, RuntimeError, FloatingPointError, OverflowError,
            KeyError, TypeError) as error:
        report['status'] = 'stopped'
        diagnostic = getattr(error, 'diagnostic', None)
        report['failure'] = dict(phase=diagnostic.get('phase', phase) if isinstance(diagnostic, dict) else phase,
                                 reason=str(error),
                                 stage_calls=None if not isinstance(diagnostic, dict) else diagnostic.get('stage_calls'),
                                 completed_stages=None if not isinstance(diagnostic, dict) else diagnostic.get('completed_stages'))
    if samples:
        report['samples'] = [item for item in samples if item is not None]
    report.update(maximum_absolute=maxima, minimum_stage_metric_certificate=minimum_A,
                  accepted_history_sha256=history.hexdigest(),
                  execution=dict(charged_values=budget.evaluations,
                                 cpu_seconds=time.process_time()-budget.started,
                                 common_stage_calls=calls, accepted_steps=accepted, last_accepted_x=last_x,
                                 source_integrations=report.get('preparation', {}).get('source_integrations_completed', 0),
                                 source_integrations_attempted=report.get('preparation', {}).get('source_integrations_attempted', 0),
                                 retries=0, raw_history_saved=False, final_state_saved=False))
    report = plain(report)
    assert_summary(report)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+'\n',
                           encoding='utf-8', newline='\n')
    print(json.dumps(dict(status=report['status'], execution=report['execution'],
                          failure=report.get('failure'))), flush=True)
    return 0 if report['status'] == 'completed_numerical_endpoint' else 1


if __name__ == '__main__':
    raise SystemExit(main())
