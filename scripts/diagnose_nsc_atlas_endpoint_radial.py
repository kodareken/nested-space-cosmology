#!/usr/bin/env python3
"""Check radial quadrature at one authenticated failed stage; no time step.

``--independent-mass-lapse`` integrates the full geometric mass, lapse and
fixed-r mass rate at that stage. It does not take another time step.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from numpy.polynomial.legendre import leggauss, legvander

import run_nsc_atlas_fixed_endpoint as endpoint
from recursive_horizons import nsc_atlas_cap_fixed_reference as fixed
from recursive_horizons import nsc_atlas_cap_variational_reference as radial
from recursive_horizons.nsc_atlas_cap_fixed_basis import exact_number
from recursive_horizons.nsc_atlas_cap_panel_jets import history_piece_indices

CASES=((32,2e-13,2e-15),(64,2e-13,2e-15),(128,2e-13,2e-15),
       (64,2e-14,1e-15))
CONFIRMATION_CASES=((32,2.5e-14,1e-15),(64,2.5e-14,1e-15))
PANEL_CASES=((32,2.5e-14,1e-15,.5),(32,2.5e-14,1e-15,.25),
             (64,2.5e-14,1e-15,.5))
ONE_SIDED_CASES=((32,2.5e-14,1e-15,None),(32,2.5e-14,1e-15,.25),
                 (64,2.5e-14,1e-15,None))


def inward_panel_point(point,span):
    """Numerical one-sided limits for a split radial integration panel.

    This diagnostic changes only points within 64 coordinate ulps of an
    endpoint. Refuse a panel too small to resolve that displacement.
    """
    lo,hi=map(float,span)
    margin=64*np.spacing(max(1.,abs(lo),abs(hi)))
    if not lo < hi or margin > 1e-6*(hi-lo):
        raise ValueError('one-sided radial endpoint is not resolved')
    if point <= lo+margin:
        return lo+margin
    if point >= hi-margin:
        return hi-margin
    return point


def panel_step_options(span, options, fraction):
    """Bound the forward step inside each existing smooth radial panel."""
    options=dict(options)
    if fraction is not None:
        if not np.isfinite(fraction) or not 0 < fraction <= 1:
            raise ValueError('radial panel fraction must be in (0, 1]')
        maximum=(float(span[1])-float(span[0]))*fraction
        if not np.isfinite(maximum) or maximum <= 0:
            raise ValueError('positive radial panel required')
        options['max_step']=min(options.get('max_step',np.inf),maximum)
        if 'first_step' in options:
            options['first_step']=min(options['first_step'],options['max_step'])
    return options


_WALL_TRACES=('f','p_trace','h_trace','q_trace','velocity','A','B')


def _safe_ratio(difference, reference):
    if reference==0.:
        return None
    return difference/abs(reference)


def physical_field_comparison(basis,reference,other,budget):
    """Galerkin L2 of one frozen stage's cap rates, plus wall-trace differences.

    ``squared_norm`` integrates the reconstructed field against the model's
    physical-z coordinate measure dz, not proper volume or energy. Diagonal Gram scaling
    stays inside that integral, so a coefficient infinity norm is not reused
    as a physical difference. This reports one stage. It does not decide
    equivalence and it does not replace the coordinate-rate gate.
    """
    left,right=reference['connection']['cap'],other['connection']['cap']
    gram=basis.gram(left['wall_coordinate'],budget=budget)
    def galerkin(coefficients):
        return math.sqrt(gram.squared_norm(np.asarray(coefficients,dtype=float),budget=budget))
    fields={}
    for key in ('chi','canonical_density'):
        original=np.asarray(reference['rates']['cap'][key],dtype=float)
        delta=original-np.asarray(other['rates']['cap'][key],dtype=float)
        reference_norm=galerkin(original)
        difference=galerkin(delta)
        fields[key]=dict(reference_galerkin_L2=reference_norm,difference_galerkin_L2=difference,
                         difference_over_reference=_safe_ratio(difference,reference_norm),
                         coordinate_infinity=float(np.max(np.abs(delta))))
    traces={}
    for key in _WALL_TRACES:
        value=float(left['wall'][key]);delta=abs(value-float(right['wall'][key]))
        traces[key]=dict(reference=value,absolute_difference=delta,
                         difference_over_reference=_safe_ratio(delta,value))
    flux=None
    if 'residual' in reference['connection'] and 'residual' in other['connection']:
        flux=abs(float(reference['connection']['residual']['conformal_wall_flux'])
                 -float(other['connection']['residual']['conformal_wall_flux']))
    return dict(gram_condition=float(gram.condition),
                wall_coordinate_absolute_difference=abs(float(left['wall_coordinate'])
                                                        -float(right['wall_coordinate'])),
                cap_rate_fields=fields,wall_traces=traces,conformal_wall_flux_absolute_difference=flux,
                scope='one frozen stage; Galerkin L2 of reconstructed cap rates on the reference wall Gram; not trajectory equivalence and not a coordinate-rate gate')


# Declared before any measured stage is accepted or rejected. The 32/64 gap is
# not one of these tolerances and is not an error bound.
FULL_MASS_ABSOLUTE_TOLERANCE=1e-10
LAPSE_ABSOLUTE_TOLERANCE=1e-10
MASS_RATE_ABSOLUTE_TOLERANCE=1e-8
INDEPENDENT_ORDERS=(32, 64)
_COMPARISON_TOLERANCES=(
    ('full_mass', FULL_MASS_ABSOLUTE_TOLERANCE),
    ('lapse', LAPSE_ABSOLUTE_TOLERANCE),
    ('fixed_radius_mass_rate', MASS_RATE_ABSOLUTE_TOLERANCE),
    ('moving_wall_mass_rate', MASS_RATE_ABSOLUTE_TOLERANCE),
    ('temporal_mass_residual', MASS_RATE_ABSOLUTE_TOLERANCE))
_GAUSS={}
INDEPENDENT_SCOPE=(
    'One authenticated frozen stage. v is the full-mass derivative at fixed areal radius: '
    'polynomial cap-coefficient rates plus background T partials times the coupled T_x, '
    'with no moving-R field coordinate term. The moving-wall rate is the separate quantity '
    'v+m_r*R_x. The temporal residual uses the actual strong-jet channel. '
    'Orders 32 and 64 are a consistency measurement, not an error bound. '
    'These tolerances do not certify continuum closure, the unsaved all-stage maximum, '
    'the trajectory, Step 2, or a parent Lambda.')


def geometric_mass_sources(r, f, p, h, df, dp, dh, alpha):
    """Full-mass densities a, s and their fixed-r directional derivatives.

    r is the areal radius. The derivatives df, dp, dh are taken at fixed r.
    Vacuum curvature sits in s; alpha multiplies only the scalar terms.
    """
    r=float(r); alpha=float(alpha)
    values=(r, float(f), float(p), float(h), float(df), float(dp), float(dh), alpha)
    if alpha<0. or r<0. or not all(math.isfinite(value) for value in values):
        raise ValueError('nonfinite geometric mass source')
    f, p, h, df, dp, dh=values[1:7]
    scale=math.sqrt(alpha/6.)
    if alpha>0. and scale==0.:
        raise ValueError('alpha scale underflows')
    if r==0.:
        if p!=0. or h!=0. or dp!=0. or dh!=0.:
            raise ValueError('mass source is singular at the centre')
        return dict(a=0., s=0., da=0., ds=0.)
    q=scale*f
    try:
        with np.errstate(over='raise', invalid='raise', divide='raise'):
            cosine, sine=float(np.cosh(q)), float(np.sinh(q))
            square=p*p+h*h
            kinetic=p*dp+h*dh
            sources=dict(a=alpha*square/(2.*r), da=alpha*kinetic/r,
                         s=1.5*r*r*cosine**4+alpha*square/4.,
                         ds=6.*r*r*scale*cosine**3*sine*df+alpha*kinetic/2.)
    except (FloatingPointError, OverflowError) as error:
        raise ValueError('nonfinite geometric mass source') from error
    if not all(math.isfinite(value) for value in sources.values()):
        raise ValueError('nonfinite geometric mass source')
    return sources


def centre_lapse_datum(f_axis, alpha):
    """Inherited centre clock, -log(cosh(sqrt(alpha/6)*f_axis))."""
    alpha=float(alpha); f_axis=float(f_axis)
    if alpha<0. or not math.isfinite(alpha) or not math.isfinite(f_axis):
        raise ValueError('nonfinite centre lapse datum')
    scale=math.sqrt(alpha/6.)
    if alpha>0. and scale==0.:
        raise ValueError('alpha scale underflows')
    q=scale*f_axis
    if q==0.:
        return 0.
    return -float(np.log1p(2.*np.sinh(q/2.)**2))


def _legendre_prefix(xi, degree):
    """Exact integrals of P_0..P_{degree-1} from -1 to each xi."""
    xi=np.asarray(xi, dtype=float)
    basis=legvander(xi, degree)
    prefix=np.empty((len(xi), degree))
    prefix[:, 0]=xi+1.
    if degree>1:
        degree_index=np.arange(1, degree)
        prefix[:, 1:]=(basis[:, 2:degree+1]-basis[:, :degree-1])/(2.*degree_index+1.)
    return prefix


def _gauss_rule(order):
    if isinstance(order, bool) or not isinstance(order, (int, np.integer)) or int(order) not in INDEPENDENT_ORDERS:
        raise ValueError('independent mass quadrature order must be 32 or 64')
    order=int(order)
    cached=_GAUSS.get(order)
    if cached is None:
        nodes, weights=leggauss(order)
        vandermonde=legvander(nodes, order-1)
        scale=(2.*np.arange(order)+1.)/2.
        prefix_nodes=_legendre_prefix(nodes, order)
        prefix_end=_legendre_prefix(np.array([1.]), order)[0]
        cached=(nodes, weights, vandermonde, scale, prefix_nodes, prefix_end)
        _GAUSS[order]=cached
    return cached


def _projected_prefix(vandermonde, weights, scale, prefix_nodes, prefix_end, values):
    coefficients=scale*(vandermonde.T@(weights*np.asarray(values, dtype=float)))
    return prefix_nodes@coefficients, float(prefix_end@coefficients)


def integrate_independent_mass(spans, coefficients, order, *, entering_mass=0., entering_rate=0.,
                               entering_lapse=0., budget=None):
    """Integrate m_r+a*m=s and v_r+a*v=ds-da*m on the caller's panels.

    ``coefficients(panel, coordinate)`` returns a, s, da, ds as derivatives
    with respect to areal radius, and ``jacobian`` = dr/d(panel coordinate).
    Those densities are multiplied by the Jacobian here. The integrator does
    not call the reference radial ODE. A budget, when supplied, is charged
    one site per coefficient sample; the callback charges its own field work.
    """
    nodes, weights, vandermonde, scale, prefix_nodes, prefix_end=_gauss_rule(order)
    mass=float(entering_mass); rate=float(entering_rate); lapse=float(entering_lapse)
    if not all(math.isfinite(value) for value in (mass, rate, lapse)):
        raise ValueError('nonfinite entering mass, rate, or lapse')
    spans=tuple(spans)
    if not spans:
        raise ValueError('at least one integration panel is required')
    samples=0
    for panel, span in enumerate(spans):
        if len(span)!=2:
            raise ValueError('each mass panel needs two endpoints')
        left, right=float(span[0]), float(span[1])
        if not (math.isfinite(left) and math.isfinite(right) and right>left):
            raise ValueError('positive radial panel required')
        half=(right-left)/2.
        incoming_mass, incoming_rate=mass, rate
        a_xi=np.empty(len(nodes)); s_xi=np.empty(len(nodes))
        da_xi=np.empty(len(nodes)); ds_xi=np.empty(len(nodes))
        for index, coordinate in enumerate(left+half*(nodes+1.)):
            row=coefficients(panel, float(coordinate))
            samples+=1
            if budget is not None:
                budget.check(1)
            jacobian=float(row['jacobian'])
            sample=tuple(float(row[key]) for key in ('a', 's', 'da', 'ds'))
            if not jacobian>0. or not all(math.isfinite(value) for value in sample):
                raise ValueError('nonfinite panel mass coefficient')
            measure=jacobian*half
            a_xi[index], s_xi[index], da_xi[index], ds_xi[index]=(value*measure for value in sample)
        factor_coefficients=scale*(vandermonde.T@(weights*a_xi))
        integrating_nodes=prefix_nodes@factor_coefficients
        integrating_end=float(prefix_end@factor_coefficients)
        try:
            with np.errstate(over='raise', invalid='raise', divide='raise'):
                weight_nodes=np.exp(integrating_nodes)
                mass_nodes, mass_end=_projected_prefix(
                    vandermonde, weights, scale, prefix_nodes, prefix_end, weight_nodes*s_xi)
                mass_nodes=np.exp(-integrating_nodes)*(incoming_mass+mass_nodes)
                mass=float(np.exp(-integrating_end)*(incoming_mass+mass_end))
                source_nodes=weight_nodes*(ds_xi-da_xi*mass_nodes)
                _, rate_end=_projected_prefix(
                    vandermonde, weights, scale, prefix_nodes, prefix_end, source_nodes)
                rate=float(np.exp(-integrating_end)*(incoming_rate+rate_end))
        except (FloatingPointError, OverflowError) as error:
            raise ValueError('panel integrating factor is not representable') from error
        lapse+=integrating_end
        if not all(math.isfinite(value) for value in (mass, rate, lapse)):
            raise ValueError('nonfinite independent mass or lapse')
    return dict(order=int(order), mass=float(mass), rate=float(rate), lapse=float(lapse), samples=samples)


def compare_independent_mass(reference, measurements, boundary):
    """Apply the predeclared absolute tolerances without relaxing them."""
    measurements=tuple(measurements)
    if tuple(int(item['order']) for item in measurements)!=INDEPENDENT_ORDERS:
        raise ValueError('independent mass comparison requires orders 32 and 64')
    compared=tuple(key for key, _tolerance in _COMPARISON_TOLERANCES)
    orders=[]; failed=[]
    for item in measurements:
        differences={}
        for key, tolerance in _COMPARISON_TOLERANCES:
            difference=float(item[key])-float(reference[key])
            differences[key]=difference
            if not math.isfinite(difference) or abs(difference)>tolerance:
                failed.append('order_%s.%s'%(int(item['order']), key))
        recorded=dict(order=int(item['order']), differences=differences)
        for key in ('samples', 'charged_values', 'metric_A', 'wall_radius', *compared):
            if key in item:
                recorded[key]=int(item[key]) if key=='samples' else float(item[key])
        orders.append(recorded)
    order_differences={key: float(measurements[1][key])-float(measurements[0][key]) for key in compared}
    return dict(
        method=('Gauss-Legendre 32 and 64 prefix integrals on the existing smooth panels; '
                'full geometric mass, including vacuum, and the fixed-r directional derivative; '
                'not the reference DOP853 mass ODE, not cap rho, and not stored gradients'),
        directional_derivative=(
            'v=m_x at fixed areal radius. Total field directions are the polynomial cap '
            'coefficient rates plus background f_T, p_T, h_T times coupled T_x. '
            'The fixed basis has no additional moving-R field term. '
            'Moving-wall mass rate is v+m_r*R_x. '
            'Temporal residual is v/T_x - alpha*A(R)^2*exp(delta(R))*p_wall*h_wall/2 '
            'on the actual strong-jet channel.'),
        channel='physical_temporal_mass_actual',
        reference={key: reference[key] for key in reference},
        orders=orders, order_differences=order_differences,
        tolerances=dict(
            full_mass_absolute=FULL_MASS_ABSOLUTE_TOLERANCE,
            lapse_absolute=LAPSE_ABSOLUTE_TOLERANCE,
            physical_directional_and_temporal_mass_rate_absolute=MASS_RATE_ABSOLUTE_TOLERANCE,
            order_difference='consistency measurement; not an error bound or a pass/fail tolerance'),
        outcome='failed' if failed else 'passed', failed_channels=failed,
        panels=int(boundary['panels']), centre_lapse=float(boundary['centre_lapse']),
        centre_field=float(boundary['centre_field']),
        wall=dict(p=float(boundary['wall_p']), h=float(boundary['wall_h']),
                  reference_p_trace=float(boundary['reference_p_trace']),
                  reference_h_trace=float(boundary['reference_h_trace']),
                  p_difference=float(boundary['wall_p'])-float(boundary['reference_p_trace']),
                  h_difference=float(boundary['wall_h'])-float(boundary['reference_h_trace'])),
        scope=INDEPENDENT_SCOPE)


def independent_stage_mass(basis, state_chi, state_density, stage, background, budget):
    """One reference stage's panels, integrated independently at orders 32 and 64."""
    cap=stage['connection']['cap']
    if 'background' not in cap:
        raise ValueError('independent mass check requires the retained background partition')
    alpha=float(cap['alpha']); radius=float(cap['radius'])
    clock_rate=float(stage['wall_rates']['T']); radius_rate=float(stage['wall_rates']['r'])
    if not clock_rate>0. or not radius>0.:
        raise ValueError('coupled clock rate and radius must be positive')
    if cap['adjoint_order']!=32:
        raise ValueError('reference stage adjoint order is not 32')
    wall_coordinate=exact_number(cap['wall_coordinate'], 'wall coordinate')
    background_time=float(cap['background']['time'])
    fields=fixed._Fields(basis, state_chi, state_density, budget)
    coefficient_rates=stage['rates']['cap']
    rate_fields=fixed._Fields(basis, coefficient_rates['chi'], coefficient_rates['canonical_density'], budget)
    breaks, owners=fixed._background_partition(basis, wall_coordinate, background_time, background)
    if float(breaks[0])!=0. or int(owners[-1])!=basis.cells-1:
        raise ValueError('mass integration must run from the centre through the wall cell')
    cut=float((wall_coordinate-basis.exact_knots[-2])/(basis.exact_knots[-1]-basis.exact_knots[-2]))
    panel_indices=[]
    for panel, cell in enumerate(owners):
        cell=int(cell)
        midpoint=float(breaks[panel]+(breaks[panel+1]-breaks[panel])/2.)
        z=fields.left[cell]+fields.width[cell]*(midpoint-cell)
        panel_indices.append(history_piece_indices(
            background.profile, [background_time-z, background_time+z], budget=budget))

    def coefficients(panel, coordinate):
        cell=int(owners[panel])
        left, right=float(breaks[panel]), float(breaks[panel+1])
        point=min(right, max(left, float(coordinate)))
        local=float(cut if panel==len(owners)-1 and point==float(breaks[-1]) else point-cell)
        departure=fields.evaluate(cell, local, budget)
        direction=rate_fields.evaluate(cell, local, budget)
        jets=background.reduced_jets(background_time, departure['z'], budget=budget,
                                     polynomial_indices=panel_indices[panel])
        total_f=departure['f']+jets['f']
        total_p=departure['p']+jets['p']
        total_h=departure['h']+jets['h']
        # Fixed-r direction: coefficient rates already differentiate the departure.
        # Background motion enters only through the coupled clock, not through R_x.
        df=direction['f']+jets['f_T']*clock_rate
        dp=direction['p']+jets['p_T']*clock_rate
        dh=direction['h']+jets['h_T']*clock_rate
        sources=geometric_mass_sources(departure['r'], total_f, total_p, total_h, df, dp, dh, alpha)
        sources.update(jacobian=fields.width[cell]*departure['A0'], r=float(departure['r']),
                       p=float(total_p), h=float(total_h))
        return sources

    axis=fields.evaluate(0, 0., budget)
    axis_jets=background.reduced_jets(background_time, 0., budget=budget)
    centre_field=float(axis['f']+axis_jets['f'])
    datum=centre_lapse_datum(centre_field, alpha)
    wall=coefficients(len(owners)-1, float(breaks[-1]))
    spans=tuple((float(left), float(right)) for left, right in zip(breaks[:-1], breaks[1:]))
    measurements=[]
    for order in INDEPENDENT_ORDERS:
        charged_before=budget.evaluations
        integrated=integrate_independent_mass(
            spans, coefficients, order, entering_mass=0., entering_rate=0.,
            entering_lapse=datum, budget=budget)
        mass=integrated['mass']; rate=integrated['rate']; lapse=integrated['lapse']
        wall_radius=float(wall['r'])
        metric_A=1.-2.*mass/wall_radius
        radial_mass=wall['s']-wall['a']*mass
        temporal=rate/clock_rate-alpha*metric_A**2*math.exp(lapse)*wall['p']*wall['h']/2.
        measurements.append(dict(
            order=order, samples=integrated['samples'], charged_values=budget.evaluations-charged_before,
            full_mass=mass, lapse=lapse, metric_A=float(metric_A), wall_radius=wall_radius,
            fixed_radius_mass_rate=rate, moving_wall_mass_rate=rate+radial_mass*radius_rate,
            temporal_mass_residual=float(temporal)))
    vacuum_mass=.5*radius**3
    cap_mass=float(cap['mass'])
    reference=dict(
        adjoint_order=32, rtol=2.5e-14, atol=1e-15, background_panel_sides=True,
        compiled_coefficients=False, radius=radius, alpha=alpha, cap_mass=cap_mass,
        vacuum_mass=vacuum_mass, full_mass=vacuum_mass+alpha*cap_mass, lapse=float(cap['lapse']),
        T_x=clock_rate, radius_rate=radius_rate,
        fixed_radius_mass_rate=alpha*(float(cap['mass_rate'])-float(cap['moving_density'])),
        moving_wall_mass_rate=float(stage['connection']['physical_mass_rate']),
        temporal_mass_residual=float(stage['connection']['residual']['physical_temporal_mass_actual']),
        temporal_channel='physical_temporal_mass_actual')
    boundary=dict(panels=int(len(owners)), centre_lapse=datum, centre_field=centre_field,
                  wall_p=float(wall['p']), wall_h=float(wall['h']),
                  reference_p_trace=float(cap['wall']['p_trace']),
                  reference_h_trace=float(cap['wall']['h_trace']))
    return compare_independent_mass(reference, measurements, boundary)


def independent_mass_lapse_main(args):
    if args.output.exists() or not args.output.parent.is_dir():
        raise ValueError('choose a new output in an existing directory')
    raw=args.receipt.read_bytes(); receipt=json.loads(raw)
    if receipt['status']!='stopped':
        raise ValueError('require a stopped evolution receipt')
    state, x=read_seed(receipt)
    commit, bindings=endpoint.frozen_sources()
    diagnostic_path=Path(__file__).resolve()
    bindings[diagnostic_path.relative_to(endpoint.ROOT).as_posix()]=hashlib.sha256(
        diagnostic_path.read_bytes()).hexdigest()
    budget=endpoint.Budget(maximum_evaluations=40_000_000, cpu_seconds=60.)
    report=dict(schema='nsc-atlas-endpoint-radial-diagnostic-v1', status='incomplete',
        source_commit=commit, source_bindings=bindings,
        input_receipt=dict(path=args.receipt.resolve().relative_to(endpoint.ROOT).as_posix(),
                           bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()),
        stage_seed=receipt['failed_stage_seed'], x=x, level=receipt['level'],
        physical_time_steps=0, comparison_kind='independent_mass_lapse',
        limits=dict(cpu_seconds=60., charged_values=40_000_000, reference_stages=1,
                    quadrature_orders=list(INDEPENDENT_ORDERS), reset=False, retry=False,
                    scope='one reference coupled stage and two independent quadrature orders on one budget'))
    try:
        endpoint.authenticate(endpoint.ROOT, budget)
        source_raw=(endpoint.ROOT/endpoint.producer.BASE_INPUT).read_bytes()
        if hashlib.sha256(source_raw).hexdigest()!=receipt['input']['sha256']:
            raise ValueError('physical source changed')
        base=json.loads(source_raw); geometry=base['geometry']; data=base['initial_data']['finite']
        rows={key: np.asarray(value, dtype=float) for key, value in data['exact_rows'].items()}
        spec=endpoint.level_spec(receipt['level'])
        basis=endpoint.FixedCapBasis(endpoint.cap_knots(receipt['level']), budget=budget)
        background, offset=endpoint.prepare_history(budget, report.setdefault('preparation', {}),
                                                   enrich_ordinary_pieces=True)
        for key in ('input_sha256', 'original_input_sha256'):
            if report['preparation'][key]!=receipt['preparation'][key]:
                raise ValueError('prepared history changed')
        snapshot=endpoint.SourceSnapshot(geometry, alpha=data['alpha'], x_nodes=rows['x'],
                                         source_rows=rows, budget=budget)
        gauge=endpoint.ReferenceWallGauge(geometry, order=64, budget=budget)
        chi=np.array(state['cap']['chi'], dtype=float, copy=True)
        density=np.array(state['cap']['canonical_density'], dtype=float, copy=True)
        stage=endpoint.coupled_fixed_stage(
            state, basis=basis, x=x, eta=np.linspace(0., 1., spec['exterior_cells']+1),
            source=snapshot.sample(x), gauge=gauge.sample(x), alpha=data['alpha'], tube=endpoint.TUBE,
            budget=budget, background=background, background_clock_offset=offset,
            rtol=2.5e-14, atol=1e-15, background_panel_sides=True, compiled_coefficients=False)
        report['independent_mass_lapse']=independent_stage_mass(basis, chi, density, stage, background, budget)
        outcome=report['independent_mass_lapse']['outcome']
        if outcome=='passed':
            report['status']='completed_one_stage_independent_mass_lapse'
        else:
            report['status']='stopped'
            report['reason']='independent mass/lapse tolerances failed: '+', '.join(
                report['independent_mass_lapse']['failed_channels'])
    except (ValueError, RuntimeError, FloatingPointError, OverflowError) as error:
        report.update(status='stopped', reason=str(error))
    report['execution']=dict(cpu_seconds=time.process_time()-budget.started,
                             charged_values=budget.evaluations)
    report=endpoint.plain(report)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+'\n',
                          encoding='utf-8', newline='\n')
    summary=dict(status=report['status'], execution=report['execution'])
    if 'reason' in report:
        summary['reason']=report['reason']
    if 'independent_mass_lapse' in report:
        summary['independent_mass_lapse']=report['independent_mass_lapse']
    print(json.dumps(summary))
    return 0 if report['status']=='completed_one_stage_independent_mass_lapse' else 1


def read_seed(receipt,root=endpoint.ROOT):
    binding=receipt['failed_stage_seed']
    path=(root/binding['path']).resolve()
    if not path.is_relative_to((root/'.build/results').resolve()):
        raise ValueError('stage seed must stay under the ignored result directory')
    raw=path.read_bytes()
    if len(raw)!=binding['bytes'] or hashlib.sha256(raw).hexdigest()!=binding['sha256']:
        raise ValueError('stage seed bytes changed')
    state={}
    with np.load(path,allow_pickle=False) as saved:
        if (saved['level'].item()!=receipt['level'] or
            saved['source_commit'].item()!=receipt['source_commit'] or
            saved['x'].item()!=binding['source_x']):
            raise ValueError('stage seed metadata changed')
        x=float(saved['x'].item())
        for key in saved.files:
            if '.' in key:
                group,name=key.split('.',1)
                state.setdefault(group,{})[name]=saved[key].copy()
    return state,x


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--confirm-forward-accuracy',action='store_true',
                        help='the admitted 32/64-point check of the selected tighter tolerance')
    mode.add_argument('--panel-step-check',action='store_true',
                      help='half/quarter-panel forward steps at the selected tolerance')
    mode.add_argument('--one-sided-panel-check',action='store_true',
                      help='test inward endpoint limits at existing radial formula breaks')
    mode.add_argument('--confirm-panel-polynomials',action='store_true',
                      help='32/64-point check of explicit panel polynomial ownership at unchanged coordinates')
    mode.add_argument('--compare-compiled-coefficients',action='store_true',
                      help='one frozen stage: reference, first compiled call, warm compiled call; no time steps')
    mode.add_argument('--independent-mass-lapse',action='store_true',
                      help='one frozen stage: independent full mass, lapse and fixed-r mass rate; no time steps')
    args=parser.parse_args(argv)
    if args.independent_mass_lapse:
        return independent_mass_lapse_main(args)
    cases=((CONFIRMATION_CASES[0],)*3 if args.compare_compiled_coefficients else
           ONE_SIDED_CASES if args.one_sided_panel_check else
           PANEL_CASES if args.panel_step_check else
           CONFIRMATION_CASES if (args.confirm_forward_accuracy or args.confirm_panel_polynomials) else CASES)
    if args.output.exists() or not args.output.parent.is_dir():
        raise ValueError('choose a new output in an existing directory')
    raw=args.receipt.read_bytes();receipt=json.loads(raw)
    if receipt['status']!='stopped':
        raise ValueError('require a stopped evolution receipt')
    state,x=read_seed(receipt)
    commit,bindings=endpoint.frozen_sources()
    diagnostic_path=Path(__file__).resolve()
    bindings[diagnostic_path.relative_to(endpoint.ROOT).as_posix()]=hashlib.sha256(
        diagnostic_path.read_bytes()).hexdigest()
    budget=endpoint.Budget(maximum_evaluations=40_000_000,cpu_seconds=60.)
    report=dict(schema='nsc-atlas-endpoint-radial-diagnostic-v1',status='incomplete',
        source_commit=commit,source_bindings=bindings,
        input_receipt=dict(path=args.receipt.resolve().relative_to(endpoint.ROOT).as_posix(),
                           bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()),
        stage_seed=receipt['failed_stage_seed'],x=x,level=receipt['level'],
        physical_time_steps=0,gradient_tolerance=2e-8,
        comparison_kind=('compiled_coefficients' if args.compare_compiled_coefficients else
                         'explicit_panel_polynomials' if args.confirm_panel_polynomials else
                         'one_sided_panel_limits' if args.one_sided_panel_check else
                         'forward_panel_step' if args.panel_step_check else
                         'selected_forward_accuracy' if args.confirm_forward_accuracy else 'error_attribution'),
        limits=dict(cpu_seconds=60.,charged_values=40_000_000,cases=len(cases)),cases=[])
    original=fixed._radial_variations
    original_solve=radial.solve_ivp
    radial_results=[];stage_results=[]
    try:
        endpoint.authenticate(endpoint.ROOT,budget)
        source_raw=(endpoint.ROOT/endpoint.producer.BASE_INPUT).read_bytes()
        if hashlib.sha256(source_raw).hexdigest()!=receipt['input']['sha256']:
            raise ValueError('physical source changed')
        base=json.loads(source_raw);geometry=base['geometry'];data=base['initial_data']['finite']
        rows={k:np.asarray(v,dtype=float) for k,v in data['exact_rows'].items()}
        spec=endpoint.level_spec(receipt['level'])
        basis=endpoint.FixedCapBasis(endpoint.cap_knots(receipt['level']),budget=budget)
        background,offset=endpoint.prepare_history(budget,report.setdefault('preparation',{}),
                                                  enrich_ordinary_pieces=True)
        for key in ('input_sha256','original_input_sha256'):
            if report['preparation'][key]!=receipt['preparation'][key]:
                raise ValueError('prepared history changed')
        snapshot=endpoint.SourceSnapshot(geometry,alpha=data['alpha'],x_nodes=rows['x'],
                                         source_rows=rows,budget=budget)
        gauge=endpoint.ReferenceWallGauge(geometry,order=64,budget=budget)
        for case_index,specification in enumerate(cases):
            order,rtol,atol=specification[:3]
            fraction=specification[3] if len(specification)==4 else None
            case=dict(adjoint_order=order,rtol=rtol,atol=atol,
                      effective_rtol=max(rtol,100*np.finfo(float).eps))
            if fraction is not None:
                case['maximum_radial_step_panel_fraction']=fraction
            panel_sides=args.confirm_panel_polynomials or args.compare_compiled_coefficients
            compiled=args.compare_compiled_coefficients and case_index>0
            if panel_sides:
                case['background_panel_sides']=True
            if args.compare_compiled_coefficients:
                case['coefficient_backend']=('reference','compiled_first','compiled_warm')[case_index]
            captured=[];started=time.process_time();started_wall=time.perf_counter();before=budget.evaluations
            stage=None
            panels=[];shifts=[]
            def solve(fun,span,initial,**kw):
                def one_sided(point,values):
                    interior=inward_panel_point(point,span)
                    if interior!=point:
                        shifts.append(abs(interior-point))
                    return fun(interior,values)
                result=original_solve(one_sided if args.one_sided_panel_check else fun,
                    span,initial,**panel_step_options(span,kw,fraction))
                budget.check(result.y.size)
                panels.append(dict(nfev=result.nfev,accepted_steps=len(result.t)-1,
                    maximum_internal_gradient=float(np.max(np.abs(result.y[2:])))))
                return result
            def routed(*a,**kw):
                kw.update(adjoint_order=order)
                result=original(*a,**kw)
                captured.append(result)
                return result
            fixed._radial_variations=routed
            radial.solve_ivp=solve
            try:
                stage=endpoint.coupled_fixed_stage(state,basis=basis,x=x,
                    eta=np.linspace(0.,1.,spec['exterior_cells']+1),source=snapshot.sample(x),
                    gauge=gauge.sample(x),alpha=data['alpha'],tube=endpoint.TUBE,budget=budget,
                    background=background,background_clock_offset=offset,rtol=rtol,atol=atol,
                    background_panel_sides=panel_sides,compiled_coefficients=compiled)
                cap=stage['connection']['cap']
                case.update(status='passed_existing_gradient_guard',mass=cap['mass'],lapse=cap['lapse'],
                    gradient_discrepancy=cap['gradient_discrepancy'],
                    projected_gradient_discrepancy=cap['projected_gradient_discrepancy'],
                    field_sites=cap['field_sites'],
                    physical_gaps=endpoint.physical_gaps(stage['connection']['residual']))
            except (ValueError,RuntimeError,FloatingPointError,OverflowError) as error:
                case.update(status='stopped',reason=str(error))
                budget.check()
            finally:
                fixed._radial_variations=original
                radial.solve_ivp=original_solve
            case['forward_panels']=dict(count=len(panels),nfev=sum(p['nfev'] for p in panels),
                accepted_steps=sum(p['accepted_steps'] for p in panels),
                maximum_internal_gradient=max((p['maximum_internal_gradient'] for p in panels),default=0.))
            if args.one_sided_panel_check:
                case['endpoint_limits']=dict(shifted_samples=len(shifts),
                    maximum_coordinate_displacement=max(shifts,default=0.),
                    coordinate_ulps=64,largest_allowed_panel_fraction=1e-6)
            radial_results.append(captured[-1] if captured else None)
            stage_results.append(stage)
            case.update(cpu_seconds=time.process_time()-started,wall_seconds=time.perf_counter()-started_wall,
                        charged_values=budget.evaluations-before)
            report['cases'].append(case)
        report['pair_differences']=[]
        pairs=((0,1),(0,2)) if (args.panel_step_check or args.one_sided_panel_check or args.compare_compiled_coefficients) else ((0,1),(1,2),(1,3))
        for i,j in pairs:
            if j>=len(radial_results):
                continue
            a,b=radial_results[i],radial_results[j]
            if a is not None and b is not None:
                report['pair_differences'].append(dict(cases=[i,j],
                    **{key:float(np.max(np.abs(a[key]-b[key]))) for key in
                       ('mass','lapse','gradient','adjoint_gradient')}))
        if args.compare_compiled_coefficients:
            comparison=[]
            for j in (1,2):
                a,b=stage_results[0],stage_results[j]
                if a is None or b is None:
                    comparison.append(dict(cases=[0,j],passed=False,reason='a stage stopped'))
                    continue
                rates={group+'.'+key:float(np.max(np.abs(value-b['rates'][group][key])))
                       for group,values in a['rates'].items() for key,value in values.items()}
                passed=all(np.allclose(value,b['rates'][group][key],rtol=1e-11,atol=1e-10)
                           for group,values in a['rates'].items() for key,value in values.items())
                scalars={key:abs(a['connection']['cap'][key]-b['connection']['cap'][key])
                         for key in ('mass','lapse')}
                passed=passed and all(value<=1e-11 for value in scalars.values())
                item=dict(cases=[0,j],passed=bool(passed),maximum_rate_differences=rates,
                          scalar_differences=scalars,
                          coordinate_gate='strict coefficient infinity norm and mass/lapse; unchanged')
                try:
                    item['physical_field_comparison']=physical_field_comparison(basis,a,b,budget)
                except (ValueError,RuntimeError,FloatingPointError,OverflowError) as error:
                    item['physical_field_comparison']=dict(status='unavailable',reason=str(error))
                comparison.append(item)
            report['compiled_equivalence']=dict(rate_atol=1e-10,rate_rtol=1e-11,
                mass_lapse_atol=1e-11,comparisons=comparison,
                passed=all(item['passed'] for item in comparison),
                scope='one frozen stage; coordinate-rate failure stays visible; Galerkin fields are not a trajectory equivalence or a replacement gate')
        report['status']='completed_one_stage_comparison'
    except (ValueError,RuntimeError,FloatingPointError,OverflowError) as error:
        report.update(status='stopped',reason=str(error))
    finally:
        fixed._radial_variations=original
        radial.solve_ivp=original_solve
    report['execution']=dict(cpu_seconds=time.process_time()-budget.started,
                              charged_values=budget.evaluations)
    report=endpoint.plain(report)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n',
                          encoding='utf-8',newline='\n')
    print(json.dumps(dict(status=report['status'],cases=report['cases'],execution=report['execution'])))
    return 0 if report['status']=='completed_one_stage_comparison' else 1


if __name__=='__main__':
    raise SystemExit(main())
