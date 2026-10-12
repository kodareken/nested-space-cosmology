"""Supplied field/preparation tangents on the existing canonical Dirac grid.

This adapter differentiates a declared numerical family, not a physical parent.
Initial columns, their tangents, incident forcing and its tangents are explicit.
The old propagator, metric generator and authenticated scientific owners are
unchanged. Time-dependent coefficients/forcing use the midpoint prescription.
"""
import numpy as np
from scipy.sparse import bmat, csr_matrix
from scipy.sparse.linalg import expm_multiply

from .nsc_transmitting_history_jets import generator_direction


def _finite_complex(value, name):
    if value is None:
        raise ValueError(f'explicit {name} required; missing data are not zero')
    result = np.asarray(value, dtype=complex)
    if not np.isfinite(result).all():
        raise ValueError(f'finite {name} required')
    return result


def _finite_real(value, name):
    if value is None:
        raise ValueError(f'explicit {name} required')
    result = np.asarray(value)
    if np.iscomplexobj(result) and np.any(result.imag != 0):
        raise ValueError(f'real {name} required')
    result = np.asarray(result.real, dtype=float)
    if not np.isfinite(result).all():
        raise ValueError(f'finite {name} required')
    return result


def restriction_covariance_tangent(F, dF, C, dC):
    """Return dF C F† + F dC F† + F C dF† without diagonal projection.

    F has trailing shape (restriction_rows, source_columns); C is a square
    source matrix. Tangents have the same trailing dimensions and may carry
    broadcast-compatible leading direction/frequency axes. dC is mandatory,
    including when its explicitly declared value is zero. No CAR, coisometry,
    source-complement allocation or physical-parent provenance is inferred.
    """
    F, dF, C, dC = (_finite_complex(value, name) for value, name in
                    ((F, 'restriction columns'), (dF, 'restriction tangents'),
                     (C, 'source covariance'), (dC, 'source covariance tangents')))
    if F.ndim < 2 or min(F.shape[-2:]) < 1:
        raise ValueError('nonempty restriction matrix dimensions required')
    sources = F.shape[-1]
    if (dF.ndim < 2 or dF.shape[-2:] != F.shape[-2:]
            or C.ndim < 2 or C.shape[-2:] != (sources, sources)
            or dC.ndim < 2 or dC.shape[-2:] != (sources, sources)):
        raise ValueError('matching restriction and source covariance matrix dimensions required')
    try:
        np.broadcast_shapes(F.shape[:-2], dF.shape[:-2], C.shape[:-2], dC.shape[:-2])
    except ValueError as error:
        raise ValueError('compatible leading restriction/tangent axes required') from error
    adjoint = lambda value: value.swapaxes(-1, -2).conj()
    return dF @ C @ adjoint(F) + F @ dC @ adjoint(F) + F @ C @ adjoint(dF)


def evolve_prepared_field_jets(owner, times, provider, initial_fields,
                               initial_tangents, incoming, *, sample_rows=None):
    """Evolve Phi'=L Phi+f and Y'=L Y+dL Phi+df on caller time nodes.

    initial_fields: complex characteristic columns (2*len(owner.x), sources).
    initial_tangents: explicit (ndir, rows, sources), including declared zeros.
    incoming(t): returns (forcing, forcing_tangents) in those same shapes.
    Both are already SAT-weighted and supported only on right inflow rows
    len(x)-1 and 2*len(x)-1; interior/left forcing is rejected.
    provider.values(t,x): raw PG metric (N,beta,q,r), shape (4,len(x)).
    provider.log_directions(t,x,metric): (ndir,4,len(x)) in
    (delta logN, delta beta, delta logq, delta logr).

    The existing owner enforces the reference metric at computational edges;
    generator_direction enforces unchanged metric edge tangents/SAT. Incident
    preparation may change through supplied f and df. No source frequency law,
    zero initial tangent or zero incident tangent is silently supplied.
    Initial columns live on the caller initial PG-time slice; they are not
    identified with the incoming rho1 KS covariance. Arbitrary finite columns
    do not define a full continuum CTP unitary.

    Optional sample_rows are unique integer characteristic row indices. Only
    those rows are retained at the initial node and after each completed step;
    no interpolation, stationary phase or Cauchy/preparation map is implied.
    """
    times = _finite_real(times, 'ordered caller time coordinates')
    if times.ndim != 1 or len(times) < 2 or np.any(np.diff(times) <= 0):
        raise ValueError('at least two strictly increasing caller time coordinates required')
    if not callable(getattr(provider, 'values', None)) or not callable(getattr(provider, 'log_directions', None)):
        raise ValueError('metric provider values and explicit log_directions required')
    if not callable(incoming):
        raise ValueError('explicit incoming forcing/tangent callable required')
    x = _finite_real(owner.x, 'owner spatial grid')
    if x.ndim != 1 or len(x) < 1:
        raise ValueError('one dimensional owner spatial grid required')
    rows = 2 * len(x)
    selected_rows = None
    if sample_rows is not None:
        selected_rows = np.asarray(sample_rows)
        if (selected_rows.ndim != 1 or selected_rows.size == 0
                or selected_rows.dtype.kind not in 'iu'
                or np.any(selected_rows < 0) or np.any(selected_rows >= rows)
                or len(np.unique(selected_rows)) != len(selected_rows)):
            raise ValueError('sample_rows must be nonempty unique integer indices within the characteristic rows')
        selected_rows = selected_rows.astype(np.intp, copy=True)
    phi = _finite_complex(initial_fields, 'initial fields')
    tangents = _finite_complex(initial_tangents, 'initial preparation tangents')
    if phi.ndim != 2 or phi.shape[0] != rows or phi.shape[1] < 1:
        raise ValueError('initial fields must have shape (2*len(x), nonzero source_count)')
    if tangents.ndim != 3 or tangents.shape[1:] != phi.shape:
        raise ValueError('initial tangents must have shape (ndir, rows, sources)')
    ndir, sources = tangents.shape[0], phi.shape[1]
    identity = np.eye(sources, dtype=complex)
    # The constant bottom block implements the affine forcing. It is a solver
    # augmentation, never a new physical state or term in the action.
    state = np.vstack((phi, *(tangents[j] for j in range(ndir)), identity))
    sampled_fields = [] if selected_rows is None else [phi[selected_rows].copy()]
    sampled_tangents = [] if selected_rows is None else [tangents[:, selected_rows].copy()]
    max_flux = max_dflux = max_speed = max_boundary = 0.
    max_auxiliary = 0.
    edge_rows = [0, len(x)-1, len(x), rows-1]
    outside_inflow = np.ones(rows, dtype=bool)
    outside_inflow[[len(x)-1, rows-1]] = False
    for left, right in zip(times[:-1], times[1:], strict=True):
        midpoint = (left + right) / 2
        metric = _finite_real(provider.values(midpoint, x), 'midpoint metric')
        if metric.shape != (4, len(x)):
            raise ValueError('metric must have shape (4, len(x))')
        L, checks = owner.generator(metric)
        directions = _finite_real(provider.log_directions(midpoint, x, metric), 'metric log directions')
        if directions.shape != (ndir, 4, len(x)):
            raise ValueError('metric direction count/shape must match initial tangents')
        supplied = incoming(midpoint)
        if not isinstance(supplied, (tuple, list)) or len(supplied) != 2:
            raise ValueError('incoming(time) must explicitly return forcing and forcing tangents')
        forcing = _finite_complex(supplied[0], 'incoming forcing')
        df = _finite_complex(supplied[1], 'incoming forcing tangents')
        if forcing.shape != phi.shape or df.shape != tangents.shape:
            raise ValueError('incoming forcing/tangent shapes must match initial columns')
        if np.any(forcing[outside_inflow] != 0) or np.any(df[:, outside_inflow] != 0):
            raise ValueError('incoming forcing and tangents must be supported only on the two right SAT inflow rows')
        if L.shape != (rows, rows):
            raise ValueError('owner generator dimension differs from characteristic columns')
        blocks = [[None] * (ndir + 2) for _ in range(ndir + 2)]
        blocks[0][0], blocks[0][-1] = L, csr_matrix(forcing)
        for index, direction in enumerate(directions):
            dL, error = generator_direction(owner, metric, direction)
            blocks[index+1][0], blocks[index+1][index+1] = dL, L
            blocks[index+1][-1] = csr_matrix(df[index])
            max_dflux = max(max_dflux, error)
        blocks[-1][-1] = csr_matrix((sources, sources), dtype=complex)
        augmented = bmat(blocks, format='csr')
        dt = right-left
        state = expm_multiply(dt*augmented, state, traceA=dt*augmented.diagonal().sum())
        if not np.isfinite(state).all():
            raise ArithmeticError('nonfinite prepared-field evolution')
        if selected_rows is not None:
            sampled_fields.append(state[selected_rows].copy())
            step_tangents = state[rows:(ndir+1)*rows].reshape(ndir, rows, sources)
            sampled_tangents.append(step_tangents[:, selected_rows].copy())
        max_auxiliary = max(max_auxiliary, float(np.max(abs(state[-sources:]-identity))))
        max_flux = max(max_flux, checks['norm_flux_algebra'])
        max_speed = max(max_speed, checks['maximum_absolute_speed'])
        if ndir:
            current = state[rows:(ndir+1)*rows].reshape(ndir, rows, sources)
            max_boundary = max(max_boundary, float(np.max(abs(current[:, edge_rows]))))
    result = {
        'evolved_field': state[:rows].copy(),
        'tangent_fields': state[rows:(ndir+1)*rows].reshape(ndir, rows, sources).copy(),
        'norm_flux_residual': float(max_flux),
        'metric_flux_tangent_residual': float(max_dflux),
        'auxiliary_identity_residual': float(max_auxiliary),
        'tangent_boundary_trace': float(max_boundary),
        'maximum_sampled_characteristic_speed': float(max_speed),
        'midpoint_steps': len(times)-1,
        'maximum_time_step': float(np.max(np.diff(times))),
        'direction_count': ndir, 'source_column_count': sources,
        'initial_and_incoming_tangents_explicit': True,
        'incoming_forcing_convention': 'already SAT-weighted, only rows n-1 and 2n-1; fixed boundary metric and zero dSAT',
        'initial_column_domain': 'caller initial PG-time slice, not an identified rho1 KS restriction or physical C0',
        'temporal_prescription': 'piecewise midpoint metric, metric tangent, forcing and forcing tangent; exponential augmented step',
        'continuum_error_bound': None, 'sampled_causal_buffer': None,
        'physical_preparation_status': 'OPEN', 'physical_history_status': 'OPEN',
        'full_spectral_covariance': None, 'physical_endpoint_family': None,
        'full_CTP_closure': False, 'full_branch_unitarity_claimed': False,
        'metric_evolution': False,
    }
    if selected_rows is not None:
        result.update(sample_times=times.copy(), sample_rows=selected_rows,
                      sampled_fields=np.stack(sampled_fields),
                      sampled_tangent_fields=np.stack(sampled_tangents).transpose(1, 0, 2, 3))
    return result
