"""Cross-platform replay controls without retained field evolution."""
import numpy as np
import pytest
from recursive_horizons import nsc_atlas_characteristic_initialization as model
from recursive_horizons.nsc_atlas_field_preparation import radius_product
from recursive_horizons.nsc_atlas_validation import source_jet_roundoff_tolerances
from recursive_horizons.nsc_atlas_validation import source_snapshot_roundoff_tolerances


def test_reference_endpoint_replay_preserves_interior_and_derivative():
    from scipy.interpolate import CubicHermiteSpline
    from recursive_horizons.nsc_atlas_validation import EndpointReplaySpline
    x=np.array([-2.5,-2.4,-2.3])
    spline=CubicHermiteSpline(x,x*x,2*x,extrapolate=False)
    replay=EndpointReplaySpline(spline)
    outside=np.array([np.nextafter(x[0],-np.inf),np.nextafter(x[-1],np.inf)])
    assert np.isnan(spline(outside)).all()
    np.testing.assert_array_equal(replay(outside),spline(x[[0,-1]]))
    np.testing.assert_array_equal(replay(outside,1),spline(x[[0,-1]],1))
    interior=np.linspace(x[0],x[-1],11)
    np.testing.assert_array_equal(replay(interior),spline(interior))
    assert replay.x is spline.x and replay.c is spline.c
    assert replay(float(outside[0]))==spline(x[0])


@pytest.mark.parametrize('point',[-2.5-1e-12,-2.3+1e-12,np.nan,np.inf])
def test_reference_endpoint_replay_rejects_larger_excursions(point):
    from scipy.interpolate import CubicHermiteSpline
    from recursive_horizons.nsc_atlas_validation import EndpointReplaySpline
    x=np.array([-2.5,-2.3])
    replay=EndpointReplaySpline(CubicHermiteSpline(x,x*x,2*x,extrapolate=False))
    with pytest.raises(ValueError,match='rounded spline endpoints'): replay(point)


@pytest.mark.parametrize('nonzero', [False, True])
def test_departure_replay_covers_reconstructed_wall_cancellation(monkeypatch, nonzero):
    import recursive_horizons.nsc_atlas_departure_diagnostics as diagnostic
    from recursive_horizons.nsc_atlas_validation import departure_history_roundoff_tolerances
    from test_nsc_atlas_departure_diagnostics import history, diagnose
    data, _, _ = history(nt=9, ne=17, nonzero=nonzero, quadratic=True)
    budget = model.Budget(maximum_evaluations=4000, cpu_seconds=5.)
    original = diagnose(data, budget)
    evaluator = diagnostic.vacuum_slice
    def shifted(U, x, **kwargs):
        row = evaluator(U, x, **kwargs)
        # One-ULP library variation propagated coherently through the analytic
        # jets; retained input states, gauge and source are not perturbed.
        direction = np.where(np.arange(len(U)) % 2, np.inf, -np.inf)
        r = np.where(U == 0, x, np.nextafter(row['r'], direction))
        L = np.log(x/r)+(x-r)/.8
        C = np.exp(L)
        p, R = C*(.8-x)/(4*x), (1-U/3.2)*C
        LU, Lx = -(1/r+1/.8)*p, 1/x+1/.8-(1/r+1/.8)*R
        row.update(r=r,L=L,C=C,r_Ubar=p,r_x=R,log_C_Ubar=LU,log_C_x=Lx,
                   r_Ubar_Ubar=LU*p,r_xx=Lx*R)
        return row
    monkeypatch.setattr(diagnostic, 'vacuum_slice', shifted)
    rebuilt = diagnose(data, budget)
    bounds = departure_history_roundoff_tolerances(rebuilt, data['source'], .2)
    for path, bound in bounds.items():
        a, b = original, rebuilt
        for key in path:
            a, b = a[key], b[key]
        assert np.all(abs(np.asarray(a)-b) <= bound), path
    assert max(np.max(abs(original[k]-rebuilt[k])) for k in ('E_Ubar','E_x')) > 2e-13
    assert max(np.max(b) for b in bounds.values()) < 1e-7
    assert {p[0] for p in bounds} == {'E_Ubar','E_x','summaries'}
    assert all(p[-1] != 'count' for p in bounds)
    budget.check()


def test_conditioned_comparison_rejects_physical_drift_and_unchanged_metadata():
    import copy
    from record_nsc_atlas_coupled_evolution import _equal
    saved = {'E_Ubar': [[0., 1.]], 'summaries': {'count': 2}, 'mass': .4}
    rebuilt = copy.deepcopy(saved)
    bounds = {('E_Ubar',): np.array([[1e-10, 2e-10]])}
    rebuilt['E_Ubar'][0][0] = 5e-11
    _equal(saved, rebuilt, 'control', roundoff=bounds)
    rebuilt['E_Ubar'][0][0] = 1e-10+1e-13
    with pytest.raises(ValueError, match='conditioned replay'):
        _equal(saved, rebuilt, 'control', roundoff=bounds)
    zero_bounds = {('E_Ubar',): np.zeros((1,2))}
    rebuilt['E_Ubar'][0][0] = 1e-15
    with pytest.raises(ValueError, match='conditioned replay'):
        _equal(saved, rebuilt, 'control', roundoff=zero_bounds)
    rebuilt['E_Ubar'][0][0] = 1e-6
    with pytest.raises(ValueError, match='conditioned replay'):
        _equal(saved, rebuilt, 'control', roundoff=bounds)
    for key, changed in [('summaries', {'count': 2.00000001}), ('mass', .40001)]:
        rebuilt = copy.deepcopy(saved)
        rebuilt[key] = changed
        with pytest.raises(ValueError, match='unconditioned'):
            _equal(saved, rebuilt, 'control', roundoff=bounds)
    for value in (np.nan, True, 0):
        rebuilt = copy.deepcopy(saved)
        rebuilt['E_Ubar'] = [[value, value]]
        with pytest.raises(ValueError, match='conditioned replay'):
            _equal(saved, rebuilt, 'control', roundoff=bounds)
    with pytest.raises(ValueError, match='inventory'):
        _equal(saved, saved, 'control', roundoff={('absent',): 1e-10})


def test_compact_source_roundoff_tracks_label_cancellation_not_signal_error():
    geometry = model.source_geometry()
    points = np.r_[geometry["x_corner"], np.linspace(*geometry["support_x"], 15)[1:-1], geometry["x_last"]]
    p, px, _ = model.source_jet(points, geometry)
    allowance_g, allowance_gx, allowance_gxx = source_jet_roundoff_tolerances(points, geometry, second=True)
    V = radius_product(points)/geometry["U0"]
    # Independent direct compact polynomial, with a one-ULP libm/null-label
    # perturbation. The analytic preparation itself has not been changed.
    shifted = np.nextafter(V, np.inf)
    z = 4*(shifted-(geometry["Va"]+geometry["Vb"])/2)/(geometry["Vb"]-geometry["Va"])
    zx = -4*points*np.exp(points/.8)/(.64*geometry["U0"]*(geometry["Vb"]-geometry["Va"]))
    active = abs(z) < 1
    g = np.where(active, (1-z*z)**4, 0.)
    gx = np.where(active, -8*z*(1-z*z)**3*zx, 0.)
    zxx = -4*np.exp(points/.8)*(1+points/.8)/(.64*geometry["U0"]*(geometry["Vb"]-geometry["Va"]))
    gxx = np.where(active, 8*(1-z*z)**2*(7*z*z-1)*zx*zx-8*z*(1-z*z)**3*zxx, 0.)
    assert np.all(abs(g-points*p) <= allowance_g+2e-13)
    assert np.all(abs(gx-(p+points*px)) <= allowance_gx+2e-13)
    pxx = model.source_jet(points, geometry)[2]
    assert np.all(abs(gxx-(2*px+points*pxx)) <= allowance_gxx+2e-13)
    assert np.all(allowance_g[[0, -1]] == 0) and np.all(allowance_gx[[0, -1]] == 0)
    assert 0 < np.max(allowance_gx) < 1e-8
    assert 0 < np.max(allowance_gxx) < 1e-6
    original = source_jet_roundoff_tolerances(points, geometry)
    assert len(original) == 2 and np.array_equal(original[1], allowance_gx)


@pytest.mark.parametrize("alpha", [0., 1e-4])
def test_snapshot_frame_residuals_keep_cancellation_scale(monkeypatch, alpha):
    import recursive_horizons.nsc_atlas_characteristic_trace as trace
    geometry = model.source_geometry()
    x = np.linspace(geometry["x_corner"], geometry["x_last"], 129)
    frame = object.__new__(trace.SourceSnapshot)
    frame._geometry, frame._alpha = geometry, alpha
    frame._bounds = dict(r_min=.2, v_max=2., F_min=.5)
    h, hx = 30*(x-.4)**3, 90*(x-.4)**2
    seed = frame._frame(x, h, hx, np.zeros_like(x))
    # An algebraic cancellation control, not a solution or fitted source.
    hxx = seed["source_ode_h_xx"]
    original = frame._frame(x, h, hx, hxx)
    assert np.all(original["source_h_xx_defect"] == 0)
    V = np.nextafter(radius_product(x)/geometry["U0"], np.inf)
    z = 4*(V-(geometry["Va"]+geometry["Vb"])/2)/(geometry["Vb"]-geometry["Va"])
    zx = -4*x*np.exp(x/.8)/(.64*geometry["U0"]*(geometry["Vb"]-geometry["Va"]))
    zxx = -4*np.exp(x/.8)*(1+x/.8)/(.64*geometry["U0"]*(geometry["Vb"]-geometry["Va"]))
    active = abs(z) < 1
    g = np.where(active, (1-z*z)**4, 0.)
    gx = np.where(active, -8*z*(1-z*z)**3*zx, 0.)
    gxx = np.where(active, 8*(1-z*z)**2*(7*z*z-1)*zx*zx-8*z*(1-z*z)**3*zxx, 0.)
    p = g/x; px = (gx-p)/x; pxx = (gxx-2*px)/x
    monkeypatch.setattr(trace, "source_jet", lambda *args: (p, px, pxx))
    changed = frame._frame(x, h, hx, hxx)
    # Also cover the separately reported smooth mixed-equation defects.
    for row in (original, changed):
        row.update(r_Ubar=np.full_like(x, .1), c_Ubar=x*x, log_C_Ubar=-x)
        rhs = trace._mixed(row, alpha)
        row["m"] = rhs["m"]
        for key, output in trace.TRANSVERSE.items():
            row[key+"_x"] = np.full_like(x, .25)
            row["mixed_"+key+"_defect"] = row[key+"_x"]-rhs[output]
    allowance = source_snapshot_roundoff_tolerances(original, geometry, alpha)
    for key, bound in allowance.items():
        assert np.all(abs(changed[key]-original[key]) <= bound+2e-13), key
    error = abs(changed["source_h_xx_defect"])
    assert np.max(error) > 2e-13  # A relative tolerance on the zero residual fails.
    assert np.max(allowance["source_h_xx_defect"]) < 1e-6
    print("snapshot conditioned replay", alpha, "max Q-defect change", float(np.max(error)),
          "max allowance", float(np.max(allowance["source_h_xx_defect"])))


def test_saved_hermite_coefficient_conditioning_covers_second_derivative():
    from scipy.interpolate import CubicHermiteSpline
    geometry = model.source_geometry()
    x = np.linspace(geometry["x_corner"], geometry["x_last"], 129)
    t = x-.4
    h, hx = 1+2*t+3*t*t+4*t*t*t, 2+6*t+12*t*t
    points = x[:-1]+np.diff(x)*.37
    original = CubicHermiteSpline(x, h, hx)
    # Alternating one-ULP coefficient perturbations avoid the exact cancellation
    # of a uniform vertical translation of every interpolation knot.
    toward = np.where(np.arange(len(x)) % 2, np.inf, -np.inf)
    changed = CubicHermiteSpline(x, np.nextafter(h, toward), np.nextafter(hx, toward))
    row = dict(x=points, h=original(points), h_x=original(points, 1), h_xx=original(points, 2),
               c=np.zeros_like(points), c_x=np.zeros_like(points), c_xx=np.zeros_like(points))
    allowance = source_snapshot_roundoff_tolerances(row, geometry, 0.,
        source_knots=dict(x=x, h=h, h_x=hx))
    for order, key in enumerate(("h", "h_x", "h_xx")):
        assert np.all(abs(changed(points, order)-row[key]) <= allowance[key]), key
    assert np.max(abs(changed(points, 2)-row["h_xx"])) > 1e-11
    assert np.max(allowance["h_xx"]) < 1e-6


def test_transverse_hermite_value_and_derivative_roundoff_is_propagated():
    import copy
    import recursive_horizons.nsc_atlas_characteristic_trace as trace
    geometry = model.source_geometry()
    alpha = 1e-4
    x = np.linspace(geometry["x_corner"], geometry["x_last"], 129)
    seed = object.__new__(trace.SourceSnapshot)
    seed._geometry, seed._alpha = geometry, alpha
    seed._bounds = dict(r_min=.2, v_max=2., F_min=.5)
    rows = seed._frame(x, 30*(x-.4)**3, 90*(x-.4)**2, 180*(x-.4))
    rU = (.8/x-1)/4
    rows.update(r_Ubar=rU, log_C_Ubar=-(1/.8+1/x)*rU, c_Ubar=-(x-.4)**2)
    original = trace.SourceSnapshot(geometry, alpha=alpha, x_nodes=x, source_rows=rows)
    changed = copy.deepcopy(original)
    for polynomial in changed._transverse.values():
        # Model one-ULP coefficient arithmetic, without changing any saved knot.
        polynomial.c[:] = np.nextafter(polynomial.c, np.inf)
    points = x[:-1]+np.diff(x)*.37
    base, replay = original.sample(points), changed.sample(points)
    limits = source_snapshot_roundoff_tolerances(base, geometry, alpha, source_knots=rows)
    changes = []
    for key in trace.TRANSVERSE:
        for field in (key, key+"_x", "mixed_"+key+"_defect"):
            error = abs(replay[field]-base[field])
            assert np.all(error <= limits[field]+2e-13), field
            changes.append(float(np.max(error)))
    assert max(changes) > 0
    assert np.all(limits["r_Ubar_x"] > 8*np.finfo(float).eps*abs(base["r_Ubar_x"]))


def _interval_diagnostic(q, knots):
    q, knots = np.asarray(q), np.asarray(knots)
    left = np.searchsorted(knots, q, side='right')-1
    right = np.searchsorted(knots, q, side='left')
    return dict(x=np.linspace(.4,.5,len(q)).tolist(),q=q.tolist(),
        left_knot_index=left.tolist(),right_knot_index=right.tolist(),
        strictly_interior_knots_per_saved_interval=np.maximum(0,right[1:]-left[:-1]-1).tolist(),
        interpretation='saved-interval knot brackets')


@pytest.mark.parametrize('toward', [-np.inf,np.inf])
def test_profile_bracket_replay_accepts_consistent_one_ulp_ties(toward):
    import copy
    from recursive_horizons.nsc_atlas_validation import profile_interval_replay_equal
    knots = [0.,.25,.5,.75,1.]
    queries = [0.,.25,.6,.75,1.]
    saved = [_interval_diagnostic(queries,knots)]
    shifted = np.nextafter(queries,toward).tolist()
    rebuilt = [_interval_diagnostic(shifted,knots)]
    before = copy.deepcopy((saved,rebuilt))
    profile_interval_replay_equal(saved,knots,rebuilt,knots)
    # A rebuilt knot can move instead of the query; counts must move with it.
    moved_knots = np.nextafter(knots,toward).tolist()
    profile_interval_replay_equal(saved,knots,[_interval_diagnostic(queries,moved_knots)],moved_knots)
    assert (saved,rebuilt) == before


@pytest.mark.parametrize('key', ['left_knot_index','right_knot_index',
                                'strictly_interior_knots_per_saved_interval'])
def test_profile_bracket_replay_rejects_unjustified_integer_changes(key):
    import copy
    from recursive_horizons.nsc_atlas_validation import profile_interval_replay_equal
    knots = [0.,.25,.5,.75,1.]
    saved = [_interval_diagnostic([.1,.3,.6,.9],knots)]
    changed = copy.deepcopy(saved)
    changed[0][key][1] += 1
    with pytest.raises(ValueError,match='integer is inconsistent'):
        profile_interval_replay_equal(saved,knots,changed,knots)
    changed[0][key][1] = float(saved[0][key][1])
    with pytest.raises(ValueError,match='integer is inconsistent'):
        profile_interval_replay_equal(saved,knots,changed,knots)


def test_profile_bracket_replay_rejects_coordinate_drift_and_multiple_knot_crossings():
    from recursive_horizons.nsc_atlas_validation import profile_interval_replay_equal
    knots = [0.,.25,.5,.75,1.]
    saved = [_interval_diagnostic([.1,.3,.6,.9],knots)]
    with pytest.raises(ValueError,match='query or metadata'):
        profile_interval_replay_equal(saved,knots,[_interval_diagnostic([.1,.30000001,.6,.9],knots)],knots)
    with pytest.raises(ValueError,match='knot coordinates'):
        profile_interval_replay_equal(saved,knots,saved,[0.,.25000001,.5,.75,1.])
    close_knots = [0.,.5,float(np.nextafter(.5,np.inf)),1.]
    a = [_interval_diagnostic([.1,np.nextafter(.5,-np.inf),.9],close_knots)]
    b = [_interval_diagnostic([.1,np.nextafter(close_knots[2],np.inf),.9],close_knots)]
    with pytest.raises(ValueError,match='not adjacent'):
        profile_interval_replay_equal(a,close_knots,b,close_knots)
    with pytest.raises(ValueError,match='more than one possible knot'):
        profile_interval_replay_equal(a,close_knots,[_interval_diagnostic([.1,.5,.9],close_knots)],close_knots)
    for invalid in ([.1,np.nan,.6,.9],[.1,.3,.2,.9]):
        with pytest.raises(ValueError,match='finite increasing'):
            profile_interval_replay_equal(saved,knots,[_interval_diagnostic(invalid,knots)],knots)
