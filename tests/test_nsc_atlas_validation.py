"""Cross-platform replay controls without retained field evolution."""
import numpy as np
import pytest
from recursive_horizons import nsc_atlas_characteristic_initialization as model
from recursive_horizons.nsc_atlas_field_preparation import radius_product
from recursive_horizons.nsc_atlas_validation import source_jet_roundoff_tolerances
from recursive_horizons.nsc_atlas_validation import source_snapshot_roundoff_tolerances


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
