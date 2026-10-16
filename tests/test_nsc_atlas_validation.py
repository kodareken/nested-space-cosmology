"""Cross-platform replay controls without retained field evolution."""
import numpy as np
from recursive_horizons import nsc_atlas_characteristic_initialization as model
from recursive_horizons.nsc_atlas_field_preparation import radius_product
from recursive_horizons.nsc_atlas_validation import source_jet_roundoff_tolerances


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


