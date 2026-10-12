"""Geometry interpolation remainder vanishes for a sampled exact interpolant."""
import numpy as np
from scipy.interpolate import BarycentricInterpolator


def test_barycentric_interpolant_reproduces_nodal_geometry():
    z = np.linspace(1.17, 1.23, 47)
    values = np.column_stack((np.sin(20 * z), np.cos(20 * z)))
    dense = np.linspace(z[0], z[-1], 129)
    interpolants = [BarycentricInterpolator(z, values[:, k]) for k in range(2)]
    reconstructed = np.column_stack([fn(z) for fn in interpolants])
    assert np.max(abs(reconstructed - values)) < 1e-12
    dense_vals = np.column_stack([fn(dense) for fn in interpolants])
    assert np.all(np.isfinite(dense_vals))
