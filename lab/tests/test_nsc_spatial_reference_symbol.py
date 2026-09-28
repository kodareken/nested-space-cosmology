"""Integration checks for the NSC spatial reference-symbol extension."""
import numpy as np
import pytest

from recursive_horizons.nsc_spatial_reference_symbol import (
    SymbolJet, star_order, supplied_KS_metric_jets, reference_projector,
)


def test_weyl_product_sign_and_noncommuting_coordinate():
    z, k = SymbolJet.variable(.2, 1), SymbolJet.variable(.7, 2)
    assert np.max(abs(star_order(z, k, 1).value-.5j*np.eye(2))) < 1e-14
    assert np.max(abs(star_order(k, z, 1).value+.5j*np.eye(2))) < 1e-14
    with pytest.raises(ValueError, match='order'):
        z.trim(0, 0).derivative(z=1)
    sigma1 = np.array([[0., 1.], [1., 0.]])
    sigma3 = np.diag([1., -1.])
    matrix = z.right(sigma1)
    assert np.max(abs((sigma3*matrix).value-.2*sigma3@sigma1)) < 1e-14
    assert np.max(abs((matrix*sigma3).value-.2*sigma1@sigma3)) < 1e-14


def test_spatial_projector_retains_scalar_correction():
    fields = supplied_KS_metric_jets(.019, .23)
    mass, angular, momentum = np.pi/2, np.sqrt(5), .2
    row = reference_projector(*fields, [momentum], [mass], [angular])
    N, beta, a, r = fields
    av, rv = a.value[0, 0, 0].real, r.value[0, 0, 0].real
    rz = r.derivative(z=1).value[0, 0, 0].real
    normal_energy = np.sqrt(mass*mass+(angular/rv)**2+(momentum/av)**2)
    expected = -.5*mass*angular*rz/(av*rv**2*normal_energy**3)
    observed = np.trace(row['orders'][1, 0])
    assert abs(observed-expected) < 3e-14
    assert abs(observed) > 1e-5
    assert max(float(v['scaled'].max()) for order in row['residuals'] for v in order.values()) < 3e-11


def test_zero_gap_is_rejected_not_regularized():
    c = SymbolJet.constant
    with pytest.raises(ValueError, match='gap|nonzero'):
        reference_projector(c(1), c(0), c(1), c(1), [0.], [0.], [0.])
