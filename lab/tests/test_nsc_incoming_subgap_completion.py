"""Mixed-channel source normalization; completed-pilot replay is added after preparation."""
import numpy as np
import pytest

from recursive_horizons.nsc_incoming_subgap_completion import recover_unsewn, mixed_vertices, MAGNITUDE
from recursive_horizons.nsc_transmitting_dirac_domain import S2
from recursive_horizons.nsc_incoming_state_moments import RADIUS


def test_recover_actual_subgap_columns_and_signed_sphere_vertex():
    v = np.array([.6, .8j]); w = np.array([.8j, .6]); R = np.exp(.7j)
    physical = np.column_stack((R*v, w, np.zeros(2)))
    assert np.max(abs(recover_unsewn(physical, np.array(R))-np.column_stack((v, w)))) < 3e-15
    assert np.max(abs(mixed_vertices(1.2, MAGNITUDE)[3]-MAGNITUDE/(2*RADIUS)*S2)) == 0
    assert np.max(abs(mixed_vertices(1.2, -MAGNITUDE)[3]+MAGNITUDE/(2*RADIUS)*S2)) == 0
    with pytest.raises(ValueError, match='unit reflection'):
        recover_unsewn(physical, np.array(.5*R))
