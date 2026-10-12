"""Independent source-measure and raw-vertex binding for the joint assembly."""
import json
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons.nsc_incoming_cauchy_jets import incoming_cauchy_jets, _raw_vertex_jets
from recursive_horizons.nsc_incoming_joint_constraints import source_action_gradient
from recursive_horizons.nsc_incoming_state_moments import incoming_group_factor, AXIAL, RADIUS
from recursive_horizons.nsc_spatial_reference_symbol import SIGMA


@pytest.mark.parametrize('group', [0, 1, 13, 14])
def test_folded_stress_measure_matches_two_sided_canonical_vertex(group):
    root = Path(__file__).resolve().parents[1]
    channels = json.loads((root/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels']
    channel = channels[group]
    domain = incoming_cauchy_jets()
    mass, ell = channel['compact_mass'], channel['angular_eigenvalue']
    signs = (1,) if ell == 0 else (-1, 1)
    # An algebraic insertion with nonzero trace tests the independent shift
    # vertex too. This is not a supplied physical stress or a state update.
    insertion = .1*SIGMA[0]-.07*SIGMA[1]+.03*SIGMA[2]+.02*np.eye(2)
    for k in (-2.3, 2.3):
        for sign in signs:
            angular = sign*ell
            raw = _raw_vertex_jets(domain.fields, np.array([k]), np.array([mass]), np.array([angular]))
            stress_vertices = [-mass*SIGMA[0]+angular/RADIUS*SIGMA[1]+k/AXIAL*SIGMA[2],
                               k/AXIAL*SIGMA[2], -k/AXIAL*np.eye(2), angular/(2*RADIUS)*SIGMA[1]]
            stress = np.array([np.trace(insertion@v).real for v in stress_vertices])*incoming_group_factor(channel, signs)
            projected = source_action_gradient(domain, stress)
            multiplicity = channel['copy_count']*channel['degeneracy']/len(signs)
            # Two folded signs times the canonical dk/(2*pi) measure.
            canonical = -2*multiplicity/(2*np.pi)*np.array([np.trace(insertion@v.value[0]).real for v in raw])
            np.testing.assert_allclose(projected, canonical, rtol=0, atol=3e-13)
