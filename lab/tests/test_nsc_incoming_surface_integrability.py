"""Nonzero first projector with zero constraint contractions and owned tails."""
import json
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons.nsc_incoming_surface_integrability import first_order_certificate,certificate,S1,S2,S3
from recursive_horizons.nsc_spatial_reference_symbol import SIGMA

ROOT=Path(__file__).resolve().parents[1]


def test_exact_first_order_contractions_keep_spatial_Weyl_terms():
    np.testing.assert_array_equal(np.array([S1.tolist(),S2.tolist(),S3.tolist()],complex),SIGMA)
    result=first_order_certificate()
    assert result['first_projector_generically_nonzero']
    assert all(set(values)=={'0'} for values in result['residuals'].values())
    for key in ('first_Weyl_idempotence','first_Weyl_transport_commutator',
                'N_symmetric_first_Moyal_trace','beta_symmetric_first_Moyal_trace'):
        assert key in result['residuals']


def test_imported_higher_powers_close_only_two_constraint_directions():
    uv=json.loads((ROOT/'results/development/nsc-reference-band-bulk.json').read_text())
    result=certificate(uv)
    assert [r['constraint_insertion_power_upper'] for r in result['imported_tail_powers']]==[-2,-3,-4]
    assert result['constraints']==['N','beta'] and result['reference_change_absolutely_integrable']
    assert not result['scope']['pressure_direction_integrability_claimed']
    assert not result['scope']['Gamma_rest_assigned']
    assert not result['scope']['spatial_domain_or_boundary_selected']
    uv['result']['recursion_power_certificate']['rows'][1]['P_power']=-2
    with pytest.raises(ValueError,match='UV power'):certificate(uv)
