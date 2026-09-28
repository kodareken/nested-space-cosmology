"""Explicit order, inherited proof owners and recurrence-free replay."""
import importlib.util
import json
from pathlib import Path

import mpmath as mp
import pytest

from recursive_horizons.nsc_incoming_retained_order24_bound import retained_order24_bound
from recursive_horizons.nsc_incoming_centered_order24 import CenteredOrder24Geometry, trimmed_defect_jets
from recursive_horizons.nsc_incoming_defect_taylor_bound import centered_coefficient_enclosure, intersect_complex
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision, _lo, _hi, _range

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT/'results/development/nsc-incoming-retained-order24-bound.json'


def test_invalid_domain_rejected_before_any_coefficients():
    for channel, kwargs in (({'index': 0}, {}), ({'index': 12}, {'intervals': 64}),
                            ({'index': 12}, {'depth': 8})):
        with pytest.raises(ValueError, match='fixed128'):
            retained_order24_bound({}, channel, **kwargs)


def test_frozen_geometry_and_proof_owners_are_reused():
    from recursive_horizons import nsc_incoming_retained_order24_bound as owner
    assert owner.trimmed_defect_jets is trimmed_defect_jets
    assert owner.centered_coefficient_enclosure is centered_coefficient_enclosure
    config = json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())['scattering_provenance']['config']
    geometry = CenteredOrder24Geometry(config)
    with _precision(40):
        assert len(geometry.cells) == 128
        assert _hi(geometry.cells[0][0]) >= _hi(geometry.horizon)
        assert _lo(geometry.cells[-1][0]) <= 1
        for i, (rho, center) in enumerate(geometry.cells):
            assert _lo(rho) <= _lo(center) <= _hi(center) <= _hi(rho)
            if i: assert _hi(rho) >= _lo(geometry.cells[i-1][0])


def test_normalized_complex_remainder_and_intersection():
    with _precision(40):
        h = _range(mp.mpf('-.25'), mp.mpf('.25')); c = mp.iv.mpc(1, -2)
        whole = [c*h**4, 4*c*h**3, 6*c*h*h, 4*c*h, c]
        value = centered_coefficient_enclosure(whole, [mp.iv.mpc(0)]*4+[c], h, 4)
        exact = c*mp.iv.mpf('.25')**4
        assert _lo(value.real) <= _lo(exact.real) <= _hi(exact.real) <= _hi(value.real)
        assert _lo(value.imag) <= _lo(exact.imag) <= _hi(exact.imag) <= _hi(value.imag)
        with pytest.raises(ArithmeticError, match='empty'):
            intersect_complex(mp.iv.mpc(1), mp.iv.mpc(2))


def test_record_replay_keeps_order24_conditional_and_cannot_prepare(monkeypatch):
    if not RESULT.exists(): pytest.skip('record is prepared after the pre-compute owner gate')
    spec = importlib.util.spec_from_file_location('retained24', ROOT/'scripts/derive_nsc_incoming_retained_order24_bound.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    def forbidden(*args, **kwargs): raise AssertionError('replay must not regenerate coefficients')
    monkeypatch.setattr(module, 'retained_order24_bound', forbidden)
    monkeypatch.setattr(module, 'local_cross_product_coefficients', forbidden)
    saved = json.loads(RESULT.read_text()); actual = module.replay(ROOT/saved['payload']['path'])
    assert actual == saved
    assert actual['group'] == 12 and actual['energy_interval'] == [40., 320.]
    assert actual['scope']['explicit_order24_source_correction_required']
    assert not actual['scope']['archived_order16_accuracy_certified']
    assert not actual['scope']['constraints_solved']
    assert actual['radial']['physical_Riccati_order'] == 24
