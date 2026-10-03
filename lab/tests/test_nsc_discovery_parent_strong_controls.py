"""Refusal and linearity checks for the read-only matched control consumer."""
import json
import numpy as np
import pytest
import assess_nsc_discovery_parent_strong_controls as consumer


def initial():
    return {**{n: np.arange(4., dtype=float) for n in consumer.EQUAL},
            'source_weights': np.array([.3, .15])}


@pytest.mark.parametrize('name', consumer.EQUAL)
def test_reuse_rejects_any_initial_mismatch(name):
    strong = initial(); frozen = initial(); frozen[name] = frozen[name].copy()
    frozen[name][1] += 1e-15
    with pytest.raises(ValueError, match='initial equality'):
        consumer.exact_initial(strong, frozen)


def test_linear_reweight_equals_direct_column_probability():
    rng = np.random.default_rng(9)
    phi = rng.normal(size=(16, 2))+1j*rng.normal(size=(16, 2))
    weights, strong_weights = np.array([.001, .002]), np.array([.3, .15])
    child = (abs(phi[:5])**2*weights).sum(axis=0)
    total = (abs(phi)**2*weights).sum(axis=0)
    result = consumer.reweight({'child': child, 'ambient_total': total}, strong_weights/weights)
    direct = float((abs(phi[:5])**2*strong_weights).sum())
    assert result['child_probability'] == pytest.approx(direct, abs=1e-14)
    assert result['child_fraction'] == pytest.approx(direct/(abs(phi)**2*strong_weights).sum())
    assert consumer.reweight({'child': [0,0], 'ambient_total': [0,0]}, [1,1])['child_fraction'] is None


def test_proper_clock_interpolation_refuses_extrapolation_or_bad_clock():
    assert consumer.interpolate([0., 1., 3.], [0., 2., 4.], 2.) == 3.
    assert consumer.interpolate([0., 1., 3.], [0., 2., 4.], 3.) == 4.
    assert consumer.interpolate([0., 1., 3.], [0., 2., 4.], 3.01) is None
    assert consumer.interpolate([0., 1., 3.], [0., 2., 4.], -.01) is None
    for clock in ([0.,0.,1.], [0.,2.,1.], [0.,float('nan'),1.]):
        with pytest.raises(ValueError):
            consumer.interpolate(clock, [0,1,2], .5)


def test_clock_anchor_uses_exact_stored_tau_and_positive_rates():
    rows = [dict(time=t, centre_rate=1.) for t in (0., .5, 1.)]
    stations = [dict(time=0., centre_tau=0.), dict(time=1., centre_tau=2.)]
    tau, gap = consumer.anchored_clocks(rows, stations)
    assert tau.tolist() == [0.,1.,2.]
    assert gap == 1.
    rows[1]['centre_rate'] = 0.
    with pytest.raises(ValueError, match='positive clock'):
        consumer.anchored_clocks(rows, stations)


def test_creation_only_readonly_writer_and_input_directory_guard(tmp_path, monkeypatch):
    path = tmp_path/'new.json'; consumer.write_record(path, {'value': 3}, [])
    assert json.loads(path.read_text()) == {'value': 3}
    assert path.stat().st_mode & 0o222 == 0
    with pytest.raises(FileExistsError):
        consumer.write_record(path, {'value': 4}, [])
    with pytest.raises(ValueError, match='sealed input'):
        consumer.write_record(tmp_path/'sealed.json', {}, [tmp_path])
    monkeypatch.setattr(consumer, 'LIMIT', 2)
    with pytest.raises(ValueError, match='64 MiB'):
        consumer.write_record(tmp_path/'large.json', {'value': 4}, [])
    assert not (tmp_path/'large.json').exists()
