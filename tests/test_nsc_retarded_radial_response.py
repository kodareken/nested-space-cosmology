"""Frame/Fourier controls and saved radial replay; no repeated solves."""
import importlib.util
from pathlib import Path
import sys
import numpy as np
import json
import pytest
import pytest

from recursive_horizons.nsc_frozen_sources import frozen_source_context
# Original producer bytes are explicit; numerical fixtures remain current local inputs.
FROZEN_PRODUCER_CONTEXT = "329ef809e466a86ed4a290ee32c19485f85add9d"


@frozen_source_context(FROZEN_PRODUCER_CONTEXT)
def test_radial_replay_rejects_errors_beyond_operand_envelopes():
    from copy import deepcopy
    record = R.check()
    _parent, _fourier, data = R.inputs()
    with np.load(ROOT / R.ARTIFACT, allow_pickle=False) as archive:
        arrays = {k: archive[k] for k in archive.files}
    meta = json.loads(arrays['metadata_json'].tobytes())
    index = meta['upstream_index']
    expected = R.archived_amplitudes(data['reference/x'], data['reference/phi'], data['source/energies'], index)
    map_bound = R.archived_map_replay_envelope(data['reference/x'], data['reference/phi'], data['source/energies'], index)
    changed = arrays['reference/upstream'].copy()
    changed[0, 0, 0] += 4 * map_bound
    with pytest.raises(ValueError, match='arithmetic envelope'):
        R.validate_archived_map(expected, changed, data['reference/x'], data['reference/phi'], data['source/energies'], index)
    label = next(iter(record['per_solve']))
    field = 'diagonal_K_inverse_antiHermitian'
    bound = R.inverse_diagnostic_replay_envelopes(arrays, label)[field]
    changed_record = deepcopy(record)
    changed_record['per_solve'][label][field] += 4 * bound
    with pytest.raises(ValueError, match='arithmetic envelope'):
        R.validate_replay_record(record, changed_record, arrays)
    changed_record = deepcopy(record)
    changed_record['status'] = 'PASS: complete physical theory'
    with pytest.raises(ValueError, match='receipt differs'):
        R.validate_replay_record(record, changed_record, arrays)

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'src'))
SPEC=importlib.util.spec_from_file_location('radial_response',ROOT/'scripts/derive_nsc_retarded_radial_response.py')
R=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(R)
from recursive_horizons.nsc_retarded_radial_response import pulse_fourier,ks_generator,local_frame_checks


def test_same_pulse_transform_has_reality_and_translation_laws():
    omega=np.array([-.7,0.,.7]);a=pulse_fourier(omega,1.2,.06,192)
    assert abs(a[0]-a[2].conjugate())<1e-15
    assert a[1].real>0 and abs(a[1].imag)==0
    shifted=pulse_fourier(omega,1.5,.06,192)
    assert np.max(abs(shifted-np.exp(1j*omega*.3)*a))<1e-15
    assert np.max(abs(a-pulse_fourier(omega,1.2,.06,96)))<3e-11
    with pytest.raises(ValueError,match='fixed96/192'):pulse_fourier(omega,1.2,.06,48)


def test_current_and_unitary_KS_generators_and_forcing_agree():
    E=np.array([-.55,-.19,.19,.55]);result=local_frame_checks((1.,1.015,1.031),E,np.pi/2,np.sqrt(5))
    assert max(result.values())<3e-11
    G=ks_generator(1.01,E,np.pi/2,np.sqrt(5))
    assert np.max(abs(G+G.swapaxes(-1,-2).conj()))==0


def test_pair_comparison_uses_every_selected_pair():
    a=np.ones((4,4,2,3),complex);b=a.copy()
    assert R.pair_comparison(a,b)['every_pair_passes']
    a[3,2,0,1]+=0.02
    result=R.pair_comparison(a,b)
    assert not result['every_pair_passes']
    assert sum(not row['passes'] for row in result['pairs'])==1


@frozen_source_context(FROZEN_PRODUCER_CONTEXT)
def test_saved_replay_has_no_propagation_and_preserves_full_response(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('radial/PDE solve called during replay')
    monkeypatch.setattr(R,'short_response',forbidden)
    monkeypatch.setattr(R.R.P,'evolve_prepared_field_jets',forbidden)
    record=R.check()
    assert len(record['control']['completed'])<=2
    assert record['scope']['physical_low_energy_preparation']=='OPEN'
    assert record['scope']['rigorous_response_error_bound'] is None
    with np.load(ROOT/record['payload']['path'],allow_pickle=False) as archive:a={k:archive[k] for k in archive.files}
    if 'fine' not in record['control']['completed']:return
    B=a['fine/response'];f=a['reference/incoming'];C=a['source/source'];P=a['source/projector']
    direct=R.F.covariance_pairs(B,f,C)
    for i in (2,3):
        expected=B[i,i]@C[i]@f[i].conj().T+f[i]@C[i]@B[i,i].conj().T
        assert np.max(abs(expected-direct[i,i]))<1e-18
        K=np.linalg.solve(f[i,:,:2].T,B[i,i,:,:2].T).T
        CS=f[i]@C[i]@f[i].conj().T
        # A conditional numerical diagnostic; the matrix solve does not
        # normalize the physical source or its archived columns.
        assert np.max(abs(expected-(K@CS+CS@K.conj().T)))<3e-11
    assert record['per_solve']['fine']['CAR_pair_residual']==R.F.maxabs(R.F.covariance_pairs(B,f,P))
    assert not record['scope']['CAR_residual_subtracted']
