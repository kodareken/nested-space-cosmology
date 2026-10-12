"""Small Fock/block identities and read-only prepared-source work forecasts."""
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.linalg import expm

from recursive_horizons import nsc_discovery_spectral_work as work
import derive_nsc_discovery_spectral_work as cli


def test_two_mode_Fock_mean_noise_and_Gaussian_reference_ratio():
    annihilator=np.array([[0.,1.],[0.,0.]],complex)
    parity=np.diag([1.,-1.])
    fields=[np.kron(annihilator,np.eye(2)),np.kron(parity,annihilator)]
    G=np.array([[.3,.7],[.7,-.3]],complex)
    operator=sum(fields[a].conj().T@fields[b]*G[a,b] for a in range(2) for b in range(2))
    angle=.23;U=expm(1j*angle*G)
    characteristics=[];noises=[];means=[]
    for lower,upper in ((0.,.6),(1.,0.),(1.,.6)):
        C=np.diag([lower,upper]);rho=np.diag([(1-lower)*(1-upper),(1-lower)*upper,lower*(1-upper),lower*upper])
        mean=np.trace(rho@operator).real
        variance=np.trace(rho@operator@operator).real-mean**2
        assert mean==pytest.approx(np.trace(C@G).real,abs=1e-14)
        assert variance==pytest.approx(np.trace(C@G@(np.eye(2)-C)@G).real,abs=1e-14)
        Fock_Z=np.trace(rho@expm(1j*angle*operator))
        determinant=np.linalg.det(np.eye(2)-C+C@U)
        assert Fock_Z==pytest.approx(determinant,abs=1e-14)
        characteristics.append(Fock_Z);noises.append(variance);means.append(mean)
    assert means[2]-means[1]==pytest.approx(means[0],abs=1e-14)
    assert abs(noises[2]-noises[1]-noises[0])>.1
    assert abs(characteristics[2]/characteristics[1]-characteristics[0])>.001


def test_finite_pulse_endpoints_area_and_stable_transform():
    settings=work.pulse_parameters(1.4)
    assert work.pulse(0.,settings)==(0.,0.)
    assert work.pulse(settings['T'],settings)==(0.,0.)
    assert quad(lambda t:work.pulse(t,settings)[0],0.,settings['T'],epsabs=1e-14)[0]==pytest.approx(0.,abs=1e-14)
    for gap in (0.,settings['omega'],1.4):
        direct=quad(lambda t:work.pulse(t,settings)[0]*np.cos(gap*t),0.,settings['T'],epsabs=1e-14)[0]+1j*quad(lambda t:work.pulse(t,settings)[0]*np.sin(gap*t),0.,settings['T'],epsabs=1e-14)[0]
        assert work.pulse_transform(gap,settings)==pytest.approx(direct,abs=1e-14)


def test_small_block_exact_work_and_passive_control_without_energy_subtraction():
    energies=np.array([[-.7,.7]])
    G=np.array([[[.3,np.sqrt(.91)],[np.sqrt(.91),-.3]]],complex)
    settings=work.pulse_parameters(1.4,amplitude=.002)
    endpoint=work.evolve_blocks(energies,G,settings,cpu_limit=5.)
    data={'energies':energies,'Gblocks':G,'positive_occupations':np.array([.6]),'degeneracy':np.ones(1),'M':4.,'Q0':1.,'N0':2.,'momenta':np.array([1.])}
    measured=work.measurement(data,endpoint);predicted=work.forecast(data,settings)
    assert endpoint['unitarity_max']<1e-10
    assert measured['work_full_C']<0.
    assert measured['work_passive_reference']>0.
    assert measured['work_vacuum_reference']>0.
    assert abs(measured['energy_work_closure_full'])<1e-12
    assert abs(measured['energy_work_closure_reference'])<1e-12
    smaller=work.pulse_parameters(1.4,amplitude=.001)
    small_endpoint=work.evolve_blocks(energies,G,smaller,cpu_limit=5.)
    small_measured=work.measurement(data,small_endpoint);small_predicted=work.forecast(data,smaller)
    large_error=abs(measured['work_full_C']-predicted['second_order_work_full_C'])/settings['amplitude']**2
    small_error=abs(small_measured['work_full_C']-small_predicted['second_order_work_full_C'])/smaller['amplitude']**2
    assert small_error<.6*large_error
    assert small_predicted['second_order_work_full_C']==pytest.approx(predicted['second_order_work_full_C']/4.,abs=1e-18)
    assert measured['CAR_min']>=-1e-12 and measured['CAR_max']<=1.+1e-10
    assert measured['mean_reference_subtraction_gap']==pytest.approx(0.,abs=1e-18)
    assert measured['energy_from_transition_probabilities_not_large_subtraction'] is True


@pytest.fixture(scope='module')
def preview():
    return work.preview()


def test_actual_prepared_source_poles_residues_and_work_preview(preview):
    assert preview['evolved'] is False
    assert preview['Q0']==pytest.approx(.4639447144311824,abs=1e-13)
    assert preview['r0']==pytest.approx(2.746344072669583,abs=1e-12)
    assert preview['M']==4. and preview['kappa']==1.
    assert preview['source_checks']['eigen_equation_max']<1e-10
    channels=preview['prediction']['channels']
    assert [row['gap'] for row in channels[:3]]==pytest.approx([1.215659930762963,2.532317371066923,4.035125235837996])
    assert [row['positive_frequency_residue_per_channel'] for row in channels[:3]]==pytest.approx([-.3130522443695124,-.43286852492854094,-.23678039659751357])
    assert all(row['degeneracy']==2 for row in channels)
    assert preview['prediction']['second_order_work_full_C']==pytest.approx(-2.269721447866569e-5)
    assert preview['prediction']['second_order_work_passive_reference']>0.
    assert abs(preview['prediction']['mean_reference_subtraction_gap'])<1e-18
    assert preview['Gaussian_reference_gap']['characteristic_ratio_minus_full_abs']>.01
    assert preview['prediction']['no_pole_evaluated'] is True
    assert preview['prediction']['physical_filled_sea_Gamma_matching']=='unresolved'
    assert preview['prediction']['autonomous_regeneration_claimed'] is False


def test_prepare_predict_are_exclusive_and_never_launch_actual_full_pulse(tmp_path,monkeypatch):
    before={str(work.INPUT.with_suffix(ext)):work.sha256(work.INPUT.with_suffix(ext)) for ext in ('.json','.npz')}
    monkeypatch.setattr(work,'evolve_blocks',lambda *args,**kwargs:pytest.fail('full nonlinear pulse ran during preparation/prediction/check'))
    prepared=work.prepare(tmp_path)
    assert prepared['source_operator_frozen_before_prediction'] is True
    assert not (tmp_path/'prediction.json').exists()
    assert work.check(tmp_path)['ok'] is True
    locked=work.predict(tmp_path)
    assert locked['locked_before_nonlinear_pulse'] is True
    assert locked['nonlinear_measurement_performed'] is False
    assert work.check(tmp_path)['ok'] is True
    assert not (tmp_path/'measurement.json').exists()
    with pytest.raises(FileExistsError):work.prepare(tmp_path)
    with pytest.raises(FileExistsError):work.predict(tmp_path)
    assert {str(work.INPUT.with_suffix(ext)):work.sha256(work.INPUT.with_suffix(ext)) for ext in ('.json','.npz')}==before
    payload=tmp_path/'prediction.npz';payload.chmod(0o644);payload.write_bytes(payload.read_bytes()+b'x')
    with pytest.raises(ValueError,match='payload binding'):work.check(tmp_path)


def test_CLI_default_is_preview_and_creation_path_is_scoped(tmp_path,monkeypatch,capsys):
    monkeypatch.setattr(work,'preview',lambda:{'mode':'preview','evolved':False})
    assert cli.main([])==0
    assert json.loads(capsys.readouterr().out)['evolved'] is False
    monkeypatch.setattr(work,'OUTPUT',tmp_path/'pulse')
    assert cli.output_path(tmp_path/'pulse'/'new')==tmp_path/'pulse'/'new'
    with pytest.raises(PermissionError):cli.output_path(tmp_path/'old')
