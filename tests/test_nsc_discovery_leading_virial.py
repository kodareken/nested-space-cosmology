"""Static topology blocker, exact adjoints and boundary-sign assessment."""
from pathlib import Path
import tempfile

import pytest

import assess_nsc_discovery_leading_virial as virial


@pytest.fixture(scope="module")
def assessed():
    return virial.assess()


def test_static_boundary_sign_and_BR_electrovac():
    assert all(virial.symbolic_checks().values())


def test_actual_positive_spectral_state_SBP_and_adjoint(assessed):
    control=assessed["finite_control"]
    assert control["nf"]==32 and control["nq"]==128
    assert not control["initial_constraints_solved"] and not control["trajectory_performed"]
    assert min(control["eigenvalues"])>0
    assert control["source_trace"]==3 and control["multiplicity"]==4
    assert control["source_energy"]==pytest.approx(control["spectral_energy"],abs=1e-10)
    assert control["source_energy"]>0 and control["magnetic_energy"]>0
    assert control["lhs"]<0
    assert abs(control["SBP_identity_indicator"])<1e-10
    assert abs(control["adjoint_pairing_indicator"])<1e-10
    assert control["rhs"]==pytest.approx(control["lhs"],abs=1e-10)
    assert not control["continuum_or_roundoff_error_certificate"]


def test_actual_e6_binding_and_narrow_claim(assessed):
    assert assessed["leading_binding"]["producing_commit"]==virial.LEADING_PRODUCER
    assert not assessed["leading_binding"]["trajectory_replayed"]
    assert assessed["assumptions"]["angular_kappa_over_r_is_not_fundamental_4D_mass"]
    assert assessed["scope"]["STATIC_class_only"]
    assert not assessed["scope"]["dynamic_or_oscillatory_regimes_excluded"]
    assert not assessed["scope"]["NSC_verdict"]
    assert assessed["status"]=="STATIC_POSITIVE_ENERGY_PERIODIC_CLASS_EXCLUDED"
    assert not assessed["evidence_written"]


def test_immutable_create_readonly_replay_and_changed_sign_rejection(assessed,monkeypatch):
    monkeypatch.setattr(virial,"verify_producer",lambda commit,hashes:commit)
    with tempfile.TemporaryDirectory(dir="/tmp",prefix="leading-virial-") as directory:
        path=Path(directory)/"record.json"
        with pytest.raises(ValueError,match="frozen"):
            virial.write_record(assessed,path)
        virial.write_record(dict(assessed,producing_commit="test-source-freeze"),path)
        before=path.read_bytes()
        assert virial.check_record(path)["bytes_written"]==0
        assert path.read_bytes()==before
        with pytest.raises(FileExistsError):
            virial.write_record(assessed,path)
        import json
        changed=json.loads(before)
        changed["equations"]["boundary"]="opposite boundary sign"
        path.write_text(json.dumps(changed))
        with pytest.raises(ValueError,match="equations"):
            virial.check_record(path)
