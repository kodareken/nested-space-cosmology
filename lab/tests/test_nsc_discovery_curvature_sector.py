"""Constrained action identities, physical pole signs and immutable replay."""
from pathlib import Path
import tempfile

import pytest
import sympy as sp

from recursive_horizons import nsc_discovery_curvature_sector as sector


@pytest.fixture(scope="module")
def assessed():
    return sector.assess()


def test_constrained_flat_variation_and_massive_branch(assessed):
    assert assessed["status"]=="CONSTRAINED_FLAT_EXTRA_POLE_CONFIRMED"
    assert all(assessed["checks"].values())
    assert assessed["checks"]["auxiliary_equation_constraint_compatibility"]
    assert assessed["checks"]["Einstein_Coulomb_separation"]
    scales=assessed["proper_scales"]
    assert scales["physical_spin2_pole_mass_square"]==pytest.approx(-26.6666666666667)
    assert scales["long_wavelength_growth_time"]==pytest.approx(.1936491673103707)
    assert scales["actual_Lorentz_Weyl_coefficient"]==-scales["C_W"]>0
    assert scales["tachyonic_extra_pole"] and scales["opposite_Einstein_pole_residue"]
    assert not assessed["evidence_written"] and assessed["producing_commit"] is None
    assert not assessed["scope"]["periodic_whole_domain_initial_data"]
    assert not assessed["scope"]["all_source_instability_proved"]
    assert not assessed["scope"]["EFT_pole_validity_established"]


def test_opposite_weyl_sign_changes_mass_but_not_extra_pole_residue():
    current=sector.proper_pole_scales(.045,-.00084)
    opposite=sector.proper_pole_scales(.045,.00084)
    assert current["physical_spin2_pole_mass_square"]<0<opposite["physical_spin2_pole_mass_square"]
    assert not opposite["tachyonic_extra_pole"]
    assert opposite["long_wavelength_growth_rate"] is None
    assert opposite["opposite_Einstein_pole_residue"]
    with pytest.raises(ValueError,match="nonzero"):
        sector.proper_pole_scales(.045,0.)


def test_variational_checks_detect_an_inconsistent_owned_radial_coupling(monkeypatch):
    original=sector.action.partial_F
    def wrong_radial(r,A,CW):
        Fr,f=original(r,A,CW)
        return sp.Rational(101,100)*Fr,f
    monkeypatch.setattr(sector.action,"partial_F",wrong_radial)
    checks=sector.symbolic_checks()
    assert not checks["flat_Euler_background"]
    assert not checks["radial_Euler_linearization"]
    assert not checks["canonical_p_Q"]


def test_creation_only_and_full_symbolic_replay(assessed,monkeypatch):
    monkeypatch.setattr(sector,"verify_producer",lambda commit,hashes:commit)
    bound=dict(assessed,producing_commit="test-frozen-source-revision")
    with tempfile.TemporaryDirectory(dir="/tmp",prefix="nsc-curvature-sector-") as directory:
        target=Path(directory)/"assessment.json"
        with pytest.raises(ValueError,match="frozen"):
            sector.write_record(assessed,target)
        sector.write_record(bound,target)
        before=target.read_bytes()
        checked=sector.check_record(target)
        assert checked["ok"] and checked["bytes_written"]==0 and checked["symbolic_recomputed"]
        assert target.read_bytes()==before
        with pytest.raises(FileExistsError):
            sector.write_record(bound,target)
        import json
        altered=json.loads(before)
        altered["proper_scales"]["physical_spin2_pole_mass_square"]*=-1
        target.write_text(json.dumps(altered))
        with pytest.raises(ValueError,match="proper_scales"):
            sector.check_record(target)
