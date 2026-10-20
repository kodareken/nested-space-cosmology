import math
import mpmath as mp
import numpy as np
import pytest

from recursive_horizons.nsc_atlas_flux_transition import FluxLaw
from recursive_horizons.nsc_atlas_flux_conversion import (
    ScalarFlux,sector_minimum,isolated_bounce,evolve_prepared,
)


def test_reachable_sector_depends_on_offset_not_transition_count():
    assert sector_minimum("1.5",".5")["minimum_vacuum"]=="0"
    assert sector_minimum("1.55",".5")["minimum_vacuum"]=="1/800"
    assert sector_minimum("2.05",".5")["minimum_vacuum"]=="1/800"
    assert sector_minimum("1.25",".5")["nearest_generations"]==[2,3]
    assert sector_minimum("-.1",".5")["nearest_generations"]==[0]
    assert sector_minimum("-.4",".5",inward_only=False)["nearest_generations"]==[-1]
    assert sector_minimum("1.5",".5",bare="-3")["minimum_vacuum"]=="-3"


@pytest.mark.parametrize("flux",[1.5,1.,.5])
def test_isolated_action_against_independent_offshell_function(flux):
    law=FluxLaw()
    row=isolated_bounce(law,flux)
    with mp.workdps(50):
        r=mp.mpf(row["radius"])
        li,lo=mp.mpf(str(row["child_Lambda"])),mp.mpf(str(row["parent_Lambda"]))
        sigma=mp.mpf("0.2")
        def bulk(L,R):
            return -R*R/2 if L==0 else (mp.power(1-L*R*R/3,mp.mpf("1.5"))-1)/L
        def action(R):
            return 2*mp.pi**2*sigma*R**3+12*mp.pi**2*(bulk(li,R)-bulk(lo,R))
        assert float(action(r))==pytest.approx(row["action"],rel=2e-14)
        assert abs(float(mp.diff(action,r)))<1e-11
    for key in ("angular_residual","inside_constraint","outside_constraint"):
        assert abs(row[key])<1e-13


def test_weak_gravity_bounce_limit():
    law=FluxLaw(kappa=1e-8)
    result=isolated_bounce(law,1.5)
    density_drop=.625
    flat=27*math.pi**2*law.tension**4/(2*density_drop**3)
    assert result["action"]==pytest.approx(flat,rel=1e-7)


def test_preparation_energy_cannot_be_created_by_relaxation():
    model=ScalarFlux()
    p=model.prepare()
    assert model.potential(p["inherited_phi"],p["child_flux"])-p["child_vacuum"]==pytest.approx(p["oscillation_energy"])
    assert p["parent_vacuum"]==pytest.approx(p["child_vacuum"]+p["oscillation_energy"]+p["released_to_transition"])
    assert p["released_to_transition"]==.0625
    assert p["actual_wall_pressure"]==pytest.approx(p["released_to_transition"])
    for field in (-.1,0.,.2):
        assert model.wall_pressure(field,1.5,.5)==pytest.approx(model.potential(field,1.5)-model.potential(field,1.))
    costly=ScalarFlux(epsilon=.01).prepare()
    assert costly["required_external_work"]==pytest.approx(.1175)
    with pytest.raises(ValueError,match="energy injection"):
        evolve_prepared(model=ScalarFlux(epsilon=.01))
    assert ScalarFlux(epsilon=0).prepare()["required_external_work"]==.125


def test_zero_charge_and_zero_mixing_controls():
    p=ScalarFlux().prepare(charge=0)
    assert p["displacement"]==0 and p["oscillation_energy"]==0 and p["released_to_transition"]==0
    p=ScalarFlux(epsilon=1).prepare()
    assert p["oscillation_energy"]==0 and p["child_vacuum"]==.5
    s,m=evolve_prepared(model=ScalarFlux(epsilon=1),efolds=.2)
    assert s.y[3,-1]==pytest.approx(m["initial_H"],abs=1e-14)


def test_proper_clock_conservation_and_independent_efold_clock():
    model=ScalarFlux()
    s,meta=evolve_prepared(rtol=2e-12,atol=2e-14)
    alt,_=evolve_prepared(clock="efolds",rtol=2e-12,atol=2e-14)
    x,v,N,H,W=s.sol(np.linspace(0,s.t[-1],257))
    energy=(v*v+model.omega**2*x*x)/2
    assert np.max(np.abs(3*H*H-model.kappa*(.125+energy)))<1e-10
    assert np.max(np.abs(energy+W-.09375))<1e-10
    assert s.t[-1]==pytest.approx(alt.y[2,-1],abs=1e-9)
    np.testing.assert_allclose(s.y[:2,-1],alt.y[:2,-1],rtol=1e-8,atol=2e-11)
    assert abs(energy[-1]*math.exp(6)/.09375-1)<.03
    assert meta["rhs_evaluations"]<30000


def test_excited_parent_transfers_state_without_resetting_energy():
    model=ScalarFlux()
    p=model.prepare(parent_displacement=.03,parent_velocity=.4)
    assert p["inherited_velocity"]==.4
    assert p["displacement"]==pytest.approx(.03-model.mu*.5/model.omega**2)
    assert p["parent_total_energy"]>p["parent_vacuum"]
    assert abs(p["energy_balance_residual"])<1e-14
    assert p["actual_wall_pressure"]==pytest.approx(p["released_to_transition"])
