import itertools
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_flux_transition import (
    FluxLaw, KottlerRegion, ChargedWall, DustFluxParent, reference_geometry, evolve_reference,
)
from recursive_horizons.nsc_atlas_parent_observables import DustParent


def test_four_form_stress_by_component_contraction():
    metric = np.diag([-1.,1.,1.,1.])
    F = np.zeros((4,)*4)
    flux = 1.5
    for p in itertools.permutations(range(4)):
        sign = (-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))
        F[p] = sign*flux
    raised3 = np.einsum("ab,cd,ef,nbdf->nace",metric,metric,metric,F)
    contraction = np.einsum("mabc,nabc->mn",F,raised3)
    square = np.einsum("mn,mn",metric,contraction)
    stress = contraction/6-metric*square/48
    np.testing.assert_allclose(stress,-flux**2*metric/2,atol=1e-14)
    assert square == -24*flux**2


def test_composition_and_parent_independent_second_difference():
    for bare,parent in [(0.,1.5),(-7.,3.2)]:
        law=FluxLaw(bare_lambda=bare)
        rows=[law.inherit(parent,j) for j in range(4)]
        assert law.inherit(rows[1]["flux"],2)==dict(rows[3],generation=2)
        for j in range(2):
            assert rows[j+2]["Lambda"]-2*rows[j+1]["Lambda"]+rows[j]["Lambda"] == pytest.approx(law.kappa*law.charge**2)


def test_shared_middle_region_and_both_unsquared_junctions():
    law,outer,inner,parent=reference_geometry()
    assert outer.inside == inner.outside
    assert parent.region == outer.outside
    for wall,r in [(outer,.8),(inner,.4)]:
        g=wall.geometry(r,0.)
        assert g["beta_outside"]>0
        for key in ("angular_residual","temporal_residual","force_residual","energy_residual"):
            assert abs(g[key])<2e-13
        bm,bp=g["beta_inside"],g["beta_outside"]
        separate=(-law.k-wall.outside.f_prime(r)/(2*bp)+wall.inside.f_prime(r)/(2*bm))/(1/bp-1/bm)
        assert separate==pytest.approx(g["acceleration"],abs=1e-12)


def test_full_shared_history_and_parent_clock():
    _,outer,inner,parent=reference_geometry()
    solution,_=evolve_reference()
    for y in solution.sol(np.linspace(0,.25,17)).T:
        assert .15<y[4]<=.4 and .8<=y[0]<1.05
        assert y[0]<parent.state(y[8])["radius"]
        for wall,j in [(outer,0),(inner,4)]:
            g=wall.geometry(y[j],y[j+1])
            assert max(abs(g[k]) for k in ("angular_residual","temporal_residual","force_residual","energy_residual"))<1e-9
    y=solution.y[:,-1]
    assert abs(y[2]-y[6])>.01  # equal wall proper times would not be the shared slice
    assert parent.static_time(y[8])==pytest.approx(y[3],abs=1e-11)


def test_cosmological_join_full_darmois_and_friedmann():
    *_,parent=reference_geometry()
    for tau in (0.,.05,.11):
        row=parent.state(tau)
        for key in ("compensation_residual","angular_residual","temporal_residual","Friedmann_residual"):
            assert abs(row[key])<1e-13
        assert row["deceleration"]<0
        assert row["vacuum_from_expansion"]==pytest.approx(parent.region.Lambda,abs=1e-14)
        assert row["dust_tidal_trace"]==pytest.approx(parent.region.Lambda-row["kappa_rho_m"]/2,abs=1e-14)
        assert row["dust_tidal_trace"]<parent.region.Lambda


def test_dust_limit_recovers_existing_parent():
    parent=DustFluxParent(KottlerRegion(.4,1e-14),4.)
    old=DustParent(.4,4.)
    for t in (0.,.3,1.):
        a,b=parent.state(t),old.state(t)
        assert a["H"]==pytest.approx(b["H"],abs=3e-14)
        # Lambda is small but nonzero: retain its O(Lambda) radius correction.
        assert a["radius"]==pytest.approx(b["radius"],abs=1e-13)
        assert a["kappa_rho_m"]==pytest.approx(8*math.pi*b["density"],abs=3e-14)


def test_tidal_observable_removes_mass_not_vacuum():
    for mass in (0.,.01,.03):
        region=KottlerRegion(mass,.5)
        for r in (.3,.6,.9):
            assert region.tidal(r)["trace"]==pytest.approx(.5,abs=2e-15)
    law=FluxLaw(charge=0.)
    inner=ChargedWall.from_inside(law,1.5,0.,.4)
    outer=ChargedWall.from_inside(law,1.5,inner.outside.mass,.8)
    traces=[region.tidal(r)["trace"] for region,r in
            zip((outer.outside,outer.inside,inner.inside),(1.1,.6,.1))]
    np.testing.assert_allclose(traces,[1.125]*3,atol=3e-15)


def test_physical_unit_scaling():
    law,outer,inner,_=reference_geometry()
    s=7.
    scaled=FluxLaw(kappa=s*s,bare_lambda=0.,charge=law.charge/s**2,tension=law.tension/s**3)
    inside=ChargedWall.from_inside(scaled,1/s**2,0.,.4*s)
    outside=ChargedWall.from_inside(scaled,1.5/s**2,inside.outside.mass,.8*s)
    for original,new,r in [(inner,inside,.4),(outer,outside,.8)]:
        assert new.outside.mass/s==pytest.approx(original.outside.mass,abs=1e-15)
        assert new.inside.Lambda*s*s==pytest.approx(original.inside.Lambda)
        assert new.acceleration(r*s)*s==pytest.approx(original.acceleration(r),abs=1e-12)


def test_changed_vacuum_cannot_be_spliced_into_fixed_action():
    law,outer,_,_=reference_geometry()
    with pytest.raises(ValueError,match="flux law"):
        ChargedWall(law,outer.outside_flux,KottlerRegion(outer.inside.mass,.51),outer.outside)
    with pytest.raises(ValueError):
        ChargedWall.from_inside(law,1.5,0.,10.)
    with pytest.raises(ValueError):
        evolve_reference(end=.5)


def test_negative_mass_trial_cannot_be_compensated_by_positive_dust():
    law,_,inner,_=reference_geometry()
    trial=ChargedWall.from_inside(law,1.5,inner.outside.mass,1.)
    assert trial.outside.mass<0
    with pytest.raises(ValueError,match="positive mass"):
        DustFluxParent(trial.outside)
