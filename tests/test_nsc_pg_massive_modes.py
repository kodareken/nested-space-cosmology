"""Focused checks of the newly assembled massive spatial mode map."""
import copy
import numpy as np
import pytest

from recursive_horizons.nsc_unruh_state import ParentDirac,outgoing_ratio_series
from recursive_horizons.nsc_massive_jost_modes import outgoing_ratio,solve_jost
from recursive_horizons.nsc_paired_horizon_preparation import PairedHorizonSeedMap,source_covariance
from recursive_horizons.nsc_pg_massive_modes import MassivePGModeResolution,open_source_projector,physical_sewing

P=PairedHorizonSeedMap(1.9006916054701435,.23832579963401956,3.973074368754331,1e-10,1e-9)

@pytest.fixture(scope='module')
def open_section():
    j=solve_jost(P.background,3.2,np.pi/2,np.sqrt(5))
    return MassivePGModeResolution(P).section(3.2,np.pi/2,np.sqrt(5),[-.7,0.,.8],[2.2,3.2],jost=j)

def test_massive_asymptotic_owner_preserves_existing_zero_mass_series():
    assert abs(outgoing_ratio(.7,0.,np.sqrt(5),200)[0]-outgoing_ratio_series(.7,np.sqrt(5),200))<3e-14

def test_fiber_rank_uses_parent_gap_not_occupation_or_small_transmission():
    energy=2.;mass=np.pi/2
    assert np.trace(open_source_projector(energy,mass))==3
    assert source_covariance(energy,P.surface_gravity,P.inheritance_ratio,mass)[2,2]==0
    _,_,future,projection=physical_sewing(energy,mass,np.sqrt(1-1e-12),1e-12)
    assert np.linalg.norm(future.conj().T@future-projection)<3e-11
    assert np.trace(open_source_projector(.2,mass))==2
    with pytest.raises(ValueError):physical_sewing(.2,mass,.8,.36)

def test_global_current_and_spectral_jump(open_section):
    assert max(open_section.residuals.values())<3e-9
    assert open_section.spectral_jump_residual()['spectral_jump_relative']<3e-9

def test_dropping_open_infinity_channel_fails_the_spatial_jump(open_section):
    broken=copy.deepcopy(open_section)
    broken.interior[:,:,2]=0.;broken.exterior[:,:,2]=0.
    assert broken.spectral_jump_residual()['spectral_jump_relative']>.1

def test_closed_channel_has_two_sources_and_stable_decaying_field():
    j=solve_jost(P.background,.35,np.pi,8.)
    assert j.transmission==0
    assert np.linalg.norm(j.field(4))<np.linalg.norm(j.field(2.1))
    section=MassivePGModeResolution(P).section(.35,np.pi,8.,[-.7,0.,.8],[2.2,3.2],jost=j)
    assert np.trace(section.source_projector)==2
    assert np.linalg.norm(section.interior[:,:,2])+np.linalg.norm(section.exterior[:,:,2])==0
    assert section.spectral_jump_residual()['spectral_jump_relative']<3e-9

def test_negative_frequency_uses_the_paired_mode_and_source_law():
    j=solve_jost(P.background,3.2,np.pi/2,-np.sqrt(5))
    owner=MassivePGModeResolution(P)
    result=owner.section(-3.2,np.pi/2,np.sqrt(5),[-.5,0.],[2.2],jost=j)
    positive=source_covariance(3.2,P.surface_gravity,P.inheritance_ratio,np.pi/2)
    assert np.linalg.norm(result.source_covariance-(np.eye(3)-positive.conj()))<3e-14
    assert result.angular==np.sqrt(5)
    assert result.energy==-3.2

def test_finite_mode_section_cannot_be_promoted_to_spatial_covariance(open_section):
    with pytest.raises(ValueError,match='global spectral completeness'):
        open_section.spatial_covariance()
    with pytest.raises(TypeError):
        MassivePGModeResolution(P).section(3.,1.,1.,[0.],[3.],seed_covariance=np.eye(2))
