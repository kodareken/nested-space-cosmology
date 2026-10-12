"""Whole-cell geometric majorants contain the same analytic history's point jets."""
from pathlib import Path
import json
import numpy as np
import pytest
pytest.importorskip('flint')
from flint import arb,ctx
from recursive_horizons.nsc_ks_uniform_axial_geometry import uniform_axial_geometry_jets
from recursive_horizons.nsc_ks_spacetime_geometry_jets import spacetime_geometry_jets
from recursive_horizons.nsc_ks_high_radius_derivatives import restore_packed
from recursive_horizons.nsc_ks_source_envelope import axial_center,computational_z_grid
from recursive_horizons.nsc_compatible_history_geometry import CompatibleIncomingMetric,CompatibleRadiusDirection
from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction
from recursive_horizons.nsc_ks_profile_identity import profile_identity

ROOT=Path(__file__).parents[1]


def inputs():
    r=json.loads((ROOT/'results/development/nsc-ks-radius-covariance-bounds.json').read_text())
    center=axial_center();w=LocalAxialFunction((0.,1.),center);u=LocalAxialFunction((0.,),center)
    f=CompatibleIncomingMetric((.001,),(CompatibleRadiusDirection(w,u,.007,.03),))
    return f,r['geometry'],computational_z_grid(2048,length=205/512)


def test_whole_axial_cell_contains_plateau_transition_and_exterior_geometry():
    f,record,grid=inputs()
    with ctx.workprec(160):
        W=[restore_packed(v) for v in record['W_fourier_l1']]
        box=uniform_axial_geometry_jets(f,arb((259,-8),(1,-20)),W,[arb(0)]*10,
            period_left=float(grid[0]),period_length=205/512,axial_profile_identity=record['axial_profile_identity'])
        for offset in (-.1,-.05,-.015,.015,.05,.1):
            point=spacetime_geometry_jets(f,259/256,axial_center()+offset,order=9)
            for t,z in ((0,0),(1,1),(1,4),(4,5),(9,0)):
                assert box.derivative('r',t,z).contains(point.derivative('r',t,z))
        assert box.history_identity==profile_identity(f)


def test_foreign_profile_or_missing_U_bounds_cannot_be_silently_used():
    f,record,grid=inputs()
    with ctx.workprec(160):
        W=[restore_packed(v) for v in record['W_fourier_l1']]
        args=dict(period_left=float(grid[0]),period_length=205/512,axial_profile_identity='foreign')
        with pytest.raises(ValueError,match='different analytic'):
            uniform_axial_geometry_jets(f,1.01,W,[0]*10,**args)
        args['axial_profile_identity']=record['axial_profile_identity']
        with pytest.raises(ValueError,match='all physical profile'):
            uniform_axial_geometry_jets(f,1.01,W,[],**args)


def test_periodic_cell_must_contain_the_original_compact_profile():
    f,record,grid=inputs()
    with ctx.workprec(160):
        W=[restore_packed(v) for v in record['W_fourier_l1']]
        with pytest.raises(ValueError,match='compact profile support'):
            uniform_axial_geometry_jets(f,1.01,W,[0]*10,period_left=axial_center()-.01,
                period_length=.02,axial_profile_identity=record['axial_profile_identity'])
