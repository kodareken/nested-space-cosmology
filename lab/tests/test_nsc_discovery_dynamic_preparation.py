"""Owned general-Q Jacobian, source, exact constraints and checkpoint checks."""
from pathlib import Path
import tempfile

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

from recursive_horizons import nsc_discovery_dynamic_preparation as prep
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons import nsc_spherical_coupling as coupling


@pytest.fixture(scope="module")
def prepared():
    with threadpool_limits(limits=1):
        return prep.build_record(execute=True,nf=64,harmonic=2,cpu_limit=20.)


def test_general_q_jacobian_is_owned_centered_difference():
    with threadpool_limits(limits=1):
        pair=prep.make_pair(32);grid=pair.grid
        Q,_=prep.declared_q(grid,harmonic=2)
        state,_=prep.spectral_state(pair,Q)
        state.r=1.4+.15*np.cos(2*np.pi*grid.xi_g/grid.length)
        fine=galerkin.prolong_state(grid,state)
        source=coupling.source_from_columns(galerkin.active_fine_system(grid,fine),fine)
        rho=source["force_L"]/grid.dx_q
        _,_,J=prep.radius_residual_jacobian(grid,state,rho)
        direction=np.cos(4*np.pi*grid.xi_g/grid.length)+.3
        images=[]
        for sign in (-1,1):
            varied=state.copy();varied.r+=sign*1e-6*direction
            images.append(prep.radius_residual_jacobian(grid,varied,rho)[0])
        np.testing.assert_allclose((images[1]-images[0])/2e-6,J@direction,rtol=2e-7,atol=2e-6)


def test_convex_seed_then_exact_constraint_and_current(prepared):
    report,arrays=prepared
    assert report["status"]=="FINITE_INITIAL_DATA"
    assert report["source_eta"]==1. and report["source_trace"]==3.
    for case in ("uniform","pattern"):
        measured=report["cases"][case]
        assert measured["converged"]
        assert measured["full_lapse_max"]<1e-6
        assert measured["full_shift_max"]<1e-10
        assert abs(measured["source_current_mean"])<1e-10
        assert measured["rho_min"]>0 and measured["r_min"]>0
        assert measured["spectral"]["commutator_max"]<1e-10
        assert not measured["holding_claim"] and not measured["evolution_performed"]
        assert np.all(arrays[case+"_nodal_chi"]==0)
        for key in ("p_Q","p_r","p_chi"):
            assert np.all(arrays[case+"_nodal_"+key]==0)
        phi=np.vstack((arrays[case+"_nodal_phi0"],arrays[case+"_nodal_phi1"]))
        C=(phi*prep.WEIGHTS)@phi.conj().T
        ev=np.linalg.eigvalsh(C)
        assert np.min(ev)>-1e-12 and np.max(ev)<1+1e-12
        np.testing.assert_allclose(ev[-6:],np.sort(prep.WEIGHTS),atol=1e-12)
    assert report["cases"]["uniform"]["declaration"]["mean_Q"]==report["cases"]["pattern"]["declaration"]["mean_Q"]
    assert not report["evidence_written"]


def test_preflight_and_checkpoint_are_explicit(prepared):
    report,arrays=prep.build_record()
    assert report["status"]=="PREPARATION_PREFLIGHT" and arrays=={}
    report,arrays=prepared
    with tempfile.TemporaryDirectory(dir="/tmp",prefix="nsc-dynamic-prep-") as directory:
        stem=Path(directory)/"checkpoint"
        prep.write_record(report,arrays,stem)
        before={p.name:p.read_bytes() for p in Path(directory).iterdir()}
        checked=prep.check_record(stem)
        assert checked["ok"] and checked["bytes_written"]==0 and not checked["solver_executed"]
        assert before=={p.name:p.read_bytes() for p in Path(directory).iterdir()}
        pair,state,_=prep.load_case(stem)
        assert pair.grid.nf==64 and state.phi0.shape==(64,6)
        with pytest.raises(FileExistsError):
            prep.write_record(report,arrays,stem)


def test_budget_and_positive_input_domain():
    with pytest.raises(ValueError,match="CPU"):
        prep.build_record(cpu_limit=31.)
    pair=prep.make_pair(32)
    with pytest.raises(ValueError,match="positive lift"):
        prep.declared_q(pair.grid,lift=0)
