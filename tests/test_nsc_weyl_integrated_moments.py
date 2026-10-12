"""Whole-cell and signed-momentum guards plus an exact zero-defect control."""
import pytest
pytest.importorskip('flint')
from flint import arb,acb,ctx
from recursive_horizons.nsc_scaled_reference_projector import scaled_reference_projector
from recursive_horizons.nsc_global_projector_bounds import global_projector_bounds
from recursive_horizons.nsc_radius_transfer_sums import _catalog_for_cutoff,_spec_key
from recursive_horizons.nsc_weyl_integrated_moments import integrated_defect_moments,matrix_nuclear_upper
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper


def setup():
    with ctx.workprec(160):
        scaled={s:scaled_reference_projector(mu=arb((1,-7),(1,-7)),sign=s,mass=1,angular=2,axial=1,radius=2) for s in (-1,1)}
        glob=global_projector_bounds(mass=1,angular=2,axial=1,radius=2)
        cat={'canonical_cutoff_lower':arb(256),'catalog':{_spec_key(s):{'total':arb(0)} for s in _catalog_for_cutoff(arb(256),160)}}
        return scaled,glob,cat


def test_nuclear_matrix_bound_keeps_the_two_singular_values():
    with ctx.workprec(160):
        assert matrix_nuclear_upper([acb(1),acb(0),acb(0),acb(1)])==2
        assert matrix_nuclear_upper([acb(2),acb(0),acb(0),acb(0)])==2


def test_both_auxiliary_choices_keep_zero_constant_defect_and_scope():
    scaled,glob,cat=setup()
    for cut in (False,True):
        r=integrated_defect_moments(scaled,glob,1,.4,cat,256,low_momentum_auxiliary=cut,uniform_axial_coverage=True)
        with ctx.workprec(160):
            assert restored_upper(r['zeroth_moment_density_upper'])==0
            assert restored_upper(r['first_moment_density_upper'])==0
        assert r['physical_local_gate']=='OPEN'
        assert not r['physical_subtraction_changed']
        assert not r['source_energy_tail_identified_with_canonical_cut']
        assert r['rho_integral_included'] is False


def test_local_box_cannot_be_relabelled_as_a_fourier_moment_bound():
    scaled,glob,cat=setup()
    with pytest.raises(ValueError,match='whole-cell'):
        integrated_defect_moments(scaled,glob,1,.4,cat,256)
    with pytest.raises(ValueError,match='both canonical'):
        integrated_defect_moments({1:scaled[1]},glob,1,.4,cat,256,uniform_axial_coverage=True)
