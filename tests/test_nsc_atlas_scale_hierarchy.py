import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_scale_hierarchy import ScaleHierarchy,depletion_integral


def test_exact_core_polynomial_recovers_original_constraint_integrator():
    model=ScaleHierarchy();row,_=model.core_data()
    original=model.core.mass_solution(model.reference.expansions[1])[0].y[0,-1]
    assert row['mass']==pytest.approx(original,rel=2e-9)
    w=model.frequency;z=model.weak_amplitude;c=row['coefficients']
    assert row['mass']==pytest.approx(z*(c[0]/w+c[1]/w**2+c[2]/w**3),rel=2e-15)


def test_resolved_reservoir_matches_independent_hamilton_geometry():
    r=ScaleHierarchy(frequency=1e9,weak_amplitude=.001).solve()
    assert abs(r['reservoir_constraint']['endpoint_mass_ratio'])<1e-9
    assert abs(r['independent_proper_distance']['endpoint_mass_ratio'])<1e-9
    assert r['reservoir_outer_radius']>r['ideal_capacity_lower_radius']>r['core']['cutoff_radius']
    assert r['parent_vacuum']==.28125
    assert r['protected_parent_angles'][0]<r['protected_parent_angles'][1]


def test_fixed_radius_capacity_obstruction_cannot_be_removed_by_units():
    m=ScaleHierarchy(frequency=8000.);r,_=m.core_data()
    available=m.reference.scalar.kappa*m.excitation*7*m.radius**3/6
    assert r['mass']>available
    # Length scaling: geometric mass scales as L, density as L^-2.
    for scale in (1e-15,1.,1e25):
        rescaled=m.reference.scalar.kappa*(m.excitation/scale**2)*7*(m.radius*scale)**3/6
        assert r['mass']*scale/rescaled==pytest.approx(r['mass']/available,rel=3e-15)


def test_reservoir_shape_dependence_remains_explicit():
    half=depletion_integral(.5)
    assert half==pytest.approx(5476736/66927861,rel=2e-15)
    assert depletion_integral(.25)<half<depletion_integral(.75)
    assert 1.5/half==pytest.approx(18.330588054636923,rel=2e-15)


def test_core_coefficients_are_independent_of_frequency_at_fixed_parent():
    lo=ScaleHierarchy();hi=ScaleHierarchy(frequency=1e8)
    a,_=lo.core_data();b,_=hi.core_data()
    np.testing.assert_allclose(a['coefficients'],b['coefficients'],rtol=2e-10,atol=2e-10)
    assert lo.reference.densities==pytest.approx(hi.reference.densities,abs=2e-15)
    assert lo.reference.expansions==pytest.approx(hi.reference.expansions,abs=2e-15)


def test_retained_hierarchy_has_checked_sources_and_finite_error_scope():
    root=Path(__file__).resolve().parents[1]
    r=json.loads((root/'results/summaries/nsc-atlas-scale-hierarchy-v1.json').read_text(encoding='utf-8'))
    assert not r['numerical_validation_errors'] and not r['budget_exhausted']
    for name,digest in r['source_sha256'].items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest,name
    assert 6000<r['fixed_ratio_obstruction']['necessary_frequency_ceiling']<6500
    exemplar=r['illustrative_physical_conversion']
    assert exemplar['assumed_Lambda_per_m2']==1e-52
    assert exemplar['core_half_mass_radius_kpc']>0
