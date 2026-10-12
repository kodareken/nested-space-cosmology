"""Bounds used by the retained-source interpolation component, with no solves."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest
pytest.importorskip('flint')
from flint import arb,ctx

spec=importlib.util.spec_from_file_location('interp_accuracy',Path(__file__).parents[1]/'scripts/derive_nsc_ks_energy_interpolation_accuracy.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_source_norm_keeps_stored_nonhermitian_roundoff_in_the_bound():
    # Row infinity norm alone would give1, below the true sqrt(3) norm.
    blocks=np.array([[[1,0,0],[1,0,0],[1,0,0]]],complex)
    with ctx.workprec(160):
        upper=module.gamma_from_blocks(blocks)
        assert upper>=arb(3).sqrt()
        assert upper<arb(3).sqrt()+arb('1e-40')


def test_control_integrals_reject_zero_amplitude_and_foreign_domain():
    import json
    radius=json.loads(module.RADIUS.read_text())
    with ctx.workprec(160):
        values,a,r=module.integrals(radius,np.pi/2,np.sqrt(5),1.03,.03)
        assert all(v>0 for v in values)
        assert a>0 and r>0
        assert values[4]>=values[2]/arb(radius['control']['history_amplitude'])
        with pytest.raises(ValueError,match='slab'):
            module.integrals(radius,1.,2.,1.2,.03)
        radius['control']['history_amplitude']=0.
        with pytest.raises(ValueError,match='single nonzero'):
            module.integrals(radius,1.,2.,1.03,.03)
