"""Exact tail accounting for the same real Fourier profile."""
import pytest
pytest.importorskip('flint')
from flint import arb,acb,acb_poly,ctx
from recursive_horizons.nsc_ks_profile_compression import compress_profiles


def test_removed_modes_enter_each_derivative_tail():
    with ctx.workprec(120):
        original={(1,0):{'coefficients':[acb(1),acb(2),acb(3),acb(2),acb(1)],
                         'tail':(arb(1)/10,arb(1)/20,arb(0))}}
        reduced=compress_profiles(original,1,1)[(1,0)]
        assert reduced['coefficients']==original[(1,0)]['coefficients'][1:4]
        assert acb_poly(reduced['coefficients'])(1)==7
        assert reduced['tail'][0]>=arb(21)/10
        assert reduced['tail'][1]>=(arb(1)/20+8*arb.pi()).upper()
        assert reduced['tail'][2]>=(2*(4*arb.pi())**2).upper()
        assert reduced['discarded_coefficients_bounded']


def test_no_uncertified_extension_and_precision_restoration():
    before=ctx.prec
    profile={(1,0):{'coefficients':[acb(0)]*3,'tail':(arb(0),arb(0))}}
    with pytest.raises(ValueError):compress_profiles(profile,1,2)
    assert ctx.prec==before
