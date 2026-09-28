"""Independent alias estimates and continuous compact-profile consistency."""
import pytest
pytest.importorskip('flint')
from flint import arb,ctx

from recursive_horizons.nsc_ks_profile_fourier_bound import alias_and_tail_bounds,enclose_profile_fourier
from recursive_horizons.nsc_ks_ball_operator import AnalyticRadiusFamily
from recursive_horizons.nsc_ks_source_envelope import computational_z_grid
from test_nsc_ks_difference_envelope import inputs


def test_alias_and_derivative_tail_majorize_direct_power_sums():
    with ctx.workprec(120):
        p,M,K=8,512,32
        alias,tails=alias_and_tail_bounds(1,1,p,M,K,2)
        factor=(1/(2*arb.pi()))**p
        actual=sum((factor/arb(abs(K+l*M))**p for l in range(-100,101) if l),arb(0))
        assert alias>=actual.upper()
        for j in range(3):
            partial=sum((2*factor*(2*arb.pi())**j*arb(k)**(j-p) for k in range(K+1,2001)),arb(0))
            assert tails[j]>=partial.upper()


def test_compact_profile_dft_retains_reality_and_monotone_bounds():
    model=AnalyticRadiusFamily(inputs()[2])
    mesh=computational_z_grid(16)
    with ctx.workprec(120):
        results=enclose_profile_fourier(model,[(1,0)],float(mesh[0]),float((mesh[1]-mesh[0])*len(mesh)),
            derivative_order=4,transition_panels=16,quadrature_points=128,retained_index=16,max_derivative=1)
        r=results[(1,0)]
        assert r.derivative_l1_upper>0 and r.alias_error_upper>0
        assert r.coefficients[16].imag.contains(0)
        for k in range(1,17):
            assert r.coefficients[16+k].overlaps(r.coefficients[16-k].conjugate())
        assert len(r.uniform_tail_bounds)==2
        assert r.derivative_cells>=36
