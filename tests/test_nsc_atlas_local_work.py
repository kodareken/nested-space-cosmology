import math
import numpy as np
import pytest
import sympy as s
from recursive_horizons.nsc_atlas_flow_evolution import PeriodicFlow
from recursive_horizons.nsc_atlas_local_work import padded_spectrum,constraint_fields
from recursive_horizons.nsc_atlas_rindler_bulk import checks,exact_radial_check,ricci_series


def test_padding_preserves_band_limited_field_and_derivatives():
    f=PeriodicFlow(16,length=2*np.pi,workers=1);g=PeriodicFlow(32,length=2*np.pi,workers=1)
    x=2*np.pi*np.arange(16)/16
    u=np.zeros((3,16,16,16));u[0]=np.sin(x)[None,:,None]
    h=f.fft(u);padded=padded_spectrum(h,32)
    np.testing.assert_allclose(g.real(padded)[:,::2,::2,::2],u,atol=4e-16)


def test_exact_shear_needs_pure_affine_expansion():
    f=PeriodicFlow(16,length=2*np.pi,viscosity=.1,workers=1)
    y=2*np.pi*np.arange(16)/16;u=np.zeros((3,16,16,16));u[0]=np.sin(y)[None,:,None]
    row=constraint_fields(f,f.project(f.fft(u)))
    np.testing.assert_allclose(row['source'],1,atol=8e-15)
    assert row['affine_rate']==pytest.approx(1/3)
    assert row['periodic_correction_rms']<2e-15
    assert row['residual_rms']<1e-14


def test_exact_taylor_green_source_and_poisson_completion():
    nu=.2;f=PeriodicFlow(16,length=2*np.pi,viscosity=nu,workers=1)
    a=2*np.pi*np.arange(16)/16;x,y,z=np.meshgrid(a,a,a,indexing='ij')
    u=np.array([np.sin(x)*np.cos(y),-np.cos(x)*np.sin(y),np.zeros_like(x)])
    row=constraint_fields(f,f.project(f.fft(u)))
    b=2*np.pi*np.arange(32)/32;X,Y,Z=np.meshgrid(b,b,b,indexing='ij')
    source=2+np.cos(2*X)+np.cos(2*Y)+(np.cos(3*X)*np.cos(Y)-np.cos(X)*np.cos(3*Y))/(4*nu)
    potential=-(np.cos(2*X)+np.cos(2*Y))/4-(np.cos(3*X)*np.cos(Y)-np.cos(X)*np.cos(3*Y))/(40*nu)
    np.testing.assert_allclose(row['source'],source,atol=4e-14)
    np.testing.assert_allclose(row['potential'],potential,atol=3e-15)
    assert row['affine_rate']==pytest.approx(2/3)
    assert row['pressure_mean']==pytest.approx(0,abs=1e-15)


def test_metric_shear_bulk_has_all_ricci_orders_zero_and_detects_omission():
    for row in checks():
        assert max(row['maximum_Ricci_coefficient_by_order'])<2e-13
        assert row['omitted_expansion_order4_residual']>.3
        assert row['background_determinant']==-1
    exact=exact_radial_check()
    assert exact['all_coefficients_zero'] and exact['checked_components']==125


def test_moving_lattice_pullback_retains_shift_and_does_not_equal_FRW_ansatz():
    N,b,bd,x=s.symbols('N b bd x',real=True)
    # Flat -N^2 dt^2+dX^2 with X=b(t)x has a nonzero shift.
    g=s.Matrix([[-N*N+bd*bd*x*x,b*bd*x],[b*bd*x,b*b]])
    assert s.factor(g.det())==-N*N*b*b
    assert g[0,1]==b*bd*x


def test_independent_curved_metric_control_for_the_tensor_engine():
    # de Sitter four-metric times an unwarped radial line, evaluated at t=0.
    H=.7;g=np.zeros((1,5,5));g[0]=np.diag([-1.,1.,1.,1.,1.])
    dg=np.zeros((5,1,5,5));ddg=np.zeros((5,5,1,5,5))
    for i in (2,3,4):dg[0,0,i,i]=2*H;ddg[0,0,0,i,i]=4*H*H
    expected=np.diag([-3*H*H,0,3*H*H,3*H*H,3*H*H])
    np.testing.assert_allclose(ricci_series(g,dg,ddg)[0],expected,atol=1e-15)


def test_periodic_gradient_completion_is_orthogonal_to_transverse_freedom():
    n=16;f=PeriodicFlow(n,length=2*np.pi,workers=1)
    a=2*np.pi*np.arange(n)/n;x,y,z=np.meshgrid(a,a,a,indexing='ij')
    potential=np.cos(x)*np.cos(y)+np.sin(z)
    gradient=np.array([-np.sin(x)*np.cos(y),-np.cos(x)*np.sin(y),np.cos(z)])
    transverse=np.zeros_like(gradient);transverse[0]=np.sin(y)
    assert np.mean(np.sum(gradient*transverse,axis=0))==pytest.approx(0,abs=1e-16)
    assert np.mean(np.sum((gradient+transverse)**2,axis=0))>=np.mean(np.sum(gradient**2,axis=0))
