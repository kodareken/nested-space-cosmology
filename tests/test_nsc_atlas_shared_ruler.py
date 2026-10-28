import numpy as np
import pytest
from scipy.integrate import quad
from scipy.linalg import solve

from recursive_horizons.nsc_atlas_shared_ruler import SharedRuler,bao_shape,marginal_score,C_OVER_100


def test_one_ruler_recovers_correlated_data_and_independent_normal_equation():
    z=np.array([.295,.51,.51,.706,.706])
    labels=['DV_over_rs','DM_over_rs','DH_over_rs','DM_over_rs','DH_over_rs']
    t=bao_shape('flat_LCDM',.33,z,labels)
    L=np.diag([.1,.2,.3,.15,.25])+np.tril(np.ones((5,5))*.01,-1)
    C=L@L.T;y=30*t+np.array([.01,-.03,.02,.04,-.02])
    answer=SharedRuler(y,C).profile(t)
    precision=solve(C,np.eye(5),assume_a='sym')
    beta=t@precision@y/(t@precision@t)
    residual=y-beta*t
    assert answer['beta']==pytest.approx(beta,rel=2e-15)
    assert answer['chi_squared']==pytest.approx(residual@precision@residual,abs=2e-13)
    assert answer['rd_h_Mpc']==pytest.approx(C_OVER_100/beta)


def test_shared_ruler_uses_one_amplitude_for_all_bins_and_all_types():
    z=np.array([.295,.51,.51,2.33,2.33])
    labels=['DV_over_rs','DM_over_rs','DH_over_rs','DH_over_rs','DM_over_rs']
    t=bao_shape('parent',np.log(1.38),z,labels)
    y=30*t
    assert SharedRuler(y,np.eye(5)*.01).profile(t)['chi_squared']<1e-22
    y[-2:]*=1.05
    assert SharedRuler(y,np.eye(5)*.01).profile(t)['chi_squared']>10


def test_marginalization_keeps_shape_uncertainty_and_normalization():
    # A continuous posterior makes an analytic Gaussian convolution.
    beta=.8;obs=1.3;noise=.4;shape_sd=.3
    nodes,weights=np.polynomial.legendre.leggauss(160)
    shape=3*nodes;qw=3*weights
    w=qw*np.exp(-shape*shape/(2*shape_sd*shape_sd));w/=w.sum()
    # mean prediction beta+shape, so the shape dependence enters b and c.
    a=np.full_like(shape,1/noise**2);b=(obs-shape)/noise**2;c=(obs-shape)**2/noise**2
    value=marginal_score(beta,np.log(w),a,b,c)
    direct=quad(lambda x:np.exp(-.5*((obs-beta-x)/noise)**2)*
                np.exp(-.5*(x/shape_sd)**2)/(np.sqrt(2*np.pi)*shape_sd),-3,3,
                epsabs=1e-13)[0]
    assert value==pytest.approx(-2*np.log(direct),abs=2e-13)
    assert abs(value-((obs-beta)/noise)**2)>.01
