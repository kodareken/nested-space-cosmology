import math
import numpy as np
import pytest
from scipy.integrate import quad
from recursive_horizons.nsc_atlas_volume_observations import VolumeBackground,FlatBackground,fit


def test_time_dependent_lightcone_against_direct_quadrature_and_density_partition():
    bg=VolumeBackground(2.33)
    No=.44;z=np.array([.1,1.,2.33]);DL,DM,DH=bg.distances(No,z,z)
    Ho=bg.state([No])[5][0]
    for i,x in enumerate(z):
        I=quad(lambda N:math.exp(-N)/bg.state([N])[5][0],No-math.log1p(x),No,epsabs=1e-11)[0]
        reference=Ho*math.exp(No)*math.sin(math.sqrt(bg.C)*I)/math.sqrt(bg.C)
        assert DM[i]==pytest.approx(reference,rel=3e-9)
    f=bg.fractions(No)
    assert sum(f[k] for k in ('Omega_Lambda','Omega_m','Omega_k','Omega_volume'))==pytest.approx(1.,abs=1e-14)
    assert f['Omega_volume']>.1
    assert f['Omega_potential']>f['Omega_Lambda']


def test_background_refinement_and_flat_limit_observer_convention():
    z=np.linspace(.03,2.33,20)
    a=VolumeBackground(2.33);b=VolumeBackground(2.33,rtol=2e-12)
    for No in (0.,.44,1.,math.log(20)):
        da=a.distances(No,z,z)[0];db=b.distances(No,z,z)[0]
        assert np.max(abs(5*np.log10(da/db)))<2e-8
    flat=FlatBackground()
    _,DM,DH=flat.distances(.3,z,z)
    assert np.allclose(DH,1/np.sqrt(.3*(1+z)**3+.7))


def test_synthetic_shared_ruler_recovers_fixed_source_epoch():
    bg=VolumeBackground(1.5);truth=.43
    z=np.linspace(.03,1.5,26);DL,_,_=bg.distances(truth,z,z)
    mag=5*np.log10(DL)+23.7;C=np.eye(len(z))*.02**2
    bz=np.array([.3,.8,1.4]);_,DM,DH=bg.distances(truth,bz,bz)
    beta=29.
    labels=['DM_over_rs','DH_over_rs','DV_over_rs']
    shape=np.array([DM[0],DH[1],(bz[2]*DM[2]**2*DH[2])**(1/3)])
    result=fit(bg,(z,z,mag,C,bz,labels,beta*shape,np.eye(3)*.02**2),order=24)
    assert result['joint_best']['theta']==pytest.approx(truth,abs=1e-7)
    assert result['joint_best']['chi_squared_total']<1e-8
    assert result['joint_best']['rd_h_Mpc']==pytest.approx(2997.92458/beta,rel=1e-7)
