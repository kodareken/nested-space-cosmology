import numpy as np
import pytest
from recursive_horizons.nsc_atlas_bao_shape import ap_prediction,ap_null_contrast


def test_ap_known_limits_and_common_ruler_cancellation():
    z=np.array([.1,.51,2.33])
    np.testing.assert_allclose(ap_prediction(z,matter=0),z,rtol=1e-14)
    np.testing.assert_allclose(ap_prediction(z,matter=1),2*((1+z)**1.5-(1+z)),rtol=1e-14)
    F=ap_prediction(z,matter=.33)
    pairs=np.array([[0,1],[2,3],[4,5]])
    y=np.array([F[0]*10,10,F[1]*8,8,F[2]*5,5])
    L=np.diag([.2,.1,.3,.2,.4,.3])+np.tril(np.ones((6,6))*.01,-1)
    C=L@L.T
    assert ap_null_contrast(y,C,pairs,F)['chi_squared']<1e-25
    changed=y.copy();changed[2]+=.4
    base=ap_null_contrast(changed,C,pairs,F)
    # Independent constrained-mean fit: unknown radial amplitude in each
    # bin, with transverse mean F times that amplitude. This must give the
    # same minimum chi-square as the linear-null projection.
    design=np.zeros((6,3))
    design[pairs[:,0],np.arange(3)]=F
    design[pairs[:,1],np.arange(3)]=1.
    inverse=np.linalg.inv(C)
    fitted=np.linalg.solve(design.T@inverse@design,design.T@inverse@changed)
    residual=changed-design@fitted
    assert base['chi_squared']==pytest.approx(residual@inverse@residual,rel=1e-13)
    scaled=ap_null_contrast(7*changed,49*C,pairs,F)
    assert scaled['chi_squared']==pytest.approx(base['chi_squared'],rel=1e-14)
    assert base['chi_squared']>0


def test_ap_refuses_duplicate_rows_and_nonfinite_data():
    with pytest.raises(ValueError):
        ap_null_contrast([1.,2.],np.eye(2),[[0,0]],[1.])
    with pytest.raises(ValueError):
        ap_null_contrast([np.nan,2.],np.eye(2),[[0,1]],[1.])
