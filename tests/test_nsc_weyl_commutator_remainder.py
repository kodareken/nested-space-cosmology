"""Independent Fourier-matrix and polynomial checks of the exact Weyl defect."""
import math
import numpy as np
import pytest
pytest.importorskip('flint')
from flint import arb,ctx
from recursive_horizons.nsc_weyl_commutator_remainder import (
    kinetic_commutator,potential_commutator,potential_taylor_commutator,
    rho_operator_defect,shifted_taylor_remainder,large_transfer_remainder,
)
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper

S1=np.array([[0,1],[1,0]],complex)
S2=np.array([[0,-1j],[1j,0]],complex)
S3=np.diag([1,-1]).astype(complex)
ZERO=np.zeros((2,2),complex)


def symbol(n,k):
    if n==0:return S3*(.2+k*k/7)+S1*k/3
    if abs(n)==1:return S1*(1+k*k/11)+S2*(.1*k)
    return ZERO


@pytest.mark.parametrize('q',[-np.pi/.4,.19,np.pi/.4])
def test_weyl_commutator_matches_independent_bloch_matrix_product(q):
    length=.4;C=.7*S3;V={0:.3*S1,2:(.02+.01j)*S2,-2:(.02-.01j)*S2}
    modes=np.arange(-8,9);wave=q+2*np.pi*modes/length
    H=np.zeros((2*len(modes),)*2,complex);P=np.zeros_like(H)
    for i,m in enumerate(modes):
        for j,n in enumerate(modes):
            H[2*i:2*i+2,2*j:2*j+2]=V.get(int(m-n),ZERO)+(C*wave[i] if i==j else ZERO)
            P[2*i:2*i+2,2*j:2*j+2]=symbol(int(m-n),(wave[i]+wave[j])/2)
    comm=H@P-P@H
    for i,m in enumerate(modes):
        for j,n in enumerate(modes):
            if abs(m)>3 or abs(n)>3:continue
            k=(wave[i]+wave[j])/2
            expected=kinetic_commutator(C,symbol(int(m-n),k),int(m-n),k,length)
            expected+=potential_commutator(V,symbol,int(m-n),k,length)
            np.testing.assert_allclose(comm[2*i:2*i+2,2*j:2*j+2],expected,atol=2e-10,rtol=3e-13)


def test_exact_shift_is_the_owned_star_taylor_for_a_polynomial():
    V={-2:.03*S1,1:.05*S3};period=2*np.pi;k=.37
    def p(n,k):return S2*k**4 if n==0 else ZERO
    def d(n,k,j):return S2*math.factorial(4)/math.factorial(4-j)*k**(4-j) if n==0 and j<=4 else ZERO
    for n in (-2,1):
        exact=potential_commutator(V,p,n,k,period)
        truncated=potential_taylor_commutator(V,d,n,k,period,4)
        np.testing.assert_allclose(exact,truncated,atol=3e-16,rtol=0)
    assert np.max(abs(potential_commutator(V,p,1,k,period)-potential_taylor_commutator(V,d,1,k,period,0)))>.01


def test_small_transfer_factor_and_both_sides_are_required():
    # p(k)=k^3*I, V=S1, k=0, shifts+/-s: commutator=-2*s^3*S1.
    s=.3
    r=shifted_taylor_remainder(s,2,1,6*np.sqrt(2),6*np.sqrt(2))
    with ctx.workprec(120):
        actual=arb(2)*arb(s)**3*arb(2).sqrt()
        bound=restored_upper(r['remainder_Frobenius_upper'])
        assert bound>=actual
        assert bound<actual*(1+arb('1e-14'))
    with pytest.raises(ValueError,match='explicit'):
        shifted_taylor_remainder(s,2,1,None,1)


def test_large_transfer_bound_preserves_missing_moments():
    r=large_transfer_remainder([2,3],[5,7],11,13)
    with ctx.workprec(120):assert restored_upper(r['remainder_Frobenius_upper'])==arb(110)
    with pytest.raises(ValueError,match='matching'):
        large_transfer_remainder([2],[3,4],1,1)
    assert r['physical_local_gate']=='OPEN'


def test_rho_defect_reproduces_covariance_duhamel_with_decreasing_rho():
    from scipy.linalg import expm
    H=.7*S1+.3*S3;P=.5*np.eye(2)+.2*S2
    C=np.array([[.4,.1j],[-.1j,.6]])
    start,end=.03,0.;U=expm(1j*H*(end-start))
    remainder=rho_operator_defect(ZERO,H@P-P@H)
    x,w=np.polynomial.legendre.leggauss(16);integral=ZERO.copy()
    for node,weight in zip(x,w):
        t=(start+end)/2+(end-start)*node/2
        step=expm(1j*H*(end-t))
        integral+=(end-start)*weight/2*(step@remainder@step.conj().T)
    direct=U@C@U.conj().T-P
    homogeneous=U@(C-P)@U.conj().T
    np.testing.assert_allclose(direct,homogeneous-integral,atol=3e-16,rtol=0)
    assert np.max(abs(direct-(homogeneous+integral)))>.001
