"""Independent canonical signs, Hessian and finite-projection controls."""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_boson_hamiltonian import (
    potential_hessian,quadratic_hamiltonian,
)


SSTAR=-.035646782552940987


def metric(N,beta,h):
    N=np.asarray(N);count=len(N);g=np.zeros((count,4,4))
    g[:,1:,1:]=h;g[:,0,1:]=np.einsum('nij,nj->ni',h,beta)
    g[:,1:,0]=g[:,0,1:]
    g[:,0,0]=-N*N+np.einsum('ni,ni->n',beta,g[:,0,1:])
    return g


def simple(nodes=2,modes=2):
    return dict(g=np.broadcast_to(np.diag([-1.,1.,1.,1.]),(nodes,4,4)).copy(),
        background=np.zeros((nodes,2)),Y=np.ones((nodes,modes)),
        DY=np.zeros((nodes,3,modes)),weights=np.ones(nodes)/nodes,
        fields=np.arange(modes)%2,sstar=SSTAR)


def test_hessian_matches_direct_symbolic_differentiation_of_the_action():
    import sympy as s
    S,chi,star=s.symbols('S chi star',real=True)
    size=star+S/s.sqrt(84)
    W=7*s.exp(-12*size)-s.Rational(15,2)*s.exp(-14*size)+s.Rational(25,32)*s.exp(-24*size)+200*s.exp(-12*size)*chi**2
    exact=s.lambdify((S,chi,star),s.hessian(W,(S,chi)),'numpy')
    for background in ([0.,0.],[.03,.02],[-.04,-.05]):
        result=potential_hessian(background,sstar=SSTAR)
        np.testing.assert_allclose(result,exact(*background,SSTAR),rtol=3e-14,atol=3e-14)
        np.testing.assert_array_equal(result,result.T)
    assert potential_hessian([0.,0.],sstar=SSTAR)[0,1]==0.
    assert potential_hessian([0.,.02],sstar=SSTAR)[0,1]<0.
    assert potential_hessian([0.,-.02],sstar=SSTAR)[0,1]>0.


def test_round_sphere_has_known_canonical_kinetic_and_mode_frequencies():
    from recursive_horizons.nsc_atlas_s3_hopf import build_hopf_frame_basis
    basis=build_hopf_frame_basis(1,(4,8));X=basis.points;volume=2*math.pi**2
    y0=np.ones(len(X))/math.sqrt(volume);y1=2*X[:,0]/math.sqrt(volume)
    r=.31;a=1.6;N=1.3
    Y=np.column_stack([y0,y1,y0,y1]);DY=np.zeros((len(X),3,4))
    derivative=2*r*np.column_stack([X[:,3],-X[:,2],X[:,1]])/math.sqrt(volume)
    DY[:,:,1]=derivative;DY[:,:,3]=derivative
    h=np.broadcast_to(a*a*np.eye(3),(len(X),3,3))
    g=metric(np.full(len(X),N),np.zeros((len(X),3)),h)
    out=quadratic_hamiltonian(g,np.zeros((len(X),2)),Y,DY,basis.weights,[0,0,1,1],sstar=SSTAR)
    masses=potential_hessian([0.,0.],sstar=SSTAR).diagonal()
    expectedK=N*a*np.array([0.,3*r*r,0.,3*r*r])+N*a**3*np.repeat(masses,2)
    np.testing.assert_allclose(out['A'],N/a**3*np.eye(4),rtol=0,atol=2e-15)
    np.testing.assert_allclose(out['K'],np.diag(expectedK),rtol=3e-15,atol=3e-13)
    np.testing.assert_array_equal(out['B'],0.)
    assert np.linalg.eigvalsh(out['G']).min()>0
    n=4;J=np.block([[np.zeros((n,n)),np.eye(n)],[-np.eye(n),np.zeros((n,n))]])
    generator=J@out['G']
    np.testing.assert_allclose(generator.T@J+J@generator,0.,atol=2e-15,rtol=0)


def test_shift_orientation_matches_direct_canonical_hamiltonian_gradient():
    z=simple(3,2);z['fields']=np.zeros(2,dtype=int)
    z['Y']=np.array([[1.,2.],[3.,-1.],[.5,1.3]])
    z['DY'][:,0,:]=[[.7,-.4],[-1.2,2.],[.9,.2]]
    beta=np.zeros((3,3));beta[:,0]=[.2,-.3,.5]
    z['g']=metric(np.ones(3),beta,np.broadcast_to(np.eye(3),(3,3,3)))
    out=quadratic_hamiltonian(**z)
    expected=np.zeros((2,2))
    for i in range(3):expected+=z['weights'][i]*np.outer(z['Y'][i],beta[i]@z['DY'][i])
    np.testing.assert_allclose(out['B'],expected,rtol=0,atol=2e-16)
    assert abs(expected[0,1]-expected[1,0])>.1
    q=np.array([.4,-.2]);p=np.array([.3,.7]);state=np.r_[q,p]
    J=np.block([[np.zeros((2,2)),np.eye(2)],[-np.eye(2),np.zeros((2,2))]])
    expected_rate=np.r_[out['A']@p+expected@q,-out['K']@q-expected.T@p]
    np.testing.assert_allclose(J@out['G']@state,expected_rate,rtol=2e-15,atol=3e-15)


def test_direct_momentum_density_projection_is_not_inverse_truncated_mass():
    z=simple(2,1);z['fields']=np.array([0]);z['Y'][:]=1.
    z['g'][1,1:,1:]=4*np.eye(3)
    out=quadratic_hamiltonian(**z)
    direct=.5*(1+1/8)
    inverse_truncated=1/(.5*(1+8))
    assert out['A'][0,0]==pytest.approx(direct,abs=1e-16)
    assert abs(out['A'][0,0]-inverse_truncated)>.3


def test_arbitrary_metric_projection_matches_pointwise_canonical_energy():
    rng=np.random.default_rng(179);z=simple(7,3);z['fields']=np.array([0,1,0])
    seed=rng.normal(scale=.2,size=(7,3,3));h=np.eye(3)+np.einsum('nik,njk->nij',seed,seed)
    N=rng.uniform(.7,1.4,7);beta=rng.normal(scale=.3,size=(7,3))
    z['g']=metric(N,beta,h);z['background']=rng.normal(scale=.02,size=(7,2))
    z['Y']=rng.normal(size=(7,3));z['DY']=rng.normal(size=(7,3,3))
    z['weights']=rng.uniform(.1,.3,size=7)
    out=quadratic_hamiltonian(**z);q=np.array([.4,-.2,.3]);p=np.array([-.7,.1,.5])
    mass=potential_hessian(z['background'],sstar=SSTAR)
    direct=0.
    for n in range(7):
        fields=np.zeros(2);momenta=np.zeros(2);gradients=np.zeros((3,2))
        for a,field in enumerate(z['fields']):
            fields[field]+=z['Y'][n,a]*q[a]
            momenta[field]+=z['Y'][n,a]*p[a]
            gradients[:,field]+=z['DY'][n,:,a]*q[a]
        determinant=math.sqrt(np.linalg.det(h[n]));hi=np.linalg.inv(h[n])
        energy=N[n]*momenta@momenta/(2*determinant)
        energy+=N[n]*determinant*(np.einsum('ij,iA,jA->',hi,gradients,gradients)+fields@mass[n]@fields)/2
        energy+=momenta@(beta[n]@gradients)
        direct+=z['weights'][n]*energy
    state=np.r_[q,p]
    assert .5*state@out['G']@state==pytest.approx(direct,rel=4e-15,abs=2e-15)
    np.testing.assert_allclose(out['lapse'],N,rtol=0,atol=4e-16)
    np.testing.assert_allclose(out['shift'],beta,rtol=0,atol=3e-16)
    np.testing.assert_array_equal(out['G'],out['G'].T)


def test_valid_large_shift_and_negative_potential_curvature_are_not_clipped():
    z=simple(2,2);z['fields']=np.zeros(2,dtype=int);z['Y']=math.sqrt(2)*np.eye(2)
    derivative=np.array([[0.,1.],[-1.,0.]])
    z['DY'][:,0,:]=z['Y']@derivative
    z['g']=metric(np.ones(2),np.tile([3.,0.,0.],(2,1)),np.broadcast_to(np.eye(3),(2,3,3)))
    out=quadratic_hamiltonian(**z)
    assert np.all(z['g'][:,0,0]>0) and np.linalg.eigvalsh(out['G']).min()<0
    z=simple(2,2);z['Y']=np.tile([1.,1.],(2,1));z['background'][:,1]=.2
    out=quadratic_hamiltonian(**z)
    assert out['minimum_nodal_mass_eigenvalue']<0
    assert np.linalg.eigvalsh(out['K']).min()<0


def test_invalid_inputs_and_invalid_slices_are_rejected():
    z=simple()
    for name,value in [('Y',np.ones((2,3))),('DY',np.zeros((2,2,2))),
            ('weights',[-1.,1.]),('weights',[0.,0.]),('fields',[0,2]),
            ('fields',[False,True]),('background',np.zeros((2,3))),
            ('sstar',True),('sstar',float('nan'))]:
        with pytest.raises(ValueError):quadratic_hamiltonian(**dict(z,**{name:value}))
    bad=z['g'].copy();bad[0,1,1]=-1.
    with pytest.raises(ValueError,match='positive definite'):quadratic_hamiltonian(**dict(z,g=bad))
    bad=z['g'].copy();bad[0,0,0]=0.
    with pytest.raises(ValueError,match='lapse'):quadratic_hamiltonian(**dict(z,g=bad))
    bad=z['g'].copy();bad[0,0,1]=.1
    with pytest.raises(ValueError,match='symmetric'):quadratic_hamiltonian(**dict(z,g=bad))
    with pytest.raises(ValueError):potential_hessian([0.,float('nan')],sstar=SSTAR)
    with pytest.raises(ValueError):potential_hessian([0.,0.],sstar=-1000)
