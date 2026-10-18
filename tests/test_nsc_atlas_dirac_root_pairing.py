"""Independent full-radial cylinder and endpoint controls; no physical grid."""
import math

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss,legvander
from scipy.special import j0,j1,jv,k0,k1

from recursive_horizons.nsc_atlas_dirac_root_pairing import (
    NullResponse,pair_cylinder_responses,logarithmic_coefficients)


def cylinder_response(U,V,mu,order,*,derivative=False):
    nodes,weights=leggauss(order)
    coefficients=[]
    for leg,L in enumerate((U,V)):
        other=(V,U)[leg]
        z=abs(mu)*np.sqrt((L-(nodes+1)*L/2)*other)
        diagonal=-mu*mu*other/2*j1(z)/z
        offdiagonal=-.5j*mu*j0(z)
        if derivative:
            f=j1(z)/z;df=-jv(2,z)/(2*z*z)
            dw=mu*mu*((L-(nodes+1)*L/2)-other)*(1 if leg==0 else -1)
            diagonal=-mu*mu/2*((1 if leg==0 else -1)*f+other*df*dw)
            offdiagonal=.25j*mu*f*dw
        value=np.stack((diagonal,offdiagonal) if leg==0 else (offdiagonal,diagonal),axis=-1)
        coefficients.append(np.einsum('s,sn,si->ni',weights,legvander(nodes,order-1),value)
                            *(2*np.arange(order)+1)[:,None]/2)
    if derivative:
        D=np.array([[mu*mu*V/4,-.5j*mu],[.5j*mu,-mu*mu*U/4]])
        return NullResponse((U,V),tuple(coefficients),D,np.diag([-1.,1.]))
    return NullResponse((U,V),tuple(coefficients),np.eye(2))


def test_one_sided_Hilbert_transform_keeps_the_endpoint_logs():
    response=NullResponse((.7,.8),(np.array([[1.,0.],[0.,1.]]),np.zeros((1,2))),np.eye(2))
    for s in (-.2,.2,.5,1.):
        z=2*s/.7-1;log=math.log(abs(s/(s-.7)))
        assert response.hilbert(0,s)==pytest.approx([log,z*log-2],abs=2e-14)
    with pytest.raises(ValueError,match='endpoint'):
        response.hilbert(0,.7)


def test_complete_radial_cylinder_pairing_matches_the_independent_Fourier_Bessel_kernel():
    eta=.008;T=.11;x=.006
    s1=np.array([[0.,1.],[1.,0.]]);s3=np.diag([1.,-1.])
    for mu in (-1/.45,1/.45):
        exact=-mu*s1*k0(abs(mu*eta))/(2*math.pi)-1j*abs(mu)*s3*k1(abs(mu*eta))/(2*math.pi)
        for order in (12,24):
            left=cylinder_response(T-x-eta/2,T+x+eta/2,mu,order)
            right=cylinder_response(T-x+eta/2,T+x-eta/2,mu,order)
            sample=pair_cylinder_responses(left,right,mu,epsabs=2e-11)
            assert sample['kernel']==pytest.approx(exact,abs=2e-9)
            assert sample['pv_quadrature_estimate']<2e-10
            assert sample['bessel_series_tail_bound']<1e-25
            assert sample['scope']=='complete_root_covariance_on_supplied_polynomial_responses'


def test_cross_leg_series_bound_and_noncoincident_domain():
    zero=(np.zeros((1,2)),np.zeros((1,2)))
    left=NullResponse((.7,.8),zero,np.eye(2));right=NullResponse((.6,.9),zero,np.eye(2))
    for terms in (12,24):
        sample=pair_cylinder_responses(left,right,4.,series_terms=terms)
        expected=np.array([[1j/(2*math.pi*(.7-.6)),-4*k0(4*math.sqrt(.7*.9))/(2*math.pi)],
            [-4*k0(4*math.sqrt(.8*.6))/(2*math.pi),1j/(2*math.pi*(.8-.9))]])
        assert np.linalg.norm(sample['kernel']-expected,ord=2)<=sample['bessel_series_tail_bound']+3e-14
    with pytest.raises(ValueError,match='noncoincident'):
        pair_cylinder_responses(left,left,4.)
    with pytest.raises(ValueError,match='unstable'):
        pair_cylinder_responses(left,right,10.)


def test_analytic_spatial_response_includes_moving_endpoint_and_delta_derivative():
    eta=.008;T=.11;x=.006;mu=1/.45
    s1=np.array([[0.,1.],[1.,0.]]);s3=np.diag([1.,-1.])
    exact=mu*mu/(2*math.pi)*(s1*k1(mu*eta)+1j*s3*(k0(mu*eta)+k1(mu*eta)/(mu*eta)))
    for order in (12,24):
        L=(T-x-eta/2,T+x+eta/2);R=(T-x+eta/2,T+x-eta/2)
        left,right=cylinder_response(*L,mu,order),cylinder_response(*R,mu,order)
        dl,dr=cylinder_response(*L,mu,order,derivative=True),cylinder_response(*R,mu,order,derivative=True)
        first=pair_cylinder_responses(dl,right,mu,epsabs=2e-11)
        second=pair_cylinder_responses(left,dr,mu,epsabs=2e-11)
        assert first['kernel']==pytest.approx(exact,abs=3e-8)
        assert second['kernel']==pytest.approx(-exact,abs=3e-8)
        assert first['pv_quadrature_estimate']<3e-10


def test_log_coefficients_match_the_independent_cylinder_Fourier_asymptotics():
    for mu in (-1/.45,1/.45):
        f=cylinder_response(.104,.116,mu,24)
        d=cylinder_response(.104,.116,mu,24,derivative=True)
        coefficient=logarithmic_coefficients(f,d)
        assert coefficient['q_alpha']==pytest.approx(mu*mu/(2*math.pi),abs=3e-10)
        assert coefficient['angular']==pytest.approx(-mu/math.pi,abs=3e-13)


def test_source_slope_can_change_the_stress_log_without_changing_value_or_mean():
    mu=1/.45;f=cylinder_response(.104,.116,mu,24)
    d=cylinder_response(.104,.116,mu,24,derivative=True)
    baseline=logarithmic_coefficients(f,d)
    altered=[g.copy() for g in f.green];epsilon=.01
    altered[0][1,0]-=epsilon;altered[0][2,0]+=epsilon
    changed=NullResponse(f.lengths,tuple(altered),f.direct)
    assert changed.value(0,.104).shape==(2,)
    assert changed.value(0,.104,derivative=1).shape==(2,)
    coefficient=logarithmic_coefficients(changed,d)
    assert changed.value(0,.104)==pytest.approx(f.value(0,.104),abs=1e-14)
    assert changed.green[0][0]==pytest.approx(f.green[0][0],abs=1e-14)
    assert coefficient['angular']==pytest.approx(baseline['angular'],abs=1e-14)
    assert coefficient['q_alpha']-baseline['q_alpha']==pytest.approx(-4*epsilon/(2*math.pi*.104),abs=1e-12)
