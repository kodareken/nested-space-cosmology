"""Extracted-map controls against independent continuum cylinder formulas."""
import math

import numpy as np
from scipy.special import j0,j1,k0,k1
from numpy.polynomial.legendre import leggauss

from recursive_horizons.nsc_atlas_dirac_retarded import (
    retarded_null_responses,retarded_response,volterra_endpoint_jets)
from recursive_horizons.nsc_atlas_dirac_root_pairing import pair_cylinder_responses


def constant_maps(U,V,mu,n,order=24):
    shape=(n,n);du=U/(n-.5);dv=V/(n-.5)
    return retarded_null_responses(np.full(shape,mu),du,dv,
        potential_u=np.zeros(shape),potential_v=np.zeros(shape),
        boundary_square_integrals=(mu*mu*V,mu*mu*U),boundary_mixed_integrals=(0.,0.),
        order=order,workers=8)


def constant_input_exact(U,V,mu):
    z=mu*math.sqrt(U*V);dz=mu*(U-V)/(2*math.sqrt(U*V))
    dp=j0(z)-j1(z)/z;factor=(1/U+1/V)/2
    value=np.array([[j0(z),-1j*math.sqrt(V/U)*j1(z)],
        [-1j*math.sqrt(U/V)*j1(z),j0(z)]])
    dx=np.array([[-j1(z)*dz,-1j*math.sqrt(V/U)*(factor*j1(z)+dp*dz)],
        [-1j*math.sqrt(U/V)*(-factor*j1(z)+dp*dz),-j1(z)*dz]])
    return value,dx


def test_extracted_field_and_analytic_spatial_maps_converge_before_any_physical_reference():
    eta=.008;T=.11;x=.006;mu=1/.45
    s1=np.array([[0.,1.],[1.,0.]]);s3=np.diag([1.,-1.])
    exact=-mu*s1*k0(mu*eta)/(2*math.pi)-1j*mu*s3*k1(mu*eta)/(2*math.pi)
    exact_dx=mu*mu/(2*math.pi)*(s1*k1(mu*eta)+1j*s3*(k0(mu*eta)+k1(mu*eta)/(mu*eta)))
    errors=[]
    for n in (32,64):
        L=(T-x-eta/2,T+x+eta/2);R=(T-x+eta/2,T+x-eta/2)
        left=constant_maps(*L,mu,n);right=constant_maps(*R,mu,n)
        row_error=[]
        for sample,point in ((left,L),(right,R)):
            target=constant_input_exact(*point,mu)
            for key,expected in zip(('field','spatial'),target):
                response=retarded_response(sample,key,24)
                reconstructed=response.direct+np.stack([
                    response.lengths[a]*response.green[a][0] for a in (0,1)],axis=-1)
                row_error.append(float(np.max(abs(reconstructed-expected))))
        kernel_errors=[]
        for order in (12,24):
            lf,rf=(retarded_response(s,'field',order) for s in (left,right))
            ld,rd=(retarded_response(s,'spatial',order) for s in (left,right))
            kernel=pair_cylinder_responses(lf,rf,mu,epsabs=2e-11)['kernel']
            dx=(pair_cylinder_responses(ld,rf,mu,epsabs=2e-11)['kernel']
                -pair_cylinder_responses(lf,rd,mu,epsabs=2e-11)['kernel'])/2
            kernel_errors.extend((np.max(abs(kernel-exact)),np.max(abs(dx-exact_dx))))
        errors.append(np.r_[row_error,kernel_errors])
        print('ENDPOINT_EXTRACTION_CONTROL',n,errors[-1].tolist(),
            'PROJECTION_L2',left['cell_projection_residual'].tolist())
    assert np.max(errors[1][:4])<1e-3
    assert np.max(errors[1][4::2])<1e-3
    assert np.max(errors[1][5::2])<.1
    assert np.all(errors[0]/errors[1]>3.)


def test_variable_potential_and_mixed_endpoint_integrals_against_exact_separable_solution():
    # ubar=int f^2 du, vbar=int g^2 dv reduces W=f(u)g(v) to unit W.
    U,V=.104,.116;p=.7;q=1/(4*p)
    f=lambda s:1+.5*s
    g=lambda s:1-.3*s
    A=lambda s:s+.5*s*s+s**3/12
    B=lambda s:s-.3*s*s+.03*s**3
    def incoming(leg,s):
        return f(s)*np.exp(-1j*p*A(s)) if leg==0 else 2*q*g(s)*np.exp(-1j*q*B(s))
    def incoming_derivative(leg,s):
        return ((.5-1j*p*f(s)**3)*np.exp(-1j*p*A(s)) if leg==0 else
                2*q*(-.3-1j*q*g(s)**3)*np.exp(-1j*q*B(s)))
    phase=np.exp(-1j*(p*A(U)+q*B(V)))
    expected=phase*np.array([f(U),2*q*g(V)])
    expected_x=phase*np.array([-.5,-.6*q])+1j*(p*f(U)**2-q*g(V)**2)*expected
    nodes,weights=leggauss(48);errors=[]
    for n in (32,64):
        du,dv=U/(n-.5),V/(n-.5)
        u,v=np.meshgrid((np.arange(n)+.5)*du,(np.arange(n)+.5)*dv,indexing='ij')
        sample=retarded_null_responses(f(u)*g(v),du,dv,potential_u=.5*g(v),potential_v=-.3*f(u),
            boundary_square_integrals=(f(U)**2*B(V),g(V)**2*A(U)),
            boundary_mixed_integrals=(.5*f(U)*B(V),-.3*g(V)*A(U)),order=24,workers=8)
        values=[]
        for kind in ('field','spatial'):
            response=retarded_response(sample,kind,24)
            value=sum(response.direct[:,a]*incoming(a,L)+response.direct_derivative[:,a]*incoming_derivative(a,L)
                +np.einsum('s,si,s->i',weights*L/2,response.value(a,(nodes+1)*L/2),incoming(a,(nodes+1)*L/2))
                for a,L in enumerate(response.lengths))
            values.append(value)
        errors.append(np.array([np.max(abs(values[0]-expected)),np.max(abs(values[1]-expected_x))]))
        print('VARIABLE_ENDPOINT_CONTROL',n,errors[-1].tolist())
    assert np.max(errors[1])<2e-5
    assert np.all(errors[0]/errors[1]>3.)


def test_order_48_uses_the_same_existing_grid_and_exact_cylinder_kernel():
    # The physical polynomial check changes representation, not field resolution.
    T,x,eta,mu=.11,.006,.008,1/.45
    left=constant_maps(T-x-eta/2,T+x+eta/2,mu,64,order=48)
    right=constant_maps(T-x+eta/2,T+x-eta/2,mu,64,order=48)
    lf,rf=(retarded_response(s,'field',48) for s in (left,right))
    ld,rd=(retarded_response(s,'spatial',48) for s in (left,right))
    kernel=pair_cylinder_responses(lf,rf,mu,epsabs=2e-11)['kernel']
    derivative=(pair_cylinder_responses(ld,rf,mu,epsabs=2e-11)['kernel']
        -pair_cylinder_responses(lf,rd,mu,epsabs=2e-11)['kernel'])/2
    s1=np.array([[0.,1.],[1.,0.]]);s3=np.diag([1.,-1.])
    exact=-mu*s1*k0(mu*eta)/(2*math.pi)-1j*mu*s3*k1(mu*eta)/(2*math.pi)
    exact_dx=mu*mu/(2*math.pi)*(s1*k1(mu*eta)+1j*s3*(k0(mu*eta)+k1(mu*eta)/(mu*eta)))
    errors=np.array([np.max(abs(kernel-exact)),np.max(abs(derivative-exact_dx))])
    print('FIXED_GRID_ORDER48_CONTROL',errors.tolist())
    assert errors[0]<1e-6
    assert errors[1]<3e-6


def test_source_endpoint_jets_against_the_independent_separable_Riemann_kernel():
    U,V=.104,.116
    f=lambda s:1+.5*s;g=lambda s:1-.3*s
    A=lambda s:s+.5*s*s+s**3/12;B=lambda s:s-.3*s*s+.03*s**3
    jets=volterra_endpoint_jets(f(U)*g(V),.5*g(V),-.3*f(U),
        (f(U)**2*B(V),g(V)**2*A(U)),(.5*f(U)*B(V),-.3*g(V)*A(U)))
    def kernel(leg,s):
        z=math.sqrt(max(0.,(A(U)-A(s))*B(V) if leg==0 else A(U)*(B(V)-B(s))))
        ratio=.5 if z==0 else j1(z)/z
        if leg==0:
            return np.array([-f(U)*f(s)*B(V)*ratio/2,-.5j*g(V)*f(s)*j0(z)])
        return np.array([-.5j*f(U)*g(s)*j0(z),-g(V)*g(s)*A(U)*ratio/2])
    for leg,L in enumerate((U,V)):
        h=2e-4;values=[kernel(leg,L-j*h) for j in range(5)]
        derivative=(25*values[0]-48*values[1]+36*values[2]-16*values[3]+3*values[4])/(12*h)
        assert np.max(abs(derivative-jets['source'][leg]))<3e-8
