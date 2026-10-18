"""Extracted-map controls against independent continuum cylinder formulas."""
import math

import numpy as np
from scipy.special import j0,j1,k0,k1
from numpy.polynomial.legendre import leggauss

from recursive_horizons.nsc_atlas_dirac_retarded import (
    retarded_null_responses,retarded_response,volterra_endpoint_jets,
    observation_moment_residual)
from recursive_horizons.nsc_atlas_dirac_retarded import (
    constrained_projection_derivative,endpoint_hilbert_moment,
    volterra_source_observation_jets)
from recursive_horizons.nsc_atlas_dirac_root_pairing import (
    NullResponse,pair_cylinder_responses,logarithmic_coefficients)


def constant_maps(U,V,mu,n,order=24):
    shape=(n,n);du=U/(n-.5);dv=V/(n-.5)
    return retarded_null_responses(np.full(shape,mu),du,dv,
        potential_u=np.zeros(shape),potential_v=np.zeros(shape),
        boundary_square_integrals=(mu*mu*V,mu*mu*U),boundary_mixed_integrals=(0.,0.),
        order=order,workers=8)


def test_projection_derivative_recovers_cubic_finite_hilbert_moment_and_mean():
    L=.8;c=np.array([L**3/4,9*L**3/20,L**3/4])[:,None]*np.array([[1.,2j]])
    zero=np.zeros((3,2),complex);v=np.array([1.,2j])
    for sign in (-1.,1.):
        d=constrained_projection_derivative(c,zero,L,L**3*v,np.zeros(2),
            3*L**2*v,sign*6*L*v,orientation=sign)
        expected=sign*np.array([0.,-3*L**2/8,3*L**2/8])[:,None]*v
        np.testing.assert_allclose(d,expected,atol=2e-15)
        np.testing.assert_allclose(d[0],0.,atol=2e-15)
        np.testing.assert_allclose(d.sum(axis=0),0.,atol=2e-15)
        np.testing.assert_allclose(endpoint_hilbert_moment(d),sign*3*L*L*v/8,atol=3e-15)


def test_projection_derivative_matches_fixed_source_differentiation_of_smooth_kernel():
    from numpy.polynomial.legendre import legvander
    p=6;L=.8;x=.15;n=np.arange(p);nodes,weights=leggauss(48);s=np.array([.07,.29,.63])
    direction=np.array([1.,.3j])
    def sample(length,position):
        values=np.exp(position+(nodes+1)*length/2)
        c=((2*n+1)/2*(legvander(nodes,p-1).T@(weights*values)))[:,None]*direction
        gamma=np.exp(position+length)*direction
        return dict(cell_field=NullResponse((length,length),(c,c),np.eye(2)),
            endpoint_values=np.stack((np.stack((gamma,gamma)),np.stack((gamma,gamma)))),
            source_endpoint_values=np.stack((gamma,gamma))),c,gamma
    base,c,gamma=sample(L,x)
    for sign in (-1.,1.):
        d=constrained_projection_derivative(c,c,L,gamma,gamma,gamma,(1+sign)*gamma,orientation=sign)
        values=legvander(2*s/L-1,p-1)@d;h=1e-4
        def shifted(delta):
            data,_,_=sample(L+sign*delta,x+delta)
            return retarded_response(data,'field',p,source_jet=True).value(0,s)
        independent=(-shifted(2*h)+8*shifted(h)-8*shifted(-h)+shifted(-2*h))/(12*h)
        np.testing.assert_allclose(values,independent,rtol=2e-8,atol=2e-8)


def test_volterra_source_observation_targets_match_separable_potential():
    U,V=.17,.21;a,b=.7,-.4
    def jets(u,v):
        W=np.exp(a*u+b*v)
        I=np.array([np.exp(2*a*u)*np.expm1(2*b*v)/(2*b),
                    np.exp(2*b*v)*np.expm1(2*a*u)/(2*a)])
        J=np.array([a,b])*I
        return W,I,J
    W,I,J=jets(U,V);K=2*np.array([a*a,b*b])*I
    actual=volterra_source_observation_jets(W,a*W,b*W,a*a*W,a*b*W,b*b*W,I,J,K)
    def source(delta):
        w,i,j=jets(U-delta,V+delta)
        return volterra_endpoint_jets(w,a*w,b*w,i,j)['source']
    h=1e-4
    independent=(-source(2*h)+8*source(h)-8*source(-h)+source(-2*h))/(12*h)
    np.testing.assert_allclose(actual,independent,rtol=2e-10,atol=2e-11)


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


def test_source_constraints_preserve_mean_and_are_the_unique_minimum_L2_correction():
    order=8;rng=np.random.default_rng(20843)
    coefficients=tuple(rng.normal(size=(order,2))+1j*rng.normal(size=(order,2)) for _ in (0,1))
    lengths=(.11,.13);jets=volterra_endpoint_jets(1.7,.2,-.3,(.32,.41),(.023,-.019))
    sample=dict(cell_field=NullResponse(lengths,coefficients,np.eye(2)),
        endpoint_values=np.stack((jets['value'],jets['spatial'])),source_endpoint_values=jets['source'])
    corrected=retarded_response(sample,'field',order,source_jet=True)
    for a,L in enumerate(lengths):
        assert np.array_equal(corrected.green[a][0],coefficients[a][0])
        assert np.max(abs(corrected.value(a,L)-jets['value'][a]))<2e-13
        assert np.max(abs(corrected.value(a,L,derivative=1)-jets['source'][a]))<2e-11
        n=np.arange(order)
        constraints=np.stack((np.eye(order)[0],np.ones(order),n*(n+1)/L))
        constraints/=np.linalg.norm(constraints,axis=1)[:,None]
        _,_,v=np.linalg.svd(constraints,full_matrices=True)
        null=v[3:].T
        delta=corrected.green[a]-coefficients[a]
        # Independent optimality condition: weighted correction is orthogonal
        # to EVERY variation preserving the three constraints.
        orthogonal=null.T@(delta*(L/(2*n+1))[:,None])
        assert np.max(abs(orthogonal))<2e-13


def test_source_compatible_reconstruction_retains_the_exact_cylinder_kernel_and_log():
    T,x,eta,mu=.11,.006,.008,1/.45
    left=constant_maps(T-x-eta/2,T+x+eta/2,mu,64,order=48)
    right=constant_maps(T-x+eta/2,T+x-eta/2,mu,64,order=48)
    s1=np.array([[0.,1.],[1.,0.]]);s3=np.diag([1.,-1.])
    exact=-mu*s1*k0(mu*eta)/(2*math.pi)-1j*mu*s3*k1(mu*eta)/(2*math.pi)
    exact_dx=mu*mu/(2*math.pi)*(s1*k1(mu*eta)+1j*s3*(k0(mu*eta)+k1(mu*eta)/(mu*eta)))
    for order in (24,48):
        lf,rf=(retarded_response(s,'field',order,source_jet=True) for s in (left,right))
        ld,rd=(retarded_response(s,'spatial',order) for s in (left,right))
        kernel=pair_cylinder_responses(lf,rf,mu,epsabs=2e-11)['kernel']
        dx=(pair_cylinder_responses(ld,rf,mu,epsabs=2e-11)['kernel']
            -pair_cylinder_responses(lf,rd,mu,epsabs=2e-11)['kernel'])/2
        print('SOURCE_COMPATIBLE_CYLINDER',order,float(np.max(abs(kernel-exact))),float(np.max(abs(dx-exact_dx))))
        assert np.max(abs(kernel-exact))<1e-6
        assert np.max(abs(dx-exact_dx))<3e-6
        for f,d in ((lf,ld),(rf,rd)):
            coefficient=logarithmic_coefficients(f,d)
            assert abs(coefficient['q_alpha']-mu*mu/(2*math.pi))<1e-10
            assert abs(coefficient['angular']+mu/math.pi)<3e-13


def test_observation_moments_keep_moving_boundary_signs_and_exact_quadrature():
    S=1.2;eta=.02;lengths=np.array([.8,.6]);sign=np.array([-1.,1.])
    vectors=(np.array([.3,.7j]),np.array([-.4j,.2]))
    greens=tuple(v[None,:] for v in vectors);zero=(np.zeros((1,2)),np.zeros((1,2)))
    left=NullResponse(tuple(lengths+sign*eta/2),greens,np.eye(2))
    right=NullResponse(tuple(lengths-sign*eta/2),greens,np.eye(2))
    dl=NullResponse(left.lengths,zero,np.zeros((2,2)))
    dr=NullResponse(right.lengths,zero,np.zeros((2,2)))
    residual=observation_moment_residual(left,dl,right,dr,separation=eta,reference_duration=S,degrees=(0,1,2))
    assert np.max(abs(residual[:,:2]))<3e-14
    for a in (0,1):
        assert np.max(abs(residual[a,2]+sign[a]*eta**2/S**2*vectors[a]))<3e-14


def test_nonconstant_probe_detects_projection_commutator_hidden_from_mean():
    def projected(L):
        coefficients=np.zeros((3,2),complex);coefficients[:,1]=L**3*np.array([.25,.45,.25])
        sample=dict(cell_field=NullResponse((.5,L),(np.zeros((3,2)),coefficients),np.eye(2)),
            endpoint_values=np.zeros((2,2,2),complex),source_endpoint_values=np.zeros((2,2),complex))
        sample['endpoint_values'][0,1,1]=L**3
        sample['source_endpoint_values'][1,1]=3*L*L
        return retarded_response(sample,'field',3,source_jet=True)
    L=.8;S=1.;zero=(np.zeros((1,2)),np.zeros((1,2)))
    for eta in (.02,.01):
        # G=s^3 has zero true observer derivative. Its constrained projection
        # changes with L, even though its mean derivative error cancels.
        left=projected(L+eta/2);right=projected(L-eta/2)
        left=NullResponse((.5-eta/2,left.lengths[1]),left.green,left.direct)
        right=NullResponse((.5+eta/2,right.lengths[1]),right.green,right.direct)
        dl=NullResponse(left.lengths,zero,np.zeros((2,2)))
        dr=NullResponse(right.lengths,zero,np.zeros((2,2)))
        residual=observation_moment_residual(left,dl,right,dr,separation=eta,reference_duration=S,degrees=(0,1))
        mean=-L*eta*eta/2
        first_monomial=-L**4/16-33*L*L*eta*eta/32-13*eta**4/256
        assert abs(residual[1,0,1]-mean)<2e-13
        assert abs(residual[1,1,1]-(2*first_monomial/S-mean))<2e-13
