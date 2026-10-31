"""Intrinsic spatial Dirac modes on the defined one-mode initial parent.

The physical isometry is left quaternion multiplication
Z=(-X4,X3,-X2,X1), so Z4=oldX1 and the right frame is preserved.
Nothing is angularly averaged to construct a different geometry.

For h=Omega² gamma, D_h(Omega^-1 f)=Omega^-2 D_gamma f. The
generalized eigenproblem is D_gamma f=lambda Omega f, with
integral Omega|f|²=1. These operators and normalized projector overlaps
do not specify an occupation, quantum state, contact source or backreaction.
"""
import math

import numpy as np
from numpy.polynomial.chebyshev import chebinterpolate,chebint,chebval
from numpy.polynomial.legendre import leggauss
from scipy.linalg import block_diag,eigh
from scipy.special import eval_gegenbauer,roots_jacobi

from .nsc_atlas_s3_spectral import _radial_normalization


def _degree(value):
    if isinstance(value,(bool,np.bool_)) or not isinstance(value,(int,np.integer)) or not 1<=value<=32:
        raise ValueError('integer degree in [1,32] required')
    return int(value)


def radial_values(l,j,u):
    """Normalized S3 radial factor; S2 spinor harmonics have unit norm."""
    return _radial_normalization(l,j)*(1-u*u)**(j/2)*eval_gegenbauer(l-j,j+1,u)


def omega_bounds(initial):
    """Analytic coefficient bounds for the declared finite psi polynomial.

    |C_l^1(u)|<=l+1 supplies both bounds. The floating coefficient and
    rounding padding are explicit; this does not enclose an unknown exact
    Einstein-data solution beyond the defined finite polynomial geometry.
    """
    labels=np.arange(0,initial.degree+1,2);c=np.asarray(initial.coefficients)
    normalization=math.sqrt(2*math.pi**2)
    center=float(c[0])/normalization
    tail=math.fsum((int(l)+1)*abs(float(v)) for l,v in zip(labels[1:],c[1:]))/normalization
    padding=128*np.finfo(float).eps*(abs(center)+tail)
    lower=center-tail-padding;upper=center+tail+padding
    if lower<=0:raise ValueError('finite polynomial coefficient bound does not establish positive psi')
    return dict(psi_lower=lower,psi_upper=upper,
                Omega_lower=math.nextafter(lower**2/initial.parent.H,0.),
                Omega_upper=math.nextafter(upper**2/initial.parent.H,math.inf),
                coefficient_roundoff_padding=padding,
                derivation='|C_l^1(u)|<=l+1; normalized constant coefficient plus/minus the absolute weighted nonconstant coefficients.',
                scope='Defined finite coefficient polynomial with floating-rounding caveat, not a certified enclosure of the exact Einstein branch.')


def radial_operator_blocks(degree,k,u,w,omega):
    """Exact round Dirac blocks and 1D Omega/Omega² multiplication.

    Total SO3 spin J=k+1/2 has multiplicity2k+2. Its two angular branches
    have orbital j=k,k+1. A constant i rephasing of the second branch makes
    the round blocks real symmetric. At each l their eigenvalues are
    l+3/2 and -l-1/2; at l=k only the positive first branch exists.
    """
    degree=_degree(degree)
    if isinstance(k,(bool,np.bool_)) or not isinstance(k,(int,np.integer)) or not 0<=k<=degree:
        raise ValueError('integer k between0 anddegree required')
    u=np.asarray(u,dtype=float);w=np.asarray(w,dtype=float);omega=np.asarray(omega,dtype=float)
    if (u.ndim!=1 or w.shape!=u.shape or omega.shape!=u.shape or
            any(np.any(~np.isfinite(v)) for v in (u,w,omega)) or
            np.any(abs(u)>=1) or np.any(w<=0) or np.any(omega<=0)):
        raise ValueError('finite interior radial nodes, positive weights and Omega required')
    first=list(range(k,degree+1));second=list(range(k+1,degree+1));n0=len(first)
    R0=np.column_stack([radial_values(l,k,u) for l in first])
    R1=np.column_stack([radial_values(l,k+1,u) for l in second]) if second else np.empty((len(u),0))
    M=block_diag(R0.T@((w*omega)[:,None]*R0),R1.T@((w*omega)[:,None]*R1))
    M2=block_diag(R0.T@((w*omega**2)[:,None]*R0),R1.T@((w*omega**2)[:,None]*R1))
    A=np.diag([k+1.5]*n0+[-k-.5]*len(second))
    for offset,l in enumerate(second):
        i=l-k;j=n0+offset;A[i,j]=A[j,i]=math.sqrt((l-k)*(l+k+2))
    return A,M,M2,R0,R1


def mode_radial_profiles(mode,u):
    """Actual unoccupied mode profiles, ordered j=k then j=k+1."""
    u=np.atleast_1d(np.asarray(u,dtype=float))
    if np.any(~np.isfinite(u)) or np.any(abs(u)>1):raise ValueError('finite u in [-1,1] required')
    L=int(mode['degree']);k=int(mode['k']);c=np.asarray(mode['coefficients'])
    n0=L-k+1
    if c.shape!=(2*(L-k)+1,):raise ValueError('mode coefficient shape does not match the angular sector')
    F=sum(v*radial_values(l,k,u) for l,v in zip(range(k,L+1),c[:n0]))
    G=sum((v*radial_values(l,k+1,u) for l,v in zip(range(k+1,L+1),c[n0:])),np.zeros_like(u))
    return F,G


def full_mode_diagnostics(A,M,M2,R0,R1,u,w,omega,eigenvalue,coefficients):
    """Positive full PDE residual, including output outside the Ritz space.

    M2=P Omega² P is not (P Omega P)². The latter subtraction can lose
    its sign at roundoff; the primary residual integrates the actual full
    radial defect squared. Its physical norm weights that defect byOmega^-1.
    """
    c=np.asarray(coefficients);lam=float(eigenvalue);n0=R0.shape[1];Ac=A@c
    F=R0@c[:n0];G=R1@c[n0:]
    rF=R0@Ac[:n0]-lam*omega*F;rG=R1@Ac[n0:]-lam*omega*G
    gamma2=float(w@(rF*rF+rG*rG));physical2=float(w@((rF*rF+rG*rG)/omega))
    density=F*F+G*G;I4=float(w@(density*density/omega))/(4*math.pi)
    identity=lam*lam*float(c@(M2-M@M)@c)
    return dict(full_PDE_residual_gamma=math.sqrt(gamma2),full_PDE_residual_physical=math.sqrt(physical2),
                residual_gamma_squared_from_Momega2=identity,positive_residual_gamma_squared=gamma2,
                Momega2_identity_absolute_gap=abs(identity-gamma2),
                relative_full_gamma_residual=math.sqrt(gamma2)/max(1e-300,abs(lam)*math.sqrt(float(c@M2@c))),
                projected_Ritz_residual=float(np.linalg.norm(A@c-lam*M@c)),
                physical_normalization=float(w@(omega*density)),
                invariant_ground_projector_I4=I4,effective_volume=1/I4,
                density_formula='rho_h=Omega^-2*(F²+G²)/(4pi), normalized complete-J projector mean')


def compact_inverse_validation(summary,bounds):
    """Isolate the global first cluster of the defined finite-psi geometry.

    K=Omega^(1/2) D_gamma^-1 Omega^(1/2). Its finite-rank approximation
    has the reciprocal generalized Ritz eigenvalues. Since D_gamma
    preserves scalar degree, ||K-K_L||<=Omega_max/(L+1.5).
    """
    first=sorted(v['eigenvalue'] for v in summary['first_spatial_modes'] if v['eigenvalue']>0)
    mu0=1/first[0];mu1=1/first[1];delta=bounds['Omega_upper']/(summary['degree']+1.5)
    gap=mu0-mu1;positive=next(v for v in summary['lowest_by_sign'] if v['sign']>0)
    lam=positive['eigenvalue'];residual=positive['full_PDE_residual_physical']
    interval=[1/(lam+residual),1/(lam-residual)] if lam>residual else None
    inside=(interval is not None and interval[0]>mu1+delta and interval[0]>=mu0-delta and interval[1]<=mu0+delta)
    return dict(inverse_operator_norm_error_upper=delta,largest_positive_inverse=mu0,next_distinct_positive_inverse=mu1,
                finite_inverse_gap=gap,twice_norm_error=2*delta,isolation_margin=gap-2*delta,
                isolated_ground_cluster=gap>2*delta,finite_ground_cluster_dimension=positive['spatial_Dirac_multiplicity'],
                residual_eigenvalue_interval=[lam-residual,lam+residual],inverse_of_residual_interval=interval,
                residual_interval_inside_isolated_global_first_cluster=bool(gap>2*delta and inside),
                proof='First omitted absolute round eigenvalueL+1.5; inverse norm error<=Omega_max/(L+1.5). Gap>2delta isolates the dimension2 first cluster; the positive full physical residual identifies its eigenvalue.',
                caveat='Analytic operator estimate with floating eigensolve/quadrature/finite-coefficient caveats, not an interval certificate or a claim about the unknown exact Einstein branch.')


def calculate_dirac_level(degree,*,initial,mass=None):
    """Return compact spectral diagnostics and raw in-memory normalized modes.

    ``initial`` is the existing NeckInitialData. ``mass`` is a separately
    declared candidate input, not selected by this spatial calculation.
    The caller owns source provenance and raw-coefficient storage.
    """
    degree=_degree(degree)
    if mass is not None and (not math.isfinite(mass) or mass<0):raise ValueError('finite nonnegative declared mass required')
    order=degree+2*initial.degree+21
    u,w=roots_jacobi(order,.5,.5);omega=initial.values(u)**2/initial.parent.H
    spectra=[];candidates=[];round_error=0.;dimension=0
    for k in range(degree+1):
        A,M,M2,R0,R1=radial_operator_blocks(degree,k,u,w,omega)
        values,vectors=eigh(A,M,check_finite=True)
        expected=sorted([l+1.5 for l in range(k,degree+1)]+[-l-.5 for l in range(k+1,degree+1)])
        round_error=max(round_error,float(np.max(abs(np.linalg.eigvalsh(A)-expected))))
        degeneracy=2*k+2;dimension+=degeneracy*len(values)
        spectra.extend(dict(eigenvalue=float(v),k=k,multiplicity=degeneracy) for v in values)
        for sign in (-1,1):
            choices=np.flatnonzero(values*sign>0)
            if not len(choices):continue
            index=int(choices[np.argmin(abs(values[choices]))]);lam=float(values[index])
            candidates.append((abs(lam),sign,k,lam,vectors[:,index],(A,M,M2,R0,R1)))
    records=[];modes=[]
    for sign in (-1,1):
        _,_,k,lam,c,pack=min(v for v in candidates if v[1]==sign)
        details=full_mode_diagnostics(*pack,u,w,omega,lam,c)
        uv,wv=roots_jacobi(2*order,.5,.5);ov=initial.values(uv)**2/initial.parent.H
        validation=full_mode_diagnostics(*radial_operator_blocks(degree,k,uv,wv,ov),uv,wv,ov,lam,c)
        details['I4_quadrature_refinement_gap']=abs(details['invariant_ground_projector_I4']-validation['invariant_ground_projector_I4'])
        details['physical_residual_quadrature_refinement_gap']=abs(details['full_PDE_residual_physical']-validation['full_PDE_residual_physical'])
        record=dict(sign=sign,k=k,spatial_Dirac_multiplicity=2*k+2,eigenvalue=lam,
                    lambda_over_H=lam/initial.parent.H,**details)
        if mass is not None:
            record.update(declared_mass=mass,energy=math.sqrt(lam*lam+mass*mass),lambda_over_energy=abs(lam)/math.sqrt(lam*lam+mass*mass))
        records.append(record);modes.append(dict(degree=degree,k=k,sign=sign,eigenvalue=lam,coefficients=c.copy()))
    summary=dict(degree=degree,spinor_dimension=dimension,expected_dimension=2*sum((l+1)**2 for l in range(degree+1)),
                 quadrature_order=order,validation_order=2*order,round_spectrum_max_error=round_error,
                 all_k_through_scalar_cutoff_retained=True,lowest_by_sign=records,
                 first_spatial_modes=sorted(spectra,key=lambda v:abs(v['eigenvalue']))[:12])
    summary['compact_inverse_validation']=compact_inverse_validation(summary,omega_bounds(initial))
    return summary,tuple(modes)


def _initial_parent_geometry_accelerations(initial,u):
    """Actual time-symmetric initial ADM jets, before any spinor source."""
    p=initial.parent;H=p.H;kappa=p.kappa
    psi=initial.values(u);psiu=initial.values(u,1);psiuu=initial.values(u,2)
    q=psi*psi/H;qu=2*psi*psiu/H;quu=2*(psiu*psiu+psi*psiuu)/H
    v=1-u*u;rr=u*q-v*qu
    radial=(-q-3*u*qu+v*quu+(qu/q)*rr)/q**3
    angular=(q*q+2*u*q*qu-v*qu*qu)/q**4
    amplitude2=initial.epsilon**2*24/(kappa*(p.mass_squared/H**2+3))
    W=p.vacuum+.5*p.mass_squared*amplitude2*u*u
    return dict(q=q,rr=rr,kb_dot=angular-radial-kappa*W,
                d_dot=(2*qu*qu-q*quu)/q**4-kappa*amplitude2/q**2,
                amplitude_squared=amplitude2)


def initial_parent_mode_pump(initial,mode,mass,orders=(220,100)):
    """Initial same-mode Hamiltonian Hessian on the actual parent geometry.

    This is a fixed INITIAL canonical energy-basis matrix element. A cubic
    offdiagonal propagator is conditional on an initial low-mode vacuum;
    it is not an instantaneous-final particle count or a UV-complete state.
    The time-dependent isotropic chart has a shift/unitary connection.
    Its same-lambda pair contribution vanishes, not all intermode terms.
    """
    if not math.isfinite(mass) or mass<0:raise ValueError('finite nonnegative declared mass required')
    if len(orders)!=2 or any(isinstance(n,(bool,np.bool_)) or not isinstance(n,(int,np.integer)) or n<16 for n in orders):
        raise ValueError('two integer quadrature orders >=16 required')
    radial_order,integral_order=map(int,orders)
    u,w=roots_jacobi(radial_order,.5,.5);F,G=mode_radial_profiles(mode,u);density=F*F+G*G
    geometry=_initial_parent_geometry_accelerations(initial,u)
    x,wx=leggauss(integral_order);inner=.5*u[:,None]*(x+1)
    Z=.5*u*np.sum(wx[None,:]*_initial_parent_geometry_accelerations(initial,inner)['d_dot'],axis=1)
    Omega_ddot=-geometry['kb_dot']*geometry['q']-geometry['rr']*Z
    lam=float(mode['eigenvalue']);energy=math.hypot(mass,lam)
    if energy==0:raise ValueError('nonzero same-mode energy required')
    lambda_ddot=-lam*float(w@(Omega_ddot*density))
    moment=float(w@(geometry['q']*density*u*u))
    mass_ddot=-.5*initial.parent.kappa*mass*initial.parent.mass_squared*geometry['amplitude_squared']*moment
    geometric=mass/energy*lambda_ddot;radion=-lam/energy*mass_ddot;total=geometric+radion
    return dict(radial_order=radial_order,integral_order=integral_order,
                mode_normalization=float(w@(geometry['q']*density)),mode_u2=moment,
                spatial_eigenvalue=lam,mass=mass,energy=energy,
                lambda_ddot=lambda_ddot,projected_mass_ddot=mass_ddot,
                geometric_pair_Hessian=geometric,radion_pair_Hessian=radion,total_pair_Hessian=total,
                initial_reference_amplitude_cubic_coefficient=total/6,
                normalized_singlet_probability_tau6_coefficient=total*total/18,
                instantaneous_final_basis_rotation_quadratic_coefficient=total/(4*energy),
                fixed_initial_beta='-i*tau³*M_pn/6+O(tau4), in the fixed initial canonical energy basis',
                instantaneous_final_beta='tau²*M_pn/(4E)+O(tau³), from rotation of the final particle basis; distinct diagnostic',
                geometry_rule='Z(u)=integral_0^u d_dot(v)dv; Omega_ddot=-kb_dot*q0-rr0*Z',
                mass_rule='m_ddot_bar=-(kappa/2)*m*mu_chi²*epsilon²*A²*<u²>, with physical Omega norm',
                chart_scope='Same-lambda pair is unchanged by the chart/spin-frame unitary; no permission to omit its shift for arbitrary intermode evolution.',
                low_mode_initial_vacuum_postulate_required=True,UV_Hadamard_state_selected=False,
                whole_field_occupation_selected=False,renormalized_stress_computed=False,baryon_photon_eta=None,
                scope='Initial same-mode fixed-reference pump only, with an added probe and conditional initial vacuum; no source backreaction or later detector prescription.')


def slice_isotropic_map(system,state,target_u,*,integral_order=128):
    """Exact slice-isomap formula with a declared numerical C2 integral.

    For h=q²dchi²+b²sin²chi dOmega² and q=b+(1-u²)c,
    J=integral_0^u(c/b)du, v=(u+tanh J)/(1+u*tanh J),
    Omega_new=b(cosh J+u sinh J). This is a spatial diffeomorphism of
    the actual supplied slice, not an average or a new zero-shift evolution.
    """
    u=np.atleast_1d(np.asarray(target_u,dtype=float))
    if u.ndim!=1 or np.any(~np.isfinite(u)) or np.any(abs(u)>=1):
        raise ValueError('finite interior isotropic target coordinates required')
    if isinstance(integral_order,bool) or not isinstance(integral_order,(int,np.integer)) or integral_order<16:
        raise ValueError('integer integral order >=16 required')
    C=system.coefficients(state)
    def values(at):
        at=np.atleast_1d(at)
        V=np.column_stack([eval_gegenbauer(l,1,at) for l in range(system.modes)])
        return C@V.T
    def ratio(at):
        b,c=values(at)[:2]
        if np.any(b<=0):raise ValueError('positive angular scale required on the supplied slice')
        return c/b
    coefficients=chebinterpolate(ratio,int(integral_order));integral=chebint(coefficients)
    offset=chebval(0.,integral);old=u.copy()
    for iteration in range(12):
        field=values(old);b,c=field[:2];J=chebval(old,integral)-offset
        g=np.cosh(J)+old*np.sinh(J)
        mapped=(old+np.tanh(J))/(1+old*np.tanh(J))
        derivative=(1+(1-old*old)*c/b)/(g*g)
        if np.any(derivative<=0):raise ValueError('positive radial scale/isomap monotonicity ended')
        residual=mapped-u
        if np.max(abs(residual))<4e-15:break
        old-=residual/derivative
        if np.any(abs(old)>=1):raise ValueError('isomap inversion left the supplied sphere chart')
    else:raise RuntimeError('bounded exact-slice isomap inversion did not converge')
    field=values(old);b,c=field[:2];J=chebval(old,integral)-offset
    g=np.cosh(J)+old*np.sinh(J);mapped=(old+np.tanh(J))/(1+old*np.tanh(J))
    derivative=(1+(1-old*old)*c/b)/(g*g);Omega=b*g;q=b+(1-old*old)*c
    angular=Omega**2*(1-mapped*mapped)/(b*b*(1-old*old))
    radial=Omega**2*derivative**2/(1-mapped*mapped)/(q*q/(1-old*old))
    return dict(original_u=old,isotropic_u=u,Omega=Omega,J=J,dv_du=derivative,
                mapping_residual=float(np.max(abs(mapped-u))),
                angular_metric_identity_relative_error=float(np.max(abs(angular-1))),
                radial_metric_identity_relative_error=float(np.max(abs(radial-1))),
                integral_order=int(integral_order),inversion_iterations=iteration,
                scope='Exact diffeomorphism formula, numerical Chebyshev integral/inversion; integration-order checks remain necessary.')


def finite_parent_mode_pump_check(system,state,proper_time,initial,mode,mass,*,radial_order=220,isomap_order=128):
    """Independent finite-geometry difference on the supplied actual slice."""
    proper_time=float(proper_time)
    if not math.isfinite(proper_time) or proper_time<=0:raise ValueError('positive finite proper time required')
    if not math.isfinite(mass) or mass<0:raise ValueError('finite nonnegative declared mass required')
    u,w=roots_jacobi(radial_order,.5,.5);F,G=mode_radial_profiles(mode,u);density=F*F+G*G
    q0=initial.values(u)**2/initial.parent.H
    S=system.jets(state,u)[0][4]
    mass_profile=mass*np.exp(-7*math.sqrt(system.kappa)*S/math.sqrt(84))
    mass_ddot=2*(float(w@(q0*density*mass_profile))-mass)/proper_time**2
    mapping=slice_isotropic_map(system,state,u,integral_order=isomap_order)
    A,M,*_=radial_operator_blocks(int(mode['degree']),int(mode['k']),u,w,mapping['Omega'])
    spectrum=eigh(A,M,eigvals_only=True)
    sign=1 if float(mode['eigenvalue'])>0 else -1
    choices=spectrum[spectrum*sign>0];lam=float(choices[np.argmin(abs(choices))])
    lam0=float(mode['eigenvalue']);energy=math.hypot(mass,lam0)
    lambda_ddot=2*(lam-lam0)/proper_time**2
    return dict(proper_time=proper_time,lambda_value=lam,lambda_ddot_finite=lambda_ddot,
                projected_mass_ddot_finite=mass_ddot,pair_Hessian_finite=(mass*lambda_ddot-lam0*mass_ddot)/energy,
                mapping_residual=mapping['mapping_residual'],
                angular_metric_identity_relative_error=mapping['angular_metric_identity_relative_error'],
                radial_metric_identity_relative_error=mapping['radial_metric_identity_relative_error'],
                radial_order=radial_order,isomap_order=isomap_order,
                scope='Actual synchronous finite slice; intrinsic eigenvalue via exact slice isomap, mass is the fixed-initial-mode matrix element. No field is physically averaged.')
