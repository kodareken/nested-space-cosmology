"""Receiver-pressure gates and an explicitly added geometric spinor candidate.

The actual parent remains the two-real-scalar action. A Dirac bundle,
fermionic statistics and a state do not follow by relabeling its scalar
harmonics. The spectral probe and local 4D EC benchmark below declare these
extra assumptions; neither is a formed child or a full 16D torsion reduction.
"""
import itertools
import math
from fractions import Fraction

import numpy as np

from .nsc_atlas_inherited_charge import harmonic_sector
from .nsc_atlas_parent_circulation import Parent


def radial_source_gate(*,radius,radius_rate,density,radial_pressure,
                       radial_spacetime_Hessian_over_radius,kappa=1.,actual_acceleration=None):
    """Necessary local force/null-stress requirements, not a matched new source.

    R_tau_tau=-(1+R_tau²+kappa R² p_r)/(2R) on the reflection sphere.
    The separate fixed-geometry null gate includes the spacetime Hessian.
    Adding a source also changes the constraints and geometry; the gate is
    not permission to hold those fixed in a subsequent Einstein solution.
    """
    values=(radius,radius_rate,density,radial_pressure,radial_spacetime_Hessian_over_radius,kappa)
    if any(not math.isfinite(float(v)) for v in values) or radius<=0 or kappa<=0:
        raise ValueError('finite data with positive radius and kappa required')
    R=float(radius);U=float(radius_rate);rho=float(density);p=float(radial_pressure)
    capacity=-kappa*R*R*p;required=1+U*U
    acceleration=-(required-capacity)/(2*R)
    nec=rho+p
    null_limit=-nec-2*radial_spacetime_Hessian_over_radius/kappa
    out=dict(radius=R,radius_rate=U,density=rho,radial_pressure=p,
        pressure_capacity=capacity,required_for_positive_acceleration=required,
        capacity_over_required=capacity/required,positive_acceleration_supported=capacity>required,
        acceleration_from_radial_Einstein_equation=acceleration,
        pressure_for_positive_acceleration_must_be_less_than=-required/(kappa*R*R),
        pressure_at_turning_point_must_be_less_than=-1/(kappa*R*R),
        extra_pressure_at_present_point_must_be_less_than=-required/(kappa*R*R)-p,
        classical_radial_NEC=nec,
        extra_total_NEC_at_fixed_geometry_must_be_less_than=null_limit,
        null_source_includes_free_and_interaction_stress=True,
        assumptions='Reflection sphere; actual state/constraints retained. A later geometry change can alter the gate. Values are in retained action units.')
    if actual_acceleration is not None:
        if not math.isfinite(float(actual_acceleration)):raise ValueError('finite actual acceleration required')
        out['actual_acceleration']=float(actual_acceleration)
        out['radial_equation_residual']=abs(acceleration-actual_acceleration)
    return out


def geometric_spinor_candidate(levels=3):
    """Quadratic massless-complex-Dirac probe on the declared round S12.

    This adds a field, not a fitted mass parameter. Both internal Dirac
    signs are retained. Gauge/mass ratios use the existing KK normalization;
    the internal U(1) is not identified with electromagnetism or baryons.
    """
    if isinstance(levels,bool) or not isinstance(levels,int) or not 1<=levels<=20:
        raise ValueError('integer levels in [1,20] required')
    scalar=harmonic_sector();d=scalar['d'];s=scalar['s_star'];a2=scalar['internal_radius_squared']
    first_mass2=d*d/(4*a2)*math.exp(-(d+2)*s)
    g2_over_g=2/(math.exp((d+2)*s)*(2*a2/(d+1)))
    modes=[dict(n=n,internal_absolute_eigenvalue=(n+d/2)/math.sqrt(a2),
                mass_squared=(n+d/2)**2/a2*math.exp(-(d+2)*s),
                multiplicity_per_internal_Dirac_sign=2**(d//2)*math.comb(n+d-1,d-1))
           for n in range(levels)]
    weights=list(itertools.product((Fraction(-1,2),Fraction(1,2)),repeat=d//2))
    ratio=Fraction(d,4);gauge_ratio=Fraction(4*(d+1),d*d);charged_ratio=gauge_ratio/4
    return dict(candidate='Added complex massless bulk Dirac probe; round internal sphere',
        existing_scalar_field_turned_into_fermions=False,added_field_postulate=True,
        round_shape_already_declared=True,round_shape_selected_by_inheritance=False,
        bulk_Dirac_mass_postulated_zero=True,fermion_statistics_imported=True,
        d=d,s_star=s,internal_scalar_curvature=scalar['internal_R_gamma'],
        internal_radius_squared=a2,internal_radius_squared_exact='44/5',
        spinor_canonical_factor='Psi16=exp(3s)*sum(psi_n tensor eta_n)',
        graded_operator='D16=exp(6s)(D4-3 slash(ds)) tensor1 + exp(-s)Gamma4,* tensorD_gamma',
        mass_law='m_n=abs(lambda_n)*exp(-7s)',
        modes=modes,lowest_mass_squared=first_mass2,lowest_mass=math.sqrt(first_mass2),
        lowest_mass_over_H=math.sqrt(first_mass2)/scalar['H'],
        logarithmic_mass_radion_derivative=-7/math.sqrt(84),
        mass_radion_derivative=-7/math.sqrt(84)*math.sqrt(first_mass2),
        Dirac_to_scalar_KK_increment_squared_ratio=float(ratio),
        Dirac_to_scalar_KK_increment_squared_ratio_exact=str(ratio),
        direct_scalar_gap_crosscheck=first_mass2/scalar['KK_mass_squared_increment'],
        gauge_squared_over_g_grav=g2_over_g,
        gauge_squared_over_kappa_Dirac_mass_squared=float(gauge_ratio),
        gauge_squared_over_kappa_Dirac_mass_squared_exact=str(gauge_ratio),
        minimal_charge_squared_over_kappa_Dirac_mass_squared=float(charged_ratio),
        minimal_charge_squared_over_kappa_Dirac_mass_squared_exact=str(charged_ratio),
        minimal_charge_coupling_squared_over_g_grav=g2_over_g/4,
        KK_gauge_charge_is_not_membrane_four_form_charge=True,
        internal_U1_identified_as_electromagnetism_or_baryons=False,
        ground_charge_multiplicities_per_Dirac_sign={
            'minus_half':sum(w[0]<0 for w in weights),'plus_half':sum(w[0]>0 for w in weights)},
        lowest_four_dimensional_Dirac_fields_both_internal_signs=2*2**(d//2),
        no_untwisted_massless_mode=True,
        Friedrich_internal_eigenvalue_squared_lower_bound=d*scalar['internal_R_gamma']/(4*(d-1)),
        symmetry='Proper rotation flipping internal n2,n13 exchanges charge weights and preserves round geometry and oriented flux.',
        invariant_state_charge_mean=0.,physical_baryon_photon_eta=None,
        g_grav_selected=False,SI_mass_selected=False,
        nonlinear_KK_and_torsion_truncation_proved=False,
        source_population_selected=False,regular_child_proved=False)


def mmp_massless_field_gate(spinor=None,magnetic_integer=1):
    """Counterfactual compatibility check, not an MMP solution of this parent.

    MMP r_e²=pi Q² G/e_min², kappa=8piG. The hypothetical gauge sector
    uses the actual smallest KK charge, so changing charge units cannot
    evade Dirac quantization. This checks the imported massless-field
    construction, not every possible massive wormhole model.
    """
    if isinstance(magnetic_integer,bool) or not isinstance(magnetic_integer,int) or not magnetic_integer:
        raise ValueError('nonzero magnetic integer required')
    spinor=spinor or geometric_spinor_candidate()
    ratio=spinor['minimal_charge_squared_over_kappa_Dirac_mass_squared'];Q=abs(magnetic_integer)
    square=Q*Q/(8*ratio)
    return dict(magnetic_integer=magnetic_integer,lowest_mass_times_extremal_radius_squared=square,
        lowest_mass_times_extremal_radius=math.sqrt(square),
        exact_squared_relation='m0²*r_e²=(18/13)*Q²',
        parametric_massless_support_limit_available=False,
        original_MMP_massless_field_assumption_satisfied=False,
        scalar_only_S3_parent_Maxwell_flux_through_separating_sphere=0.,
        missing_current_ingredients=['Maxwell field and magnetic topology/source','massless charged spinor sector','physical return state/loop'],
        scope='Additional operator/domain gate even if an external magnetic RN sector were provided. Not a no-go for other massive or twisted-bundle sources.')


def killing_ground_contact_block(r,signs,*,internal_dimension=12):
    """Exact ordered-index bare EC flavor tensor in Killing trivialization.

    eta_plus(v,n)=v/sqrt(Vgamma), eta_minus=Gamma(n)v/sqrt(Vgamma).
    Cq sums ALL ordered ambient Clifford q-indices. Signs name the left
    and right fields of two bilinears. Chiral mass rotations are separate.
    """
    if (isinstance(r,bool) or not isinstance(r,int) or not 0<=r<=3 or
            isinstance(internal_dimension,bool) or not isinstance(internal_dimension,int) or
            internal_dimension<2 or internal_dimension%2):
        raise ValueError('r in [0,3] and positive even internal dimension>=2 required')
    if len(signs)!=4 or any(v not in (-1,1) or isinstance(v,bool) for v in signs):
        raise ValueError('four Killing signs +/-1 required')
    a,b,c,d=signs;N=internal_dimension+1
    diagonal1=a==b;diagonal2=c==d
    if diagonal1!=diagonal2:rank=0;coefficient=Fraction(0)
    elif diagonal1:
        rank=(0,2,2,4)[r]
        coefficient=(Fraction(1),Fraction(-a*c,N),Fraction(N-2,N),Fraction(-a*c,N))[r]
    else:
        rank=(1,1,3,3)[r]
        coefficient=(Fraction(1,N),Fraction(-a*c*(N-1),N),Fraction(1,N),Fraction(-a*c*(N-3),N))[r]
    if rank>N:coefficient=Fraction(0)
    return dict(r=r,Killing_signs=list(signs),ambient_Clifford_rank=rank,
                ordered_tensor_coefficient=str(coefficient),coefficient=float(coefficient),
                external_three_form_weight=math.comb(3,r)*(-1)**r)


def complete_ground_contact_vertex():
    """All bare ground-sign couplings; not an occupied or conserved truncation."""
    return dict(added_action='Minimal 16D Einstein-Cartan plus complex massless Dirac on the declared round internal sphere',
        universal_contact='L16=-kappa16*J_ABC*J^ABC/32; J=barPsi Gamma_[ABC] Psi',
        normalization='Psi16=exp(3s)*sum psi4*eta; quartic/volume warps cancel, kappa16/Vgamma=kappa4',
        ambient_tensor='Cq[v,w;x,z]=sum_(I1...Iq ordered) Gamma_[I1...Iq][v,w]*Gamma_[I1...Iq][x,z]',
        bare_bilinear='B_r=(barpsi_(a,v) gamma_[mu1...mu(3-r)] gamma5^r psi_(b,w))',
        reduced_action='L4=-kappa4/32*sum_r choose(3,r)(-1)^r*B_r^(av,bw)*B_r^(cx,dz)*K_r^(ab,cd)[vw;xz]',
        blocks=[killing_ground_contact_block(r,signs) for r in range(4) for signs in itertools.product((1,-1),repeat=4)],
        positive_mass_rotation='psi_a=exp(-i*a*pi*gamma5/4)*chi_a; transform every bilinear with its own left/right signs',
        one_sign_positive_mass_action='-kappa4/32*[6A²-(3/13)sum_IJ T_IJ,munu T_IJ^munu+(33/13)sum_IJ V_IJ,mu V_IJ^mu-(1/13)sum_IJKL S_IJKL²]',
        all_internal_sums_ordered=True,
        ground_vertex_projection_derived=True,independent_contact_parameter_added=False,
        nonlinear_KK_truncation_proved=False,occupied_state_selected=False,
        radial_stress_computed=False,regular_child_proved=False,
        scope='Bare tree-level matrix-element projection of the added action, both Killing signs. Higher modes, state transport, composite renormalization and backreaction remain separate.')


def signed_gaussian_coincidence(covariance,*,density=1.):
    """Complete a declared paired Gaussian equal-time covariance difference.

    Input C is traceless Hermitian2x2 Delta<psi psi^dagger> in the positive
    mass tau basis, after its actual initial energy rotation. The opposite
    spatial branch is tau3*C*tau3, with equal local mode density. Both
    complete internal Killing-sign populations are IDENTICAL in this
    positive-mass convention. This restricted preparation is not selected
    by C/P symmetry or by the geometry alone.

    The result is the vacuum-Wick coincidence difference, not kinetic
    stress, a full renormalized source, or a positive covariance by itself.
    """
    C=np.asarray(covariance,dtype=complex)
    if (C.shape!=(2,2) or np.any(~np.isfinite(C)) or
            not math.isfinite(density) or density<0):
        raise ValueError('finite2x2 covariance difference and nonnegative density required')
    scale=max(1.,float(np.linalg.norm(C)))
    if np.max(abs(C-C.conj().T))>64*np.finfo(float).eps*scale or abs(np.trace(C))>64*np.finfo(float).eps*scale:
        raise ValueError('Hermitian traceless covariance difference required')
    tau3=np.diag([1.,-1.]);eye=np.eye(4,dtype=complex)
    gamma5=np.kron(np.array([[0.,1.],[1.,0.]]),np.eye(2))
    external=density*np.kron(C@tau3+tau3@C,np.eye(2))
    rotations=[(eye-1j*a*gamma5)/math.sqrt(2) for a in (-1,1)]
    complete=sum(R@external@R for R in rotations)
    return dict(external_normal_covariance=external,complete_normal_covariance=complete,
                external_norm=float(np.linalg.norm(external)),complete_norm=float(np.linalg.norm(complete)),
                scalar_coefficient=float((2*density*C[0,0]).real),
                scope='Identically prepared complete internal signs and paired spatial branch; fixed-reference Gaussian difference only.')


def scalar_coherent_gaussian_control(c=.1,*,D=16,density=1.):
    """Local finite-CAR counterexample, not a selected global parent state.

    K_ref=I/2 and K=I/2+c*gamma0 are valid mixed Gaussian covariances
    for |c|<=1/2. Their relative coincidence is c*density*I, outside the
    identical-positive-mass internal populations in the previous function.
    The local scalar Clifford component is C/P even, but chiral coherence
    is nonzero. No global C/P/Hadamard preparation follows from this
    local matrix; the half-filled reference is not a vacuum. Density is
    only a profile normalization in this algebraic control, not a ruler.
    """
    if (isinstance(D,bool) or not isinstance(D,int) or D<4 or D>16 or D%2 or
            not math.isfinite(c) or abs(c)>.5 or not math.isfinite(density) or density<0):
        raise ValueError('evenD in[4,16], |c|<=1/2 and nonnegative finite density required')
    ns=2**(D//2);eye=np.eye(ns);half=np.eye(ns//2)
    gamma0=np.kron(np.diag([1.,-1.]),half)
    grading=np.kron(np.array([[0.,1.],[1.,0.]]),half)
    reference=.5*eye;covariance=reference+c*gamma0
    normal=c*density*eye
    xi=covariance@grading-grading@covariance
    rank3_trace=ns*D*(D-1)*(D-2)
    wick=-rank3_trace*(c*density)**2
    return dict(covariance=covariance,reference_covariance=reference,
                gamma0=gamma0,chirality=grading,normal_covariance=normal,chiral_commutator=xi,
                covariance_eigenvalues=np.repeat([.5-c,.5+c],ns//2),
                chiral_commutator_norm=float(np.linalg.norm(xi)),
                Clifford_dimension=D,spinor_dimension=ns,
                rank3_trace_coefficient=rank3_trace,literal_Wick_three_form=wick,
                contact_lagrangian_over_kappa=-wick/32,
                contact_coefficient_per_squared_scalar_coincidence=rank3_trace/32,
                reference_is_vacuum=False,global_parent_state_selected=False,UV_completion_selected=False,
                radial_NEC_computed=False,
                scope='Local normalized finite-CAR scalar/chiral-coherence witness, not a global C/P invariant Hadamard state or an inherited occupation.')


def parent_dirac_source_operator(initial,degree,*,mass,k=0,quadrature_order=None):
    """Full initial radial Dirac Hamiltonian Hessian in a fixed proper chart.

    No same-mode closure or time-dependent isotropic-chart shift is used.
    Return a compact diagnostic and in-memory matrices. This is a finite
    initial operator, not a selected occupied state or renormalized stress.
    All spatial modes within the specified k/cutoff block are retained.
    """
    from scipy.linalg import eigh
    from scipy.special import roots_jacobi,eval_gegenbauer
    from .nsc_atlas_s3_spectral import _radial_normalization
    from .nsc_atlas_valve_dirac_modes import (
        radial_operator_blocks,_initial_parent_geometry_accelerations,initial_parent_mode_pump)
    if (isinstance(degree,(bool,np.bool_)) or not isinstance(degree,(int,np.integer)) or
            not 2<=degree<=32 or isinstance(k,(bool,np.bool_)) or
            not isinstance(k,(int,np.integer)) or not 0<=k<degree or
            not math.isfinite(mass) or mass<0):
        raise ValueError('degree2..32, integer0<=k<degree and finite nonnegative mass required')
    order=2*degree+2*initial.degree+35 if quadrature_order is None else quadrature_order
    if isinstance(order,(bool,np.bool_)) or not isinstance(order,(int,np.integer)) or order<2*degree+5:
        raise ValueError('integer quadrature order>=2degree+5 required')
    u,w=roots_jacobi(order,.5,.5);v=1-u*u;t=np.sqrt(v)
    psi=initial.values(u);q=psi*psi/initial.parent.H
    qu=2*psi*initial.values(u,1)/initial.parent.H
    A,M,*_=radial_operator_blocks(degree,k,u,w,q)
    labels=(list(range(k,degree+1)),list(range(k+1,degree+1)))
    n0=len(labels[0]);n=len(A)
    F=np.zeros((len(u),n));G=np.zeros_like(F);Fu=np.zeros_like(F);Gu=np.zeros_like(F)
    for branch,(j,ls,values,deriv) in enumerate(((k,labels[0],F,Fu),(k+1,labels[1],G,Gu))):
        offset=0 if branch==0 else n0
        for h,l in enumerate(ls):
            factor=_radial_normalization(l,j);C=eval_gegenbauer(l-j,j+1,u)
            values[:,offset+h]=factor*v**(j/2)*C
            deriv[:,offset+h]=factor*((-j*u*v**(j/2-1)*C if j else 0)+
                (2*(j+1)*v**(j/2)*eval_gegenbauer(l-j-1,j+2,u) if l>j else 0))
    # The polar radial basis is O=[[c,-s],[-s,-c]], a reflection.
    # The old right-frame convention F-i G sigma.n is unchanged.
    c=np.sqrt((1+u)/2);s=np.sqrt((1-u)/2)
    P0=c[:,None]*F-s[:,None]*G;P1=-s[:,None]*F-c[:,None]*G
    dP0=-.5*s[:,None]*F-.5*c[:,None]*G-t[:,None]*(c[:,None]*Fu-s[:,None]*Gu)
    dP1=-.5*c[:,None]*F+.5*s[:,None]*G+t[:,None]*(s[:,None]*Fu+c[:,None]*Gu)
    scale=np.sqrt(q)*t;logder=-t*qu/(2*q)+u/t
    g0=scale[:,None]*P0;g1=scale[:,None]*P1
    dg0=scale[:,None]*dP0+logder[:,None]*g0;dg1=scale[:,None]*dP1+logder[:,None]*g1
    wc=w/v
    gram=g0.T@(wc[:,None]*g0)+g1.T@(wc[:,None]*g1)
    def radial(B):
        weighted=(wc*B)[:,None]
        return .5*(-g0.T@(weighted*dg1)+g1.T@(weighted*dg0)+
                   dg0.T@(weighted*g1)-dg1.T@(weighted*g0))
    def angular(B):
        weighted=(wc*B)[:,None]
        return g0.T@(weighted*g1)+g1.T@(weighted*g0)
    original=radial(1/q)+angular(-(k+1)/(q*t))
    jets=_initial_parent_geometry_accelerations(initial,u)
    ka_dot=jets['kb_dot']+v*jets['d_dot']
    Ddd=radial(ka_dot/q)+angular(-(k+1)*jets['kb_dot']/(q*t))
    md=-.5*initial.parent.kappa*mass*initial.parent.mass_squared*jets['amplitude_squared']*u*u
    Mdd=g0.T@((wc*md)[:,None]*g0)+g1.T@((wc*md)[:,None]*g1)
    lam,C=eigh(A,M);Dp=C.T@Ddd@C;Mp=C.T@Mdd@C
    energy=np.hypot(mass,lam)
    aa=np.sqrt((energy+mass)/(2*energy));bb=lam/np.sqrt(2*energy*(energy+mass))
    pair=(aa[:,None]*aa[None,:]-bb[:,None]*bb[None,:])*Dp-(aa[:,None]*bb[None,:]+bb[:,None]*aa[None,:])*Mp
    choices=np.flatnonzero(lam>0);index=int(choices[np.argmin(lam[choices])])
    mode=dict(degree=int(degree),k=int(k),sign=1,eigenvalue=float(lam[index]),coefficients=C[:,index])
    independent=initial_parent_mode_pump(initial,mode,mass)
    column=float(np.sum(pair[:,index]**2));diagonal=float(pair[index,index]**2)
    strong=sorted([dict(spatial_eigenvalue=float(lam[j]),pair_Hessian=float(pair[j,index]))
                   for j in range(len(lam)) if j!=index],key=lambda r:abs(r['pair_Hessian']),reverse=True)[:5]
    summary=dict(degree=int(degree),k=int(k),angular_multiplicity=2*k+2,quadrature_order=int(order),
        spatial_block_dimension=n,ground_spatial_eigenvalue=float(lam[index]),
        Gram_identity_error=float(np.max(abs(gram-M))),
        original_Dirac_projection_error=float(np.max(abs(original-A))),
        Hessian_hermiticity_error=float(np.max(abs(Ddd-Ddd.T))),
        same_mode_pair_Hessian=float(pair[index,index]),
        same_mode_independent_isotropic_HF_value=independent['total_pair_Hessian'],
        same_mode_independent_difference=abs(float(pair[index,index])-independent['total_pair_Hessian']),
        ground_column_squared_coefficient_norm=column,
        same_mode_squared_coefficient_fraction=diagonal/column if column else None,
        other_mode_squared_coefficient_fraction=1-diagonal/column if column else None,
        strongest_other_mode_coefficients=strong,
        scope='Finite initial k-block Hamiltonian Hessian in the fixed proper chart; fractions are reference-amplitude coefficient norms, not energy fractions, physical particle counts or renormalized NEC.',
        whole_field_state_selected=False,renormalized_stress_computed=False,child_formed=False)
    matrices=dict(spatial_eigenvalues=lam,mass_hessian=Mp,spatial_Dirac_hessian=Dp,
                  positive_negative_pair_hessian=pair,spatial_eigenvectors=C,ground_index=index)
    return summary,matrices


def dirac_clifford_control():
    """Literal four-component Clifford/Wick normalization in mostly-plus units."""
    zero=np.zeros((2,2),complex);I=np.eye(2,dtype=complex)
    pauli=(np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1.,-1.]).astype(complex))
    gamma=[np.block([[I,zero],[zero,-I]])]+[np.block([[zero,sigma],[-sigma,zero]]) for sigma in pauli]
    eta=np.diag([-1.,1.,1.,1.]);eye=np.eye(4,dtype=complex)
    defect=max(np.max(abs(gamma[a]@gamma[b]+gamma[b]@gamma[a]+2*eta[a,b]*eye)) for a in range(4) for b in range(4))
    gamma5=1j*gamma[0]@gamma[1]@gamma[2]@gamma[3]
    axial=[g@gamma5 for g in gamma]
    trace=sum(eta[a,a]*np.trace(axial[a]@axial[a]) for a in range(4))
    return dict(metric='-+++',Clifford_relation='{gamma^a,gamma^b}=-2 eta^ab',
                Clifford_defect=float(defect),axial_pair_trace=float(trace.real),
                imaginary_trace_residual=float(abs(trace.imag)),
                literal_Wick_A_squared_over_m_squared_sigma_squared=-float(trace.real)/16)


def neutral_thermal_ec_benchmark(spinor=None):
    """Explicit 4D EC thermal-difference benchmark; not a total 16D source.

    No physical temperature, occupation or gravitational normalization is
    selected. hat(rho)=kappa_nat*Lunit²*rho_phys and g=kappa_nat/Lunit².
    The free and interaction NEC are both included. The literal operator
    convention differs from a half-weighted Wick-diagram convention.
    """
    spinor=spinor or geometric_spinor_candidate();m2=spinor['lowest_mass_squared']
    return dict(state='Local free unpolarized C-neutral Gaussian, mu=0; vacuum Wick thermal-difference prescription',
        state_selected_by_parent=False,background='Local thermal/adiabatic benchmark, not the evolving neck state',
        action='L_contact=-(3*kappa_nat/16)*A_mu*A^mu in the declared -+++ Clifford convention',
        literal_contraction='A²_TT=-m² sigma_T², sigma_T=-T²/6 at T>>m',
        contact_pressure_hat_over_g_squared_T_hat_fourth=m2/192,
        contact_density_hat_over_g_squared_T_hat_fourth=m2/64,
        contact_NEC_hat_over_g_squared_T_hat_fourth=m2/48,
        free_NEC_hat_over_g_T_hat_fourth=7*math.pi**2/45,
        contact_to_free_NEC_ratio_over_g=15*m2/(112*math.pi**2),
        reference_half_Wick_ratio_over_g=15*m2/(224*math.pi**2),
        contact_NEC_sign='positive',total_NEC_sign_in_controlled_regime='positive',
        supports_minimum_neck_turnaround=False,
        controlled_regime='m_hat << T_hat and g*T_hat² <<1; hence g*m_hat² <<1',
        quantum_gravitational_normalization='g=kappa_nat/Lunit² remains free; arithmetic kappa=1 is not a choice g=1',
        g_selected=False,temperature_selected=False,
        C_even_contact_selects_eta=False,current_scalar_parent_fermionic_contact=0.,
        benchmark_is_not_the_full_16d_contact=True,
        required_next_source='N(τ)=k^mu*k^nu <T_Dirac+interaction_mu_nu>_(inherited state,ren), including kinetic stress and backreaction.',
        missing_state_data='Actual inherited covariance and connected four-point function, with renormalization and all retained torsion channels.')


def parent_mode_axial_kernel(*,mass,spatial_eigenvalue,density_overlap):
    """Bare external-axial contact difference on an actual k=0 parent mode.

    The mode field is sqrt(rho(x))*(Q_tau tensor U_spin(x))*psi_rest,
    Q=[[A,-B],[B,A]], A²=(E+m)/(2E). The singlet/triplet states both
    have N=2 and identical quadratic covariance. Their local quartic
    difference is [8+4(m/E)²]*rho². This is one torsion channel, not a
    total 16D interaction, pressure or claimed inherited occupation.
    """
    if any(not math.isfinite(float(v)) for v in (mass,spatial_eigenvalue,density_overlap)) or mass<0 or density_overlap<=0:
        raise ValueError('finite mass/eigenvalue/overlap with mass>=0 and overlap>0 required')
    energy=math.hypot(mass,spatial_eigenvalue)
    if energy==0:raise ValueError('positive mode energy required')
    fraction=(mass/energy)**2;quartic=8+4*fraction
    return dict(mass=mass,spatial_eigenvalue=spatial_eigenvalue,energy=energy,
        density_overlap=density_overlap,effective_volume=1/density_overlap,
        mass_squared_over_energy_squared=fraction,
        local_quartic_gap_over_density_squared=quartic,
        integrated_quartic_gap=quartic*density_overlap,
        triplet_minus_singlet_energy_hat_over_g=3*quartic*density_overlap/16,
        singlet_contact_energy_hat_over_g=3*(-8-2*fraction)*density_overlap/16,
        isotropic_triplet_contact_energy_hat_over_g=3*(2*fraction)*density_overlap/16,
        normalization='Delta(E*Lunit)=g*(3/4)*(2+m_hat²/E_hat²)*I4_hat; g=kappa_nat/Lunit² remains free.',
        source_state_selected=False,full_torsion_channels_matched=False,
        physical_pressure_or_NEC_derived=False,
        scope='One explicit initial-parent angular doublet and one external axial contact channel; finite-mode state comparison, not the total interacting ground state.')


def axial_pair_correlation_witness():
    """Exact finite CAR states with equal quadratic data and different contact.

    One orbital/rest-frame witness only: no actual curved-space mode,
    inherited state, pressure/NEC or full16D truncation is inferred.
    """
    import time
    import sympy as s
    from sympy.physics.matrices import mgamma

    def zero(matrix):return all(s.simplify(x)==0 for x in matrix)
    def scalar(value):
        value=s.simplify(value)
        if value.is_real is not True:raise ArithmeticError('nonreal expectation')
        return str(value)
    started=time.perf_counter();I=s.eye(2);Z=s.diag(1,-1);ann=s.Matrix([[0,1],[0,0]])
    car=[]
    for i in range(4):
        car.append(s.kronecker_product(*[Z if j<i else ann if j==i else I for j in range(4)]))
    create=[a.H for a in car];identity=s.eye(16);vac=s.zeros(16,1);vac[0]=1
    car_ok=all(zero(car[i]*create[j]+create[j]*car[i]-(identity if i==j else s.zeros(16)))
               and zero(car[i]*car[j]+car[j]*car[i]) for i in range(4) for j in range(4))
    assert car_ok
    singlet=(create[0]*create[3]-create[1]*create[2])*vac/s.sqrt(2)
    triplets=[create[0]*create[2]*vac,
              (create[0]*create[3]+create[1]*create[2])*vac/s.sqrt(2),
              create[1]*create[3]*vac]
    rhoS=singlet*singlet.H;rhoT=sum((v*v.H for v in triplets),s.zeros(16))/3
    rhoG=identity/16
    number=sum((create[i]*car[i] for i in range(4)),s.zeros(16))
    charge=sum((create[i]*car[i]*(1 if i<2 else -1) for i in range(4)),s.zeros(16))
    pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1)]
    spin=[]
    for sigma in pauli:
        spin.append(sum((create[base+i]*sigma[i,j]*car[base+j]/2
                         for base in (0,2) for i in range(2) for j in range(2)),s.zeros(16)))
    spin_square=sum((J*J for J in spin),s.zeros(16))

    # Second quantization of a<->b, with vacuum invariant and fermionic signs.
    permutation=(2,3,0,1);C=s.zeros(16)
    for column in range(16):
        occupied=[i for i in range(4) if (column>>(3-i))&1]
        mapped=vac
        for i in reversed(occupied):mapped=create[permutation[i]]*mapped
        C[:,column]=mapped
    assert zero(C.H*C-identity)
    assert all(zero(C*car[i]*C.H-car[permutation[i]]) for i in range(4))

    # Overall V^-1/2 is omitted in matrices and restored in all readouts.
    # psi=(a_up,a_down,-b_down^dagger,b_up^dagger)/sqrt(V).
    field_terms=[(1,('a',0)),(1,('a',1)),(-1,('c',3)),(1,('c',2))]
    field=[coefficient*(car[index] if kind=='a' else create[index])
           for coefficient,(kind,index) in field_terms]
    gamma5=mgamma(5);gamma0=mgamma(0);current=[];current_terms=[]
    for mu in range(4):
        matrix=gamma0*mgamma(mu)*gamma5;A=s.zeros(16);terms=[]
        for i in range(4):
            for j in range(4):
                if matrix[i,j]==0:continue
                A+=matrix[i,j]*field[i].H*field[j]
                ci,(ki,ii)=field_terms[i];cj,(kj,ij)=field_terms[j]
                adjoint_i=('c' if ki=='a' else 'a',ii)
                terms.append((ci*cj*matrix[i,j],(adjoint_i,(kj,ij))))
        current.append(A);current_terms.append(terms)
    AA=-current[0]**2+sum((A*A for A in current[1:]),s.zeros(16))
    assert all(zero(current[i+1]-2*spin[i]) for i in range(3))

    def fully_normal_ordered(tokens):
        tokens=list(tokens);sign=1
        # All creation operators move to the left; contractions are discarded
        # across the FULL four-operator product, not each current separately.
        for i in range(len(tokens)):
            for j in range(len(tokens)-1):
                if tokens[j][0]=='a' and tokens[j+1][0]=='c':
                    tokens[j],tokens[j+1]=tokens[j+1],tokens[j];sign=-sign
        product=sign*identity
        for kind,index in tokens:product=product*(car[index] if kind=='a' else create[index])
        return product

    normal=s.zeros(16)
    for mu in range(4):
        metric_sign=-1 if mu==0 else 1
        for c1,t1 in current_terms[mu]:
            for c2,t2 in current_terms[mu]:
                normal+=metric_sign*c1*c2*fully_normal_ordered(t1+t2)
    normal_order_identity=zero(normal-AA-2*identity+4*number)
    assert normal_order_identity

    def expectation(rho,operator):return s.simplify(s.trace(rho*operator))
    def inspect(rho):
        covariance=s.Matrix(4,4,lambda i,j:expectation(rho,create[i]*car[j]))
        anomalous=s.Matrix(4,4,lambda i,j:expectation(rho,car[i]*car[j]))
        assert covariance==s.eye(4)/2 and anomalous==s.zeros(4)
        raw=expectation(rho,AA);ren=expectation(rho,normal)
        return dict(trace=scalar(s.trace(rho)),N=scalar(expectation(rho,number)),Q=scalar(expectation(rho,charge)),
            spin_squared=scalar(expectation(rho,spin_square)),
            one_body_covariance='diag(1/2,1/2,1/2,1/2)',anomalous_covariance='zero',
            raw_AA_times_V_squared=scalar(raw),full_vacuum_Wick_AA_times_V_squared=scalar(ren),
            raw_contact_H_times_V_over_kappa=scalar(s.Rational(3,16)*raw),
            normal_contact_H_times_V_over_kappa=scalar(s.Rational(3,16)*ren),
            normal_contact_density_times_V_squared_over_kappa=scalar(s.Rational(3,16)*ren),
            free_vacuum_subtracted_energy_over_m=scalar(expectation(rho,number)),
            rotation_invariant=all(zero(rho*J-J*rho) for J in spin),
            C_invariant=zero(C*rho*C.H-rho),charge_commuting=zero(rho*charge-charge*rho))
    records={name:inspect(rho) for name,rho in (('singlet',rhoS),('isotropic_triplet',rhoT),('matching_gaussian',rhoG))}
    for rho in (rhoS,rhoT):
        assert zero(rho*number-2*rho)
        assert zero(rho*charge)
    assert all(zero(C*A*C.H-A) for A in current)
    assert all(zero(C*J-J*C) for J in spin)
    assert zero(C*charge*C.H+charge)
    sector=rhoS+3*rhoT
    assert zero(sector*sector-sector) and s.trace(sector)==4
    assert zero(AA*sector-sector*AA)
    assert zero(sector*AA*sector+4*rhoS-24*rhoT)
    difference=s.simplify(expectation(rhoT,AA)-expectation(rhoS,AA))
    assert difference==12
    out=dict(schema='nsc-ec-two-body-connected-correlation-witness-v1',
        CAR_exact=car_ok,CAR_modes=['a_up','a_down','b_up','b_down'],Hilbert_dimension=16,
        field='psi=(a_up,a_down,-b_down^dagger,b_up^dagger)/sqrt(V)',
        metric='-+++; Clifford {gamma,gamma}=-2*eta; gamma5^2=1',
        minimal_EC='L_contact=-3*kappa*A_a*A^a/16; H_contact=+3*kappa*integral A_a*A^a dV/16',
        full_normal_order_identity=':A_a A^a:=A_a A^a+(2-4N)/V^2',
        full_normal_order_identity_exact=normal_order_identity,
        exact_fixed_N2_Q0_spectrum={'singlet':-4,'triplet_eigenvalue':8,'triplet_multiplicity':3},
        current_identity='A_i=2 J_i/V; A0=-sqrt(2)*(S_pair^dagger+S_pair)/V',
        states=records,
        triplet_minus_singlet_AA_times_V_squared=scalar(difference),
        triplet_minus_singlet_contact_H_times_V_over_kappa=scalar(s.Rational(3,16)*difference),
        connected_quartic_departures_from_matching_Gaussian={
            name:scalar(expectation(rho,AA)-expectation(rhoG,AA)) for name,rho in (('singlet',rhoS),('isotropic_triplet',rhoT))},
        common_lower_order_subtraction_cancels='Both states have identical expectations of every constant and quadratic polynomial; their quartic difference survives any common such subtraction.',
        scope=['Finite four-mode single-orbital rest-frame algebra; not an actual S3 Dirac mode or Hadamard continuum state.',
               'New 4D fermion extension; not derived from existing bosonic fields or a full 16D torsion projection.',
               'Singlet minimizes contact energy in this fixed-N=2,Q=0 sector; no pressure, NEC, continuum ground state or bounce is asserted.',
               'Both states have exactly zero charge; no matter-antimatter asymmetry is generated.'],
        elapsed_seconds=time.perf_counter()-started)
    return out
