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
        full_16d_contact_reduction_performed=False,
        required_next_source='N(τ)=k^mu*k^nu <T_Dirac+interaction_mu_nu>_(inherited state,ren), including kinetic stress and backreaction.',
        missing_state_data='Actual inherited covariance and connected four-point function, with renormalization and all retained torsion channels.')


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
