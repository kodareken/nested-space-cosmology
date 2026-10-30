"""Internal harmonic charge and its conditional inheritance across a horizon.

Round S12 is an explicit shape completion of the existing volume reduction.
The quadratic charged modes and Euclidean preparation are not a solved
higher-dimensional formation history or a baryon-sector identification.
"""
import math
from .nsc_atlas_parent_scale_selector import VolumeParent


def harmonic_sector(omega0=20.):
    model=VolumeParent();c=model.coefficients
    minimum=model.vacuum_roots()[1];s=minimum['s'];d=model.d
    curvature=2*float(c['B'])
    radius2=d*(d-1)/curvature
    eigenvalue=d/radius2
    base=omega0**2*math.exp(-d*s)
    shift=eigenvalue*math.exp(-(d+2)*s)
    mass2=base+shift;H=math.sqrt(minimum['vacuum']/3)
    return dict(d=d,ancestors=model.ancestor_steps,omega0=omega0,
        s_star=s,Lambda=minimum['vacuum'],H=H,
        internal_R_gamma=curvature,internal_radius_squared=radius2,
        internal_radius_squared_exact='44/5',ell1_eigenvalue=eigenvalue,
        ell1_eigenvalue_exact='15/11',ell1_real_degeneracy=d+1,
        Cartan_complex_pairs=d//2,Cartan_neutral_modes=1,
        zero_mode_mass_squared=base,KK_mass_squared_increment=shift,
        charged_mass_squared=mass2,charged_mass=math.sqrt(mass2),
        charged_mass_over_H=math.sqrt(mass2)/H,
        temperature=H/(2*math.pi),
        shape_selected_by_action=False,
        scope='Round internal sphere with curvature fixed by the retained B coefficient; linear harmonics on the effective late-vacuum background.')


def pair_transport(n_plus,n_minus,squeeze):
    """Diagonal input number state control for a neutral Bogoliubov pair."""
    if min(n_plus,n_minus,squeeze)<0:
        raise ValueError('nonnegative populations and finite squeeze required')
    a2=math.cosh(squeeze)**2;b2=math.sinh(squeeze)**2
    return dict(plus=a2*n_plus+b2*(n_minus+1),
                minus=a2*n_minus+b2*(n_plus+1),
                conserved_charge=n_plus-n_minus)


def thermal_log_ratio(gap,H):
    if gap<=0 or H<=0:
        raise ValueError('positive central-clock gap and H required')
    return -2*math.pi*gap/H


def calculation():
    h=harmonic_sector();M=h['charged_mass'];H=h['H']
    log_ratio=thermal_log_ratio(M,H)
    log_variance=math.log(2)+log_ratio-2*math.log1p(-math.exp(log_ratio))
    return dict(schema='nsc-atlas-inherited-charge-v1',
        sector=h,
        derivation=dict(
            metric='ds_D^2=exp(-d s) ds_4E^2+exp(2s) gamma_ab dy^a dy^b',
            curvature='V_curvature=-R_gamma exp[-(d+2)s]/(2 kappa4); B=15/2 implies R_gamma=15',
            normalized_harmonics='Y_i=sqrt((d+1)/V_gamma)*n_i, integral Y_i Y_j=delta_ij',
            field='deltaPhi=x Y_1+y Y_2=psi Y_plus^*+psi^* Y_plus; psi=(x+i y)/sqrt(2)',
            eigenvalue='-Delta_gamma Y_ell=ell(ell+d-1)/a0^2 Y_ell',
            mass='M_ell^2=omega0^2 exp(-d s)+ell(ell+d-1)/a0^2 exp[-(d+2)s]',
            action='S_2=-integral sqrt(-g) [|D psi|^2+M_1^2 |psi|^2]',
            charge='Internal rotational charges +/-1; j_mu=x D_mu y-y D_mu x',
            meaning='Geometric Kaluza-Klein charge, not yet baryon number; gauge fields from mixed metric components are required when charged modes backreact.',
            closure='Quadratic test sector only; finite excitation generally sources additional metric and harmonic modes.'),
        orientation_symmetry=dict(
            map='P=diag(1,-1,1,...,1,-1) on R^13; flip n2 and n13',
            determinant=1,charge_generator='P J12 P^-1=-J12',
            flux_preserved=True,
            consequence='A proper internal rotation conjugates the charge while preserving the round metric, homogeneous scalar and oriented top-form background. Nonzero background flux does not by itself select a charge sign.'),
        horizon=dict(
            radial_equation="f u''+(f'-2i omega)u'-(f'/r+M^2+l(l+1)/r^2)u=0",
            simple_root_indices='p=0 and p=2i omega/f_prime(rh)',
            regular_branch='Ingoing branch has arbitrary complex amplitude; regularity supplies no internal phase bias.',
            charges_have_same_operator=True,
            mode_frequency_selected=False,
            scope='Smooth neutral spherical horizon, constant stationary internal radius. A future collapse history was not specified.'),
        euclidean_preparation=dict(
            condition='Analytic continuation of the positive Gaussian operator on S4 at the late de Sitter vacuum, with internal-rotation invariant boundary preparation.',
            eigenvalues='H^2 L(L+3)+M1^2, L=0,1,...',
            lowest_eigenvalue=h['charged_mass_squared'],
            cross_covariance_xy=0.,equal_component_covariances=True,
            phase_distribution='Uniform on [0,2pi) for the isotropic complex Gaussian; no coherent ellipse offset.',
            coherent_phase_offset=None,mean_charge=0.,chemical_potential=0.,
            selection_status='Unique regular Euclidean Gaussian inverse under the stated preparation; horizon regularity alone does not select this state.',
            whole_universe_wavefunction_derived=False,
            de_Sitter_state_is_not_automatically_Kottler_state=True),
        inheritance=dict(
            pair_map='a_out=alpha a_in+beta* b_in^dagger; b_out=alpha b_in+beta* a_in^dagger',
            identity='N_plus_out-N_minus_out=(|alpha|^2-|beta|^2)(N_plus_in-N_minus_in)=Q_in',
            arbitrary_depth='A charge-conjugation-equivariant channel takes an invariant input state to an invariant state at every depth; regional <Q_n>=0 when the internal conjugation preserves the region.',
            regional_flux='Regional charge may change by current through the boundary; it is not set to zero merely because total charge is zero.',
            neutral_controls=[dict(squeeze=x,**pair_transport(0.,0.,x)) for x in (.1,1.,2.)],
            charged_control=dict(squeeze=1.,**pair_transport(2.,1.,1.)),
            controls_are_observation_fits=False),
        thermal_control=dict(
            comparison_gap=M,log_excitation_to_deexcitation_ratio=log_ratio,
            log10_ratio=log_ratio/math.log(10),
            excitation_to_deexcitation_ratio=math.exp(log_ratio),
            charge_variance_formula='Var(Q_mode)=2 n(1+n), n=1/(exp(2pi omega/H)-1)',
            local_charge_distribution='P(Q=k)=(1-z)/(1+z)*z^abs(k), z=exp(-2pi omega/H)',
            nonzero_charge_probability=2*math.exp(log_ratio)/(1+math.exp(log_ratio)),
            log_charge_variance=log_variance,
            single_stationary_mode_charge_rms=math.exp(log_variance/2),
            equal_for_both_charges=True,
            scope='KMS detailed balance at the selected illustrative gap M1. This is neither total abundance nor baryon/photon eta; static frequencies form a continuum.'),
        imbalance=dict(
            invariant_state_charge_asymmetry=0.,eta_if_identified_with_baryons_and_positive_photons=0.,
            target=6.1e-10,conditional_miss_percent=-100.,
            physical_baryon_photon_eta=None,
            missing='An inherited charge-asymmetric preparation/interaction and a baryon/photon identification. The geometric charge pair alone provides neither.',
            fluctuation_scope='Local charge fluctuations need region, mode weights and photon population before a stochastic baryon/photon ratio can be computed. Zero expectation is not identically zero in every local measurement.',
            no_target_used_in_field_or_state=True),
        nonuniqueness=dict(
            construction='Smooth coherent test-field states related by charge conjugation have equal external stress and opposite charge; convex mixtures vary mean charge while preserving that stress.',
            mean_charge='(2p-1)Q_coherent',
            caveat='Fixed-background test-field statement. Full gauge/mixed-metric constraints must be restored for backreacting charged solutions; equal external stress is not identical full higher-dimensional initial data.',
            continuation='Cauchy data and their inherited current, not a coordinate zero, select the nonsymmetric state.'),
        dispositions=[
            'Homogeneous neutral truncation: no internal rotation modes; it cannot identify the offset with a particle charge.',
            'Round-sphere harmonic completion: derives a degenerate complex pair and quantized internal charges without fitting a coupling or mass.',
            'Invariant Euclidean preparation: zero mean charge, not a tiny nonzero phase; neutral horizon transport preserves this result.',
            'Biased or internally rotating preparation: can carry charge, but its amplitude/current and baryon conversion remain new state data.'],
        requested_nonzero_eta_derived=False)
