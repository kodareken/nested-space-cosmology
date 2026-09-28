"""Covariance error transport in absolute Fourier coefficients and k moments.

The Fourier-diagonal Hermitian Dirac evolution is treated exactly as a
left/right unitary on each symbol coefficient. Only the spatially varying
potential enters the growth bound. No physical state or action is replaced.
"""
from flint import arb,ctx
from .nsc_bloch_trace_bound import _norm_upper,_positive_interval,local_stress_from_traces
from .nsc_ks_ball_trajectory import exact_upper,restored_upper


def propagate_covariance_moments(initial_zero,initial_first,residual_zero,residual_first,
                                 potential_zero_integral,potential_first_integral,*,bits=160):
    """M0 <= e^(2K0) A0; M1 <= e^(2K0)(A1+K1 A0).

    Mj=sum_l integral |k|^j ||D_l(k)||_1 dk/(2pi), j=0,1.
    Kj=integral sum_l |omega_l|^j ||delta V_l||_op |d rho|.
    A_j=initial moment + time-integrated true-defect moment. None is unknown.
    """
    with ctx.workprec(bits):
        i0,i1,r0,r1,v0,v1=[_norm_upper(v,n,bits=bits) for v,n in (
            (initial_zero,'initial D_up zeroth moment'),(initial_first,'initial D_up first moment'),
            (residual_zero,'time-integrated R zeroth moment'),(residual_first,'time-integrated R first moment'),
            (potential_zero_integral,'potential zeroth integral'),(potential_first_integral,'potential first integral'))]
        a0,a1=i0+r0,i1+r1;growth=(2*v0).exp()
        m0=growth*a0;m1=growth*(a1+v1*a0)
        return {'zeroth_moment_upper':exact_upper(m0),'first_moment_upper':exact_upper(m1),
            'growth_upper':exact_upper(growth),'moment_shift_coefficient':exact_upper(v1),
            'matrix_norm':'nuclear per2x2 Fourier coefficient; potential operator norm',
            'momentum_weights':('1','abs(k)'),'momentum_measure':'dk/(2pi)',
            'half_transfer_factor_retained':True,'Fourier_diagonal_evolution_unitary':True,
            'initial_D_up_retained':True,'time_integrated_R_retained':True,
            'physical_local_gate':'OPEN','source_law_changed':False,'new_action_term':False}


def pure_radius_potential_integrals(q_fourier_zero,q_fourier_first,*,absolute_angular,
                                    axial_lower,rho_up,rho_sigma=1,bits=160):
    """Uniform q=1/r_g-1/r_ref Fourier moments -> true potential integral bounds."""
    with ctx.workprec(bits):
        q0=_norm_upper(q_fourier_zero,'reciprocal zeroth Fourier moment',bits=bits)
        q1=_norm_upper(q_fourier_first,'reciprocal first Fourier moment',bits=bits)
        ell=_norm_upper(absolute_angular,'absolute angular',bits=bits)
        a=_positive_interval(axial_lower,'axial lower',bits=bits)
        up=_positive_interval(rho_up,'upstream coordinate',bits=bits)
        end=_positive_interval(rho_sigma,'incoming coordinate',bits=bits)
        if not up>=end:raise ValueError('decreasing preparation coordinate required')
        factor=(up-end)*ell/a
        return {'potential_zero_integral_upper':exact_upper(factor*q0),
            'potential_first_integral_upper':exact_upper(factor*q1),
            'homogeneous_mass_counted_as_growth':False,'potential_moment_measure':'sum of cell Fourier coefficients',
            'clock':'absolute d rho; H_rho potential=(-m S1+ell/r S2)/a',
            'physical_local_gate':'OPEN'}


def covariance_matter_error(moment_record,*,mass,absolute_angular,axial_lower,radius_lower,
                            multiplicity,continuum_embedding_established,bits=160):
    """Uniform coincidence/current bounds from M0/M1; no additional k powers."""
    if continuum_embedding_established is not True:
        raise ValueError('continuum local embedding must be established separately')
    with ctx.workprec(bits):
        stress=local_stress_from_traces(
            moment_record['zeroth_moment_upper'],moment_record['first_moment_upper'],mass=mass,
            absolute_angular=absolute_angular,axial_lower=axial_lower,radius_lower=radius_lower,
            multiplicity=multiplicity,bits=bits)
        return {**stress,'uniform_on_local_interval':True,'physical_local_gate':'OPEN',
            'covariance_coincidence_uses_zeroth_moment':True,'current_uses_first_moment':True,
            'numerical_operator_error_included':False,'source_accuracy_included':False}
