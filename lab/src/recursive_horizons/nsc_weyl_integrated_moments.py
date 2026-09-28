"""Momentum-integrated defect bounds from certified weighted Fourier sums.

This bounds high canonical-k moments of an auxiliary Weyl defect. It is not
an identification with the missing input-frequency tail of the source sum.
"""
from math import factorial
from flint import arb,ctx
from .nsc_bloch_trace_bound import _positive_interval,_norm_upper
from .nsc_ks_ball_trajectory import exact_upper
from .nsc_scaled_reference_projector import momentum_D


def matrix_nuclear_upper(v):
    frob=sum((a.abs_upper()**2 for a in v),arb(0))
    determinant=(v[0]*v[3]-v[1]*v[2]).abs_upper()
    return (frob+2*determinant).sqrt().upper()


def matrix_fourier_sum_upper(jet,period):
    """A0(matrix symbol)<=sup nuclear + L^2/12 sup nuclear(dz^2)."""
    return (matrix_nuclear_upper(jet.value)+period**2/12*matrix_nuclear_upper(jet.derivative(z=2).value)).upper()


def integrated_defect_moments(scaled_by_sign,global_bounds,inverse_axial,period_length,
                              transfer_catalog,canonical_cutoff,*,low_momentum_auxiliary=False,
                              uniform_axial_coverage=False,bits=160):
    """High-|k| R_rho moment densities before the rho integral and channel factor.

    Optional Q_aux=P0+chi_K sum(P1..P4), with chi_K=0 below K/4 and1 above
    K/2. This only changes the error-analysis split C=Q_aux+D. Physical
    subtraction remains unchanged; Q_aux-P_ref must be retained on reconstruction.
    """
    if uniform_axial_coverage is not True:raise ValueError('whole-cell axial derivative bounds required')
    if set(scaled_by_sign)!= {-1,1}:raise ValueError('both canonical momentum signs required')
    with ctx.workprec(bits):
        K=_positive_interval(canonical_cutoff,'canonical momentum split',bits=bits)
        L=_positive_interval(period_length,'numerical period',bits=bits)
        inva=_positive_interval(inverse_axial,'inverse axial scale',bits=bits)
        if global_bounds.formal_order!=4:raise ValueError('fourth-order global symbol bounds required')
        catalog=transfer_catalog['catalog']
        if transfer_catalog['canonical_cutoff_lower']!=K.lower():raise ValueError('Fourier weights belong to another canonical split')
        def weight(name):return _norm_upper(catalog[name]['total'],'weighted reciprocal Fourier sum',bits=bits)
        rows=[]
        for sign in (-1,1):
            jets=scaled_by_sign[sign]
            if jets.sign!=sign or jets.formal_order!=4 or jets.physical<7 or jets.momentum<5:
                raise ValueError('correctly signed projector jets with retained true-defect derivatives required')
            if jets.history_identity!=global_bounds.history_identity:
                raise ValueError('local and global projector histories differ')
            if not global_bounds.mass.contains(jets.mass) or not global_bounds.angular.contains(jets.angular):
                raise ValueError('global symbol bound does not cover the local channel')
            needed=(4 if low_momentum_auxiliary else 2)/K
            if not jets.mu.contains(arb(0)) or not jets.mu.contains(needed):
                raise ValueError('mu box must cover all shifted or cutoff-supported momenta')
            ell=jets.angular.abs_upper();B4=jets.B[4]
            leading=[]
            for z in (0,2):
                T=B4.derivative(t=1,z=z).value;Z=B4.derivative(z=z+1).value
                value=[T[0]+inva*Z[0],T[1],T[2],T[3]-inva*Z[3]]
                leading.append(matrix_nuclear_upper(value))
            A5=leading[0]+L**2/12*leading[1]
            Bnorm=[matrix_fourier_sum_upper(B,L) for B in jets.B]
            global_norm=[]
            for j in range(5):
                if j==0:
                    # P0 is a rank-one orthogonal projector at real k.
                    bound=1+L**2/6*global_bounds.derivative_bound(0,z=2)
                elif low_momentum_auxiliary:
                    bound=needed**(j+1)*Bnorm[j]
                else:
                    bound=2*global_bounds.derivative_bound(j)+L**2/6*global_bounds.derivative_bound(j,z=2)
                global_norm.append(bound.upper())
            for s in (0,1):
                kinetic=A5/((4-s)*K.lower()**(4-s));small=arb(0);shifted=arb(0);taylor=arb(0)
                for j in range(5):
                    n=5-j
                    core=momentum_D(jets.B[j],j+1,n,sign,jets.mu)
                    core_norm=matrix_fourier_sum_upper(core,L)
                    small+=arb(128)*ell*core_norm*weight(f'small_n{n}_s{s}')/factorial(n)
                    shifted+=2*ell*global_norm[j]*weight(f'shift_s{s}')
                    for n in range(5-j):
                        if j==0 and n==0:
                            taylor+=2*ell*weight(f'shift_s{s}')
                        core=momentum_D(jets.B[j],j+1,n,sign,jets.mu)
                        core_norm=matrix_fourier_sum_upper(core,L)
                        p=j+1+n
                        taylor+=2*ell*core_norm*weight(f'taylor_n{n}_p{p}_s{s}')/factorial(n)
                factor=inva.abs_upper()/(2*arb.pi())
                total=factor*(kinetic+small+shifted+taylor)
                rows.append({'sign_k':sign,'moment_order':s,
                    'kinetic_time_upper':exact_upper(factor*kinetic),
                    'small_transfer_upper':exact_upper(factor*small),
                    'large_shifted_upper':exact_upper(factor*shifted),
                    'large_Taylor_upper':exact_upper(factor*taylor),
                    'total_upper':exact_upper(total)})
        from .nsc_ks_ball_trajectory import restored_upper
        sums=[sum((restored_upper(r['total_upper']) for r in rows if r['moment_order']==s),arb(0)) for s in (0,1)]
        return {'schema':'NSC-WEYL-INTEGRATED-MOMENTS-v1','rows':rows,
            'zeroth_moment_density_upper':exact_upper(sums[0]),'first_moment_density_upper':exact_upper(sums[1]),
            'momentum_measure_included':'dk/(2pi), both signs exactly once',
            'channel_multiplicity_included':False,'rho_integral_included':False,
            'auxiliary':('P0+chi_K*(P1+P2+P3+P4)' if low_momentum_auxiliary else 'P0+P1+P2+P3+P4'),
            'auxiliary_cutoff':({'zero_through':'K/4','one_from':'K/2','smooth_even_between':True} if low_momentum_auxiliary else None),
            'physical_subtraction_changed':False,'auxiliary_reconstruction_correction_required':low_momentum_auxiliary,
            'low_canonical_momenta_need_separate_evaluation':True,
            'initial_covariance_error_bound':None,'source_energy_tail_identified_with_canonical_cut':False,
            'physical_local_gate':'OPEN'}
