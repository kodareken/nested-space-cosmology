"""Quantitative small-transfer coefficients of the actual fourth-order defect.

The bound is pointwise on one geometry box and all canonical |k|>=K.
Large transfer, full box coverage and time integration remain separate.
"""
from math import comb,factorial
from flint import arb,ctx
from .nsc_ks_ball_trajectory import exact_upper,restored_upper
from .nsc_bloch_trace_bound import _norm_upper,_positive_interval
from .nsc_scaled_reference_projector import ScaledProjectorJets,momentum_D


def frobenius_upper(matrix):
    return sum((entry.abs_upper()**2 for entry in matrix),arb(0)).sqrt().upper()


def small_transfer_coefficients(jets,inverse_axial,reciprocal_fourier_moments,
                                 canonical_cutoff,*,axial_order=2,bits=160):
    """||d_z^q R_T,small||_F <= C5_q/|k|^5+C6_q/|k|^6, q<=axial_order.

    R_T,small includes P4_T plus the kinetic first-star term and potential
    transfers |omega|<=|k|. No input source-energy cutoff is asserted here.
    """
    if not isinstance(jets,ScaledProjectorJets) or jets.formal_order!=4:
        raise TypeError('owned fourth-order scaled projector jets required')
    if isinstance(axial_order,bool) or int(axial_order)!=axial_order or not 0<=axial_order<=4:
        raise ValueError('axial derivative order in0..4 required')
    if jets.physical<5+axial_order or jets.momentum<5:
        raise ValueError('projector jets do not retain the required true-defect derivatives')
    if len(reciprocal_fourier_moments)<6+axial_order:
        raise ValueError('reciprocal Fourier moments through5+axial_order required')
    with ctx.workprec(max(bits,jets.bits)):
        cutoff=_positive_interval(canonical_cutoff,'canonical momentum split',bits=bits)
        inv_a=_positive_interval(inverse_axial,'inverse axial scale on geometry box',bits=bits)
        if not jets.mu.contains(arb(0)) or not jets.mu.contains(2/cutoff):
            raise ValueError('mu box must cover[0,2/K] for both shifted momentum segments')
        moments=[_norm_upper(v,'reciprocal Fourier moment',bits=bits) for v in reciprocal_fourier_moments]
        ell=jets.angular.abs_upper()
        B4=jets.B[4];rows=[]
        for q in range(axial_order+1):
            temporal=B4.derivative(t=1,z=q).value
            spatial=B4.derivative(z=q+1).value
            # (1/2){sigma3/a,B4_z} keeps exactly its diagonal entries.
            leading=[temporal[0]+inv_a*spatial[0],temporal[1],temporal[2],temporal[3]-inv_a*spatial[3]]
            C5=frobenius_upper(leading);C6=arb(0)
            for j in range(5):
                n=5-j
                core=momentum_D(jets.B[j],j+1,n,jets.sign,jets.mu)
                for h in range(q+1):
                    derivative=frobenius_upper(core.derivative(z=q-h).value)
                    # Two shifted segments and |k+/-omega/2|>=|k|/2 give2*2^6.
                    C6+=arb(128*comb(q,h))*ell*moments[n+h]*derivative/(2**n*factorial(n))
            rows.append({'axial_order':q,'C5_T_upper':exact_upper(C5),'C6_T_upper':exact_upper(C6),
                'C5_rho_upper':exact_upper(inv_a.abs_upper()*C5),
                'C6_rho_upper':exact_upper(inv_a.abs_upper()*C6)})
        return {'schema':'NSC-WEYL-SMALL-TRANSFER-COEFFICIENTS-v1','rows':rows,
            'canonical_cutoff_lower':exact_upper(cutoff.lower()),'sign_k':jets.sign,
            'history_identity':jets.history_identity,'region':'|k|>=K, |omega|<=|k|',
            'clock_conversion':'R_rho=-R_T/a; a is z-independent',
            'norm':'pointwise2x2 Frobenius on the supplied geometry box',
            'large_transfer_bound':None,'whole_history_coverage':False,'time_integral_bound':None,
            'input_energy_tail_identified_with_canonical_k':False,'physical_local_gate':'OPEN'}


def evaluate_coefficient_row(row,absolute_momentum,*,clock='rho',bits=160):
    if clock not in ('rho','T'):raise ValueError('owned clock rho or T required')
    with ctx.workprec(bits):
        k=_positive_interval(absolute_momentum,'absolute canonical momentum',bits=bits).lower()
        return exact_upper(restored_upper(row['C5_'+clock+'_upper'])/k**5+
                           restored_upper(row['C6_'+clock+'_upper'])/k**6)
