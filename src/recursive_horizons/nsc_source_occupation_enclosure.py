"""Vacuum-distance bound for the unchanged three-channel affine source.

The horizon 2x2 block of C-P_vac has eigenvalues +/-sqrt(f). The third
entry is the existing incoming occupation n, including its declared gap.
A normalized sewing and homogeneous unitary transport cannot enlarge its
operator norm. No reflection phase is set to zero or discarded.
"""
from flint import arb,ctx
from .nsc_vacuum_source_remainder import _real


def occupation_vacuum_distance(energy,mass,kappa,omega,*,bits=192):
    with ctx.workprec(bits):
        E,m,k,w=(_real(v,n) for v,n in zip((energy,mass,kappa,omega),
                                         ('energy','mass','kappa','omega')))
        if not all(v.is_finite() for v in (E,m,k,w)) or not E>0 or not m>=0 or not k>0 or not w>0:
            raise ValueError('positive source scales/energy and nonnegative mass required')
        f=1/(1+(2*arb.pi()*E/k).exp())
        gap=bool(m>0 and E<w*m)
        incoming=arb(0) if gap else 1/(1+(2*arb.pi()*E/(w*k)).exp())
        bound=f.sqrt().max(incoming).upper()
        return {'operator_distance_upper':bound,'horizon_occupation_upper':f.upper(),
                'incoming_occupation_upper':incoming.upper(),'incoming_gap_used':gap,
                'coherence_retained':True,'signed_source_norm_equal':True,
                'requires_normalized_sewing':True,'source_occupation_law_changed':False,
                'field_or_quadrature_error_included':False,'physical_local_gate':'OPEN'}
