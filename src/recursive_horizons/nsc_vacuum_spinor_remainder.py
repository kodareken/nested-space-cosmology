"""Phase-fixed upstream vacuum spinor from the owned Bloch remainder.

This supplies an initial column bound only. It does not propagate A_j, bound
changed-history L0 A_M, or remove the finite-occupation source remainder.
"""
from flint import arb,arb_series
from .nsc_massive_jost_phase_bound import _series_context,_polynomial
from .nsc_vacuum_source_remainder import VacuumSourceExpansion,_real


def spinor_remainder(model,rho,cutoff,*,cells=64):
    """Return degree M-1 spinor coefficients and ||u-p|| <= C/E^M.

    At this one homogeneous upstream slice choose u0 real positive. The
    per-energy phase is held fixed during subsequent spacetime transport.
    The exact vacuum Bloch vector has unit norm by its skew evolution.
    """
    if not isinstance(model,VacuumSourceExpansion):
        raise TypeError('owned vacuum expansion required')
    M=model.order
    with _series_context(model.bits,M):
        E=_real(cutoff,'energy cutoff')
        if not E.is_finite() or not E>0:raise ValueError('positive cutoff required')
        x0=1/E.lower();delta=model.target_distance(rho)
        _,ps,zs=model.jets(delta)
        transverse=[[arb(0)]+[delta.sqrt()*ps[j][k][0] for j in range(1,M+1)] for k in (0,1)]
        longitudinal=[zs[j][0] for j in range(M+1)]
        C=model.remainder_constant(rho,cells=cells)
        def chart(x):
            z=_polynomial(longitudinal,x,M)
            a=((1+z)/2).sqrt()
            return (a,*(_polynomial(w,x,M)/(2*a) for w in transverse))
        origin=chart(arb_series([0,1],prec=M+1))
        domain=arb(0).union(x0)
        zrange=_polynomial(longitudinal,arb_series([domain],prec=M+1),M)[0]
        # The segment from n_M to the exact n stays in this north chart.
        min_square=((1+zrange.lower()-C*x0**M)/2).lower()
        if not min_square>0:raise ArithmeticError('vacuum north chart not separated from its pole')
        amin=min_square.sqrt().lower()
        wmax=sum((_polynomial(w,arb_series([domain],prec=M+1),M)[0].abs_upper()**2
                  for w in transverse),arb(0)).sqrt().upper()
        # Jacobian bound: |da|<=|dn|/(4a),
        # |dv|<=|dn|/(2a)+|w_M||dn|/(8a^3), using an exact denominator identity.
        lipschitz=(3/(4*amin)+wmax/(8*amin**3)).upper()
        whole=chart(arb_series([domain,1],prec=M+1))
        taylor=sum((v[M].abs_upper()**2 for v in whole),arb(0)).sqrt().upper()
        total=(C*lipschitz+taylor).upper()
        if not total.is_finite():raise ArithmeticError('nonfinite spinor remainder')
        return {'coefficients':tuple(tuple(v[j] for j in range(M)) for v in origin),
                'order':M,'cutoff':E,'bloch_remainder_constant':C,
                'chart_major_lower':amin,'conversion_lipschitz_upper':lipschitz,
                'chart_taylor_constant':taylor,'spinor_remainder_constant':total,
                'upstream_envelope_z_derivatives_zero':True,
                'finite_occupation_error_included':False,
                'changed_history_remainder':None,'physical_local_gate':'OPEN'}


def truncated_upstream_envelope(result,degree,period_length):
    """Initial H2 triple for the SAME phase-fixed polynomial through degree.

    Norm convention is (||e||L2, ||e_z||L2, ||e_zz||L2) for the envelope
    on a finite period. This does not bound derivatives of the physical carrier.
    The caller must use the returned coefficients as its upstream recurrence data.
    """
    M=result['order']
    if type(degree) is not int or not 0<=degree<M:
        raise ValueError('truncation degree must be below the proved remainder order')
    with _series_context(192,M):
        L=_real(period_length,'period length')
        if not L.is_finite() or not L>0:raise ValueError('positive finite period required')
        E=result['cutoff'].lower();coeff=result['coefficients']
        tail=result['spinor_remainder_constant']/E**(M-degree)
        for j in range(degree+1,M):
            tail+=sum((row[j].abs_upper()**2 for row in coeff),arb(0)).sqrt()/E**(j-degree)
        tail=tail.upper()
        return {'coefficients':tuple(row[:degree+1] for row in coeff),
                'inverse_energy_order':degree,'pointwise_scaled_error_upper':tail,
                'scaled_envelope_H2_triple':((L.sqrt()*tail).upper(),arb(0),arb(0)),
                'norm_domain_length':L,'carrier_derivatives_included':False,
                'requires_same_upstream_coefficients':True,
                'changed_history_remainder':None,'physical_local_gate':'OPEN'}
