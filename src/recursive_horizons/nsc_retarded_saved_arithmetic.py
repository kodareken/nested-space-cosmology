"""Finite saved-array diagnostic checks, with no physical error budget.

Every allowance uses the uncancelled operands of the stated finite operation.
Linear-system checks test backward residuals rather than equality of two host
LAPACK solutions. The spatial checks separate binary64 operation allowances
from the configured chart integration indicators; the latter are unproved
accuracy settings. Neither component closes an evolution or matching gate.
"""
import numpy as np

from .nsc_ks_saved_value_arithmetic import _Ball, _window, _gamma, _up


def _finite(value):
    a = np.asarray(value)
    if a.dtype.kind not in 'fciu' or not np.isfinite(a).all():
        raise ValueError('finite numeric operands required')
    return a


def require_linear_system(matrix, rhs, current, saved, *, label):
    """Both solutions must satisfy the same sealed square system.

    A complex length-n dot uses at most eight real operations per term;
    eight additional operations cover the final subtraction and modulus.
    The residual cap is relative to |A||X|+|B|, never a cancellation result.
    This is a backward-residual check, not a forward/physical error bound.
    """
    A, B, X, Y = map(_finite, (matrix, rhs, current, saved))
    n = A.shape[-1]
    if A.shape[-2] != n or X.shape != Y.shape or (A @ X).shape != B.shape:
        raise ValueError(label + ': linear-system layout changed')
    results = []
    for solution in (X, Y):
        residual = abs(A @ solution - B)
        cap = _up(_gamma(8*n+8)*(abs(A) @ abs(solution)+abs(B)))
        if np.any(residual > cap):
            raise ValueError(label + ': sealed linear-system residual failed')
        results.append(float(np.max(residual)))
    return {'scope': 'finite backward residual only', 'residuals': results,
            'physical_error_bound': None}


def require_matrix_product(current, saved, left, right, *, label):
    A, B, X, Y = map(_finite, (left, right, current, saved))
    if X.shape != Y.shape or (A @ B).shape != X.shape:
        raise ValueError(label + ': product layout changed')
    # Charge both finite contractions and their complex modulus comparison.
    cap = _up(2*_gamma(8*A.shape[-1]+8)*(abs(A) @ abs(B)))
    if np.any(abs(X-Y) > cap):
        raise ValueError(label + ': finite product replay failed')
    return float(np.max(cap))


def require_singular_invariants(matrix, values):
    """Check the stored three singular values against Gram invariants.

    This checks numerical polynomial residuals, rather than declaring a
    LAPACK SVD error bound. The operation count covers a 3x3 Gram product,
    its six-term determinant, and the three elementary symmetric sums.
    """
    A = _finite(matrix)
    s = _finite(values)
    if (A.shape != (3, 3) or s.shape != (3,) or np.iscomplexobj(A) or np.iscomplexobj(s)
            or np.any(s <= 0) or np.any(np.diff(s) > 0)):
        raise ValueError('selected singular-value layout changed')
    G = A.T @ A
    lam = s*s
    determinant_terms = np.array([G[0,0]*G[1,1]*G[2,2], G[0,1]*G[1,2]*G[2,0],
        G[0,2]*G[1,0]*G[2,1], -G[0,2]*G[1,1]*G[2,0],
        -G[0,1]*G[1,0]*G[2,2], -G[0,0]*G[1,2]*G[2,1]])
    expected = np.array([np.trace(G),
        G[0,0]*G[1,1]+G[0,0]*G[2,2]+G[1,1]*G[2,2]
        -G[0,1]*G[1,0]-G[0,2]*G[2,0]-G[1,2]*G[2,1],
        determinant_terms.sum()])
    actual = np.array([lam.sum(), lam[0]*lam[1]+lam[0]*lam[2]+lam[1]*lam[2], np.prod(lam)])
    mass = np.array([np.sum(abs(A)**2)+lam.sum(),
        np.sum(abs(G))**2+lam.sum()**2,
        np.sum(abs(determinant_terms))+np.prod(lam)])
    # 64*n^3 real operations is a conservative count for the independently
    # evaluated fixed 3x3 invariants, including Gram-operand accumulation.
    if np.any(abs(expected-actual) > _up(_gamma(64*3**3)*mass)):
        raise ValueError('selected singular-value invariant residual failed')


def require_selected_fit(current, saved, current_coefficients, saved_coefficients):
    if current.keys() != saved.keys():
        raise ValueError('selected calibration metadata keys changed')
    for key in current:
        if key not in ('singular_values', 'linear_system_residual') and current[key] != saved[key]:
            raise ValueError('selected calibration metadata changed: '+key)
    A = np.array(saved['matrix']); b = np.array(saved['target'])
    require_linear_system(A, b, current_coefficients, saved_coefficients, label='selected real calibration')
    require_singular_invariants(A, current['singular_values'])
    require_singular_invariants(A, saved['singular_values'])
    for fit, c in ((current, current_coefficients), (saved, saved_coefficients)):
        observed = float(np.max(abs(A @ c-b)))
        cap = float(np.max(_up(_gamma(8*3+8)*(abs(A) @ abs(c)+abs(b)))))
        if not np.isfinite(fit['linear_system_residual']) or abs(fit['linear_system_residual']-observed) > cap:
            raise ValueError('selected calibration residual diagnostic changed')


def require_diagnostic_receipt(actual, saved, allowances):
    """Only enumerated floating diagnostics may differ; typed metadata exact."""
    differences = []
    def visit(a, b, path):
        if type(a) is not type(b):
            raise ValueError('diagnostic receipt type changed: '+str(path))
        if isinstance(a, dict):
            if a.keys() != b.keys():
                raise ValueError('diagnostic receipt keys changed: '+str(path))
            for key in a: visit(a[key], b[key], path+(key,))
        elif isinstance(a, list):
            if len(a) != len(b): raise ValueError('diagnostic receipt length changed')
            for index, (x, y) in enumerate(zip(a, b)): visit(x, y, path+(index,))
        elif a != b:
            cap = allowances.get(path)
            if (type(a) is not float or cap is None or not np.isfinite(a)
                    or not np.isfinite(b) or not np.isfinite(cap) or cap < 0 or abs(a-b) > cap):
                raise ValueError('diagnostic receipt changed beyond finite operation policy: '+str(path))
            differences.append({'path': list(path), 'abs_difference': abs(a-b), 'allowance': float(cap)})
    visit(actual, saved, ())
    return differences


def _sqrt(value):
    lower = value.v-value.e
    if np.any(lower <= 0): raise ValueError('diagnostic square-root indicator reaches zero')
    root = np.sqrt(value.v)
    return _Ball(root, _up(value.e/(root+np.sqrt(lower))+4*np.finfo(float).eps*root))


def _norm(value, error=0.):
    a = _finite(value)
    e = np.broadcast_to(error, a.shape)
    mass = abs(a.real)+abs(a.imag)
    # Complex modulus/square: hypot plus multiplication, eight operations.
    square = abs(a)**2
    radius = 2*abs(a)*e+e*e+_gamma(8)*mass*mass
    total = _Ball(np.sum(square), _up(np.sum(radius)+_gamma(a.size)*np.sum(square)))
    return _sqrt(total)


def _difference(a, b):
    a, b = _finite(a), _finite(b)
    return a-b, _up(_gamma(2)*(abs(a)+abs(b)))


def localization_allowances(arrays):
    p = ('localization',)
    result = {}
    n = len(arrays['401_64/x'])
    coarse = arrays['401_64/field_tangent'].reshape(2,n,12)
    fine = arrays['801_64/field_tangent'].reshape(2,801,12)[:,::2]
    difference, de = _difference(coarse,fine)
    dt, dte = _difference(arrays['401_64/trace_tangent'], arrays['801_64/trace_tangent'])
    cd, cde = _difference(arrays['401_64/field'],arrays['401_64/stationary_field'])
    fd, fde = _difference(arrays['801_64/field'],arrays['801_64/stationary_field'])
    bd, bde = _difference(cd.reshape(2,n,12),fd.reshape(2,801,12)[:,::2])
    bde = _up(bde+cde.reshape(2,n,12)+fde.reshape(2,801,12)[:,::2])
    ct, cte = _difference(arrays['401_64/trace'],arrays['401_64/stationary_trace'])
    ft, fte = _difference(arrays['801_64/trace'],arrays['801_64/stationary_trace'])
    btd, btde = _difference(ct,ft); btde = _up(btde+cte+fte)
    for name,a,ae,b,be in [('field',difference,de,bd,bde),('trace',dt,dte,btd,btde)]:
        numerator = float(abs(np.vdot(a,b)))
        radius = _up(np.sum(abs(a)*be+abs(b)*ae+ae*be)+_gamma(8*a.size+8)*np.sum(abs(a)*abs(b)))
        correlation = _Ball(numerator,radius)/(_norm(a,ae)*_norm(b,be))
        result[p+('baseline_discriminator',name+'_discrepancy_normalized_correlation')] = float(_up(2*correlation.e))
    # Every recorded stationary-drift maximum uses the same finite subtraction
    # and complex-modulus policy. Charge the original operands before their
    # cancellation, and both evaluations; metadata and locations stay exact.
    for case in ('401_64', '801_64', '801_128'):
        for component in ('field', 'trace'):
            drift, error = _difference(arrays[case+'/'+component],
                                       arrays[case+'/stationary_'+component])
            result[p+('per_case',case,'stationary_'+component+'_drift','maximum')] = float(
                _up(2*np.max(error+_gamma(8)*abs(drift))))
    kc, kf = arrays['401_64/covariance_tangent'],arrays['801_64/covariance_tangent']
    delta, error = _difference(kc,kf)
    ratio = _norm(delta,error)/_norm(kf)
    result[p+('spatial_covariance','Frobenius_relative')] = float(_up(2*ratio.e))
    roi = np.flatnonzero((arrays['401_64/x']>.8)&(arrays['401_64/x']<1.12))
    delta, error = _difference(coarse[0,roi],fine[0,roi])
    ratio = _norm(delta,error)/_norm(fine[0,roi])
    result[p+('translation_discriminator','unshifted','relative_shape_difference')] = float(_up(2*ratio.e))
    return result


def _bump(z, center, halfwidth):
    u = (z-center)/halfwidth
    active = abs(u.v)<1
    result = _Ball(np.zeros_like(u.v))
    if np.any(abs(abs(u.v)-1) <= u.e):
        raise ValueError('diagnostic axial cutoff indicator crosses a branch')
    if not np.any(active): return result
    v = _Ball(u.v[active],u.e[active])
    exponent = 1-1/(1-v*v)
    value = np.exp(exponent.v)
    error = _up(np.exp(exponent.v+exponent.e)*exponent.e+4*np.finfo(float).eps*value)
    result.v[active],result.e[active] = value,error
    return result


def local_forcing_allowances(x801,phi801,data,support,pilot,*,indicators):
    """Charge the unchanged radius-only sparse two-matvec/DFT diagnostics.

    The finite DFT uses its defining n-term sum with complex products and
    elementary twiddle operations. This is an accumulation policy, not a
    proof of a particular FFT library implementation. Chart epsabs/epsrel
    settings are accounted separately as unproved coordinate indicators.
    """
    from .nsc_transmitting_history_jets import generator_direction
    allowances = {}
    direction = pilot.CompatibleRadiusDirection(pilot.axial_bump(support['axial_center']),lambda z,n:0.,.007,.03)
    grids = {}
    for n in (401,801):
        stride = 800//(n-1);x=x801[::stride]
        phi=phi801.reshape(2,801,12)[:,::stride].reshape(2*n,12)
        owner=pilot.FourthOrderModePropagator(x,data['mass'].item(),data['angular'].item())
        provider=pilot.CachedZeroRadiusProvider(x,direction)
        grids[n]=(x,phi,owner,provider)
    for index,t in enumerate((.11,.15,.19)):
        samples={}
        for n,(x,phi,owner,provider) in grids.items():
            metric=provider.values(t,x);d=provider.log_directions(t,x,metric)[0]
            coordinates=np.array([pilot.chart_coordinates(float(rho)) for rho in x])
            anchor=pilot.chart_coordinates(1.)[0]
            ne=_gamma(2)*(abs(coordinates[:,0])+abs(anchor))
            se=_gamma(2)*abs(coordinates[:,1])
            if indicators:
                ne+=2e-13*(np.maximum(1,abs(coordinates[:,0]))+max(1,abs(anchor)))
                se+=2e-13*np.maximum(1,abs(coordinates[:,1]))
            normal=_Ball(provider.normal,_up(ne));space=_Ball(t+provider.shift,_up(se+_gamma(2)*(abs(t)+abs(provider.shift))))
            radius=_sqrt(1+_Ball(x)**2)
            raw=_window(normal,.007,.03,0)[0]*normal*_bump(space,support['axial_center'],.06)
            dlog=raw/radius;basis=dlog/radius
            if np.any(abs(dlog.v-d[3])>dlog.e):
                raise ValueError('current radius direction differs from finite operation model')
            L,_=owner.generator(metric);dL,_=generator_direction(owner,metric,d)
            D=abs(owner.D)
            beta=metric[1]
            # The geometry uses arctan/sqrt; include four ulps each and
            # uncancelled beta_squared arithmetic, propagated through sqrt.
            angle=np.pi/2-np.arctan(x)
            angle_error=_gamma(8)*(abs(np.pi/2)+abs(np.arctan(x)))
            squared_error=_gamma(8)*3*((1+x*x)*abs(angle)+abs(x))+3*(1+x*x)*angle_error
            beta_error=_up(squared_error/(2*beta)+4*np.finfo(float).eps*abs(beta))
            speed_error=beta_error+_gamma(2)*(1+abs(beta))
            from scipy.sparse import diags,block_diag,bmat
            transport_error=.5*(diags(speed_error)@D+D@diags(speed_error))
            transport_error+=_gamma(16)*abs(L[:n,:n])
            potential_error=_up(abs(owner.angular)*(radius.e/radius.v**2)+_gamma(8)*(abs(owner.mass)+abs(owner.angular)/radius.v))
            Le=block_diag([transport_error,transport_error],format='csr')+bmat([[None,diags(potential_error)],[diags(potential_error),None]],format='csr')
            de=_up(abs(owner.angular)*(dlog.e/radius.v+abs(dlog.v)*radius.e/radius.v**2)+_gamma(8)*abs(owner.angular*d[3]/metric[3]))
            dLe=bmat([[None,diags(de)],[diags(de),None]],format='csr')
            phase=np.exp(-1j*np.repeat(data['energies'],3)*t)
            phased=phi*phase[None,:]
            pe=_gamma(16)*abs(phi)
            F=dL@phased
            Fe=_up(abs(dL)@pe+dLe@(abs(phased)+pe)+_gamma(8*max(np.diff(dL.indptr))+8)*(abs(dL)@abs(phased)))
            LF=L@F
            LFe=_up(abs(L)@Fe+Le@(abs(F)+Fe)+_gamma(8*max(np.diff(L.indptr))+8)*(abs(L)@abs(F)))
            samples[n]=(F,Fe,LF,LFe)
            p=('local_discriminator','times',index,'grids',str(n))
            allowances[p+('F_maximum',)]=float(_up(2*np.max(Fe+_gamma(8)*abs(F))))
            allowances[p+('LF_maximum',)]=float(_up(2*np.max(LFe+_gamma(8)*abs(LF))))
            spectrum=np.fft.fft(d[3]/metric[3])
            fft_error=_up(np.sum(basis.e)+_gamma(16*n)*np.sum(abs(basis.v)+basis.e))
            power=abs(spectrum)**2
            power_error=_up(2*abs(spectrum)*fft_error+fft_error**2+_gamma(8)*power)
            total=_Ball(np.sum(power),_up(np.sum(power_error)+_gamma(n)*np.sum(power)))
            theta=2*np.pi*np.fft.fftfreq(n);group=(4*np.cos(theta)-np.cos(2*theta))/3
            group_error=_gamma(32)*(4+1)/3
            if np.any(abs(abs(group-1)-.25)<=group_error) or np.any(abs(group)<=group_error):
                raise ValueError('DFT group mask differs within its operation allowance')
            for key,mask in [('geometry_power_fraction_group_error_above_0_25',abs(group-1)>.25),('geometry_power_fraction_reversed_group_factor',group<0)]:
                numerator=_Ball(np.sum(power[mask]),_up(np.sum(power_error[mask])+_gamma(int(np.sum(mask)))*np.sum(power[mask])))
                ratio=numerator/total
                allowances[p+(key,)]=float(_up(2*ratio.e))
        common=lambda a:a.reshape(2,801,12)[:,::2].reshape(802,12)
        _,_,lfc,lfe=samples[401];_,_,lff,lffe=samples[801]
        lff,lffe=common(lff),common(lffe)
        delta=lfc-lff;error=_up(lfe+lffe+_gamma(2)*(abs(lfc)+abs(lff)))
        numerator=_Ball(float(np.max(abs(delta))),float(np.max(error+_gamma(8)*abs(delta))))
        denominator=_Ball(float(np.max(abs(lff))),float(np.max(lffe+_gamma(8)*abs(lff))))
        p=('local_discriminator','times',index)
        allowances[p+('LF_common_node_difference',)]=float(_up(2*numerator.e))
        allowances[p+('LF_fine_common_norm',)]=float(_up(2*denominator.e))
        allowances[p+('LF_relative_difference',)]=float(_up(2*(numerator/denominator).e))
    return allowances


def spatial_diagnostic_policy(arrays,x,phi,data,support,pilot):
    arithmetic={**localization_allowances(arrays),**local_forcing_allowances(x,phi,data,support,pilot,indicators=False)}
    combined={**localization_allowances(arrays),**local_forcing_allowances(x,phi,data,support,pilot,indicators=True)}
    return {'allowances':combined,'arithmetic':arithmetic,
        'configured_chart_indicators':{'epsabs':2e-13,'epsrel':2e-13,'proved_error_bound':False},
        'scope':'saved finite diagnostic comparison only; no field/radial propagation',
        'physical_error_bound':None,'continuum_error_bound':None,'matching_status':'OPEN'}
