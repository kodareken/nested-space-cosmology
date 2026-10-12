"""Continuous time/space bound of the polynomial KS residual via convolution.

The Fourier profile table encloses continuum coefficients. Bernstein
convexity covers the entire time cell, and coefficient absolute sums cover
the spatial period. Time-coefficient and radius remainders are added below.
"""
from flint import acb, acb_poly, arb, ctx


def _complex_l1(value):
    return (abs(value.real)+abs(value.imag)).upper()


def _polynomial_l1(poly, mode_origin):
    zeroth, first = arb(0), arb(0)
    for j,value in enumerate(poly.coeffs()):
        upper = _complex_l1(value)
        zeroth += upper
        first += abs(mode_origin+j)*upper
    return zeroth.upper(), first.upper()


def profile_sup_bounds(profiles, length):
    omega = 2*arb.pi()/arb(length)
    result = {}
    for key, profile in profiles.items():
        coefficients, tail = profile['coefficients'], profile['tail']
        K = (len(coefficients)-1)//2
        norms = [arb(0),arb(0)]
        for j,value in enumerate(coefficients):
            bound = _complex_l1(value)
            norms[0] += bound
            norms[1] += abs(j-K)*omega*bound
        result[key] = ((norms[0]+tail[0]).upper(),(norms[1]+tail[1]).upper())
    return result


def polynomial_residual_bounds(residual, profiles):
    """Two row-l2 upper bounds for S and S_z on the whole xi/z cell.

    The finite convolution uses an l1 source-column majorant, which is
    conservative for row l2. Profile tail products use row-l2 moments.
    """
    if set(profiles)!=set(residual.terms):
        raise ValueError('every spatial potential profile is required')
    field=residual.field
    with ctx.workprec(residual.bits):
        sizes={(len(row['coefficients'])-1)//2 for row in profiles.values()}
        if len(sizes)>1:
            raise ValueError('one common retained profile band is required')
        K=next(iter(sizes),0)
        if any(len(row['coefficients'])!=2*K+1 for row in profiles.values()):
            raise ValueError('odd complete Fourier profile bands required')
        order=sorted(range(field.spatial_count),key=lambda j:field.modes[j])
        mode_origin=min(field.modes)-K
        profile_polys={key:acb_poly(row['coefficients']) for key,row in profiles.items()}
        omega=2*arb.pi()/field.length
        moment0={key:residual._moment(data,0) for key,data in residual.terms.items()}
        moment1={key:residual._moment(data,1) for key,data in residual.terms.items()}
        maximum=[arb(0),arb(0)]
        for q in range(residual.degree+1):
            for spin in range(2):
                bound0,bound1=arb(0),arb(0)
                for source in range(field.source_count):
                    free=acb_poly([residual.free[q][spin][source][j] for j in order]).left_shift(K)
                    for key,profile in profile_polys.items():
                        part=acb_poly([residual.terms[key][q][spin][source][j] for j in order])
                        free += profile*part
                    row0,row1=_polynomial_l1(free,mode_origin)
                    bound0 += row0
                    bound1 += omega*row1
                for key,row in profiles.items():
                    t0,t1=map(arb,row['tail'][:2])
                    if not t0>=0 or not t1>=0:raise ValueError('nonnegative profile tails required')
                    bound0 += t0*moment0[key][q][spin]
                    bound1 += t1*moment0[key][q][spin]+t0*moment1[key][q][spin]
                maximum[0]=arb.max(maximum[0],bound0.upper())
                maximum[1]=arb.max(maximum[1],bound1.upper())
        return tuple(v.upper() for v in maximum)


def operator_remainder_bounds(field, time_errors, profile_sup, radius_tail,
                               mass, angular, energies):
    """Continuous normalized residual from time Taylor and radius remainders."""
    with ctx.workprec(field.bits):
        moments=field.derivative_moment_bounds(2)
        X=[]
        for j in range(3):
            X.append(arb.max(*(sum((moments[p][j][spin] for p in range(8)),arb(0)).upper()
                               for spin in range(2))))
        h=abs(arb(field.rho_end)-arb(field.rho_start))
        m,ell=abs(arb(float(mass))),abs(arb(float(angular)))
        E=max(abs(float(e)) for e in energies)
        a2,a,ar=(arb(time_errors[k]) for k in ('inv_a2','inv_a','inv_ar'))
        error0=a2*(X[1]+E*X[0])+(m*a+ell*ar)*X[0]
        error1=a2*(X[2]+E*X[1])+(m*a+ell*ar)*X[1]
        for key,error in time_errors.items():
            if isinstance(key,tuple):
                Z0,Z1=profile_sup[key]
                error0 += error*Z0*X[0]
                error1 += error*(Z1*X[0]+Z0*X[1])
        v0,v1=map(arb,radius_tail[:2])
        error0 += v0*X[0]
        error1 += v1*X[0]+v0*X[1]
        return ((h*error0).upper(),(h*error1).upper())
