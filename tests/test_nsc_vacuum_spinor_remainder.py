import numpy as np
import pytest
from flint import arb,ctx
from scipy.integrate import solve_ivp
from recursive_horizons.nsc_vacuum_source_remainder import VacuumSourceExpansion
from recursive_horizons.nsc_vacuum_source_correction import VacuumCorrection
from recursive_horizons.nsc_vacuum_spinor_remainder import spinor_remainder

HINT=1.9006916054701435
UP=1.0300000000000002

@pytest.fixture(scope='module')
def result():
    with ctx.workprec(192):
        m=VacuumSourceExpansion(HINT,arb.pi()/2,arb(5).sqrt(),order=8)
        return m,spinor_remainder(m,UP,32,cells=64)


def test_initial_coefficients_and_normalization(result):
    m,r=result
    with ctx.workprec(192):
        a,re,im=r['coefficients'];delta=m.target_distance(UP)
        _,ps,zs=m.jets(delta)
        assert a[0].contains(1) and a[1].contains(0)
        assert re[0].contains(0) and im[0].contains(0)
        assert (re[1]-delta.sqrt()*ps[1][0][0]/2).contains(0)
        assert (im[1]-delta.sqrt()*ps[1][1][0]/2).contains(0)
        assert (a[2]+(re[1]**2+im[1]**2)/2).contains(0)
        for j in range(1,r['order']):
            assert sum((v[k]*v[j-k] for v in (a,re,im) for k in range(j+1)),arb(0)).contains(0)


def test_independent_bloch_ode_and_phase_fixed_reconstruction(result):
    m,r=result
    with ctx.workprec(192):
        end=float(m.target_distance(UP).log().mid())
        for E in (32.,48.):
            flow=VacuumCorrection(m,E)
            n0=np.array([float(v.mid()) for v in flow.expansion_at_log_distance(-20.)])
            # Full equation, not the driven correction or spinor polynomial.
            from recursive_horizons.nsc_subgap_source_covariance import BlochSource
            run=solve_ivp(lambda y,n:BlochSource.rhs_numeric(flow,y,n),(-20.,end),n0,
                          method='DOP853',rtol=2e-13,atol=2e-15,max_step=.005)
            assert run.success
            nx,ny,nz=run.y[:,-1];u0=np.sqrt((1+nz)/2)
            numerical=np.array([u0,(nx+1j*ny)/(2*u0)])
            coeff=r['coefficients'];v=[sum(float(c.mid())/E**j for j,c in enumerate(row)) for row in coeff]
            polynomial=np.array([v[0],v[1]+1j*v[2]])
            assert np.linalg.norm(numerical-polynomial)<float(r['spinor_remainder_constant']/arb(E)**r['order'])
        assert r['changed_history_remainder'] is None
        assert r['finite_occupation_error_included'] is False


def test_chart_failure_and_context(result):
    m,r=result;before=ctx.prec,ctx.cap
    with pytest.raises(ValueError):spinor_remainder(m,UP,None)
    with pytest.raises(ArithmeticError,match='chart'):spinor_remainder(m,UP,.1,cells=2)
    assert (ctx.prec,ctx.cap)==before


def test_signed_vacuum_partner_requires_the_complement():
    # Negative occupied vacuum is the complemented signed covariance, not
    # sigma3-conjugation of the positive occupied column by itself.
    a=np.sqrt(.7);v=np.sqrt(.3)*np.exp(.6j)
    positive=np.array([a,v]);negative=np.array([v,a]);s3=np.diag([1.,-1.])
    Q=np.outer(positive,positive.conj())
    target=np.eye(2)-s3@Q.conj()@s3
    assert np.max(abs(np.outer(negative,negative.conj())-target))<1e-15
    wrong=s3@positive.conj()
    assert np.linalg.norm(np.outer(wrong,wrong.conj())-target)>1
    # The swap preserves the approximation-error norm, including complex data.
    error=np.array([1e-5+2e-5j,-3e-5+4e-5j])
    assert np.linalg.norm(error)==np.linalg.norm(error[::-1])


def test_same_polynomial_truncation_and_envelope_norm(result):
    from recursive_horizons.nsc_vacuum_spinor_remainder import truncated_upstream_envelope
    m,r=result
    with ctx.workprec(192):
        reduced=truncated_upstream_envelope(r,4,arb(3))
        assert reduced['coefficients']==tuple(row[:5] for row in r['coefficients'])
        C=reduced['pointwise_scaled_error_upper']
        for E in (32,160,320):
            omitted=[sum((row[j]/arb(E)**j for j in range(5,r['order'])),arb(0)) for row in r['coefficients']]
            norm=sum((v.abs_upper()**2 for v in omitted),arb(0)).sqrt()
            assert norm+r['spinor_remainder_constant']/arb(E)**r['order']<=C/arb(E)**4
        h0,h1,h2=reduced['scaled_envelope_H2_triple']
        assert h0>=arb(3).sqrt()*C and h1==h2==0
        assert reduced['carrier_derivatives_included'] is False
        assert reduced['changed_history_remainder'] is None
        with pytest.raises(ValueError):truncated_upstream_envelope(r,8,1)
        with pytest.raises(ValueError):truncated_upstream_envelope(r,4,None)


def test_original_channel_record_keeps_transport_and_source_gaps():
    import runpy
    from pathlib import Path
    driver=runpy.run_path(str(Path(__file__).resolve().parents[1]/'scripts/derive_nsc_vacuum_spinor_remainder.py'))
    result=driver['calculate']()
    assert len(result['cases'])==8
    assert {(c['group'],c['angular_sign'],c['cutoff']) for c in result['cases']}=={
        (g,s,e) for g in (1,14) for s in (-1,1) for e in (160,320)}
    assert all(c['spatial_norm_domain_bound'] is None for c in result['cases'])
    assert result['changed_history_C_M'] is None
    assert result['finite_occupation_source_error_included'] is False
    assert result['physical_local_gate']=='OPEN'


def test_upstream_minor_coefficients_match_the_dirac_recurrence(result):
    # Independent metric formulas and the L0 recurrence, rather than reading
    # the same Bloch coefficients back into their own normalization check.
    m,result=result
    with ctx.workprec(192):
        rho=arb(UP);radius=(1+rho*rho).sqrt()
        axial=(3*((1+rho*rho)*(1/rho).atan()-rho)-1).sqrt()
        axial_prime=3*(rho*(1/rho).atan()-1)/axial
        _,real,imag=result['coefficients']
        first_real=axial*m.mass/2;first_imag=-axial*m.angular/(2*radius)
        second_real=axial**2*m.angular*(-axial_prime/radius+axial*rho/radius**3)/4
        second_imag=-axial**2*axial_prime*m.mass/4
        for got,expected in ((real[1],first_real),(imag[1],first_imag),
                             (real[2],second_real),(imag[2],second_imag)):
            assert (got-expected).contains(0)
        # Reversing the Dirac i-sign must not satisfy the same relation.
        assert not (imag[2]+second_imag).contains(0)
