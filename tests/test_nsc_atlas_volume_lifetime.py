import math
import mpmath as mp
import pytest
import sympy as s
from recursive_horizons.nsc_atlas_parent_scale_selector import VolumeParent
from recursive_horizons.nsc_atlas_volume_lifetime import calculation,pole_shot


def test_completed_square_exclusion_identity():
    q,p,V,Vp,Vpp=s.symbols('q p V Vp Vpp',real=True)
    qprime=-(p*p+V)/3-q*q;pprime=Vp-3*q*p
    Z=q/10+s.Rational(13,20)*p
    original=Vpp+s.Rational(4,3)*(p*p+V)
    transformed=s.expand(original+s.diff(Z,q)*qprime+s.diff(Z,p)*pprime+5*q*Z-Z*Z)
    wanted=Vpp+s.Rational(13,10)*V+s.Rational(13,20)*Vp+s.Rational(39,100)*(q+s.Rational(3,2)*p)**2
    assert s.simplify(transformed-wanted)==0


def test_five_coefficients_and_four_global_exclusions():
    rows=calculation(search=False)['rows']
    assert [r['regular_compact_CDL_excluded'] for r in rows]==[False,True,True,True,True]
    assert all(float(r['exact_exclusion_certificate']['lower_bound'])>0 for r in rows[1:])
    assert all(rows[i]['HM_exponent_times_g']>rows[i+1]['HM_exponent_times_g'] for i in range(4))
    assert rows[0]['HM_exponent_times_g']==pytest.approx(403.7817410357,rel=1e-10)
    # Independent high precision evaluation of the explicit global minimum.
    with mp.workdps(60):
        for row in rows:
            d,j=row['d'],row['j'];model=VolumeParent(d,j);c=model.coefficients
            val=lambda x:mp.mpf(x.numerator)/x.denominator
            K=val(c['kinetic']);p=[d/mp.sqrt(K),(d+2)/mp.sqrt(K),2*d/mp.sqrt(K)]
            a,b,t=[val(c[k])*(w*w-mp.mpf(13)*w/20+mp.mpf(13)/10) for k,w in zip(('A','B','T'),p)]
            x=(2*b/(d*t))**(mp.mpf(2)/(d-2));v=a-(d-2)*b*x/d
            assert float(v)==pytest.approx(row['CDL_exclusion']['reduced_global_minimum'],rel=2e-10)


def test_equator_constraint_and_constant_HM_control():
    model=VolumeParent();top=model.vacuum_roots()[0]
    shot=pole_shot(model,top['s'])
    assert shot['s']==pytest.approx(top['s'],abs=1e-10)
    assert abs(shot['derivative'])<1e-10
    assert shot['radius']==pytest.approx(math.sqrt(3/top['vacuum']),rel=1e-9)
    assert abs(shot['constraint_residual'])<2e-9
