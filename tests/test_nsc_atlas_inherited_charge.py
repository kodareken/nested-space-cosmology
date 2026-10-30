import math
import numpy as np
import pytest
import sympy as s
from recursive_horizons.nsc_atlas_inherited_charge import harmonic_sector,pair_transport,calculation


def test_warped_action_factors_and_fixed_curvature_radius():
    d,v=s.symbols('d v',positive=True)
    volume=s.exp(-2*d*v)*s.exp(d*v)
    assert s.simplify(volume*s.exp(d*v))==1
    assert s.simplify(volume*s.exp(-2*v)-s.exp(-(d+2)*v))==0
    assert s.Rational(12*11,15)==s.Rational(44,5)
    assert 12/s.Rational(44,5)==s.Rational(15,11)
    row=harmonic_sector()
    assert row['charged_mass_squared']>row['zero_mode_mass_squared']>0


def test_harmonic_pair_comes_from_one_real_field_and_proper_rotation():
    n1,n2,x,y=s.symbols('n1 n2 x y',real=True)
    Y=(n1+s.I*n2)/s.sqrt(2);psi=(x+s.I*y)/s.sqrt(2)
    assert s.expand(psi*s.conjugate(Y)+s.conjugate(psi)*Y)==x*n1+y*n2
    K=lambda f:-n2*s.diff(f,n1)+n1*s.diff(f,n2)
    assert s.simplify(K(Y)-s.I*Y)==0
    # Ambient sphere formula: Delta_S P = projected Hessian - d n.grad P.
    assert s.simplify(-12*(n1*s.diff(Y,n1)+n2*s.diff(Y,n2))+12*Y)==0
    P=s.eye(13);P[1,1]=-1;P[12,12]=-1
    J=s.zeros(13);J[0,1]=-1;J[1,0]=1
    assert P.det()==1 and P*J*P.T==-J


def test_horizon_radial_equation_from_covariant_wave_operator():
    r,v,w,M,l=s.symbols('r v w M l',real=True)
    f=s.Function('f')(r);u=s.Function('u')(r)
    field=s.exp(-s.I*w*v)*u/r
    box=2*s.diff(field,v,r)+2*s.diff(field,v)/r+f*s.diff(field,r,2)
    box+=(2*f/r+s.diff(f,r))*s.diff(field,r)-l*(l+1)*field/r**2
    expected=f*s.diff(u,r,2)+(s.diff(f,r)-2*s.I*w)*s.diff(u,r)
    expected-=(s.diff(f,r)/r+M*M+l*(l+1)/r**2)*u
    assert s.simplify(r*s.exp(s.I*w*v)*(box-M*M*field)-expected)==0
    p,k=s.symbols('p k')
    indicial=k*p*(p-1)+(k-2*s.I*w)*p
    assert s.factor(indicial)==p*(k*p-2*s.I*w)


def test_euclidean_kernel_selects_covariance_not_an_ellipse_phase():
    row=harmonic_sector()
    for L in (0,1,2,10,100):
        eig=row['H']**2*L*(L+3)+row['charged_mass_squared']
        covariance=np.linalg.inv(np.eye(2)*eig)
        assert covariance[0,1]==0
        assert covariance[0,0]==covariance[1,1]>0
    assert calculation()['euclidean_preparation']['coherent_phase_offset'] is None


def test_neutral_bogoliubov_charge_identity_with_nonzero_input():
    for q in (.1,1.,2.):
        neutral=pair_transport(0,0,q)
        assert neutral['plus']==neutral['minus']
        charged=pair_transport(2,1,q)
        assert charged['plus']-charged['minus']==pytest.approx(1,abs=4e-14)
    # Polynomial coefficients of operator number difference, including pair cross terms.
    A,B,np_,nm,one,cross=s.symbols('A B np nm one cross')
    outplus=A*np_+B*(nm+one)+cross
    outminus=A*nm+B*(np_+one)+cross
    assert s.expand(outplus-outminus-(A-B)*(np_-nm))==0


def test_state_conjugation_preserves_external_quadratic_stress():
    # y -> -y flips j while each external stress bilinear is unchanged.
    x,y,dx,dy,ex,ey=s.symbols('x y dx dy ex ey',real=True)
    flip={y:-y,dy:-dy,ey:-ey}
    J=x*dy-y*dx;stress=dx*ex+dy*ey
    assert s.expand(J.subs(flip,simultaneous=True)+J)==0
    assert s.expand(stress.subs(flip,simultaneous=True)-stress)==0


def test_no_thermal_factor_is_relabelled_eta_or_an_observed_scale():
    r=calculation()
    assert r['thermal_control']['log_excitation_to_deexcitation_ratio']<-500
    assert r['thermal_control']['equal_for_both_charges'] is True
    assert r['imbalance']['conditional_miss_percent']==-100
    assert r['imbalance']['physical_baryon_photon_eta'] is None
    assert r['sector']['shape_selected_by_action'] is False
    assert r['requested_nonzero_eta_derived'] is False


def test_local_thermal_fluctuations_survive_zero_mean_charge():
    q=s.symbols('q',positive=True)
    mean=q/(1-q);second=q*(1+q)/(1-q)**2
    assert s.simplify(2*(second-mean**2)-2*mean*(1+mean))==0
    prefactor=(1-q)/(1+q)
    assert s.simplify(prefactor*(1+2*q/(1-q))-1)==0
    assert s.simplify(prefactor*2*q*(1+q)/(1-q)**3-2*q/(1-q)**2)==0
    r=calculation()['thermal_control']
    assert r['single_stationary_mode_charge_rms']>0
    assert math.log(r['single_stationary_mode_charge_rms'])==pytest.approx(r['log_charge_variance']/2)
