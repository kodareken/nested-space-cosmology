"""Inherited curvature audit and an explicitly imposed flat-volume control."""
import math
from fractions import Fraction as F
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from .nsc_atlas_parent_scale_selector import VolumeParent,evolve
from .nsc_atlas_inherited_matter import exact_initial_geometry
from .nsc_atlas_volume_observations import VolumeBackground


class FlatVolumeModel(VolumeParent):
    def hubble(self,N,s,v,E):
        return math.sqrt((self.K*v*v/2+self.potential(s)+E)/3)


class ClosedVolumeComparison(VolumeBackground):
    name='closed_volume_12_2'
    def fractions(self,theta):
        result=super().fractions(theta)
        result['Omega_dark_energy']=result['Omega_Lambda']+result['Omega_volume']
        return result


class FlatVolumeControl(ClosedVolumeComparison):
    """C=0 is imposed before fitting; no inherited flatness is asserted."""
    name='imposed_flat_volume_12_2'
    def __init__(self,zmax,*,rtol=2e-11):
        self.model=FlatVolumeModel();self.C=0.
        self.Vstar=self.model.vacuum_roots()[1]['vacuum']
        self.minimum_N=-math.log1p(zmax)-1e-8;self.calls=0
        def rhs(N,y):
            self.calls+=1
            if self.calls>30000:raise RuntimeError('flat-control allowance exhausted')
            s,v,J=y;dy=self.model.averaged_rhs(0.,[N,s,v]);H=dy[0]
            return v/H,dy[2]/H,math.exp(-N)/H
        self.past=solve_ivp(rhs,(0,self.minimum_N),[0,0,0],method='DOP853',rtol=rtol,atol=rtol/100,dense_output=True)
        self.future=solve_ivp(rhs,(0,self.bounds[1]),[0,0,0],method='DOP853',rtol=rtol,atol=rtol/100,dense_output=True)
        if not self.past.success or not self.future.success:raise RuntimeError('flat-control domain unresolved')

    def distances(self,theta,z,hel):
        z=np.asarray(z);ob=self.state([theta]);em=self.state(theta-np.log1p(z))
        H=ob[5][0];delta=ob[2][0]-em[2]
        if np.any(delta<=0):raise ValueError('nonpositive light-cone distance')
        DM=H*math.exp(theta)*delta
        return (1+np.asarray(hel))*DM,DM,H/em[5]


def exact_energy_obstruction():
    """A coarse rational certificate; no precision campaign is needed."""
    def exp_bounds(x):
        lo=sum(x**n/F(math.factorial(n)) for n in range(9))
        return lo,lo+x**9/F(math.factorial(9))/(1-x/10)
    e1=exp_bounds(F(1));e76=exp_bounds(F(7,6));e2=exp_bounds(F(2))
    potential_lower=7*e1[0]-F(15,2)*e76[1]+F(25,32)*e2[0]
    assert potential_lower>F(21,32)
    C=exact_initial_geometry()['curvature']
    ratio=3*C/F(3,8)/exp_bounds(F(1,2))[1]
    return dict(field_lower_bound='-1/12',potential_at_bound_lower=str(potential_lower),
        curvature_over_matter_lower_exact=str(ratio),curvature_over_matter_lower=float(ratio),
        matter_upper_at_curvature_0_002=float(F(1,500)/ratio))


def curvature_audit():
    model=VolumeParent();_,solution=evolve(endpoint=6.)
    C=float(exact_initial_geometry()['curvature']);star=model.vacuum_roots()[1]
    def fractions(N):
        _,s,v=solution.sol(N)[:3];E=model.matter(N,s);H=model.hubble(N,s,v,E)
        return dict(N=N,Omega_k=-C*math.exp(-2*N)/(H*H),Omega_m=E/(3*H*H),
            Omega_Lambda=star['vacuum']/(3*H*H),
            Omega_volume=(model.K*v*v/2+model.potential(s)-star['vacuum'])/(3*H*H),s=s)
    threshold=brentq(lambda N:fractions(N)['Omega_k']+.002,1.,6.,xtol=1e-12)
    lower_s=brentq(lambda s:model.potential(s)-21/32,-.3,star['s'],xtol=1e-14)
    ratio_lower=3*C/.375*math.exp(6*lower_s)
    return dict(schema='nsc-atlas-curvature-mechanisms-v1',
        route_dispositions=[
            'Open CDL: an open slicing fixes the sign k=-1, not the observed small magnitude; the compatible nucleation state, slow-roll duration and exit history are missing.',
            'Inherited oscillatory matter: dust domination increases |Omega_k| approximately as a; late vacuum acceleration dilutes it but also removes the required matter fraction.',
            'Nested similarity: a->r*a and H->H/r preserve Omega_k; a dimensionless curvature/state transfer rule and an observer level are not supplied.'],
        threshold_crossing=fractions(threshold),threshold_is_selected=False,
        declared_N6_endpoint=fractions(6.),
        reference_observer=fractions(.43923013327691346),
        energy_bound=dict(initial_total_density='21/32',volume_field_lower_bound=lower_s,
            absolute_curvature_over_matter_lower_bound=ratio_lower,
            implication='For any expanding N>=0 in this fixed source, |Omega_k|/Omega_m=(3C/E0)*exp(N+6s) exceeds this bound. Small curvature together with substantial matter requires changed preparation or new physics.',
            matter_upper_if_curvature_below_0_002=.002/ratio_lower),
        exact_coarse_energy_obstruction=exact_energy_obstruction(),
        curvature_derivative='d ln|Omega_k|/dN=-2-2*Hdot/H^2=2*q_deceleration',
        inherited_observed_small_curvature_derived=False,
        flat_control='C=0 is an explicitly imposed diagnostic initial condition, never labelled inherited or fitted.')
