"""Retain original Hermite data and evaluate regular background derivatives.

The existing CharacteristicProfile owns the numerical interpolation. This
adapter retains its input values for exact-C2 bounds and supplies the z/T
jets needed by the same constrained mass calculation. No source is fitted.
"""
import numpy as np

from .nsc_atlas_background_cap import _budget
from .nsc_atlas_characteristic_background import CharacteristicProfile
from .nsc_atlas_cap_mass_derivatives import third_moments
from .nsc_atlas_cap_geometric_reference import _positive_gap
from .nsc_atlas_polar_cap import _scalar

SMALL_Z = .01


class RetainedCapBackground:
    """The same finite C2 profile plus its original binary Hermite inputs."""

    def __init__(self, labels, values, first, second, *, onset=None,interior_jets=None):
        self.profile = CharacteristicProfile(labels,values,first,second,onset=onset,interior_jets=interior_jets)
        self._onset=self.profile.onset
        self._inputs = {name:np.array(value,dtype=float,copy=True) for name,value in
                        zip(('labels','values','first','second'),(labels,values,first,second))}
        if self.profile.interior_jets is not None:
            self._inputs.update({'interior_'+key:value for key,value in self.profile.interior_jets.items()})
        for value in self._inputs.values():
            value.setflags(write=False)

    @property
    def onset(self):
        """Original binary onset coefficient and piece index, before conversion."""
        return None if self._onset is None else self._onset.copy()

    @property
    def inputs(self):
        """Copies of the supplied data, before interpolation-coefficient rounding."""
        return {key:value.copy() for key,value in self._inputs.items()}

    @property
    def interior_jets(self):
        return self.profile.interior_jets

    def reduced_jets(self, T, z, *, budget):
        """One bounded point: f,p,h, their T/z partials, chi,k and k_T.

        p=r*Pi, h=r*Phi and z=atanh(r). Near the centre, knot-split
        polynomial moments avoid subtraction of almost equal field values.
        The gamma series uses the existing bounded remainder implementation;
        numerical interpolation/rounding errors are not interval-enclosed.
        """
        budget = _budget(budget)
        T,z = _scalar(T,'T'),_scalar(z,'z',nonnegative=True)
        self.profile._domain(np.array([T-z,T+z]))
        r,A0 = float(np.tanh(z)),float(1/np.cosh(z)**2)
        if z and not 0 < r < 1:
            raise ValueError('background radius is not representable inside (0,1)')
        budget.check(24)
        if z == 0:
            jet = self.profile.jet(T,budget=budget)
            f,fT = 2*float(jet['F_prime']),2*float(jet['F_second'])
            return dict(r=0.,A0=1.,chi=0.,k=0.,k_T=0.,chi_z=f,k_z=fT,
                        f=f,p=0.,h=0.,f_T=fT,p_T=0.,h_T=0.,f_z=0.,p_z=fT,h_z=0.)
        with np.errstate(over='raise',invalid='raise',divide='raise'):
            if z <= SMALL_Z:
                I1,I2,J2 = self.profile._split_integrals(T,z,budget)
                J0,J1,_ = third_moments(self.profile,T,z,budget=budget)
                gamma = _positive_gap(z)
                chi,k,b = z*I1,z*I2,z*J0
                f,p,h = z*I1/r,z*I2/A0,z*(gamma*I1+J2)/A0
                fT,pT,hT = z*I2/r,z*J0/A0,z*(gamma*I2+J1)/A0
                az = I2+z*J1
            else:
                jet = self.profile.jet(np.array([T-z,T+z]),budget=budget)
                fm,fp = jet['F']; pm,pp = jet['F_prime']; sm,sp = jet['F_second']
                chi,k,b,az = float(fp-fm),float(pp-pm),float(sp-sm),float(sp+sm)
                f,p = chi/r,k/A0
                h = (pp+pm)/A0-f
                fT,pT,hT = k/r,b/A0,az/A0-k/r
            fz = A0*h/r
            pz = az/A0+2*r*p
            hz = b/A0+2*r*(f+h)-fz
        result = dict(r=r,A0=A0,chi=chi,k=k,k_T=b,chi_z=A0*(f+h),k_z=az,
                      f=f,p=p,h=h,f_T=fT,p_T=pT,h_T=hT,f_z=fz,p_z=pz,h_z=hz)
        if not all(np.isfinite(value) for value in result.values()):
            raise ValueError('nonfinite background derivative')
        budget.check()
        return result
