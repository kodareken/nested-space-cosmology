"""Exact projection oracle for coalesced quiet cuts; certificate untouched."""
from fractions import Fraction as Q
from math import comb

import numpy as np
from numpy.polynomial import Polynomial
import pytest
import sympy as s
from flint import arb,ctx

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_cap_background_jets import RetainedCapBackground,SMALL_Z
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
import recursive_horizons.nsc_atlas_cap_fixed_reference as reference
from recursive_horizons.nsc_atlas_cap_background_bound import _background_controls,_piece_ranges
from recursive_horizons.nsc_atlas_cap_variational_reference import _exact_ball

CONTROLS=ControlBudget(maximum_evaluations=3000000,cpu_seconds=60.)


def test_quiet_cut_removal_keeps_exact_projection_and_original_certificate(monkeypatch):
    with CONTROLS.active() as budget:
        knots=np.array([-1.,-.75,-.5,0.,.5]);A=Q(1,64);B=Q(1,512)
        poly=float(A)*Polynomial([.5,1])**4+float(B)*Polynomial([.5,1])**7
        jets=[poly.deriv(j)(knots) for j in range(3)]
        for row in jets: row[:3]=0.
        bg=RetainedCapBackground(knots,*jets,onset=dict(index=2,quartic_coefficient=float(A)))
        saved=bg.inputs;basis=FixedCapBasis([0,1],budget=budget)
        Z=Q(1,2);T=-.6;time=Q(T);front=Q(-.5)-time
        full=sorted({Q(0),Q(SMALL_Z),Z,*[abs(Q(float(k))-time) for k in knots if 0<abs(Q(float(k))-time)<Z]})
        original=np.array(list(map(float,full)))
        reduced,owners=reference._background_partition(basis,Z,T,bg)
        removed=float(abs(Q(-.75)-time))
        assert removed in original and removed not in reduced
        assert float(front) in reduced and float(abs(Q(-1)-time)) in reduced
        assert len(reduced)==len(original)-1
        assert np.all(owners==0)
        zero=np.zeros(basis.dimension)
        opts=dict(alpha=.002,T_x=.8,radius_rate=.001,incoming_w_plus=.01,
                  background=bg,background_time=T,budget=budget)
        small=reference.fixed_cap_reference(basis,zero,zero,Z,**opts)
        monkeypatch.setattr(reference,'_background_partition',
                            lambda *args:(original.copy(),np.zeros(len(original)-1,dtype=int)))
        full_row=reference.fixed_cap_reference(basis,zero,zero,Z,**opts)
        assert small['background']['radial_intervals']==full_row['background']['radial_intervals']-1
        for key in ('mass','lapse','mass_rate','explicit_time_remainder'):
            assert small[key]==pytest.approx(full_row[key],rel=0.,abs=2e-12)
        assert small['gradient']['T']==pytest.approx(full_row['gradient']['T'],rel=0.,abs=2e-12)
        for row in (small,full_row): assert abs(row['remainder_identity_defect'])<2e-12
        assert small['metric_certificate']==full_row['metric_certificate']
        # Independent exact integration of phi(z)*F'(T+z) and phi(z)*F''(T+z).
        # The other characteristic foot is quiet on this whole radial interval.
        shapes=basis.exact_shapes(0);n=basis.dimension
        M=s.Matrix(n,n,lambda i,j:sum(shapes[p][i]*shapes[q][j]*Z**(p+q+1)/Q(p+q+1)
                                     for p in range(6) for q in range(6)))
        for name,terms in [('projected_time_chi',((3,4*A),(6,7*B))),
                           ('projected_time_density',((2,12*A),(5,42*B)))]:
            loads=[]
            for j in range(n):
                value=sum(shapes[p][j]*c*comb(power,k)*(-front)**(power-k)*
                          (Z**(p+k+1)-front**(p+k+1))/Q(p+k+1)
                          for p in range(6) for power,c in terms for k in range(power+1))
                loads.append(value)
            budget.check(2048)
            exact=np.array(M.inv()*s.Matrix(loads),dtype=float).ravel()
            for row in (small,full_row):
                np.testing.assert_allclose(row['background'][name],exact,rtol=0.,atol=2e-11)
        with ctx.workprec(160):
            kb=list(map(_exact_ball,knots));v,d,e=([_exact_ball(x) for x in row] for row in jets)
            controls=_background_controls(bg,kb,v,d,e)
            bounds=_piece_ranges(kb,controls,arb(-.75),arb(-.5),budget)
            # Even a quiet query ending at the onset sees the full active hull.
            for bound,coefficients in zip(bounds,controls[2]):
                assert bound==max(abs(c).upper() for c in coefficients)
        assert bg.onset==dict(index=2,quartic_coefficient=float(A))
        for name,value in saved.items():np.testing.assert_array_equal(bg.inputs[name],value)
        print('quiet partition:',len(original)-1,'to',len(reduced)-1,'radial panels; certificate',
              small['metric_certificate']['interval_count'],'unchanged')
    print('quiet partition controls:',CONTROLS.budget.evaluations,'charged values;',
          CONTROLS.active_cpu_seconds,'active CPU seconds')
