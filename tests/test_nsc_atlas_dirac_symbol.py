"""Independent metric, junction and transport checks of the radar UV input."""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_dirac_radar import WallRadarChart
from recursive_horizons.nsc_atlas_dirac_compact_source import CompactMassChart
from recursive_horizons.nsc_atlas_dirac_symbol import (
    radar_metric_jets,radar_reference_projector,PAULI_CYCLE,second_order_symbol_jump)
from recursive_horizons.nsc_spatial_reference_symbol import INDEX,SIGMA


def scalar(jet,t=0,x=0):
    return (math.factorial(t)*math.factorial(x)*jet.data[INDEX[t,x,0],0,0,0]).real


def test_fourth_radar_jets_obey_independent_vacuum_curvature_and_null_constraints():
    chart=CompactMassChart(WallRadarChart(),1e-4)
    S=chart.duration
    for U,V,side in ((.8*S,.7*S,'cap'),(.3*S,.6*S,'exterior')):
        jets=radar_metric_jets(chart,U,V,region=side)
        direct=chart.fields(U,V,region=side)
        N,r=jets['lapse'],jets['radius'];C=N*N;L=C.log()
        for jet,prefix in ((N,'lapse'),(r,'radius')):
            assert scalar(jet)==pytest.approx(direct[prefix],abs=3e-13)
            assert scalar(jet,1)==pytest.approx(direct[prefix+'_t'],abs=3e-12)
            assert scalar(jet,0,1)==pytest.approx(direct[prefix+'_x'],abs=3e-12)
        mass=chart.central_mass if side=='cap' else chart.parent_mass
        f=1-2*mass/r-r*r if side=='cap' else 1-2*mass/r
        fp=2*mass/(r*r)-2*r if side=='cap' else 2*mass/(r*r)
        fpp=-4*mass/(r*r*r)-2 if side=='cap' else -4*mass/(r*r*r)
        rt,rx=r.derivative(t=1),r.derivative(z=1)
        rtt,rtx,rxx=r.derivative(t=2),r.derivative(t=1,z=1),r.derivative(z=2)
        Lt,Lx=L.derivative(t=1),L.derivative(z=1)
        residuals=[(L.derivative(t=2)-L.derivative(z=2))/C+fpp,
            rtt-rxx+C*fp,(rx*rx-rt*rt)/C-f,
            (rtt-2*rtx+rxx-(Lt-Lx)*(rt-rx))/4,
            (rtt+2*rtx+rxx-(Lt+Lx)*(rt+rx))/4]
        # Differentiating these curvature equations twice controls the actual
        # order-four geometry, independently of the optical Taylor recipe.
        for residual in residuals:
            for t in range(3):
                for x in range(3-t):
                    assert abs(scalar(residual,t,x))<2e-9


def test_wall_second_jet_retains_the_curvature_jump_at_both_source_masses():
    reference=WallRadarChart()
    for mass in (0.,1e-4):
        chart=CompactMassChart(reference,mass);T=.4*chart.duration
        cap=radar_metric_jets(chart,T,T,region='cap')['potential_per_k']
        ext=radar_metric_jets(chart,T,T,region='exterior')['potential_per_k']
        R=float(chart.radius_at_clock(T))
        for t in range(5):
            assert scalar(ext,t)==pytest.approx(scalar(cap,t),rel=5e-13,abs=2e-10)
        for t in range(4):
            assert scalar(ext,t,1)==pytest.approx(scalar(cap,t,1),rel=5e-13,abs=2e-10)
        predicted=-6*(chart.parent_mass-mass)/R**4
        assert scalar(ext,0,2)-scalar(cap,0,2)==pytest.approx(predicted,rel=2e-12,abs=2e-10)


def test_actual_fourth_order_symbol_frame_time_correction_parity_and_joint_scaling():
    chart=CompactMassChart(WallRadarChart(),1e-4);S=chart.duration
    p=2.3*np.array([1.,1.,2.,2.]);k=np.array([1,-1,2,-2])
    row=radar_reference_projector(chart,.8*S,.7*S,p,k)
    orders=row['orders']
    for old,new in ((1,0),(2,1),(0,2)):
        assert PAULI_CYCLE@SIGMA[old]@PAULI_CYCLE.conj().T==pytest.approx(SIGMA[new],abs=2e-15)
    F=row['metric_jets']['potential_per_k']
    W=k*scalar(F);Wt=k*scalar(F,1)
    omega=np.sqrt(p*p+W*W)
    expected=(np.eye(2)-(p[:,None,None]*SIGMA[1]+W[:,None,None]*SIGMA[0])/omega[:,None,None])/2
    assert orders[0]==pytest.approx(expected,abs=3e-14)
    # Solve i*Cdot=[H,C] at first order directly: this temporal term
    # is absent from the bare instantaneous negative projector.
    correction=(-p*Wt/(4*omega**3))[:,None,None]*SIGMA[2]
    assert orders[1]==pytest.approx(correction,abs=3e-13)
    for n in range(5):
        assert orders[n,1]==pytest.approx(SIGMA[1]@orders[n,0]@SIGMA[1],abs=3e-12)
        assert orders[n,2]==pytest.approx(orders[n,0]/2**n,abs=3e-12)
        assert orders[n,3]==pytest.approx(orders[n,1]/2**n,abs=3e-12)
    assert max(float(v['scaled'].max()) for order in row['residuals'] for v in order.values())<3e-10
    with pytest.raises(ValueError,match='integer'):
        radar_reference_projector(chart,.8*S,.7*S,[1.],[.5])


def test_second_order_symbol_jump_agrees_with_independent_Pauli_Sylvester_identity():
    chart=CompactMassChart(WallRadarChart(),1e-4);T=.4*chart.duration
    p=np.array([2.3,-2.3]);k=np.ones(2,dtype=int)
    cap=radar_reference_projector(chart,T,T,p,k,region='cap')
    ext=radar_reference_projector(chart,T,T,p,k,region='exterior')
    R=float(chart.radius_at_clock(T))
    J=-6*(chart.parent_mass-chart.central_mass)/R**4
    expected=second_order_symbol_jump(p,1/R,J)
    assert np.linalg.norm(expected)>.01
    assert ext['orders'][:2]==pytest.approx(cap['orders'][:2],abs=3e-12)
    assert ext['orders'][2]-cap['orders'][2]==pytest.approx(expected,rel=3e-10,abs=3e-12)
    # The closed identity comes from the two different defects, not from
    # merely differencing the full recursive projector calculation.
    W=1/R;omega=np.sqrt(p*p+W*W)
    nW=np.stack((p*p/omega**3,-W*p/omega**3,np.zeros_like(p)),axis=-1)
    npp=np.stack((W*(2*p*p-W*W)/omega**5,-3*W*W*p/omega**5,np.zeros_like(p)),axis=-1)
    P=cap['orders'][0];Q=np.eye(2)-P
    H=p[:,None,None]*SIGMA[1]+W*SIGMA[0]
    cW=-np.einsum('pi,ijk->pjk',nW,SIGMA)/2
    cpp=-np.einsum('pi,ijk->pjk',npp,SIGMA)/2
    G=-J*(cW@cpp+cpp@cW)/8
    F=J*(SIGMA[0]@cpp-cpp@SIGMA[0])/8
    independent=-P@G@P+Q@G@Q+(H@F-F@H)/(4*omega**2)[:,None,None]
    assert expected==pytest.approx(independent,abs=3e-14)
