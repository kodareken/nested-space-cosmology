"""Independent propagated-derivative and bounded-field controls."""
import numpy as np
import pytest

from recursive_horizons.nsc_atlas_dirac_mean import ExteriorPacket,MeanModeGrid,PreparedDiracSource
from recursive_horizons.nsc_atlas_dirac_scattering import null_transfer
from recursive_horizons.nsc_atlas_dirac_radar import WallRadarChart


def test_plane_wave_values_and_analytic_characteristic_derivatives():
    errors = []
    W,p = 1.7,1.1
    q = W*W/(4*p)
    for n in (8,16,32):
        h = 1/n
        t = (np.arange(n)+.5)*h
        plus = np.exp(-1j*p*t)
        minus = 2*q/W*np.exp(-1j*q*t)
        grid = MeanModeGrid(np.full((n,n),W),h,h,potential_u=np.zeros((n,n)),
            potential_v=np.zeros((n,n)),incoming_plus=plus,incoming_minus=minus,
            plus_u=-1j*p*plus,minus_v=-1j*q*minus)
        # A fixed envelope of subcell positions controls bounded queries.
        # One query whose fractional cell position changes under refinement
        # is not by itself a reliable convergence-order estimator.
        coordinates = np.concatenate([(np.arange(n//4,3*n//4-1)+.5+f)*h
                                      for f in (0.,.25,.5,.75)])
        u,v = np.meshgrid(coordinates,coordinates,indexing='ij')
        expected = np.exp(-1j*(p*u+q*v))[...,None]*np.array([1.,2*q/W])
        values = grid.evaluate(u,v)
        du = grid.evaluate(u,v,derivative_u=1)
        dv = grid.evaluate(u,v,derivative_v=1)
        errors.append(max(np.max(np.abs(values-expected)),
                          np.max(np.abs(du+1j*p*expected)),
                          np.max(np.abs(dv+1j*q*expected))))
        assert grid.flux_error < 32*n*np.finfo(float).eps*grid.incoming_flux
    assert errors[-1] < 2e-3
    assert 3.2 < errors[0]/errors[1] < 4.8
    assert 3.2 < errors[1]/errors[2] < 4.8


def test_zero_energy_exact_solution_across_C1_wall():
    def W(x):
        return .7+3.2*np.maximum(x,0)**2
    def amplitude(x):
        return np.exp(-.7*x-3.2*np.maximum(x,0)**3/3)
    errors = []
    for n in (8,16,32):
        h = 1/n
        t = (np.arange(n)+.5)*h
        x = (t[None,:]-t[:,None])/2
        plus = amplitude(-t/2)/np.sqrt(2)
        minus = -1j*amplitude(t/2)/np.sqrt(2)
        grid = MeanModeGrid(W(x),h,h,potential_u=-3.2*np.maximum(x,0),
            potential_v=3.2*np.maximum(x,0),incoming_plus=plus,incoming_minus=minus,
            plus_u=W(-t/2)*plus/2,minus_v=-W(t/2)*minus/2)
        u,v = np.array([.37,.63,.5]),np.array([.63,.37,.5])
        x0 = (v-u)/2
        expected = amplitude(x0)[:,None]*np.array([1.,-1j])/np.sqrt(2)
        errors.append(max(np.max(np.abs(grid.evaluate(u,v)-expected)),
            np.max(np.abs(grid.evaluate(u,v,derivative_u=1)-W(x0)[:,None]*expected/2)),
            np.max(np.abs(grid.evaluate(u,v,derivative_v=1)+W(x0)[:,None]*expected/2))))
    assert errors[-1] < 2e-3
    assert 3.1 < errors[0]/errors[1] < 5.
    assert 3.1 < errors[1]/errors[2] < 5.


def test_same_CAR_transfer_flux_weights_and_bounded_evaluation():
    rng = np.random.default_rng(4026)
    w = rng.uniform(-1,1,size=(6,7))
    plus = rng.normal(size=6)+1j*rng.normal(size=6)
    minus = rng.normal(size=7)+1j*rng.normal(size=7)
    grid = MeanModeGrid(w,.1,.13,potential_u=np.zeros_like(w),
        potential_v=np.zeros_like(w),incoming_plus=plus,incoming_minus=minus,
        plus_u=np.zeros(6),minus_v=np.zeros(7))
    transfer,_ = null_transfer(w,.1,.13)
    direct = transfer@np.r_[np.sqrt(.1)*plus,np.sqrt(.13)*minus]
    assert np.max(np.abs(direct-grid.outgoing)) < 32*np.finfo(float).eps*max(1,np.linalg.norm(direct))
    assert grid.flux_error < 32*13*np.finfo(float).eps*grid.incoming_flux
    for u,v in ((0.,.4),(.3,.95),(.7,.3)):
        with pytest.raises(ValueError):
            grid.evaluate(u,v)


def test_packet_continuum_normalization_and_derivative():
    S=.1256873632901346
    packet = ExteriorPacket(S)
    nodes,weights = np.polynomial.legendre.leggauss(128)
    values,derivatives = packet.data((.6+.3*nodes)*S)
    assert .3*S*np.dot(weights,np.abs(values)**2) == pytest.approx(1.,abs=3e-14)
    point=.62*S
    h=2e-6*S
    samples = [packet.data(point+a*h)[0] for a in (-2,-1,1,2)]
    fd = (samples[0]-8*samples[1]+8*samples[2]-samples[3])/(12*h)
    assert fd == pytest.approx(packet.data(point)[1],rel=3e-10,abs=2e-7)
    quiet,quiet_derivative = packet.data(np.array([0.,.29*S,.91*S,S]))
    assert np.array_equal(quiet,np.zeros(4))
    assert np.array_equal(quiet_derivative,np.zeros(4))


def test_source_interface_causal_zero_incoming_energy_and_transmission():
    # Interface/known-limit check only: no clock or mass target is evaluated.
    chart = WallRadarChart();S = chart.duration
    source = PreparedDiracSource(chart,16)
    quiet = source.evaluate(.4*S,.2*S)
    for value in quiet['orthonormal'].values():
        assert value == 0.
    incoming = source.evaluate(.1*S,.6*S)['orthonormal']
    assert incoming['energy_density']>0
    assert incoming['outward_energy_flux']<0
    cap = source.evaluate(.6*S,.6*S,region='cap')['orthonormal']
    exterior = source.evaluate(.6*S,.6*S,region='exterior')['orthonormal']
    for key in cap:
        assert cap[key] == pytest.approx(exterior[key],rel=2e-10,abs=2e-12)
    assert -incoming['energy_density']+incoming['radial_pressure']+2*incoming['angular_pressure'] == pytest.approx(
        0.,abs=2e-11*incoming['energy_density'])
