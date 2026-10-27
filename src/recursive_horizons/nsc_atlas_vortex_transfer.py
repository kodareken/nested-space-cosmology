"""Exact initial Reynolds transfer for compact smooth Navier--Stokes data.

This is a finite initial-rate theorem for unforced flow, inspired by the
zero-mean pulse mechanism in OpenAI's 2026 forced-blowup paper. It does not
construct a blowup or a sustained cascade. Coordinates are cylindrical;
the rotational mean averages cylindrical components with their basis carried
by the rotation group. Density, reference length and speed are one.
"""
from __future__ import annotations

from functools import lru_cache
import math

import numpy as np


KAPPA = -2.0


def bump(s):
    """C-infinity compact bump and its first two derivatives; peak value one."""
    s = np.asarray(s, dtype=float)
    values = [np.zeros_like(s) for _ in range(3)]
    keep = np.abs(s) < 1
    x = s[keep]
    d = 1-x*x
    g = np.exp(1-1/d)
    first_log = -2*x/d**2
    second_log = -2/d**2-8*x*x/d**3
    values[0][keep] = g
    values[1][keep] = g*first_log
    values[2][keep] = g*(first_log**2+second_log)
    return tuple(values)


def profile(r):
    a, ap, app = bump(2*(np.asarray(r)-1))
    return a, 2*ap, 4*app


def jets(r, z, mode, kappa=KAPPA):
    """Complex velocity amplitudes and derivatives including the phase.

    u=Re(W exp(i phase)), phase=m theta+kappa*m log(r).
    The real vector potential has A_r=0, A_theta=b chi cos(phase)/m,
    A_z=a chi cos(phase)/m, b=-a/2. Taking its curl gives W exactly.
    """
    if not isinstance(mode, (int, np.integer)) or mode <= 0:
        raise ValueError("a positive integer angular mode is required")
    r, z = np.broadcast_arrays(np.asarray(r, dtype=float), np.asarray(z, dtype=float))
    if np.any(r <= 0) or not np.all(np.isfinite(r)) or not np.all(np.isfinite(z)):
        raise ValueError("finite coordinates with r>0 required in this cylindrical chart")
    a, ap, app = profile(r)
    b, bp, bpp = -a/2, -ap/2, -app/2
    c, cp, cpp = bump(z)
    m, k = mode, float(kappa)
    W = np.array([-b*cp/m+1j*a*c/r,
                  -ap*c/m-1j*k*a*c/r,
                  (bp+b/r)*c/m+1j*k*b*c/r])
    dr = np.array([-bp*cp/m+1j*(ap/r-a/r**2)*c,
                   -app*c/m-1j*k*(ap/r-a/r**2)*c,
                   (bpp+bp/r-b/r**2)*c/m+1j*k*(bp/r-b/r**2)*c])
    dr += 1j*k*m/r*W
    dz = np.array([-b*cpp/m+1j*a*cp/r,
                   -ap*cp/m-1j*k*a*cp/r,
                   (bp+b/r)*cp/m+1j*k*b*cp/r])
    return W, dr, 1j*m*W, dz


def reynolds(r, z, mode):
    W = jets(r, z, mode)[0]
    return .5*np.real(np.einsum("i...,j...->ij...", W, np.conj(W)))


def stress_divergence(r, z, mode):
    W, dr, _, dz = jets(r, z, mode)
    conjugate = np.conj(W)
    R = .5*np.real(np.einsum("i...,j...->ij...", W, conjugate))
    Rr = .5*np.real(np.einsum("i...,j...->ij...", dr, conjugate)
                     +np.einsum("i...,j...->ij...", W, np.conj(dr)))
    Rz = .5*np.real(np.einsum("i...,j...->ij...", dz, conjugate)
                     +np.einsum("i...,j...->ij...", W, np.conj(dz)))
    return np.array([Rr[0, 0]+Rz[0, 2]+(R[0, 0]-R[1, 1])/r,
                     Rr[0, 1]+Rz[1, 2]+2*R[0, 1]/r,
                     Rr[0, 2]+Rz[2, 2]+R[0, 2]/r])


def mean_theta_rate(r, mode):
    """Actual unforced NS mean tangential acceleration at t=0,z=0.

    The rotational mean velocity is zero, hence so is the mean viscous term;
    the axisymmetric mean pressure has no theta gradient. Thus only this
    nonlinear Reynolds force remains, for every constant viscosity nu>0.
    """
    return -stress_divergence(r, 0., mode)[1]


def exact_theta_rate(r):
    a, ap, _ = profile(r)
    return -2*a*ap/np.asarray(r)**2


def angular_advection_average(r, z, mode, samples=None):
    """Independent angular sampling of the original cylindrical convection."""
    samples = samples or 4*mode+1
    theta = np.arange(samples)*(2*math.pi/samples)
    phase = np.exp(1j*(mode*theta+KAPPA*mode*math.log(r)))
    W, dr, dtheta, dz = jets(r, z, mode)
    U, Ur, Ut, Uz = [np.real(j[:, None]*phase) for j in (W, dr, dtheta, dz)]
    convection = (U[0]*Ur[1]+U[1]*Ut[1]/r+U[2]*Uz[1]+U[0]*U[1]/r)
    return -float(np.mean(convection))


def _rule(order, lo, hi):
    x, w = np.polynomial.legendre.leggauss(order)
    return (lo+hi)/2+(hi-lo)*x/2, (hi-lo)*w/2


@lru_cache(maxsize=8)
def moments(order=128):
    r, wr = _rule(order, .5, 1.5)
    z, wz = _rule(order, -1., 1.)
    a, ap, app = profile(r)
    b, bp, bpp = -a/2, -ap/2, -app/2
    c, cp, cpp = bump(z)
    I, Ip, Ipp = [float(wz @ (g*g)) for g in (c, cp, cpp)]
    k = KAPPA
    integrate = lambda g: float(wr @ g)
    E0 = math.pi/2*I*integrate(((1+k*k)*a*a+k*k*b*b)/r)
    E2 = math.pi/2*(Ip*integrate(r*b*b)+I*integrate(r*(ap*ap+(bp+b/r)**2)))
    Da, Db = app+ap/r, bpp+bp/r-b/r**2
    D2 = math.pi*(1+k*k)*I*integrate(((1+k*k)*a*a+k*k*b*b)/r**3)
    D0 = math.pi*(Ip*integrate(((2*k*k+1)*b*b+(1+k*k)*a*a)/r)
                     +I*integrate(((bp+b/r)**2+4*k*k*(ap*ap+bp*bp)
                                   -2*k*k*b*Db-2*(1+k*k)*a*Da)/r))
    Dminus2 = math.pi*(Ipp*integrate(r*b*b)+Ip*integrate(r*(ap*ap-2*b*Db))
                        +I*integrate(r*(Db*Db+Da*Da)))
    return dict(Ichi=I, Ichi_prime=Ip, Ichi_second=Ipp,
                E0=E0, E2=E2, D2=D2, D0=D0, Dminus2=Dminus2)


def energy(mode, order=128):
    q = moments(order)
    return q["E0"]+q["E2"]/mode**2


def dissipation(mode, viscosity=.001, order=128):
    if not math.isfinite(viscosity) or viscosity <= 0:
        raise ValueError("positive finite viscosity required")
    q = moments(order)
    return viscosity*(q["D2"]*mode**2+q["D0"]+q["Dminus2"]/mode**2)


def radial_fluxes(r, order=128):
    a = profile(r)[0]
    b = -a/2
    I = moments(order)["Ichi"]
    return dict(angular=-math.pi*KAPPA*a*a*I,
                axial=math.pi*KAPPA*a*b*I/r)


def direct_integrals(mode, order=128, viscosity=.001):
    """Direct 2D integrations of energy, curl loss and full gradient loss."""
    rv, wr = _rule(order, .5, 1.5)
    zv, wz = _rule(order, -1., 1.)
    r, z = rv[:, None], zv[None, :]
    W, dr, dt, dz = jets(r, z, mode)
    curl = np.array([dt[2]/r-dz[1], dz[0]-dr[2], dr[1]+W[1]/r-dt[0]/r])
    angular_gradient = np.array([(dt[0]-W[1])/r, (dt[1]+W[0])/r, dt[2]/r])
    weight = wr[:, None]*wz[None, :]*r
    norm = lambda v: np.sum(np.abs(v)**2, axis=0)
    E = math.pi/2*np.sum(weight*norm(W))
    Dcurl = viscosity*math.pi*np.sum(weight*norm(curl))
    Dgradient = viscosity*math.pi*np.sum(weight*(norm(dr)+norm(dz)+norm(angular_gradient)))
    return dict(energy=float(E), curl_loss=float(Dcurl), gradient_loss=float(Dgradient),
                quadrature_points=order*order)
