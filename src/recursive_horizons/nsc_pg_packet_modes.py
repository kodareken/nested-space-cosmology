"""Packet projections of computed massive PG mode fields.

Only short spatial continuations of already evaluated mode columns are used.
There is no horizon/scattering restart, metric step, or seed covariance input.
The complete bulk remains in the spectral representation; packets are probes.
"""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .nsc_common_time_bulk_split import CommonTimeBulkSplit
from .nsc_lorentzian import geometry
from .nsc_pg_lll_preparation import radial_packet
from .nsc_transmitting_dirac_domain import I2,S2,S3

PACKETS=('parent','child','child_bulk','exterior_bulk')
SIGNED_PACKET_MAP=np.kron(np.eye(4),S3)

def project_computed_modes(energy,mass,angular,mode_at_zero,mode_at_three,*,rtol=2e-12,atol=2e-14):
    """Return integral J^dagger Phi_E d rho, before the dE/(2pi) measure.

    Inputs are columns of the authenticated PG mode resolution, not a
    covariance matrix. The eight probes contain both spin components and
    therefore close under the existing signed-energy antiunitary map.
    """
    if energy<=0 or mass<0 or np.shape(mode_at_zero)!=(2,3) or np.shape(mode_at_three)!=(2,3):
        raise ValueError('positive-energy global spatial mode columns required')
    return _project_columns(energy,mass,angular,mode_at_zero,mode_at_three,rtol,atol)


def project_horizon_bases(z,mass,angular,basis_at_zero,basis_at_three,*,rtol=2e-12,atol=2e-14):
    """Analytic finite-horizon basis projection for a frequency contour.

    These are two fundamental columns, before physical horizon/infinity
    sewing. Complex z is an integration variable, not a quantum state or
    a complex metric/time history.
    """
    if np.real(z)<=0 or not np.isfinite(z) or mass<0 or np.shape(basis_at_zero)!=(2,2) or np.shape(basis_at_three)!=(2,2):
        raise ValueError('two analytic PG horizon-basis columns required')
    return _project_columns(z,mass,angular,basis_at_zero,basis_at_three,rtol,atol)


def _project_columns(energy,mass,angular,mode_at_zero,mode_at_three,rtol,atol):
    columns=np.shape(mode_at_zero)[1]
    if not np.isfinite(mode_at_zero).all() or not np.isfinite(mode_at_three).all():
        raise ValueError('spatial mode columns must be finite')
    output=np.zeros((8,columns),complex);evaluations=0
    for start,end,initial,active in ((0.,1.,mode_at_zero,(0,)),(0.,-1.,mode_at_zero,(1,2)),(3.,4.,mode_at_three,(3,))):
        direction=np.sign(end-start)
        def rhs(rho,state):
            beta,bp,_=geometry(rho)
            potential=CommonTimeBulkSplit.local_potential(1.,np.sqrt(1+rho*rho),mass,angular)
            generator=np.linalg.solve(S2-beta*I2,1j*(energy*I2-potential)+.5*bp*I2)
            field=state[:2*columns].reshape(2,columns)
            derivatives=[direction*float(radial_packet(np.array(rho),PACKETS[k])[0])*field for k in active]
            return np.r_[(generator@field).ravel(),np.array(derivatives).ravel()]
        y=np.r_[np.asarray(initial,complex).ravel(),np.zeros(len(active)*2*columns,complex)]
        run=solve_ivp(rhs,(start,end),y,method='DOP853',rtol=rtol,atol=atol)
        if not run.success:raise ArithmeticError(run.message)
        projections=run.y[2*columns:,-1].reshape(len(active),2,columns)
        for k,value in zip(active,projections):output[2*k:2*k+2]=value
        evaluations+=run.nfev
    return output,evaluations


@dataclass
class SpectralWindowProjection:
    """An evaluated window and its still unassigned complement.

    Contact reconstruction is conditional on the continuum spectral identity;
    the remainder is never replaced by chosen occupation numbers.
    """
    gram: np.ndarray
    occupied: np.ndarray
    centered: np.ndarray
    energy_maximum: float

    @property
    def contact_plus_window(self):return .5*np.eye(len(self.gram))+self.centered

    def missing_gram(self):return np.eye(len(self.gram))-self.gram

    def conditional_tail_bound(self,quadrature_error):
        """Exact completeness + positive spectral measure; error input required.

        For ||C_source-P/2||<=1/2, the omitted centered covariance has norm
        <= (||I-G_window||+quadrature_error)/2. A grid comparison alone is
        measured convergence, not automatically a rigorous quadrature bound.
        """
        if quadrature_error<0:raise ValueError('nonnegative declared quadrature error required')
        return .5*(float(np.linalg.norm(self.missing_gram(),2))+quadrature_error)

    def as_complete_covariance(self):
        raise ValueError('a finite spectral window leaves a covariance tail; it is not full C1b')


def assemble_signed_window(energies,weights,plus,minus,source_plus,projectors):
    """Paired signed-energy projection for one fixed angular sign.

    `minus` denotes the positive-energy *opposite angular* field. It is not
    an algebraic completion of a stored X margin. Zero angular sectors use
    the same physical radial family and the locked charge-conjugate bundle.
    """
    energies=np.asarray(energies);weights=np.asarray(weights)
    if not len(energies) or np.any(energies<=0) or np.any(weights<=0):
        raise ValueError('declared positive-energy quadrature required')
    gram=np.zeros((8,8),complex);occupied=gram.copy();centered=gram.copy()
    for w,fp,fm,source,P in zip(weights,plus,minus,source_plus,projectors,strict=True):
        negative=SIGNED_PACKET_MAP@fm.conj()
        g=fp@P@fp.conj().T+negative@P@negative.conj().T
        c=fp@source@fp.conj().T+negative@(P-source.conj())@negative.conj().T
        gram+=w/(2*np.pi)*g;occupied+=w/(2*np.pi)*c
        centered+=w/(2*np.pi)*(c-.5*g)
    return SpectralWindowProjection(gram,occupied,centered,float(max(energies)))
