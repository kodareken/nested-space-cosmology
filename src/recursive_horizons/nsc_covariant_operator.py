"""Static four-dimensional covariant proper-time Dirac modulus and variation.

Finite spatial Fourier regulator and angular/frequency quadrature are removed
by independent refinements. The physical proper-time cutoff remains explicit.
No invariant finite completion, compensator or metric saddle is assumed.
"""
from dataclasses import dataclass
from functools import cached_property
from math import pi, sqrt

import numpy as np
from scipy.linalg import eigh
from scipy.special import exp1, erfc
from numpy.polynomial.legendre import leggauss
from .nsc_regulated import OperatorConventions, RegulatedOperator


@dataclass(frozen=True)
class CovariantStaticMetric:
    length: float
    lapse: np.ndarray
    radial_scale: np.ndarray
    sphere_radius: np.ndarray
    eta: float = .5

    def __post_init__(self):
        if not np.isfinite(self.length) or self.length <= 0 or self.eta not in (0.,.5):
            raise ValueError('positive length and P/AP spin structure required')
        sizes=[]
        for name in ('lapse','radial_scale','sphere_radius'):
            if np.iscomplexobj(getattr(self,name)):raise ValueError('metric arrays must be real')
            value=np.array(getattr(self,name),dtype=float,copy=True)
            if value.ndim!=1 or len(value)<9 or not np.isfinite(value).all() or min(value)<=0:
                raise ValueError('finite positive metric arrays with at least nine samples required')
            value.setflags(write=False);object.__setattr__(self,name,value);sizes.append(len(value))
        if len(set(sizes))!=1:raise ValueError('metric array sizes differ')
        if (self.eta==0. and sizes[0]%2!=1) or (self.eta==.5 and sizes[0]%2!=0):
            raise ValueError('use odd P or even AP grid to retain symmetric signed momenta')

    @property
    def points(self):return len(self.lapse)
    @property
    def spacing(self):return self.length/self.points
    @property
    def x(self):return (np.arange(self.points)-self.points//2)*self.length/self.points

    @cached_property
    def momenta(self):
        indices=np.arange(self.points)-self.points//2
        return 2*pi*(indices+self.eta)/self.length

    @cached_property
    def momentum_matrix(self):
        fourier=np.exp(1j*self.x[:,None]*self.momenta[None,:])/sqrt(self.points)
        matrix=(fourier*self.momenta[None,:])@fourier.conj().T
        return (matrix+matrix.conj().T)/2


def smooth_metric(points=32, radius_parameter=2., throat_radius=1., eta=.5, general=False):
    length=2*radius_parameter
    x=(np.arange(points)-points//2)*length/points
    wave=pi/length
    radius=np.sqrt(throat_radius**2+(np.sin(wave*x)/wave)**2)
    angle=2*pi*x/length
    lapse=np.exp(.08*np.cos(angle)+.03*np.sin(2*angle)) if general else np.ones(points)
    radial=np.exp(.06*np.sin(angle)-.02*np.cos(2*angle)) if general else np.ones(points)
    return CovariantStaticMetric(length,lapse,radial,radius,eta)


def cylinder_metric(points=32,length=4.,radius=1.,eta=.5):
    return CovariantStaticMetric(length,np.ones(points),np.ones(points),np.full(points,radius),eta)


SIGMA1=np.array([[0,1],[1,0]],complex)
SIGMA2=np.array([[0,-1j],[1j,0]],complex)
SIGMA3=np.diag([1.,-1.])


def radial_momentum(metric):
    inverse=1/metric.radial_scale
    p=metric.momentum_matrix
    return .5*(inverse[:,None]*p+p*inverse[None,:])


def euclidean_operator(metric,omega,kappa=1):
    """One angular-sign block after χ=r sqrt(Nq) ψ; full weight 4κ.

    Dω=ρ3 ω/N+ρ1 p_q−ρ2 κ/r. The paired ±κ 4-spinor operator
    anticommutes physical Γ5=−η1ρ2. Static real metric and symmetric Fourier
    modes give equal modulus traces for the two angular signs and ±ω.
    """
    if not np.isfinite(omega) or isinstance(kappa,bool) or not isinstance(kappa,(int,np.integer)) or kappa==0:
        raise ValueError('finite frequency and nonzero integer angular label required')
    return (np.kron(SIGMA3,np.diag(omega/metric.lapse))
            +np.kron(SIGMA1,radial_momentum(metric))
            -np.kron(SIGMA2,np.diag(kappa/metric.sphere_radius)))


def frequency_extent(metric,cutoff):
    """Finite-operator spectral lower estimate at the quadrature endpoint.

    Dω² >= ω²/Nmax²−|ω| ||[p_q,1/N]||, since the remaining spatial
    square is positive. Choose endpoint where this lower estimate is 64Λ².
    Refinements, not this endpoint value alone, control the omitted integral.
    """
    pq=radial_momentum(metric);ni=1/metric.lapse
    c=float(np.linalg.norm(pq*ni[None,:]-ni[:,None]*pq,2))
    nmax=float(max(metric.lapse))
    return .5*nmax*nmax*(c+sqrt(c*c+256*cutoff*cutoff/(nmax*nmax)))


def frequency_energy_tail_bound(metric,cutoff,angular_max,extent):
    """Analytic upper bound on the omitted finite-angular energy integral.

    E1(u)<=exp(-u)/u and the quadratic lower eigenvalue bound is convex.
    This bound concerns energy, not derivatives or infinite angular/spatial
    truncations. Displayed evaluation is binary64, not interval-certified.
    """
    pq=radial_momentum(metric);ni=1/metric.lapse
    c=float(np.linalg.norm(pq*ni[None,:]-ni[:,None]*pq,2))
    a=float(max(metric.lapse))**2
    q=extent**2/a-c*extent;slope=2*extent/a-c
    if q<=0 or slope<=0:raise ValueError('frequency endpoint has no positive spectral bound')
    sum_kappa=angular_max*(angular_max+1)/2
    return float(4*metric.points*sum_kappa/pi*cutoff**4/(q*slope)*np.exp(-q/cutoff**2))


def frequency_density(metric,omega,kappa,cutoff,gradients=True):
    matrix=euclidean_operator(metric,omega,kappa)
    values,vectors=eigh(matrix,driver='evd',check_finite=False)
    if min(abs(values))<1e-10:
        raise ValueError('unresolved zero mode requires a separate determinant/domain prescription')
    value=.5*float(np.sum(exp1((values/cutoff)**2)))
    if not gradients:return value,None
    derivative=(vectors*(-np.exp(-(values/cutoff)**2)/values)[None,:])@vectors.conj().T
    blocks=derivative.reshape(2,metric.points,2,metric.points)
    trace=lambda gamma:np.einsum('ba,aibj->ij',gamma,blocks)
    g_n=-omega*np.diag(trace(SIGMA3)).real/metric.lapse**2
    t=trace(SIGMA1);p=metric.momentum_matrix
    g_q=-.5*np.diag(t@p+p@t).real/metric.radial_scale**2
    # Angular term is -rho2*kappa/r, hence positive radius derivative.
    g_r=kappa*np.diag(trace(SIGMA2)).real/metric.sphere_radius**2
    return value,np.array([g_n,g_q,g_r])


def cutoff_response(metric,cutoff=2.,angular_max=8,frequency_points=48,extent_factor=1.,gradients=True):
    if not np.isfinite(cutoff) or cutoff<=0 or isinstance(angular_max,bool) or not isinstance(angular_max,int) or angular_max<1:
        raise ValueError('positive cutoff and angular maximum required')
    if not isinstance(frequency_points,int) or frequency_points<8 or extent_factor<1:
        raise ValueError('resolved quadrature and extent factor >=1 required')
    extent=frequency_extent(metric,cutoff)*extent_factor
    nodes,weights=leggauss(frequency_points)
    frequencies=(nodes+1)*extent/2;weights=weights*extent/2
    total=0.;gradient=np.zeros((3,metric.points));sectors=[]
    for kappa in range(1,angular_max+1):
        energy=0.;sector_gradient=np.zeros_like(gradient)
        for omega,weight in zip(frequencies,weights):
            value,g=frequency_density(metric,omega,kappa,cutoff,gradients)
            factor=4*kappa*weight/pi
            energy+=factor*value
            if gradients:sector_gradient+=factor*g
        total+=energy;gradient+=sector_gradient
        sectors.append({'kappa':kappa,'energy':float(energy)})
    result={'energy':float(total),'angular_sectors':sectors,'frequency_extent':extent,
            'finite_angular_frequency_energy_tail_bound':frequency_energy_tail_bound(metric,cutoff,angular_max,extent)}
    if gradients:
        n,q,r,h=metric.lapse,metric.radial_scale,metric.sphere_radius,metric.spacing
        rho=gradient[0]/(4*pi*q*r*r*h)
        px=-gradient[1]/(4*pi*n*r*r*h)
        pt=-gradient[2]/(8*pi*n*q*r*h)
        result.update({'metric_gradients':gradient,'rho':rho,'p_x':px,'p_perp':pt,
                       'radial_null':rho+px,'trace':rho-px-2*pt,
                       'lapse_homogeneity_residual':float(np.dot(n,gradient[0])-total)})
    return result


def ultrastatic_reference(metric,cutoff=2.,angular_max=8):
    if not np.all(metric.lapse==1):raise ValueError('reference assumes unit constant lapse')
    total=0.
    for kappa in range(1,angular_max+1):
        values=np.linalg.eigvalsh(euclidean_operator(metric,0.,kappa))
        positive=values[values>0]
        total+=4*kappa*np.sum(cutoff/np.sqrt(pi)*np.exp(-(positive/cutoff)**2)
                              -positive*erfc(positive/cutoff))
    return float(total)


def cylinder_zero_winding_energy(length=4.,radius=1.,cutoff=2.,angular_max=64):
    k=np.arange(1,angular_max+1,dtype=float);m2=(k/radius)**2
    return float(length/pi*np.sum(k*(cutoff**2*np.exp(-m2/cutoff**2)-m2*exp1(m2/cutoff**2))))


def directional_variation(metric,log_direction,cutoff=2.,angular_max=16,frequency_points=48):
    """First/second metric variations from the existing noncommuting calculus.

    N,q,r -> exp(epsilon*direction)*(N,q,r). Only the finite matrix
    variations are reused; its additive rank normalization is NOT integrated.
    """
    direction=np.asarray(log_direction,dtype=float)
    if direction.shape!=(3,metric.points) or not np.isfinite(direction).all():
        raise ValueError('three finite log-metric directions required')
    n,a,b=direction;p=metric.momentum_matrix
    pq1=-.5*((a/metric.radial_scale)[:,None]*p+p*(a/metric.radial_scale)[None,:])
    pq2=.5*((a*a/metric.radial_scale)[:,None]*p+p*(a*a/metric.radial_scale)[None,:])
    extent=frequency_extent(metric,cutoff)
    nodes,weights=leggauss(frequency_points);frequencies=(nodes+1)*extent/2;weights=weights*extent/2
    gradient=hessian=0.
    for kappa in range(1,angular_max+1):
        w=kappa/metric.sphere_radius
        for omega,weight in zip(frequencies,weights):
            f=omega/metric.lapse
            first=(-np.kron(SIGMA3,np.diag(f*n))+np.kron(SIGMA1,pq1)
                   +np.kron(SIGMA2,np.diag(w*b)))
            second=(np.kron(SIGMA3,np.diag(f*n*n))+np.kron(SIGMA1,pq2)
                    -np.kron(SIGMA2,np.diag(w*b*b)))
            op=RegulatedOperator(euclidean_operator(metric,omega,kappa),OperatorConventions(cutoff=cutoff))
            g,h=op.variation(first,second)
            factor=4*kappa*weight/pi;gradient+=factor*g;hessian+=factor*h
    return {'first':float(gradient),'second':float(hessian)}


def heat_trace(metric,proper_time,angular_max=24,frequency_points=48):
    if not np.isfinite(proper_time) or proper_time<=0:raise ValueError('positive proper time required')
    cutoff=1/sqrt(proper_time)
    extent=frequency_extent(metric,cutoff)
    nodes,weights=leggauss(frequency_points);frequencies=(nodes+1)*extent/2;weights=weights*extent/2
    total=0.
    for kappa in range(1,angular_max+1):
        for omega,weight in zip(frequencies,weights):
            values=np.linalg.eigvalsh(euclidean_operator(metric,omega,kappa))
            total+=4*kappa*weight/pi*np.sum(np.exp(-proper_time*values**2))
    return float(total)


def heat_coefficients(metric):
    """Integrated a0,a2,a4 per coordinate time in the covariant Dirac convention.

    Reuse independently derived curvature integrals on the same even grid;
    R_E=-R_L, while squared invariants are unchanged. The closed-cell
    Laplacian integrates to zero; no boundary coefficient is suppressed.
    """
    from .nsc_finite_terms import finite_response
    from .nsc_shape_response import StaticAxialMetric
    if metric.points%2:raise ValueError('curvature control uses the existing even metric grid')
    existing=StaticAxialMetric(metric.length,metric.lapse,metric.radial_scale,metric.sphere_radius)
    f=finite_response(existing)['energy']
    return {'a0':float(4*f[0]),'a2':float(f[1]/3),
            'a4':float((-18*f[2]+11*f[4])/360),
            'integrated_E4':float(f[4]),'integrated_boxR':float(f[5])}


def ultraviolet_subtraction(metric,cutoff,normalization=1.):
    """Known divergent local functional, not a fitted finite completion."""
    from .nsc_finite_terms import finite_response
    from .nsc_shape_response import StaticAxialMetric
    if normalization<=0 or cutoff<=0:raise ValueError('positive independent scales required')
    existing=StaticAxialMetric(metric.length,metric.lapse,metric.radial_scale,metric.sphere_radius)
    local=finite_response(existing)
    logarithm=np.log(cutoff**2/normalization**2)
    # R_E=-R_L. The integrated box term vanishes on this closed domain.
    coefficients=np.array([2*cutoff**4,cutoff**2/3,-logarithm/20,0.,11*logarithm/360,0.])/(32*pi*pi)
    energy=float(np.dot(coefficients,local['energy']))
    gradient=metric.spacing*np.einsum('k,kij->ij',coefficients,local['gradients'])
    return {'energy':energy,'metric_gradients':gradient,
            'coefficients_in_Lorentzian_finite_term_basis':coefficients,
            'local_normalization_log_derivative':-heat_coefficients(metric)['a4']/(16*pi*pi),
            'remainder_normalization_log_derivative':heat_coefficients(metric)['a4']/(16*pi*pi)}


def subtracted_response(metric,cutoff=2.,normalization=1.,angular_max=16,frequency_points=48):
    """Compute an explicitly defined subtraction remainder and its variation.

    This supplies the determinant part for matching. Setting any remaining
    invariant/compensator terms to zero is NOT inferred by this decomposition.
    """
    raw=cutoff_response(metric,cutoff,angular_max,frequency_points)
    local=ultraviolet_subtraction(metric,cutoff,normalization)
    gradient=raw['metric_gradients']-local['metric_gradients']
    n,q,r,h=metric.lapse,metric.radial_scale,metric.sphere_radius,metric.spacing
    rho=gradient[0]/(4*pi*q*r*r*h)
    px=-gradient[1]/(4*pi*n*r*r*h)
    pt=-gradient[2]/(8*pi*n*q*r*h)
    return {'raw_energy':raw['energy'],'subtracted_local_energy':local['energy'],
            'energy':raw['energy']-local['energy'],'metric_gradients':gradient,
            'rho':rho,'p_x':px,'p_perp':pt,'radial_null':rho+px,
            'trace':rho-px-2*pt,
            'normalization_log_derivative':local['remainder_normalization_log_derivative'],
            'cutoff':cutoff,'normalization':normalization,
            'local_coefficients':local['coefficients_in_Lorentzian_finite_term_basis']}


def spatial_cutoff_proxy(metric,cutoff=2.,angular_max=16):
    """Deliberately different prescription: fixed cutoff on the canonical H.

    It agrees at N=1 but not under lapse variation. Retained only as a control
    of the previously exposed spatial-versus-spacetime regulator distinction.
    """
    total=0.;root=np.sqrt(np.r_[metric.lapse,metric.lapse]);beta=np.kron(SIGMA3,np.eye(metric.points))
    for kappa in range(1,angular_max+1):
        h0=-1j*beta@euclidean_operator(metric,0.,kappa)
        h=root[:,None]*h0*root[None,:]
        values=np.linalg.eigvalsh(h);positive=values[values>0]
        total+=4*kappa*np.sum(cutoff/np.sqrt(pi)*np.exp(-(positive/cutoff)**2)
                              -positive*erfc(positive/cutoff))
    return float(total)


def cutoff_weyl_trace(metric,sigma,cutoff=2.,angular_max=16,frequency_points=48):
    """Direct raw cutoff trace for the continuum Weyl Ward comparison."""
    sigma=np.asarray(sigma,dtype=float)
    if sigma.shape!=(metric.points,):raise ValueError('Weyl weight must match the grid')
    extent=frequency_extent(metric,cutoff)
    nodes,weights=leggauss(frequency_points);frequencies=(nodes+1)*extent/2;weights=weights*extent/2
    total=0.;diagonal=np.r_[sigma,sigma]
    for kappa in range(1,angular_max+1):
        for omega,weight in zip(frequencies,weights):
            values,vectors=eigh(euclidean_operator(metric,omega,kappa),driver='evd',check_finite=False)
            measured=np.sum(abs(vectors)**2*diagonal[:,None],axis=0)
            total+=4*kappa*weight/pi*np.dot(measured,np.exp(-(values/cutoff)**2))
    return float(total)


def scalar_derivative(values,length):
    """Real periodic metric derivative; odd grids have no Nyquist mode to drop."""
    values=np.asarray(values,dtype=float)
    frequency=2*np.pi*np.fft.fftfreq(values.size,length/values.size)
    symbol=1j*frequency
    if values.size%2==0:symbol[values.size//2]=0.
    return np.fft.ifft(np.fft.fft(values)*symbol).real


def ward_summary(metric,response):
    dx=lambda values:scalar_derivative(values,metric.length)
    conservation=(dx(response['p_x'])+dx(metric.lapse)/metric.lapse*response['radial_null']
                  +2*dx(metric.sphere_radius)/metric.sphere_radius*(response['p_x']-response['p_perp']))
    return {'maximum_radial_conservation_residual':float(max(abs(conservation))),
            'rms_radial_conservation_residual':float(np.sqrt(np.mean(conservation**2))),
            'lapse_homogeneity_residual':response['lapse_homogeneity_residual']}


def local_trace_coefficient(metric):
    """a4/(16pi²) in the declared +--- energy-to-stress convention.

    R_E=-R_L while Box_E R_E=Box_L R_L for this static Wick rotation.
    Finite R² terms can alter the removable Laplacian part; none is inferred.
    """
    from .nsc_finite_terms import finite_response
    from .nsc_shape_response import StaticAxialMetric
    existing=StaticAxialMetric(metric.length,metric.lapse,metric.radial_scale,metric.sphere_radius)
    f=finite_response(existing)
    volume=metric.lapse*metric.radial_scale*metric.sphere_radius**2
    return (-18*f['densities'][2]+11*f['densities'][4]-12*f['densities'][5])/(360*16*pi*pi*volume)
