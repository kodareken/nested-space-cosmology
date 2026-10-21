"""One-sided charged-wall scalar reconstruction and finite-time energy balance.

These are necessary local identities for the declared Einstein/scalar/wall
action, not sufficient conditions for a global bulk solution or proof of
the existence of the excited-parent evolution.
"""
from dataclasses import dataclass, field
import math

from .nsc_atlas_scalar_preparation import ScalarPreparation
from .nsc_atlas_flux_conversion import ScalarFlux


@dataclass(frozen=True)
class WallGeometry:
    """Child-side orthonormal geometry; dots use wall proper time.

    n points child to parent. Extrinsic mixed eigenvalues are A,B,B.
    B=n(R)/R; A=n.acceleration. Normal derivatives are essential inputs.
    """
    radius: float
    radius_dot: float
    radius_ddot: float
    A: float
    B: float
    A_dot: float
    B_dot: float

    def __post_init__(self):
        if not all(math.isfinite(v) for v in vars(self).values()) or self.radius<=0:
            raise ValueError('finite one-sided geometry and positive radius required')

    def projections(self, kappa):
        if not math.isfinite(kappa) or kappa<=0:
            raise ValueError('positive finite Einstein coupling required')
        R,v,d,A,B,Ad,Bd=(self.radius,self.radius_dot,self.radius_ddot,
                         self.A,self.B,self.A_dot,self.B_dot)
        H=v/R
        Gnn=-2*d/R-(1+v*v)/(R*R)+2*A*B+B*B
        Gnu=2*(H*(A-B)-Bd)
        mass=R*(1+v*v-B*B*R*R)/2
        mass_dot=v*(1+v*v)/2-1.5*B*B*R*R*v+R*v*d-B*R**3*Bd
        return dict(radial_pressure=Gnn/kappa,inward_flux=Gnu/kappa,
                    G_nn=Gnn,G_nu=Gnu,Misner_Sharp_mass=mass,
                    geometric_mass_dot=mass_dot,
                    flux_contribution=kappa*R**3*B*(Gnu/kappa)/2,
                    pressure_work=-kappa*R*R*v*(Gnn/kappa)/2)


@dataclass(frozen=True)
class ScalarWallObserver:
    prepared: ScalarPreparation = field(default_factory=ScalarPreparation)

    @property
    def force_scale(self):
        S=self.prepared.scalar.mu*self.prepared.charge
        if S==0:
            raise ValueError('mu*e=0: this wall pressure cannot identify the scalar')
        return S

    @property
    def ground_pressure(self):
        m=self.prepared
        return m.charge*(m.scalar.epsilon*m.child_flux+m.charge/2)

    @property
    def tidal_quadratic_coefficient(self):
        m=self.prepared.scalar
        return m.kappa*m.omega**2/(2*self.force_scale**2)

    def reconstruct(self, geometry, *, measured_Ruu=None):
        m=self.prepared
        S=self.force_scale
        g=geometry.projections(m.scalar.kappa)
        pressure=m.tension*(geometry.A+2*geometry.B-1.5*m.k)
        pressure_dot=m.tension*(geometry.A_dot+2*geometry.B_dot)
        displacement=(pressure-self.ground_pressure)/S
        phi=m.scalar.minimum(m.child_flux)+displacement
        Pi=pressure_dot/S
        U=m.child_vacuum+m.scalar.omega**2*displacement**2/2
        Q=2*(g['radial_pressure']+U)-Pi*Pi
        F=g['inward_flux']
        signed_normal=F/Pi if Pi!=0 else None
        predicted_Ruu=m.scalar.kappa*(Pi*Pi-U)
        vacuum=m.scalar.kappa*m.child_vacuum
        R,v=geometry.radius,geometry.radius_dot
        out=dict(pressure=pressure,pressure_dot=pressure_dot,
            scalar_value=phi,child_displacement=displacement,scalar_time_derivative=Pi,
            normal_square=Q,normal_magnitude=math.sqrt(Q) if Q>=0 else None,
            signed_normal=signed_normal,
            normal_sign_identified=Pi!=0,
            normal_flux_sensitivity=1/abs(Pi) if Pi!=0 else None,
            scalar_stress_residual=F*F-Pi*Pi*Q,
            potential=U,predicted_Ruu=predicted_Ruu,
            predicted_tidal_corrected=predicted_Ruu-m.scalar.kappa*pressure_dot**2/S**2,
            vacuum_subtracted_mass=g['Misner_Sharp_mass']-vacuum*R**3/6,
            vacuum_subtracted_mass_dot=g['geometric_mass_dot']-vacuum*R*R*v/2,
            excitation_pressure_work=-m.scalar.kappa*R*R*v*(g['radial_pressure']+m.child_vacuum)/2,
            **g)
        if measured_Ruu is not None:
            if not math.isfinite(measured_Ruu):raise ValueError('finite independent Ricci projection required')
            out.update(measured_Ruu=measured_Ruu,tidal_residual=measured_Ruu-predicted_Ruu,
                measured_tidal_corrected=measured_Ruu-m.scalar.kappa*pressure_dot**2/S**2,
                inferred_vacuum=m.scalar.kappa*(Pi*Pi-m.scalar.omega**2*displacement**2/2)-measured_Ruu)
        return out

    def calibrate_tidal_pair(self, first, second):
        """Infer Pc and Lambda from two samples with other couplings fixed.

        Samples supply pressure, its proper-time derivative and an independent
        Ruu. The model's existing parent flux and bare density are not used to
        fit the pair. Its fixed bare term is then tested, not silently refitted.
        """
        C=self.tidal_quadratic_coefficient
        kap=self.prepared.scalar.kappa
        S=self.force_scale
        p1,p2=first['pressure'],second['pressure']
        if not all(math.isfinite(s[key]) for s in (first,second)
                   for key in ('pressure','pressure_dot','Ruu')) or p1==p2:
            raise ValueError('finite samples with distinct pressures required')
        y1=first['Ruu']-kap*first['pressure_dot']**2/S**2
        y2=second['Ruu']-kap*second['pressure_dot']**2/S**2
        z1,z2=y1+C*p1*p1,y2+C*p2*p2
        slope=(z2-z1)/(p2-p1)
        Pc=slope/(2*C)
        vacuum=-(z1-slope*p1)-C*Pc*Pc
        eps=self.prepared.scalar.epsilon
        e=self.prepared.charge
        qc=(Pc/e-e/2)/eps if eps!=0 else None
        fixed_bare_vacuum=kap*(self.prepared.scalar.bare_density+eps*qc*qc/2) if qc is not None else None
        return dict(ground_pressure=Pc,vacuum=vacuum,inferred_child_flux=qc,
                    fixed_action_vacuum_residual=vacuum-fixed_bare_vacuum if qc is not None else None,
                    calibration_pressure_separation=abs(p2-p1))

    def predict_tidal(self, pressure, pressure_dot, calibration):
        return (self.prepared.scalar.kappa*pressure_dot**2/self.force_scale**2
                -calibration['vacuum']
                -self.tidal_quadratic_coefficient*(pressure-calibration['ground_pressure'])**2)

    def divided_difference(self, samples):
        if len(samples)!=3:raise ValueError('exactly three independent samples required')
        p=[s['pressure'] for s in samples]
        if len(set(p))!=3:raise ValueError('distinct pressures required')
        y=[s['Ruu']-self.prepared.scalar.kappa*s['pressure_dot']**2/self.force_scale**2 for s in samples]
        return sum(y[i]/((p[i]-p[(i+1)%3])*(p[i]-p[(i+2)%3])) for i in range(3))


def regular_wall_geometry(prepared, solution, eta):
    """Geometry-only jet of the previously solved O(3,1) wall continuation."""
    if not 0<=eta<=1:raise ValueError('reference interval is 0<=eta<=1')
    R0=float(solution.y[0,-1])
    B=float(solution.y[1,-1]/R0)
    return WallGeometry(R0*math.cosh(eta),math.sinh(eta),math.cosh(eta)/R0,B,B,0.,0.)


def regular_bulk_Ruu(solution):
    """Independent one-sided Ricci from the metric profile and its derivatives."""
    a,ap=map(float,solution.y[:2,-1])
    app=float(solution.sol(1.,1)[1]/solution.p[0])
    return (a*app+2*ap*ap-2)/(a*a)


def infer_mixing_and_charge(tidal_C, vacuum_second_difference, *, kappa=1.):
    """Potential-level identification from local response and three ground vacua.

    Requires common couplings and equal oriented charge steps. It does not
    establish that the three regions have a compatible shared spacetime.
    """
    if (not all(map(math.isfinite,(tidal_C,vacuum_second_difference,kappa)))
            or min(tidal_C,kappa)<=0 or vacuum_second_difference<0):
        raise ValueError('positive C,kappa and nonnegative vacuum second difference required')
    x=vacuum_second_difference/kappa
    z=kappa/(2*tidal_C)
    return dict(epsilon=x/(x+z),charge_magnitude=math.sqrt(x+z),charge_squared=x+z)


def equivalent_one_wall_example():
    """Two declared coupling choices with identical two-sector reduced dynamics.

    The third charge sector distinguishes them. Absolute scalar values are
    translated; this equivalence does not include independent charge probes,
    additional matter couplings or quantum preparation probabilities.
    """
    first=ScalarPreparation()
    second=ScalarPreparation(scalar=ScalarFlux(epsilon=25/28),
        parent_flux=3*math.sqrt(7)/10,charge=math.sqrt(7)/2,tension=.2)
    return first,second,3*math.sqrt(3)/100
