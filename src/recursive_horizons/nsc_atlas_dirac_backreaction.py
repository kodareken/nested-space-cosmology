"""First-order spherical metric from a supplied conserved finite mean source.

Transport the nonvacuum Misner-Sharp excess and outgoing-EF lapse, use the
existing signed wall-mass equation, then differentiate the null-ray equation
analytically to obtain the proper-wall radar clock. All vacuum terms and
background geometry are those of the current quiet reference.

This computes an explicitly supplied finite state difference only. It is
the mean response required for the later mixed metric/source/tadpole solve;
it is not the full vacuum polarization or Psi/Phi.
"""
from __future__ import annotations

import time

import numpy as np
from scipy.integrate import solve_ivp

from .nsc_atlas_cap_wall_response import reference_coefficients,wall_response_coefficients


def mass_one_form(fields,stress,*,gravity_coupling=1.):
    """d(nonvacuum mass); subtract Lambda*r^3/6 at the ACTUAL areal radius."""
    r,C = fields['radius'],fields['conformal_factor']
    ru,rv = (fields['radius_t']-fields['radius_x'])/2,(
        fields['radius_t']+fields['radius_x'])/2
    weight = gravity_coupling*r*r/C
    return dict(u=weight*(stress['uv']*ru-stress['uu']*rv),
                v=weight*(stress['uv']*rv-stress['vv']*ru))


def ef_lapse_v(fields,stress,*,gravity_coupling=1.):
    rv = (fields['radius_t']+fields['radius_x'])/2
    return gravity_coupling*fields['radius']*stress['vv']/(2*rv)


class MeanBackreaction:
    """Physical nested wall, or smooth-cap control with a supplied worldtube.

    A smooth cap has no membrane. Its optional shared worldtube is explicitly
    driven along the already computed nested areal-radius history. It then
    compares geometry at the same observer motion; it is not an alternative
    material-wall dynamics or independently calibrated cosmological model.
    """

    def __init__(self,source,end,*,gravity_coupling=1.,initial_masses=(0.,0.),
                 initial_displacement=0.,shared_worldtube=None,
                 wall_rtol=2e-8,wall_atol=2e-10,wall_cpu_limit=None):
        self.source,self.chart = source,source.chart
        self.gravity_coupling = gravity_coupling
        self.start,self.end = float(source.start),float(end)
        self.shared_worldtube = shared_worldtube
        if not 0<=self.start<self.end<=self.chart.duration:
            raise ValueError('wall history must lie in the declared proper-clock interval')
        if self.start>0 and (any(initial_masses) or initial_displacement):
            raise ValueError('a quiet prefix requires zero corner perturbations')
        if shared_worldtube is not None and source.region!='cap':
            raise ValueError('the smooth control must use the same continued cap')
        self.initial = np.array([*initial_masses,initial_displacement],dtype=float)
        self.wall_rhs_calls = 0
        deadline = None if wall_cpu_limit is None else time.process_time()+wall_cpu_limit
        def rhs(t,y):
            if deadline is not None and time.process_time()>deadline:
                raise RuntimeError('mean wall history exhausted its setup CPU allowance')
            self.wall_rhs_calls += 1
            g = self._geometry(t)
            r = float(g['x'])
            mean = self.source.evaluate(t,t)['orthonormal']
            j,pr = float(mean['outward_energy_flux']),float(mean['radial_pressure'])
            mu_dot = -gravity_coupling*r*r*(g['beta_plus']*j+g['v']*pr)/2
            nu_dot = -gravity_coupling*r*r*(g['beta_minus']*j+g['v']*pr)/2
            if self.shared_worldtube is None:
                rate = self._response(g,y[2],y[0],y[1])['displacement_rate']
            else:
                rate = self.shared_worldtube.wall_state(t)['displacement_rate']
                mu_dot = nu_dot
            return [mu_dot,nu_dot,rate]
        solution = solve_ivp(rhs,(self.start,self.end),self.initial,method='DOP853',
            rtol=wall_rtol,atol=wall_atol,max_step=self.chart.duration/50,dense_output=True)
        if not solution.success:
            raise RuntimeError(solution.message)
        self.wall_solution = solution.sol
        self.wall_nfev = solution.nfev
        self._quadrature_rules = {}
        self.source_evaluations = 0
        self._cpu_deadline = None

    def _geometry(self,t):
        if hasattr(self.chart,'wall_geometry'):
            return self.chart.wall_geometry(t)
        return reference_coefficients(float(self.chart.radius_at_clock(t)))

    def _response(self,g,rho,mu,nu):
        if not hasattr(self.chart,'wall_geometry'):
            return wall_response_coefficients(float(g['x']),rho,mu,nu)
        r,v,bp=g['x'],g['v'],g['beta_plus']
        d=(mu-nu)/(g['k_wall']*r*r)
        return dict(displacement_rate=(g['radial_acceleration']*rho+bp*d+mu/r)/v,
            delta_beta_minus=g['beta_minus_x']*rho+d,
            delta_beta_plus=g['beta_plus_x']*rho+d)

    def wall_state(self,t):
        t = float(t)
        if not 0<=t<=self.end:
            raise ValueError('wall query exits the computed history')
        y = self.initial if t<=self.start else self.wall_solution(t)
        g = self._geometry(t)
        R = float(g['x'])
        if self.shared_worldtube is None:
            response = self._response(g,y[2],y[0],y[1])
            rho,rate = float(y[2]),float(response['displacement_rate'])
            beta_changes = (float(response['delta_beta_minus']),
                            float(response['delta_beta_plus']))
        else:
            supplied = self.shared_worldtube.wall_state(t)
            rho,rate = supplied['displacement'],supplied['displacement_rate']
            fr = g.get('cap_f_prime',-2*R)
            db = (g['v']*rate+fr*rho/2-float(y[1])/R)/g['beta_minus']
            beta_changes = (float(db),float(db))
        return dict(radius=R,displacement=rho,displacement_rate=rate,
            exterior_mass=float(y[0]),cap_mass=float(y[1]),delta_beta=beta_changes)

    def _source(self,U,V,*,region=None):
        self._guard()
        self.source_evaluations += int(np.broadcast_arrays(U,V)[0].size)
        return self.source.evaluate(U,V,region=region)

    def _guard(self):
        if self._cpu_deadline is not None and time.process_time()>self._cpu_deadline:
            raise RuntimeError('mean clock exhausted its case CPU allowance')

    def _integrate(self,lower,upper,callback,order):
        if lower==upper:
            return np.zeros(2)
        sign = 1 if upper>lower else -1
        a,b = min(lower,upper),max(lower,upper)
        knots = np.asarray(self.source.breakpoints)
        inside = knots[(knots>a)&(knots<b)]
        if a<self.source.start<b:
            inside = np.r_[inside,self.source.start]
        bounds = np.r_[a,np.unique(inside),b]
        if order not in self._quadrature_rules:
            self._quadrature_rules[order] = np.polynomial.legendre.leggauss(order)
        nodes,weights = self._quadrature_rules[order]
        mid,half = (bounds[:-1]+bounds[1:])/2,np.diff(bounds)/2
        query = (mid[:,None]+half[:,None]*nodes).ravel()
        quadrature_weights = (half[:,None]*weights).ravel()
        result = callback(query)
        return sign*np.einsum('p,...p->...',quadrature_weights,result)

    def ef_fields(self,U,V,*,order=8,region=None):
        """Nonvacuum mass and delta log(a), each integrated on one null side."""
        self._guard()
        U,V = float(U),float(V)
        wall = self.wall_state(U)
        fields = self.source.points(U,V,region=region)
        actual_cap = bool(fields['is_cap'])
        side = 0 if actual_cap else 1
        g = self._geometry(U)
        beta = float(g['beta_minus'] if actual_cap else g['beta_plus'])
        mass = wall['cap_mass'] if actual_cap else wall['exterior_mass']
        psi = -(wall['delta_beta'][side]+wall['displacement_rate'])/(beta+g['v'])
        def integrand(v):
            row = self._source(U,v,region=region)
            dm = mass_one_form(row['fields'],row['null'],gravity_coupling=self.gravity_coupling)
            da = ef_lapse_v(row['fields'],row['null'],gravity_coupling=self.gravity_coupling)
            return np.array([dm['v'],da])
        addition = self._integrate(U,V,integrand,order)
        return fields,float(mass+addition[0]),float(psi+addition[1])

    def clock_response(self,U,V,*,order=8,ray_rtol=2e-8,ray_atol=2e-10,cpu_limit=None):
        """Analytic radial-null sensitivity and fixed proper-wall radar clock."""
        U,V = float(U),float(V)
        started = time.process_time()
        self._cpu_deadline = None if cpu_limit is None else started+cpu_limit
        region = ('cap' if self.source.region=='cap' or U>=V else 'exterior')
        wall = self.wall_state(V)
        fields,mass,psi = self.ef_fields(V,V,order=order,region=region)
        r,C = float(fields['radius']),float(fields['conformal_factor'])
        ru,rv = [(float(fields['radius_t'])+s*float(fields['radius_x']))/2 for s in (-1,1)]
        a = C/(2*rv)
        fp = -4*float(fields['radius_uv'])/C
        su = -a*fp*wall['displacement']/2+a*mass/r-a*float(fields['static_f'])*psi/2
        initial = [wall['displacement'],wall['displacement_rate']-su]
        wall_q = float(psi+initial[1]/rv)
        def rhs(u,y):
            fields,m,psi = self.ef_fields(u,V,order=order,region=region)
            r,C,f = (float(fields[k]) for k in ('radius','conformal_factor','static_f'))
            rv = (float(fields['radius_t'])+float(fields['radius_x']))/2
            a = C/(2*rv)
            fp = -4*float(fields['radius_uv'])/C
            fpp = -4*float(fields['log_C_uv'])/C
            row = self._source(u,V,region=region)
            mv = float(mass_one_form(fields,row['null'],gravity_coupling=self.gravity_coupling)['v'])
            psiv = float(ef_lapse_v(fields,row['null'],gravity_coupling=self.gravity_coupling))
            b = -a*fp/2
            g = a*m/r-a*f*psi/2
            bv = -a*fpp*rv/2
            gv = a*(mv/r-m*rv/r**2)-a*(fp*rv*psi+f*psiv)/2
            return [b*y[0]+g,b*y[1]+bv*y[0]+gv]
        if U==V:
            endpoint = np.array(initial)
            calls = 0
        else:
            solution = solve_ivp(rhs,(V,U),initial,method='DOP853',rtol=ray_rtol,
                atol=ray_atol,max_step=self.chart.duration/40)
            if not solution.success:
                raise RuntimeError(solution.message)
            endpoint,calls = solution.y[:,-1],solution.nfev
        fields,m,psi = self.ef_fields(U,V,order=order,region=region)
        rv = (float(fields['radius_t'])+float(fields['radius_x']))/2
        q = psi+endpoint[1]/rv
        def other_path(u):
            row = self._source(u,V,region=region)
            dm = mass_one_form(row['fields'],row['null'],gravity_coupling=self.gravity_coupling)
            return np.array([dm['u'],np.zeros_like(dm['u'])])
        mass_at_V = wall['cap_mass'] if bool(fields['is_cap']) else wall['exterior_mass']
        alternate = mass_at_V+self._integrate(V,U,other_path,order)[0]
        return dict(radius_change=float(endpoint[0]),log_radar_C_change=float(q),
            frequency_fractional_change=float(-q/2),
            received_over_wall_frequency_change=float(-q/(2*np.sqrt(fields['conformal_factor']))),
            background_received_over_wall_frequency=float(1/np.sqrt(fields['conformal_factor'])),
            mass=float(m),mass_other_path=float(alternate),mass_path_gap=float(m-alternate),
            wall_clock_log_C_at_corner=wall_q,
            ray_rhs_calls=int(calls),source_points=self.source_evaluations,
            process_CPU_seconds=time.process_time()-started)
