"""Short-domain propagation from derived high-energy boundary mode data.

The asymptotic expansion is used only outside/above the trapped probe cell;
the same Dirac operator is then integrated exactly through the probe supports.
No seed covariance or metric Cauchy data is assigned.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import make_interp_spline

from .nsc_pg_high_energy import branch_envelope,local_series,branch_riccati_residual,riccati_coefficients
from .nsc_pg_lll_preparation import LLLFlowAtlas,radial_packet
from .nsc_lorentzian import geometry
from .nsc_pg_retarded_packets import _parts
from .nsc_transmitting_dirac_domain import MODE_TO_CURRENT


class FastVacuumPacketProjector:
    """Reuse boundary/phase coefficients across energy quadrature nodes."""
    def __init__(self,mass,angular,preparation,*,order=16,minimum_energy=40.,rtol=2e-13,atol=2e-15):
        self.mass=mass;self.angular=angular;self.order=order;self.minimum_energy=minimum_energy
        self.rtol=rtol;self.atol=atol
        rho=np.linspace(1.,4.,193);coeff,_=riccati_coefficients(rho,mass,angular,order)
        self.at_one=coeff[0];self.at_three=coeff[128]
        U=angular/np.sqrt(1+rho*rho)+1j*mass
        primitive=make_interp_spline(rho,U[:,None]*coeff.conj(),k=7).antiderivative()
        self.phase_coeff=primitive(3.)-primitive(1.)
        self.beta_one=float(geometry(1.)[0]);self.beta_three=float(geometry(3.)[0])
        atlas=LLLFlowAtlas(preparation.horizon_rho,preparation.surface_gravity)
        self.clock_difference=float(atlas.coordinate(-1,np.array(3.))-atlas.coordinate(-1,np.array(1.)))

    def project(self,energy):
        if energy<self.minimum_energy:raise ValueError('frequency below the declared boundary approximation range')
        powers=(1/(2*energy))**np.arange(1,self.order+1)
        S1=self.at_one@powers;S3=self.at_three@powers;b1=self.beta_one;b3=self.beta_three
        def vector(b,S,branch):
            y=np.array([1.,1j*(1-b)*S]) if branch==1 else np.array([-1j*(1+b)*S.conjugate(),1.])
            current=(1-b)*abs(y[0])**2+(-1-b)*abs(y[1])**2
            return MODE_TO_CURRENT@y/np.sqrt(abs(current))
        partner=vector(b1,S1,1);outgoing=vector(b3,S3,1);incoming_one=vector(b1,S1,-1)
        q1=-1j*(1+b1)*S1.conjugate();j1=(1-b1)*abs(q1)**2-1-b1
        q3=-1j*(1+b3)*S3.conjugate()
        incoming_three=(MODE_TO_CURRENT@np.array([q3,1.]))/np.sqrt(abs(j1))*np.sqrt((1+b1)/(1+b3))
        incoming_three*=np.exp(1j*(energy*self.clock_difference+self.phase_coeff@powers))
        fields=np.column_stack((np.zeros(2),partner,incoming_one));output=np.zeros((8,3),complex)
        def derivative(rho,field):
            beta,bp,A=geometry(rho);lr=self.angular/np.sqrt(1+rho*rho)
            d0=(1j*(energy+lr)+.5*bp)*field[0]+1j*self.mass*field[1]
            d1=1j*self.mass*field[0]+(1j*(energy-lr)+.5*bp)*field[1]
            return np.array([(beta*d0-1j*d1)/A,(1j*d0+beta*d1)/A])
        for start,end,active in ((1.,0.,((0,'parent'),)),(0.,-1.,((2,'child'),(4,'child_bulk')))):
            initial=np.r_[fields.ravel(),np.zeros(6*len(active),complex)]
            def rhs(rho,y):
                field=y[:6].reshape(2,3)
                integ=[-float(radial_packet(np.array(rho),kind)[0])*field for _,kind in active]
                return np.r_[derivative(rho,field).ravel(),np.array(integ).ravel()]
            run=solve_ivp(rhs,(start,end),initial,method='DOP853',rtol=self.rtol,atol=self.atol)
            if not run.success:raise ArithmeticError(run.message)
            fields=run.y[:6,-1].reshape(2,3)
            for (index,_),value in zip(active,run.y[6:,-1].reshape(len(active),2,3)):output[index:index+2]=value
        initial=np.r_[np.column_stack((outgoing,np.zeros(2),incoming_three)).ravel(),np.zeros(6,complex)]
        def rhs_outer(rho,y):
            field=y[:6].reshape(2,3)
            return np.r_[derivative(rho,field).ravel(),(float(radial_packet(np.array(rho),'exterior_bulk')[0])*field).ravel()]
        run=solve_ivp(rhs_outer,(3.,4.),initial,method='DOP853',rtol=self.rtol,atol=self.atol)
        if not run.success:raise ArithmeticError(run.message)
        output[6:]=run.y[6:,-1].reshape(2,3)
        return output

    def project_many(self,energies):
        """Vectorized independent mode equations; no frequency coupling added."""
        from threadpoolctl import threadpool_limits
        E=np.asarray(energies,float);count=len(E)
        if E.ndim!=1 or np.any(E<self.minimum_energy):raise ValueError('energies below the numerical boundary range')
        powers=(1/(2*E[:,None]))**np.arange(1,self.order+1)
        S1=powers@self.at_one;S3=powers@self.at_three;b1=self.beta_one;b3=self.beta_three
        def vectors(b,S,branch):
            y=np.stack((np.ones_like(S),1j*(1-b)*S),axis=1) if branch==1 else np.stack((-1j*(1+b)*S.conj(),np.ones_like(S)),axis=1)
            current=(1-b)*abs(y[:,0])**2+(-1-b)*abs(y[:,1])**2
            return (y@MODE_TO_CURRENT.T)/np.sqrt(abs(current))[:,None]
        partner=vectors(b1,S1,1);outgoing=vectors(b3,S3,1);incoming_one=vectors(b1,S1,-1)
        q1=-1j*(1+b1)*S1.conj();q3=-1j*(1+b3)*S3.conj();j1=(1-b1)*abs(q1)**2-1-b1
        incoming_three=(np.stack((q3,np.ones_like(q3)),axis=1)@MODE_TO_CURRENT.T)/np.sqrt(abs(j1))[:,None]*np.sqrt((1+b1)/(1+b3))
        incoming_three*=np.exp(1j*(E*self.clock_difference+powers@self.phase_coeff))[:,None]
        fields=np.stack((np.zeros_like(partner),partner,incoming_one),axis=2);output=np.zeros((count,8,3),complex)
        def derivative(rho,field):
            beta,bp,A=geometry(rho);lr=self.angular/np.sqrt(1+rho*rho)
            d0=(1j*(E+lr)+.5*bp)[:,None]*field[:,0]+1j*self.mass*field[:,1]
            d1=1j*self.mass*field[:,0]+(1j*(E-lr)+.5*bp)[:,None]*field[:,1]
            return np.stack(((beta*d0-1j*d1)/A,(1j*d0+beta*d1)/A),axis=1)
        with threadpool_limits(limits=1):
            for start,end,active in ((1.,0.,((0,'parent'),)),(0.,-1.,((2,'child'),(4,'child_bulk')))):
                initial=np.r_[fields.ravel(),np.zeros(count*6*len(active),complex)]
                def rhs(rho,y):
                    field=y[:count*6].reshape(count,2,3)
                    integ=np.stack([-float(radial_packet(np.array(rho),kind)[0])*field for _,kind in active],axis=1)
                    return np.r_[derivative(rho,field).ravel(),integ.ravel()]
                run=solve_ivp(rhs,(start,end),initial,method='DOP853',rtol=self.rtol/2,atol=self.atol/2,t_eval=[end])
                if not run.success:raise ArithmeticError(run.message)
                fields=run.y[:count*6,-1].reshape(count,2,3)
                integ=run.y[count*6:,-1].reshape(count,len(active),2,3)
                for k,(index,_) in enumerate(active):output[:,index:index+2]=integ[:,k]
            initial=np.r_[np.stack((outgoing,np.zeros_like(outgoing),incoming_three),axis=2).ravel(),np.zeros(count*6,complex)]
            def rhs_outer(rho,y):
                field=y[:count*6].reshape(count,2,3)
                return np.r_[derivative(rho,field).ravel(),(float(radial_packet(np.array(rho),'exterior_bulk')[0])*field).ravel()]
            run=solve_ivp(rhs_outer,(3.,4.),initial,method='DOP853',rtol=self.rtol/2,atol=self.atol/2,t_eval=[4.])
            if not run.success:raise ArithmeticError(run.message)
            output[:,6:]=run.y[count*6:,-1].reshape(count,2,3)
        return output


def _normalized_boundary(rho,energy,mass,angular,branch,order):
    S,_=local_series([rho],energy,mass,angular,order);b=float(geometry(rho)[0])
    if branch==1:vector=np.array([1.,1j*(1-b)*S[0]],complex)
    else:vector=np.array([-1j*(1+b)*S[0].conjugate(),1.],complex)
    current=(1-b)*abs(vector[0])**2+(-1-b)*abs(vector[1])**2
    return MODE_TO_CURRENT@vector/np.sqrt(abs(current))


def fast_vacuum_packet_modes(energy,mass,angular,preparation,*,order=16,minimum_energy=40.,rtol=2e-13,atol=2e-15):
    """Return the three tail source columns and their boundary residuals.

    Warm source and beyond-all-orders reflection corrections require the
    acceptance record's explicit bound/control; this routine is not full C1b.
    """
    if energy<minimum_energy:raise ValueError('frequency is below this numerical boundary approximation range')
    atlas=LLLFlowAtlas(preparation.horizon_rho,preparation.surface_gravity)
    bridge=branch_envelope(np.linspace(1.,4.,193),energy,mass,angular,-1,order=order,reference=1.)
    incoming_one=bridge.envelope[0]*np.exp(1j*energy*atlas.coordinate(-1,np.array(1.)))
    incoming_three=bridge.envelope[128]*np.exp(1j*energy*atlas.coordinate(-1,np.array(3.)))
    partner=_normalized_boundary(1.,energy,mass,angular,1,order)
    outgoing=_normalized_boundary(3.,energy,mass,angular,1,order)
    fields=np.column_stack((np.zeros(2),partner,incoming_one));output=np.zeros((8,3),complex)
    for start,end,active in ((1.,0.,((0,'parent'),)),(0.,-1.,((2,'child'),(4,'child_bulk')))):
        initial=np.r_[fields.ravel(),np.zeros(6*len(active),complex)]
        def rhs(rho,state):
            A,_=_parts(rho,energy,mass,angular);field=state[:6].reshape(2,3)
            integrals=[-float(radial_packet(np.array(rho),kind)[0])*field for _,kind in active]
            return np.r_[(A@field).ravel(),np.array(integrals).ravel()]
        run=solve_ivp(rhs,(start,end),initial,method='DOP853',rtol=rtol,atol=atol)
        if not run.success:raise ArithmeticError(run.message)
        fields=run.y[:6,-1].reshape(2,3)
        for (index,_),value in zip(active,run.y[6:,-1].reshape(len(active),2,3)):output[index:index+2]=value
    initial=np.r_[np.column_stack((outgoing,np.zeros(2),incoming_three)).ravel(),np.zeros(6,complex)]
    def outer_rhs(rho,state):
        A,_=_parts(rho,energy,mass,angular);field=state[:6].reshape(2,3)
        return np.r_[(A@field).ravel(),(float(radial_packet(np.array(rho),'exterior_bulk')[0])*field).ravel()]
    run=solve_ivp(outer_rhs,(3.,4.),initial,method='DOP853',rtol=rtol,atol=atol)
    if not run.success:raise ArithmeticError(run.message)
    output[6:]=run.y[6:,-1].reshape(2,3)
    residuals={'incoming_bridge_equation':bridge.equation_residual,'incoming_bridge_current':bridge.current_residual,
               'outgoing_inner_boundary_equation':float(branch_riccati_residual([1.],energy,mass,angular,1,order)[0]),
               'outgoing_exterior_boundary_equation':float(branch_riccati_residual([3.],energy,mass,angular,1,order)[0])}
    return output,residuals
