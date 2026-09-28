"""Vectorized independent projections of already computed PG mode fields."""
import numpy as np
from scipy.integrate import solve_ivp
from threadpoolctl import threadpool_limits

from .nsc_pg_retarded_packets import _parts
from .nsc_pg_lll_preparation import radial_packet


def project_mode_family(energies,mass,angular,modes_at_zero,modes_at_three,*,rtol=1e-13,atol=1e-15):
    E=np.asarray(energies,float);zero=np.asarray(modes_at_zero,complex);three=np.asarray(modes_at_three,complex)
    count=len(E)
    if E.ndim!=1 or np.any(E<=0) or zero.shape!=(count,2,3) or three.shape!=zero.shape:
        raise ValueError('one physical positive-energy mode family required')
    if not all(np.isfinite(a).all() for a in (E,zero,three)):raise ValueError('finite mode fields required')
    answer=np.zeros((count,8,3),complex)
    with threadpool_limits(limits=1):
        for start,end,initial,active in ((0.,1.,zero,((0,'parent'),)),
                                        (0.,-1.,zero,((2,'child'),(4,'child_bulk'))),
                                        (3.,4.,three,((6,'exterior_bulk'),))):
            orientation=np.sign(end-start)
            def rhs(rho,state):
                field=state[:count*6].reshape(count,2,3)
                A,B=_parts(rho,0.,mass,angular)
                derivative=np.einsum('ij,njk->nik',A,field)+E[:,None,None]*np.einsum('ij,njk->nik',B,field)
                integral=np.stack([orientation*float(radial_packet(np.array(rho),kind)[0])*field for _,kind in active],axis=1)
                return np.r_[derivative.ravel(),integral.ravel()]
            state=np.r_[initial.ravel(),np.zeros(count*6*len(active),complex)]
            run=solve_ivp(rhs,(start,end),state,method='DOP853',rtol=rtol,atol=atol,t_eval=[end])
            if not run.success:raise ArithmeticError(run.message)
            integrals=run.y[count*6:,-1].reshape(count,len(active),2,3)
            for i,(index,_) in enumerate(active):answer[:,index:index+2]=integrals[:,i]
    return answer
