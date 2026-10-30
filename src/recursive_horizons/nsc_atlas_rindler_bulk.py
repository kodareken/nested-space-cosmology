"""Direct finite-order Ricci check of the published Rindler reconstruction.

An exact decaying shear solution verifies the local expansion through order
four, including the required affine velocity correction. This verifier is
independent of the spectral constraint calculation on the retained vortex.
"""
from functools import lru_cache
import numpy as np
import sympy as s

ORDER=4
DIM=5


def ricci_series(g,dg,ddg):
    """Covariant Ricci tensor from metric jets; coefficients of epsilon^0..4."""
    P=len(g);B=np.zeros_like(g);exact=g.dtype==object
    half=s.Rational(1,2) if exact else .5
    B[0]=np.array(s.Matrix(g[0]).inv().tolist(),dtype=object) if exact else np.linalg.inv(g[0])
    for n in range(1,P):
        B[n]=-B[0]@sum((g[k]@B[n-k] for k in range(1,n+1)),np.zeros((DIM,DIM),dtype=g.dtype))
    C=np.zeros((P,DIM,DIM,DIM),dtype=g.dtype);dC=np.zeros((DIM,P,DIM,DIM,DIM),dtype=g.dtype)
    for a in range(DIM):
        for b in range(DIM):
            for c in range(DIM):
                C[:,a,b,c]=dg[b,:,a,c]+dg[c,:,a,b]-dg[a,:,b,c]
                for e in range(DIM):dC[e,:,a,b,c]=ddg[e,b,:,a,c]+ddg[e,c,:,a,b]-ddg[e,a,:,b,c]
    dB=np.zeros_like(dg)
    for e in range(DIM):
        for n in range(P):
            for i in range(n+1):
                for j in range(n-i+1):dB[e,n]-=B[i]@dg[e,j]@B[n-i-j]
    G=np.zeros_like(C);dG=np.zeros_like(dC)
    for n in range(P):
        for i in range(n+1):
            G[n]+=half*np.tensordot(B[i],C[n-i],axes=(1,0))
            for e in range(DIM):
                dG[e,n]+=half*(np.tensordot(dB[e,i],C[n-i],axes=(1,0))+np.tensordot(B[i],dC[e,n-i],axes=(1,0)))
    R=np.zeros_like(g)
    for n in range(P):
        for a in range(DIM):
            for b in range(DIM):
                R[n,a,b]=sum(dG[c,n,c,a,b]-dG[b,n,c,a,c] for c in range(DIM))
        for i in range(n+1):
            R[n]+=np.tensordot(np.trace(G[i],axis1=0,axis2=1),G[n-i],axes=(0,0))-np.tensordot(G[i],G[n-i],axes=([0,2],[2,0])).T
    return R


@lru_cache(maxsize=2)
def shear_metric_evaluator(corrected=True):
    """Compile metric jets; no imported tensor package or field equation residual."""
    eps=s.symbols('eps');T,r,X,Y,Z,nu,phase=s.symbols('T r X Y Z nu phase',real=True)
    coords=(T,r,X,Y,Z);weights=(2,0,1,1,1);space=(X,Y,Z)
    A=s.exp(-nu*T);u=[A*s.sin(Y+phase),s.Integer(0),s.Integer(0)]
    w=[A*A*q/3 if corrected else s.Integer(0) for q in space]
    v=[eps*u[i]+eps**3*w[i] for i in range(3)]
    lap=[eps**2*sum(s.diff(vi,q,2) for q in space) for vi in v]
    sigma=s.Matrix(3,3,lambda i,j:eps*(s.diff(v[j],space[i])+s.diff(v[i],space[j])))
    omega=s.Matrix(3,3,lambda i,j:eps*(s.diff(v[j],space[i])-s.diff(v[i],space[j])))
    v2=sum(vi*vi for vi in v);vlap=sum(v[i]*lap[i] for i in range(3))
    sig2=sum(q*q for q in sigma);omcontract=s.trace(omega*omega)
    gradv2=[eps*s.diff(v2,q) for q in space]
    lapsigma=s.Matrix(3,3,lambda i,j:eps**2*sum(s.diff(sigma[i,j],q,2) for q in space))
    g=s.zeros(5);g[0,0]=-r+(1-r/nu)*v2;g[0,1]=g[1,0]=1+v2/(2*nu)
    for i in range(3):
        g[0,i+2]=g[i+2,0]=-(1-r/nu)*v[i]+(r-nu)/(2*nu)*(2*v2*v[i]/nu-(r+nu)*lap[i])
        g[1,i+2]=g[i+2,1]=-v[i]/nu
        for j in range(3):
            sijomega=sum((sigma[k,i]*omega[j,k]+sigma[k,j]*omega[i,k])/2 for k in range(3))
            term=v[i]*v[j]*v2/nu**2-(sigma*sigma)[i,j]/2+(r-5*nu)/(4*nu)*(omega*omega)[i,j]
            term+=sijomega-(r+nu)/(2*nu)*(v[i]*lap[j]+v[j]*lap[i])+(r+5*nu)*lapsigma[i,j]/4
            term-=(v[i]*gradv2[j]+v[j]*gradv2[i])/(2*nu)+sigma[i,j]*v2/(2*nu)
            g[i+2,j+2]=int(i==j)+(1-r/nu)*v[i]*v[j]/nu+(1-r/nu)*term
    F=9*v2*v2/(8*nu)-2*nu*vlap-nu*sig2/2
    g[0,0]+=-(r-nu)**3*omcontract/(8*nu**2)+(r-nu)**2*(8*vlap+sig2)/(8*nu)-(r-nu)*F/nu
    coefficients=[g.applyfunc(lambda expr:s.expand(expr).coeff(eps,n)) for n in range(ORDER+1)]
    specs=[];expressions=[]
    for n in range(ORDER+1):
        for a in range(DIM):
            for b in range(a,DIM):
                f=coefficients[n][a,b]
                if f==0:continue
                specs.append(('g',n,a,b));expressions.append(f)
                for i in range(DIM):
                    if n+weights[i]>ORDER:continue
                    fi=s.diff(f,coords[i])
                    if fi!=0:specs.append(('dg',i,n+weights[i],a,b));expressions.append(fi)
                    for j in range(i,DIM):
                        if n+weights[i]+weights[j]>ORDER:continue
                        fij=s.diff(fi,coords[j])
                        if fij!=0:specs.append(('ddg',i,j,n+weights[i]+weights[j],a,b));expressions.append(fij)
    evaluator=s.lambdify((T,r,X,Y,Z,nu,phase),expressions,'numpy',cse=True)
    def evaluate(T0,r0,X0,Y0,Z0,nu0=1.,phase0=.4,*,exact=False):
        dtype=object if exact else float
        arrays=dict(g=np.zeros((5,5,5),dtype=dtype),dg=np.zeros((5,5,5,5),dtype=dtype),ddg=np.zeros((5,5,5,5,5),dtype=dtype))
        values=[e.subs(dict(zip((T,r,X,Y,Z,nu,phase),(T0,r0,X0,Y0,Z0,nu0,phase0))),simultaneous=True) for e in expressions] if exact else evaluator(T0,r0,X0,Y0,Z0,nu0,phase0)
        for spec,val in zip(specs,values):
            key,*indices=spec;arr=arrays[key];arr[tuple(indices)]=val
            swap=indices.copy();swap[-1],swap[-2]=swap[-2],swap[-1];arr[tuple(swap)]=val
            if key=='ddg':
                swap0=indices.copy();swap0[0],swap0[1]=swap0[1],swap0[0];arr[tuple(swap0)]=val
                swap0[-1],swap0[-2]=swap0[-2],swap0[-1];arr[tuple(swap0)]=val
        return arrays
    return evaluate


def exact_radial_check():
    r=s.symbols('r',real=True)
    jets=shear_metric_evaluator(True)(0,r,s.Rational(1,5),0,-s.Rational(1,7),1,s.pi/6,exact=True)
    R=ricci_series(**jets)
    nonzero=[(index,str(s.factor(value))) for index,value in np.ndenumerate(R) if s.factor(value)!=0]
    return dict(all_coefficients_zero=not nonzero,checked_components=R.size,
        radial_coordinate_symbolic=True,nonzero=nonzero,
        scope='Exact rational/trigonometric shear-flow jet at T=0, Y=0, phase=pi/6, X=1/5, Z=-1/7, nu=1; all radial coordinates. General profiles use the attributed reconstruction.')


def checks():
    rows=[]
    for nu in (1.,.2):
        for fraction in (0.,.5,1.):
            args=(.1,fraction*nu,.2,.3,-.15,nu,.4)
            good=shear_metric_evaluator(True)(*args);bad=shear_metric_evaluator(False)(*args)
            R=ricci_series(**good);Rb=ricci_series(**bad)
            rows.append(dict(nu=nu,radial_fraction=fraction,
                maximum_Ricci_coefficient_by_order=np.max(np.abs(R),axis=(1,2)).tolist(),
                omitted_expansion_order4_residual=float(np.max(np.abs(Rb[4]))),
                background_determinant=float(np.linalg.det(good['g'][0]))))
    return rows
