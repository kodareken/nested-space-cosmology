"""Finite-radius canonical massless reference stress in the conformal child.

Development owner: full S2 channels as massive 2D Dirac fields. This is
not the finite-Lambda compact complement or a parent-selected state.
"""
from functools import lru_cache
from math import atan, ceil, pi, sqrt
import ctypes
import hashlib
from pathlib import Path
import subprocess
import tempfile

import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss
import sympy as sp


@lru_cache(maxsize=1)
def _extension_functions():
    q,d=sp.symbols("q d",real=True)
    t=2*(sp.tan(q)+d)/d
    f,g=sp.exp(-1/t),sp.exp(-1/(1-t))
    base=3*(sp.pi-q)+3*sp.sin(q)*sp.cos(q)-sp.sin(q)**2
    w=3*sp.pi+f/(f+g)*(base-3*sp.pi)
    return sp.lambdify((q,d),[sp.diff(w,q,n) for n in range(5)],"numpy",cse=True)


def profile_jets(q,width=.1,flat=False):
    """W and its first four derivatives in q=atan(u); d eta/dq=-1/sqrt(W)."""
    if flat or q<=-atan(width):return np.array([3*pi,0.,0.,0.,0.])
    if q < -atan(width/2):return np.asarray(_extension_functions()(q,width),dtype=float)
    s,c=np.sin(2*q),np.cos(2*q)
    return np.array([3*(pi-q)+1.5*s-np.sin(q)**2,
                     -3+3*c-s,-6*s-2*c,-12*c+4*s,24*s+8*c])


@lru_cache(maxsize=1)
def _adiabatic_functions():
    p,m,H,H1,H2,H3=sp.symbols("p m H H1 H2 H3",real=True)
    omega=sp.sqrt(p*p+m*m)
    derivative=lambda v: -H*p*sp.diff(v,p)+H1*sp.diff(v,H)+H2*sp.diff(v,H1)+H3*sp.diff(v,H2)
    angle_dot=m*p*H/omega**2
    C1=-angle_dot/(2*omega)
    B2=sp.simplify(derivative(C1)/(2*omega))
    A2=-C1*C1/2
    C3=sp.simplify(-(derivative(B2)+angle_dot*A2)/(2*omega))
    A4=sp.simplify(-C1*C3-(A2*A2+B2*B2)/2)
    B4=sp.simplify(derivative(C3)/(2*omega))
    energy2=sp.simplify(-omega*A2)
    energy4=sp.simplify(-omega*A4)
    pressure2=sp.simplify(-p*(A2*p-B2*m)/omega)
    pressure4=sp.simplify(-p*(A4*p-B4*m)/omega)
    return sp.lambdify((p,m,H,H1,H2,H3),(energy2,pressure2,energy4,pressure4),"numpy",cse=True)


@lru_cache(maxsize=1)
def _expansion_functions():
    w=sp.symbols("w0:5",real=True)
    derivative=lambda v:-sp.sqrt(w[0])*sum(sp.diff(v,w[i])*w[i+1] for i in range(4))
    H=-w[1]/(2*sp.sqrt(w[0]))
    H1=sp.simplify(derivative(H));H2=sp.simplify(derivative(H1));H3=sp.simplify(derivative(H2))
    return sp.lambdify(w,(H,H1,H2,H3),"numpy",cse=True)


def adiabatic_coefficients(p,m,jets):
    return _adiabatic_functions()(p,m,*_expansion_functions()(*jets))


def curvature_squared_tensor(jets):
    w,w1,w2,w3,w4=jets
    energy=.5*w2*w2-w1*w3
    mixed_z=energy-2*w*w4
    return energy,mixed_z


def static_cylinder_density(normalization=1.,radius=1.):
    """Proper-time-renormalized full 4D Dirac cylinder; joint zeta normalization."""
    zeta_prime=float(mp.diff(lambda x:mp.zeta(x),-3))
    return -(2*zeta_prime+(np.log(normalization**2*radius**2)+1-float(mp.euler))/120)/(4*pi*pi*radius**4)


def _unitary_step(positive,negative,hy,hz,step):
    norm=np.sqrt(hy*hy+hz*hz)
    angle=step*norm
    cosine=np.cos(angle)
    sine_over=np.sinc(angle/pi)*step
    z=sine_over*hz;y=sine_over*hy
    return (cosine-1j*z)*positive-y*negative, y*positive+(cosine+1j*z)*negative


def smooth_integration_window(values,begin,end):
    """A numerical C-infinity taper of the convergent remainder, not a regulator."""
    if not 0<=begin<end:raise ValueError("ordered nonnegative window edges required")
    x=np.asarray(values,dtype=float);out=np.ones_like(x)
    out[x>=end]=0.;inside=(x>begin)&(x<end)
    t=(x[inside]-begin)/(end-begin)
    first=np.exp(-1/t);second=np.exp(-1/(1-t))
    out[inside]=second/(first+second)
    return out


@lru_cache(maxsize=1)
def _native_transport():
    source=Path(__file__).with_name("_angular_transport.cpp")
    flags=("-O3","-std=c++17","-shared","-fPIC","-pthread")
    digest=hashlib.sha256(source.read_bytes()+" ".join(flags).encode()).hexdigest()[:16]
    library=Path(tempfile.gettempdir())/("nsc-angular-"+digest+".so")
    if not library.exists():
        subprocess.run(["c++",*flags,str(source),"-o",str(library)],check=True,capture_output=True)
    loaded=ctypes.CDLL(str(library));function=loaded.angular_transport
    pointer=ctypes.POINTER(ctypes.c_double)
    function.argtypes=[ctypes.c_int]*3+[pointer,pointer,ctypes.c_int]+[pointer]*7+[ctypes.c_int]
    function.restype=ctypes.c_int
    return loaded,function


def angular_reference(points=(.1,.2,.5,1.,pi/2),angular_max=24,
                      momentum_points=96,momentum_ratio=60.,steps=6000,
                      width=.1,normalization=1.,flat=False,momentum_max=None,
                      backend="python",workers=8,angular_start=1,momentum_min=0.,
                      momentum_windows=(),integrand_directory=None):
    """Evolve the specified reference, sum subtracted modes and restore local terms.

    Four-dimensional angular degeneracy is 4*kappa; the +/-k integral is
    performed using parity. Convergence must be assessed before accepting
    a finite-radius value. No asymptotic-mode tail is silently extrapolated.
    """
    points=tuple(sorted(points))
    start=-atan(width)
    if not points or points[0]<start or angular_max<1 or momentum_points<4 or steps<4 or not 1<=angular_start<=angular_max:
        raise ValueError("ordered reference region and resolved positive mode counts required")
    nodes,weights=leggauss(momentum_points)
    masses=np.arange(angular_start,angular_max+1,dtype=float)[:,None]
    mode_count=len(masses)
    extent=momentum_ratio if momentum_max is None else momentum_max
    if not 0<=momentum_min<extent:raise ValueError("ordered nonnegative momentum integration band required")
    grid=momentum_min+(nodes+1)*(extent-momentum_min)/2;weights=weights*(extent-momentum_min)/2
    if momentum_windows and momentum_max is None:raise ValueError("windows require an absolute momentum grid")
    y=grid if momentum_max is None else grid[None,:]/masses
    positive=np.zeros((mode_count,momentum_points),complex)
    negative=np.ones_like(positive)
    output=[]
    step_bound=(points[-1]-start)/steps
    earlier=start
    a=(3-2*sqrt(3))/12;b=(3+2*sqrt(3))/12
    c1=.5-sqrt(3)/6;c2=.5+sqrt(3)/6

    def generator(q):
        jets=profile_jets(q,width,flat);w,w1=jets[:2]
        hz=-masses*np.sqrt(w+y*y)/w
        hy=-y*w1/(4*np.sqrt(w)*(w+y*y))
        return hy,hz

    for point in points:
        count=max(1,ceil((point-earlier)/step_bound));step=(point-earlier)/count
        if backend=="native":
            first=np.array([profile_jets(earlier+(i+c1)*step,width,flat)[:2] for i in range(count)])
            second=np.array([profile_jets(earlier+(i+c2)*step,width,flat)[:2] for i in range(count)])
            arrays=[np.ascontiguousarray(v,dtype=float) for v in
                    (masses[:,0],grid,np.full(count,step),first[:,0],first[:,1],second[:,0],second[:,1])]
            ptr=lambda v:v.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            _,transport=_native_transport()
            status=transport(mode_count,momentum_points,count,ptr(arrays[0]),ptr(arrays[1]),int(momentum_max is None),
                             *(ptr(v) for v in arrays[2:]),ptr(positive),ptr(negative),workers)
            if status:raise RuntimeError("native angular transport rejected its dimensions")
        elif backend=="python":
            for i in range(count):
                q=earlier+i*step
                y1,z1=generator(q+c1*step);y2,z2=generator(q+c2*step)
                positive,negative=_unitary_step(positive,negative,b*y1+a*y2,b*z1+a*z2,step)
                positive,negative=_unitary_step(positive,negative,a*y1+b*y2,a*z1+b*z2,step)
        else:raise ValueError("backend must be python or native")
        jets=profile_jets(point,width,flat);s=np.sqrt(jets[0])
        p=y/s;omega=np.sqrt(1+p*p)
        occupancy=abs(positive)**2
        tangent=-2*np.real(positive.conj()*negative)
        delta_energy=2*masses*omega*occupancy
        delta_pressure=2*masses*p*p/omega*occupancy+masses*p/omega*tangent
        E2,P2,E4,P4=adiabatic_coefficients(p,1.,jets)
        energy_remainder=delta_energy-E2/masses-E4/masses**3
        pressure_remainder=delta_pressure-P2/masses-P4/masses**3
        factor=masses[:,0]**(2 if momentum_max is None else 1)/(pi*pi*s)
        energy_terms=factor*(energy_remainder@weights)
        pressure_terms=factor*(pressure_remainder@weights)
        windowed=[]
        for low,high in momentum_windows:
            weighted=weights*smooth_integration_window(grid,low,high)
            windowed.append({"begin":low,"end":high,
                             "angular_energy_remainders":(factor*(energy_remainder@weighted)).tolist(),
                             "angular_pressure_remainders":(factor*(pressure_remainder@weighted)).tolist()})
        if integrand_directory is not None:
            folder=Path(integrand_directory);folder.mkdir(parents=True,exist_ok=True)
            target=folder/("q_"+format(point,".12g")+".npz")
            if target.exists():raise FileExistsError("refusing to overwrite scientific integrands")
            np.savez_compressed(target,energy=energy_remainder,pressure=pressure_remainder,
                                factor=factor,grid=grid,weights=weights,masses=masses[:,0],q=point,W_jets=jets)
        henergy,hz=curvature_squared_tensor(jets)
        harmonic=np.log(normalization)+float(mp.euler)/2
        rho_static=static_cylinder_density(normalization)
        rho_local=harmonic*henergy/(480*pi*pi)
        pz_local=-harmonic*hz/(480*pi*pi)
        rho=rho_static+rho_local+float(sum(energy_terms))
        pz=-rho_static+pz_local+float(sum(pressure_terms))
        w,w1,w2,w3,w4=jets;R2=-w2;boxR=-w*w4-w1*w3
        anomaly=(-(R2-2)**2/60-11*R2/90-boxR/30)/(16*pi*pi)
        output.append({"q":point,"rho_coordinate":-1/np.tan(point),"sphere_radius":1/np.sin(point),
                       "bar_rho":rho,"bar_p_parallel":pz,"bar_p_sphere":(rho-pz-anomaly)/2,
                       "bar_trace":anomaly,"static_density":rho_static,
                       "restored_curvature_rho":rho_local,"restored_curvature_p_parallel":pz_local,
                       "angular_energy_remainders":energy_terms.tolist(),
                       "angular_pressure_remainders":pressure_terms.tolist(),
                       "windowed_momentum_remainders":windowed,
                       "unitarity_defect":float(np.max(abs(abs(positive)**2+abs(negative)**2-1))),
                       "maximum_positive_frequency_occupation":float(np.max(occupancy)),
                       "W_jets":jets.tolist()})
        earlier=point
    return {"parameters":{"angular_max":angular_max,"momentum_points":momentum_points,
                           "momentum_ratio":momentum_ratio,"steps":steps,"width":width,
                           "momentum_max":momentum_max,
                           "momentum_min":momentum_min,
                           "momentum_windows":[list(v) for v in momentum_windows],
                           "angular_start":angular_start,"backend":backend,"workers":workers,
                           "normalization":normalization,"flat_control":flat},"rows":output}


@lru_cache(maxsize=1)
def _metric_jet_functions():
    w=sp.symbols("w0:5",positive=True)
    derivative=lambda v:sum(sp.diff(v,w[i])*w[i+1] for i in range(4))
    n=[1/sp.sqrt(w[0])];a=[sp.sqrt(w[0])]
    for _ in range(4):n.append(derivative(n[-1]));a.append(derivative(a[-1]))
    return sp.lambdify(w,(*n,*a),"numpy",cse=True)


@lru_cache(maxsize=1)
def _conformal_source_functions():
    """Bulk variations of the integrated anomaly, with sigma held fixed.

    Sphere radius is fixed ONLY for the two base-metric variations. The
    sphere pressure is recovered from the exact 4D trace after variation.
    Variations have compact support, so the displayed integration by parts
    does not discard a physical endpoint stress.
    """
    n=sp.symbols("n0:5",positive=True);a=sp.symbols("a0:5",positive=True)
    s=sp.symbols("s0:5",real=True)
    derivative=lambda v:sum(sp.diff(v,row[i])*row[i+1] for row in (n,a,s) for i in range(4))
    Ha=a[1]/(n[0]*a[0])
    acceleration=(a[2]-a[1]*n[1]/n[0])/(n[0]**2*a[0])
    sd=s[1]/n[0];sdd=(s[2]-s[1]*n[1]/n[0])/n[0]**2
    U=sdd+Ha*sd;v=sd*sd
    euler_average=8*(acceleration+U/2+(acceleration*v+2*Ha*sd*sdd)/3+(U*v+2*v*sdd)/4)
    R0=-2*(acceleration+1)
    alpha,beta,gamma=-sp.Rational(1,20),sp.Rational(11,360),-sp.Rational(1,30)
    density=n[0]*a[0]/(16*sp.pi**2)*(alpha*s[0]*sp.Rational(4,3)*(acceleration+1)**2
        +beta*s[0]*euler_average+gamma*(R0*(U+v)-3*(U+v)**2))
    def euler(row):
        return sp.diff(density,row[0])-derivative(sp.diff(density,row[1]))+derivative(derivative(sp.diff(density,row[2])))
    rho=euler(n)/a[0]
    pressure=-euler(a)/n[0]
    return sp.lambdify((*n,*a,*s),(rho,pressure),"numpy",cse=True)


@lru_cache(maxsize=1)
def _physical_trace_function():
    q=sp.Symbol("q",real=True);w=sp.symbols("w0:5",real=True)
    derivative=lambda v:sp.diff(v,q)+sum(sp.diff(v,w[i])*w[i+1] for i in range(4))
    s1=-sp.cot(q);s2=1/sp.sin(q)**2
    A=w[2]/2+w[1]*s1+w[0]*s2
    B=w[1]*s1/2+w[0]*s1*s1
    C=w[0]*s2+w[1]*s1/2
    D=1+w[0]*s1*s1
    scalar=-2*sp.sin(q)**2*(A+2*B+2*C+D)
    weyl=sp.sin(q)**4*(w[2]+2)**2/3
    euler=8*sp.sin(q)**4*(A*D+2*B*C)
    box=sp.sin(q)**2*(w[0]*derivative(derivative(scalar))+(w[1]-2*w[0]*sp.cot(q))*derivative(scalar))
    trace=(-weyl/20+11*euler/360-box/30)/(16*sp.pi**2)
    return sp.lambdify((q,*w),(scalar,trace),"numpy",cse=True)


def physical_source(row):
    """Conformal transport of the FULL computed reference tensor at finite q."""
    q=row["q"]
    w=np.asarray(row["W_jets"])
    if not 0<q<pi or w[0]<=0:raise ValueError("finite timelike interior point with W>0 required")
    s,c=np.sin(q),np.cos(q);cot=c/s
    sigma=(-np.log(s),-cot,1/s**2,-2*cot/s**2,4*cot*cot/s**2+2/s**4)
    jets=_metric_jet_functions()(*w)
    anomaly_rho,anomaly_pz=_conformal_source_functions()(*jets,*sigma)
    rho=s**4*(row["bar_rho"]+anomaly_rho)
    pz=s**4*(row["bar_p_parallel"]+anomaly_pz)
    scalar,trace=_physical_trace_function()(q,*w)
    return {"q":q,"rho_coordinate":row["rho_coordinate"],"sphere_radius":row["sphere_radius"],
            "rho":float(rho),"p_parallel":float(pz),"p_sphere":float((rho-pz-trace)/2),
            "radial_null":float(rho+pz),"trace":float(trace),"R_L":float(scalar),
            "transported_bar_rho":float(s**4*row["bar_rho"]),
            "local_anomaly_rho":float(s**4*anomaly_rho),
            "H_parallel":float(np.sqrt(w[0])*c-w[1]*s/(2*np.sqrt(w[0]))),
            "H_sphere":float(np.sqrt(w[0])*c)}
