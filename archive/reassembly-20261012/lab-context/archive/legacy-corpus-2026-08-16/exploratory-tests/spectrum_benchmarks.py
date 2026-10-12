import numpy as np

# ============================================================
# THE STIFF-BOUNCE SPECTRUM + THE SETTING OF c
# BM2: dust+spin FLRW bounce (eps_eff = V^-1 - V^-2), direction-resolved
# curvature-perturbation spectrum with the zeta-matching across the
# singular points. RESULT (cited): n_s = 2.94 (blue) - the stiff bounce's
# horizon opens and closes inside the turn (aH peaks at 0.397 at V=4,
# dies at the bounce), so its own patch genuinely blues; the observed
# 0.9649 must come from the parent's outer dust epoch (the two-region
# junction). The earlier BM1 (spliced-bounce benchmark) was a degenerate
# design and is removed; the pipeline certification is dS_benchmark.py.
# ============================================================
def deriv(t, y, S):
    b, V, ha, hb = y
    R = 1.0 + S
    eps = R * V**(-1.0) - S * V**(-2.0)
    P = -S * V**(-2.0)
    db  = hb * b
    dV  = (ha + 2.0*hb) * V
    dhb = -(3.0*hb*hb + 1.0/b**2 + P)/2.0
    dha = -(ha*ha + hb*hb + ha*hb + P) - dhb
    return np.array([db, dV, dha, dhb])

def rk4seg(f, y0, t0, t1, hstep, S, stop=None):
    n = int(abs(t1-t0)/hstep) + 1
    h = (t1-t0)/n
    ts=[t0]; ys=[np.array(y0,dtype=float)]; t=t0; y=np.array(y0,dtype=float)
    for _ in range(n):
        if stop is not None and stop(t,y): break
        k1=f(t,y,S); k2=f(t+h/2,y+h*k1/2,S); k3=f(t+h/2,y+h*k2/2,S); k4=f(t+h,y+h*k3,S)
        y=y+h*(k1+2*k2+2*k3+k4)/6; t+=h
        ts.append(t); ys.append(y.copy())
    return np.array(ts), np.array(ys)

S = 2.0
t1a, y1a = rk4seg(deriv, [1.,1.,0.,0.], 0.0, -1.5, 4.0e-5, S, stop=lambda t,y: y[0]<0.05 or y[1]<0.2)
t2a, y2a = rk4seg(deriv, y1a[-1], t1a[-1], -2.65, 4.0e-5, S, stop=lambda t,y: y[0]<0.05 or y[1]<0.2)
tb = np.concatenate([t2a[::-1][:-1], t1a[::-1][:-1]]); yb = np.concatenate([y2a[::-1][:-1], y1a[::-1][:-1]])
t3a, y3a = rk4seg(deriv, [1.,1.,0.,0.], 0.0, 1.5, 4.0e-5, S, stop=lambda t,y: y[0]<0.05 or y[1]<0.2)
t4a, y4a = rk4seg(deriv, y3a[-1], t3a[-1], 2.5, 4.0e-5, S, stop=lambda t,y: y[0]<0.05 or y[1]<0.2)
tf = np.concatenate([t3a,t4a[1:]]); yf=np.concatenate([y3a,y4a[1:]])
tau = np.concatenate([tb, t1b[1:] if False else tf]) if False else np.concatenate([tb, t3a, t4a[1:]])
Y = np.concatenate([yb, y3a, y4a[1:]])
b_all=Y[:,0]; Vall=Y[:,1]; ha=Y[:,2]; hb=Y[:,3]
a_all = Vall/b_all**2
R = 1.0+S
epP = R*Vall**-1 - 2.0*S*Vall**-2
Nt = len(tau)
nec_past=None; nec_fut=None
for j in range(1,Nt):
    if epP[j-1] > 0 and epP[j] <= 0 and nec_past is None and tau[j] < 0: nec_past=tau[j]
    if epP[j-1] <= 0 and epP[j] > 0 and nec_fut is None and tau[j] > 0: nec_fut=tau[j-1]
z2A = a_all**2*epP/ha**2
z2B = b_all**2*epP/hb**2

def build_zppz(z2):
    zpp=np.zeros(Nt); h=tau[1]-tau[0]
    for j in range(Nt):
        t=tau[j]
        if abs(t)<0.10:
            x=t
            zpp[j]=(2.0/x**2+2.0) if abs(x)>1e-10 else 1e8
        elif abs(t-nec_past)<0.10 or abs(t-nec_fut)<0.10:
            tn=nec_past if t<0 else nec_fut
            x=t-tn
            if abs(x)<1e-10: zpp[j]=-1e8
            else:
                jz=np.argmin(np.abs(tau-tn))
                z2p=(z2[jz+1]-z2[jz-1])/(tau[jz+1]-tau[jz-1])
                z2pp=(z2[jz+1]-2*z2[jz]+z2[jz-1])/((tau[jz+1]-tau[jz-1])/2)**2
                c2c1=z2pp/z2p
                zpp[j]=-1.0/(4.0*x*x)-c2c1/(2.0*x)
        else:
            if j<3:
                z2p=(-3*z2[0]+4*z2[1]-z2[2])/(2*h); z2pp=(2*z2[0]-5*z2[1]+4*z2[2]-z2[3])/h**2
            elif j>Nt-4:
                z2p=(3*z2[j]-4*z2[j-1]+z2[j-2])/(2*h); z2pp=(2*z2[j]-5*z2[j-1]+4*z2[j-2]-z2[j-3])/h**2
            else:
                z2p=(z2[j+1]-z2[j-1])/(2*h); z2pp=(z2[j+1]-2*z2[j]+z2[j-1])/h**2
            zpp[j]=z2pp/(2*z2[j])-z2p*z2p/(4*z2[j]**2)
    return zpp

zppA=build_zppz(z2A); zppB=build_zppz(z2B)

def solve(k, dirn):
    zpp = zppA if dirn=='A' else zppB
    ai = a_all if dirn=='A' else b_all
    z2 = z2A if dirn=='A' else z2B
    ka = k*k/ai**2
    jstart = np.argmax((ka >= 50.0*np.abs(zpp)) & (tau <= tau[np.argmin(np.abs(tau - (nec_past - 0.05)))]))
    if jstart < 2: jstart = 2
    vr = 1.0/np.sqrt(2.0*k); vi=0.0
    wr = 0.0; wi = -k/ai[jstart]*vr
    jmatch = np.argmin(np.abs(tau - (nec_past - 0.05)))
    for j in range(jstart, jmatch):
        h = tau[j+1]-tau[j]
        f = zpp[j] - ka[j]
        k1r=wr; k1i=wi; k1wr=f*vr; k1wi=f*vi
        k2r=wr+h*k1wr/2; k2i=wi+h*k1wi/2; k2wr=f*(vr+h*k1r/2); k2wi=f*(vi+h*k1i/2)
        k3r=wr+h*k2wr/2; k3i=wi+h*k2wi/2; k3wr=f*(vr+h*k2r/2); k3wi=f*(vi+h*k2i/2)
        k4r=wr+h*k3wr; k4i=wi+h*k3wi; k4wr=f*(vr+h*k3r); k4wi=f*(vi+h*k3i)
        vr+=h*(k1r+2*k2r+2*k3r+k4r)/6; vi+=h*(k1i+2*k2i+2*k3i+k4i)/6
        wr+=h*(k1wr+2*k2wr+2*k3wr+k4wr)/6; wi+=h*(k1wi+2*k2wi+2*k3wi+k4wi)/6
    return (vr*vr+vi*vi)/abs(z2[jmatch])

kk=[4,5,6,8,10,14,20,28,40]
P=[]
for k in kk:
    amp=solve(k,'A')
    P.append(k**3*amp)
ns = 1 + np.polyfit(np.log(kk), np.log(P), 1)[0]
print(f"stiff-bounce direction-A spectrum: n_s = {ns:.4f}   [cited: 2.94]")

print("\nHorizon profile: aH(V) = sqrt(1/3) sqrt(V^-1/3 - V^-4/3):")
for V in [100, 25, 9, 4, 2.5, 1.5, 1.1]:
    f = V**(-1/3) - V**(-4/3)
    print(f"  V={V:6.1f}: aH = {np.sqrt(1/3)*np.sqrt(f):.4f}")

print("""
=== THE SETTING OF c (cited: the-consistency-audit.md §4) ===
T_bb = 1.152e32 K (the bounce, P4)
T_CMB = 2.7255 +/- 0.0006 K (COBE-FIRAS)
gradient braking factor: T_bb/T_CMB = %.3e = e^%.2f  (72.8 e-folds of braking)
horizon condition: sqrt(2GM/r) = c  -  where the local gradient equals c
finite c <=> horizons form <=> pockets seal <=> nested gradients exist
""" % (1.152e32/2.7255, np.log(1.152e32/2.7255)))
