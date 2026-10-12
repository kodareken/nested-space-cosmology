import numpy as np, time

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
t0 = time.time()
y0 = [1.,1.,0.,0.]
t1a, y1a = rk4seg(deriv, y0, 0.0, -1.45, 4.0e-5, S, stop=lambda t,y: y[0]<0.05 or y[1]<0.2)
t2a, y2a = rk4seg(deriv, y1a[-1], t1a[-1], -2.95, 4.0e-5, S, stop=lambda t,y: y[0]<0.05 or y[1]<0.2)
tb = np.concatenate([t2a[::-1][:-1], t1a[::-1][:-1]])
yb = np.concatenate([y2a[::-1][:-1], y1a[::-1][:-1]])
t3a, y3a = rk4seg(deriv, y0, 0.0, 1.45, 4.0e-5, S, stop=lambda t,y: y[0]<0.05 or y[1]<0.2)
t4a, y4a = rk4seg(deriv, y3a[-1], t3a[-1], 2.95, 4.0e-5, S, stop=lambda t,y: y[0]<0.05 or y[1]<0.2)
tf = np.concatenate([t3a, t4a[1:]])
yf = np.concatenate([y3a, y4a[1:]])
# fine core
t1b, y1b = rk4seg(deriv, y1a[-1], t1a[-1], 0.0, 4.0e-6, S)
t3b, y3b = rk4seg(deriv, y0, 0.0, 1.45, 4.0e-6, S)
tau = np.concatenate([tb, t1b[1:], t3b[1:], tf[1:]])
Y = np.concatenate([yb, y1b[1:], y3b[1:], yf[1:]])
print(f"trajectory: {len(tau)} pts over [{tau[0]:.2f},{tau[-1]:.2f}], {time.time()-t0:.1f}s")

b_all=Y[:,0]; Vall=Y[:,1]; ha=Y[:,2]; hb=Y[:,3]
a_all = Vall/b_all**2
R = 1.0+S
epP = R*Vall**-1 - 2.0*S*Vall**-2
Nt = len(tau)

nec_past=None; nec_fut=None
for j in range(1,Nt):
    if epP[j-1] > 0 and epP[j] <= 0 and nec_past is None and tau[j] < 0: nec_past = tau[j]
    if epP[j-1] <= 0 and epP[j] > 0 and nec_fut is None and tau[j] > 0: nec_fut = tau[j-1]
print(f"NEC zeros: {nec_past:.4f}, {nec_fut:.4f}")

z2A = a_all**2 * epP/ha**2
z2B = b_all**2 * epP/hb**2

def build_zppz(z2):
    zpp = np.zeros(Nt)
    h = tau[1]-tau[0]
    for j in range(Nt):
        t = tau[j]
        if abs(t) < 0.10:
            x = t
            zpp[j] = (2.0/x**2 + 2.0) if abs(x) > 1e-10 else 1e8
        elif abs(t - nec_past) < 0.10 or abs(t - nec_fut) < 0.10:
            tn = nec_past if t < 0 else nec_fut
            x = t - tn
            if abs(x) < 1e-10:
                zpp[j] = -1e8
            else:
                jz = np.argmin(np.abs(tau - tn))
                z2p = (z2[jz+1] - z2[jz-1])/(tau[jz+1]-tau[jz-1])
                z2pp = (z2[jz+1] - 2*z2[jz] + z2[jz-1])/((tau[jz+1]-tau[jz-1])/2)**2
                c2c1 = z2pp/z2p
                zpp[j] = -1.0/(4.0*x*x) - c2c1/(2.0*x)
        else:
            if j < 3:
                z2p = (-3*z2[0]+4*z2[1]-z2[2])/(2*h)
                z2pp = (2*z2[0]-5*z2[1]+4*z2[2]-z2[3])/h**2
            elif j > Nt-4:
                z2p = (3*z2[j]-4*z2[j-1]+z2[j-2])/(2*h)
                z2pp = (2*z2[j]-5*z2[j-1]+4*z2[j-2]-z2[j-3])/h**2
            else:
                z2p = (z2[j+1]-z2[j-1])/(2*h)
                z2pp = (z2[j+1]-2*z2[j]+z2[j-1])/h**2
            zpp[j] = z2pp/(2*z2[j]) - z2p*z2p/(4*z2[j]**2)
    return zpp

zppA = build_zppz(z2A)
zppB = build_zppz(z2B)
print(f"zpp built, {time.time()-t0:.1f}s; max|zpp| A={np.abs(zppA).max():.1e} B={np.abs(zppB).max():.1e}")

# cumulative conformal-ish time for BD phase: eta_i = int dt/a_i
etaA = np.cumsum(np.gradient(tau)/a_all)
etaB = np.cumsum(np.gradient(tau)/b_all)

def solve(k, dirn):
    zpp = zppA if dirn=='A' else zppB
    ai = a_all if dirn=='A' else b_all
    eta = etaA if dirn=='A' else etaB
    # start at first index where k^2/a^2 >= 50|zpp| (deep sub-horizon), min 2
    ka = k*k/ai**2
    jstart = np.argmax(ka >= 50.0*np.abs(zpp))
    if jstart < 2: jstart = 2
    ph = -k*eta[jstart]
    vr = np.cos(ph)/np.sqrt(2.0*k); vi = np.sin(ph)/np.sqrt(2.0*k)
    wr = 0.0; wi = -k/ai[jstart]*vi + 0.0
    wr2 = -k/ai[jstart]*vi; wi2 = k/ai[jstart]*vr
    wr, wi = wr2, wi2
    for j in range(jstart, Nt-1):
        h = tau[j+1]-tau[j]
        f = zpp[j] - ka[j]
        k1r = wr;          k1i = wi;          k1wr = f*vr;          k1wi = f*vi
        k2r = wr+h*k1wr/2; k2i = wi+h*k1wi/2; k2wr = f*(vr+h*k1r/2); k2wi = f*(vi+h*k1i/2)
        k3r = wr+h*k2wr/2; k3i = wi+h*k2wi/2; k3wr = f*(vr+h*k2r/2); k3wi = f*(vi+h*k2i/2)
        k4r = wr+h*k3wr;   k4i = wi+h*k3wi;   k4wr = f*(vr+h*k3r);   k4wi = f*(vi+h*k3i)
        vr += h*(k1r+2*k2r+2*k3r+k4r)/6
        vi += h*(k1i+2*k2i+2*k3i+k4i)/6
        wr += h*(k1wr+2*k2wr+2*k3wr+k4wr)/6
        wi += h*(k1wi+2*k2wi+2*k3wi+k4wi)/6
    jend = Nt//2 + np.argmax(Vall[Nt//2:])
    z2end = (a_all[jend]**2 if dirn=='A' else b_all[jend]**2) * epP[jend] / (ha[jend]**2 if dirn=='A' else hb[jend]**2)
    zend = np.sqrt(abs(z2end))
    return (vr*vr + vi*vi)/(zend*zend)

kk = [1.5, 2, 3, 4, 6, 8, 12, 16, 24]
for dirn in ['A','B']:
    P = []
    t2 = time.time()
    for k in kk:
        amp = solve(k, dirn)
        P.append(k**3*amp)
        print(f"  dir {dirn}, k={k:4.1f}: k^3|zeta|^2 = {k**3*amp:.6e}   ({time.time()-t2:.0f}s)")
    print(f"  n_s({dirn}) = {1 + np.polyfit(np.log(kk), np.log(P), 1)[0]:.4f}")
