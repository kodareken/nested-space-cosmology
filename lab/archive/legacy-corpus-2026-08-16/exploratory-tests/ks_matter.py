import numpy as np

# KS + MATTER (w=0) + Weyssenhoff spin fluid
# eps_eff = R V^-1 - S V^-2,  P_eff = 0 - S V^-2  (dust + stiff spin, spin p_s = eps_s)
# constraint: 2 ha hb + hb^2 + 1/b^2 = eps_eff
# bounce: ha=hb=0, b=1, V=1,  Rbar = 1+Sbar,  Sbar > 1 required

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

def constraint(y, S):
    b, V, ha, hb = y
    R = 1.0 + S
    eps = R * V**(-1.0) - S * V**(-2.0)
    return 2.0*ha*hb + hb*hb + 1.0/b**2 - eps

def rk4(f, y0, t0, t1, n, S, stop=None):
    ts = [t0]; ys = [np.array(y0, dtype=float)]
    h = (t1-t0)/n; t = t0; y = np.array(y0, dtype=float)
    for _ in range(n):
        if stop is not None and stop(t, y): break
        k1 = f(t, y, S); k2 = f(t+h/2, y+h*k1/2, S)
        k3 = f(t+h/2, y+h*k2/2, S); k4 = f(t+h, y+h*k3, S)
        y = y + h*(k1+2*k2+2*k3+k4)/6; t += h
        ts.append(t); ys.append(y.copy())
    return np.array(ts), np.array(ys)

for S in [1.2, 1.5, 2.0, 3.0, 5.0]:
    tb, yb = rk4(deriv, [1.,1.,0.,0.], 0.0, -30.0, 300000, S, stop=lambda t,y: y[0]<1e-3 or y[1]<1e-3)
    tf, yf = rk4(deriv, [1.,1.,0.,0.], 0.0,  30.0, 300000, S, stop=lambda t,y: y[0]<1e-3 or y[1]<1e-3)
    Cb = np.array([constraint(y,S) for y in yb]); Cf = np.array([constraint(y,S) for y in yf])
    Cmax = max(np.max(np.abs(Cb[yb[:,1]<2.0])), np.max(np.abs(Cf[yf[:,1]<2.0])))
    imax = np.argmax(yb[:,1])
    sh = (yf[:,2]-yf[:,3])/np.sqrt(3); th = (yf[:,2]+2*yf[:,3])/3
    print(f"S={S:4.1f}: Cmax={Cmax:.2e}  pancake tau={tb[-1]:.4f}  V_past_max={yb[imax,1]:.4f}  "
          f"shear/|th| late={abs(sh[-1]/th[-1]):.3f}  eps_S/eps_m at bounce={S/(1+S):.3f}")
