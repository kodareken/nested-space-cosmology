import numpy as np

# KS background with Weyssenhoff spin fluid + radiation
# metric: ds^2 = -dT^2 + A(T)^2 dR^2 + B(T)^2 dOmega^2
# dimensionless: tau = T/B_b, b = B/B_b, Vbar = V/V_b (V = A B^2), h_a = B_b H_A, h_b = B_b H_B
# bounce: h_a = h_b = 0, b = 1, Vbar = 1, 1/B_b^2 = kappa eps_bb
# Sbar = spin energy / eps_bb at bounce,  Rbar = 1 + Sbar
# constraint: 2 h_a h_b + h_b^2 + 1/b^2 = Rbar V^-4/3 - Sbar V^-2

def deriv(t, y, S):
    b, V, ha, hb = y
    R = 1.0 + S
    P = R * V**(-4.0/3.0)/3.0 - S * V**(-2.0)
    db  = hb * b
    dV  = (ha + 2.0*hb) * V
    dhb = -(3.0*hb*hb + 1.0/b**2 + P)/2.0
    dha = -(ha*ha + hb*hb + ha*hb + P) - dhb
    return np.array([db, dV, dha, dhb])

def constraint(y, S):
    b, V, ha, hb = y
    R = 1.0 + S
    eps = R * V**(-4.0/3.0) - S * V**(-2.0)
    return 2.0*ha*hb + hb*hb + 1.0/b**2 - eps

def rk4(f, y0, t0, t1, n, S, stop=None):
    h = (t1 - t0)/n
    t = t0; y = np.array(y0, dtype=float)
    ts = [t0]; ys = [np.array(y0, dtype=float)]
    for _ in range(n):
        if stop is not None and stop(t, y):
            break
        k1 = f(t, y, S)
        k2 = f(t + h/2, y + h*k1/2, S)
        k3 = f(t + h/2, y + h*k2/2, S)
        k4 = f(t + h, y + h*k3, S)
        y = y + h*(k1 + 2*k2 + 2*k3 + k4)/6
        t += h
        ts.append(t); ys.append(y.copy())
    return np.array(ts), np.array(ys)

for S in [2.5, 3.0, 5.0, 10.0, 100.0]:
    y0 = [1.0, 1.0, 0.0, 0.0]
    back = rk4(deriv, y0, 0.0, -40.0, 400000, S,
               stop=lambda t, y: y[0] < 1e-3 or y[1] < 1e-3)
    fwd  = rk4(deriv, y0, 0.0,  40.0, 400000, S,
               stop=lambda t, y: y[0] < 1e-3 or y[1] < 1e-3)

    tb, yb = back[:,0], back[:,1:]
    tf, yf = fwd[:,0],  fwd[:,1:]

    Cback = np.array([constraint(y, S) for y in yb])
    Cfwd  = np.array([constraint(y, S) for y in yf])
    Cmax = max(np.max(np.abs(Cback)), np.max(np.abs(Cfwd)))

    shearb = (yb[:,2] - yb[:,3])/np.sqrt(3.0)
    thb = (yb[:,2] + 2.0*yb[:,3])/3.0
    # find max Vbar backward (past turnaround)
    imax = np.argmax(yb[:,1])
    print(f"S={S:8.1f}  Cmax={Cmax:.2e}")
    print(f"  BACKWARD: pancake at tau={tb[-1]:.4f}, V_past_max={yb[imax,1]:.4f} at tau={tb[imax]:.4f}")
    print(f"    past h_a={yb[-1,2]:.4f} h_b={yb[-1,3]:.4f}  shear/|th| last={abs(shearb[-1]/thb[-1]):.4f}")
    print(f"  FORWARD:  ends tau={tf[-1]:.4f}, b={yf[-1,0]:.4e} V={yf[-1,1]:.4e} h_a={yf[-1,2]:.4e} h_b={yf[-1,3]:.4e}")
    shf = (yf[:,2]-yf[:,3])/np.sqrt(3.0); thf = (yf[:,2]+2.0*yf[:,3])/3.0
    print(f"    fwd shear/|th| last={abs(shf[-1]/thf[-1]):.4f}   (1=Kasner-like)")
    print()
