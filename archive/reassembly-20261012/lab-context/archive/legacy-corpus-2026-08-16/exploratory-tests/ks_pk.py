import numpy as np

# NOTE: RUNG 0 below (the FLRW MS rung) is superseded by dS_benchmark.py,
# which certifies the pipeline against the exact de Sitter answer.
# The CITED result of this file is RUNG 1: the KS direction-resolved
# test-scalar spectra (n_s = 2.77 R-family, 2.81 Omega-family, converged).
# ============================================================
# RUNG 0 (superseded): FLRW+spin symmetric bounce:
# evolve each mode through the bounce, extract FROZEN post-bounce zeta_k,
# fit n_s.  Analytic expectation: n_s = 3 (blue).
# Background (from earlier work): eta = sqrt(3)*acosh(y), y=a/a_b,
# z^2 = 2(2-3y^-2)/(1-y^-2), bounce at y=1, symmetric in eta.
# ============================================================
y0, y1 = 1e6, 1.000
Ngrid = 400000
y_grid = np.geomspace(y0, y1, Ngrid)                    # contracting side
eta_cont = np.sqrt(3.0)*(np.arccosh(y0) - np.arccosh(y_grid))  # 0..eta_b
z2 = 2.0*(2.0 - 3.0*y_grid**-2)/(1.0 - y_grid**-2)
z = np.sqrt(z2)

def d2(f, x):
    d2f = np.empty_like(f); h = np.diff(x); hh = h[:-1] + h[1:]
    d2f[1:-1] = 2*(f[2:]*h[:-1] - f[1:-1]*hh + f[:-2]*h[1:])/(h[:-1]*h[1:]*hh)
    d2f[0] = d2f[1]; d2f[-1] = d2f[-2]
    return d2f

zzpp_cont = d2(z, eta_cont)/z
eta_b = eta_cont[-1]
# symmetric full grid: -eta_b .. +eta_b
Nu = 600000
eta = np.linspace(-eta_b, eta_b, Nu)
zz = np.empty(Nu)
zz[:Nu//2] = np.interp(-eta[Nu//2-1::-1], eta_cont, zzpp_cont)   # eta<0 side
zz[Nu//2:] = np.interp(eta[Nu//2:], eta_cont, zzpp_cont)         # eta>0 side
# z on full grid for extracting zeta = v/z
zz2 = np.empty(Nu)
zz2[:Nu//2] = np.interp(-eta[Nu//2-1::-1], eta_cont, z2)
zz2[Nu//2:] = np.interp(eta[Nu//2:], eta_cont, z2)
z_full = np.sqrt(zz2)
h = eta[1] - eta[0]
print(f"FLRW: eta_b={eta_b:.4f}, z''/z deep past={zz[0]:.3e}")

def solve_flrw(k):
    vr = 1.0/np.sqrt(2*k); vi = 0.0
    for j in range(Nu-1):
        f = zz[j] - k*k
        k1r = vi;            k1i = f*vr
        k2r = vi + h*k1i/2;  k2i = f*(vr + h*k1r/2)
        k3r = vi + h*k2i/2;  k3i = f*(vr + h*k2r/2)
        k4r = vi + h*k3i;    k4i = f*(vr + h*k3r)
        vr += h*(k1r + 2*k2r + 2*k3r + k4r)/6
        vi += h*(k1i + 2*k2i + 2*k3i + k4i)/6
    return vr, vi

ks = [0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0]
PR = []
print("\nFLRW+spin bounce, frozen post-bounce P(k):")
for k in ks:
    vr, vi = solve_flrw(k)
    Pk = k**3*(vr*vr + vi*vi)/(2*np.pi**2*z_full[-1]**2)
    PR.append(Pk)
    print(f"  k={k:5.1f}: P_R = {Pk:.6e}")
ns = 1 + np.polyfit(np.log(ks), np.log(PR), 1)[0]
print(f"\nFitted n_s = {ns:.4f}  (analytic expectation for radiation+stiff: 3.0)")

# ============================================================
# RUNG 1: KS direction-resolved TEST-SCALAR spectrum
# metric: ds^2 = -dtau^2 + a(tau)^2 dR^2 + b(tau)^2 dOmega^2
# wave equation: f'' + (h_a+2h_b) f' + [k^2/a^2 + l(l+1)/b^2] f = 0
# vacuum: WKB adiabatic at the deep-past side, evolve through the bounce,
# extract frozen late-time amplitude. Families:
#   R-family: l=0, k>0   (radial modes)
#   O-family: k=0, l>0   (sphere modes)
# ============================================================
exec(open('/tmp/ks_background.py').read().split('for S in')[0])

S = 2.5
y0 = [1.0, 1.0, 0.0, 0.0]
# full symmetric trajectory: from pancake (tau<0) through bounce to pancake (tau>0)
tb, yb = rk4(deriv, y0, 0.0, -2.63, 263000, S, stop=lambda t, y: y[0] < 1e-3 or y[1] < 1e-3)
tf, yf = rk4(deriv, y0, 0.0,  2.63, 263000, S, stop=lambda t, y: y[0] < 1e-3 or y[1] < 1e-3)
tau = np.concatenate([tb[::-1][:-1], tf])
Y = np.concatenate([yb[::-1][:-1], yf])
b_all = Y[:,0]; Vall = Y[:,1]; ha = Y[:,2]; hb = Y[:,3]
a_all = Vall/b_all**2
Nt = len(tau)

def solve_ks(k, l):
    # WKB vacuum ICs at first grid point (deep past side)
    w2 = k*k/a_all[0]**2 + l*(l+1.0)/b_all[0]**2
    w = np.sqrt(w2)
    fr = 1.0/np.sqrt(2*w); fi = 0.0
    d = np.zeros_like(tau)
    # state: f, f'; integrate with the ODE f'' = -(ha+2hb) f' - w2(t) f
    dpr = -w*fi; dpi = w*fr   # f' = -i w f at start
    for j in range(Nt-1):
        dt = tau[j+1] - tau[j]
        c = ha[j] + 2.0*hb[j]
        w2j = k*k/a_all[j]**2 + l*(l+1.0)/b_all[j]**2
        def F(r, i, pr, pi):
            return pr, pi, -c*pr - w2j*r, -c*pi - w2j*i
        k1 = F(fr, fi, dpr, dpi)
        k2 = F(fr+dt*k1[0]/2, fi+dt*k1[1]/2, dpr+dt*k1[2]/2, dpi+dt*k1[3]/2)
        k3 = F(fr+dt*k2[0]/2, fi+dt*k2[1]/2, dpr+dt*k2[2]/2, dpi+dt*k2[3]/2)
        k4 = F(fr+dt*k3[0], fi+dt*k3[1], dpr+dt*k3[2], dpi+dt*k3[3])
        fr += dt*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6
        fi += dt*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6
        dpr += dt*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6
        dpi += dt*(k1[3]+2*k2[3]+2*k3[3]+k4[3])/6
    return fr, fi

print("\nKS test-scalar, R-family (l=0, radial k), frozen late-time:")
kk = [2, 4, 8, 16, 32, 64, 128]
Pk = []
for k in kk:
    fr, fi = solve_ks(k, 0)
    amp = fr*fr + fi*fi
    Pk.append(k**3*amp)
    print(f"  k={k:5.1f}: |f|^2*k^3 = {k**3*amp:.6e}")
print(f"  fitted n_s(R-family) = {1 + np.polyfit(np.log(kk), np.log(Pk), 1)[0]:.4f}")

print("\nKS test-scalar, Omega-family (k=0, sphere l):")
ll = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48]
Pl = []
for l in ll:
    fr, fi = solve_ks(0, l)
    amp = fr*fr + fi*fi
    keff = l + 0.5
    Pl.append(keff**3*amp)
    print(f"  l={l:3d}: |f|^2*k_eff^3 = {keff**3*amp:.6e}")
print(f"  fitted n_s(Omega-family) = {1 + np.polyfit(np.log(ll), np.log(Pl), 1)[0]:.4f}")
