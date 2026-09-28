import numpy as np

# ============================================================
# CERTIFICATION BENCHMARK: de Sitter vacuum
# v'' + (k^2 - 2/eta^2) v = 0,  z = z0/eta (a = -1/(H eta))
# EXACT ANSWER: BD mode v = e^{-ik eta}(1 - i/(k eta))/sqrt(2k)
#  -> |zeta_k|^2 = 1/(2 k^3 z0^2)  ->  P = k^3|zeta|^2 = const  ->  n_s = 1 EXACTLY
# This certifies the integration + extraction + fit procedure.
# ============================================================
eta0 = -1000.0; eta1 = -1e-3
N = 400000
eta = np.linspace(eta0, eta1, N)
h = eta[1]-eta[0]
zpp = 2.0/eta**2

def solve_dS(k):
    vr = np.cos(k*eta0)/np.sqrt(2.0*k)
    vi = np.sin(k*eta0)/np.sqrt(2.0*k)
    wr = -k*vi; wi = k*vr
    for j in range(N-1):
        f = zpp[j] - k*k
        k1r=wr; k1i=wi; k1wr=f*vr; k1wi=f*vi
        k2r=wr+h*k1wr/2; k2i=wi+h*k1wi/2; k2wr=f*(vr+h*k1r/2); k2wi=f*(vi+h*k1i/2)
        k3r=wr+h*k2wr/2; k3i=wi+h*k2wi/2; k3wr=f*(vr+h*k2r/2); k3wi=f*(vi+h*k2i/2)
        k4r=wr+h*k3wr; k4i=wi+h*k3wi; k4wr=f*(vr+h*k3r); k4wi=f*(vi+h*k3i)
        vr+=h*(k1r+2*k2r+2*k3r+k4r)/6; vi+=h*(k1i+2*k2i+2*k3i+k4i)/6
        wr+=h*(k1wr+2*k2wr+2*k3wr+k4wr)/6; wi+=h*(k1wi+2*k2wi+2*k3wi+k4wi)/6
    zeta_sq = eta1**2*(vr*vr + vi*vi)     # z = z0/eta, z0=1
    return zeta_sq

kk = [0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8]
P = []
for k in kk:
    zsq = solve_dS(k)
    P.append(k**3*zsq)
    print(f"  k={k:5.1f}: k^3|zeta|^2 = {k**3*zsq:.6e}")
ns = 1 + np.polyfit(np.log(kk), np.log(P), 1)[0]
print(f"CERTIFICATION RESULT: n_s = {ns:.6f}   [exact: 1.000000]")
print(f"pipeline certified: {'YES' if abs(ns-1.0) < 0.01 else 'NO - bug must be found'}")
