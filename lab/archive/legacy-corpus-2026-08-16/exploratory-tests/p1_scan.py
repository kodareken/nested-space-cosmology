import numpy as np

# P1 fertility scan: semi-analytic scaling model anchored on Adams (2008)
# constants varied: a = alpha/alpha_obs, r = (m_e/m_p)/obs, g = G/G_obs, l = Lambda/Lambda_obs
# observed point: (1,1,1,1)

# --- stellar model ---
def stellar(a, r, g):
    # hydrogen-burning minimum (degeneracy vs fusion): Adams M_min ~ a^-3/2 r^3/4 g^-3/2
    M_min = 0.08 * a**-1.5 * r**0.75 * g**-1.5
    # maximum star (radiation pressure / PISN): ~150 Msun, mild dependence on g
    M_max = 150.0 * g**-0.25
    # CO core mass of star M: calibrated so M=25 Msun -> 2.2 Msun core at g=1
    def core(M):
        return 2.2 * (M/25.0)**1.4
    # NS max mass ~ Chandrasekhar: M_NS = 2.2 g^-3/2
    M_NS = 2.2 * g**-1.5
    # BH threshold: core > M_NS  ->  M_BH = 25 g^(-1.5/1.4)
    M_BH = 25.0 * g**(-1.5/1.4)
    # PISN core threshold ~ m_e dependent: 40 r^-1
    M_PISN = 40.0 * r**-1.0
    # main-sequence lifetime: t_MS(M) = 10 Gyr a^-1 (M/1)^-2.5
    def t_MS(M):
        return 10.0 * a**-1.0 * M**-2.5
    return M_min, M_max, core, M_NS, M_BH, M_PISN, t_MS

def fertility(a, r, g, l):
    M_min, M_max, core, M_NS, M_BH, M_PISN, t_MS = stellar(a, r, g)
    if M_min >= M_max:
        return 0.0
    # --- cosmology module ---
    t_L = 43.0 / np.sqrt(l)                 # Gyr; Lambda-dominated age
    t_form = 0.15 * g**-0.5                 # Gyr; first-star timescale ~ t_ff
    if t_L <= t_form:
        return 0.0
    # BH formation requires the massive stars to die before t_L
    if t_MS(M_BH) > (t_L - t_form):
        return 0.0
    # baryon mass available: M_b ~ rho_Lambda (c t_L)^3 ~ l^-1/2
    M_b = 1.0 * l**-0.5
    N_star = M_b / M_min                    # max star number (all baryons in stars)
    # IMF: Salpeter M^-2.35 over [M_min, M_max], normalized
    # fraction of stellar mass in stars that die as BH:
    # f = int_{M_BH}^{M_up} M^-1.35 dM / int_{M_min}^{M_max} M^-1.35 dM
    # BH for core in (M_NS, M_PISN): M in (M_BH, M_PISN_star) where core(M_PISN_star)=M_PISN
    M_PISN_star = 25.0 * (M_PISN/2.2)**(1.0/1.4)
    M_BH_up = min(M_PISN_star, M_max)
    if M_BH >= M_BH_up:
        return 0.0
    p = -0.35
    f_BH = (M_BH_up**p - M_BH**p) / (M_max**p - M_min**p)
    if f_BH <= 0:
        return 0.0
    # stellar conversion efficiency + lifetime window
    eff = min(1.0, (t_L - t_form)/t_form)
    return N_star * f_BH * eff

# --- calibration: our universe ---
F_obs = fertility(1.0, 1.0, 1.0, 1.0)
print(f"F_obs = {F_obs:.4e}  (relative units; N_star~M_b/M_min, f_BH fraction)")

# --- 4D scan ---
grid = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0]
N = len(grid)
F = np.zeros((N, N, N, N))
for i, a in enumerate(grid):
    for j, r in enumerate(grid):
        for k, g in enumerate(grid):
            for m, l in enumerate(grid):
                F[i, j, k, m] = fertility(a, r, g, l)
print(f"grid {N}^4 = {N**4} points;  F range: [{F.min():.3e}, {F.max():.3e}]")
n_nonzero = np.count_nonzero(F)
print(f"nonzero points: {n_nonzero} ({100*n_nonzero/N**4:.1f}%)")

# --- local maxima ---
maxima = []
for i in range(1, N-1):
    for j in range(1, N-1):
        for k in range(1, N-1):
            for m in range(1, N-1):
                v = F[i, j, k, m]
                if v <= 0:
                    continue
                nb = [F[i+1,j,k,m], F[i-1,j,k,m], F[i,j+1,k,m], F[i,j-1,k,m],
                      F[i,j,k+1,m], F[i,j,k-1,m], F[i,j,k,m+1], F[i,j,k,m-1]]
                if v >= max(nb):
                    maxima.append((v, grid[i], grid[j], grid[k], grid[m]))
maxima.sort(reverse=True)
print(f"\nlocal maxima found: {len(maxima)}")
for v, a, r, g, l in maxima[:10]:
    print(f"  F={v:.3e} at (a,r,g,l)=({a:.1f},{r:.1f},{g:.1f},{l:.1f})")

# --- p-value of observed point ---
pval = np.count_nonzero(F >= F_obs) / N**4
print(f"\np-value: fraction of landscape with F >= F_obs = {pval:.4f}")

# --- local structure around (1,1,1,1) ---
i0 = grid.index(1.0)
print(f"\nlocal structure around observed point (a,r,g,l)=(1,1,1,1):")
print(f"  F values: {F[i0,i0,i0,i0]:.3e} (obs)")
for d, name in [(1,'a'), (2,'r'), (3,'g'), (4,'l')]:
    ii = [i0]*4; ii[d] = i0+1; Fp = F[tuple(ii)]
    ii = [i0]*4; ii[d] = i0-1; Fm = F[tuple(ii)]
    print(f"  {name}: up={Fp:.3e}, down={Fm:.3e}")
# is observed point a local max?
nb = [F[i0+1,i0,i0,i0], F[i0-1,i0,i0,i0], F[i0,i0+1,i0,i0], F[i0,i0-1,i0,i0],
      F[i0,i0,i0+1,i0], F[i0,i0,i0-1,i0], F[i0,i0,i0,i0+1], F[i0,i0,i0,i0-1]]
print(f"  neighbors: {['%.3e'%v for v in nb]}")
print(f"  observed is local max: {F_obs >= max(nb)}")
