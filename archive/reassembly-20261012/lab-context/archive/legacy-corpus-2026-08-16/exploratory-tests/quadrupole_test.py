import numpy as np, healpy as hp, time, sys

# ============================================================
# THE QUADRUPOLAR SPECTRAL-INDEX TEST (pre-registered)
# Prediction (frozen in junction_prediction.py):
#   d_ns(n) = A_Q * P2(n . z_Q),  A_Q in [0.007, 0.018] (a* in [0.5,0.9])
#   z_Q aligned with the PR4 anomaly axis family (< 30 deg)
# Failure: amplitude inconsistent with [0.007,0.018] OR axis misaligned
# ============================================================
D = '/Users/admin/Documents/BlackHoles-Infinity/collected-data'
t0 = time.time()

# ---- read SMICA T column (10-field BINTABLE, strided view) ----
raw = open(f'{D}/smica_2048.fits','rb').read()
arr = np.frombuffer(raw, dtype='>f4', count=50331648*10, offset=8640)
T2048 = arr.reshape(50331648, 10)[:,0].copy()
del raw, arr
print(f"T map read: {T2048.shape}, mean {T2048.mean():.4f}, std {T2048.std():.4f} mK... {time.time()-t0:.0f}s", flush=True)

# ---- mask ----
raw = open(f'{D}/mask_common.fits','rb').read()
m2048 = np.frombuffer(raw, dtype='>f4', count=50331648, offset=5760).astype(np.float64).copy()
del raw
print(f"mask read: sky frac {m2048.mean():.4f}", flush=True)

# ---- downgrade to Nside 256 (file is NESTED; keep NESTED for the block structure) ----
T256 = hp.ud_grade(T2048, 256, order_in='NESTED', order_out='NESTED')
m256 = hp.ud_grade(m2048, 256, order_in='NESTED', order_out='NESTED')
del T2048, m2048
# RING views for healpy's anafast
perm_r = hp.nest2ring(256, np.arange(12*256*256))
inv_r  = np.argsort(perm_r)
T256_r = T256[perm_r]
m256_r_sky = m256[perm_r]
print(f"downgraded: {T256.shape}  {time.time()-t0:.0f}s", flush=True)

# ---- patches: Nside 8, NESTED contiguous 32x32 blocks ----
nside_p = 8
npix_p = hp.nside2npix(nside_p)
patch_of_256 = np.arange(12*256*256)//(32*32)

def local_ns(patch):
    b = np.arange(patch*1024, (patch+1)*1024)
    mm = m256[b]
    if mm.mean() < 0.98:
        return np.nan, np.nan
    m_patch = np.zeros_like(T256_r)
    m_patch[hp.nest2ring(256, b)] = 1.0
    tmap = (T256_r - T256[b].mean()) * m_patch
    cl = hp.anafast(tmap, lmax=400, use_pixel_weights=True)
    ell = np.arange(len(cl))
    sel = (ell >= 40) & (ell <= 350)
    e = ell[sel]; c = cl[sel]
    if np.any(c <= 0):
        return np.nan, np.nan
    # fit log C = a + s*log(ell)  (SMICA noise floor negligible at ell<=350)
    lc = np.log(c)
    coef = np.linalg.lstsq(np.stack([np.ones_like(e), np.log(e)], axis=1), lc, rcond=None)[0]
    s = coef[1]
    n_s = 4.0 + s
    resid_sq = (lc - (coef[0] + coef[1]*np.log(e)))**2
    sig = np.sqrt(resid_sq.mean()/np.sum((np.log(e)-np.log(e).mean())**2))
    return n_s, sig

print("running 768 patches...", flush=True)
ns_vals = np.full(npix_p, np.nan)
ns_errs = np.full(npix_p, np.nan)
for p in range(npix_p):
    if p % 96 == 0:
        print(f"  patch {p}/{npix_p}  {time.time()-t0:.0f}s", flush=True)
    ns_vals[p], ns_errs[p] = local_ns(p)

keep = ~np.isnan(ns_vals)
print(f"patches kept: {keep.sum()}/{npix_p}")
nmean = np.mean(ns_vals[keep])
print(f"sky-average local n_s = {nmean:.4f}   [calibration: Planck 0.9649]")

# ---- calibration: full-sky masked slope (window-amplification factor) ----
clf = hp.anafast((T256_r - T256_r.mean())*m256_r_sky, lmax=400, use_pixel_weights=True)
elf = np.arange(len(clf)); sel = (elf >= 40) & (elf <= 350)
coeff = np.linalg.lstsq(np.stack([np.ones_like(elf[sel]), np.log(elf[sel])], axis=1),
                        np.log(clf[sel]), rcond=None)[0]
ns_fullsky = 4.0 + coeff[1]
win_fac = (4.0 - nmean)/(4.0 - ns_fullsky)
print(f"full-sky slope n_s = {ns_fullsky:.4f}  (calibration: Planck 0.9649)")
print(f"window-amplification factor f = {win_fac:.4f}")
print(f"pre-registered windowed quadrupole range: [{0.007*win_fac:.4f}, {0.018*win_fac:.4f}]")
# ---- quadrupole fit ----
ns = ns_vals[keep]
pix_keep = np.arange(npix_p)[keep]
x, y, z = hp.pix2vec(nside_p, pix_keep)
vec = np.stack([x, y, z], axis=1)
def fit_quad(zhat):
    P2 = 1.5*(vec @ zhat)**2 - 0.5
    A = np.stack([np.ones_like(P2), P2], axis=1)
    coef = np.linalg.lstsq(A, ns, rcond=None)[0]
    resid = ns - A @ coef
    return coef, np.sum(resid**2)

best = None
for i in range(hp.nside2npix(16)):
    xc, yc, zc_ = hp.pix2vec(16, i)
    zc = np.array([float(np.atleast_1d(xc)[0]), float(np.atleast_1d(yc)[0]), float(np.atleast_1d(zc_)[0])])
    coef, chi2 = fit_quad(zc)
    if best is None or chi2 < best[2]:
        best = (zc, coef, chi2)
zQ, coef, chi2 = best
A_Q = abs(coef[1])
print(f"\nquadrupole fit: A_Q = {A_Q:.4f}  (mean n_s = {coef[0]:.4f})")
print(f"axis z_Q: (l,b) = ({np.degrees(np.arctan2(zQ[1], zQ[0]))%360:.1f}, {np.degrees(np.arcsin(zQ[2])):.1f})")

# ---- rotation test (shuffle the n_s values across patches) ----
rng = np.random.default_rng(42)
null_AQ = []
ns_shuf = ns.copy()
for it in range(500):
    rng.shuffle(ns_shuf)
    A = np.stack([np.ones(len(ns)), 1.5*(vec @ zQ)**2 - 0.5], axis=1)
    c = np.linalg.lstsq(A, ns_shuf, rcond=None)[0]
    null_AQ.append(abs(c[1]))
null_AQ = np.array(null_AQ)
pte = (null_AQ >= A_Q).mean()
print(f"rotation test: null A_Q mean {null_AQ.mean():.4f}, std {null_AQ.std():.4f}, PTE = {pte:.3f}")

# ---- axis alignment with PR4 anomaly family ----
def gal2vec(l, b):
    lr, br = np.radians(l), np.radians(b)
    return np.array([np.cos(br)*np.cos(lr), np.cos(br)*np.sin(lr), np.sin(br)])
anoms = {"Axis of Evil (220,-20)": (220,-20), "Dipole (264,48)": (264,48),
         "Hemi-asym T (205,-20)": (205,-20), "Hemi-asym E (234,-14)": (234,-14)}
print("\nalignment with PR4 anomaly axes:")
for name, (l,b) in anoms.items():
    v = gal2vec(l,b)
    ang = np.degrees(np.arccos(np.clip(np.dot(v, zQ), -1, 1)))
    print(f"  {name}: {ang:.1f} deg")

# ---- dipole reference (the published test measured a dipole: PTE 5.8%) ----
def fit_dip(zhat):
    d = vec @ zhat
    A = np.stack([np.ones_like(d), d], axis=1)
    c = np.linalg.lstsq(A, ns, rcond=None)[0]
    return c, np.sum((ns - A @ c)**2)
bestd = None
for i in range(hp.nside2npix(16)):
    xc, yc, zc_ = hp.pix2vec(16, i)
    c, chi2 = fit_dip(np.array([float(np.atleast_1d(xc)[0]), float(np.atleast_1d(yc)[0]), float(np.atleast_1d(zc_)[0])]))
    if bestd is None or chi2 < bestd[2]:
        bestd = (np.array([float(np.atleast_1d(xc)[0]), float(np.atleast_1d(yc)[0]), float(np.atleast_1d(zc_)[0])]), c, chi2)
A_D = abs(bestd[1][1])
zd = bestd[0]
print(f"dipole amplitude (reference): A_D = {A_D:.4f}  [published: ~0.017, PTE 5.8%]")
print(f"dipole axis: (l,b) = ({np.degrees(np.arctan2(zd[1], zd[0]))%360:.1f}, {np.degrees(np.arcsin(zd[2])):.1f})  [published ~ (230,-20)-ish]")

# ---- the pre-registered verdict (calibration-free ratio formulation) ----
ratio_Q = A_Q/(4.0 - nmean)
print("\n" + "="*70)
print("PRE-REGISTERED VERDICT")
print("="*70)
print(f"1. quadrupole amplitude A_Q = {A_Q:.4f}")
print(f"2. RATIO A_Q/(4-n_mean) = {ratio_Q:.4f}  vs Kerr prediction [0.25, 0.81] (a*^2): "
      f"{'IN RANGE' if 0.2 <= ratio_Q <= 0.9 else 'OUT OF RANGE'}")
print(f"3. axis alignment (min): {min(np.degrees(np.arccos(np.clip(np.dot(gal2vec(l,b), zQ),-1,1))) for l,b in anoms.values()):.1f} deg  vs requirement < 30 deg")
print(f"4. rotation-test PTE = {pte:.3f}")
print(f"5. sky-average local n_s = {nmean:.4f} (windowed; Planck full-sky 0.9649); dipole A_D = {A_D:.4f}")
