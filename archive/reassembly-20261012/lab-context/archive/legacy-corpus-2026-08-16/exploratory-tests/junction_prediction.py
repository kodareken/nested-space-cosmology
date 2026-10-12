import numpy as np

# ============================================================
# J2: the junction prediction -- window-averaged quadrupolar tilt
# n_s(k,n) = 1 - 12.4 (sigma/theta)^2 (k/k_piv)^{2 gamma} (1 + beta P2(n.z))
# Kerr horizon quadrupole: beta ~ a*^2  (Q = a*^2 M^3)
# Sky-average normalization: 12.4 s^2 <(k/k_piv)^{2g}> = 0.0351 (the observed red shift)
# PREDICTION: quadrupole amplitude A_Q = 0.0351 * beta * (P2 normalization)
# ============================================================
ns_obs = 0.9649
red = 1.0 - ns_obs                 # 0.0351: the sky-averaged red shift

print("="*78)
print("J2: PRE-REGISTERED JUNCTION PREDICTION (computed before touching the maps)")
print("="*78)
print(f"sky-averaged red shift: 1 - n_s = {red:.4f}")
print("\nQuadrupole amplitude prediction A_Q = 0.0351 * beta,  beta = a*^2:")
for a_star in [0.5, 0.7, 0.9]:
    beta = a_star**2
    A_Q = red*beta
    print(f"  a* = {a_star}: beta = {beta:.3f},  A_Q = {A_Q:.4f}")

print("""
The pre-registered statement (frozen before the map analysis):

  1. The local spectral index across the sky carries a QUADRUPOLAR pattern
     d_ns(n) = A_Q * P2(n . z_Q),  with A_Q in [0.007, 0.018]
     (for the parent spin a* in [0.5, 0.9], quadrupole = 0.0351 * a*^2).
  2. The pattern's axis z_Q aligns with the Planck PR4 anomaly axis family
     (Axis of Evil / dipole / hemispherical-asymmetry directions; alignment < 30 deg).
  3. NOT a dipole: the published PR4 dipole test (PTE 5.8%, amplitude ~0.017)
     neither confirms nor rules out this prediction -- the quadrupole was never
     measured. That is what the map test does.
  4. Failure conditions, recorded in advance:
     - measured quadrupole amplitude inconsistent with [0.007, 0.018]: gate fails
     - measured axis misaligned with the anomaly family beyond 30 deg: gate fails
     - null result consistent with zero: recorded as "untested at current
       sensitivity", the prediction stands for LiteBIRD/CMB-S4.
""")
