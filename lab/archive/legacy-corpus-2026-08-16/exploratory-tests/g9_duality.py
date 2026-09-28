import numpy as np

# ============ CHAIN NUMBERS ============
# per-level horizon ratio q = 1/sqrt(Omega_Lambda), verified at our level (-0.00 sigma)
Om = 0.6889; sOm = 0.0056
q = 1.0/np.sqrt(Om)
sq = 0.5*Om**-1.5*sOm
# descendant stack (children, grandchildren, ...): geometric series, ratio 1/q < 1
Sd = q/(q-1.0)
sSd = sq/(q-1.0)**2
# ancestor pull-back stack (mother, grandmother, ...): 1/q + 1/q^2 + ...
Sa = 1.0/(q-1.0)
sSa = sq/(q-1.0)**2
print("=== The nested chain ===")
print(f"per-level ratio q = 1/sqrt(Omega_L) = {q:.5f} +/- {sq:.5f}")
print(f"descendant stack   S_d = q/(q-1)      = {Sd:.4f} +/- {sSd:.4f}   (children down: CONVERGENT)")
print(f"ancestor pull-back S_a = 1/(q-1)      = {Sa:.4f} +/- {sSa:.4f}   (mothers up: finite pull-back sum)")
print(f"maturity factor pi/2 * 1/sqrt(Omega_L) = {np.pi/2*q:.4f}")

# ============ G9: THE COLLAPSE RESIDUAL (neutron star -> black hole) ============
G = 6.6743e-11; c = 2.99792458e8; Msun = 1.98892e30
print("\n=== G9: the two-sided event, neutron-star collapse ===")
for M in [1.4, 2.5]:
    Mkg = M*Msun
    rs = 2*G*Mkg/c**2
    R = 12.0e3                       # neutron-star radius
    t_light = R/c
    t_ff = (np.pi/2)*np.sqrt(R**3/(2*G*Mkg))     # free-fall from R to center
    t_over = (R-rs)/c                           # overhang crossing at light speed (upper-bound speed)
    t_ff2rs = (np.pi/2)*np.sqrt(R**3/(8*G*Mkg)) # free-fall to r=rs (standard t_ff to center/2... approximate)
    print(f"M={M} Msun: r_s={rs/1e3:.2f} km, R={R/1e3:.0f} km, shrink 12km->{rs/1e3:.1f}km (overhang { (R-rs)/1e3:.1f} km)")
    print(f"  t_light = {t_light*1e6:.0f} us   t_ff(center) = {t_ff*1e3:.3f} ms   t_overhang/c = {t_over*1e6:.0f} us")
print("\nObserved (literature): dynamical NS->BH collapse ~ 0.1-10 ms;")
print("GW170817 post-merger remnant collapsed ~1.74 s later (hypermassive NS survival, differential rotation)")

# residual analysis vs chain candidates
print("\n=== Residual analysis ===")
tff_14 = (np.pi/2)*np.sqrt((12e3)**3/(2*G*1.4*Msun))
cands = {"1/sqrt(Om)": q, "S_d (desc. stack)": Sd, "S_a (anc. pull-back)": Sa,
         "maturity pi/2q": np.pi/2*q, "2pi": 2*np.pi, "e": np.e}
obs_ms = 1.0  # representative observed dynamical collapse, ms
res = obs_ms*1e-3/tff_14
print(f"naive free-fall: {tff_14*1e3:.3f} ms;  observed ~{obs_ms} ms -> residual factor {res:.2f}")
print("nearest chain candidates:")
for n, v in cands.items():
    print(f"  {n:22s} = {v:.3f}   (ratio residual/{v} = {res/v:.3f})")

# ============ DEFLATION-INFLATION DUALITY ============
print("\n=== Deflation-Inflation duality (one event, two clocks) ===")
Bb = 0.69e-15          # bounce radius (P4: 0.69 fm)
Tbb = 1.152e32         # K, bounce temperature (P4)
t_throat = Bb/c
print(f"child-side throat crossing time: B_b/c = {t_throat:.2e} s  (T_bb = {Tbb:.2e} K)")
print(f"parent-side deflation: 12 km -> {2*G*1.4*Msun/c**2/1e3:.1f} km in {tff_14*1e3:.2f} ms")
print(f"ratio t_deflation/t_inflation = {tff_14/t_throat:.3e}")
print(f"ratio R_NS/B_b = {12e3/Bb:.3e}")
print(f"locking check: (t_def/t_infl)/(R_NS/B_b) = {(tff_14/t_throat)/(12e3/Bb):.2f}  (O(1) = same geometry)")
RH = c/(67.4e3/3.0857e22)
print(f"\nchild total expansion factor: R_H/B_b = {RH/Bb:.2e}")
print(f"parent-side total compression factor: R_NS/B_b = {12e3/Bb:.2e}")
print(f"asymmetry (expansion/compression) = {RH/12e3:.2e}  <- the gradient/entropy transfer")
print(f"\nS_d * R_H = {Sd:.2f} * R_H  = finite total length of the descendant string")
print(f"ancestor tower: DIVERGENT (each level 1.2048x larger) - the 'Infinity' of BlackHoles-Infinity")
