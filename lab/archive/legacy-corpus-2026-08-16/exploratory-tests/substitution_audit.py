import numpy as np

# ============================================================
# SUBSTITUTION AUDIT: all closures recomputed with DESI DR2 inputs
# (independent of the Planck combination they were derived with)
# Sources: collected-data/desi-dr2-planck-pr4-cosmology-constraints.md
# ============================================================
G = 6.67430e-11; c = 2.99792458e8; hbar = 1.054571817e-34; kB = 1.380649e-23
lp = np.sqrt(hbar*G/c**3)
kmMpc = 3.08568e19   # km per Mpc

def H_of(v): return v/3.08568e19

print("="*78)
print("SUBSTITUTION AUDIT: closures with DESI DR2 (non-Planck) inputs")
print("="*78)

sets = {
 "DESI BAO only":          dict(Om=0.2975, sOm=0.0086, Ol=0.7025, sOl=0.0086, H=None, sH=None),
 "DESI BAO+BBN":          dict(Om=0.2977, sOm=0.0086, Ol=0.7023, sOl=0.0086, H=68.51,  sH=0.58),
 "DESI BAO+BBN+theta*":   dict(Om=0.2967, sOm=0.0045, Ol=0.7033, sOl=0.0045, H=68.45,  sH=0.47),
 "DESI BAO+CMB":          dict(Om=0.3027, sOm=0.0036, Ol=0.6973, sOl=0.0036, H=68.17,  sH=0.28),
 "DESI FS(DR1)+BAO(DR2)": dict(Om=0.3035, sOm=0.0085, Ol=0.6965, sOl=0.0085, H=68.76,  sH=0.59),
}
# Planck reference (the corpus combination)
planck = dict(Om=0.3111, sOm=0.0056, Ol=0.6889, sOl=0.0056, H=67.4, sH=0.5, t0=13.797, st0=0.023)

print("\nPlanck reference:  1/sqrt(Om_L) = {:.4f} +/- {:.4f};  T_GH/T_H = {:.4f};  t0 = {:.2f} Gyr".format(
    1/np.sqrt(planck['Ol']), 0.5*planck['sOl']/planck['Ol']**1.5, 2/np.sqrt(planck['Ol']), planck['t0']))

print("\n| Dataset | 1/sqrt(Om_L) [M_parent/M_H] | T_GH/T_H = 2/sqrt(Om_L) | t0 (Gyr) | sigma vs Planck 1.2048 |")
print("|---|---|---|---|---|")
for name, d in sets.items():
    q = 1/np.sqrt(d['Ol']); sq = 0.5*d['sOl']/d['Ol']**1.5
    tr = 2/np.sqrt(d['Ol']); str_ = d['sOl']/d['Ol']**1.5
    # age (flat LCDM, using this dataset's H0 when available)
    t0s = None
    if d['H']:
        x = np.sqrt(d['Ol']/d['Om'])
        t0s = 2/(3*H_of(d['H'])*np.sqrt(d['Ol']))*np.arcsinh(x)
        t0s /= 3.15576e16   # Gyr
    sig = abs(q - 1.204819)/(np.sqrt(sq**2 + 0.0049**2))
    t0str = f"{t0s:.3f}" if t0s else "H0-independent"
    print(f"| {name} | {q:.4f} +/- {sq:.4f} | {tr:.4f} +/- {str_:.4f} | {t0str} | {sig:.2f} |")

print("""
H0-INDEPENDENT CLOSURES (no Hubble input at all):
  M_parent/M_H = 1/sqrt(Omega_L)          -- from Omega_L alone
  T_GH/T_H     = 2/sqrt(Omega_L)          -- from Omega_L alone
  S_parent = S_dS = pi*3/(Lambda*l_P^2)   -- from Lambda alone
  -> the identities are structurally H0-free; the dataset test isolates Omega_L.

c-BOUND with DESI errors:
  BAO only:      |c_p/c_c - 1| <= 2*0.5*0.0086/0.7025 = {:.3f}% (1 sigma)
  BAO+CMB:       |c_p/c_c - 1| <= 2*0.5*0.0036/0.6973 = {:.3f}% (1 sigma)
""".format(2*0.5*0.0086/0.7025*100, 2*0.5*0.0036/0.6973*100))
