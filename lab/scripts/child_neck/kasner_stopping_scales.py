"""Where does the stretch stop? Three scales for the shear-produced pairs in a Kasner collapse.

Inputs: kasner_shear_production.json (C = n tau^3 and eps tau^4 per u_c, from the exact two-level reduction),
eck_collapse_consistency_check.json (the count N_need that the paper's N ~ a^-3 assumption requires).
Units: x = c tau / l_P; energy densities in eps_P = hbar c / l_P^4; hbar = c = 1 inside the formulas.

Test-field pair gas:   n_pair = C / tau^3,   eps = e / tau^4  (e = eps tau^4 coefficient).
Thread (Schwarzschild interior near the crunch, Kasner):  V = 6 pi r_g^2 (c tau) * L,  L = thread length / (2M).
  N(x)         = 6 pi C L (r_g / l_P)^2 / x^2
  shear scale  : theta^2 / 3 = sigma^2 = 1 / (3 tau^2) for vacuum Kasner; backreaction when kappa eps = sigma^2:
                 kappa = 8 pi l_P^2  ->  x_br^2 = 24 pi e
  torsion      : eps_eff = eps - alpha n_f^2, alpha = kappa (hbar c)^2 / 32, n_f = 2 n_pair (fermions + antifermions)
                 -> alpha n_f^2 / eps = (x_t / x)^2,  x_t^2 = pi C^2 / e
  count scale  : N(x_N) = N_need
Stage 2 (isotropised radiation after x_br; N conserved; eps ~ a^-4, n ~ a^-3):
  a_b / a_br = x_t / x_br,  eps_b = eps_br (x_br / x_t)^4,  E_b = E_br x_br / x_t.
All of it is order-of-magnitude: single Kasner epoch, thread length set by hand, counting cutoff u_c in {0.5, 1, 2}.
"""
import json
import math

G = 6.67430e-11
c = 2.99792458e8
hbar = 1.054571817e-34
Msun = 1.98847e30
mn = 1.67492750e-27
lP = math.sqrt(hbar * G / c**3)
epsP = hbar * c / lP**4
M_B_PAPER_SUN = 1.7e31        # recorded earlier from Eq. 9 scaling (prefactor not re-derived here)


def scales(m_bh_sun=10, species=1.0, L=1.0, uc="1.0"):
    k = json.load(open("kasner_shear_production.json"))
    eck = json.load(open("eck_collapse_consistency_check.json"))[str(m_bh_sun)]
    N_need = eck["N_needed_at_bounce"]
    rg = 2 * G * m_bh_sun * Msun / c**2
    R = rg / lP
    S_BH = math.pi * R**2
    C = species * k["C_pairs_density_coeff"][uc]
    e = species * k["eps_tau4_coeff"][uc]
    x_br = math.sqrt(24 * math.pi * e)
    x_t = C * math.sqrt(math.pi / e)
    x_N = math.sqrt(6 * math.pi * C * L * R**2 / N_need)
    N = lambda x: 6 * math.pi * C * L * R**2 / x**2
    V = lambda x: 6 * math.pi * L * rg**2 * x * lP
    E_br = (e / x_br**4) * epsP * V(x_br)
    out = {
        "M_BH_solar": m_bh_sun, "species_factor": species, "thread_length_factor": L, "u_c": uc,
        "C": C, "e": e, "pair_energy_times_tau": e / C,
        "x_br": x_br, "x_t": x_t, "x_N": x_N,
        "N_at_planck_length": N(1.0), "N_at_x_br": N(x_br), "N_need": N_need, "S_BH": S_BH,
        "N_br_over_N_need": N(x_br) / N_need, "N_br_over_S_BH": N(x_br) / S_BH,
        "N_need_over_S_BH": N_need / S_BH,
        "thread_length_factor_needed": N_need / N(x_br) * L,
        "eps_br_over_epsP": e / x_br**4,
        "E_br_solar": E_br / (Msun * c**2),
        "E_stage1_paper_solar": N_need * mn / Msun,
    }
    if x_t < x_br:
        out.update({
            "a_b_over_a_br": x_t / x_br,
            "eps_b_over_epsP": (e / x_br**4) * (x_br / x_t) ** 4,
            "E_b_solar": E_br / (Msun * c**2) * x_br / x_t,
        })
        out["E_b_over_M_b_paper"] = out["E_b_solar"] / M_B_PAPER_SUN
        out["E_b_over_E_stage1_paper"] = out["E_b_solar"] / out["E_stage1_paper_solar"]
    return out


if __name__ == "__main__":
    res = {}
    for sp in (1.0, 22.5):
        for uc in ("0.5", "1.0", "2.0"):
            res[f"species{sp}_uc{uc}"] = scales(species=sp, uc=uc)
    # scaling in M_BH at fixed u_c = 1 (x scales do not depend on M_BH)
    res["M_scaling"] = {str(m): scales(m_bh_sun=m) for m in (3, 10, 30)}
    json.dump(res, open("kasner_stopping_scales.json", "w"), indent=1)
    for key, r in res.items():
        if key == "M_scaling":
            continue
        print(key, {k: (float(f"{v:.3g}") if isinstance(v, float) else v) for k, v in r.items()
                    if k in ("x_br", "x_t", "x_N", "N_at_x_br", "N_br_over_N_need", "N_br_over_S_BH",
                             "thread_length_factor_needed", "E_br_solar", "eps_b_over_epsP", "E_b_solar",
                             "E_b_over_M_b_paper", "E_b_over_E_stage1_paper")})
