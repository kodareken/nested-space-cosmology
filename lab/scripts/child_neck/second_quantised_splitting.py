"""Particle-antiparticle level splitting from the torsion contact terms in second quantisation (mean field).

Effective action of Freidel, Minic, Takeuchi (hep-th/0507253, their Eq. 23):
    S_int = -(3/2) pi G (gamma^2/(gamma^2+1)) Int [ A.A + (2 alpha/gamma) A.V - alpha^2 V.V ],
    V^mu = psibar gamma^mu psi, A^mu = psibar gamma^mu gamma5 psi, alpha = non-minimal coupling, gamma = Immirzi parameter.
Energy density of a gas at rest (no net particle current):  H = kappa' [ h^2 - s^2 + (2 alpha/gamma) h n_V - alpha^2 n_V^2 ],
kappa' = (3 kappa / 16) gamma^2/(gamma^2+1), kappa = 8 pi G.  h = A^0 = N_R - N_L = net helicity density of all quanta
(particles and antiparticles), n_V = V^0 = N_R + N_L = net number, s = spin density (A spatial).
Mean-field single-quantum shift = dH/dn_quantum for each of the four helicity-and-charge types:
    helicity +, particle   (R field):       dN_R = +1
    helicity +, antiparticle (L field):     dN_L = -1
    helicity -, particle   (L field):       dN_L = +1
    helicity -, antiparticle (R field):     dN_R = -1
Checks the claim: the spin-spin (A.A) term shifts a particle and an antiparticle of the same physical helicity equally;
the splitting comes from the vector-axial term (needs net helicity h) and from the V.V term (needs net number n_V).
Then: eta = 34.5 (beta/2) chi_h with beta = 2 gamma alpha / (gamma^2+1), chi_h = h / n; inheritance bound chi_h <= N0 / N_f.
Run after torsion_splitting_asymmetry.py. Writes second_quantised_splitting.json.
"""
import json

import sympy as sp

x = sp.symbols("N_R N_L", real=True)
NR, NL = x
a, g, kp = sp.symbols("alpha gamma kappa_p", positive=True)
s2 = sp.symbols("s2", real=True)            # |s|^2, spin density squared; depends on both helicities and charges equally (C-even)
h = NR - NL
nV = NR + NL
H = kp * (h**2 - s2 + (2 * a / g) * h * nV - a**2 * nV**2)

shift = {
    "helicity+ particle": sp.diff(H, NR),
    "helicity+ antiparticle": -sp.diff(H, NL),
    "helicity- particle": sp.diff(H, NL),
    "helicity- antiparticle": -sp.diff(H, NR),
}
# s^2 piece: the same quantum (same physical spin) shifts it equally whether particle or antiparticle, so treat ds2 as a common term.


def parts(expr):
    return {"A.A (h^2)": sp.simplify(expr.coeff(kp, 1).subs(a, 0).subs(s2, 0)),
            "V.A cross": sp.simplify(sp.expand(expr).coeff(a, 1)),
            "V.V (alpha^2)": sp.simplify(sp.expand(expr).coeff(a, 2))}


out = {"shifts": {k: str(sp.simplify(v)) for k, v in shift.items()}}
# splitting between particle and antiparticle of the same helicity (+)
split_plus = sp.simplify(shift["helicity+ particle"] - shift["helicity+ antiparticle"])
split_minus = sp.simplify(shift["helicity- particle"] - shift["helicity- antiparticle"])
out["splitting_helicity_plus"] = str(sp.expand(split_plus))
out["splitting_helicity_minus"] = str(sp.expand(split_minus))
# set n_V = 0 (no net matter) to isolate the cross term
sp0 = sp.simplify(split_plus.subs(NL, -NR + 0))        # N_L = -N_R means n_V = 0, h = 2 N_R
out["splitting_helicity_plus_at_nV0"] = str(sp.simplify(sp.expand(sp0)))
# A.A only (alpha = 0)
out["A.A only: splitting_helicity_plus"] = str(sp.simplify(split_plus.subs(a, 0)))


def main():
    t = json.load(open("torsion_splitting_asymmetry.json"))
    k = t["eta_per_polarisation"]                        # eta = 34.5 p
    Mstar = json.load(open("child_lifetime_entropy.json"))["M_star_kg"]
    Msun = 1.98847e30
    res = {"symbolic": out, "eta_per_effective_polarisation": k}
    res["chi_h_needed_for_beta_1"] = {e: 2 * float(e) / k for e in ("6e-10", "9e-8")}
    res["inheritance_bound_chi_h"] = {str(m): Mstar / (m * Msun) for m in (3, 10, 30)}
    res["beta_needed_from_inherited_helicity"] = {
        f"eta={e}": {str(m): 2 * (float(e) / k) / (Mstar / (m * Msun)) for m in (3, 10, 30)} for e in ("6e-10", "9e-8")}
    res["parent_mass_kg_for_beta_1"] = {e: Mstar / (2 * float(e) / k) for e in ("6e-10", "9e-8")}
    res["FMT_experimental_bound_on_beta"] = 1e30
    # amplifier: V.V energy / A.A energy for a spin-averaged gas (A.A average = n^2/6 -> alpha_HD n^2 with alpha_HD = kappa/32)
    res["VV_over_AA"] = "6 (gamma^2/(gamma^2+1)) alpha^2 eta^2: about 1e-18 at eta = 1e-9, order one only for eta ~ 0.4/alpha"
    json.dump(res, open("second_quantised_splitting.json", "w"), indent=1)
    return res


if __name__ == "__main__":
    print(json.dumps(main(), indent=1))
