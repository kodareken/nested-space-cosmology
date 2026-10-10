"""Which discrete maps of the Kruskal geometry preserve the Dirac action, and what they do to
charge density, axial (helicity) density, chirality and the vector-axial (V.A) torsion term.

Maps of the (T, X) plane of the Kruskal form  ds^2 = f(s)^2 (dT^2 - dX^2) - R(s)^2 dOmega^2,
s = X^2 - T^2  (f, R arbitrary even functions of s):
  Sigma_X  : (T, X) -> ( T, -X)   exchanges the two exteriors (I <-> III), interior fixed
  Sigma_T  : (T, X) -> (-T,  X)   exchanges black-hole and white-hole interiors (II <-> IV), exteriors fixed
  Sigma_XT : (T, X) -> (-T, -X)   exchanges both pairs (the wedge reflection)
For each map every spinor lift psi'(x) = S psi(Rx) (linear) or B conj(psi(Rx)) (antilinear), S, B one of the
16 gamma-matrix products, is tested numerically for  D psi' = kappa * S (D psi)(Rx).
Convention: signature (+,-,-,-), Dirac representation, (i gamma^a e_a^mu (d_mu + Gamma_mu) - m) psi = 0.
"""
import itertools, json
import numpy as np
import sympy as sp

I2 = sp.eye(2); Z2 = sp.zeros(2)
SX = sp.Matrix([[0, 1], [1, 0]]); SY = sp.Matrix([[0, -sp.I], [sp.I, 0]]); SZ = sp.Matrix([[1, 0], [0, -1]])
def blk(a, b, c, d): return sp.Matrix(sp.BlockMatrix([[a, b], [c, d]]))
G = [blk(I2, Z2, Z2, -I2)] + [blk(Z2, s, -s, Z2) for s in (SX, SY, SZ)]
ETA = sp.diag(1, -1, -1, -1)
G5 = sp.I * G[0] * G[1] * G[2] * G[3]
Gn = [np.array(g.tolist(), dtype=complex) for g in G]
G5n = np.array(G5.tolist(), dtype=complex)
ETAn = np.diag([1., -1., -1., -1.])
G0n = Gn[0]

# 16 gamma-matrix products, labelled by the subset of indices
BASIS = {}
for k in range(5):
    for sub in itertools.combinations(range(4), k):
        M = np.eye(4, dtype=complex)
        for i in sub: M = M @ Gn[i]
        BASIS["".join(str(i) for i in sub) or "1"] = M


def dirac_data(e, coords):
    """e[a, mu]: tetrad e^a_mu (sympy). Returns inverse tetrad e_a^mu and spin connection Gamma_mu."""
    g = sp.simplify(e.T * ETA * e)
    gi = g.inv()
    n = 4
    Chr = [[[sum(gi[l, k] * (sp.diff(g[k, m], coords[nn]) + sp.diff(g[k, nn], coords[m]) - sp.diff(g[m, nn], coords[k]))
                 for k in range(n)) / 2 for nn in range(n)] for m in range(n)] for l in range(n)]
    einv = e.inv()                       # einv[mu, a]  (e_a^mu)
    elow = ETA * e                       # e_{b nu} = eta_bc e^c_nu
    Gam = []
    for m in range(n):
        Mu = sp.zeros(4, 4)
        for a in range(n):
            for b in range(n):
                if a == b: continue
                w = sum(einv[nu, a] * (sp.diff(elow[b, nu], coords[m]) - sum(Chr[l][m][nu] * elow[b, l] for l in range(n)))
                        for nu in range(n))
                Mu += sp.Rational(1, 8) * w * (G[a] * G[b] - G[b] * G[a])
        Gam.append(Mu)
    return einv, Gam


def make_psi(coords, seed=1):
    rng = np.random.default_rng(seed)
    t, x, y, z = coords
    comps = []
    for k in range(4):
        c = rng.normal(size=6) + 1j * rng.normal(size=6)
        expr = (c[0] + c[1] * t + c[2] * x * y + c[3] * sp.cos(z) + c[4] * sp.sin(y) * t + c[5] * x * x) \
            * sp.exp(sp.I * (0.4 * x + 0.2 * z))
        comps.append(sp.nsimplify(0) + expr)
    return sp.Matrix(comps)


class Geometry:
    def __init__(self, e, coords, seed=1):
        self.coords = coords
        self.einv, self.Gam = dirac_data(e, coords)
        self.psi = make_psi(coords, seed)
        self.dpsi = [self.psi.diff(c) for c in coords]
        f = lambda m: sp.lambdify(coords, m, 'numpy')
        self.fpsi = f(self.psi); self.fdpsi = [f(d) for d in self.dpsi]
        self.feinv = f(self.einv); self.fGam = [f(m) for m in self.Gam]

    def ev(self, fn, p):
        return np.array(fn(*p), dtype=complex)

    def D(self, psi, dpsi, p):
        einv = self.ev(self.feinv, p).reshape(4, 4)
        Gam = [self.ev(g, p).reshape(4, 4) for g in self.fGam]
        out = np.zeros(4, dtype=complex)
        for a in range(4):
            for mu in range(4):
                if abs(einv[mu, a]) < 1e-300: continue
                out += 1j * einv[mu, a] * (Gn[a] @ (dpsi[mu] + Gam[mu] @ psi))
        return out

    def D_original(self, p):
        psi = self.ev(self.fpsi, p).reshape(4)
        dpsi = [self.ev(d, p).reshape(4) for d in self.fdpsi]
        return self.D(psi, dpsi, p)

    def test_lift(self, r, lift, anti, pts):
        """D psi'(x) = kappa * S (D psi)(Rx)?  returns (kappa, worst residual)."""
        kap = []; res = 0.0
        for p in pts:
            q = [r[i] * p[i] for i in range(4)]
            phiR = self.ev(self.fpsi, q).reshape(4)
            dphiR = [r[i] * self.ev(self.fdpsi[i], q).reshape(4) for i in range(4)]
            if anti:
                phiR = np.conj(phiR); dphiR = [np.conj(d) for d in dphiR]
            lhs = self.D(lift @ phiR, [lift @ d for d in dphiR], p)
            Dq = self.D_original(q)
            rhs = lift @ (np.conj(Dq) if anti else Dq)
            k = np.vdot(rhs, lhs) / np.vdot(rhs, rhs)
            kap.append(k); res = max(res, np.linalg.norm(lhs - k * rhs) / np.linalg.norm(lhs))
        kap = np.array(kap)
        return kap, res

    def search(self, r, pts, tol=1e-9):
        found = []
        for name, M in BASIS.items():
            for anti in (False, True):
                kap, res = self.test_lift(r, M, anti, pts)
                if res < tol and np.ptp(kap) < 1e-8:
                    k = kap[0]
                    found.append(dict(lift=name, antilinear=anti, kappa=complex(np.round(k, 8)), residual=float(res)))
        return found


def bilinear_signs(lift, anti):
    """Classical (commuting-spinor) transformation of psibar Gamma psi under psi -> lift psi (or lift conj psi):
    returns diagonal factors for V^a, A^a, S, P."""
    def eff(Gam):
        if anti:
            return G0n @ (lift.conj().T @ G0n @ Gam @ lift).T
        return G0n @ lift.conj().T @ G0n @ Gam @ lift
    def coef(X, basis, norm):
        return [np.trace(norm(b) @ X) / np.trace(norm(b) @ b) for b in basis]
    # vector
    Mv = []; Na = []
    for a in range(4):
        X = eff(Gn[a]); Y = eff(Gn[a] @ G5n)
        cv = [np.trace(Gn[b] @ X) / 4 * ETAn[b, b] for b in range(4)]
        ca = [-np.trace(Gn[b] @ G5n @ Y) / 4 * ETAn[b, b] for b in range(4)]
        # X should be c_a gamma^a only
        Mv.append(cv[a]); Na.append(ca[a])
        assert max(abs(cv[b]) for b in range(4) if b != a) < 1e-9 and max(abs(ca[b]) for b in range(4) if b != a) < 1e-9
    s = np.trace(eff(np.eye(4, dtype=complex))) / 4
    p = np.trace(G5n @ eff(1j * G5n)) / 4 / 1j * 1j  # coefficient of i*gamma5 in eff(i g5)
    p = np.trace((1j * G5n) @ eff(1j * G5n)) / np.trace((1j * G5n) @ (1j * G5n))
    return np.real_if_close(Mv), np.real_if_close(Na), np.real_if_close(s), np.real_if_close(p)


def classify(lift_name, anti, r):
    lift = BASIS[lift_name]
    Mv, Na, s, p = bilinear_signs(lift, anti)
    flipE = (r[0] * (-1 if anti else 1)) == -1          # sign of frequency after the map
    eps = -1 if flipE else 1                            # second-quantised sign for particle<->antiparticle
    Mq = eps * np.real(Mv); Nq = eps * np.real(Na)
    chir = np.linalg.norm(lift @ G5n - G5n @ lift) < 1e-9     # commutes with gamma5
    chir_sign = +1 if chir else -1
    VA = np.real(Mv * Na)  # product of classical diagonal factors, same for every a if V.A is covariant
    return dict(E_flip=bool(flipE), V_signs=[int(round(x)) for x in Mq], A_signs=[int(round(x)) for x in Nq],
                S_sign=int(round(eps * np.real(s))), P_sign=int(round(eps * np.real(p))),
                chirality_preserved=bool(chir), VA_sign=[int(round(x)) for x in VA])


# ----------------------------------------------------------------------------------------------
KMAPS = {"Sigma_X": (1, -1, 1, 1), "Sigma_T": (-1, 1, 1, 1), "Sigma_XT": (-1, -1, 1, 1)}
FLAT_EXPECT = {   # standard C, P, T, CPT parity table (psibar psi, i psibar g5 psi, V^mu, A^mu), from memory
    "P":   dict(lift="0",  anti=False, r=(1, -1, -1, -1), S=1, P=-1, V=[1, -1, -1, -1], A=[-1, 1, 1, 1]),
    "T":   dict(lift="13", anti=True,  r=(-1, 1, 1, 1),   S=1, P=-1, V=[1, -1, -1, -1], A=[1, -1, -1, -1]),
    "C":   dict(lift="2",  anti=True,  r=(1, 1, 1, 1),    S=1, P=1,  V=[-1, -1, -1, -1], A=[1, 1, 1, 1]),
    "CPT": dict(lift="0123", anti=False, r=(-1, -1, -1, -1), S=1, P=1, V=[-1, -1, -1, -1], A=[-1, -1, -1, -1]),
}


def main(out="swap_symmetry.json"):
    T, X, th, ph = sp.symbols('T X theta phi', real=True)
    coords = [T, X, th, ph]
    s = X ** 2 - T ** 2
    f = sp.exp(-s / 3) / sp.sqrt(1 + s ** 2 / 5)
    Rr = sp.sqrt(2 + s + s ** 2 / 7)
    rng = np.random.default_rng(5)
    pts = [[rng.uniform(-.7, .7), rng.uniform(-.7, .7), rng.uniform(.5, 2.6), rng.uniform(0, 6)] for _ in range(4)]
    res = {"flat_validation": {}, "kruskal": {}, "negative_control": {}}
    # flat space: the lifts reproduce the standard C, P, T, CPT table
    t, x, y, z = sp.symbols('t x y z', real=True)
    gF = Geometry(sp.eye(4), [t, x, y, z], seed=3)
    ptsF = [list(rng.uniform(-.8, .8, 4)) for _ in range(4)]
    for k, ex in FLAT_EXPECT.items():
        found = [d for d in gF.search(ex["r"], ptsF) if abs(d["kappa"] - 1) < 1e-6]
        hit = [d for d in found if d["lift"] == ex["lift"] and d["antilinear"] == ex["anti"]]
        assert hit, (k, found)
        cl = classify(ex["lift"], ex["anti"], ex["r"])
        ok = (cl["S_sign"] == ex["S"] and cl["P_sign"] == ex["P"] and cl["V_signs"] == ex["V"] and cl["A_signs"] == ex["A"])
        res["flat_validation"][k] = dict(lift=ex["lift"], antilinear=ex["anti"], classification=cl, matches_standard_table=bool(ok))
        assert ok, (k, cl)
    # Kruskal-form geometry
    geo = Geometry(sp.diag(f, f, Rr, Rr * sp.sin(th)), coords)
    for name, r in KMAPS.items():
        lifts = []
        for d in geo.search(r, pts):
            entry = dict(lift=d["lift"], antilinear=d["anti"] if "anti" in d else d["antilinear"], kappa=float(d["kappa"].real),
                         residual=d["residual"])
            if abs(d["kappa"] - 1) < 1e-6:
                entry.update(classify(d["lift"], entry["antilinear"], r))
                entry["preserves_VA"] = all(v == 1 for v in entry["VA_sign"])
            lifts.append(entry)
        res["kruskal"][name] = dict(r=list(r), lifts=lifts)
    # negative control: break evenness
    fb = f * (1 + T / 4); Rb = Rr * (1 + X / 5)
    gb = Geometry(sp.diag(fb, fb, Rb, Rb * sp.sin(th)), coords)
    for name, r in KMAPS.items():
        res["negative_control"][name] = len(gb.search(r, pts))
    return res


if __name__ == "__main__":
    import os
    res = main()
    # net baryon number per quantum if the child carries only the image of the parent's baryons
    Msun, mn = 1.98847e30, 1.67492750e-27
    N0 = 10 * Msun / mn
    N_need = json.load(open("child_energy_floor.json"))["N_need"]   # written by child_energy_floor.py
    eta = N0 / N_need
    res["derived"] = dict(
        N0_parent_baryons=N0, N_need_child_quanta=N_need, eta_if_child_holds_only_parent_baryons=eta,
        shortfall_factor_for_eta_6e_minus_10=6e-10 / eta,
        note="eta = N0/N_need is the net baryon number per quantum if the child carries only the image of the parent's N0 baryons; sign + (T-type) or - (CP/CRT-type)")
    json.dump(res, open("swap_symmetry.json", "w"), indent=1, default=float)
