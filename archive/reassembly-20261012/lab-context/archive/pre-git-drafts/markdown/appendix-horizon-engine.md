# Appendix: The Horizon Engine — Compactness Equations

*Standalone derivation showing that a critical-density universe with Hubble radius \(R = c/H\) naturally has compactness \(C = 1\), giving rise to relations for mass, entropy, dark energy, and cosmic time.*

---

## The Compactness Identity

Start with a sphere of radius \(R\) and average mass–energy density \(\rho\). The mass inside is:

\[
M = \rho \cdot \frac{4}{3}\pi R^3.
\]

Take \(\rho\) to be the critical density from the Friedmann equations:

\[
\rho_c = \frac{3H^2}{8\pi G},
\tag{1}
\]

where \(H\) is the Hubble parameter and \(G\) is Newton's constant.

Define the **compactness ratio**:

\[
C = \frac{2GM}{c^2 R}.
\]

Substituting \(M = \rho_c \cdot \frac{4}{3}\pi R^3\) and \(\rho_c\) from (1):

\[
C = \frac{2G}{c^2 R} \cdot \frac{3H^2}{8\pi G} \cdot \frac{4}{3}\pi R^3
= \frac{H^2 R^2}{c^2}.
\tag{2}
\]

Now choose the sphere's radius equal to the **Hubble radius**:

\[
R = R_H = \frac{c}{H}.
\tag{3}
\]

Then \(C = H^2(c/H)^2 / c^2 = 1\). **A flat critical-density region whose size is the Hubble radius has the same compactness as a Schwarzschild black hole.**

---

## Derived Relations

### Mass Scale

\[
M_H = \rho_c \cdot \frac{4}{3}\pi R_H^3
     = \frac{c^3}{2GH}.
\tag{4}
\]

For \(H_0 \approx 70 \; \mathrm{km\,s^{-1}\,Mpc^{-1}} \approx 2.27 \times 10^{-18} \; \mathrm{s^{-1}}\), this gives \(M_H \sim 10^{53} \; \mathrm{kg}\).

### Cosmic Time

\[
t \approx \frac{R_H}{c} = \frac{1}{H}.
\tag{5}
\]

For the present Hubble parameter this yields \(\sim 13.8\) billion years, matching the observed age when integrated with the full Friedmann equation (radiation + matter + dark energy).

### Horizon Entropy

\[
S_H = \frac{k_B \cdot 4\pi R_H^2}{4\ell_P^2}
     = \frac{k_B \pi R_H^2}{\ell_P^2}.
\tag{6}
\]

With \(R_H \sim 1.4 \times 10^{26} \; \mathrm{m}\) and \(\ell_P \sim 1.6 \times 10^{-35} \; \mathrm{m}\), this gives \(S_H/k_B \sim 10^{122}\) — the observable universe contains a finite information budget of staggering size.

### Dark Energy (Vacuum Energy)

\[
S_{\rm dS} = \frac{3\pi}{\Lambda \ell_P^2}.
\tag{7}
\]

If we equate the child universe's de Sitter entropy to the parent black hole's horizon entropy:

\[
\Lambda = \frac{3}{R_H^2}.
\tag{8}
\]

For the observed \(\Lambda \sim 10^{-52} \; \mathrm{m^{-2}}\), this implies a parent horizon \(R_H \sim 10^{26} \; \mathrm{m}\) — consistent with the cosmic scale. A **gargantuan parent horizon produces a tiny cosmological constant in the child**.

---

## Summary Table

| Quantity | Expression | Today's approximate value |
|----------|-----------|---------------------------|
| Hubble radius \(R_H\) | \(c/H\) | \(1.4 \times 10^{26} \; \mathrm{m}\) |
| Mass \(M_H\) | \(c^3/(2GH)\) | \(\sim 10^{53} \; \mathrm{kg}\) |
| Cosmic time \(t\) | \(1/H\) | \(13.8 \times 10^9 \; \mathrm{yr}\) |
| Horizon entropy \(S_H/k_B\) | \(\pi R_H^2 / \ell_P^2\) | \(\sim 10^{122}\) |
| Compactness \(C\) | \(2GM_H/(c^2 R_H)\) | 1 |
| Cosmological constant \(\Lambda\) | \(3/R_H^2\) (if \(S_{\rm BH} = S_{\rm dS}\)) | \(\sim 10^{-52} \; \mathrm{m^{-2}}\) |

---

*The companion interactive HTML visualization (`presentation.html`) covers the full narrative including the bounce cosmology, parent-child transition kernel, and the eternal engine cycle.*
