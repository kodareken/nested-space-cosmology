# Nested-Space Cosmology: The Geometry of a Universe

**Douglas Ek** · Foundation Edition (`v0.27.0`) · [Zenodo Record (`10.5281/zenodo.23890258`)](https://zenodo.org/records/23890258) · [Illustrated PDF](paper/nested-space-cosmology.pdf) · [Manuscript Source](paper/nested-space-cosmology.md)

---

## Why We Built This: The Three Broken Textbook Premises

When students spend ten years memorizing textbooks, they are asked on day one to accept three unwritten premises:
1. **That a radius can only be an unsigned positive number ($r \ge 0$).**
2. **That a circle or sphere is an infinitely smooth continuum with an infinite number of sides ($\infty$ facets).**
3. **That the universe is an isolated, closed box.**

Once you swallow those three premises, you are mathematically forced to invent patches for point singularities, repulsive dark-energy fluids, and infinite curvature.

**In reality, all three premises are false:**

1. **Signed Radius (`+ -> 0 -> -`, Axiom A2):**
   Seen edge-on along a diameter, a closed circle $z(\tau) = R_0 e^{i\omega\tau}$ is the harmonic pendulum $R(\tau) = R_0\cos(\omega\tau)$, which swings smoothly across `0` into the negative direction $R < 0$ (`-`). Restricting the radius of space to an unsigned size $r = |R(\tau)| \ge 0$ folds the negative half of the swing onto the positive half and manufactures an artificial velocity cusp $\Delta\dot r = +2R_0\omega$ and a point singularity at $R = 0$. Keeping the signed coordinate $R(\tau)$ removes the cusp ($\Delta\dot R = 0$, $\ddot R(\tau_0) = 0$), flips the 3D orientation Jacobian $\det F = (R/R_0)^3$ from `+` to `-` (`(-) x (-) = (+)` across two generations), and turns inward parent implosion ($\dot R < 0$) into outward child expansion ($\dot r_- = -\dot R > 0$) with finite bounce curvature $R_{\mathrm{Ricci}}(0) = 6\omega^2$.
2. **The Finite-Facet Sphere Between Absolute Zero and the Planck Length (Axiom A3 & Proposition 9.2b):**
   We sit in a nested room where **Absolute Zero is the floor wavelength distance** ($\lambda_{\max} = 2R_k$, a wave as wide as the room itself) and **the Planck length is the ceiling resolution limit** ($\lambda_{\min} = \ell_P = \sqrt{\hbar G/c^3}$). The spherical shape taken by particles, planets, and stars is a force-versus-force equilibrium across that finite frequency band. A sphere does **not** have an infinite number of sides: a polyhedral sphere with $F < \infty$ facets has excess isoperimetric ratio $A_F^3/(36\pi V_F^2) = 1 + \pi\sqrt{3}/F + O(F^{-2}) > 1$. A star distributes the enclosing pressure it faces evenly across the frequency range the room allows until its surface facets hit the room's Planck ceiling ($\mathcal N_{1\mathrm{D}} = \lambda_n/\ell_P = M_P/m_n = 1.29942\times 10^{19}$, $F_\star = (\lambda_n/\ell_P)^2 = 1.688\times 10^{38}$ Planck-resolution facets, $N_\star = (\lambda_n/\ell_P)^3 = 2.194\times 10^{57}$ baryons, where the geometric mean of the stellar radius and the Planck length equals the baryon wavelength $\sqrt{R_\star\ell_P} = \lambda_n$). At that exact resolution ceiling, the enclosing room catches the star, forces spaghettification ($\chi_{\mathrm{spaghetti}} = \tfrac{1}{2}\ln 2$, compactness $R_S/R_{\mathrm{NS}} = 1/2 = \theta^\star(\chi_{\mathrm{spaghetti}})$ at $R_{\mathrm{NS}} = 4\lambda_n^2/\ell_P = 10.916\ \mathrm{km}$), and collapses it inward across `0` into `-` above the Theorem 2(3) `(1,1,1)` hard-core barrier $M_{\mathrm{TOV}} \in [2.0414, 2.1454]\ M_\odot$, inflating an internal child space. **Because the Planck length $\ell_P$ is the smallest resolution step of the room ($\mathcal R(0) = R_{\mathrm{core}} = \ell_P > 0$, Proposition 1.2), a zero-volume point singularity does not exist:** the Planck ceiling is the exact pen tip where parent compression finishes and the 1D pen (**Axiom A1**) starts drawing out the new child room (`-`), deflating in frequency (stretching in wavelength) from the Planck ceiling back down toward the new room's Absolute-Zero floor $\lambda_{\max,k+1} = 2R_{k+1}$ (*setting the room again* with invariant dynamic range $\mathcal B_{k+1} = \mathcal B_k$).
3. **Open-Valve Nested Circulation Instead of a Closed Box (Axiom A4 & Theorem 1):**
   Every local room sits between an outer birth `0` valve (inflow from the parent room) and inner black-hole `0` valves (outflow into child rooms). When a constant expansion stretch density $\rho_{\mathrm{DE}} = 3P = \mathrm{const}$ is forced into a closed-box continuity equation $dE + p_{\mathrm{eff}} dV = 0$, it creates the illusion of a mysterious repulsive fluid with $w_{\mathrm{eff}} = -1$. In reality, each individual 1D filament has ordinary positive elastic stretch tension ($E_i \propto R^{-1} \propto V^{-1/3}$, $w_{\mathrm{closed}} = +1/3$), while the open `0` valve ($\partial E/\partial\theta = 0 \iff \mu_+ = \mu_-$) continuously transfers winding across the horizon ($\left.dE_{\mathrm{DE}}^\star\right\vert_{\mathrm{valve}} = +\frac{4}{3}\rho_{\mathrm{DE}} dV_-$), maintaining constant stretch density per unit volume:
   $$dE_{\mathrm{DE}}^\star = \left(-\frac{1}{3} + \frac{4}{3}\right)\rho_{\mathrm{DE}} dV_- = +\rho_{\mathrm{DE}} dV_- \quad\Rightarrow\quad w_{\mathrm{DE}} = -1.$$

---

## Five Openly Declared Geometric Axioms (`A0–A4`)

In Nested-Space Cosmology, **Axioms A0–A4 are the fundamental generating law of space and time**, while 3+1D continuum field equations and general relativity are coarse-grained local dictionary readouts (Propositions 1.1–1.2b) of what an interior observer measures before the single-horizon approximation breaks down (Part I, `NSC-P006`).

- **A0 — Path of Least Resistance:** Geometry evolves by minimizing a single intrinsic cost functional
  $$E = Q_{\mathrm{stretch}} + PV = \sum_{i=1}^3 \frac{m_i}{\lambda_i} + PV_0\lambda_1\lambda_2\lambda_3,$$
  combining the Cauchy-Schwarz 1D generator stretch cost $Q_{\mathrm{stretch}}(A) = \tfrac{1}{2}\mathrm{tr}[(A^T A)^{-1/2} C^T C]$ with enclosing compression work $PV$. No fluid incompressibility or borrowed surface tension is assumed.
- **A1 — 1D Seed + Inherited Offset = Time:** Parent spherical convergence ($\lambda_1 = \lambda_2 = e^{-\chi}$) spaghettifies infalling parcels into a 1D radial filament above $\chi_{\mathrm{spaghetti}} = \tfrac{1}{2}\ln 2$. Enclosing compression buckles the 1D seed into a helix ($v_\perp^\star = \sqrt{1 - e^{-2\chi}}$), driving repeated offset accumulation ($1\mathrm{D} \to 2\mathrm{D\ circle} \to 3\mathrm{D\ sphere}$) and generating the effective metric $g_{\mu\nu}$ (Proposition 1.1b).
- **A2 — Signed-Zero Crossing (`+ -> 0 -> -`):** When inward compression exceeds outward resistance, the signed radius $R(\tau) = R_0\cos(\omega\tau)$ crosses smoothly through `0` into $R < 0$ (`-`). In $d=3$ dimensions, $\det[(R/R_0)I_3] = (R/R_0)^3$ flips sign (`+1 -> 0 -> -1`) and parent implosion ($\dot R < 0$) becomes outward child expansion ($\dot r_- = -\dot R > 0$, $r_- = -R$) without a velocity cusp or point singularity.
- **A3 — Moving Spectral Octave (`[Planck Ceiling, Absolute-Zero Floor]`):** Each room $k$ at relative compression $\chi$ has a finite wavelength band $[\lambda_{\min}(\chi), \lambda_{\max}(\chi)] = e^{-\chi}[\lambda_{\min,0}, \lambda_{\max,0}]$ bounded by its shortest resolution wavelength before `0` crosses into `-` (Planck ceiling $\lambda_{\min} = \ell_P$) and its longest wavelength as wide as the room itself (Absolute-Zero floor $\lambda_{\max} = 2R_k$), with invariant dynamic range $\mathcal B = \lambda_{\max}/\lambda_{\min}$.
- **A4 — Two-Valve Nested Circulation:** Each `0` horizon is a one-way valve from `+` to `-`. Every local room $k$ sits between an outer birth `0` valve (inflow from parent $k-1$) and inner black-hole `0` valves (outflow to children $k+1$) on a closed vortex-diffusion ring (`O <-> o`).

---

## Deductive Core: Parameter-Free Theorems & Observational Readout (`v0.27.0`)

All seven core theorems and propositions are derived in closed form in [Section 9 of the Main Manuscript](paper/nested-space-cosmology.md) and verified to machine precision in [`nsc_atlas_offset_generator.py`](src/recursive_horizons/nsc_atlas_offset_generator.py):

| Theorem / Proposition | Exact Closed-Form Result | Physical & Observational Consequence |
| --- | --- | --- |
| **Prop. 1.1–1.2b** (*Dictionary Projections*) | $R_{\mathrm{Ricci}}(0) = 6\omega^2$, $K(0) = 12\omega^4$, $\rho_{\mathrm{core}}(0) = \frac{3\omega^2}{8\pi G}$ | Smooth orthonormal bounce through $R=0$; **A1** 3-axis stretch tensor generates effective metric $g_{\mu\nu}$ |
| **Theorem 1** (*5-Variable Cost & Open-Valve EoS*) | $\sqrt{\theta^\star} = \frac{\sqrt{16e^{2\chi}+17}-3}{4e^\chi}$, $9DVB = 4\left(V+\frac{B}{2}\right)^3$ | Exact algebraic budget surface for all $\chi \ge 0$; open-valve balance yields $w_{\mathrm{DE}} = -1$, $w_{\mathrm{DM}} = 0$, $w_b = 0$ |
| **Theorem 2** (*Exponential-Wall Bypasses & `(1,1,1)` Triad*) | $E^\star(n_1,n_2,n_3) = 4(PV_0)^{1/4}m_0^{3/4}\sqrt{n_1 n_2 n_3}$ | `(1,1,1)` ($N=3$) is the unique 3D volumetric ground state ($E^\star=4$); two `(1,1,1)` triads have attraction $-2m_0/d$ + hard-core barrier $\Delta E_{\mathrm{core}} = 8(\sqrt{2}-1) \approx +3.3137$ |
| **Prop. 9.2b** (*Finite-Facet Sphere & Stellar Collapse*) | $\sqrt{R_\star \ell_P} = \lambda_n$, $R_{\mathrm{NS}} = \frac{4\lambda_n^2}{\ell_P}$, $M_{\star,0} = \frac{M_P^3}{m_n^2}$ | $F_\star = (\lambda_n/\ell_P)^2 = 1.688\times 10^{38}$ facets; $R_{\mathrm{NS}} = 10.916\ \mathrm{km}$ ($R_S/R_{\mathrm{NS}} = 1/2 = \theta^\star(\chi_{\mathrm{spaghetti}})$); TOV black-hole collapse threshold $M_{\mathrm{TOV}} \in [2.0414, 2.1454]\ M_\odot$ |
| **Theorem 3** (*Gauss-Divergence $3:1$ Virial Split*) | $E_{\mathrm{surr}}^\star / W_{\mathrm{press}}^\star = 3$ ($75$% exterior, $25$% cavity), $\Delta E_{\mathrm{bind}}^\star = -m_0/d$ | Exact $3:1$ exterior-to-cavity partition in $d=3$; potential-flow added mass $M_{\mathrm{eff}} = \frac{3}{2}M_{\mathrm{cav}}$ clusters with baryonic cores |
| **Theorem 4** (*Surface-Area Conservation in Fission*) | $\sum_i 4\pi(R_{n_i}^\star)^2 = 4\pi(R_N^\star)^2$, $\Delta E > 0$ | Same-sign bubble fission conserves horizon area while releasing positive kinetic energy ($+29.29$% for $2\to 1+1$) |
| **Theorem 5** (*Whole-Exterior Quadrupole Cubic*) | $p^3 - 3p - 3 = 0 \Rightarrow p_+ = \phi^{2/3} + \phi^{-2/3}$ | Exact causal neck-instability growth rate $p_+ = 2.1038034027$ ($\phi = \frac{3+\sqrt{5}}{2}$) without matching-radius cutoff |
| **Theorem 6** (*Harmonicity $1\to 2\to 4\to 8$ Cascade*) | $\partial_{\zeta_x}v + \partial_{\zeta_y}v + \partial_{\zeta_z}v = 0$, ratios $-1/2$ and $-2$ | Laplace trace zero forces division sequentially along all 3 orthogonal axes $\mathbf e_z \to \mathbf e_x \to \mathbf e_y$ |
| **Theorem 7 & Prop. 9.7b** (*$\chi_{\mathrm{freeze}}=3$, WHIM, $H_0$, $S_8$*) | $\chi_{\mathrm{freeze}} = 3\chi_1 = 3$, $\Omega_{\mathrm{DE}} = 70.88\text{--}71.14$%, $\Omega_{\mathrm{DM}} = 23.81\text{--}23.88$% | Baryon baseline $\Omega_b = e^{-3} = 4.9787$% + late WHIM $\Delta\Omega_{\mathrm{WHIM}} = 0.3252$%; $H_0 = 67.03 \pm 0.22\ \mathrm{km\ s^{-1}\ Mpc^{-1}}$; $S_8 = 0.776 \pm 0.010$ via $3.114$% uncondensed horizon slip |

*(Post-hoc numerical coincidences such as $6\pi^5$, $K=7$, and ad-hoc CMB exponents are explicitly quarantined from the deductive theorems in Section 9.8 of the manuscript.)*

---

## Active Desk & Navigation

- **[Main Manuscript (`v0.27.0` Markdown)](paper/nested-space-cosmology.md)** and **[Illustrated PDF (`v0.27.0`)](paper/nested-space-cosmology.pdf)**: The complete 30-page Foundation Edition monograph and 348-page preserved archive.
- **[Zenodo Publication Record (`10.5281/zenodo.23890258`)](https://zenodo.org/records/23890258)**: Permanent DOI archive of Foundation Edition `v0.27.0`.
- **[Geometric Foundation (`THEORY.md`)](THEORY.md)**: Openly declared axioms (`A0–A4`), the three broken textbook premises, derived theorems, and quarantined coincidences.
- **[Offset & Winding Construction (`docs/nsc-atlas-offset-generator.md`)](docs/nsc-atlas-offset-generator.md)**: Detailed derivations of the 5-variable minimization, coupled parent-receiver evolution, and multipole cascade.
- **[Current Result (`docs/current-result.md`)](docs/current-result.md)**: Summary of active mathematical connections and verification status.
- **[Paper Guide (`docs/papers.md`)](docs/papers.md)** and **[Code & Verification Map (`docs/active-code-map.md`)](docs/active-code-map.md)**.
- **[Working Method (`docs/instructions.md`)](docs/instructions.md)** and **[Retained Library (`docs/research-index-history-20261102.md`)](docs/research-index-history-20261102.md#retained-library-accessed-deliberately)**.

---

## Reproducing the Calculations

Run the complete test suite and non-local-to-local atlas readouts using the isolated project runtime:

```text
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_offset_generator.py
python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --atlas-only
python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --evolution-only
python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --last-echo-only
python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --split-onset-only
python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --split-waveform-only
python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --cascade-inheritance-only
python scripts/run_python.py scripts/check_publication.py
```

Default `rg` searches show the active foundation files. Use `rg --no-ignore 'pattern' <specific-path>` to inspect retained historical work in the [research library](docs/research-index-history-20261102.md#retained-library-accessed-deliberately).

