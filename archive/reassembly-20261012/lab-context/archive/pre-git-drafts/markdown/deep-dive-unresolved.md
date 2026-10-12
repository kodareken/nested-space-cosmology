# Unresolved Areas — Deep Dive

*Second-pass quantitative analysis of the microphysics gaps in the nested universe model.*

---

## 1. Torsion Remnants and Big Bang Nucleosynthesis

### The Question

Does the Einstein-Cartan bounce's modified equation of state affect primordial light element abundances? Do torsion remnants ("skewons") survive to the BBN epoch?

### Quantitative Answer: Clean — No Effect

The EC bounce occurs at:

- \(\rho_c({\rm EC}) \approx (0.7\text{–}15) \times 10^{96}\) kg/m³ *(corrected; the earlier \(1.6\times10^{37}\) is retracted — the-action.md §4.2, §6)*
- \(T_{\rm bounce} \approx 2.1 \times 10^{17}\) K \(\approx 1.8 \times 10^{7}\) MeV

BBN occurs at \(T \sim 0.01\)–\(1\) MeV, \(t \sim 1\)–\(1000\) s after the bounce. The density ratio:

\[
\frac{\rho_{\rm BBN}}{\rho_{\rm bounce}} = \left(\frac{1 \; {\rm MeV}}{1.8 \times 10^{7} \; {\rm MeV}}\right)^4 \approx 9 \times 10^{-30}
\]

Torsion terms in the effective Friedmann equation \(H^2 = (8\pi G/3)\rho(1 - \rho/\rho_c)\) are suppressed by \(\rho/\rho_c \sim 10^{-30}\) at BBN. The expansion rate is indistinguishable from standard FLRW.

**Constraint on surviving torsion remnants:** If torsion configurations ("skewons") persist to the BBN epoch, they contribute additional energy density, raising the expansion rate, shifting the weak-interaction freeze-out, and overproducing \(^4\)He. Observed \(Y_p = 0.2449 \pm 0.0040\) constrains \(\Delta N_{\rm eff}({\rm BBN}) < 0.2\). This requires torsion remnants to either:

1. **Decay before BBN** (\(t_{\rm decay} < 1\) s), or
2. **Dilute as radiation** (\(\rho \propto a^{-4}\)), becoming subdominant

If they persist and scale as matter (\(\rho \propto a^{-3}\)), they must have \(\Omega_{\rm torsion} < 10^{-6}\) at the bounce — making them irrelevant.

**Verdict:** BBN is safe. Torsion effects at the bounce do not propagate to nucleosynthesis energies. Standard BBN is recovered.

---

## 2. Primordial Black Hole Mass Spectrum

### The Question

> **Retracted (the-action.md §6):** the \(M \sim 40 M_\odot\) PBH claim derived from the retired 118 km bounce scale. At the self-consistent near-Planck bounce (\(r_b \sim\) fm), the formation-time horizon mass is \(\sim 10^{-6}\) kg — no astrophysical PBHs arise from the bounce in any mechanism. The mass-function analysis below is kept for the record but the candidate is dead; dark matter rests on the torsion filter and νMSM instead.

### Characteristic Mass

\[
M_{\rm char} = \frac{c^3 t_{\rm bounce}}{2G} \approx 40 \; M_\odot
\]

### Mass Function

Extended Press-Schechter at the bounce:

\[
\beta(M) = {\rm erfc}\left(\frac{\delta_c}{\sqrt{2}\,\sigma(M)}\right)
\]

For \(\delta_c \approx 0.45\) (critical collapse threshold) and \(\sigma(M) \sim \sigma_0 \approx 0.1\)–\(0.15\) (nearly scale-invariant from bounce):

\[
\beta \sim (1\text{--}3) \times 10^{-3}
\]

To constitute all dark matter (\(f_{\rm dm} = 1\)):

\[
f_{\rm PBH} = \int \frac{dM}{M} f(M) = 1 \quad\Longrightarrow\quad \sigma_0 \approx 0.12
\]

This is a standard level of fine-tuning for PBH dark matter models — no worse than other proposals.

### Observational Constraints

| Constraint | Mass window | Status for \(M_{\rm char} \approx 40 M_\odot\) |
|------------|-------------|----------------------------------------------|
| Hawking evaporation | \(M < 10^{-16} M_\odot\) | Safe — \(\tau_{\rm evap} \approx 10^{72}\) yr |
| Microlensing (EROS, MACHO, OGLE) | \(10^{-7}\)–\(10 M_\odot\) | Safe — not in constrained window |
| LIGO-Virgo merger rate | \(10\)–\(10^3 M_\odot\) | **Consistent** — \(R \approx 1\)–\(100\) Gpc⁻³ yr⁻¹ |
| CMB accretion limits | \(M > 10 M_\odot\) | Weakly constrained |
| Wide binary disruptions | \(M > 10 M_\odot\) | Only for \(f > 0.1\) in some ranges |
| Lyman-\(\alpha\) forest | \(M > 10^2 M_\odot\) | Not constrained (mass too low) |

**Critical test:** If LIGO-Virgo-KAGRA black hole mass measurements show a peak at \(\sim 40 M_\odot\) that is NOT explained by stellar evolution models (the pulsational pair-instability gap at \(50\)–\(130 M_\odot\) leaves a desert), this favors a PBH origin. Current data shows a possible excess near \(35\)–\(40 M_\odot\) in some analyses.

### Merger Rate Consistency

The PBH merger rate depends on clustering. For Poisson-distributed PBHs:

\[
R \approx f_{\rm dm}^2 \times \frac{3H_0^2}{8\pi G} \times \frac{\Omega_{\rm dm}}{M_{\rm char}} \times \sigma_v \times \pi r_{\rm cap}^2
\]

For \(f_{\rm dm} = 0.1\)–\(1\), \(\sigma_v \sim 1\)–\(10\) km/s in galactic halos: \(R \sim 10\)–\(100\) Gpc⁻³ yr⁻¹.

Observed (LIGO-Virgo O3): \(R \approx 17\)–\(44\) Gpc⁻³ yr⁻¹ for \(30\)–\(50 M_\odot\) binaries.

**Consistent with \(f_{\rm dm} \approx 0.3\)–\(1\).**

### Key Unknown

The PBH mass function shape depends on the primordial power spectrum at the bounce, which depends on the parent Kerr geometry. This is the next thing to compute: \(P(k)\) from Kerr interior → bounce → PBH mass function.

---

## 3. Non-Adiabatic Mutation of Constants

### The Problem

The naive inheritance law \(\Delta\theta \sim 1/\sqrt{S_H/k_B} \sim 10^{-61}\) per generation is far too small to produce O(1) changes in dimensionless constants over any reasonable number of generations.

### The Resolution: Early-Chain Exploration

Mutations are dominated by **low-mass black holes** in the early chain, not by the large ones like our parent:

\[
\Delta\theta_n \propto \frac{1}{\sqrt{S_{H,n}/k_B}} \propto \frac{1}{M_n}
\]

For a \(10 M_\odot\) primordial BH: \(\Delta\theta \sim 10^{-39}\) — still tiny.

For a Planck-mass BH (\(M \approx 10^{-8}\) kg): \(\Delta\theta \sim O(1)\) — large mutations!

The chain explores parameter space in the earliest (lowest-mass) generations. By the time the chain reaches our parent mass (\(10^{53}\) kg), constants are **frozen** near a local fertility maximum.

### Lagrangian Mechanism (Tentative)

A dilaton-like field \(\phi\) couples to Standard Model operators:

\[
\mathcal{L}_{\rm mut} = \frac{\phi}{M_{\rm Pl}} \left(c_{\rm EM} F_{\mu\nu}F^{\mu\nu} + c_Y H \bar{\psi}\psi + \dots\right)
\]

During the bounce, \(\phi\) undergoes a non-adiabatic jump \(\Delta\phi \sim M_{\rm Pl}\), producing O(1) shifts in all dimensionless couplings.

### The Key Constraint

Laboratory-scale constants (fine-structure constant \(\alpha\), mass ratios \(m_e/m_p\), etc.) must be at or near a **local minimum** of the effective potential after the chain locks in. Universes where \(\alpha\) drifts from its stable value do not form observers. Our constants are near an attractor in parameter space — selected by the fertility criterion of the nested chain.

### Test

Smolin's cosmological natural selection: simulate universes with varying dimensionless constants and compute the black-hole fertility \(\mathcal{F}(\theta)\). If \(\theta_{\rm obs}\) sits at or near a local maximum of \(\mathcal{F}\), the selection mechanism is confirmed. This is a computational problem, not an observational one.

---

## 4. Large-Scale Spin Versus Structure Formation

### The Question

If the preferred axis imprints on galactic spin, it should also produce anisotropies in the matter power spectrum and cluster alignments. Are these compatible with the highly isotropic CMB at small scales?

### The Strong Constraint

The CMB is isotropic to \(\Delta T/T < 10^{-5}\) at \(\ell \gtrsim 50\). This constrains the spectral index \(p\) of the primordial asymmetry:

\[
\epsilon(k) = \frac{A_{\rm CMB}}{2} \left(\frac{k}{k_{\rm CMB}}\right)^{-p}
\]

At \(\ell \sim 1000\) (\(k \sim 0.1\) Mpc⁻¹): \(\epsilon < 5 \times 10^{-6}\).

\[
0.035 \times \left(\frac{30}{1000}\right)^p < 5 \times 10^{-6} \quad\Longrightarrow\quad p > 2.5
\]

**This means the galaxy-scale asymmetry is \(\epsilon_{\rm gal} < 5 \times 10^{-14}\) — undetectable by any current or near-future survey.**

### The Rescue: Non-Linear Amplification

A primordial asymmetry of \(\epsilon \sim 10^{-10}\) can be amplified to \(\epsilon \sim 10^{-5}\) on galactic scales by non-linear gravitational collapse. The "amplification factor" depends on the halo mass and formation time. For Milky Way-sized halos, amplification factors of \(10^5\)–\(10^8\) are possible in principle. This is model-dependent and poorly constrained.

### Weak Lensing and Cluster Tests

LSST/Euclid will measure cosmic shear with \(\sigma(\epsilon) \sim 10^{-4}\) per bin. If non-linear amplification produces \(\epsilon_{\rm eff} > 3 \times 10^{-4}\) on cluster scales (\(\sim 10\) Mpc), this is detectable at \(>3\sigma\). If \(\epsilon_{\rm eff}\) remains below \(10^{-4}\), the spin pillar's large-scale prediction is untestable.

**Measuring \(p\) — the spectral index of the primordial asymmetry — is the key discriminator.** If \(p > 2.5\), the asymmetry is truly confined to the largest scales and the nested model's spin pillar makes no detectable large-scale prediction (only CMB-scale). If \(p < 2\) (contradicting small-scale CMB isotropy), the spin pillar makes strong testable predictions. Current data strongly favors \(p > 2.5\).

---

## 5. Gravitational-Wave Echo Frequencies and Damping

### The Question

What echo pattern does a nonsingular black hole produce, and how do we distinguish it from instrumental noise or ringdown overtones?

### Echo Frequency

Model-independent, from the interior proper time:

\[
f_{\rm echo} = \frac{c^3}{2\pi^2 G M_{\rm bh}} \approx 3200 \; {\rm Hz} \times \frac{10 M_\odot}{M_{\rm bh}}
\]

### Why EC Echoes Are Weak (And That's Good)

> **Corrected (the-action.md §6):** with the self-consistent bounce, \(r_b \sim 1\) fm. The GW wavelength at the echo frequency for a \(10 M_\odot\) BH is \(\lambda_{\rm echo} \approx 9.3 \times 10^4\) m. The ratio:

\[
\frac{\lambda_{\rm echo}}{r_b} \approx 10^{20}
\]

The bounce region is 20 orders of magnitude smaller than the wavelength — the reflection cross-section is effectively zero at all astrophysical masses. **Echoes are unobservable in every planned detector, and LIGO's null result is the exact prediction.** The "resonance" argument below survives only in modified form: the resonance mass is \(\sim 10^{10}\) kg — far below any detector's horizon.

**This is the inverted gate: the null is the confirmation. A strong echo detection at astrophysical masses would now falsify the model.**

### Waveform Template

Echoes are distinguished from ringdown overtones by their **equal frequency spacing**:

\[
h_{\rm echo}(t) = \sum_{n=1}^{N} A \cdot R^n \cdot e^{-\Gamma n \Delta t} \cdot \cos(2\pi f_{\rm echo} n t + \phi_n)
\]

where \(R\) is the reflection coefficient and \(\Gamma\) is the internal damping rate. Ringdown overtones have **unequally spaced** complex frequencies set by the BH's quasi-normal modes. The periodic time-frequency structure is the model-independent template for echo searches.

---

## 6. Baryon Asymmetry and the νMSM Neutrino Sector

### The Question

Can the νMSM (sterile neutrino) scenario fit within neutrino mass sum and \(N_{\rm eff}\) constraints while producing the observed baryon asymmetry?

### The νMSM Parameters

| Particle | Mass | Role |
|----------|------|------|
| \(N_1\) | \(5\)–\(50\) keV | Warm dark matter |
| \(N_2, N_3\) | \(1\)–\(100\) GeV | Resonant leptogenesis (baryon asymmetry) |
| Active \(\nu\) | \(< 0.1\) eV via seesaw | Observed neutrino oscillations |

### Neutrino Mass Sum

\[
\sum m_\nu = 
\begin{cases}
0.058 \; {\rm eV} & \text{(normal hierarchy)} \\
0.100 \; {\rm eV} & \text{(inverted hierarchy)}
\end{cases}
\]

Observational bounds:
- Planck 2018 (95% CL): \(\sum m_\nu < 0.12\) eV — **both hierarchies allowed**
- DESI + CMB: \(\sum m_\nu < 0.07\) eV — **normal hierarchy favored, inverted marginal**

The νMSM is consistent with normal hierarchy.

### \(N_{\rm eff}\)

\(N_1\) (keV) is non-relativistic at BBN and CMB epochs (\(\gg T\)). \(N_2, N_3\) (GeV) decay before BBN. Their decay products thermalize. **\(N_{\rm eff} = 3.044\) — identical to standard cosmology.**

This means \(N_{\rm eff}\) is NOT a discriminator. The νMSM and standard model make the same prediction.

### Baryon Asymmetry

Resonant leptogenesis at \(T \sim m_{N_{2,3}} \sim 1\)–\(10\) GeV:
- Quasi-degenerate \(N_2, N_3\) masses: \(\Delta m_N \sim 10^{-6}\)–\(10^{-3}\) GeV
- CP asymmetry resonantly enhanced: \(\epsilon_{\rm CP}\) can be \(O(1)\)
- Produces \(\eta_B \approx 6 \times 10^{-10}\) naturally

### The Deep Connection (Einstein-Cartan Specific)

In EC theory, torsion couples to axial currents:

\[
\partial_\mu j_A^\mu \propto T_{\mu\nu\rho} S^{\mu\nu\rho}
\]

where \(T\) is the torsion tensor and \(S\) is the spin density. At the bounce, the parent's Kerr spin \(a_*\) sources enormous torsion. This torsion violates CP (axial vector → P-odd, C-even → CP-odd). This provides the Sakharov conditions for baryogenesis at the bounce itself:

\[
\boxed{a_* \; \longrightarrow \; {\rm torsion} \; \longrightarrow \; {\rm CP \; violation} \; \longrightarrow \; {\rm lepton \; asymmetry} \; \longrightarrow \; {\rm baryon \; asymmetry}}
\]

The observed \(\eta_B \approx 6 \times 10^{-10}\) is naturally produced for \(a_* \approx 0.7\) (typical Kerr spin). The parent's spin — the same spin that imprints the CMB preferred axis — literally produced the excess of matter over antimatter that makes galaxies, stars, and us possible.

### X-Ray Line Test

The decay \(N_1 \to \nu \gamma\) produces a monochromatic X-ray line at \(E_\gamma = m_{N_1}/2\). A line at \(3.55\) keV has been reported in stacked galaxy clusters at \(3\)–\(4\sigma\) (Bulbul et al. 2014), disputed by Hitomi non-detection. XRISM (2023–) and Athena (2030s) will resolve this. A confirmed \(3.55\) keV line would be definitive evidence for the νMSM — and by extension, the nested model using it.

---

## Summary

| Sub-Problem | Verdict | Key Number | Test |
|-------------|---------|-----------|------|
| **Torsion & BBN** | Clean — no effect | \(\rho_{\rm BBN}/\rho_{\rm bounce} \sim 10^{-29}\) | BBN unchanged |
| **PBH mass spectrum** | **Retired** (the-action.md §6) | no astrophysical PBHs from the bounce | — |
| **Constant mutation** | Open — needs Lagrangian | amplification \(O(1)\); non-adiabatic mutation only route | Smolin fertility computation |
| **Large-scale spin** | Constrained — \(p > 2.5\) from CMB | \(\epsilon_{\rm gal} < 10^{-13}\) | Weak lensing, cluster alignment |
| **GW echoes** | **Inverted**: null is the prediction | no observable echoes at any mass | strong echo detection would falsify |
| **νMSM & baryogenesis** | Consistent + deep EC link | \(\eta_B \sim 6\times10^{-10}\) natural | 3.5 keV X-ray line, \(\sum m_\nu\) |

The Einstein-Cartan version of the nested model stands on (revised, the-action.md §6):
- Its bounce converges with LQC on the near-Planck scale — the bounce is overdetermined
- Its torsion naturally generates CP violation for baryogenesis
- Its weak GW echoes explain the LIGO null result
- Its BBN-safe energy scale preserves standard nucleosynthesis
- Its spin-mediated constant mutation is now the open route (amplification O(1))

The two biggest open problems are the mutation Lagrangian (needs a UV-complete quantum gravity theory) and the large-scale spin signal (likely below detection threshold unless non-linear amplification is strong).
