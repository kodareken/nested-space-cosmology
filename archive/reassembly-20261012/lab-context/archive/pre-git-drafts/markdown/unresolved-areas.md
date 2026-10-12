# Unresolved Areas: A Quantitative Tackle

*Five open problems in the nested universe model, each addressed with a concrete proposal, numerical prediction, and falsifiability path.*

---

## 1. Quantum Gravity and the Bounce Mechanism

### The Problem

The model requires collapse to halt at finite density and transition to expansion. Standard GR predicts a singularity. Which quantum gravity theory delivers the bounce, and does it conserve entropy?

### Three Candidate Mechanisms

| | LQC | Einstein-Cartan | Asymptotic Safety |
|---|---|---|---|
| **Mechanism** | Discrete area spectrum → max curvature | Fermion spin → spacetime torsion → repulsion | Running G(k) → G→0 at high energy |
| **ρ_c** | 2.1 × 10⁹⁶ kg/m³ | **(0.7–15) × 10⁹⁶ kg/m³** *(corrected; the earlier 10³⁷ is retracted — the-action.md §4.2, §6)* | ~10⁹⁵ kg/m³ |
| **r_b for M=5.6×10²² M_⊙** | 2.3 fm | **0.7–2.3 fm** | ~Planck scale |
| **Entropy conservation?** | Yes (difference eq.) | Yes (spin conserved) | Yes (RG unitary) |
| **Mathematical maturity** | High (symmetry-reduced only) | Medium | Medium |
| **Testable signatures** | Possible echoes | **Torsion → DM filter, spin alignment** | Running G in CMB |

### Recommendation: Favor Einstein-Cartan (revised)

> **Retraction:** the original case for EC rested on a "macroscopic bounce" at 118 km / 10³⁷ kg/m³. That branch violates the equation of state at its own bounce and is retired (the-action.md §6). The revised case for EC does not need a macroscopic bounce:

1. **Torsion is sourced by spin — the mechanism is built into the fermions.** No extra fields, no free parameters (the Hehl–Datta coefficient 3κ/16 is fixed).
2. **The spin-torsion bounce converges with LQC**: (0.7–15)ρ_Pl vs 0.41ρ_Pl — two independent mechanisms, one scale. The near-Planck bounce is overdetermined, whichever mechanism delivers it.
3. **The EC-specific observables are the spin-dependent ones**: the torsion filter (no spin-0 superpartners in the child — collider-falsifiable), spin-mediated baryogenesis, and the parent's preferred axis.
4. **BBN is safe in both mechanisms** (deep-dive §1); the thermal history of bounce-physics.md (T_bounce ≈ 10³² K, ~73 e-folds) was already built on the near-Planck temperature and is now internally consistent for EC as well.

The key challenge is computing the full ECSK field equations for a realistic collapse-to-bounce-to-expansion scenario. This is an active research area (Popławski 2010–2023, Gaspár & Hidalgo 2022).

### Falsifiability

If gravitational-wave echoes are detected with \(f_{\rm echo} \propto 1/M_{\rm bh}\) (see §5), and the implied bounce scale is \(\sim 100\) km rather than Planckian, ECSK is favored over LQC.

---

## 2. Dark Matter and Baryogenesis

### The Problem

ΛCDM requires cold dark matter as a separate ingredient with \(\Omega_{\rm DM}/\Omega_b \approx 5.4\). The nested model must either reproduce this or provide an alternative. Baryogenesis (matter-antimatter asymmetry) also requires a mechanism.

### Three Candidate Mechanisms

#### Candidate A: Primordial Black Holes from Bounce Fluctuations — RETRACTED

> **Retracted (the-action.md §6):** the \(M \sim 40 M_\odot\) PBH prediction was derived from the now-retired 118 km bounce scale (\(t_{\rm bounce} = r_b/c = 3.9\times10^{-4}\) s). At the self-consistent near-Planck bounce (\(r_b \sim\) fm), the horizon mass at formation is \(c^3 t_b/2G \sim 10^{-6}\) kg — no astrophysical PBHs are produced at the bounce. The LIGO 30–50 M_⊙ black hole population is therefore stellar or formed later by accretion/mergers; it is not a bounce signature. Candidate B (torsion remnants) and Candidate C (νMSM) are unaffected and remain the dark-matter mechanisms.

**Mechanism:** Density fluctuations at the bounce (\(\delta\rho/\rho\)) exceeding the collapse threshold \(\delta_c \approx 0.3\)–\(0.5\) form primordial black holes. The characteristic PBH mass is set by the horizon mass at formation:

\[
M_{\rm PBH} \approx \frac{c^3 t_{\rm bounce}}{2G}
\]

where \(t_{\rm bounce} \approx r_b/c\) is the light-crossing time of the bounce region.

*(The mechanics and table below are kept for the record; the candidate is superseded by the retraction above — at the self-consistent bounce, \(M_{\rm PBH} \sim 10^{-6}\) kg for every mechanism, so no bounce mechanism produces astrophysical dark matter PBHs.)*

| Bounce theory | \(t_{\rm bounce}\) (s) | \(M_{\rm PBH}\) (\(M_\odot\)) | Hawking lifetime (yr) | Viable DM? |
|---|---|---|---|---|
| LQC | \(7.8 \times 10^{-24}\) | \(8 \times 10^{-19}\) | \(10^{-40}\) | ✗ (evaporates instantly) |
| Einstein–Cartan (self-consistent) | \(~10^{-23}\) | \(\sim 10^{-12}\) | \(10^{-30}\) | ✗ (evaporates instantly) |

**Key result:** EC bounce produces PBHs of \(\sim 40 M_\odot\) — solar-mass to intermediate-mass. These:
- Survive Hawking evaporation (\(t_{\rm evap} \gg t_0\))
- Interact only gravitationally
- Can constitute all of dark matter if the PBH mass function is broad enough
- Are consistent with LIGO-Virgo binary black hole merger rates (many observed BHs are \(\sim 30\)–\(50 M_\odot\))

**Test:** The PBH mass function from EC bounce predicts a peak near \(40 M_\odot\). Gravitational-wave observations of BH mass distributions can test this. LIGO-Virgo-KAGRA has already observed BHs in the \(30\)–\(50 M_\odot\) range — these could BE primordial black holes from the bounce rather than stellar remnants.

#### Candidate B: Torsion Remnants (EC-Specific)

**Mechanism:** In ECSK, torsion is sourced by fermion spin density. At the bounce, spin density is maximal. As the universe expands, spin density drops, but torsion configurations can "freeze in" — analogous to magnetic flux freezing in a conducting plasma.

These torsion configurations:
- Have mass from gravitational self-energy (\(\sim 10^{20}\) kg each)
- Interact only gravitationally
- Are topologically stable (Skyrmion-like)
- Are produced only in EC bounce — a unique prediction

**Status:** Requires detailed EC field theory calculations to establish stability. Currently speculative. If confirmed, would be a "smoking gun" for EC bounce.

#### Candidate C: νMSM — Sterile Neutrino Dark Matter + Leptogenesis

**Mechanism:** The bounce triggers inflation. The inflaton decays into Standard Model particles. If right-handed neutrinos exist (required for neutrino masses via seesaw), the inflaton decay asymmetry (CP violation) produces both:
- **Baryon asymmetry** (via leptogenesis at \(T \sim 10^{10}\) GeV, \(N_2, N_3\) decays)
- **Dark matter** (sterile neutrino \(N_1\) at keV scale, produced non-thermally)

This is the νMSM (Shaposhnikov et al.), independently motivated by neutrino physics. It fits naturally into the nested model as the "microphysics module."

**Parameters:** \(m_{N_1} \sim 5\) keV, \(m_{N_{2,3}} \sim 1\)–\(100\) GeV. Mixing angle \(\theta_1 \sim 10^{-6}\) with active neutrinos.

**Test:** X-ray line at \(E = m_{N_1}/2 \approx 2.5\) keV from \(N_1 \to \nu \gamma\) decay. Marginal detections exist (3.5 keV line in galaxy clusters — controversial). Future X-ray missions (XRISM, Athena) will confirm or rule out.

### Combined Picture

The candidates are not mutually exclusive:
- **νMSM** handles baryogenesis + provides warm DM component
- **Torsion remnants** provide an additional exotic DM component (if they exist)
- ~~EC PBHs~~ — **retired** (Candidate A retraction above; the bounce produces no astrophysical PBHs at the self-consistent scale)

The observed \(\Omega_{\rm DM}/\Omega_b \approx 5.4\) could be a combination. This is an advantage: the nested model naturally produces multiple DM components from the bounce physics.

---

## 3. Dimensionless Constant Inheritance

### The Problem

If black holes spawn new universes with slightly mutated constants (Smolin 1992), what determines the mutation amplitude? The original proposal (\(\Delta\theta \propto 1/\sqrt{S_H}\)) is dimensional analysis. We need a physical mechanism.

### Proposal: Spin-Dependent Mutation

\[
\theta_{n+1} = \theta_n + \frac{A(a_*)}{\sqrt{S_{H,n}/k_B}} \; \xi
\]

where \(A(a_*) = A_0 \cdot a_*\) and \(\xi \sim \mathcal{N}(0,1)\) is a random vector in parameter space.

**Physical motivation:** The "violence" of the bounce should scale with the deviation from spherical symmetry. A Schwarzschild black hole (\(a_* = 0\)) is perfectly symmetric — the bounce is maximally "clean," producing minimal constant drift. A maximally spinning Kerr black hole (\(a_* \to 1\)) has extreme frame-dragging, an ergosphere, and a highly distorted interior — the bounce is "dirty," producing larger mutations.

**Mutation step size for our parent:**

\[
\Delta\theta \approx \frac{a_*}{\sqrt{S_H/k_B}} \approx \frac{0.7}{10^{61}} \approx 7 \times 10^{-62}
\]

This is minuscule per generation. Over many generations (\(N \sim 10^2\)–\(10^3\)):

\[
\theta_N = \theta_0 + \sum_{n=1}^{N} \Delta\theta_n \approx \theta_0 \pm \sqrt{N} \times 7\times10^{-62} \approx \theta_0 \pm 2\times10^{-60}
\]

**Problem:** The drift is too small to shift dimensionless constants by O(1) over a reasonable number of generations. If our chain is \(\sim 760\) generations deep (from \(M_{\rm Planck}\) to \(M_{\rm parent}\)), the accumulated RMS drift is only \(\sim 2\times10^{-60}\) — far too small.

**Resolution:** The mutation step must be larger than \(1/\sqrt{S_H}\). Two possibilities:

1. **Mutation dominated by low-mass, high-spin parents:** Lower-mass black holes have smaller \(S_H\) (by \(M^2\)), so their \(\Delta\theta\) is much larger. A \(10 M_\odot\) black hole has \(S_H \sim 10^{79} k_B\), giving \(\Delta\theta \sim 0.7/10^{39.5} \sim 2\times10^{-40}\) per generation — still tiny.

2. **Mutation at the bounce is non-adiabatic:** The \(1/\sqrt{S_H}\) scaling assumes the mutation is smooth. If the bounce involves a violent, non-adiabatic phase transition, the mutation could be much larger — potentially O(1) per generation in some parameters. For example, if the metric signature change at the bounce (space ↔ time inversion) induces a discontinuity in the effective action, constants could jump by amounts set by the curvature scale at the bounce.

**Revised proposal:** 

\[
A(a_*, M) = A_0 \cdot a_*^{\alpha} \cdot \left(\frac{M_{\rm Planck}}{M}\right)^{\beta}
\]

where \(\beta\) determines the mass scaling. For \(\beta = 1\): low-mass BHs produce larger mutations. For \(\beta = 1\), a \(10 M_\odot\) parent produces \(\Delta\theta \sim 0.7 \times 10^{-8}/10^{53} \sim ...\) — still tiny in absolute terms, but accumulated over many generations in the early (low-mass) chain, it could be significant.

**The computational test:** Simulate a chain of \(10^3\)–\(10^6\) generations, starting from \(M \sim M_{\rm Planck}\) with random mutations at each bounce. Track the distribution of dimensionless constants. If the distribution converges to a narrow range around values that maximize black-hole formation, Smolin's natural selection is confirmed.

---

## 4. Galactic Spin Alignment

### The Problem

The model predicts a preferred axis (parent Kerr spin). The CMB shows this axis at low \(\ell\) (3.1σ). Should galaxies also show spin alignment along this axis?

### Quantitative Prediction

The CMB hemispherical power asymmetry has amplitude \(A_{\rm CMB} \approx 0.07\) at \(\ell \sim 10\)–\(30\). This asymmetry extends to smaller scales (higher \(k\)) with a spectral index \(p\):

\[
A(k) = A_{\rm primordial} \times \left(\frac{k}{k_{\rm bounce}}\right)^{-p}
\]

For galaxy-scale modes (\(k_{\rm gal} \approx 1\) Mpc⁻¹ \(\approx 10^3 \times k_{\rm CMB, low}\)):

| \(p\) | \(A_{\rm gal}\) | Detectability with \(10^9\) galaxies |
|-------|----------------|--------------------------------------|
| 1 (nearly scale-invariant) | \(7 \times 10^{-5}\) | **4–5σ** ✓ |
| 2 | \(7 \times 10^{-8}\) | 0.04σ ✗ |
| 3 | \(7 \times 10^{-11}\) | hopeless ✗ |

The key unknown is \(p\) — the spectral index of the primordial asymmetry. Standard ΛCDM predicts \(p > 2\) (the asymmetry is only a low-\(\ell\) fluke and vanishes on small scales). The nested model is consistent with \(p \approx 1\) (the preferred axis imprints on all scales, diluted evenly by inflation).

**LSST, Euclid, and Roman will measure \(\sim 10^9\) galaxy spins. If they detect alignment at 4–5σ along the CMB Axis of Evil, this confirms the nested model's spin prediction. If they detect nothing (\(p > 1.5\) confirmed), the spin pillar is weakened.**

---

## 5. Gravitational-Wave Echoes from Nonsingular Interiors

### The Problem

If black holes bounce rather than singularly, the bounce produces a secondary gravitational-wave pulse that reflects between the horizon's potential barrier and the bounce region, creating echoes.

### Quantitative Echo Prediction

The echo delay is the round-trip proper time from horizon to bounce:

\[
\Delta t_{\rm echo} \approx 2\tau_{\rm collapse} = \frac{2\pi G M_{\rm bh}}{c^3}
\]

| Black hole mass | \(\Delta t_{\rm echo}\) | \(f_{\rm echo}\) | Detector |
|---|---|---|---|
| 10 \(M_\odot\) (stellar) | 0.31 ms | 3.2 kHz | LIGO, Virgo, KAGRA, ET |
| \(10^4 M_\odot\) (intermediate) | 0.31 s | 3.2 Hz | LIGO, ET |
| \(10^9 M_\odot\) (supermassive) | 8.6 hr | \(3.2 \times 10^{-5}\) Hz | LISA |

**Key prediction:** \(f_{\rm echo} \propto 1/M_{\rm bh}\). This scaling is model-independent — it follows from the Schwarzschild/Kerr interior geometry. If echoes are detected across multiple black hole masses with this exact scaling, it confirms the nonsingular interior model regardless of which quantum gravity theory produces the bounce.

### Current Constraints

LIGO-Virgo O1–O3 searches: **no significant echoes detected.** This places an upper bound on the echo amplitude relative to the main merger signal. Two interpretations:

1. **Echo efficiency is low:** Only a small fraction of the collapse energy goes into the echo pulse. This is plausible — the bounce is a quantum process, and the energy conversion efficiency may be \(\ll 1\%\).

2. **Echoes are absorbed:** The potential barrier near the horizon may not be perfectly reflective. In the EC theory, torsion modifies the effective potential, potentially damping echoes.

### Next-Generation Tests

- **Einstein Telescope (2030s):** 10× more sensitive than LIGO, probing \(f_{\rm echo} \sim 1\)–\(10^4\) Hz. Will detect echoes from stellar-mass BHs even if the efficiency is \(<1\%\).
- **LISA (2030s):** Probing \(f_{\rm echo} \sim 10^{-5}\)–\(10^{-1}\) Hz. Will detect echoes from supermassive BHs (\(10^4\)–\(10^7 M_\odot\)). Echo periods of minutes to hours — easily resolvable.

If ET and LISA both find no echoes across all mass ranges, the nonsingular interior model is severely constrained — either the bounce is perfectly absorbing (no echo) or the bounce doesn't happen (singularity survives).

---

## Summary: The Five Gaps and Their Status

| Gap | Best Current Proposal | Promising? | Falsifiable? | Timescale |
|-----|----------------------|-----------|-------------|-----------|
| **Bounce mechanism** | Einstein-Cartan (ECS) | Yes — macroscopic, testable | GW echoes, PBH DM mass spectrum | 5–15 yr |
| **Dark matter** | EC PBHs (~40 \(M_\odot\)) + νMSM (sterile ν, 5 keV) | Yes — EC PBHs naturally | GW BH mass distribution, X-ray line | 5–10 yr |
| **Constant inheritance** | \(A(a_*, M) \propto a_* M^{-\beta}\) | Needs computation | Smolin fertility maximum | 5–20 yr |
| **Galactic spin alignment** | \(A_{\rm gal} \sim 7\times10^{-5}\) if \(p\)≈1 | Conditional on \(p\) | LSST/Euclid/Roman | 5–10 yr |
| **GW echoes** | \(f_{\rm echo} \propto 1/M_{\rm bh}\) | Yes — model-independent | ET, LISA | 10–15 yr |

The nested model is incomplete — but so is every model of the early universe. The difference is that the nested model's gaps are *specific* — each has a quantitative proposal, a numerical prediction, and a defined experiment that will confirm or falsify it. This is how a research programme advances: not by being complete, but by being testable.
