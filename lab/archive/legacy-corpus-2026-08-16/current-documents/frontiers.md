# Frontiers: The Open Problems, Merged and Current

*One document, current statuses. This replaces the former `unresolved-areas.md`, `deep-dive-unresolved.md`, and the earlier `frontiers.md` (all archived). Every open problem with its executed result, its surviving route, and its referee.*

---

## 1. The Bounce Mechanism

**Status: verified numerically (P4), mechanism identified.**

The three-candidate comparison, updated:

| Candidate | Bounce density | Radius | Status |
|---|---|---|---|
| Loop quantum gravity | ~0.41 ρ_Pl | ~2.3 fm | Consistent, independent mechanism |
| Einstein–Cartan (spin-torsion) | (0.7–15) ρ_Pl | 0.7–2.3 fm | **Verified end-to-end** (RK4 integration: bounce at a_min exactly, Friedmann consistency 10⁻¹⁴, focusing reversed, SEC violated as required, thermal content T_bb = 1.152×10³² K reproducing Popławski) |
| Asymptotic safety | ~ρ_Pl | fm-scale | Open |

The two independent mechanisms converge on the same fm-scale bounce. The retired branches (macroscopic bounce, 118 km, 1.6×10³⁷ kg/m³) are documented in the-action.md §6. The open piece: the full 3D numerical-relativity treatment (the 1D RK4 is done).

## 2. Dark Matter and Baryogenesis

**Status: the torsion filter is the mechanism; the candidate priority: sterile neutrino first, neutralino second.**

- **The torsion filter** (spin-½ survives, spin-0 suppressed): alive; ATLAS/CMS 2024–2025: squarks ≳1.6–1.85 TeV, gluinos ≳2.0–2.3 TeV excluded, all null, exactly as predicted. Kill clause: any squark at HL-LHC/FCC.
- **BBN safety:** the bounce's torsion corrections are suppressed by (ρ_BBN/ρ_bounce) ~ 10⁻²⁹; ΔN_eff < 0.2; clean.
- **νMSM module:** the sterile neutrino (~3.55 keV; X-ray line test with XRISM/Athena) and resonant leptogenesis via the EC axial-current/torsion Sakharov chain — the parent's spin drives the CP violation. Note: the old "a* = 0.7 produces η_B naturally" statement is updated — the derived parent spin is **a* = 0.39 ± 0.04** (the-consistency-audit §3b); the baryogenesis parameters inherit this.
- **Torsion remnants** (Skyrmion-like, ~10²⁰ kg): speculative, retained as a candidate.

## 3. Gravitational-Wave Echoes: The Inverted Gate

**Status: the null IS the prediction.**

The bounce wall is femtometer-scale; the echo wavelength-to-wall ratio is λ_echo/r_b ≈ 10²⁰ — echoes are unobservable at all masses. LIGO–Virgo–KAGRA's null results are the expected outcome; a strong echo detection at any mass would falsify the bounce scale. The equal-frequency-spacing template remains as the discriminator if anything is ever found. (The old "echoes at f ∝ 1/M, testable" framing is retired.)

## 4. Constant Inheritance and Fertility (P1)

**Status: EXECUTED — UNSUPPORTED as computed. The MESA-class scan is the open door.**

The 12⁴ = 20,736-point calibrated scan (Adams-anchored stellar scalings: M_min ∝ α^(−3/2)(m_e/m_p)^(3/4)G^(−3/2), Chandrasekhar scaling, CO-core mapping, PISN window, Salpeter IMF, Λ-lifetime and star-formation gates) over (α, m_e/m_p, G, Λ):

- The observed constants sit at fertility fraction **0.64** (below-median), not at a local maximum.
- No local maximum exists at (1,1,1,1); the global max sits at the landscape corner.
- The old mutation-amplification claims (5.6×10²⁹, ~700-generation convergence, "attractors") are retired — the amplification at the self-consistent bounce is O(1).

The surviving route: the MESA-class scan with nuclear-stability and cooling-chemistry gates — the toy model's absence of those gates is the honest limitation of the negative. The simulation *machinery* (2D log-space landscape, per-generation mutation loop, sensitivity analysis) remains as a reusable scaffold for any future mutation law.

## 5. The Spin Pillar

**Status: the CMB axis family persists (a posteriori); the galaxy-spin support is downgraded.**

- **CMB anomalies (PR4/NPIPE):** all four persist — quadrupole–octupole alignment, hemispherical asymmetry, parity asymmetry, cold spot; individual 2–2.5σ, a posteriori. The combined 6.5σ headline is retired.
- **The executed quadrupolar-tilt test (Planck SMICA, healpy pipeline):** the pre-registered prediction (A_Q/(1−n_s) = a*², axis in the anomaly family) — verdict **NOT CONFIRMED** (measured ratio 0.154; marginal ~2σ; the raw-sky quadrupole tracked the galactic mask — identified and mitigated). LiteBIRD/CMB-S4 polarization is the referee.
- **The derived spin:** a* = 0.39 ± 0.04 from the measured tilt ratio — consistent with the spin ladder (the parent cannot spin like the objects inside it; spin-per-mass drops ~10²² up the chain; the two signatures — tilt quadrupole and hemispherical asymmetry — agree within the scale-dependent factor).
- **Galactic spin alignment:** downgraded — the local Taylor & Jagannathan (2016) alignment is contested (2023 reanalysis: consistent with isotropy); the global radio-axis directions point toward Virgo, not the CMB family. The tidal-torque amplification arithmetic (10³–3.7×10⁴) and the p > 2.5 isotropy bound remain the quantitative frame: the galactic signal is likely undetectable; the CMB is the spin's observable.

## 6. The Primordial Spectrum (P3)

**Status: components computed with a certified pipeline; the two-region junction is the open computation.**

- Pipeline certified: de Sitter benchmark n_s = 1.0002 vs the exact 1.000000.
- KS patch alone: blue (n_s = 2.77–3.0, direction-resolved, converged) — the observable modes crossed in the parent's outer dust collapse, not in the patch.
- Two-region structure: dust crossing (n_s = 1 baseline) → KS-throat transmission → child leaf.
- The junction mechanism: n_s = 1 − 12.4(σ/θ)²; the observed 0.9649 ⟺ σ/θ ≈ 5.3% at crossing.
- The amplitude anchor: A_s = 2.1×10⁻⁹ ⟺ S_eff = 4.76×10⁸ at 1.23×10⁴ ℓ_P — an integrated quantity no single horizon supplies; the open target.
- Open finding: the measured n_s–A_s sky correlation (r = −0.35, −7.3σ) with a sign opposite the naive prediction — either the junction's convention or a window systematic; the full-likelihood treatment decides.

## 7. The Gates, Current

| Gate | Status |
|---|---|
| G1/G2 — w = −1 | ALIVE, stressed: DESI DR2 2.8–4.2σ, below the 5σ kill clause; Euclid Oct 2026 is the referee |
| G3/G7 — echo null | ALIVE (inverted gate; null expected) |
| G4/G5 — CMB anomaly persistence | ALIVE, instrumented (LiteBIRD/CMB-S4) |
| G6 — 3.55 keV line | ALIVE (XRISM/Athena) |
| G8 — no squarks | ALIVE (all nulls to date) |
| G9 — neutron-star balance | standard physics anchor |
| P1 — fertility | EXECUTED, UNSUPPORTED as computed; MESA-class open |
| P2 — chain simulation | scaffold only; mutation Lagrangian open |
| P3 — spectrum | two-region junction open; anchors fixed |
| P4 — EC collapse | 1D verified; 3D open |
| P5 — wall cross-check | heuristic only; full semiclassical derivation open |

---

*The archives: the superseded documents live in `archive/` — the pre-retraction books, the six-proofs structure, the process reports, the early analyses. The current core is the index's reading order.*
