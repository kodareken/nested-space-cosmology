# BlackHoles-Infinity

*Every black hole is a universe. Our universe is a black hole. Dark energy, dark matter, the constants, the arrow of time — all follow from three established laws of physics, one assumption, and zero free parameters.*

This repository is a self-contained cosmological framework: **energy is the only fabric, thermodynamics is the only instrument.** It claims no new forces, no new particles, no new constants. It merges general relativity's horizons, quantum mechanics' quanta, and thermodynamics' walls into one chain — and tests that chain against public data with pre-registered, reproducible computations.

---

## The theory in one paragraph

Space is a gradient chamber, and gradients form pockets. When a pocket collapses, the inward geometry flow eventually equals the causal speed — that equality is a black hole horizon. Inside, the roles of space and time swap, and the collapsing region is a whole contracting universe. The collapse takes finite time (πGM/c³ — 27.45 Gyr for our parent), and at the bottom it cannot crush to zero: resolution itself has a wall (Δx ≥ √2 ℓ_P — probing smaller enlarges the wall). The collapse halts, bounces (the Einstein–Cartan spin-torsion bounce, numerically verified), and expands as a new universe, bound by the entropy of the parent horizon. Everything "dark" is inheritance: the parent's horizon geometry becomes our Λ; the parent's fermions become our dark matter; the parent's spin becomes our preferred axis. The chain repeats; the walls are warm; the engine never stops.

**Three laws + one assumption:** General Relativity (1915) · unitary quantum mechanics · horizon thermodynamics (1973–77) · the saturation boundary condition (the bounce is a one-shot unitary entropy transfer).

---

## The core results

| Result | Status |
|--------|--------|
| Λ = 3/r_s² from entropy conservation | Derived; matches Planck's Λ to the published precision |
| M_parent/M_H = 1/√Ω_Λ = 1.2048 | −0.00σ against Planck 2018 |
| T_GH/T_H = 2/√Ω_Λ = 2.4097 | 4 digits |
| S_parent = S_dS = 3.2885×10¹²² k_B | exact |
| C = 1.0000000000 (compactness) | exact |
| c_p = c_c across the bounce | ≤ 0.813% (1σ) |
| The bounce at (0.7–15)ρ_Pl, converging with LQC | numerically verified (P4) |
| w = −1 exactly | DESI tension 2.8–4.2σ, below the 5σ kill clause |
| No squarks at any energy | all ATLAS/CMS nulls to date |
| The substitution audit | all closures hold under five independent DESI datasets (0.75–2.0σ) |
| The spectrum | two-region junction in progress; the pre-registered quadrupole test: NOT CONFIRMED at current sensitivity |

**The honest ledger:** the corpus states what holds, what failed as computed (the P1 fertility scan: unsupported), what is unconfirmed (the quadrupolar tilt), and what is parked outside the proof (the ancestor infinity, the simulation reading). Every claim carries its kill condition.

---

## Repository structure

```
├── README.md                  ← you are here
├── index.md                   ← the router: reading order, statuses, retracted numbers
├── the-action.md              ← the flagship derivation (ECKS action, Λ, theorems, gates)
├── the-engine.md              ← the condensed theory + proof plan
├── empirical-validation.md    ← the mesh: every closure against public data
├── the-consistency-audit.md   ← instruments, benchmarks, substitution audit, forward tests
├── the-adversarial-ledger.md  ← every attack a physicist can land, with its status
├── why-we-see-what-we-see.md  ← the observation ledger + possibility proof + meaning
├── energy-fabric.md           ← the axioms, the entropy chain, the ladders, the pixel
├── gradient-pressure-picture.md ← the fluid/gradient language
├── holographic-chain.md       ← the Λ = 3/r_s² derivation
├── c-test.md                  ← the causal-speed conservation
├── frontiers.md               ← the open problems, merged and current
├── finishing-blow.md          ← the anomaly scorecard across models
├── einstein-epilogue.md       ← the historical and philosophical close
├── argument-plan.md           ← the presentation doctrine and attack map
├── tests/                     ← 13 reproducible scripts (see tests/README.md)
├── collected-data/            ← public data (see its README)
└── archive/                   ← superseded material — never quoted
```

## Reproducing the computations

```bash
cd tests
python3 consistency_audit.py      # numpy only
python3 dS_benchmark.py           # numpy only
python3 substitution_audit.py     # numpy only
# the forward test (needs healpy + the Planck maps):
python3 -m venv cmbvenv && cmbvenv/bin/pip install healpy numpy
cmbvenv/bin/python quadrupole_test.py
```

## Reading for reviewers

- **Physicist, maximum rigor:** `the-action.md` → `empirical-validation.md` → `the-consistency-audit.md` → `the-adversarial-ledger.md`.
- **General reader:** `the-engine.md` → `gradient-pressure-picture.md` → `einstein-epilogue.md`.

## Retracted numbers (do not quote)

The 1.6×10³⁷ kg/m³ / 118 km bounce · the 40 M_⊙ PBHs · the 13 M_⊙ echoes · the 5.6×10²⁹ mutation amplification · the combined 6.5σ anomaly headline · the 2–3σ galaxy-spin claim. All retired with reasons in `index.md`.

---

*The model does not ask to be believed. It asks to be killed — and says exactly how: w ≠ −1 at 5σ, a squark at HL-LHC/FCC, vanishing CMB anomalies, a null 3.55 keV line, or a junction computation that misses its anchors. Until then, it holds as a coherent, falsifiable, zero-parameter explanation of everything we see.*
