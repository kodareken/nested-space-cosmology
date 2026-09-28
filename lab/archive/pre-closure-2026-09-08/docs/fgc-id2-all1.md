# FGC-1-ID2-ALL1: all-amplitude, all-case static initial admission

**Authors:** Douglas Ek & ChatGPT 5.6 Sol
**Status:** machine-reproduced pre-calibration input ledger; no trajectory

```text
all_amplitude_all_case_static_initial_admission_completed = true
fresh_GR0_dynamic_calibration_authorized = false
classical_spherical_diagnostic_authorized = false
retained_EFT_evolution_authorized = false
physical_transition_claim_authorized = false
```

## Question

PROTO4 inherits seven ordered GR-0 calibration amplitudes. The amplitude that
first produces a controlled trapped interval is later copied into five frozen
FGC-QR holdout perturbations. PROTO3 also requires the constraint, centre,
finite-mass, compactness, and no-initial-trapping premises of those expanded
inputs to be checked before GR-0 dynamic eligibility.

The selected amplitude is not known before calibration. Therefore checking
only five slices after selection would be outcome-dependent. The smallest
outcome-neutral static ledger is

\[
  7\ \text{amplitudes}\times 5\ \text{FGC-QR modifiers}=35
\]

FGC-QR slices, together with all seven GR-0 calibration inputs. This artifact
constructs that Cartesian product before reading a trajectory.

## Frozen inputs

The amplitude order remains

\[
  2,\;\frac52,\;3,\;\frac72,\;4,\;\frac92,\;5.
\]

For every amplitude, the five FGC-QR modifiers are:

| Case | Matter-width factor | Regulator-seed factor | Exact matter centre buffer |
|---|---:|---:|---:|
| `FGCQR-CENTRAL` | \(1\) | \(1\) | \(5L_0/2\) |
| `FGCQR-SEED-HALF` | \(1\) | \(1/2\) | \(5L_0/2\) |
| `FGCQR-SEED-DOUBLE` | \(1\) | \(2\) | \(5L_0/2\) |
| `FGCQR-WIDTH-SEVEN-EIGHTHS` | \(7/8\) | \(1\) | \(41L_0/16\) |
| `FGCQR-WIDTH-NINE-EIGHTHS` | \(9/8\) | \(1\) | \(39L_0/16\) |

The regulator has its separately frozen half-width \(L_0/2\). Thus the
combined exact-vacuum centre buffer in the narrowed matter case is still
\(5L_0/2\); the table's \(41L_0/16\) is the matter-profile buffer required by
the future five-row holdout manifest.

Every input is constructed independently with RK4 and SSPRK3 radial constraint
integration on the full nested grids

\[
  N_r\in\{1025,2049,4097\},\qquad 0\le r\le128,
\]

using the unredefined ACT1/VAR1 physical constraints and the frozen
gauge-compatible initial time derivatives.

## Complete-grid construction

For FGC-QR the support solve is embedded between exact continuations:

\[
  k=0,\quad\lambda=1
  \qquad\text{inside the centre vacuum},
\]

and

\[
  k=\frac{J}{r^3},\qquad
  \lambda^{-2}=1-\frac{2M}{r}+\frac{J^2}{r^4}
  \qquad\text{outside the compact support}.
\]

The generated grid carries the canonical state

\[
 U=(u,p,q),\qquad
 u=(\alpha,v,\lambda,R,\phi,\chi),
\]

with exact centre parity, \(R(0)=0\), and
\(\partial_rR(0)=\lambda(0)=1\). The certificate evaluates the specialized
Hamiltonian and radial-momentum residuals at every positive-radius grid point,
not merely at the pulse peak. Each complete \((u,p,q)\) array is hashed as
little-endian binary64 with its branch, shape, and field-block labels. There
are 252 distinct state hashes:

\[
  (7\ \mathrm{GR0}+35\ \mathrm{FGCQR})
  \times2\ \mathrm{methods}\times3\ \mathrm{grids}.
\]

## Static decision rule

An FGC-QR slice passes this PROTO3 static contract only when both methods and
all three resolutions establish:

- complete-grid physical-constraint residual below \(10^{-10}\);
- regular centre plus exact inner and outer scalar-vacuum buffers;
- finite positive mass and positive analytic exterior metric factor;
- a finite positive constraint Jacobian and effective Planck ratio at least
  \(1/2\);
- no initially trapped sphere;
- peak initial compactness in the inclusive frozen interval
  \([1/10,3/4]\);
- convergent constraint profiles with observed RK4 order at least \(3\) and
  SSPRK3 order at least \(5/2\);
- fine-grid cross-method mass and compactness agreement within the frozen
  tolerances; and
- every individual PROTO4 weighted spectral budget below its absolute field,
  derivative, and RMS-scale threshold.

The GR-0 input additionally must pass HLT2's two-method, three-grid nested
weighted spectral rule. An amplitude becomes eligible for fresh dynamic GR-0
calibration only when that GR-0 input and all five corresponding FGC-QR static
slices pass.

## Result

The ledger leaves exactly two amplitudes eligible, in the original order:

\[
  \boxed{\frac52,\;3}.
\]

Two amplitude/case combinations fail the frozen compactness window:

- amplitude \(2\), width \(9/8\): the converged peak is about \(0.0954\),
  below \(0.1\);
- amplitude \(5\), width \(7/8\): the converged peak is about \(0.794\),
  above \(0.75\).

Two further narrowed FGC-QR inputs miss an individual spectral budget:

- at amplitude \(4\), the coarsest primary-grid
  `lambda_minus_1` derivative-tail fraction is approximately
  \(1.0356\times10^{-3}\), slightly above \(1/1024\);
- the amplitude \(9/2\), width \(7/8\) input also misses an individual
  weighted budget.

The GR-0 amplitudes \(7/2\), \(9/2\), and \(5\) independently fail the
literal PROTO4 nested-tail admission on their static calibration inputs. No
threshold or ordering is changed in response. Those amplitudes are skipped
before dynamics; they are not failed collapse experiments.

In total, 31 of 35 FGC-QR slices pass the declared static premises. This is a
real narrowing of the executable experiment, but it says nothing about whether
the regulator activates or defocuses.

## Disclosed nested-tail diagnostic

Apart from the two narrowed high-amplitude inputs named above, the individual
FGC-QR weighted field, derivative, and RMS-scale budgets pass. Only 17 of the
35 cases pass the additional literal nested-tail ratio on both methods. In the
failures, an already extremely small high-frequency tail can
decrease sharply from coarse to medium and then decrease too slowly—or rise
slightly—between medium and fine after interpolation and roundoff dominate.
For example, a tail can remain many orders below its absolute admission ceiling
while its ratio exceeds \(1/4\).

PROTO3's future static manifest has exactly five booleans: constraint
compatibility, regular centre, finite mass, compactness window, and no initial
trapping. It does not name the nested FGC-QR spectral ratio. ID2 therefore
serializes the 17/35 result but does not retroactively add it to that static
decision. HLT2 remains unchanged. Before a later FGC-QR trajectory can be
accepted, the runtime monitor must still resolve when and how its nested
spectral comparison applies, without using an outcome to move a threshold.

## Scientific boundary

`FGC-1-ID2-ALL1` proves that the repaired protocol has a nonempty, fully
enumerated static route into fresh GR-0 calibration. It does **not** run that
calibration, form a trapped sphere dynamically, activate \(\phi\), calculate a
Raychaudhuri margin, resolve `PRO4-HLD1`, supply the missing SGB-L branch-owned
health definition, or convert constraint residuals into Raychaudhuri-observable
error units.

The two excluded inputs are failures of the frozen input window, not failures
of FGC-QR, Finite Gradient Closure, a black-hole transition, or a general
gradient mechanism. Retained-EFT evolution and every physical-transition claim
remain false.

## Reproduction

```bash
python3 scripts/reproduce_fgc_id2_all1.py \
  --config configs/fgc/fgc-1-id2-all1.toml \
  --output results/fgc-1-id2-all1.json
```

The authoritative machine record is
[`results/fgc-1-id2-all1.json`](../results/fgc-1-id2-all1.json). The source
configuration, predecessor records, complete-grid implementation, reproducer,
and this derivation are SHA-256 bound inside it.
