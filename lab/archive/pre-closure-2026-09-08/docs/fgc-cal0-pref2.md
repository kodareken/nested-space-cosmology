# FGC-1-CAL0-PREF2: PROTO3 pre-holdout calibration-contract diagnosis

## Decision

`FGC-1-CAL0-PREF2` rejects **FGC-2-SF1-PROTO3 as an executable resolved
holdout contract**, before any FGC-QR or SGB-L evolution outcome is opened.
The obstruction is in the test premises, not in the FGC-QR field equations.

Two independently checkable defects cause the stop.

1. PROTO3 defines the nonzero initial matter pulse with exact compact support,
   but HLT1 declares the largest discrete Fourier bin above a relative
   `2^-40` amplitude floor to be the physical support and stops whenever that
   bin lies in the top eighth of the grid. The declared pulse itself occupies
   that top band at `t=0` on all three frozen grids. Amplitude scaling cannot
   repair the failure because it leaves every relative Fourier ratio
   unchanged.
2. PROTO3 says that every RUN1 health and constraint stop applies to its GR-0
   amplitude calibration, while HLT1, DOM4, and HYP2 explicitly own an
   FGC-QR target-branch envelope. The FGC deformation monitors vanish on
   GR-0, its acceleration solve is a unique affine Einstein root, and a raw
   coordinate-time acceleration cap is not an invariant criterion for whether
   a GR collapse baseline forms a trapped sphere. No branch-applicability map
   says which target-only predicates are meaningful during calibration.

The machine decision is:

```text
PROTO3_pre_holdout_numerical_contract_obstruction_verified = true
PROTO3_resolved_holdout_manifest_authorized = false
PROTO4_premise_revision_required = true
```

This is not a collapse result, a failed defocusing result, or a rejection of
Finite Gradient Closure.

## The compact-support versus last-bin problem

PROTO3 freezes

```text
r chi(r) = A B((r-12)/2),
B(x) = exp(1-1/(1-x^2)) for |x|<1 and zero otherwise.
```

The profile is smooth and compactly supported, but compact support does not
mean exact band limitation. Its Fourier transform has a decaying tail. HLT1's
current estimator identifies support with the final bin whose amplitude is
larger than

```text
max(1e-30, 2^-40 times the spectrum peak).
```

That rule is useful as an injected alias detector, but it is too literal as a
physical-admission rule for this input. The canonical reproducer samples the
actual `chi` field over the frozen measurement interval, applies HLT1's own
compact window and `windowed_spectral_support`, and obtains top-eighth
occupation at every declared resolution. The nonzero-amplitude rescaling
control gives the same support bins for `A=2` and `A=5`.

The reproducer also reports, without using it to excuse PROTO3, the fraction
of total spectral energy and derivative-weighted spectral energy in the top
eighth. Those quantities converge rapidly even while the last-bin Boolean
remains true. This shows exactly why a successor should bound unresolved
spectral weight and require cross-resolution convergence instead of demanding
that every Fourier coefficient above a near-roundoff floor disappear.

For the first declared amplitude (`A=2`), the canonical values are:

| Full-domain points | HLT1 support/Nyquist bin | Top-eighth field power | Top-eighth derivative-weighted power |
|---:|---:|---:|---:|
| 1025 | `96/96` | `8.054e-8` | `6.347e-5` |
| 2049 | `192/192` | `1.587e-10` | `4.793e-7` |
| 4097 | `384/384` | `1.387e-13` | `1.741e-9` |

The `A=5` control gives the same bins and the same fractions to floating-point
roundoff. Thus refinement makes the top-band weight smaller by orders of
magnitude while the Boolean rejection becomes no less absolute.

NUM1 remains correct within its scope: it verifies that the estimator detects
a low sinusoid and an injected near-Nyquist sinusoid. It never claims that the
declared compact pulse passes the estimator. CAL0 closes that missing
composition.

## Calibration is not target-branch containment

PROTO3 uses GR-0 to choose the first matter amplitude that demonstrably forms
a trapped sphere before the final time. That calibration must be numerically
converged, constraint controlled, centre regular, Lorentzian, causally
isolated from the boundary, and independently reproduced. It must not silently
inherit conditions whose only mathematical owner is the FGC-QR local
weak-coupling box.

In particular:

- the all-covector coefficient-deformation inequalities in HYP2 describe the
  FGC-QR principal operator relative to Einstein-scalar theory;
- DOM4's acceleration box is a sufficient local container around sampled
  FGC-QR initial-slice roots, not an invariant theorem about all GR collapse;
- the GR-0 REF1 acceleration system is affine with a unique root whenever its
  kinetic block is nonsingular; and
- coordinate accelerations depend on the chosen time and gauge variables,
  whereas a trapped sphere is defined by the signs of physical metric-null
  expansions.

A successor protocol therefore needs an explicit applicability ledger. GR-0
calibration must retain the universal numerical and geometric stops. FGC-QR
holdout execution must retain its branch, multidirectional-health, activation,
scale, and observable stops. Separating those scopes does not waive a target
premise; it prevents a target-only local coordinate box from defining the
control experiment.

## Required successor shape

The diagnosis permits only a premise-based successor frozen before any
FGC-QR output exists. It must preserve the action, scalar profiles, amplitude
order, held-out cases, physical equations, null observable, and nonclaims. It
must add:

- a branch-applicability table for every stop;
- a spectral-tail energy/error budget with nested-resolution convergence;
- a deterministic trajectory constraint scale and common-event convergence
  rule; and
- a new output namespace and resolved-manifest identity.

The existing PROTO3 file and its `7/8` RUN1 result remain immutable historical
evidence. A successor does not retroactively turn PROTO3 into a passing test.

## Reproduction and nonclaims

Run:

```bash
python3 scripts/reproduce_fgc_cal0_pref2.py \
  --output results/fgc-1-cal0-pref2.json
python3 -m unittest tests/test_fgc_cal0_pref2_reproduction.py -v
```

The artifact binds the active PROTO3, HLT1, DOM4, NUM1, the estimator and
profile implementations, and the GR-0 direct-source cross-check by SHA-256.
It reads no path under the FGC-QR or SGB-L output namespaces. It derives no
collapse, activation, defocusing, transition, singularity resolution, child
domain, dark sector, or varying locally measured speed of light.
