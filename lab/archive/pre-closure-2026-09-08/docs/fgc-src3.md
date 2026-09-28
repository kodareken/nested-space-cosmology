# FGC-1-SRC3: reference-covariant GR-0 source evaluator

`FGC-1-SRC3` closes the numerical-instrument obligation opened by
[SRC2/PREF11](fgc-src2-pref11.md). It derives and freezes a binary64 evaluation
order for the unchanged GR-0 REF1 equations that does not form the physical
and flat-spherical connection derivatives separately before subtracting them.

The machine decisions are:

```text
SRC3_reference_covariant_binary64_evaluator_derived = true
SRC3_exact_reference_bitwise_preserved = true
SRC3_compact_and_independent_exact_oracle_controls_passed = true
SRC3_captured_full_grid_source_gate_passed = true
SRC3_runtime_source_instrument_authorized_for_a_separately_frozen_GR0_protocol = true
amplitude_three_spectral_veto_cleared = false
PROTO12_frozen = false
fresh_GR0_calibration_completed = false
classical_spherical_diagnostic_authorized = false
FGCQR_holdout_execution_authorized = false
retained_EFT_evolution_authorized = false
physical_transition_claim_authorized = false
```

This is an instrument result. It does not turn an existing rejected trajectory
into an accepted one and it contains no new trajectory.

## The diagnosed cancellation

The legacy evaluator constructs

```text
Gamma(g),  partial Gamma(g),  Gamma_bar,  partial Gamma_bar
```

and only then forms their differences. In spherical coordinates, the physical
and reference terms separately contain contributions proportional to `1/r`
and `1/r^2`, even arbitrarily close to Minkowski space. PREF11 showed that at
the first positive-radius point of the captured source call the exact dyadic
six-by-six system is nonsingular, but the binary64 residual loses enough of
the small difference to falsely reject a correctly rounded exact root.

SRC3 changes the arithmetic graph, not the equations. Define

```text
h_ab = g_ab - gbar_ab
```

and use the flat spherical reference derivative to form the tensorial
connection difference directly:

```text
C^a_bc = Gamma^a_bc - Gammabar^a_bc
       = 1/2 g^ad (̅nabla_b h_dc + ̅nabla_c h_db - ̅nabla_d h_bc).
```

The implementation differentiates this identity before rounding to obtain
`partial_e C^a_bc`. Since the reference is flat, Ricci is then assembled as

```text
R_bd = partial_a C^a_db - partial_d C^a_ab
     + Gammabar^a_ae C^e_db + C^a_ae Gammabar^e_db + C^a_ae C^e_db
     - Gammabar^a_de C^e_ab - C^a_de Gammabar^e_ab - C^a_de C^e_ab.
```

The modified-harmonic constraint and its derivative are likewise assembled
from `C` and `partial C`. This avoids subtracting two already-rounded
`1/r^2` connection derivatives. The action, metric equations, scalar
equations, flat spherical reference, auxiliary cones, independent fields,
affine acceleration branch, and six REF1 rows are unchanged.

## Frozen proof burden

The source implementation is
[`src3_reference_balanced_source.py`](../src/recursive_horizons/fgc/evolution/src3_reference_balanced_source.py).
Its frozen controls are literal binary64 values in
[`fgc-1-src3-controls.json`](../configs/fgc/fgc-1-src3-controls.json), and its
scope and thresholds are fixed by
[`fgc-1-src3.toml`](../configs/fgc/fgc-1-src3.toml).

The proof has four layers.

### 1. Exact reference

At seven frozen dyadic radii from `1/64` through `16`, exact spherical
Minkowski data produces bitwise-zero physical metric rows, gauge rows, scalar
rows, Hamiltonian and momentum projections, Ricci invariants, acceleration,
and verified residual. This is equality, not a near-reference tolerance or an
equality-only runtime bypass.

### 2. PREF11 compact point

At the captured `r=1/64` lower jet, the SRC3 affine root has complete
binary64 residual

```text
1.7165613297609217e-27
```

and exact-dyadic residual

```text
7.642196077734746e-27.
```

Its four nonzero components lie within four binary64 ULPs of the independently
solved, correctly rounded exact root. The legacy complete evaluator still
reports a value above `1e-12` at the SRC3 root, so the control retains the
original failure rather than weakening its threshold.

### 3. Independent nontrivial exact controls

Two dyadic, non-Minkowski lower jets were selected outside the CAP1 capture.
For each one, the exact `Fraction` oracle evaluates the zero acceleration and
all six unit-acceleration seeds. The SRC3 binary64 residual agrees with those
exact values after scale normalization within `16 epsilon_64`; its affine root
also passes the unchanged strict `1e-12` residual gate. These controls prevent
the captured point from becoming a one-fixture special case.

### 4. Captured full grid

When the optional hash-bound CAP1 raw bundle is present, the certificate
recomputes all 8,192 positive-radius systems. The unchanged legacy evaluator
returns

```text
1.0659145473163184e-12 > 1e-12,
```

while SRC3 returns

```text
3.296668493746324e-14 < 1e-12.
```

The largest SRC3 residual occurs in row `3` at point `872`, `r=13.640625`,
not at the centre-adjacent cancellation wall. The maximum infinity-norm
kinetic condition estimate is unchanged at `2859227.196975944`, far below the
frozen `1e10` limit. Canonical hashes bind the full acceleration and residual
arrays. A clean clone can reproduce the exact compact and independent controls
without the ignored raw bundle; if any part of that bundle is present, all
four CAP1 files and their frozen hashes are mandatory.

## Prospective boundary

The existing CAP1 fixture is disclosed development evidence. No fresh GR-0
trajectory was inspected while defining SRC3. The evaluator may be consumed
only by a separately frozen successor protocol and fresh namespace after the
canonical SRC3 certificate passes.

SRC3 does not adjudicate amplitude `3`'s independent direct coarse-to-medium
`phi/Lambda` derivative-tail veto. The next gate is therefore
**FGC-1-RSP1**, a prospective resolution-spectrum study. It may test whether
that veto reflects under-resolution, a discretization-specific diagnostic, or
a genuine lack of convergence, but it may not alter the source evaluator or
tune the existing threshold after seeing a successor outcome. Only after both
SRC3 and that independent gate pass may PROTO12 be considered.

## Scientific boundary

SRC3 establishes that one precisely identified GR-0 source wall was caused by
binary64 evaluation order and supplies a bounded replacement instrument. It
does not establish an eligible GR-0 collapse case, trapped-region evolution,
SGB-L or FGC-QR dynamics, regulator activation, affine-null defocusing,
singularity resolution, a daughter domain, a dark-sector mechanism, variable
locally measured light speed, retained-EFT validity, or a physical model of
nature. A later model failure cannot be declared an arithmetic artifact merely
because this earlier arithmetic artifact was real; each later obstruction must
be localized independently.

## Reproduction

```bash
python3 scripts/reproduce_fgc_src3.py --check
```

The command performs no evolution. With the optional CAP1 bundle it also
recomputes and hashes the captured full grid; without it, the command verifies
the immutable bundle binding and reproduces the exact compact and independent
controls.
