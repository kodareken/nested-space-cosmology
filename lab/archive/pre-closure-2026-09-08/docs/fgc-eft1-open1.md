# FGC-1-EFT1-OPEN1: retained-EFT open-run authorization audit

**FGC-1-EFT1-OPEN1** defines the conjunction that must pass before the
repository may call a nonlinear calculation an FGC-QR retained-effective-field-
theory evolution. It then evaluates that conjunction against the exact current
artifacts. The result is deliberately negative: the present evidence does not
authorize a run.

This is a useful result rather than a placeholder. It prevents the exact local
component inequalities in EFT0-LED1, the compact radial hyperbolicity theorem
in DOM3-UHYP1, the frozen annular main-system boundary estimate in BND1-MD1,
or the conditional boundary-free gauge theorem in CON3-CAU1 from being
silently promoted into a claim none of them establishes.

## The required hierarchy

For a retained EFT on a spacetime run domain \(D\), a declared cutoff
\(\Lambda\) is only the beginning. The calculation needs invariant bounds on
every physical scale that enters the derivative expansion. Schematically,

\[
 \epsilon_D=\sup_D\max\left(
 \frac{|\phi|}{\Lambda},
 \frac{|\nabla\phi|^{1/2}}{\Lambda},
 \frac{|\nabla\chi|^{1/2}}{\Lambda},
 \frac{|\mathrm{Riemann}|^{1/2}}{\Lambda},
 \frac{\omega_{\rm proper}}{\Lambda},\ldots
 \right)<\epsilon_\star<1,
\]

with each schematic norm replaced by a dimensionally correct invariant in the
implemented convention. A finite coordinate-component box at one radius does
not bound Fourier support or proper frequency: arbitrarily rapid oscillations
can have arbitrarily small amplitude. Nor does a short representative list of
omitted operators bound the complete truncation error.

The audit therefore requires all of the following, simultaneously:

1. a declared cutoff/matching scale not inferred from the numerical fixture;
2. the operator basis allowed by the frozen symmetries through a declared
   target order, with matched or bounded coefficients;
3. covariant field-amplitude, gradient, curvature, and proper-frequency bounds
   on a nonzero spacetime run domain;
4. a uniform bound on the **sum** of omitted contributions below a declared
   tolerance;
5. multidirectional full-system hyperbolicity and a physical/gauge
   constraint-complete local-existence theorem or IBVP;
6. regular-centre or declared-boundary data; and
7. a bootstrap monitor that stops the evolution before any premise fails.

The authorization rule is an all-of conjunction. One false predicate stops
the run; passing radial or amplitude subtests cannot compensate for a missing
frequency, remainder, constraint, or existence bound.

## What passes now

- EFT0-LED1 declares \(\Lambda=16\) as an external assumption, explicitly not
  a value inferred from the FGC-QR fixture, and its four exact local
  component-amplitude controls pass.
- DOM3-UHYP1 proves compact nonflat **radial** strong hyperbolicity on the
  implicit REF1 branch graph.
- BND1-MD1 proves uniform frozen radial main-system maximal dissipation and
  compatibility with the zero-speed kinematic reduction subsystem.
- CON3-CAU1 proves the conditional boundary-free uniqueness implication that
  preserves the metric-derived gauge vector on any sufficiently smooth full
  solution with compatible Cauchy data.

These are retained as positive evidence. They are not erased by the failed
authorization.

## What blocks the run

EFT0 explicitly has only a representative omitted basis, no proper-frequency
bound, no covariant full-jet/curvature run box, and no omitted-remainder sum.
DOM3 is radial rather than multidirectional. BND1 does not supply the incoming
metric-derived gauge or physical-constraint map, a regular centre, or a
quasilinear IBVP. CON3 supplies conditional gauge uniqueness, not a nonzero-
width compatible hypersurface or solution existence. No evolution exists on
which a bootstrap invariant could yet be proved.

Consequently the machine record sets

```text
retained_EFT_open_run_envelope_passed = false
evolution_authorized = false
```

and carries a stop mask into COL1 and DEF1. The appropriate next scientific
work is to construct the missing evidence, not to tune collapse data until a
desired bounce appears.

Reproduce with:

```bash
python3 scripts/reproduce_fgc_eft1_open1.py
```
