# Coupled evolution stops at the first same-state source evaluation

The source-selected neck passes the ADM constraints when the previously
integrated tensor is held fixed.  Joint evolution asks a stronger question:
does the same Gaussian covariance produce that tensor on the replacement
geometry, and what pressures does it produce one step later?

## The geometry change acts on the source immediately

The replacement radius is

$$
r_\star=0.945328394434129
$$

instead of the source calculation's unit radius.  Consequently the angular
term in every retained Dirac Hamiltonian changes by

$$
\frac{1/r_\star}{1/r_0}-1
=0.0578334533139646.
$$

The existing metric-force identity gives the fixed-covariance radius vertex

$$
\boxed{
\left.\frac{\partial\rho}{\partial\log r}\right|_C
=-2(\rho_L+p_{\perp,L})
=-0.460687999769215.
}
$$

It is nonzero.  The old tensor therefore cannot be copied unchanged onto the
new radius and called a fresh CTP evaluation.  Its linearized finite-shift
diagnostic is $+0.02590054110$; that number is not substituted as the new
density.

## Why the requested covariance cannot yet be advanced

The charged angular and positive-compact result records retain mode labels,
quadrature metadata, integrated stress rows and covariance eigenvalue bounds.
They do not persist the complex Gaussian covariance for every
compact/angular/frequency mode.  Four integrated stress moments cannot recover
those matrices: future evolution depends on

$$
\dot C_j=-i[H_j[g],C_j],
$$

and many covariances with the same current vertex traces have different
commutators and later pressures.

The Landau transformation recorded for the stress tensor also does not by
itself define the spinor two-point function on the new Cauchy slice.  The
canonical half-density gives the correct common Hilbert-space representation,
but its mode covariance, slice map and fourth-order reference must be carried
explicitly.  Replaying the old source generator alone would reconstruct the
old unit-radius slice, not this missing map.

## Recorded break and next executable contract

The trajectory contains the accepted fixed-tensor $T=0$ row and stops before
the first metric step.  Its fixed-tensor constraints are at
$2.2\times10^{-16}$ and both null signs are negative.  The break occurs at

$$
\boxed{T_{\rm break}=0}
$$

during same-covariance stress re-evaluation, before a constraint or null-sign
evolution test exists.  Advancing with constant $w$, frozen pressures or a
four-moment surrogate would introduce a new source law and is not used.

The exact next owner is `ModeResolvedCauchyState`:

1. persist frequency weights and complex covariance matrices for every
   retained angular and compact channel;
2. transport them in the canonical half-density basis to the source-selected
   Landau normal;
3. carry the existing fourth-order reference and local induced allocation on
   that same metric history;
4. pass this state to `CausalCommonFunctional` and the ADM stepper.

This is the shallowest path that produces later pressures from the state
rather than inferring an equation of state from the four neck moments.

The focused reproducer is

```sh
python3 scripts/derive_nsc_coupled_ctp_metric_evolution.py --check
```

It authenticates the source, constraint and CTP-owner records.  It runs no
angular, compact, metric, MMP or historical generator.
