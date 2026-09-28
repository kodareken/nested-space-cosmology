# Bounded normalization-scale flow of the covariant Dirac remainder

The already calculated static Dirac remainder depends on the subtraction
mass `M` only through the logarithmic heat term. This note derives that
dependence, the three metric Euler–Lagrange densities that go with it, and
the finite local coefficient shifts that cancel a change `M → s M` for the
massless torsionless operator. Unspecified initial coefficients stay
`None` or symbols; they are not replaced by a selected zero action. The
scale ratio `s` is not a recursive inheritance factor. No independently
weighted gravitational or scalar term is added, and no throat is fitted.

The [record](../results/nsc-12-normalization-flow.json) is generated and
checked by `scripts/check_nsc_normalization_flow.py`. The owner is
`src/recursive_horizons/nsc_normalization_flow.py`. It reads the existing
covariant subtraction and finite-term libraries; it does not replace them.

## Claim

### Local Interpretation

Inside one room, changing the subtraction mass at fixed cutoff and fixed
geometry changes only the finite local part of the already defined
remainder. Physical predictions of a completed common action cannot depend
on that label.

### Non-Local Interpretation

A parent/child inheritance ratio `Omega` is a different object: it maps
solved rooms through the recursive operator. The present `s = M'/M` only
relabels the free-Dirac subtraction. The two ratios are not identified.
This logarithm also does not supply a physical Newton or cosmological
beta function.

### Simplest Relational Analogue

A logarithm changes by an additive constant when its argument is
rescaled. The heat density `a4` is that constant for this operator. On a
smooth closed cell the continuum bulk piece of `a4` is `C^2`.

### Constructive Mathematics

The identities below are exact. Numerical controls compare
determinant-source changes to the derived coefficient increments. Finite-code
cancellation keeps the discrete Euler column of the existing owner; grid
refinement shows that column is product-rule error.

### Inherited Invariant

The combination `E_sub(M) + c(M) · F(M_ref)` is invariant under `M → s M`
after the derived shift of `c`. Continuum bulk running on this cell is
`C^2`. The `R^2` coefficient does not run.

### Observable or Logical Consequence

Unmatched `M` changes the remainder energy and every `N,q,r` gradient.
Matched discrete shifts cancel those changes. A `C^2`-only shift already
cancels the continuum bulk; its coarse-grid leftover vanishes at 64 and
128 points.

### Next Intervention

The concurrent Weyl compensator, the remaining finite matching data, and
the recursive `Omega` map are separate owners. This flow does not supply
them, and it does not infer `G` or `Λ` running.

## Remainder and exact `M`-derivative

The covariant source defines

\[
E_{\rm sub}(M)
=
E_\Lambda
-
\frac{
A_0\Lambda^4/2+A_2\Lambda^2+A_4\log(\Lambda^2/M^2)
}{32\pi^2}.
\]

`E_Lambda` is independent of `M`. Differentiating at fixed geometry and
cutoff gives

\[
\frac{d E_{\rm sub}}{d\log M}
=
\frac{A_4}{16\pi^2}.
\]

The local heat density of the massless torsionless Dirac operator, in the
Lorentzian finite-term basis with `R_E=-R_L` and `Box_E R_E=Box_L R_L`, is

\[
a_4=\frac{-18 C^2+11 E_4-12\Box\mathcal R}{360}.
\]

`R^2` is absent. The integrated coefficient is the same linear combination
of the six finite-term energies `F_i`. The Euler–Lagrange densities of
`A_4` are that combination of the existing `e_{if}`, `f=N,q,r`, so

\[
\frac{d}{d\log M}\frac{\delta E_{\rm sub}}{\delta f}
=
\frac{1}{16\pi^2}\frac{\delta A_4}{\delta f}.
\]

Nodal gradients use the covariant-source convention: the inner product
against `δf` is an energy variation.

## Two masses

The finite local basis is

\[
I_i=(M_{\rm ref}^4,\,M_{\rm ref}^2\mathcal R,\,C^2,\,\mathcal R^2,\,E_4,\,\Box\mathcal R).
\]

`M_ref` is held fixed while the determinant normalization `M` varies.
`M_ref^4` and `M_ref^2 R` therefore do not run in this basis: the
subtraction keeps cutoff powers `Lambda^4` and `Lambda^2`, not `M^4` or
`M^2`. That scheme fact is not a physical Newton or cosmological-constant
beta function. Changing `M_ref` rescales those two functionals and is a
different operation from changing `M`.

If the same basis is rewritten with `M` itself,
`\hat c_{M^4}=c_{M^4}(M_{\rm ref}/M)^4` and
`\hat c_{M^2R}=c_{M^2R}(M_{\rm ref}/M)^2`. The equivalent beta functions
are then the classical

\[
\beta_{M^4}=-4\hat c_{M^4},\qquad
\beta_{M^2R}=-2\hat c_{M^2R},
\]

together with the same `a4` anomaly as in the `M_ref` basis. Those `-4`
and `-2` factors are unit rewriting, not running of `G` or `Λ`.

## Finite coefficient shifts

Write `F=\sum_i c_i F_i` for arbitrary initial `c_i`. Unspecified `c_i`
remain `None` with symbolic names `(c_{M^4},c_{M^2R},c_{C^2},c_{R^2},c_{E_4},c_{\Box R})`.
They are not replaced by the zero vector. Under `M→sM` at fixed `M_ref`,

\[
\Delta E_{\rm sub}
=
\frac{A_4\log s}{16\pi^2}.
\]

The compensating increments are `Δc_i=β_i\log s` with

\[
\begin{aligned}
\beta_{M^4}&=0,\\
\beta_{M^2R}&=0,\\
\beta_{C^2}&=\frac{1}{320\pi^2},\\
\beta_{R^2}&=0,\\
\beta_{E_4}&=-\frac{11}{5760\pi^2},\\
\beta_{\Box R}&=\frac{1}{480\pi^2}.
\end{aligned}
\]

These are constant. The flow is linear in `\log M` along the single
direction `a4`. It has no finite fixed point along that direction, and it
is not a recursive `Omega` fixed point. The `R^2` coefficient is marginal
and stationary.

## Continuum bulk running versus discrete Euler leftover

On a smooth closed fixed-topology four-dimensional cell the Euler density
is a total derivative and the Euler functional derivative vanishes
identically. The box term is likewise a total derivative, with vanishing
bulk Euler–Lagrange density in the finite-term owner. Therefore the
continuum bulk running of this logarithm is `C^2` alone:

\[
\Delta c^{\rm bulk}(s)=\bigl(0,0,\log s/(320\pi^2),0,0,0\bigr).
\]

The finite-term owner still differentiates its discrete Euler column as a
control. Finite-code cancellation therefore keeps `Δc_{E_4}` so that the
discrete subtraction owner cancels field-by-field. The difference between
that full discrete shift and the continuum `C^2`-only shift is
product-rule error, not topology-generated bulk stress. On the nonconstant
`R=2` cell at `s=2` the `C^2`-only gradient leftover is about `1.6\times10^{-10}`
at 32 points and falls below `10^{-12}` at 64 and 128 points. Local
fourth-derivative roundoff of the discrete Euler column may increase
slightly from 64 to 128; that is the known finite-term owner artefact and
does not restore continuum Euler stress.

Composition and inversion of the full increment are exact:

\[
\Delta c(s_1)+\Delta c(s_2)=\Delta c(s_1 s_2),\qquad
\Delta c(s)+\Delta c(1/s)=0.
\]

## Numerical controls

The determinant-source change is the change of the existing local
subtraction; `E_Lambda` cancels. That change is compared to the derived
`Δc·F` and `Δc·e_f`, not used to select `c_i`. Controls:

- cylinder `L=4`, `a=1`; smooth `R=2`; nonconstant `N,q,r` on the same
  `R=2` profile;
- scale ratios `1/2`, `2`, `e`, `10`, `1/π`;
- unspecified coefficients (`None`) and one explicit unfitted control
  vector;
- centered finite-difference `d/d\log M` of the subtraction, independent
  of the beta formula;
- composition `2,3` versus `6` and the inverse pair `2,1/2`;
- unmatched `M` without `Δc`, which changes energy and every metric
  gradient;
- 32/64/128 refinement of the discrete Euler leftover versus `C^2`-only
  bulk running;
- a tiny existing remainder evaluation showing `E_Lambda` independent of
  `M`.

On the cylinder the analytic value is `d E_sub/d\log M=-1/(15π)`.

## Comparison tolerance

Primitive identities, matched cancellations, and ordinary numerical fields
use `atol=2e-12` and `rtol=2e-12`. Centered finite-difference residuals of
the local subtraction,

`energy_minus_analytic` and `maximum_gradient_minus_analytic`,

are cancellation-sensitive (observed subtraction error about `10^{-11}` at
step `10^{-5}`) and are compared at `atol=2e-10`. Keys, booleans, strings
and hashes remain exact. The comparator is tested on both classes.

## Scope

The Weyl compensator is a concurrent derivation in the main laboratory and
is not this `M`-flow. Four independent finite bulk matching conditions,
state, link, and a self-sourced geometry remain open. The calculation does
not construct a complete ultraviolet functional, does not infer physical
Newton or cosmological running, and does not close a recursive fixed point.

Run:

```sh
python3 scripts/check_nsc_normalization_flow.py --check
PYTHONPATH=src:. python3 -m unittest discover -s tests -p test_nsc_normalization_flow.py -v
```

Existing outputs cannot be overwritten by `--output`. `--output` and
`--check` are mutually exclusive.
