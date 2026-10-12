# Regular normal-data branch of the compatible constraint chart

Reuse the general closed `cu(Ha,Hr)` expression, its unchanged certified
baseline `A0`, and the constant d/e certificate in `c967361`. The only new
coefficient is `A1=partial_w cu(Ha,Hr+w/r)`. No Cu, Cv, principal, reference
quadrature or physical-source producer is rerun.

The missing connection is the domain of the lapse-eliminated spatial ODE.
Differentiate the existing formula, enclose its two affine root locations,
and classify the connected regular w-component containing0. Stop after
this coefficient/branch certificate. These roots are **normal-data values**
`w=delta(r_T)`, not spatial lengths, durations, profiles or physical initial
conditions. No periodicity, boundary data or global solution is chosen.

## Affine extension from the existing coefficient

Let `M2=m^2+L^2`, `L=ell/r` and
`K=-32*pi*C_W/3+(h_q+log r)/(30*pi)`. Direct differentiation gives

\[
A(w)=A_0+A_1w,
\quad A_{1,\mathrm{local}}=a\left(K-\frac1{60\pi}\right),
\]
\[
A_{1,\mathrm{ref},g}=
-\frac{D aL^2(L^2+5m^2)}{120\pi r^2M^6}.
\]

To reuse existing directed data, evaluate the identical expressions
`A1_local=d_local/r-a/(30*pi)` and
`A1_ref,g=(d_ref,g/r)*(1+4*m^2/M2)`. Zero-angular reference slopes vanish
exactly by the formula, without thresholding numerical residuals.

The already-proved general principal identities extend as

\[
B(w)=-\frac{2A(w)+(H_r+w/r)d}{a^2},\qquad e=-d/a^2,
\]
\[
\Delta(w)=A(w)e-B(w)d
=\frac d{a^2}\{R_0+R_1w\},
\quad R_0=A_0+H_rd,\quad R_1=A_1+d/r.
\]

The independent quadratic coefficient has the exact consistency identity
`A1+d/r=-a^2*C`. This is checked symbolically using the frozen reference C
formula and the independently certified local C expression `-2K/a`. The
[separate C certificate](nsc-incoming-surface-coefficients.md) is then read
for an interval cross-check, without using C to compute these roots or
evaluating its D/F producer. This also gives the exact identity
`partial_w Delta=-d*C`; the momentum first-integral reduction is owned
separately.

## Two distinct algebraic obstructions

The lapse pivot root is `w_A=-A0/A1`; the principal root is `w_D=-R0/R1`.
The [record](../results/development/nsc-incoming-surface-regular-branch.json)
stores directed root intervals, signs and their strict ordering. The exact
connected regular branch is expressed using its root symbol; the record
also gives a conservative numerical subinterval valid for every value in
the input enclosures. It never uses the upper endpoint of an uncertain root
interval as though the whole preceding interval were certified regular.

At A=0 the lapse equation cannot be solved for U by that pivot. If Delta
remains nonzero, the original principal matrix is nevertheless nonsingular
and another algebraic pivot is available. At Delta=0 its principal rank
really drops in this chosen family; a lower-order compatibility analysis is
needed, not a declaration of non-existence.

Neither root changes the fixed positive intrinsic N,a,r or the normal Dirac
gaps. Finite normal-data values are not a spacetime singularity. No claim
about a solution reaching or crossing either value follows from this chart
certificate. The fixed rho0 seam, C0, action terms and physical scales stay
unchanged.

```sh
python3 scripts/derive_nsc_incoming_surface_regular_branch.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_surface_regular_branch.py
```
