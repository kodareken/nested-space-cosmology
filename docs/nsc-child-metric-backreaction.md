# The locked child backreaction gate stops at its Cauchy surface

This gate asks whether the completed charged angular plus compact CTP tensor
can evolve the stored Bronnikov child metric from the exact neck jet.  It keeps

$$
q=4,\qquad \Omega=3.973074368754331,\qquad
\zeta=15.785319939652627,\qquad V_{\rm full}=0
$$

and the Wilsonian Einstein coefficient

$$
A=0.04501936182826115
$$

fixed.  The existing ADM/spherical equations, Bronnikov chart and CTP tensor
are imported.  No gravity equation or source is recomputed.

## The constraints are the first evolution gate

For the locked unit-radius minimum, the existing two-derivative source map
requires

$$
\rho_{\rm required}=2A=0.0900387236565223,
\qquad T_{01,{\rm required}}=0.
$$

The completed CTP tensor instead has

$$
\rho_{\rm CTP}=0.100748362886786,
\qquad T_{01,{\rm CTP}}=0.00122387015580509.
$$

Using the imported homogeneous ADM normalization, the two initial residuals
are

$$
\boxed{
\mathcal C_H
=H_\perp^2+2H_\parallel H_\perp+\frac1{r^2}
-\frac{\rho}{2A}
=-0.118944813912720,
}
$$

$$
\boxed{
\mathcal C_M=-\frac{T_{01}}{2A}
=-0.0135927088490713.
}
$$

Both exceed the declared $2\times10^{-10}$ tolerance at $T=0$.  The
Hamiltonian mismatch corresponds to

$$
\Delta\rho=\rho_{\rm required}-\rho_{\rm CTP}
=-0.0107096392302639,
$$

and homogeneous momentum closure requires

$$
\Delta T_{01}=-0.00122387015580509.
$$

The initial null components remain negative,

$$
T_{++}=-0.246591992070731,
\qquad T_{--}=-0.251487472693951.
$$

Their signs therefore pass at the surface, but sign alone does not place the
metric and source on the same constraint surface.

## Why no time steps are reported

Hamiltonian and momentum constraints are conditions on the initial Cauchy
slice.  Evolution cannot repair an initial violation; it only propagates it.
Consequently the trajectory record contains the complete $T=0$ row and stops
before replaying the costly CTP covariance.  This is the requested gate's
mathematically determined failure point, rather than an incomplete numerical
run.

Two diagnostics show why an automatic adjustment would change the problem.
Matching the density alone would require

$$
A=0.0503741814433931,
$$

an $11.8944813913\%$ coefficient change.  Keeping $A$ fixed and changing only
the transverse expansion would require

$$
H_\perp(0)=0.0377389178484633,
$$

so $\dot r(0)\ne0$ and the initial slice would no longer be the stored minimum;
the momentum constraint would still fail.  Neither change is applied.

The next owner is singular and explicit: a constraint-complete same-action
neck assembly must supply
$(\Delta\rho,\Delta T_{01})$ or derive a different initial geometry from the
same source.  Until that dependency changes, a longer background projection
or CTP Hessian cannot start from this branch.

The focused reproducer is

```sh
python3 scripts/derive_nsc_child_metric_backreaction.py --check
```

It reads and authenticates the completed tensor, scale, ADM and spherical
records.  It runs no angular, compact, MMP, metric or historical generator.
