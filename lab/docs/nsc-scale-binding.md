# Binding the charged throat radius to the inherited spectral scale

This calculation joins three equations that previously remained separate: the
relational law fixes the homogeneous coefficient to (V_{\rm full}=0), the
compact spectral weight supplies Wilsonian (A) and (C), and the imported
Maldacena–Milekhin–Popov (MMP) relation fixes the charged throat radius.  It
does not rerun the MMP geometry or reuse the matching cutoff as a fitted scale.

## Keep the charged light field explicit

The MMP source is the charged massless lowest-Landau sector.  It must therefore
remain in the causal theory rather than also being integrated into the local
coefficients.  MMP specify that the running gauge coupling is evaluated at the
magnetic scale

$$
\mu_B=\frac{\sqrt{|q|}}{r_e}.
$$

Using the already verified compact/light allocation, the coefficients above
that scale are

$$
A(\Lambda,\mu_B)=\frac{Q_1[h_\Lambda]-\mu_B^2}{6(4\pi)^2},
\qquad
C(\Lambda,\mu_B)=\frac{H_{\Lambda,\mu_B}(0)}{3(4\pi)^2}.
$$

The explicit light field supplies its Casimir stress and causal response.  It
is counted once.  The relational projector acts after the same-action assembly
and gives (V_{\rm full}=0); it does not remove (A), (C), the Casimir term,
or the boundary response.

## One radius equation becomes one scale equation

Set the throat radius as the room ruler, (r_e=L_\star=1).  The parent
resolution is (1/r_e), so the child cutoff gives

$$
\Omega=\Lambda r_e,
\qquad \zeta=(\Lambda L_\star)^2=\Omega^2.
$$

The imported MMP normalization

$$
r_e^2=\frac{q^2C}{4A}
$$

then becomes the dimensionless root equation

$$
\boxed{F_q(\Omega)=\frac{q^2C(\Omega,\sqrt{|q|})}
{4A(\Omega,\sqrt{|q|})}-1=0.}
$$

For each calculation (q) is fixed as an integer before the continuous root
is solved.  No observational quantity sets the bracket or the result.

## The first cutoff-resolved development branch

The adjacent fixed sectors (q=2,3,4) give one radius root each above their
magnetic matching scales.  The first two roots lie below the stored child's
asymptotic curvature scale (H_{\rm child}L_\star=\sqrt{3\pi}), so the same
finite cutoff cannot resolve the geometry to which it would be applied.  The
first checked integer passing that independent condition is (q=4):

$$
\Omega=3.973074368754324,
\qquad
\zeta=15.78531993965257,
$$

$$
A L_\star^2=0.04501936182826088,
\qquad
C=0.011254840457065235,
\qquad
V_{\rm full}L_\star^4=0.
$$

The MMP radius residual is (1.33\times10^{-15}), and
(H_{\rm child}/\Lambda=0.7726963653).  Increasing the compact basis from 26
to 34 modes changes (Omega) by (8.64\times10^{-9}).  The conditional old
scalar-gap expression is real on this branch, but it is not applied to the
massless charged MMP channel.

At this value the stored smooth-chain endpoint control has
(\Omega q_{\rm transfer}<1), so that particular carrier remains in its
limit-point regime.  This is a compatibility check, not an identification of
the smooth periodic carrier with the charged throat.

## What this closes and what it leaves

The calculation supplies a concrete fixed-(q) geometry/scale candidate and
removes the arbitrary light/heavy split from its coefficient equation.  It
does not select magnetic flux from the complete action; (q=4) is the first
adjacent development branch satisfying the stated cutoff condition, and it is
not in MMP's parametrically large-(q) control regime.  The actual charged
recursive boundary state and the four-component CTP metric source still have
to accept the same branch before it becomes a `PhysicalSolution`.

The focused reproducer is

```sh
python3 scripts/derive_nsc_scale_binding.py --check
```

It evaluates only this new coefficient/root binding and one resolution
increase.  It does not rerun the MMP Einstein equations, the angular vacuum
source, the historical scale scans, or any collapse calculation.

