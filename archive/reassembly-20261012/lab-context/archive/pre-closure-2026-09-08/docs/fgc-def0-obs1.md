# FGC-1-DEF0-OBS1: exact metric-null observable contract

**FGC-1-DEF0-OBS1** defines and tests the physical-metric observable required
by the later DEF1 gate. It is a geometry contract with exact point controls,
not a collapse or defocusing result.

For

\[
ds^2=-\alpha^2dt^2+\Lambda^2(dr+vdt)^2+R^2d\Omega^2,
\]

the future unit normal and outward radial unit vector are

\[
n^a=\alpha^{-1}(\partial_t-v\partial_r)^a,
\qquad s^a=\Lambda^{-1}(\partial_r)^a.
\]

The unscaled radial null frame

\[
\ell_\pm^a=n^a\pm s^a
\]

is checked directly against the physical metric:

\[
g(\ell_\pm,\ell_\pm)=0,
\qquad g(\ell_+,\ell_-)=-2.
\]

Its round-sphere expansions are

\[
\vartheta_\pm=\frac{2}{R}\ell_\pm(R).
\]

The observable accepts an exact positive rescaling \(k^a=\nu\ell_\pm^a\)
and rejects it unless

\[
g(k,k)=0,
\qquad k^b\nabla_bk^a=0
\]

hold exactly. It then evaluates the expansion along that affine generator in
two independent ways. The trajectory-jet derivative is compared with the
complete spherical, hypersurface-orthogonal Raychaudhuri right-hand side

\[
\frac{d\theta}{d\lambda}
=-\frac12\theta^2-\sigma_{ab}\sigma^{ab}
-R_{ab}k^ak^b.
\]

For a radial spherical congruence the screen shear vanishes by symmetry; the
selected radial null hypersurface has zero twist. Those are declared
specializations, not non-spherical robustness claims. The output stores the
expansion, shear, twist, and Ricci terms separately so that
\(-R_{kk}>0\) cannot be mislabeled as defocusing when the complete right-hand
side remains nonpositive.

## Exact controls

1. Spherical Minkowski space gives a normal round sphere and the exact
   outgoing focusing identity
   \(d\theta/d\lambda=-1/8\) at \(R=4\).
2. A contracting flat de Sitter chart supplies a trapped round sphere while
   the affine radial congruence still focuses; de Sitter has \(R_{kk}=0\).
3. A synthetic smooth local warped geometry has \(\theta=0\) and
   \(-R_{kk}=1/2\), so the complete right-hand side is positive. This verifies
   the sign-detection path. It is not an ACT1 field-equation solution and is
   not evidence for the FGC mechanism.

## What DEF1 still requires

A true DEF1 record must consume a constraint-compatible, HYP1- and EFT-valid
COL1 evolution; integrate an affinely normalized null generator through a
finite trapped interval; demonstrate a positive lower error bound for the
complete right-hand side at multiple resolutions; and show constraint,
affine, null, interpolation, gauge, and boundary errors are smaller than that
margin. OBS1 supplies none of those data. It derives no transition surface,
physical wall, or singularity resolution.

Reproduce with:

```bash
python3 scripts/reproduce_fgc_def0_obs1.py
```
