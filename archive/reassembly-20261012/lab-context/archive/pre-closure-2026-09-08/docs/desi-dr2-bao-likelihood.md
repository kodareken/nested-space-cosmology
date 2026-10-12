# DESI DR2 BAO-only constant-w profile

This artifact is a deliberately narrow empirical calculation: a real Gaussian
likelihood for the public 13-observable DESI DR2 BAO distance vector. It tests
only a flat late-time constant-`w` expansion history with a free positive
`H0*r_d` scale. It is not a DESI+CMB or DESI+supernova analysis, not a
full-shape calculation, and not evidence that dark energy has an external or
boundary origin.

The raw mean vector, covariance, official chain configuration, marginal
summaries, and best-fit record are held in `collected-data/desi-dr2-bao/`.
[`scripts/reproduce_bao.py`](../scripts/reproduce_bao.py) records their sizes
and SHA-256 hashes in
[`results/desi-dr2-bao-profile.json`](../results/desi-dr2-bao-profile.json). The upstream likelihood configuration is
pinned by the DESI-provided `b7b8a36e9bccb063081f811f323cada21ab5fbdd` revision identifier recorded in that
JSON; the files must remain a matched set.

## Model and likelihood

For each point with redshift `z <= 2.33`, this limited calculation uses

\[
E^2(z)=\Omega_m(1+z)^3+(1-\Omega_m)(1+z)^{3(1+w)}.
\]

Radiation is deliberately omitted: it is a sub-percent late-time
approximation over this redshift range, not an assertion that radiation is
absent from cosmology. The line-of-sight and transverse distances are

\[
I(z)=\int_0^z\frac{dz'}{E(z')},\qquad
\frac{D_M}{r_d}=A I(z),\qquad
\frac{D_H}{r_d}=\frac{A}{E(z)},
\]

with `A=c/(H0*r_d)`. The isotropic observable is

\[
\frac{D_V}{r_d}=\left[z\left(\frac{D_M}{r_d}\right)^2
\frac{D_H}{r_d}\right]^{1/3}.
\]

At fixed `(Omega_m,w)`, `A` is minimized analytically using the supplied full
covariance,

\[
A_*=\frac{f^T C^{-1}d}{f^T C^{-1}f},\qquad
\chi^2=(d-A_*f)^T C^{-1}(d-A_*f),
\]

and must be positive. The remaining nuisance-free profile is searched with
deterministic bounded golden searches, fixed in this artifact at
`Omega_m in [0.05,0.7]` and `w in [-2,0]`. Composite Simpson integration uses
256 even intervals per distance integral; the generated result repeats the
fit at 128 intervals and records the difference. One-parameter `w` profile
intervals use the conventional asymptotic thresholds `Delta chi2=1` (68%) and
`3.841458820694124` (95%). These choices are transparent and fixed for the
artifact, but they are not a prospective preregistration made before DESI
published the target data.

Run it with:

```bash
python3 scripts/reproduce_bao.py
python3 -m unittest tests/test_bao.py -v
```

## What the result can and cannot say

The calculation can reject a pre-fixed *constant-w background subclass* when
its profile likelihood is poor. It provides a stronger test than comparing a
proposal to a rounded posterior mean, because it uses the thirteen distances
and their correlations. It still cannot identify the physical origin of an
effective component: an ordinary fluid, a scalar field, an interaction, and a
boundary-inspired model can share the same background `H(z)`.

It does not infer `H0` itself. BAO alone constrains `H0*r_d`; inferring `H0`
requires an independently modeled/calibrated sound horizon and early-universe
physics. It also cannot apply the separately conserved DEB-1 relation
`w=-1+qs/3` to a flux-fed `Q != 0` component. Such a model requires a
covariant flux/junction law and, for an origin claim, derived perturbation,
lensing, and growth predictions tested against independent data.

Primary context: [DESI DR2 Results II: BAO Measurements and Cosmological
Constraints](https://arxiv.org/abs/2503.14738). The public calculation here is
not a replacement for the collaboration's full CAMB implementation, priors,
posterior sampling, or complete analysis.
