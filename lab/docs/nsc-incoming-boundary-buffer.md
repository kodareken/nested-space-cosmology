# Incoming value control with a longer computational buffer

The exact harmonic source-phase control reduced the reference matter drift
by about590 times. Its remaining transient peaks at the incoming slice after
the short original numerical right-boundary travel time. Reuse the original
Dirac/SBP/SAT operator and source; extend only the numerical edge from rho1.2
to1.8. This tests that identified numerical error, not a different state law.

## Preparation and phase

Preserve all3201 authenticated reference-grid nodes and columns byte-for-byte;
the original801 grid is their exact spin-wise nested subset. Append only
rho>1.2 using one dense DOP853 solve of the owned canonical radial generator

$$
G_E=\frac{i}{a_0}\left[-m\sigma_1+\frac{\ell}{r}\sigma_2-
\frac{E}{a_0}\sigma_3\right].
$$

The initial amplitude is the actual three-fiber KS restriction at1.2. The
inverse includes the fixed trace/half-density map and exp(-i E S(rho)); no
normalization, source occupation, reflection phase or scattering is changed.
Both new grids sample this same continuation. The extended coordinates retain
the original prefix bytes, including the original linspace endpoint.

For the reference geometry, reuse the harmonic augmentation

$$
\Phi'=L_0\Phi+BQ,\qquad Q'=-i\operatorname{diag}(E)Q,
\qquad Q(0)=I.
$$

At rho1 the incoming coordinate is z=tau+S(1). Restrict the actual fields and
their PDE derivatives: F=R Phi and F_z=R(L_0 Phi+BQ). The coherent matter
contraction uses the original source matrix, weights once and multiplicity12.
It is compared with the stationary analytic reference without subtracting
the measured drift. This four-energy group14_1 check does not cover the full
source and does not constrain the incoming covariance to be C0.

## Bounded experiment and decision

Run exactly two reference evolutions,951 and3801 nodes on[-2,1.8], with65
samples on[0,.3]. Continue the new .6 radial interval once (rtol2e-13,
atol2e-15,max_step1/2000). Total CPU cap60seconds. Refuse a new run when its
output exists; --check is array-only. The fine-grid raw N,beta drift target
is3e-11; the coarse run measures the spatial change, without requiring that
coarse data meet the fine target. Check phase/flux/frame/source-fiber algebra
and exact inherited prefixes. A timeout or failed target is OPEN; no extra
case starts automatically. Record both attempted and completed field runs.

For rho>0, rho*arctan(1/rho)<1 in the owned beta derivative, so beta decreases.
On this reference interval the incoming fast-characteristic
travel time is at least .8/(1+beta(1))>.3. This is a continuum buffer argument.
It does not make the finite SBP grid exactly causal, nor bound its continuum
field error. The remaining low/subgap, all-family, changed-history and
continuum uncertainties stay OPEN regardless of this numerical control.

The next experiment uses the nonzero prepared history and its retarded
tangent with these same buffered initial columns. No global matching,
extended stationarity or metric timestep follows from this control.

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_boundary_buffer.py
python3 scripts/derive_nsc_incoming_boundary_buffer.py --run
python3 scripts/derive_nsc_incoming_boundary_buffer.py --check
```
