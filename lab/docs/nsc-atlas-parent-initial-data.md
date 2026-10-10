# A shared spatial parent from Brill-Lindquist data

The atlas now has one consistent two-puncture spatial metric instead of two
unrelated local Schwarzschild charts. It is imported initial-data mathematics.
It supplies neither child formation nor a solved parent-child junction.

Owner: [adapter](../src/recursive_horizons/nsc_atlas_parent_initial_data.py),
[independent controls](../tests/test_nsc_atlas_parent_initial_data.py).
The existing [atlas](nsc-cosmological-atlas.md) keeps its original sealed record;
a successor producer may embed `parent_record()` from this adapter.

## Borrowed construction and reading boundary

[Pook-Kolb et al., arXiv:2308.10906v1](https://arxiv.org/html/2308.10906v1),
section III, equations (17)-(21) and the surrounding initial-data/MOTS
paragraphs were read directly on 10 October 2026. The original attribution is
[Brill and Lindquist, Phys. Rev. 131, 471 (1963)](https://doi.org/10.1103/PhysRev.131.471);
that original article was not separately read in this bounded review. The
quasi-local-mass analysis and numerical merger evolution in the newer article
are not imported or reproduced here.

With $G=c=1$, positive bare masses $m_i$ and distinct punctures $x_i$,

$$
h_{ab}=\psi^4\delta_{ab},\qquad K_{ab}=0,\qquad
\psi(x)=1+\sum_i\frac{m_i}{2|x-x_i|}.
$$

The vacuum Hamiltonian constraint is
$R[h]=-8\psi^{-5}\Delta_{\rm flat}\psi=0$ away from punctures; the momentum
constraint vanishes for $K_{ab}=0$. Spatial scalar curvature zero does not
mean a flat spatial metric or vanishing tidal curvature. The adapter evaluates
the analytic first and second derivatives of the imported conformal factor.

The outer asymptotic end has $M_\infty=\sum_i m_i$. The end at puncture $i$
has $M_i=m_i[1+\sum_{j\ne i}m_j/(2d_{ij})]$. These are asymptotic-end charges,
not individual apparent-horizon or other quasi-local masses. For the recorded
inputs $m_A=1$, $m_B=1.5$, $d=12$, the three charges are $2.5$, $1.0625$ and
$1.5625$. Their difference $M_\infty-M_A-M_B=-0.125$ uses this explicitly
stated convention. Coordinate separation is not proper distance.

The spatial topology is $\mathbb R^3$ minus two punctures, with three
asymptotically flat ends. Those additional ends are not identified with newly
formed expanding children, and an initial spatial connection supplies no
traversable causal route. No apparent-horizon/minimal-surface search is run;
the adapter does not claim to exclude a common horizon. A finite-separation
two-hole setting cannot be assigned one global spherical symmetry.

## What a local spherical approximation loses

On a coordinate sphere of radius $r<d$ about puncture $i$, orient
$\mu=\cos\theta$ toward the other puncture $j$. The companion contribution is

$$
\frac{m_j}{2\sqrt{d^2+r^2-2dr\mu}}
=\frac{m_j}{2d}\sum_{\ell=0}^{\infty}(r/d)^\ell P_\ell(\mu).
$$

This convergent analytic identity is imported elementary multipole mathematics,
not an infinite numerical work queue. The implemented profile uses the exact
finite-distance expression at 17 angular samples. The monopole comparison is
$\psi_0=1+m_i/(2r)+m_j/(2d)$. The dipole is a coordinate gradient; the
quadrupole coefficient in $\psi$ is $m_jr^2/(2d^3)$. Neither is labelled a
measured four-dimensional tide: no observer or lapse has been supplied.

The sphere's induced metric is
$q_{AB}=r^2\psi(\mu)^4\Omega_{AB}$, so its relative angular-metric correction
is $(\psi/\psi_0)^4-1$. Its intrinsic Gaussian curvature is

$$
\mathcal K=\frac{1-2\Delta_{S^2}\log\psi}{r^2\psi^4},\qquad
\Delta_{S^2}g=(1-\mu^2)g''-2\mu g'.
$$

The adapter records angular profiles of both quantities and the maximum
relative conformal-factor correction; it does not invent a threshold for an
"acceptable" junction. The local sphere radii are the respective bare masses,
$r_A=1$, $r_B=1.5$. They are coordinate probes, not computed horizons or child
interfaces. Their angular variation prevents interpreting the monopole
approximation as an exact Schwarzschild chart at finite separation. A genuine
attachment would need compatible induced geometry and extrinsic curvature,
plus a physical source and an evolution that defines its timelike history.

## Finite admission and verification

This is the integration exception under the existing
[admission rule](../../docs/prior-art-and-open-claim.md#compute-admission-rule):
known formulas evaluated to improve the map, not a research novelty claim.

| Required item | Recorded decision |
|---|---|
| Borrowed | Published conformal initial data and asymptotic charges above; ordinary conformal and multipole identities |
| Missing connection | The first atlas's two local holes lacked a shared spatial parent; child interface and spacetime evolution remain absent |
| Decision and unlock | Supply a concrete nonspherical parent and show what the local spherical approximation drops, enabling a later specified interface test |
| Domain | One time-symmetric vacuum slice, two positive unequal masses at finite separation; four off-puncture controls and two 17-point angular profiles |
| Budget and stop | At most one CPU minute for evaluation and focused tests; no PDE or horizon search. Reject singular samples and local-sphere expansion domains with $r\ge d$; failed constraint/unit controls block this adapter |
| Independent check | Finite-difference gradients/Hessian and Laplacian, outer ADM flux quadrature, puncture-end inversion, harmonic sphere mean-value identity and common unit rescaling |

All masses and coordinates scale together under a unit change. Dimensionless
angular profiles remain unchanged, gradients scale inversely with length,
Hessians with its square, and ADM masses scale as lengths. Fixed-slice data do
not set lapse, elapsed proper time, inherited constants or branching rates.

The API returns NumPy field arrays for finite batches. `parent_record()` is
deterministic and JSON-safe, carries source/version and missing conditions,
and does not write files. The atlas producer owns output versioning and hashes.
