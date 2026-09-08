# Nested-Space closure verification — 7 September 2026

This follow-up checks the mathematics behind the curated `v0.1.0` snapshot.
The compact JSON and this note are imported as the 59th public record. The
original 58 JSON records and their generators are retained byte-for-byte.
It preserves the constructive objective: one common operator must determine
the particle, outside-response, and parent/child sectors. It records one new
calculation, two normalization issues, the strongest surviving identities,
and the precise equations still needed for physical closure.

Source checkpoint: `ff2cf2722b966589b98a61accdbb6cee819a58c7`.
The curated snapshot is `c0e4fee30c4cbd53d67ef67c9e1c971b675ab13d`.
Its successful CI establishes reproducibility under its comparison policy.
It does not independently derive the implemented equations. All historical
JSON and generator bytes are retained unchanged by this follow-up.

## Result of the calculation

The ZETA1 mass term uses inconsistent units. Its definitions are

$$
\lambda_j=L_\star^2\,\operatorname{eig}_j(D_{\rm spatial}^2),\qquad
\zeta=\Lambda^2L_\star^2,\qquad
\mu^2=\frac{|\Phi|^2}{\Lambda^2}=1-\frac{a}{\zeta},\qquad a=\frac{3\pi}{2}.
$$

Consequently the same declared additive-gap operator requires

$$
\boxed{q_j=\frac{\lambda_j}{\zeta}+\mu^2
       =1+\frac{\lambda_j-a}{\zeta}.}
$$

The historical runners instead used `(lambda_j + mu_squared)/zeta`.
That expression implements a physical mass squared smaller by `1/zeta`
than the declared mass squared. This is a conversion error, independent of
the interpretation of the nested-space hypothesis. The source is explicit
in `run_nsc_zeta1_regulated_determinant.py::_mode_values`,
`run_nsc_zeta1_lowest_mode.py`, and the compact-direction successors.

The new [calculation](../scripts/check_nsc_scale_closure.py) independently
assembles the same radial potentials and compact warp. It computes exact
symbolic unit identities, reconstructs the old derivative as a control,
then changes only the mass conversion. It retains radius 30, 750 radial
half-intervals, angular sectors 1–12, 32 compact intervals, and the historical
factorization and angular weights. Those approximations remain explicit.

At the old warped candidate `zeta = 4.748389947082489`:

| Quantity | Historical argument | Argument with consistent units |
|---|---:|---:|
| Determinant logarithmic derivative | approximately 0.0120152 | approximately -24.31995 |
| Derivative after the inherited cutoff subtraction | approximately -12.32458 | approximately -36.58293 |

The independent old-argument reconstruction differs from the stored derivative
by `1.53e-8`; separately generated half-grids and a literal principal
compression differ at floating-point rounding level. With radial spacings
`0.08`, `0.04`, and `0.02`, the corrected subtracted derivative at this point is
`-36.58813190`, `-36.58292833`, and `-36.58162543`. Its sign is resolved.

The corrected finite determinant has a diagnostic root near
`zeta = 6.09675392014`. At that root, the subtracted derivative is approximately
`-14.15469`. It therefore does not solve the adopted subtracted scale equation.
No physical scale is inferred from this new diagnostic root.

The [raw follow-up result](../results/nsc-2-zeta1-unit-closure-check.json)
contains all 257 scan points, the three radial resolutions, exact algebra,
source hashes, and the noncommuting block-matrix control.

## A sign result stronger than the scan

For the inherited prescription, let

$$
\Gamma_{\det,j}=\tfrac12 E_1(q_j),\qquad
q_j=1+(\lambda_j-a)/\zeta.
$$

At fixed spatial geometry, subtracting the same pure-cutoff term used in the
historical calculation gives

$$
\boxed{\frac{d\Gamma_{{\rm sub},j}}{d\log\zeta}
       =-\frac{e^{-q_j}}{2q_j}.}
$$

This identity is verified symbolically, and the determinant derivative is
checked against finite differences. The disconnected radial matrix is the
principal submatrix obtained by deleting the throat node from the joined
matrix. Let its eigenvalues be `nu_j`, with joined eigenvalues `lambda_j`.
Hermitian eigenvalue interlacing gives

$$
\lambda_j\le\nu_j\le\lambda_{j+1}.
$$

Since `f(lambda)=exp(-q(lambda))/q(lambda)` is positive and decreasing,

$$
\operatorname{Tr}f(H_{\rm joined})-
\operatorname{Tr}f(H_{\rm disconnected})>0.
$$

Every angular and compact degeneracy is positive. Thus the total subtracted
derivative is strictly negative for every `zeta > 3*pi/2` in this finite,
fixed-geometry, additive-gap family. The numerical scan corroborates the sign;
it is not the basis for extending the statement between sample points.

**Consequence:** refining this family or selecting another determinant-only
minimum cannot close its scale equation. A further attempt needs a derived
change in the physical operator, its geometry/link dependence on scale, or
the anomaly prescription. Merely naming the omitted recursive tail does not
prove that it supplies that change. This result does not cover the full
noncompact operator, a general recursive geometry, or the hypothesis as a whole.

## Two further connections that must be made explicit

### Raw determinant and regulated determinant

For the very same relative-sheet symbol,

$$
D_\delta=\begin{pmatrix}pe^{-\delta}&\Phi\\\Phi&-pe^{\delta}\end{pmatrix},
\qquad \det D_\delta=-(p^2+\Phi^2),
$$

the unregulated finite determinant is independent of `delta`. But the
proper-time functional actually used later in ZETA1 satisfies

$$
\boxed{
\left.\frac{\partial^2}{\partial\delta^2}
\left[\frac12\sum_\pm E_1(E_\pm(\delta)^2/\Lambda^2)\right]
\right|_{\delta=0}
=\frac{4p^2}{\Lambda^2}e^{-(p^2+\Phi^2)/\Lambda^2}>0
\quad (p\ne0).
}
$$

The new script verifies this equality exactly. Therefore raw determinant
invariance cannot justify setting the regulated relative Hessian to zero.
The existing full-exponential calculation remains evidence about its stated
heat-trace Hessian. Identifying it with the complete physical anomaly or
covariance requires a common regulator, field measure, and derived variation.
The primary anomaly paper distinguishes its normalization scale from the
spectral cutoff explicitly; it is not a licence to interchange regularizations.
[Andrianov, Kurkov and Lizzi, Sections 3–5](https://arxiv.org/html/1106.3263v1)

### Normalization of the recursive link

At a fixed common dimensional energy, `x=E/Lambda_parent` is seen as
`x/Omega` by a child with `Lambda_child=Omega*Lambda_parent`. For a first-order
inverse response, the child's dimensionless block in parent units is
`Omega*Gamma_child(x/Omega)`. If `b=B_dim/Lambda_parent`, block inversion gives

$$
\boxed{\Gamma_p(x)=K_p(x)-\frac1\Omega\,
 b\,\Gamma_c(x/\Omega)^{-1}b^\dagger.}
$$

Using `b_sym=B_dim/sqrt(Lambda_parent*Lambda_child)` absorbs the factor and
recovers the displayed unweighted form. Both are legitimate conventions;
they cannot silently share the same numerical link. The new script verifies
the scalar normalization and a noncommuting complex-energy block-matrix
example. It also separates a transformed mode, for which `E_child/Lambda_child`
is unchanged, from probing both rooms at one fixed dimensional energy.
These statements do not choose `Omega` or specify a clock map in a general
time-dependent geometry.

## Strong points and discovery status

| Result | What is supported | What can be claimed about originality |
|---|---|---|
| One `Phi` gives the mass gap and visible self-energy | Exact block algebra: `D^2=(p^2+Phi^2)I` and `G_pp^-1=E-p-Phi^2/(E+p)`; both rechecked independently here | A sound realization of the desired connection. Two-sheet Dirac/Higgs geometry and Schur elimination predate this project. |
| One off-diagonal field also couples two metrics | The published doubled spectral action supplies this relationship | Imported result; explicit in Bochniak–Sitarz. |
| Regular black-universe geometry | A finite mathematical benchmark with an expanding interior exists in the cited action | Imported construction; no new existence priority or observational ancestry claim. |
| Explicit warped carrier and its `1015/144` coefficient relation | Algebraic statements for the selected ansatz are recorded; the critical interface has unresolved physical stability | A model-specific candidate contribution. A literature search and a reproducible value alone establish neither uniqueness nor a new healthy solution. |
| Relative heat kernel at finite momentum | Positive sampled Euclidean Hessian values are computed for the declared truncation | Potential numerical study, subject to the regulator distinction above. Priority over all prior work is unestablished. |
| New unit correction and finite-family sign result | Exact normalization plus a scope-specific monotonicity result, with numerical corroboration | New findings in this repository. They are not evidence of a new law of nature or of new-to-world priority. |

The scalar equation `Gamma=K-b^2/Gamma` is also established surface-Green-function
mathematics: a semi-infinite chain gives the same recursion and square root.
Its value here is the proposed physical identification of its operators, which
remains to be solved. [Nemec, Tomanek and Cuniberti, Appendix A.3](https://arxiv.org/pdf/0711.1088)

The two-metric action and its shared field are explicit in
[Bochniak and Sitarz (2022), Section 2](https://arxiv.org/html/2201.03839v1).
Their FLRW reduction already gives `Lambda_e=12*(Lambda^2/c-kappa*abs(Phi)^2)`
and `alpha=12*kappa*abs(Phi)^2`.
[Bochniak and Sitarz (2021), Eq. II.33](https://arxiv.org/pdf/2012.06401)

The regular black-universe benchmark and its phantom-source assumptions are
documented in [Bronnikov, Dehnen and Melnikov (2007)](https://arxiv.org/pdf/gr-qc/0611022).
The critical EGB phenomenon also has existing literature; the critical vacuum
alone supplies no healthy graviton proof.
[Fan, Chen and Lu (2016)](https://arxiv.org/pdf/1606.02728)

This was a targeted comparison with the closest primary sources, not an
exhaustive priority review. Absence of a search hit cannot certify novelty.

## Remaining physical closure, in dependency order

1. **Specify one regulated curved operator and one field measure.** Apply the
   corrected mass units and explicit parent/child normalization. Derive its
   scale and relative-field variations from the same functional. The physical
   stationarity equation must be derived before another root is interpreted.
2. **Compute the throat response from that operator.** The current domain
   artifact computes radial partner identities and a two-component flux
   cancellation example. `N_parent(E)`, `N_child(E)`, the APS projectors,
   lapse/shift additions, and `Phi_throat=R(E)` are recorded as symbolic
   declarations, not evaluated boundary maps. The cited BFK theorem assumes
   compactness, a product collar, and spectral conditions; its direct
   applicability to this noncompact warped geometry is unproved.
   [Lee (2003), Theorems 1.1 and 1.5](https://arxiv.org/html/math/0304347v2)
3. **Solve the normalized recursive boundary problem and stationary scale
   together.** Demonstrate that the same operator produces its own boundary
   data with a controlled tail error. Only that result may determine `Omega`,
   `zeta`, and the physical gap. Independent link fitting is excluded by the
   model's own one-parameter-set objective.
4. **Establish physical propagation and nonlinear continuation.** Use the
   gauge-invariant metric/field Hessian, retarded boundary conditions, conserved
   charges, and full stress tensor. A positive two-component correlator does
   not establish the graviton's health or a healthy collapse trajectory.
5. **Make a shared, unfitted prediction.** The existing coefficient audit
   imposes three constraints on four Skyrme coefficients and leaves six other
   coefficients unclosed. Its rank calculation is not a computation of the
   full renormalization-group relevant directions. Particle masses, nuclear
   binding, lensing, expansion, and the parent/child energy account must still
   be computed from the same specified parameter vector and compared with
   data not used to choose it.

The existing 58-step record is a collection of useful exact reductions and
finite diagnostics. It is not yet the continuous, coupled physical solution
needed to say the full mathematics and physics work. The strongest remaining
target is still a quantitative relation forced by one operator across two
physical scales.

## Reproduce this follow-up

```bash
python3 scripts/check_nsc_scale_closure.py --check
```

The command compares every exact field and every numeric field in the new
record; it does not apply the public snapshot's selective numeric comparison.
It does not launch a campaign or modify any historical output. The script
contains no imported runner decisions, no observed particle/cosmological fit,
and no action retuning. The source anchor and input SHA-256 values are checked
before and after the calculation.

## Consequence for the curated PDF

The `v0.1.0` PDF is a historical presentation of the previous checkpoint. Its
ZETA1 scale numbers must now be described as diagnostics of an incorrectly
normalized additive-gap proxy. Its raw determinant argument must be scoped
to the raw determinant, and the throat identification must remain a proposed
operator relation until the boundary maps are solved. Those corrections are
required before the PDF can support stronger physical or novelty claims.
