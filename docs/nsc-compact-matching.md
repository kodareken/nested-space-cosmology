# Match the compact Dirac determinant to its light field

The calculated warped determinant now has a finite matching map to one
canonical four-dimensional Dirac field. The same spectral weight determines
its complementary vacuum, Einstein, gauge and curvature-squared coefficients.
They are no longer unspecified pieces of this particular free determinant.
Contributions from the complete measure, compensator, physical state and
other interactions remain distinct unknowns.

This is an assembly step using the [warped-source result](nsc-warped-source.md),
the [known Dirac heat coefficients](nsc-covariant-source.md), and the existing
[normalized state functional](nsc-influence.md). It does not repeat those
calculations or introduce a second gravitational action. The new quantities
are two spectral moments and the correctly normalized low-argument limit.

## Separate a field, not a five-dimensional energy interval

For the existing paired chiral compact domain, write

\[
\Gamma_{5,\Lambda}=\tfrac12\operatorname{Tr}_4 h_\Lambda(D_4^2),\qquad
\Gamma_{\mathrm{light},\nu}
 =\tfrac12\operatorname{Tr}_4 E_1(D_4^2/\nu^2),
\]
\[
H_\nu(y)=h_\Lambda(y)-E_1(y/\nu^2),\qquad
\Gamma_{H,\nu}=\tfrac12\operatorname{Tr}_4 H_\nu(D_4^2).
\]

Thus Gamma5=Gamma_light+Gamma_H exactly within the stated Euclidean
prescription. The light functional belongs to the one canonical massless
4D Dirac field already supplied by the opposite-parity compact pair. Its
subtraction also moves its high-frequency contribution into Gamma_H.
This is not a sharp energy projector on the warped 5D singular values,
nor is nu a selected particle mass, a measured coupling or a determinant
normalization mass. Lambda remains the fixed cutoff of the original
functional. A change of this partition cannot change its total source.

## An exact proper-length limit fixes the finite gauge coefficient

Let sigma(Y) be the already specified Gaussian and define

\[
J=\frac1\ell\int_{-\ell/2}^{\ell/2}e^{s\sigma(Y)}\,dY,
\qquad \ell_{\rm proper}=\ell J.
\]

At zero tangential eigenvalue, the proper coordinate
z=int exp(s sigma) dY makes the compact sandwich Dirac operator a flat
first-order derivative on this proper interval. Its modulus has one zero
mode and pairs of nonzero eigenvalues squared (p pi/ell_proper)^2.

The zero singular mode has unit norm with respect to exp(s sigma)dY.
The matrix element of a small tangential lambda between the normalized
left/right zero modes is lambda ell/ell_proper=lambda/J. Consequently

\[
\epsilon_0(y)=y/J^2+O(y^2),\quad y=\lambda^2,
\]
\[
\boxed{H_\nu(0)=\log\frac{\Lambda^2J^2}{\nu^2}
 +2\sum_{p\geq1}E_1\!\left[
       \left(\frac{p\pi}{\ell_{\rm proper}\Lambda}\right)^2\right].}
\]

The warp-gradient mixing retained in the previous calculation is essential
to the small-y slope. A diagonal-only expectation of the quadratic form
would give a different ratio of warp integrals. The new record compares
the proper-length spectrum and this limit with that existing full matrix.
The light canonical mass is still zero: this is a regulator-norm matching
factor, not a generated fermion mass or a changed classical KK tower.

The logarithms cancel in H. An actual discrete spacetime zero mode still
requires its separate determinant/state prescription. The preceding
noncompact R2 density has no L2 kernel contribution; that condition is not
silently transferred to a compact spacetime box.

## Apply the existing heat coefficients once

Define the four-dimensional Mellin quantities

\[
Q_1[H]=\int_0^\infty H(y)\,dy,\qquad
Q_2[H]=\int_0^\infty yH(y)\,dy,\qquad Q_0[H]=H(0).
\]

They have dimensions length^-2, length^-4 and one, respectively. They are
not the half-integer 5D proper-time moments in the earlier leading UV match.
The light field subtracts nu^2 from Q1[h] and nu^4/2 from Q2[h].

Reuse the trace and Mellin identities in
[Codello–Percacci–Rahmede, Appendix A, equations A10 and A14–A15](https://arxiv.org/pdf/0805.2909v5), and
[Vassilevich, equations4.26–4.28](https://arxiv.org/html/hep-th/0306138v3).
For one canonical 4D Dirac field the retained coefficients are
a0=4, a2=-R_E/3 and
a4=(-18 C_E^2+11 E4-12 box_E R_E)/360+(2/3)F^2.
The factor one-half in Gamma_H then gives the local Euclidean density

\[
\mathcal L_{H,E}=V_D-A_D R_E+C_{F,D}F^2
 +C_{W,D}C_E^2+C_{E,D}E_4+C_{\Box,D}\Box_E R_E+\cdots,
\]
\[
V_D=\frac{2Q_2}{(4\pi)^2},\qquad
A_D=\frac{Q_1}{6(4\pi)^2},\qquad
C_{F,D}=\frac{Q_0}{3(4\pi)^2},
\]
\[
C_{W,D}=-\frac{Q_0}{40(4\pi)^2},\quad
C_{E,D}=\frac{11Q_0}{720(4\pi)^2},\quad
C_{\Box,D}=-\frac{Q_0}{60(4\pi)^2}.
\]

There is no independent R_E^2 term in this free Dirac a4 combination.
That derived zero does not set the complete theory's allowed finite R^2
coefficient to zero. These are contributions in the declared scheme, not
predictions of the complete Newton constant or gauge coupling.

The compact warp is retained in h, including its finite effects. Truncating
the resulting g4 action at a4 still requires small external curvatures,
derivatives and gauge fields compared with the eliminated scales. This
condition has not been established for the imposed NSC neck. Use the full
spectral remainder wherever that local expansion is not controlled.

At Lambda=2, ell=2, s=1 and nu=1, in natural units:

| Quantity | Calculated value |
|---|---:|
| J | 0.99411998706 |
| H_nu(0) | 2.29160403721 |
| Q1[H] | 5.21331560062 |
| Q2[H] | 13.62328050962 |
| V_D | 0.17254086329 |
| A_D | 0.00550228445 |
| C_F,D | 0.00483725070 |
| C_W,D | -0.00036279380 |

The values identify the numerical development record and its chosen units;
they do not select a physical cutoff, compact size or matching scale. From
26 to 34 compact basis levels Q1 changes by about 1.79e-8 and Q2 by 7.04e-8.
These are observed Galerkin convergence measures. Separate momentum-tail
bounds do not cover retained-eigenvalue errors. Independent analytic flat
moments check the new integral and its normalization without rerunning the
old compact-spectrum generator.

## Matching-scale changes cancel; they cannot be fitted

For fixed physical operator and cutoff,

\[
\partial_\nu Q_1[H]=-2\nu,\quad
\partial_\nu Q_2[H]=-2\nu^3,\quad
\partial_\nu Q_0[H]=-2/\nu,
\]
\[
\partial_\nu\Gamma_{H,\nu}
 =-\frac1\nu\operatorname{Tr}_4e^{-D_4^2/\nu^2}
 =-\partial_\nu\Gamma_{\mathrm{light},\nu}.
\]

The record checks the finite vacuum, Einstein and gauge changes at nu=0.5,
1 and1.4. In the gauge sector this compares changes in the light functional;
it does not assign its massless, nonlocal response an arbitrary finite
local constant. A ratio formed from V_D,A_D,C_F,D alone is generally
partition-dependent. Varying nu to make such a ratio fit a geometry would
omit the compensating light contribution.

## Connect the source owner and identify the remaining equation

The existing static metric owner uses the energy basis
(M_b^4,M_b^2 R_L,C^2,R_L^2,E4,box_L R_L), with S=-int F dt.
Its already verified conventions give R_E=-R_L; the squared invariants
and static box_E R_E=box_L R_L agree. The coefficients to feed that owner
are therefore

\[
(V_D/M_b^4,\ A_D/M_b^2,\ C_{W,D},\ 0,\ C_{E,D},\ C_{\Box,D}).
\]

M_b just normalizes the basis. It is distinct from nu and from a determinant
normalization mass. The sign/units conversion is checked without recomputing
the metric basis. The full invariant, including its required boundary
variation, must be counted once when moving the induced Einstein term to
the geometric side of the source equations.

The known light-field in-in functional carries the actual state, occupations
and correlations. The Euclidean identity above does not by itself produce
the heavy sector's causal response, cross-correlations or initial state.
The remaining assembly is

\[
\Gamma_{\rm one}^{\rm CTP}
 =\Gamma_{\rm light}^{\rm CTP}
  +\Gamma_H^{\rm CTP}+\Gamma_{\rm completion}^{\rm CTP}.
\]

Gamma_completion denotes the contribution still required by the common
measure, compensator, transmitting interaction and other fields; it is
not an independently adjustable gravitational action. Its real finite
terms are not fixed by Z[g,g]=1, as established by the earlier normalized
state calculation. The exact spectral complement is now available to
match that completion. Repeating raw Dirac coefficient scans will not
determine it. No self-sourced geometry or cosmological energy-transfer
prediction follows from this coefficient map alone.

## Reproduction

`scripts/check_nsc_compact_matching.py --check` authenticates its dependencies
and compares every field of `results/development/compact-matching.json`.
The existing comparison policy uses absolute tolerance 3e-9 and relative
tolerance 3e-8, with exact keys, structure, non-float types and hashes.
New output paths refuse overwrite. Only the new matching calculation runs;
prior scientific generators remain stored, authenticated inputs.
