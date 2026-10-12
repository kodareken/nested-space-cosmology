# FGC-1-CON4-PHY1: physical, gauge, and reduction constraint closure

## Result

**FGC-1-CON4-PHY1** closes the physical, metric-defined gauge, and
one-dimensional kinematic-reduction constraint system for the declared
spherical REF1 formulation, conditional on a sufficiently smooth solution and
restricted to its boundary-free domain of dependence.

The result is structural rather than empirical. It does not infer constraint
propagation from a small residual in one numerical run. It derives:

1. the Hamiltonian and radial-momentum projections of the unredefined
   ACT1/VAR1 equation that must be imposed as
   physical initial constraints;
2. why those projections contain no coordinate-time accelerations;
3. an exact invertible map between those physical constraints and the normal
   derivative of the metric-defined REF1 gauge vector;
4. the composition with the already certified homogeneous gauge subsidiary
   equation and conditional Cauchy uniqueness theorem;
5. the complete six-field FO1 kinematic-reduction subsidiary identity; and
6. executable dimensionless cancellation monitors for every constraint
   family.

This advances **FGC-1-RUN1-SYM1** by one predicate. It does not supply a
regular centre, a compatible nonzero-width initial-data family, an incoming
constraint-preserving boundary map, a quasilinear existence theorem, a
collapse trajectory, or retained-EFT authorization.

## 1. Physical constraints are taken from the unredefined equations

Write the spherical base metric as

\[
g_{AB}=
\begin{pmatrix}
A & B\\
B & C
\end{pmatrix},
\qquad C>0,
\qquad AC-B^2<0.
\]

Define

\[
s=\frac{B}{C},
\qquad
\ell^2=-A+\frac{B^2}{C}>0,
\qquad
w^A=(1,-s).
\]

Here \(w^A\) is the unnormalized normal to a coordinate-time slice and
\(\ell\) is its lapse. For the unredefined ACT1/VAR1 metric residual
\(E_{ab}\), CON4 uses

\[
\mathcal H = E_{ab}w^aw^b
              =E_{tt}-2sE_{tr}+s^2E_{rr},
\]

\[
\mathcal M = E_{ab}w^a(\partial_r)^b
              =E_{tr}-sE_{rr}.
\]

Their unit-normal versions differ only by strictly positive factors:

\[
\mathcal H_{\rm unit}=\frac{\mathcal H}{\ell^2},
\qquad
\mathcal M_{\rm unit}=\frac{\mathcal M}{\ell\sqrt C}.
\]

The zero sets are therefore identical. No REF1 gauge-extension term is
included in these definitions.

The exact evaluator decomposes each projection into the signed contributions

\[
F G_{ab},\qquad
(g_{ab}\Box-\nabla_a\nabla_b)F,\qquad
-T^{(\phi)}_{ab},\qquad
-T^{(\chi)}_{ab},\qquad
8P_{acbd}\nabla^{cd}f.
\]

The machine certificate verifies that the five projected terms reassemble
\(\mathcal H\) and \(\mathcal M\) exactly at every control state. This is the
same unredefined tensor evaluator used by RED1 and VAR1.

## 2. Why these are constraints rather than acceleration equations

A legitimate initial-data equation cannot depend on the six coordinate-time
accelerations that the REF1 evolution rows solve for. The absence is derived
sector by sector.

### Einstein sector

The normal projections are the Gauss--Codazzi identities,

\[
G_{nn}=\frac12\left({}^{(3)}R+K^2-K_{ij}K^{ij}\right),
\]

\[
G_{ni}=D_jK^j{}_i-D_iK.
\]

They contain intrinsic slice curvature, extrinsic curvature, and its spatial
derivatives, but no second normal derivative of the metric.

### Nonminimal \(F\) sector

The double-normal scalar Hessian cancels in the Hamiltonian projection:

\[
n^an^b(g_{ab}\Box F-\nabla_a\nabla_bF)
=-h^{ij}\nabla_i\nabla_jF.
\]

The momentum projection is mixed normal--spatial,

\[
n^ah_i{}^b(g_{ab}\Box F-\nabla_a\nabla_bF)
=-n^ah_i{}^b\nabla_a\nabla_bF,
\]

so it also contains no double-normal scalar derivative.

### Matter stresses

Both scalar stress tensors depend on field values and first derivatives only.
They contain no coordinate-time accelerations.

### Scalar--Gauss--Bonnet response

The important identity is the double-dual representation

\[
P_{acbd}\propto
\epsilon_{acEF}\epsilon_{bdGH}R^{EFGH}.
\]

For the Hamiltonian projection, both exposed indices are normal. Nonzero
Levi-Civita support forces the remaining curvature indices and the two Hessian
indices to be spatial. Thus neither factor contains a double-normal
derivative.

For the radial-momentum projection, one exposed index is normal and one is
spatial. Nonzero support permits at most one normal index in either the
curvature or scalar-Hessian factor. Those are Codazzi or mixed derivatives,
not coordinate-time accelerations.

The certificate enumerates the finite four-dimensional Levi-Civita support:
36 nonzero abstract terms for the Hamiltonian pattern and 36 for one momentum
direction. Every Hamiltonian curvature/Hessian index is spatial; every
momentum term contains at most one normal index.

As an independent implementation regression, exact first-tangent automatic
differentiation is then propagated through the complete unredefined tensor
evaluator at three distinct rational FGC-QR two-jets, including two shifted,
nonflat states. All 36 derivatives—six acceleration directions times two
constraints times three states—vanish exactly. Those samples test the code;
the covariant decomposition above is the general proof.

## 3. Exact physical-to-gauge map on the REF1 shell

REF1 evolves

\[
E_{ab}+X_{ab}=0,
\]

where

\[
X^{\mu\nu}
=F\,\hat P_{\alpha}{}^{\beta\mu\nu}
\nabla_\beta C^\alpha
\]

and \(C^\alpha[g]\) is the metric-defined modified-harmonic gauge vector.

If \(C^\alpha=0\) as a function on the initial slice, all spatial covariant
derivatives of \(C^\alpha\) vanish there. In spherical symmetry the active
components are \(C^t,C^r\). Put

\[
x=\nabla_tC^t,
\qquad
y=\nabla_tC^r,
\qquad
q>1
\]

with \(q\) the declared hat-cone normal factor. Direct projection of
\(E_{ab}=-X_{ab}\) gives

\[
\begin{pmatrix}
\mathcal H\\
\mathcal M
\end{pmatrix}
=
\frac{Fq}{2}
\begin{pmatrix}
\ell^2 & 0\\
-B & -C
\end{pmatrix}
\begin{pmatrix}
x\\
y
\end{pmatrix}.
\]

Its determinant is

\[
\det=-\frac{F^2q^2\ell^2C}{4}<0
\]

whenever

\[
F>0,\qquad q>1,\qquad \ell^2>0,\qquad C>0.
\]

Therefore, on the complete REF1 metric shell and on a slice where \(C=0\),

\[
\mathcal H=\mathcal M=0
\quad\Longleftrightarrow\quad
n^a\nabla_aC^\mu=0.
\]

The implementation obtains the two matrix columns by inserting exact basis
values for \(\nabla_tC^t\) and \(\nabla_tC^r\) into the existing REF1 extension
evaluator. It independently constructs the analytic matrix above and requires
entry-by-entry equality. The equality and strict determinant sign pass at all
three exact controls, including nonzero shifts where the lower-left matrix
entry is nonzero.

This upgrades COMP1's one-point compatibility observation to a general
spherical algebraic map on the stated open metric domain. It still does not
construct a compatible hypersurface.

## 4. Conditional propagation chain

The result composes two already sealed predecessors rather than replacing
them.

**CON2-MPROP1** derives the metric-defined identity

\[
L_{\hat g}(C[g])
=\nabla_aE_{\rm REF1}^{ab}
 +E_\phi\nabla^b\phi
 +E_\chi\nabla^b\chi.
\]

On the complete REF1 metric and both unmodified scalar equations this becomes
the homogeneous, normally hyperbolic system

\[
L_{\hat g}(C[g])=0.
\]

**CON3-CAU1** records the conditional boundary-free uniqueness consequence:
if \(F>0\), \(\hat g\) is Lorentzian, the initial slice is
\(\hat g\)-spacelike, and both \(C\) and its normal derivative vanish, then
\(C=0\) throughout the boundary-free domain of dependence of any sufficiently
smooth full solution.

CON4 supplies the missing physical handoff:

```text
C = 0 on the initial slice
        -> spatial nabla C = 0
H = M = 0 and the REF1 metric equations
        -> normal nabla C = 0 by the invertible 2x2 map
CON2 homogeneous subsidiary equation + CON3 uniqueness
        -> C = 0 before a boundary enters
C = 0
        -> X = 0
        -> REF1 recovers the unredefined ACT1/VAR1 metric equation
        -> H = M = 0 remains true
```

This is a conditional physical-constraint propagation theorem. It is not a
solution-existence theorem and not an initial-boundary value theorem.

## 5. Kinematic reduction constraints

For each of the six FO1 fields, define

\[
D=\partial_tu-p,
\qquad
C_R=q-\partial_ru,
\qquad
K=\partial_tq-\partial_rp.
\]

Commutation of differentiable \(t,r\) partials gives the exact off-shell
identity

\[
\partial_tC_R+\partial_rD-K=0.
\]

Thus \(D=K=0\) implies \(\partial_tC_R=0\). There is no independent spatial
curl constraint in one dimension. The six reduction subsidiary fields have
zero radial principal speed and supply no incoming boundary characteristic.

The certificate exercises both an on-shell control with \(D=C_R=K=0\) and an
off-shell control with all three residual families nonzero. The subsidiary
identity vanishes exactly in both cases.

## 6. Constraint characteristics

The active spherical gauge subsidiary components share the hat-metric
principal polynomial

\[
\hat g^{ab}\xi_a\xi_b=0,
\qquad
\xi_a=(-c,1).
\]

At each exact control the two roots are real, nonzero, and have signs
\((-1,+1)\). Each of \(C^t,C^r\) therefore contributes one incoming hat-cone
branch at each end of the annular control domain. The physical
\(\mathcal H,\mathcal M\) pair is algebraic normal gauge data on the REF1 shell,
not a separate propagating characteristic family. The six reduction
constraints have zero speed.

CON4 classifies these families but intentionally does not translate the
incoming gauge fields into a stable nonlinear boundary operator for the main
variables. That remains the distinct **FGC-1-BND2-CP1** obligation.

## 7. Canonical monitors

For every residual represented as signed terms \(r=\sum_i r_i\), CON4 defines

\[
\mathcal N_{\rm rel}(r)
=\frac{|\sum_i r_i|}{\sum_i|r_i|}.
\]

The monitor is dimensionless, invariant under a common rescaling of all terms,
and bounded by one. Every record retains:

- the raw residual;
- every signed contribution;
- the absolute-term denominator;
- the normalized value; and
- whether the exact zero-scale convention was used.

If every term is exactly zero, the denominator vanishes and the monitor is
defined to be zero with an explicit flag. This avoids an arbitrary
dimensionful numerical floor.

The same executable contract is used for:

- the five physical contributions to \(\mathcal H\) and \(\mathcal M\);
- every contraction term in \(C^\mu\);
- every partial and connection contribution to \(\nabla_aC^\mu\); and
- every \(D,C_R,K\) and reduction-subsidiary identity.

This relative cancellation is not an absolute truncation-error norm. HLT1 now
adds physical scales and strict transactional thresholds; NUM1 must still add
discretization estimates and validate the integrated numerical stop decision.
Smallness in one run is never used as proof of propagation.

## 8. What is now closed

The formulation now distinguishes, in executable form:

| Role | Equations or fields |
|---|---|
| Solved by future evolution | Six complete REF1 metric/scalar rows plus the FO1 kinematic equations |
| Imposed on initial data | \(C^t=C^r=0\), unredefined \(\mathcal H=\mathcal M=0\), and six radial reduction constraints |
| Propagated analytically | Gauge constraints by CON2/CON3, physical constraints by the invertible REF1-shell map, and reduction constraints by the FO1 identity |
| Monitored numerically | Raw and normalized physical, gauge, gauge-derivative, and reduction residuals |

The gate

```text
physical_gauge_reduction_constraint_system_closed = true
```

means exactly this conditional, boundary-free spherical closure. It does not
authorize a run by itself.

## 9. Machine reproduction

```bash
python3 scripts/reproduce_fgc_con4_phy1.py \
  --output results/fgc-1-con4-phy1.json
```

The canonical result is bound to:

- ACT1 and VAR1;
- COMP1's exact compatible local datum;
- CON2's metric-derived gauge and reduction subsidiary identities;
- CON3's conditional boundary-free gauge theorem;
- the frozen SF1 protocol and RUN1 scope; and
- the exact implementation and this derivation document.

No FGC-QR collapse holdout is read or executed.

## Nonclaims

CON4-PHY1 does not derive or authorize:

- a regular-centre formulation;
- a compatible finite-mass, nonzero-width initial-data family;
- a constraint-preserving outer-boundary map or full IBVP;
- existence of a smooth nonlinear REF1 evolution;
- a classical collapse calculation;
- retained-Wilsonian-EFT validity;
- metric-null defocusing or singularity resolution;
- a child domain or topology change;
- a dark-sector mechanism; or
- a varying locally measured speed of light.
