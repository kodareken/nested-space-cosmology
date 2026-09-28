# FGC-1-HYP1-DOM1-QIFT1: quantified full-dimensional implicit branch box

**FGC-1-HYP1-DOM1-QIFT1** closes the first quantified-domain obligation left
by [FO1-RC1](fgc-hyp1-fo1-rc1.md). It proves, using exact rational closed
intervals, that the complete nonlinear REF1 residual has one and only one
six-component coordinate-time acceleration in a declared acceleration box for
every point of a nonzero-width, thirty-dimensional mixed/spatial two-jet
parameter box around the IMP1 flat FGC-QR root.

This is a real continuous-box theorem. It is not a grid, random sample,
floating-point observation, closed-form nonlinear solution, hyperbolicity
theorem, constraint theorem, EFT-validity claim, or authorization to evolve.

## 1. Frozen local relation

The six base fields remain

\[
u=(h_{tt},h_{tr},h_{rr},R,\phi,\chi).
\]

FO1-RC1 defines

\[
p=\partial_tu,
\qquad
q=\partial_ru,
\]

and rewrites the complete REF1 physical rows as

\[
\mathcal R(a;z)=0,
\qquad
a=p_t,
\qquad
z=(u,p,q,p_r,q_r).
\]

The order of the thirty independent parameter coordinates is frozen as

\[
(u^I,p^I,q^I,p_r^I,q_r^I),
\qquad I=1,\ldots,6,
\]

in the existing base-field order. The six acceleration coordinates are
\(a^I=p_t^I\). The coordinate radius is fixed at \(r=4\); it is a coordinate
parameter of this local chart, not one of the thirty state-box axes.

The exact center is the IMP1 root

\[
u_0=(-1,0,1,4,0,0),
\qquad
p_0=0,
\qquad
q_0=(0,0,0,1,0,0),
\]

\[
(p_r)_0=(q_r)_0=0,
\qquad
a_0=0.
\]

The nonzero action parameters are retained:

\[
M_{\rm Pl}=2,
\quad \mu=3,
\quad g_4=\frac12,
\quad \beta=-\frac14,
\quad \eta=\frac12,
\quad \alpha=0.
\]

At the center, the exact complete REF1 residual vanishes and

\[
J_0=D_a\mathcal R(a_0;z_0),
\qquad
\det J_0=-165888.
\]

## 2. Declared boxes

Every one of the thirty parameter axes varies independently over the same
closed rational interval:

\[
Z=
\left\{
z:\|z-z_0\|_\infty\le \rho_z
\right\},
\qquad
\rho_z=\frac1{65536}.
\]

The acceleration enclosure is

\[
A=
\left\{
a:\|a-a_0\|_\infty\le \rho_a
\right\},
\qquad
\rho_a=\frac1{128}.
\]

Thus \(Z\) is a compact box with nonempty open interior in the complete
thirty-dimensional local jet space used by the implicit residual. It is not a
one-dimensional path or a lower-dimensional parameter slice.

These coordinate widths are a mathematical local-domain certificate in the
frozen component norm. They are not claimed to be invariant physical length,
energy, or curvature scales.

## 3. Exact interval lift

The implementation introduces closed intervals

\[
[x]=[\underline x,\overline x],
\qquad
\underline x,\overline x\in\mathbb Q,
\]

and evaluates every arithmetic operation with exact `Fraction` endpoints.
There is no conversion to binary floating point and therefore no hidden
rounding direction to assume. Addition, subtraction, multiplication,
reciprocal, division, and integer powers return closed interval enclosures.
A reciprocal or division fails whenever its denominator interval contains
zero.

The interval type has no implicit truth value and no ordering operators. Every
sign-dependent branch in the inherited REF1 path must establish its sign over
the whole interval. An ambiguous metric signature, time covector, effective
Planck coefficient, auxiliary determinant, inverse pivot, or reciprocal fails
closed.

The complete existing REF1 tensor evaluator is used directly. QIFT1 does not
maintain a second curvature or field-equation implementation. A separate
interval first-tangent type propagates one exact interval directional
derivative. Six seeded evaluations enclose every entry of

\[
D_a\mathcal R(A,Z)\subseteq[J_a].
\]

The exact IMP1 center matrix is checked entrywise to lie inside this interval
Jacobian.

## 4. Uniform contraction certificate

Let

\[
B=J_0^{-1}
\]

be the exact rational inverse already fixed by IMP1, and define, for each
\(z\in Z\),

\[
T_z(a)=a-B\mathcal R(a;z).
\]

QIFT1 computes the interval matrix

\[
[E]=I-B[J_a]
\]

and the componentwise enclosure

\[
[K]
=
-B[\mathcal R(a_0;Z)] + [E]\,(A-a_0).
\]

The executable gate verifies both strict inequalities

\[
\|[E]\|_\infty<1,
\tag{QIFT-1}
\]

and

\[
[K]\subset\operatorname{int}(A).
\tag{QIFT-2}
\]

Equation (QIFT-1) makes every \(T_z\) a contraction on \(A\), uniformly in
\(z\). Equation (QIFT-2) makes \(A\) invariant under every such map. The
Banach fixed-point theorem then gives

\[
\boxed{
\forall z\in Z,
\quad
\exists!\,a(z)\in A
\quad\text{such that}\quad
\mathcal R(a(z);z)=0.
}
\]

It also proves convergence of the fixed-preconditioner iteration from every
starting point in \(A\), and the uniform inverse estimate

\[
\|(D_a\mathcal R)^{-1}\|_\infty
\le
\frac{\|B\|_\infty}{1-\|[E]\|_\infty}.
\]

The machine-readable result stores the complete exact interval residual,
interval Jacobian, preconditioned matrix, contraction bound, image box, and
six strict inclusion margins. For readable scale, its exact values satisfy

\[
\|[E]\|_\infty<\frac1{250},
\qquad
[K]\subset
\left[-\frac1{512},\frac1{512}\right]^6
\subset
\left(-\frac1{128},\frac1{128}\right)^6.
\]

These simpler rational comparisons are summaries of the stored exact values,
not replacements for them.

## 5. Regular-domain ledger

The same interval path proves strict margins throughout \(Z\) for:

- positive areal radius;
- negative two-dimensional base-metric determinant;
- a timelike coordinate-time covector, \(g^{tt}<0\);
- positive effective Planck coefficient
  \(F=M_{\rm Pl}^2+\beta\phi^2\);
- Lorentzian \(\widetilde g\) and \(\widehat g\) auxiliary inverse metrics;
- positive separation from the inner edge of the flat spherical reference
  annulus.

The exact lower margins are serialized. They obey the simpler strict bounds

\[
R>3,
\qquad
-\det h>\frac{99}{100},
\qquad
-g^{tt}>\frac{99}{100},
\qquad
F>3,
\]

\[
-\det\widetilde g^{ab}>\frac1{65},
\qquad
-\det\widehat g^{ab}>\frac1{29},
\qquad
r-r_{\min}=\frac72.
\]

These inequalities certify this local algebraic chart. They are not a
uniform physical-characteristic or symmetrizer theorem.

## 6. What this closes

IMP1 used the ordinary implicit-function theorem at one exact root, so its
neighborhood existed but had no published radius. FO1-RC1 then exposed the
complete local first-order implicit relation but did not solve it away from
that root.

QIFT1 now proves all of the following for the frozen REF1 formulation:

- the parameter domain has thirty independently varying, nonzero-width axes;
- the complete nonlinear residual is enclosed on that continuous domain;
- the complete nonlinear acceleration Jacobian is enclosed on
  \(A\times Z\);
- every declared denominator and sign branch stays regular;
- one exact fixed preconditioner gives a contraction uniformly on the box;
- every parameter point has exactly one acceleration root in \(A\);
- the acceleration Jacobian has a uniform inverse-norm bound on the box.

The theorem encloses the branch \(a(z)\); it does not claim an elementary
closed-form expression for it.

## 7. Why hyperbolicity remains separate

The flat root is an excellent center for the implicit solve, but it is not a
free open-domain hyperbolicity proof. At the inactive point

\[
\phi=0,
\qquad
\nabla\phi=0,
\qquad
\nabla\nabla\phi=0,
\]

the regulator quadratic coincides with the physical metric quadratic. MHG1
proved that this repeated root is semisimple at the exact point. Continuity of
a repeated real eigenvalue does not by itself prove that every nearby matrix
has a real, uniformly conditioned eigenbasis.

The next characteristic gate must therefore use the QIFT1 solved-branch
enclosure and either:

1. construct an explicit uniformly bounded eigenbasis on a nonflat subbox
   where the physical and regulator cones have a strict separation margin; or
2. derive an explicit positive full-system symmetrizer with uniform upper and
   lower bounds.

Pointwise discriminants or a finite fixture list are not substitutes for that
gate.

The downstream [FGC-1-HYP1-CON1-COMP1](fgc-hyp1-con1-comp1.md) artifact now
chooses one activated point strictly inside this box, with nonzero `phi` and
`partial_r phi`. It proves exact local physical/gauge/reduction compatibility
at QIFT1's unique root and a distinct regulator cone at that point. It does
not turn QIFT1 into a compatible initial hypersurface, propagation theorem, or
uniform hyperbolicity box.

## 8. Why “retained EFT” remains false

The present action has fixed numerical couplings, but the repository has not
yet declared an EFT cutoff, dimensional normalization, omitted-operator
basis, or bounds on their Wilson coefficients. A small local field/jet box
and a regular principal system cannot supply those missing physical inputs by
definition.

QIFT1 is therefore a **quantified regular local model branch box**. It is not
called an EFT-validity domain. A later validity ledger must state the cutoff,
normalizations, retained operators, omitted-operator assumptions, and the
dimensionless margins that make the truncation meaningful—or explicitly
treat ACT1 only as a classical candidate model rather than an EFT truncation.

## 9. Fail-closed nonclaims

QIFT1 does **not** establish:

- an exact closed-form nonlinear acceleration map;
- a complete quasilinear first-order evolution system;
- metric-derived gauge-constraint propagation;
- a complete reduction-constraint subsystem;
- Hamiltonian, momentum, or other physical constraints and their propagation;
- constraint-compatible initial or boundary data;
- a uniform radial strong-hyperbolicity box or full-system symmetrizer;
- an EFT cutoff, omitted-operator envelope, or retained-EFT domain;
- a regular center, horizon boundary system, or well-posed IBVP;
- an evolution, independent-`chi` collapse, restoring interval,
  metric-null affine defocusing, singularity resolution, or child spacetime.

## 10. Reproduction

Run

```bash
python3 scripts/reproduce_fgc_hyp1_dom1_qift1.py \
  --config configs/fgc/fgc-1-hyp1-dom1-qift1.toml \
  --output results/fgc-1-hyp1-dom1-qift1.json
```

The canonical record binds the configuration, complete predecessor config
chain, exact interval implementation, inherited REF1 implementation, this
derivation document, full interval certificate, and retained nonclaims by
SHA-256.
