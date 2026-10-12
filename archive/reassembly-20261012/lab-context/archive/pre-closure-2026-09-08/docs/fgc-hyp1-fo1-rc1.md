# FGC-1-HYP1-FO1-RC1: exact local first-order lift and radial reduction identity

**FGC-1-HYP1-FO1-RC1** closes a narrow formulation prerequisite left open by
[IMP1](fgc-hyp1-mhg-implicit.md) and
[PROP1](fgc-hyp1-mhg-propagation.md). It rewrites the complete six-equation
REF1 spherical two-jet residual as an exact local first-order implicit relation
for eighteen variables, differentiates the complete implicit branch at the
IMP1 flat root, and derives the kinematic radial reduction-constraint identity
for the declared coordinate-derivative reduction.

The result is deliberately bounded. By itself it is not a nonlinear acceleration
solver, a quantified implicit-function neighborhood, a complete
gauge/physical-constraint formulation, an initial-boundary value problem, an
open-domain hyperbolicity theorem, or an evolution authorization.

The machine-readable authority is
[`results/fgc-1-hyp1-fo1-rc1.json`](../results/fgc-1-hyp1-fo1-rc1.json),
generated from
[`configs/fgc/fgc-1-hyp1-fo1-rc1.toml`](../configs/fgc/fgc-1-hyp1-fo1-rc1.toml)
by
[`scripts/reproduce_fgc_hyp1_fo1_rc1.py`](../scripts/reproduce_fgc_hyp1_fo1_rc1.py).

## 1. Why this gate exists

RED1 evaluates the unredefined ACT1/VAR1 equations on exact spherical
two-jets. MHG1 supplies the pointwise gauge-fixed second-order principal
symbol and its standard twelve-variable principal reduction. REF1 adds the
complete lower-order reference-connection gauge residual. IMP1 then proves
that the complete REF1 residual can be solved locally for the six
coordinate-time accelerations at one exact flat FGC-QR root. PROP1 derives a
conditional lower-order operator on an independently supplied gauge-vector
two-jet.

None of those artifacts had yet written the complete REF1 residual in a
declared first-order variable set. In particular, MHG1's twelve-by-twelve
matrix is a **principal** reduction only. It omits the base fields and the
complete lower-order residual, so it cannot be relabeled as a complete
first-order evolution system.

FO1-RC1 performs the smallest exact lift that current evidence supports. It
does not alter the action, redefine the physical equations, or substitute a
principal approximation for REF1.

## 2. Frozen fields and first-order state

The base fields remain exactly the RED1/REF1 variables

\[
u^I=(h_{tt},h_{tr},h_{rr},R,\phi,\chi),
\qquad I=1,\ldots,6.
\]

Define coordinate-derivative variables

\[
p^I:=\partial_t u^I,
\qquad
q^I:=\partial_r u^I,
\]

and freeze the eighteen-component state

\[
\mathsf U=(u^I,p^I,q^I).
\]

This is a base-metric variable set. It is **not** silently relabeled as
generalized ADM variables \((\alpha,s,\lambda,R,\phi,\chi)\). Such a change
would require the full nonlinear field and jet chain rule already discussed
in the IMP1 owner document.

At one local first jet of \(\mathsf U\), FO1-RC1 stores

\[
(u,p,q;u_t,p_t,q_t;u_r,p_r,q_r).
\]

Keeping \(u_t\) and \(u_r\) distinct from \(p\) and \(q\) makes the
definition and reduction constraints observable rather than assumed.

## 3. Exact two-jet recomposition

The complete REF1 evaluator requires the local two-jet

\[
(u,u_t,u_r,u_{tt},u_{tr},u_{rr}).
\]

On the declared first-order variables the exact substitution is

\[
u_t=p,
\qquad
u_r=q,
\qquad
u_{tt}=p_t,
\qquad
u_{tr}=p_r=q_t,
\qquad
u_{rr}=q_r.
\]

No third jet is needed to write or evaluate this relation. The implementation
performs a lossless round trip between every ordinary-rational `SphericalState`
two-jet and its first-order representation. Both the flat IMP1 root and the
activated RED1 control return entry by entry to their original two-jets.

## 4. The exact eighteen-row implicit relation

With REF1 equation order

\[
\mathcal R^A_{\rm REF1}=
(E^{\rm MHG}_{tt},E^{\rm MHG}_{tr},E^{\rm MHG}_{rr},
E^{\rm MHG}_{\theta\theta},E_\phi,E_\chi),
\]

FO1-RC1 freezes the rows

\[
\boxed{
\begin{aligned}
\mathcal D_u^I
  &:=\partial_tu^I-p^I=0,\\
\mathcal D_R^A
  &:=\mathcal R^A_{\rm REF1}
    (u,p,q,\partial_t p,\partial_r p,\partial_r q;r)=0,\\
\mathcal D_q^I
  &:=\partial_tq^I-\partial_rp^I=0.
\end{aligned}}
\]

The six physical rows are the **complete** unredefined ACT1/VAR1 plus REF1
reference-gauge residual. FO1-RC1 delegates to the existing tensor evaluator;
it does not reconstruct those rows from a truncated principal matrix.

The word “implicit” matters. Gauss--Bonnet curvature--Hessian products make the
residual generally nonlinear in second jets. Consequently

\[
\mathcal R_{\rm REF1}
(u,p,q,p_t,p_r,q_r;r)=0
\]

is a first-order fully nonlinear differential relation. It cannot yet be
advertised globally as

\[
\partial_t\mathsf U+A(\mathsf U)\partial_r\mathsf U=S(\mathsf U)
\]

or used as a numerical source evaluator.

## 5. Kinematic radial reduction constraint

For this exact coordinate-derivative reduction define

\[
\mathcal C_r^I:=q^I-\partial_ru^I.
\]

Assuming sufficiently smooth fields and commuting coordinate partials,

\[
\begin{aligned}
\partial_t\mathcal C_r^I
&=\partial_tq^I-\partial_r\partial_tu^I\\
&=\partial_tq^I-\partial_rp^I\\
&=\mathcal D_q^I.
\end{aligned}
\]

Therefore the declared rows \(\mathcal D_u=0\) and \(\mathcal D_q=0\) imply

\[
\boxed{\partial_t\mathcal C_r^I=0.}
\]

This is an exact **kinematic auxiliary-constraint identity for this chosen
reduction**. It is not a universal property of every first-order reduction.
Advection terms or additions proportional to reduction constraints generally
produce a different homogeneous constraint subsystem.

It also does not imply propagation of any of the following:

- the metric-defined REF1 gauge vector \(C^\mu[g]\);
- Hamiltonian, momentum, or other physical initial constraints;
- mixed auxiliary constraints introduced by a different reduction;
- regular-center or outer-boundary compatibility conditions.

The standard distinction between derivative-reduction constraints and the
physical/gauge constraints is explicit in Gundlach and Martín-García's
parameterized first-order reductions and their auxiliary-constraint subsystem
(arXiv:[gr-qc/0506037](https://arxiv.org/abs/gr-qc/0506037), Sec. III,
Eqs. 15--21). FO1-RC1 uses that distinction as scope discipline; the elementary
identity above is derived directly for this project rather than attributed to
the paper.

## 6. Exact complete residual Jacobian

At fixed reference radius and parameters, freeze the physical argument order

\[
z_{\rm all}=(u,p,q,p_t,p_r,q_r).
\]

The implementation seeds every one of the 36 slots with the exact tangent
number \(x+\epsilon\dot x\), \(\epsilon^2=0\), and evaluates the complete REF1
code path unchanged. It publishes the six blocks

\[
J_u,
J_p,
J_q,
J_{p_t},
J_{p_r},
J_{q_r}.
\]

Every seeded primal reproduces the unseeded complete residual. A second pass
seeds all 36 directions simultaneously with distinct rational weights and
agrees exactly with the corresponding weighted sum of stored columns. This is
a column-assembly and linear-superposition check, not a finite-difference
approximation.

FO1-RC1 extends the exact tangent path through field values and first
derivatives. The auxiliary inverse-metric constructor preserves tangents while
performing its Lorentzian-domain check on the exact primal matrix. Ordinary
rational behavior remains unchanged.

## 7. What IMP1 gives at the flat root

Write

\[
a:=p_t,
\qquad
z:=(u,p,q,p_r,q_r),
\qquad
\mathcal R(a,z)=0.
\]

IMP1 established at its exact flat FGC-QR reference-gauge root

\[
\mathcal R(a_0,z_0)=0,
\qquad
J_a:=\left.\frac{\partial\mathcal R}{\partial a}\right|_0,
\qquad
\det J_a=-165888\ne0.
\]

The finite-dimensional implicit function theorem gives an unquantified smooth
local map

\[
a=\mathcal A(z).
\]

FO1-RC1 does not pretend to evaluate \(\mathcal A\) nonlinearly. It computes
its **complete derivative at the root**:

\[
\boxed{
D\mathcal A|_0=-J_a^{-1}D_z\mathcal R|_0.
}
\]

The generated certificate publishes all \(6\times30\) rational entries and
checks the exact implicit-derivative identity

\[
J_aD\mathcal A|_0+D_z\mathcal R|_0=0.
\]

This is the complete linearized local source/derivative map, not itself a
nonlinear source function or quantified open neighborhood. The downstream
[FGC-1-HYP1-DOM1-QIFT1](fgc-hyp1-dom1-qift1.md) certificate now quantifies the
frozen REF1 implicit branch on one full-dimensional rational box.

## 8. Linearized first-order matrices

Let

\[
D\mathcal A|_0=
(A_u,A_p,A_q,A_{p_r},A_{q_r}).
\]

The exact linearized system can be written

\[
\partial_t\delta\mathsf U+
\mathsf A_0\partial_r\delta\mathsf U
=\mathsf S_0\delta\mathsf U,
\]

with frozen state order \((u,p,q)\),

\[
\mathsf A_0=
\begin{pmatrix}
0&0&0\\
0&-A_{p_r}&-A_{q_r}\\
0&-I&0
\end{pmatrix},
\qquad
\mathsf S_0=
\begin{pmatrix}
0&I&0\\
A_u&A_p&A_q\\
0&0&0
\end{pmatrix}.
\]

Both eighteen-by-eighteen matrices are serialized exactly. The lower
twelve-by-twelve \((p,q)\) radial principal block agrees entry by entry with
MHG1's independently assembled standard first-order principal matrix at the
same flat root. This regression fixes the signs and the `dtr` ordering.

Kovács and Reall's Appendix A gives the relevant generic distinction between a
first-order system, a solved time form, and a smooth uniformly bounded
symmetrizer (arXiv:[2003.08398](https://arxiv.org/abs/2003.08398),
Eqs. 119--132). Their result does not supply the ACT1/REF1 lower-order matrices
or prove this project's open-domain hyperbolicity.

## 9. Activated off-shell control

The second control is RED1's exact `FGCQR_activated_generic` fixture at
\(r=9/2\). It has nonzero regulator value, gradient, and Hessian and exercises
the coupled ACT1/REF1 code away from flat space.

FO1-RC1 verifies:

- exact first-order-to-two-jet reconstruction;
- exact equality with the direct complete REF1 residual;
- zero declared kinematic/reduction rows for the packed genuine two-jet;
- a nonzero complete residual, proving that the control is off shell and
  nontrivial;
- all 36 exact residual-derivative columns and their weighted-superposition
  check;
- the observed pointwise acceleration-block determinant and rank.

The fixture is explicitly serialized as
`off_shell_control_not_a_solution`. A nonsingular derivative block at an
off-shell point does not create an IFT solution branch there, because the
required residual root has not been found.

## 10. Exact pass conditions

The generated gate requires all of the following:

1. the frozen state and equation orders each contain exactly eighteen entries;
2. both control two-jets round-trip exactly through the first-order variables;
3. both complete REF1 residuals reproduce entry by entry;
4. all declared kinematic rows and reduction constraints vanish on the packed
   controls;
5. the kinematic radial reduction-constraint time derivative vanishes on the
   declared rows;
6. all 36 tangent primals reproduce the direct residual on both controls;
7. both weighted tangent passes equal the stored-column superpositions;
8. the flat acceleration block equals IMP1, has rank six and determinant
   `-165888`;
9. the complete root branch derivative satisfies the implicit derivative
   identity;
10. the flat \((p,q)\) radial principal matrix equals MHG1 entry by entry;
11. the activated control is explicitly off shell and has a nonzero residual;
12. no nonlinear acceleration solve is performed or claimed.

The configuration parser also rejects unknown keys, noncanonical rational
encodings, repository traversal, predecessor drift, reference-annulus drift,
auxiliary-cone drift, reordered fields/equations, a promoted activated
solution status, or any promoted open gate.

## 11. What is now established

FO1-RC1 establishes, for the frozen ACT1/VAR1/REF1 formulation:

- an exact local eighteen-variable first-order implicit/DAE lift;
- exact reconstruction of the complete REF1 physical residual from those
  variables;
- the exact kinematic radial reduction-constraint identity for the declared
  coordinate derivative rows;
- the complete exact derivative of the IMP1 acceleration branch at its flat
  root;
- exact eighteen-by-eighteen linearized radial-principal and lower-order source
  matrices at that root;
- exact agreement of the linearized \((p,q)\) principal block with MHG1;
- a nontrivial activated off-shell coefficient control.

These are deterministic-computation (`D-C`) results for the declared fixtures
and conventions.

## 12. What remains open

FO1-RC1 does **not** establish:

- an explicit nonlinear map \(p_t=\mathcal A(u,p,q,p_r,q_r;r)\);
- a solved acceleration at the activated fixture or any nonflat solution;
- an FO1-RC1-native quantified neighborhood (the downstream QIFT1 record
  supplies a frozen REF1 branch box);
- affine dependence on radial derivatives or a complete quasilinear system;
- metric-derived REF1 gauge-constraint propagation;
- a complete reduction-constraint subsystem for a production formulation;
- Hamiltonian, momentum, or other physical constraints and their propagation;
- compatible initial data, regular-center data, outer-boundary conditions, or
  an initial-boundary value problem;
- a smooth uniform symmetrizer on an open retained-EFT domain;
- an evolution, independent-`chi` collapse, restoring interval, metric-null
  affine defocusing, singularity resolution, or child spacetime.

QIFT1 now closes the quantified regular implicit-branch-box threshold. The
downstream [FGC-1-HYP1-CON1-COMP1](fgc-hyp1-con1-comp1.md) result now closes
one activated metric-defined gauge/reduction/physical-compatible local root.
The next formulation threshold remains direct constraint propagation and a
compatible initial/boundary system, followed by a uniform nonflat
hyperbolicity and retained-validity certificate.

## 13. Reproduction

Run

```bash
python3 scripts/reproduce_fgc_hyp1_fo1_rc1.py \
  --output results/fgc-1-hyp1-fo1-rc1.json
python3 -m unittest -v \
  tests.test_fgc_modified_harmonic_first_order \
  tests.test_fgc_hyp1_fo1_rc1_reproduction
```

The reproducer source-binds ACT1, VAR1, RED1, SYM1, MHG1, REF1, IMP1, PROP1,
this configuration, the owner document, and the complete implementation path.
