# FGC literature ledger: gravity-transition wedge

## Purpose and current boundary

This is a **precedent and novelty ledger**, not an assertion that Finite
Gradient Closure (FGC) is realised in nature.  The approved first wedge is
narrow:

> In a covariant spherical model, can collapse of an independent matter field
> `chi` drive an invariant-defined threshold that activates a regulator field
> `phi`, while the physical metric has a healthy initial-value principal
> symbol and an explicitly computed affine Raychaudhuri margin?

The target is *not* a prescribed bounce metric, a generic claim that curvature
caps itself, an identification of `phi` with known particles, or a proof that
black holes produce child universes.  The matter field, regulator, physical
metric, characteristic cones, and initial data must be stated in the action.

For an affinely parametrised, hypersurface-orthogonal null congruence of the
**matter-coupled physical metric** `g_ab`, the required local diagnostic is

```text
d theta/d lambda = - theta^2/2 - sigma_ab sigma^ab - R_ab[g] k^a k^b.
```

Thus a proposed restoring interval must report the full evaluated quantity
`R_ab k^a k^b` (and its uncertainty/convergence), rather than calling a
positive potential, a regular metric ansatz, or a modified field equation
“defocusing.”  The statement is local; an integrated trapped-to-anti-trapped
transition and a global continuation are additional problems.  In a modified
gravity theory, metric null rays, scalar characteristics, and tensor
characteristics need not coincide.

## Direct scalarization / strong-field precedent

The nearest established family is scalar--Gauss--Bonnet (sGB) scalarization,
for example

```text
S = integral sqrt(-g) [M_Pl^2 R/2 - (nabla phi)^2/2 - V(phi)
                       + lambda^2 f(phi) G] d^4x,
G = R^2 - 4 R_ab R^ab + R_abcd R^abcd.
```

With `f'(0)=0`, the GR branch can remain present.  Linearisation about a
specified background can give an invariant-dependent effective mass,
`m_eff^2 = V''(0) - lambda^2 f''(0) G`, and hence a tachyonic bifurcation.
This establishes that curvature invariants can trigger nonlinear scalar
branches; it does **not** establish a finite-gradient regulator, regular
collapse core, or a negative physical-metric `R_ab k^a k^b`.

| Primary source | What is already covered | Remaining FGC gap |
| --- | --- | --- |
| [Silva et al. (2018), arXiv:1711.02080](https://arxiv.org/abs/1711.02080) | Curvature/GB-triggered spontaneous scalarization of black holes and compact stars. | The proposed wedge must use an **independent collapsing `chi` source**, a regulator `phi`, and calculate the collapse and Raychaudhuri diagnostics rather than only a vacuum scalarized branch. |
| [Silva et al. (2018), arXiv:1812.05590](https://arxiv.org/abs/1812.05590) | The nonlinear coupling controls quenching and radial stability; a quadratic choice alone need not be stable. | A regulator claim needs the analogous nonlinear saturation and perturbation calculation. |
| [Papallo (2017), arXiv:1710.10155](https://arxiv.org/abs/1710.10155) | Einstein--dilaton--GB fails strong hyperbolicity in any generalized-harmonic gauge of the class analyzed there on a generic weak-field background; more general Horndeski terms share the obstruction established there. | “Second order” and that generalized-harmonic class cannot be used as health certificates; a distinct explicit formulation and an action-specific principal analysis are required. |
| [Kovacs & Reall (2020), arXiv:2003.08398](https://arxiv.org/abs/2003.08398) | A modified-harmonic formulation separates physical, pure-gauge, and gauge-constraint cones and is strongly hyperbolic for the weakly coupled one-scalar Horndeski system studied there. | Use the theorem as a formulation architecture and weak-coupling comparator, not as proof that ACT1's activated regime or independent-`chi` extension lies in its open neighborhood; derive ACT1's dynamic-`F`, two-scalar system directly. |
| [East & Ripley (2021), arXiv:2011.03547](https://arxiv.org/abs/2011.03547) | Full nonlinear shift-symmetric ESGB evolutions implement the modified-harmonic architecture and remain hyperbolic in the exterior for couplings comparable to the range in which the spherical problem remains well posed. | Numerical success in a neighboring action does not transfer the retained health domain. ACT1 must close its own lower-order constraints, open-domain symmetrizer/EFT bounds, and independent-matter collapse. |
| [Ripley & Pretorius (2019), arXiv:1903.07543](https://arxiv.org/abs/1903.07543) | sGB collapse can evolve from hyperbolic to elliptic character. | Demonstrate that the proposed activation does not lose hyperbolicity before or during the claimed restoring interval. |
| [Hegade, Ripley & Yunes (2023), *Where and why does Einstein-Scalar-Gauss-Bonnet theory break down?*](https://arxiv.org/abs/2211.08477) | A gauge-covariant hyperbolicity diagnostic applied to collapse in two ESGB theories.  It separates null-convergence violation from loss of hyperbolicity and finds only a restricted predictive coupling region for the studied data. | Hyperbolicity and defocusing must be reported as two independent gates.  A negative `R_ab k^a k^b` is not evidence that the evolution remains predictive. |
| [Thaalba et al. (2023), *Spherical collapse in scalar-Gauss-Bonnet gravity: taming ill-posedness with a Ricci coupling*](https://arxiv.org/abs/2306.01695) | A quadratic Ricci coupling can mitigate hyperbolicity loss in spherical collapse of a scalar cloud and can yield scalarized black holes in a bounded mass range. | This is direct precedent for the lead action.  The proposed wedge must add genuinely independent collapsing matter `chi`, the full affine Raychaudhuri margin, and evolution that stops rather than passes beyond loss of predictivity. |
| [Thaalba et al. (2024), *Hyperbolicity in scalar-Gauss-Bonnet gravity: a gauge invariant study for spherical evolution*](https://arxiv.org/abs/2410.16264) | Gauge-invariant tracking ties hyperbolicity loss to physical degrees of freedom and shows that derivative-dependent field redefinitions, including the studied disformal example, can alter the initial-value formulation and apparent equation character; this is not an indiscriminate statement about every invertible linear transformation. | Freeze the evolved variables and physical metric, report the unredefined system's principal symbol, and do not treat a formulation-dependent discriminant as a theory-wide certificate. |
| [Thaalba et al. (2024), *The dynamics of spherically symmetric black holes in scalar-Gauss-Bonnet gravity with a Ricci coupling*](https://arxiv.org/abs/2409.11398) | Nonlinear black-hole evolutions use excision; the study relates an interior elliptic region to a finite-area singularity and extracts scalar modes where possible. | Excision of a region after the original equations become elliptic cannot establish the proposed healthy restoring interval there.  Any retained FGC interval must precede that boundary in the original system. |
| [Ye, Liu & Zhang (2026), *Spontaneous scalarization and dynamical evolution of black holes in scalar-Gauss-Bonnet gravity*](https://doi.org/10.1140/epjc/s10052-025-15272-w) | Ricci/GB-coupled nonlinear evolutions already report scalarization thresholds, Misner--Sharp energy transport, and null-convergence-condition violation, while dynamically excising elliptic regions. | Neither a threshold, Misner--Sharp ledger, nor null-convergence violation is new by itself.  The remaining distinction is independent `chi` collapse plus a converged affine defocusing margin wholly inside the original system's predictive domain. |

[FGC-1-HYP1-SYM1](fgc-hyp1-symbol.md) closes the physical quotient burden for
the frozen local fixtures. [FGC-1-HYP1-MHG1](fgc-hyp1-modified-harmonic.md)
then applies the Kovacs--Reall architecture directly to ACT1's unredefined
VAR1 symbol: it retains the independent `chi`, dynamic `F(phi)`, and derived
regulator factor; adds two separated auxiliary cones; and proves pointwise
complete real eigenbases for the full radial principal matrices. This is
source-constrained but not literature-derived: the ACT1 coefficients and
factorization are recalculated exactly. The implementation also stays in the
smooth VAR1 representation at `X=0` without invoking a standard-basis rewrite
or any `log X` coefficient. The downstream
[FGC-1-HYP1-MHG2-REF1](fgc-hyp1-mhg-reference.md) gate now freezes a tensorial
flat-spherical reference gauge, evaluates the full gauge extension, and
recovers MHG1's principal gauge block entry by entry. The downstream
[FGC-1-HYP1-MHG3-IMP1](fgc-hyp1-mhg-implicit.md) gate then uses exact forward
tangents on the complete REF1 residual to prove one nonsingular
coordinate-time acceleration Jacobian and finite-dimensional IFT branch at an
exact flat inactive FGC-QR vacuum/reference-gauge root, with the nonzero action
couplings retained. The downstream
[FGC-1-HYP1-MHG4-PROP1](fgc-hyp1-mhg-propagation.md) gate then evaluates the
lower-order `nabla(F hat_P nabla C)` operator on independent spherical
`C^mu`, homogeneous only after the ACT1 Noether handoff is combined with both
unmodified scalar equations; it uses two exact routes, extracts every spherical coefficient,
and regresses the principal columns to MHG1. It is a project-specific exact
operator certificate, not a claim that the cited literature or PROP1 proves
metric-derived constraint propagation. Kovacs--Reall motivate the
separation between kinetic invertibility and the later open weak-coupling
hyperbolicity proof; they do not establish this ACT1/REF1 Jacobian or PROP1's
lower-order coefficients. PROP1 remains upstream of direct metric-derived
gauge/reduction-constraint propagation, a quantified open-domain health/EFT proof, collapse, and the affine
Raychaudhuri margin that define the publication wedge.

The downstream [FGC-1-HYP1-DOM1-QIFT1](fgc-hyp1-dom1-qift1.md) result is a
self-contained project computation: exact rational interval evaluation of the
complete REF1 residual proves a unique acceleration root throughout one
declared full-dimensional local box. It needs no new literature citation.
Kovacs--Reall remain motivation for separating a local algebraic branch from a
weak-coupling hyperbolicity theorem; they neither prove QIFT1's ACT1/REF1 box
nor convert it into a constraint, nonflat-hyperbolicity, EFT, or evolution
result.

The further [FGC-1-HYP1-CON1-COMP1](fgc-hyp1-con1-comp1.md) result is likewise
a project-specific exact local theorem. At one nonflat activated mixed/spatial
two-jet, two selected spatial jets solve
`H_ww=0` and `nabla_r C^r=0`; the frozen choices and exact coefficient
identities give `M=C^mu=nabla_r C^t=0`. Composing that datum with QIFT1's
unique acceleration root and the nonsingular normal gauge-extension map then
forces the remaining gauge derivatives and recovers the unredefined metric
and both scalar equations at that point. Kovacs--Reall supply the
formulation architecture, Gundlach and Martín-García supply the standard
distinction between gauge and first-order reduction constraints, and East--
Ripley supply a neighboring compatible-data/evolution comparator; none of
those works proves COMP1's ACT1 datum. Conversely, one compatible local
two-jet is not a constraint-propagation theorem, a regular-center or
asymptotically admissible initial slice, a uniform open-domain symmetrizer, a
retained-EFT validity envelope, or an evolution result.

These comparisons were resolved against the linked primary records on
2026-08-21.  They remain a **refresh obligation** because versions, nearby
models, and priority can change before release.

There is also a convention-sensitive sign boundary. Equation (1) of Thaalba
et al. (2023) contains a negative coefficient for its published positive
Ricci parameter in the `R phi^2` term. At the level of signs, that corresponds
to **negative** `beta` in the present
`F(phi)=M_Pl^2+beta phi^2` convention; coefficient magnitudes also depend on
the overall and scalar normalizations. FGC-1-ACT1 therefore permits signed
nonzero `beta` and uses a negative algebra fixture. It does not import the
paper's hyperbolicity result into this differently normalized action. The
effective-Planck boundary `F(phi)>0` and the full principal symbol must be
checked anew.

## Exact overlap, and the remaining publishable gap

Already saturated by the literature:

- curvature-dependent scalar instabilities and bifurcating hairy branches;
- nonlinear saturation/stability as coupling-dependent questions;
- mixed-type or non-strongly-hyperbolic failure of apparently benign
  higher-curvature scalar--tensor equations;
- numerical strong-field and waveform phenomenology in nearby theories.

The remaining, potentially publishable, result is narrower and must be shown
from one action:

1. `chi` is an independently specified matter-collapse sector, with finite
   spherical initial data and a GR recovery limit;
2. an invariant such as `I[g, chi]` controls a regulator-sector threshold,
   e.g. `m_phi,eff^2(I)` changes sign or an auxiliary second-order equation
   activates, without a coordinate-density trigger;
3. the regulator nonlinearly saturates at a finite scale and all relevant
   spherical characteristics stay real, finite, and with a declared physical
   cone;
4. the calculation reports the complete affine Raychaudhuri right-hand side,
   including shear and `R_ab[g] k^a k^b`, before making any focusing or
   defocusing statement; and
5. the code reproduces the branch under refinement and produces one
   GR-baseline discriminator (for example a threshold, transient waveform,
   compactness response, or exterior multipole).

This is a proposed *combination and test package*, not a claim of broad
theoretical priority.  It becomes novel only if its action, healthy domain,
and discriminator differ concretely from the cited models.

## Declared fallbacks and comparators

**Quadratic Palatini/Ricci-based gravity — fallback regular-core route.**  A
restricted independent-connection action can be written

```text
S = (1/2 kappa^2) integral sqrt(-g)
    [R(g,Gamma) + ell^2(a R^2 + b R_(ab) R^(ab))] d^4x + S_m[g,chi].
```

In the specified projectively invariant Ricci-based subclass, the connection
is auxiliary after the matter relation is solved and the metric equations are
second order; this does not license a claim about arbitrary metric-affine
models.  Static charged cores/wormholes are known, sometimes conditionally
regular and sometimes geodesically complete despite a curvature divergence.
It is a regular-core comparator, not current evidence for an `chi -> phi`
phase transition or generic collapse completion.  Primary sources:
[Olmo & Rubiera-Garcia (2011)](https://arxiv.org/abs/1112.0475),
[Olmo & Rubiera-Garcia (2013)](https://arxiv.org/abs/1301.2430), and
[Olmo, Rubiera-Garcia & Sanchez-Puente (2015)](https://arxiv.org/abs/1508.03272).

**Extended mimetic limiting curvature — comparator only.**  A representative
form is `S = integral sqrt(-g)[M_Pl^2 R/2 + L(X+1)+F(Box phi)] d^4x` with a
mimetic clock constraint.  It is useful direct precedent for a model-specific
limiting-curvature black-hole interior, but its preferred foliation and known
generic gradient/extra-mode hazards preclude treating it as the healthy FGC
default.  Primary sources: [Ben Achour et al. (2017)](https://arxiv.org/abs/1712.03876),
[Langlois et al. (2017)](https://arxiv.org/abs/1708.02951), and
[Frolov, Markov & Mukhanov (1990)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.41.383).

## Approved and disallowed public wording

**Use:** “We study a toy covariant `chi`--`phi` model in which a specified
invariant can activate a regulator sector during spherical collapse.  We test
the physical characteristics and the full affine Raychaudhuri source before
asking whether a finite transition exists.”

**Do not use:** “first finite-gradient theory,” “curvature automatically
prevents singularities,” “a positive regulator term proves defocusing,”
“healthy because the equations are second order,” “black holes create child
universes,” “waves become particles,” or any claim that a static core, an
effective stress tensor, or a coordinate threshold establishes the required
collapse result.

**First-order reduction scope.** FGC-1-HYP1-FO1-RC1 uses the standard
auxiliary-constraint distinction exemplified by Gundlach and Martín-García,
[arXiv:gr-qc/0506037](https://arxiv.org/abs/gr-qc/0506037), Sec. III,
Eqs. 15--21. That paper is precedent for separating derivative-reduction
constraints from physical constraints and for reduction-choice-dependent
homogeneous subsystems; it does not derive ACT1/REF1, the FO1 coefficient
matrices, or this project's metric-gauge and physical-constraint propagation.

## Release-time novelty refresh

Within seven days of a public release, and again on the release day:

1. resolve the Hegade/Ripley/Yunes, Thaalba, and Ye/Liu/Zhang entries to exact
   version-of-record metadata and read their abstract, action, principal-part
   treatment, collapse setup, and observables;
2. search INSPIRE/arXiv and the relevant journal databases for `scalarization`,
   `Gauss-Bonnet`, `hyperbolicity`, `spherical collapse`, `limiting curvature`,
   and the exact new action terms; date-stamp the query;
3. update this table with an explicit overlap/no-overlap decision for the
   action, activation invariant, health proof, nonlinear solution, and
   observable; and
4. remove every priority adjective unless the comparison supports a narrowly
   worded one.  The release should cite the closest positive and negative
   precedents, including hyperbolicity failures, rather than presenting them
   as omissions by prior work.
