# FGC-1-SGB1-CTL1 — branch-owned SGB-L controls (prospective)

This is an implementation map, not a completed certificate. The accepted
PLAN.md permits this outcome-blind Wave-0 work alongside numerical recovery.
No SGB-L trajectory, holdout or candidate execution is authorized here.

## Required complete owner

SGB1-CTL1 must eventually bind matched initial data, a branch-owned source
and runtime, independent source/Jacobian checks, and kinetic, cone, continuity,
regular-center, constraint and health evidence. Weak-field, constant-coupling,
Schwarzschild and established SGB controls remain required. FGC-QR health
evidence cannot be borrowed. The annular constraint kernel and the
pointwise source/Jacobian algebra below close two implementation
dependencies; they are not the full owner.

The action is the existing linear branch `f(phi)=alpha_gb*phi`, with
`beta=eta=0`, not a new force or a modification chosen from an outcome.
`AGENTS.md` remains unchanged creative context. The action, RED1 residuals,
and ID1 slice—not that creative document—own the calculations.

## Annular constraint kernel and structural proof

Use unit lapse, zero shift, `R=r` on the initial slice, and
`K^r_r=-2k`, `K^theta_theta=K^varphi_varphi=k`. This means
`R_t=-r*k`; it does not mean the time derivative vanishes.
For `r>0`, `L=lambda>0`, define

```text
A = L_r/(L^3*r) - 2*k^2
B = k^2 + (1-L^-2)/r^2
C = (k_r+3*k/r)/L
X = phi_rr/L^2 - L_r*phi_r/L^3 - 2*k*Pi_phi
Y = k*Pi_phi + phi_r/(L^2*r)
Z = (Pi_phi_r-2*k*phi_r)/L.
```

These follow from the actual warped-product curvature and scalar Hessian.
RED1's double dual gives `P0101=-B`, `P0202=P0303=-A`,
`P0110=B`, `P0212=P0313=-C`, and `G00=2A+B`, `G01=2C`.
Mixed Hessian indices are raised with signature `-+++`. Thus the physical
projections (`H=E00`, `M=Etr=L*E_orth01`) are

```text
rho = (Pi_phi^2+Pi_chi^2)/2 + (phi_r^2+chi_r^2)/(2*L^2)
      + mu^2*phi^2/2 + g4*phi^4/4
H = Mpl^2*(2*A+B) - rho - 8*alpha_gb*(B*X+2*A*Y)
M = 2*Mpl^2*L*C - Pi_phi*phi_r - Pi_chi*chi_r
    - 8*alpha_gb*L*(B*Z+2*C*Y).
```

`A` and `X` are affine in `L_r` and never multiply each other. Only `C`
contains `k_r`, linearly. This proves `H=H0+H_L*L_r`,
`M=M0+M_k*k_r`, with vanishing cross and higher-degree coefficients.
Neither scalar momentum nor the independent chi stress is set to zero.
The plus `8 alpha_gb P Hess(phi)` term is on RED1's residual left-hand
side; it is consistent with VAR1's negative right-hand-side GB stress.

Sources: `spherical_reduction.py` (`warped_2plus2_curvature`, `_p`,
`residuals`), `modified_harmonic_constraints.py`
(`physical_constraint_projections`), and `initial_data_family.py`'s ID1
slice and gauge-compatible metric time derivative. Preserve those files.

## File/code implementation map

First bounded slice (annular constraint kernel; committed, leave unchanged):

- `src/recursive_horizons/fgc/sgb1_ctl1_constraints.py`: exact rational
  annular inputs, the four coefficients, complete H/M reconstruction, and
  an optional exact two-derivative solve with a typed singular-Jacobian stop.
- `tests/test_fgc_sgb1_ctl1_constraints.py`: independent direct RED1
  projections, a symbolic component/degree identity, mixed scalar stresses,
  nonzero momenta, GR/flat limits and exact singular controls.

Public inputs use `planck_mass`, not its square; the formulas square it
explicitly. `alpha_gb=0` is permitted only as the algebraic GR-limit control,
not as an SGB-L production branch. Fractions and built-in integers are
accepted; float/bool aliases and invalid annular metrics are refused.

The algebraic derivative solve establishes only local constraint zeros when
both diagonal coefficients are nonzero. It must not invent a quantitative
Jacobian floor, claim a neighborhood/existence theorem, or evaluate `1/r`
at the center. No radial ODE integration, actual initial-data family, source
trajectory or production-state codec is part of this first slice.

Second bounded slice (pointwise source/Jacobian algebra):

- `src/recursive_horizons/fgc/sgb1_ctl1_source.py`: exact annular source-point
  inputs, seven-evaluation `R0` and six-by-six acceleration Jacobian,
  complete affine reconstruction, exact determinant, and an optional exact
  solve with a separately named singular-source-Jacobian stop.
- `tests/test_fgc_sgb1_ctl1_source.py`: independent 28-node quadratic-support
  controls, holdout reconstruction, actual residual re-evaluation at the
  solved root, unredefined six-row determinant, coupling-limit and
  Schwarzschild controls, malformed-input refusal, and unforgeable records.

The source Jacobian is not the two-dimensional annular constraint Jacobian.

Third bounded slice (matched family and initial center):

- `src/recursive_horizons/fgc/sgb1_ctl1_family.py`: lambda-free compact
  profiles, affine H/M constraint ODEs, named radial constraint integrators,
  exact Minkowski buffer, analytic exterior formula with floating evaluation,
  sampled-versus-continuous trapping language, and a closed health gate.
- `src/recursive_horizons/fgc/sgb1_ctl1_center.py`: SGB-L initial-center
  Laurent/jet validator, empty-buffer series, first-grid and elementary-flatness
  controls.
- `tests/test_fgc_sgb1_ctl1_family.py`, `tests/test_fgc_sgb1_ctl1_center.py`.

Fourth bounded slice (branch-owned point principal adapter):

- `src/recursive_horizons/fgc/sgb1_ctl1_principal.py`: maps a complete SGB-L
  spherical two-jet to the universal covariant ACT1 principal background with
  `F=Mpl^2`, `F'=0`, `f'=alpha_gb`, and
  `nabla nabla f=alpha_gb*nabla nabla phi`.
- `tests/test_fgc_sgb1_ctl1_principal.py`: compares the full covariant radial
  restriction against the independently assembled exact spherical SGB-L MHG
  symbol on a nonflat, nonzero-shift, two-scalar solved source point.

This is code reuse of branch-neutral tensor algebra, not reuse of an FGC-QR
health certificate. It is one point comparison, not an interval, cone,
symmetrizer, source-runtime or branch-health result.

Fifth bounded slice (honest point coefficient/directional facts):

- `src/recursive_horizons/fgc/sgb1_ctl1_principal.py`:
  `sgbl_principal_point_facts` retains the universal 12-field coefficient
  tensors, hashes, radial/angular evaluated symbols and binary64 spectra
  at orthonormal `(1,0,0)`, `(0,1,0)` and `(3/5,4/5,0)`. Health, pass, cone
  and interval fields are absent as constructor data and are explicitly
  false as derived properties. The FGC-QR spherical builder and
  weak-coupling pass bits are not imported.
- `tests/test_fgc_sgb1_ctl1_principal.py`: radial restriction still matches
  the exact spherical SGB-L symbol; radial and angular spatial tensors
  differ; sampled spectra do not set health true.

Sixth bounded slice (declared-box ADM source-Jacobian invertibility):

- `src/recursive_horizons/fgc/sgb1_ctl1_interval_health.py`: encloses
  `J=d(complete MHG rows)/d(alpha_tt,shift_tt,lambda_tt,R_tt,phi_tt,chi_tt)`
  on a caller-declared lower-jet box with arithmetic/resource caps. The
  inverse is a Neumann certificate only when a predeclared `rho_∞<1`
  condition proves it. `rho>=1` wrapping is typed
  `interval_inconclusive`. An exactly singular center keeps
  `singular_source_jacobian`. This is not SRC1's metric-`dtt` Jacobian,
  not a frozen physical production width, and not a trajectory.
- `tests/test_fgc_sgb1_ctl1_interval_health.py`: exact singleton reduction,
  a proved `rho<1` enclosure with sampled containment and inverse residual
  bound, a wide wrapping inconclusive, distinct exact-singular ownership,
  and injected chart/sign/resource failures.

Interval invertibility on a declared box is an instrument. It does not
open `SGBL_branch_owned_and_healthy`. All-direction cone/symmetrizer
enclosure and the evolving centre after matter arrives remain separate
gates. A compact result and live execution target are not created from
this source, family, principal or interval slice.

Seventh bounded slice (branch-owned source/runtime):

- `src/recursive_horizons/fgc/sgb1_ctl1_runtime.py`: floating and small-grid
  evaluation of the same complete MHG six-row residual, ADM 6x6 Jacobian
  and inverse (or typed singular/inconclusive), method-owned SBP
  reconstruction of `(p_r,q_r)`, both numerical-engine time methods, and a
  finite-state capture/restore transaction.
- `tests/test_fgc_sgb1_ctl1_runtime.py`: exact nonflat/nonzero-shift/
  two-scalar floating comparison, injected chart/action/sign controls,
  missing-jet refusal, both time methods on a 9-point annulus, and
  rollback fingerprints for source/CFL/health stops.

`EvolutionState` stores only `(u,p,q)`. The source algebra needs the
complete lower two-jet `(u,p,q,p_r,q_r)`. Those last two slots are never
filled with zeros and are never borrowed from an FGC-QR reference
projection. A caller must supply the arrays or bind the method-owned SBP
operator that defines `p_r=D_r p` and `q_r=D_r q`. An `EvolutionState`
without that completion is `missing_lower_two_jets`, not a runnable
smaller source.

RK4/SSPRK3 in this module are the numerical-engine time methods with
their method-owned spatial operators. The family's identically named
radial constraint integrators remain a different owner. Grids are
synthetic and annular (`r>1/2`, at most 17 points). This slice does not
qualify the evolving centre, enclose the cone, or set
`SGBL_branch_owned_and_healthy`.

## Highest-time source structure

The source problem is separate from the annular constraint solve. Holding
field values, first derivatives, mixed derivatives and radial second
derivatives fixed, consider the six ADM coordinate-time accelerations.
The ADM-to-metric second-derivative map is affine on that subspace.
Write `ell=dt`, `v=g^{-1}ell`. The acceleration-dependent Riemann variation
contains one `ell` in each antisymmetric index pair. In a normal orthonormal
frame this is the arbitrary symmetric block `delta R_0i0j=A_ij`.

Two such curvature variations give zero in the Gauss--Bonnet quadratic
contraction: the generalized antisymmetrization repeats `ell`. Likewise
`delta P_acbd v^c v^d=0` for the double dual. The only potential quadratic
metric/scalar cross is therefore zero:

```text
delta P_acbd * delta Hess(f)^{cd}
  = alpha_gb * delta(phi_tt) * delta P_acbd * v^c * v^d = 0.
```

Einstein, scalar-wave and modified-harmonic gauge terms are affine at fixed
lower jets; scalar stress and the potential have no accelerations. Thus
these curvature identities remove the possible quadratic highest-time
terms in the SGB-L equations. They do not say that solving the source is
linear in the lower jets, or that its acceleration Jacobian is nonsingular.

`tests/test_fgc_sgb1_ctl1_time_affinity.py` checks the two curvature identities
as exact formal polynomials in all six independent symmetric `A_ij` atoms
using RED1's actual contractions. This is not an interpolation over numerical
samples. It also pulls the tensor into a nonzero-shift chart and verifies
that the twice-raised time-covector contraction vanishes while the naive
coordinate `P_theta0theta0` need not. Tensor covariance with timelike `dt`
connects the normal-frame polynomial identity to general ADM coordinates.

The pointwise source algebra below uses that identity. Quantitative
invertibility, source-domain qualification, regular-center continuation,
source runtime and branch kinetic/cone health remain separate obligations.

## Source acceleration algebra

Holding the lower jet `z` fixed, the complete six-row MHG residual of the
linear branch is affine in the six ADM coordinate-time accelerations:

```text
R(a; z) = R0(z) + J(z) a,
a = (alpha_tt, shift_tt, lambda_tt, R_tt, phi_tt, chi_tt).
```

The map `z |-> a` remains nonlinear. Invertibility of `J(z)` is a separate
question from affinity: a nonzero exact determinant at one point is not a
neighborhood theorem, a conditioning floor, or kinetic health.

The reference instrument is `modified_harmonic_full_residuals` with the
frozen flat spherical annulus `r>1/2`, `tilde=4`, `hat=9`, and
`beta=eta=0`. `R0` and `J` are the exact constant and six unit-basis
evaluations of that residual; the highest-time identity makes those seven
calls exact rather than a fit. A complete reconstruction at arbitrary
accelerations is `R0+J a`. An optional exact solve inverts `J` when
`det J != 0` and then re-evaluates the actual full residual at the algebraic
root, requiring the exact zero vector. The typed stop
`SGBLSourceJacobianSolveStop` (`singular_source_jacobian`) is distinct
from the annular `SGBLConstraintSolveStop`.

Lower-jet inputs are Fraction or built-in int, or exact `Jet2` components
with `dtt=0`. Public five-tuples are ordered `(value, dt, dr, dtr, drr)`.
A supplied nonzero `dtt` is refused rather than dropped. Lapse, ADM
`lambda`, areal radius `R`, coordinate radius, `Mpl`, `mu` and `g4` are
positive. `alpha_gb` is nonzero on the SGB-L branch; `alpha_gb=0` is only
an algebraic coupling-limit control.

Independent 28-node quadratic-support controls uniquely determine a
quadratic in six accelerations (constant, six linear, twenty-one quadratic
monomials). All twenty-one quadratic coefficients vanish in all six rows
on the exact fixtures. A holdout vector that is not a support node checks
reconstruction rather than the affine arithmetic `J a = -R0` alone.

Exact fixture controls, with `dtt=0` in every lower jet:

- A: coordinate radius `5/2`, `Mpl=2`, `mu=3`, `g4=1/2`,
  `alpha_gb=-1/4`; `alpha=(2,1/3,-1/2,1/5,-1/7)`,
  `shift=(1/3,-2/5,1/4,-1/6,2/9)`, `lambda=(3/2,1/8,-1/3,1/7,-1/4)`,
  `R=(4,-1/5,3/2,1/9,-2/5)`, `phi=(1/2,2/3,-4/5,1/4,-1/8)`,
  `chi=(-3/4,1/2,3/5,-2/7,1/6)`. `det J_A = 6519803496633/6553600000`.
- B: coordinate radius `3`, `Mpl=3/2`, `mu=2/3`, `g4=5/7`,
  `alpha_gb=2/5`; `alpha=(5/2,-2/7,3/4,-1/8,2/11)`,
  `shift=(-2/5,1/6,-3/7,2/9,-1/5)`, `lambda=(4/3,-1/4,2/5,-1/9,3/8)`,
  `R=(7/2,3/8,-5/4,2/7,1/10)`, `phi=(-2/3,-5/8,1/2,-3/11,2/13)`,
  `chi=(4/5,-1/3,-2/9,3/10,-1/12)`.
  `det J_B = 57291897169490646528/7593631744384765625`.
- The unredefined six-row acceleration Jacobian has determinant `0` on both
  fixtures; without the MHG rows the source map is not a six-by-six solve.
- `alpha_gb=0` on A gives `det J = 2187/2`. Nonzero `phi` means this is not
  the canonical `phi=0` GR-0 collapse control.
- Holdout `a=(2,-1,1/2,3,-2,1)` reconstructs the actual residual exactly.
- Flat vacuum with `alpha_gb != 0` is an exact zero source; Schwarzschild
  with `phi=0` and `alpha_gb != 0` has unredefined scalar source
  `alpha_gb*GB` and is not a full SGB-L solution.
- A constructed off-shell point with flat metric, `R=r=5/2`, `phi_r=-5`,
  zero other scalar data, `Mpl=2`, `mu=g4=1` and `alpha_gb=-1/4` has an
  exactly singular full source Jacobian. It exercises the actual
  `singular_source_jacobian` refusal, not just an arbitrary test matrix.
  It is neither a constraint-satisfying solution nor a physical singularity.

Public coefficient records store only `R0` and `J`. Construction and
`replace` cannot attach coefficient-identity or solved-residual claims;
the determinant is a derived property of the stored matrix.

This source slice is still not a source runtime, interval invertibility
certificate, or branch kinetic/cone health owner. Matched family and
initial-center controls are a later bounded slice in this owner, not a
promotion of the pointwise Jacobian. It does not open a compact result
or an execution target.

## First-slice proof and failure forecast

Require exact agreement with the unredefined RED1 projections, symbolic
affine/diagonal support, preserved chi and Pi_phi terms, and deterministic
typed refusal at a singular Jacobian. `r=5/2`, `L=1`, `k=Pi_phi=0`,
`alpha_gb=-1/4`, `phi_r=-5` is a singular control for `Mpl=2`, not `Mpl=1`.
Run focused exact tests and Ruff; make no commit before coordinator review.

Likely failures are a convention/sign mismatch (compare actual RED1
components, not an isolated Codazzi mnemonic), omitted scalar stress or
momentum terms (independent projections with all such inputs nonzero), and
overgeneralizing annular nondegeneracy to health claims. Family and
initial-center controls are implemented as a later bounded slice; they
do not set branch-owned health. Even a complete SGB1-CTL1 owner would
not itself open the holdout or establish defocusing.

## Source-slice proof and remaining obligations

Require the seven-evaluation Jacobian to reconstruct the actual MHG
residual at a holdout acceleration, all twenty-one quadratic coefficients
to vanish in all six rows on both fixtures, the solved root to make the
re-evaluated residual exactly zero, and the declared exact determinants.
Do not treat `J a = -R0` arithmetic as a substitute for residual
re-evaluation. Do not resurrect the disproved claim that `alpha_gb*GB`
makes the residual quadratic in accelerations. Do not borrow FGC-QR
health or source-domain evidence.

Likely failures are silently dropping a supplied `dtt` from a lower jet,
confusing the six-by-six source Jacobian with the two-by-two annular
constraint Jacobian, treating `alpha_gb=0` on a nonzero-`phi` point as
the canonical GR-0 collapse control, or converting a pointwise exact
`det J != 0` into a neighborhood, interval-invertibility or kinetic-health
claim. Interval invertibility, source runtime and branch kinetic/cone
health remain unqualified. This slice creates no compact result and no
execution target.

## Matched family and initial-center slice

Third bounded slice (matched initial family and initial-center controls):

- `src/recursive_horizons/fgc/sgb1_ctl1_family.py`: declared compact SGB-L
  profiles, affine H/M radial constraint ODEs, named RK4/SSPRK3 constraint
  integrators on a small fixed mesh, exact interior Minkowski buffer,
  initial GR vacuum outer continuation, independent tensor/FD checks, and
  a closed health gate.
- `src/recursive_horizons/fgc/sgb1_ctl1_center.py`: SGB-L initial-center
  profile/series validator using neutral `LaurentSeries`/`SeriesJet2`,
  exact empty Minkowski buffer, nontrivial regular first-grid control, and
  elementary-flatness/odd-scalar failure controls.
- `tests/test_fgc_sgb1_ctl1_family.py`, `tests/test_fgc_sgb1_ctl1_center.py`.

The action and seed are the declared linear branch

```text
alpha_gb = -1/4,  beta = eta = 0,  Mpl = 2,  mu = 3,  g4 = 1/2,
A_phi = 1/131072,  center = 12,  half-width = 2,  outer = 128.
```

`A_chi` is a declared family input with nominal value 3. It is not selected
here and is not an ID1-box vertex. Profiles use `compact_bump_with_derivatives`
(or the same checked formulas). Live `compact_family_fields` takes no
`lambda`; on the unit-lapse/zero-shift slice

```text
Pi_chi = chi_t = A_chi B_r / r
```

is fixed independently of the solved geometry and is not rematched after the
branch solve. This module does not call `FGCQRActionParameters` or
`solve_initial_data`, and does not borrow an ID1 family or health certificate.

The two radial ODEs are the branch affine pair already proved in the annular
kernel. RK4 and SSPRK3 are radial constraint integrators, not a time method
and not a trajectory. Meshes larger than the small fixed budget raise a typed
`resource_limit` stop. An exact `Fraction` zero diagonal is
`singular_constraint_jacobian`. A floating `0.0` diagonal is
`uncertified_floating_jacobian` unless the exact kernel on the same finite
binary-rational inputs is itself singular. There is no conditioning floor.
Nonpositive metric is separately typed. Defining-RHS zeros are not the
constraint check: stored `lambda(r)` and `k(r)` are reconstructed by finite
differences and compared through the affine H/M formulas, and exact rational
points are compared to independent RED1 tensor projections. A known vacuum
constraint solution, compared to the analytic exterior formula, is the
radial-order control; a cross-method mass tolerance is not.

The interior of the compact support is the exact Minkowski buffer
`lambda=1`, `k=0`. Outside the support the analytic GR vacuum formula

```text
k = J / r^3,
lambda^{-2} = 1 - 2M/r + J^2/r^4
```

is evaluated in floating arithmetic. The formula is not the rounded
evaluation. That continuation satisfies the initial constraints when the
scalars vanish; it is not a static full SGB-L solution, because
`alpha_gb*GB` can source `phi` acceleration there. Sampled-node compactness
`<1` on at most 64 steps is `sampled_nodes_untrapped` only.
`continuous_no_initial_trapped_sphere` remains unqualified. No SGB-L pass is
inferred from an FGC-QR family-box certificate.

Public family records store nodes and parameters. Forbidden scope claims
(`SGBL_branch_owned_and_healthy`, full static SGB-L, continuous no-trap,
exact floating exterior) are derived properties and cannot be attached by
construction or `replace`. Finite/shape/mesh/derived compactness-mass
inconsistencies are refused. Each node's stored scalar profiles and defining
RHS metadata are recomputed from its declared family and geometry on
construction; this prevents relabelling foreign scalar data as the empty
buffer. That consistency check is not the independent finite-difference
constraint test or a certificate of the integrated trajectory.

The initial-center owner reuses Laurent/jet arithmetic, not
`validate_regular_profile` or `regular_center_series_certificate`. Empty
Minkowski is a coefficient/geometry fact, never a `profile_id` prefix. A
nontrivial even regular profile is only a first-grid/elementary-flatness
instrument. Cached `formal_series` must reproduce the complete six-equation
exact vector of the profile; incomplete, foreign, or altered caches are
refused. The source algebra's frozen annulus minimum `1/2` does not license
an epsilon-at-zero shortcut: pointwise RED1 requires `r>0`, and `r=0` is a
Laurent chart. The empty initial-center proof does not qualify the
evolving centre after matter arrives.

This slice does not set `SGBL_branch_owned_and_healthy`. Missing health
closes that holdout gate; it is not itself a proof of nonhyperbolicity.
No production-grid family sweep or trajectory is owned here.

## Family-slice proof and remaining obligations

Require lambda-free profile matching, `Pi_chi=chi_t=A_chi B_r/r` before and
after each branch solve, exact vacuum tensor H/M zeros, preserved chi
stresses and momenta, positive metric, typed exact-singular versus
uncertified-floating Jacobian stops, sampled-versus-continuous trapping
language, independent FD reconstruction that is not the defining RHS, a
known vacuum-constraint refinement control, emptiness from coefficients
rather than names, refusal of incomplete/foreign/altered first-grid
caches, an exact empty initial-center series, a nontrivial regular
first-grid control, and explicit elementary-flatness/pole failures.

Likely failures are rematching `Pi_chi` with solved `lambda`, treating
defining-RHS zeros as a tensor check, evaluating `1/r` at the origin or
borrowing the source annulus minimum as a centre chart, calling a floating
exterior evaluation an exact continuation or a static SGB-L solution,
labelling a renamed nonempty profile empty, accepting a partial
`formal_series` cache, promoting sampled-node untrappedness to a continuous
no-trap gate, treating a floating zero diagonal as the exact kernel
singularity, or setting branch-owned health true because a conservative
envelope was not violated.

The declared-box ADM Jacobian instrument can prove invertibility on a
caller-chosen lower-jet box, or return typed wrapping-inconclusive. That
is not all-direction cone health and does not set
`SGBL_branch_owned_and_healthy`. A synthetic source runtime with rollback
is implemented; all-direction cone/symmetrizer enclosure, the evolving
centre after matter arrives, and continuous no-initial-trap enclosure
remain unqualified. This owner creates no compact result and no
execution target.

## Principal-fact and interval-slice proof and remaining obligations

Require the 12-field adapter to consume a correctly built SGB-L
`CovariantPrincipalBackground` with `F=Mpl^2`, `F'=0`, `f'=alpha_gb` and
`Hess(f)=alpha_gb Hess(phi)`; retain coefficient tensors and hashes; match
the existing exact radial spherical symbol; show that the angular spatial
symbol differs; and keep every health/pass/cone/interval field false.
Binary64 directional spectra remain diagnostics with uncertainty visible.

Require the interval owner to evaluate the ADM six-by-six source Jacobian,
not SRC1 metric-`dtt`. An exact singleton box must reproduce the frozen
source `J` and determinant. A nonzero declared box enclosure must contain
independent sampled Jacobians. A proved `rho<1` case must verify the
inverse residual bound. A wide wrapping box is `interval_inconclusive`.
An exact singular center remains `singular_source_jacobian`. Injected
chart, sign and resource failures must not be relabelled as that
singularity. Do not import the FGC-QR spherical builder, QIFT1, SRC1
solve, weak-coupling pass bits or an old calibration runner.

Likely failures, and how the tests distinguish them:

1. Chart error: seeding metric `dtt` (`h_tt,h_tr,h_rr,…`) produces a
   different six-by-six from ADM `J`; requesting that chart is
   `interval_chart_error`. Lorentzian/positivity loss on a too-wide box
   is `interval_chart_domain`.
2. Interval wrapping: the enclosure still contains sampled exact
   Jacobians, but `rho>=1`, so invertibility is `interval_inconclusive`
   rather than a scientific nonpass or a fitted floor.
3. Genuine singularity: the exact center has `det J=0` and raises the
   existing `singular_source_jacobian` before any interval inverse is
   attempted.

Remaining gates after the interval slice were source runtime with
rollback, all-direction cone/symmetrizer enclosure, evolving centre
after matter arrives, continuous no-initial-trap enclosure, and
`SGBL_branch_owned_and_healthy`. The source-runtime slice below closes
only the synthetic owner; the other gates remain open. Missing health
still closes holdout; it does not prove nonhyperbolicity.

## Source-runtime proof and remaining obligations

Require the floating residual, Jacobian and algebraic root to match the
exact `sgbl_source_state`/`sgbl_source_solve` witness on the nonflat,
nonzero-shift, two-scalar fixtures; refuse metric-`dtt`, `beta!=0`/`eta!=0`
and a sign-flipped `alpha_gb`; keep exact `det J=0` as
`singular_source_jacobian`; refuse `EvolutionState`-only evaluation
without declared `(p_r,q_r)` or a method-owned operator; recover those
jets from the bound SBP operator on the polynomial embedding; run both
time methods on the small annulus; and restore state, monitor, causal,
tracers and counters on source, CFL and health stops. Point invertibility
is not branch health.

The grid adapter returns a floating affine candidate together with absolute
and coefficient-scaled residual witnesses. This slice deliberately records
`source_solve_admission_qualified=false`: it introduces no fitted absolute
residual threshold and does not promote a merely finite residual to a
production source solve. Failure of a binary64 grid inversion is
`source_inversion_inconclusive`; `singular_source_jacobian` remains reserved
for the exact point owner with exact `det J=0`. Public coefficient/grid arrays
are immutable. After a synthetic step, the initial declared `(p_r,q_r)` arrays
are cleared because current lower jets are reconstructed by the bound SBP
operator; restoring the captured initial member restores the original
declarations exactly.

Likely failures, and how the tests distinguish them:

1. Silent zero jets: an `EvolutionState` without `(p_r,q_r)` or an
   operator raises `missing_lower_two_jets`. Explicit zero arrays produce
   a different acceleration from the declared mixed/second-radial jets.
2. Relabelled source: AST import controls refuse SRC1
   `solve_accelerations`, GR0 vectorized source, `FGCQRActionParameters`
   and old `run_fgc*` runners.
3. Chart or action mix-up: `metric_dtt` is `source_chart_error`; nonzero
   `beta`/`eta` is `action_identity_error`; flipping `alpha_gb` changes
   the residual relative to the exact fixture.
4. Constraint-integrator confusion: family `RK4`/`SSPRK3` names are not
   the numerical-engine time methods and cannot construct a member.

Remaining gates: all-direction cone/symmetrizer enclosure, the evolving
centre after matter arrives, continuous no-initial-trap enclosure, and
`SGBL_branch_owned_and_healthy`. The synthetic runtime is not a
trajectory, production family run or physical authorization. Missing
health still closes holdout; it does not prove nonhyperbolicity.

## Prospective admission, initial-health, cone and continuity instruments

Later outcome-blind slices now live in the following focused owners:

- [source admission](fgc-sgb1-ctl1-admission.md) implements a declared-box
  parametric Krawczyk uniqueness certificate for the affine six-row source;
  it uses no PROTO4 residual floor and does not set production admission;
- [initial health](fgc-sgb1-ctl1-initial-health.md) encloses the actual radial
  H/M constraint ODE by a fixed-resource interval-Picard construction and
  evaluates continuous compactness directly on each graph box;
- [cone/symmetrizer](fgc-sgb1-ctl1-cone.md) owns exact interval coefficient
  tensors and an all-direction reference/perturbation instrument without
  importing HYP2 or FGC-QR pass bits;
- [constraint continuity](fgc-sgb1-ctl1-continuity.md) rechecks the scalar,
  affine-source, H/M and hat-wave identities on SGB-L data;
- [no-trap refinement 1](fgc-sgb1-ctl1-trap-refinement.md) and
  [refinement 2](fgc-sgb1-ctl1-trap-refinement2.md) preserve the exact
  depth-8/10/12 and depth-14/16/18 product-box resource ladders;
- [Taylor no-trap enclosure](fgc-sgb1-ctl1-trap-taylor.md) preserves a
  separately frozen order-4 interval-jet route, also without a nominal pass;
- [matched-family principal feeder](fgc-sgb1-ctl1-family-principal.md)
  maps authenticated family two-jets/source witnesses into interval
  orthonormal Riemann/Hess(phi) boxes and names missing whole-cell owners;
- [whole-cell source admission](fgc-sgb1-ctl1-cell-admission.md) owns the
  nonuniform product-box Krawczyk contract and retains its family-cell
  wrapping nonpass;
- [point-local symmetrizer](fgc-sgb1-ctl1-local-symmetrizer.md) proves the
  exact ESF 24-mode owner while retaining the nonflat missing-mode result;
- [trap barrier](fgc-sgb1-ctl1-trap-barrier.md) binds the exact failure of a
  whole-domain `D=1-C` Nagumo proof without claiming the orbit traps;
- [SOL1 orbit continuation](fgc-sgb1-ctl1-sol1.md) freezes a 16-step,
  `rho=1/8`, depth-0 Picard--Lindelof tube from the authenticated `(C,k)`
  prefix; nominal `A_chi=3` is interval-inconclusive because the prefix
  endpoint is wider than the first tube;
- [SOL1-FRZ1](fgc-sgb1-ctl1-sol1-frz1.md) compactly freezes that exact policy,
  nominal obstruction, and control-only branch;
- [SOL1-PREF1](fgc-sgb1-ctl1-sol1-pref1.md) independently authenticates the
  freeze and reconstructs the unchanged nominal/control outcomes without
  promoting aggregate branch health;
- `src/recursive_horizons/fgc/sgb1_ctl1_controls.py` owns weak-field,
  constant-coupling/topological, Schwarzschild and established
  decoupling-limit controls.

These remain prospective instruments. The nominal `A_chi=3` first-order
Picard compactness enclosure remains inconclusive because interval wrapping
makes its direct upper bound slightly exceed one; this is not a trapped-sphere
result. Two prospectively finite resource ladders through depth 18 and a
separately frozen order-4 Taylor route preserve gap-free prefixes and never
tile the full support with `C<1`; no informal deeper/order retry is licensed.
The proved cone is still the curvature-free
SGB-L/ESF reference, not the nonflat matched-family box. The family feeder
passes the exact buffer; whole-cell acceleration boxes and the local nonflat
symmetrizer remain typed inconclusive. The exact Nagumo witness rejects only
the trajectory-free invariant-barrier proof, not SGB-L or the actual orbit.
SOL1 preserves the same family and refuses before its first continuation
step; it does not retune the tube or turn the `A_chi=1/8` control into nominal
data. Unique remaining-support continuation therefore requires a separately
frozen matched-data redesign with a tighter authenticated endpoint.
Constraint continuity is a conditional boundary-free propagation
theorem; smooth existence, compatible data, regular center, boundary map and
IBVP remain explicit premises.
The evolving-center contract remains unexecuted until a later branch-owned
trajectory exists. Accordingly `SGBL_branch_owned_and_healthy`, holdout,
execution and physical claims remain false. SOL1-FRZ1/PREF1 seal only the
scoped continuation inconclusive.
