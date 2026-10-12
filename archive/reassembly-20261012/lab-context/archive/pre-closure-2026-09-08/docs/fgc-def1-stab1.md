# FGC-1-DEF1-STAB1: bounded outcome-blind error map and exact mass algebra

**FGC-1-DEF1-STAB1** is the Wave 0 constraint-to-Raychaudhuri stability/error
map. It owns the saved `RaychaudhuriErrorBudget` and `ActivationAssessment`
interfaces, fail-closed componentwise `Q`-error assembly, exact covariant
Misner–Sharp algebra, typed `mass_flux_inconsistency` mapping, activation
arithmetic with an explicit case reference, and independent trappedness,
complete-`Q`, affine, and control-margin helpers.

It is an algebra, input-contract, and outcome-blind provider-conversion slice.
The prospective qualification module
`src/recursive_horizons/fgc/def1_stab1_qualification.py` adds conversion-
ownership contracts required before the later FRZ1/PREF1 freeze: the exact
RED1/constraint/source-to-MHG2 six-row identity, the FO1-to-ADM/26-input
geometry-slot map, the IMP1 18-channel-to-`Q` Lipschitz debit, honest owner
provenance, and a fourteen-component route coverage matrix. Those contracts
do not freeze PREF1, do not read a candidate or control trajectory, do
not certify a global PDE error, and do not emit a complete physical DEF1 gate.
The compact [FRZ1 freeze](fgc-def1-stab1-frz1.md) binds them without setting
`DEF1_error_map_passed`. Independent
[STAB1-PREF1](fgc-def1-stab1-pref1.md) now binds candidate-blind pre-holdout
map readiness only.

The separate [geometric sensitivity instrument](fgc-def1-stab1-geometry.md)
computes whole-box Q sensitivities from complete metric two-jets and a radial
tangent. The provider module converts each of the fourteen source-specific
errors into those 26-input radii and/or a justified additive Q debit, then
evaluates the Jacobian on the single joint box. Input provenance stays
visible. Independent PREF1 now binds candidate-blind map readiness; the
nine trajectory DEF1 booleans and actual error radii remain later.

The saved 27 August interface is the frozen field list. Current campaign
routing stays in [PLAN.md](../PLAN.md) and the
[research roadmap](research-roadmap.md). This document does not replace
either.

## Frozen interfaces

```text
RaychaudhuriErrorBudget
  spatial_temporal
  physical_constraint
  gauge_constraint
  reduction_constraint
  initial_data
  nonlinear_source
  affine
  nullness
  trajectory_alignment
  interpolation
  extraction
  boundary
  conservation
  arithmetic
  total_upper_bound

ActivationAssessment
  raw_factor
  lower_bound
  matched_control_upper_bound
  threshold_passed
  control_dominance_passed
```

Every stored error-budget field is a finite nonnegative float obtained from
exact rational arithmetic and outward upper conversion. Direct construction
is a checked contract: `total_upper_bound` must dominate the exact
`Fraction` sum of the stored binary64 component values, not merely each
component and not a floating accumulation that can round down. Missing,
negative, nonfinite, or non-dominating totals refuse.
Convenience integer inputs must convert exactly to binary64 and are stored
as floats; silently rounding a supplied bound is refused. An activation lower
bound may not exceed its raw factor.

## Componentwise `Q` error

Complete `Q` is the spherical twist-free Raychaudhuri right-hand side already
defined by DEF0:

```text
Q = -theta^2/2 - sigma_ab sigma^ab - R_ab k^a k^b.
```

This slice does not estimate `E_Q` from a measured positive `Q`. Each frozen
component requires an explicit sensitivity enclosure and an explicit
nonnegative input-error enclosure. The component contribution is the exact
product of the sensitivity absolute-upper bound and the input-error upper
bound. The total is the exact sum, then converted outward to binary64.

Missing, extra, duplicate-owned, negative, nonfinite, boolean, or malformed
premises refuse. They are not replaced by zero. An explicit supplied zero
enclosure is allowed and is distinct from a missing premise. Overlap
rejection means duplicate ownership of the same frozen component. Correlated
or common physical sources across different components remain allowed under
the conservative sum; the sum does not assume independence.

Premise status labels are caller declarations, not independently
authenticated enclosure certificates:

```text
proven_enclosure
conditional_premise
```

Allowed sources are `supplied_certified_enclosure`, `exact_algebra`,
`richardson`, `imp1_admission_debit`, and `declared_conditional_enclosure`.
Those names record how the caller classified the premise. This slice does
not authenticate that a `supplied_certified_enclosure` is a true continuum
certificate. Exact-control zero helpers are algebraic test inputs, not
generic certified scientific evidence.

Richardson may appear only on `spatial_temporal`. IMP1 admission debit may
appear only on `spatial_temporal` or `arithmetic`. Both are conditional.
Neither is a global PDE-error theorem, and the assembly record keeps

```text
all_components_declared_proven_enclosures  (derived from status declarations)
global_pde_error_certified = false
used_measured_q = false
richardson_treated_as_global_pde_error = false
imp1_admission_debit_treated_as_global_pde_error = false
```

even when those conditional debits are present in the sum. Assembled records
validate exact component types, order, nonnegativity, exact sum, matching
stored bounds, ordered status/source inventories, and derived declaration
flags. Scope flags are strict bools. A later owner may tighten a declared
enclosure; it may not relabel Richardson or the IMP1 admission debit as a
universal solution-error bound. This algebra slice does not set
`DEF1_error_map_passed`. The provider assembly also leaves that flag false
and records the remaining conversions.

## Activation and GR-0 phi policy

Activation uses one explicit positive case reference `S_ref`:

```text
A = max_measure |phi| / S_ref
lower = A - E_A
```

`S_ref` must be strictly positive. `0/0` is refused. The two GR-0 policies are
different and must be declared:

```text
canonical_gr0_phi_zero
historical_planted_calibration_seed
```

Canonical GR-0 has `phi = 0`, so the GR-0 numerator is zero and the
denominator remains `S_ref`. The historical planted calibration seed (for
example `1/131072` in earlier GR-0 campaigns) is not substituted silently,
is not used as a hidden denominator, and cannot be attached to the canonical
policy. Historical policy requires that seed to be supplied as a positive
declared value; the ratio still uses `S_ref`.

Threshold and dominance remain independent:

```text
A - E_A >= 2
A - E_A > max(A_GR0 + E_A,GR0, A_SGBL + E_A,SGBL)
```

Exact rationals never weaken `>= 2` or strict dominance. The frozen
five-field DTO stores outward floats and decides `threshold_passed` and
`control_dominance_passed` from those saved bounds, so a True
strict-comparison flag cannot survive after the stored interval meets or
overlaps. Direct construction and `dataclasses.replace` attacks are
checked against the same saved-bound rule.

## Exact Misner–Sharp algebra

On the Lorentzian 2D base,

```text
m = R (1 - h^{ab} R_a R_b) / 2
d_a m = (R^2 / 2) (G_a^b - delta_a^b G_c^c) R_b
```

with the trace over the 2D base only. If `E = F G - T_eff` and `F > 0`, the
equation defect is

```text
(R^2 / (2 F)) (E_a^b - delta_a^b E_c^c) R_b.
```

Hamiltonian and momentum residuals alone are insufficient for the full
spacetime mass-flux vector. In an orthonormal frame

```text
P[E]_n = -E_ss N(R) + E_ns S(R)
P[E]_s = -E_ns N(R) + E_nn S(R)
```

`E_ss` is required for the normal/time component, not generally for the
radial component. Public APIs that accept an inverse metric require a
symmetric Lorentzian matrix with determinant `< 0`. That symmetry test is
not applied to mixed tensors. Incomplete, non-Lorentzian, nonsymmetric
inverse, nonpositive-`F`, or nonpositive-`R` inputs refuse.

Independent exact controls used here:

- spherical Minkowski: `m = 0`, `d_a m = 0`, `E = 0`;
- contracting flat de Sitter: `Lambda = 3`, `m = Lambda R^3 / 6 = 4` at the
  DEF0 chart point, `d_a m = (Lambda/2) R^2 R_a`, and dropping the 2D trace
  adjustment fails that identity;
- injected `E` with a compensating `T_eff` restores the algebraic
  `dm = flux + defect` split. That checker is an `E = F G - T_eff` algebra
  identity on supplied mixed tensors. It is not an independent geometric or
  trajectory conservation check.

The ledger test is

```text
|Delta m_MS - F_integrated| <= E_M.
```

Equality passes. Exceeding the enclosure is a typed veto. The protocol token
remains `mass_flux_inconsistency`. The interpretive class remains
`stopped_mass_flux_inconsistency`. This slice does not invent a new stop
token and does not edit the PROTO19 stop list.

## Independent helpers and nine booleans

Helpers preserve the existing numerical gates and do not collapse them:

```text
theta_+ < -4 E_theta+
theta_- < -4 E_theta-
min_I Q > 4 E_Q > 0
eight affine samples
Delta lambda / L0 >= 1/64
|direct - assembled| <= 1/100000000
affine residual <= 1/10000000000
three nested resolutions
two methods
p >= 3/2
```

Trappedness and complete-`Q` inequalities are strict. The activation
threshold is closed at 2. Control dominance is strict. Affine samples must be
strictly increasing; overlapping or unordered samples refuse.

The nine DEF1 booleans stay independently visible:

```text
resolved_activation
control_dominance
resolved_trapped_interval
resolved_complete_defocusing
direct_Raychaudhuri_agreement
all_health_constraints_scales_valid
mass_flux_ledger_valid
three_resolution_convergence
two_method_agreement
```

A resolution-count or method-count helper is not three-resolution
convergence or two-method agreement. Missing health or convergence cannot be
inferred from a passed margin. An incomplete boolean mapping refuses.

Exact DEF0 geometry remains a control, not a DEF1 result: Minkowski is
normal and focusing; the contracting de Sitter chart is trapped with
`R_kk = 0` and negative complete `Q`; the synthetic positive-`Q` point is
marginal, not trapped.

## Proven scope

This slice proves and tests only:

- frozen dataclass field contracts and finite nonnegative bound storage;
- fail-closed assembly of a linearized componentwise `Q` upper bound from
  supplied enclosures, independent of measured `Q`;
- distinction between declared proven/conditional labels and authenticated
  enclosure certificates;
- refusal to treat Richardson or IMP1 admission debit as global PDE error;
- exact 2D Misner–Sharp mass, gradient, and algebraic `E = F G - T_eff`
  split on supplied tensors;
- typed mapping onto the existing mass-flux protocol token;
- activation arithmetic with explicit `S_ref` and a declared GR-0 phi policy;
- independent margin helpers and an explicit nine-boolean record;
- exact Minkowski, trace-sensitive de Sitter, synthetic-`Q`, injected-defect,
  invalid-input, absent-premise, and inequality-edge controls;
- fourteen explicit source-to-geometry or additive-`Q` conversions with
  declaration, provenance, context, and unit;
- one joint-box Lipschitz evaluation and componentwise
  `sum_i L_i epsilon_{c,i} + a_c` assembly;
- refusal of missing C-to-jet, Gronwall, interpolation/extraction derivative
  bounds, H/M-only or C-only residuals, residual-times-step affine bounds,
  and causality-only zero boundary debit;
- preservation of the mass-flux veto and of `DEF1_error_map_passed = false`.

It does not prove a continuum PDE truncation theorem, a COL1 trajectory
enclosure, a physical DEF1 classification, ROB1 robustness, or retained-EFT
validity.

## Provider conversions and joint-box assembly

`src/recursive_horizons/fgc/def1_stab1_providers.py` owns the fourteen
explicit conversions. Focused tests are
`tests/test_fgc_def1_stab1_providers.py`. The frozen algebra and geometry
instruments remain `def1_stab1.py` and `def1_geometry_error.py`. Conversion
ownership before FRZ1/PREF1 is
`src/recursive_horizons/fgc/def1_stab1_qualification.py`, tested by
`tests/test_fgc_def1_stab1_qualification.py`. That module does not promote
the catalog, freeze config/result/reproducer artifacts, or set any DEF1
boolean.

Each provider returns a `ProviderRecord` with:

- complete 26-input nonnegative radii in `INPUT_NAMES` order, including
  explicit conversion zeros on unused slots;
- a nonnegative additive complete-`Q` debit;
- `status`, `source`, a shared evaluation `context` and `unit`, component
  `provenance`, and a canonical `conversion` string;
- component-specific guards where zero additional debit is claimed.

Duplicate radius names, aliases, and out-of-order inventories refuse rather
than collapsing through a dict. Missing context, unit, or provenance
refuses. Float, bool, and nonfinite aliases refuse. An explicit supplied
zero is distinct from a missing premise. Registered zero-only
Richardson/IMP1 stubs refuse. `context` and `unit` are a declared
evaluation-context/normalization identity shared by every assembled
record; provenance may vary. Those strings remain declarations, not
authenticated production provenance.

### Joint box, not fourteen small boxes

Let `epsilon_{c,i}` be provider `c`'s radius on geometry input `i`, and
`a_c` its additive debit. The assembly constructs one box

```text
epsilon_i = sum_c epsilon_{c,i}
Box = [x_i - epsilon_i, x_i + epsilon_i]
L_i = sup_Box |partial Q / partial x_i|
E_c = sum_i L_i * epsilon_{c,i} + a_c
E_Q = sum_c E_c
```

`L_i` is the existing whole-box interval Jacobian from
[the geometry instrument](fgc-def1-stab1-geometry.md). It is evaluated once
on the joint box. Independent smaller boxes underbound mixed partials: a
Lipschitz constant that depends on another component's coordinates is larger
on the joint box than on a component's own slice. The mixed-box control with
Minkowski nominal, `R.dr` radius `1/8` and `k.r` radius `1/4` has exact
`|ΔQ| = 1001/2048`, separate-box sum `29/64`, and joint bound `315/512`.
An `R.dr`-only radius `1/8` has honest bound `9/64`.

The checked DTO stores the nominal point, shared evaluation context and
unit, joint radii, the reconstructed box whose widths are those radii,
Lipschitz enclosures, per-component `radii_term + additive`, and the frozen
`AssembledQErrorBudget`. Outer source/status labels are tied to that inner
inventory. `replace` of `joint_box`, Jacobian, contributions, or inner
budget is rederived from the stored nominal and records; cached flags and
attacker-supplied Jacobians are not trusted. The factory derives those
fields once; `replace` and direct construction rederive. `QComponentPremise.input_error`
on that budget is the already-converted Q contribution, with sensitivity
`1`.

The bound does not use measured `Q` or its sign. Correlated sources remain
allowed under the conservative sum.

### Conversion ownership and status

| Component | Supported nonzero conversion | Refuses |
|---|---|---|
| `spatial_temporal` | declared complete radii and/or additive debit; Richardson or IMP1 stay `conditional_premise` | proven Richardson/IMP1; zero-only Richardson/IMP1 stubs; global PDE labels |
| `physical_constraint` | complete 4×4 componentwise `sum \|E_ab\| \|k^a\| \|k^b\|/F`, enclosed over the joint `k` box | H/M-only scalars; missing `F`; 2×2 projectors; `k.t<=0`; nonzero angular `k`; foreign tangent |
| `gauge_constraint` | complete gauge-extension componentwise residual debit enclosed over the joint `k` box, or certified C-to-jet inverse bound times the full `C` vector | `C` alone; missing inverse and missing extension; foreign tangent |
| `reduction_constraint` | explicit eight metric first-derivative discrepancies times a supplied operator bound | missing slot; missing operator; nonzero discrepancy with bound `0` |
| `initial_data` | explicit complete 26-input radii | missing/extra slots |
| `nonlinear_source` | certified inverse-`J` bound times the full six-residual, assigned to metric `.dtt` radii | partial residual; missing inverse bound; nonzero residual with bound `0` |
| `affine` | transport defect times a supplied finite-interval Gronwall factor, assigned to `k` radii | residual times step; missing factor; nonzero defect with factor `0` |
| `nullness` | `g(k,k)` residual times a supplied null-to-tangent bound | missing bound; nonzero residual with bound `0` |
| `trajectory_alignment` | explicit complete 26-input radii | missing/extra slots |
| `interpolation` | Lipschitz remainder: enclosed first derivatives times sample spacing, plus explicit remainder radii for uncovered jet slots | missing derivative bounds; second-order remainder without second derivatives |
| `extraction` | `dQ/dλ` bound times the maximum unsampled affine gap | missing derivative bound; one sample; positive sample values as a continuous-`Q` certificate |
| `boundary` | supplied additional debit, or explicit zero with no-influence, coverage, and an independent numerical-boundary guard | physical causality alone as exact zero numerical boundary error |
| `conservation` | supplied additional debit inside a valid mass-flux ledger, or explicit zero with coverage | invalid ledger (`mass_flux_inconsistency`); zero additional without coverage |
| `arithmetic` | explicit additive debit, or nonzero IMP1 local admission debit | IMP1 zero-only stubs; IMP1 treated as a PDE theorem |

The physical and gauge additive debits are componentwise magnitude
enclosures of a residual, not a global solution-error theorem and not
`|E_ab k^a k^b|/F`. The signed contraction `E_ab k^a k^b` remains an exact
diagnostic; exact cancellation there is legitimate and is not an
implementation defect. Entrywise error bounds cannot cancel. The supplied
tangent must match the nominal `(k.t, k.r)` and be future-directed and
radial; the debit is then enclosed over the joint `k` box with a positive
`F` lower bound. A point contraction is not a whole-box bound.
Hamiltonian/momentum projections or the
gauge vector `C` do not determine that debit. Source inversion uses
the independent six-equation residual order from the spherical reduction,
not a floating Jacobian. Affine transport requires an explicit finite
interval together with a validated Gronwall factor; the factor is not
replaced by `defect * Δλ`. Interpolation may use a first-derivative
Lipschitz bound; it may not advertise a second-order remainder without
nonempty enclosed second derivatives covering the claimed slots. Empty
second-derivative maps with remainder zeros are not a second-order proof.
Extraction covers the open gaps between
strictly increasing samples; positivity of sampled `Q` is ignored.

Zero *additional* boundary or conservation debit is a guarded conversion,
not an absent provider. The existing mass-flux veto remains: an invalid
ledger refuses both the conservation provider and joint assembly, and still
maps to `mass_flux_inconsistency`. Frozen thresholds and the nine final
DEF1 booleans are unchanged and are not inferred from this map. If an
independent `Def1BooleanRecord` is supplied, it is stored only as a
caller-owned record.

Richardson may appear only on `spatial_temporal`. IMP1 admission debit may
appear only on `spatial_temporal` or `arithmetic`. The assembled DTO keeps

```text
def1_error_map_passed = false
used_measured_q = false
global_pde_error_certified = false
richardson_treated_as_global_pde_error = false
imp1_admission_debit_treated_as_global_pde_error = false
```

and copies the concrete remaining-work list. A compact certificate is not
emitted.

## Conversion-ownership qualification contracts

The qualification module is candidate-blind. It reconstructs live owner
constants rather than inventing maps.

1. **RED1/MHG2 six-row identity.** `INDEPENDENT_EQUATION_ORDER` maps onto
   `MHG2_FULL_EQUATION_ORDER` by the exact pairs
   `(metric_tt, metric_tt_mhg)`, `(metric_tr, metric_tr_mhg)`,
   `(metric_rr, metric_rr_mhg)`, `(metric_theta_theta, metric_theta_theta_mhg)`,
   `(scalar_phi, scalar_phi)`, `(scalar_chi, scalar_chi)`. MHG1, REF1,
   MHG3-IMP1, FO1 physical rows, SGB1 source, and SRC1-NL1 already publish
   that MHG2 order. CON4 unredefined residual *names* use the RED1
   vocabulary. Hamiltonian/momentum projections, gauge `H`/`C` vectors, FO1
   kinematic reduction rows, and regular-centre equation names are typed
   non-identities. Omissions, dict aliases, and reordering refuse.
2. **FO1-to-ADM/26-input geometry slots.** FO1 group-to-jet identity is
   `(u,p,q,p_t,p_r,q_r) -> (value,dt,dr,dtt,dtr,drr)`. FO1-RC1 remains a
   BASE-metric inventory. ADM-to-BASE is owned by
   `state_from_generalized_adm_pg_fixture`. The inverse BASE-to-ADM two-jet
   map is now owned as exact `Jet2` calculus: given complete
   `h_tt,h_tr,h_rr` and caller-supplied positive roots that square exactly
   to the radicands,
   `lambda=sqrt(h_rr)` via `Jet2.compose` with `q'=1/(2 lambda)` and
   `q''=-1/(4 lambda^3)`, `shift=h_tr/h_rr`, and
   `alpha=sqrt(h_tr^2/h_rr-h_tt)` with the same compose coefficients.
   The inverse requires a Lorentzian positive branch and proves exact
   roundtrip of all three BASE metric two-jets through the ADM fixture.
   `areal_radius` is a documented name identity for geometry `R`. Together
   with a caller-owned tangent, that inverse packs the complete 26-input
   inventory. First-order `(u,p,q)` state still cannot fill second jets.
   No 26-input slot has a mathematically justified implicit zero; missing
   roots, missing jet slots, or missing tangent refuse rather than
   substitute zero. `phi`/`chi` are not geometry inputs. IMP1 `v` is a
   documented name identity for `shift`, not a silent alias in either
   inventory. ADM/source-named FO1 packing still identity-fills the 24
   metric two-jet slots without using the BASE inverse.
3. **IMP1 18-channel-to-`Q`.** The channel order is the live TDG6/IMP1
   `u/p/q` of `(alpha, v, lambda, R, phi, chi)`. Numeric debit is the exact
   sum of caller-supplied nonnegative per-channel Lipschitz times debit.
   A universal numeric factor is refused. The record stays
   `conditional_premise` / `imp1_admission_debit` and cannot set
   `global_pde_error_certified`, `DEF1_error_map_passed`, or measured-`Q`
   flags.
4. **Honest provenance.** Records bind owner artifact id, repo-relative
   path, live file SHA-256, and conversion-identity digest. They
   authenticate those owner bytes. They cannot claim future COL1 values.
5. **Fourteen-component coverage matrix.** Every current provider route has
   an executable positive control and an injected-failure refusal.
   Trajectory values and the nine DEF1 booleans remain unevaluated.

C238 remains a partial-instrument claim. FRZ1 freezes the candidate-blind
contract. Independent [STAB1-PREF1](fgc-def1-stab1-pref1.md) now binds
candidate-blind pre-holdout map readiness (`DEF1_error_map_passed` with
`map_readiness_only`). Real COL1 input-error radii and the nine trajectory
DEF1 booleans stay later application. A universal PDE theorem is not claimed.

`src/recursive_horizons/fgc/def1_stab1_freeze_contract.py` aggregates these
live conversion owners into the deterministic candidate-blind payload bound
by [FGC-1-DEF1-STAB1-FRZ1](fgc-def1-stab1-frz1.md). The compact freeze
directly exercises the conversion, coverage, algebra, and margin
controls and keeps `DEF1_error_map_passed`, trajectory, ROB1, holdout, and
physical flags false. Independent PREF1 reconstructs those controls from the
low-level owners rather than importing the freeze contract.

## Remaining gate-closing work

Provider qualification is outcome-blind and belongs before holdout. The
conversion-ownership contracts above are the map layer. Actual COL1 data
binding is later. `PLAN.md` requires `DEF1_error_map_passed` before holdout,
together with `GR0_case_eligible` and `SGBL_branch_owned_and_healthy`.
Routing that gate to Wave 5 after COL1 would create a dependency cycle.
Independent [STAB1-PREF1](fgc-def1-stab1-pref1.md) now sets
`DEF1_error_map_passed=true` with `map_readiness_only=true`. That is not a
trajectory result and does not add a universal PDE theorem. FRZ1 remains a
compact nonexecuting instrument freeze.

Distinguish two later jobs:

1. Remaining outcome-blind STAB1 map qualification before holdout: keep
   Richardson/IMP1 as honest conditional numerical estimators under frozen
   hypotheses, supply C-to-jet or complete extension residuals from the
   gauge owner when available, and keep every conversion explicit. A new
   universal nonlinear PDE theorem is not a prerequisite. Do not read
   holdout or candidate trajectories.
2. Later application/binding of that qualified map to COL1 trajectories,
   then Wave 5 `FGC-1-DEF1-PREF1` classification. That application is not
   the pre-holdout map gate and must not create a post-COL1 dependency
   cycle.

| Component or boolean | Remaining owner |
|---|---|
| `spatial_temporal` declared debit | still open: Richardson and IMP1 are implemented as conditional local estimators under frozen hypotheses; a universal PDE theorem is not required to qualify this map |
| `physical_constraint` COL1 residual | conversion exists; later binding of the complete `E_ab` from the constraint monitor |
| `gauge_constraint` C-to-jet inverse | conversion exists only if the caller supplies the inverse or the complete extension residual; the inverse is not derived here |
| `reduction_constraint` COL1 discrepancies | conversion exists; later binding of FO1 first-derivative residuals and operator bounds |
| `initial_data` family residual | conversion exists; later binding of explicit complete ID radii |
| `nonlinear_source` certified inverse-`J` | conversion exists; later binding of the full residual and a certified inverse bound |
| `affine` Gronwall factor | conversion exists; the factor remains a supplied finite-interval premise, not a derived transport theorem |
| `nullness` | conversion exists; later binding of `g(k,k)` and a tangent bound |
| `trajectory_alignment`, `interpolation`, `extraction` | conversions exist as supplied-radii / Lipschitz / gap-remainder maps; COL1 trajectory reading is later |
| `boundary` numerical debit | conversion exists; causality is not a zero-error theorem; later binding of coverage and independent guards |
| `conservation` `Delta m` and `F_integrated` | ledger veto and additional-debit conversion exist; later evolution binding using this algebra |
| `arithmetic` binary64 path debit | IMP1 local admission debit is implemented as conditional arithmetic, not a PDE bound |
| health / scale / constraint validity | HLT; later COL1 application |
| three-resolution convergence and two-method agreement | later COL1 / DEF1-PREF1 application after trajectories exist |
| `DEF1_error_map_passed` | independent STAB1-PREF1 now sets this true with `map_readiness_only`; it is not a trajectory result and does not evaluate the nine DEF1 booleans |

Reproduce the algebra slice, geometry instrument, provider layer, and
conversion-ownership contracts with:

```bash
/opt/homebrew/bin/python3.14 -I -B tests/test_fgc_def1_stab1.py
/opt/homebrew/bin/python3.14 -I -B tests/test_fgc_def1_geometry_error.py
/opt/homebrew/bin/python3.14 -I -B tests/test_fgc_def1_stab1_providers.py
/opt/homebrew/bin/python3.14 -I -B tests/test_fgc_def1_stab1_qualification.py
/opt/homebrew/bin/python3.14 -I -B tests/test_fgc_def1_stab1_freeze_contract.py
/opt/homebrew/bin/python3.14 -I -B tests/test_fgc_def1_stab1_frz1_certificate.py
/opt/homebrew/bin/python3.14 -I -B tests/test_fgc_def1_stab1_pref1_binder.py
```
