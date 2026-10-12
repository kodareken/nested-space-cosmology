# FGC-1-HLT3-MON3: PROTO5 GR-0 runtime authorization

## Decision

`FGC-1-HLT3-MON3` closes the pre-trajectory runtime-composition boundary
left by `FGC-1-PRO5-FRZ1`. Its two positive decisions are deliberately narrow:

```text
PROTO5_successor_runtime_compositor_implemented = true
PROTO5_fresh_GR0_dynamic_calibration_authorized = true
```

The artifact authorizes only the fresh GR-0 amplitude-calibration campaign
defined by `FGC-1-CAL2-RUN1-PLAN`. It does not complete that campaign,
select an amplitude, resolve the future holdout manifest, run SGB-L or FGC-QR,
or answer the mechanism question.

The authorization is evaluated before a physical trajectory is advanced.
Its physical inputs are the already declared `t=0` ID2 slices, reconstructed
from the immutable PROTO5 checkpoint. Its time-dependent controls use only a
synthetic zero right-hand side and deliberate stop injections.

## Immutable lineage and fresh namespace

The reproducer reads PROTO5/PRO5-FRZ1, CAL1, ID2, HLT2, and the CAL1 runtime
adapter directly from immutable Git commit
`cf48cf729ff1c3a6effce6ecd3fdbeaf8d9645bf`. Every blob is checked against the
hash frozen in the HLT3 configuration. This prevents a later successor file
from silently changing the premises that authorize the run.

The two PROTO5 output roots are inspected without being created:

```text
runs/fgc-2-sf1/proto5/calibration
runs/fgc-2-sf1/proto5/holdout
```

Both must be absent or empty. The calibration runner must refuse an unrelated
pre-existing campaign, refuse overwrite, and bind any legitimate resume to
the frozen plan and campaign-manifest hashes. The holdout path remains closed.

## A real GR-0 monitor, not padded FGC-QR data

HLT2 proved that the historical 26-stop HLT1 record has three ownership
classes. GR-0 owns exactly nine universal runtime stops:

```text
nonpositive_lapse
nonpositive_radial_metric
nonpositive_areal_radius_away_from_center
hat_cone_not_Lorentzian
newton_residual_limit
newton_iteration_limit
newton_residual_not_monotonic
kinetic_condition_limit
boundary_causal_buffer
```

HLT3 implements a branch-specific transaction containing exactly those nine.
It does not populate the eleven FGC-QR-only fields with synthetic safe values.
The GR-0 linear acceleration solve exposes its residual, iterative-refinement
count and monotonicity, and kinetic condition number; those actual diagnostics
feed the historically named solver stops.

The hat-cone test is calculated from the actual ADM stage metric. For hat
normal factor `q=9`, its radial-base conditions are monitored through

```text
-hat_g^tt = q / alpha^2 > 0,
-det(hat_g_base) = q / (alpha^2 lambda^2) > 0,
1 / lambda^2 > 0.
```

Positive lapse, radial metric, and areal radius away from the regular centre
remain independent stops.

Every internal Runge--Kutta stage and the explicit candidate endpoint enters
one atomic transaction. A scientific failure latches the first premise and
advances neither the accepted state nor the causal ledger. BND2's coordinate-
speed envelope is debited over the proposed interval before acceptance.

A Courant violation has a different meaning. It rejects the numerical
proposal and requires a smaller retry, but it does not become evidence against
GR-0 or FGC. The synthetic controls demonstrate healthy acceptance, late-stage
rollback, immutable first-failure ownership, boundary rollback, and this
retryable CFL distinction for the two frozen integrator implementations.

## Common-event numerical admission

The six numerical decisions replaced by HLT2 are not evaluated one grid at a
time. At every predeclared common output event, each method supplies all three
nested resolutions.

The constraint compositor evaluates the Hamiltonian, radial momentum, two
metric-defined gauge constraints, and six reduction constraints directly from
the unredefined GR-0 residual. PROTO5's method-owned raw guards apply:

```text
RK4:    coarse < 1/50, fine < 1/1000
SSPRK3: coarse < 1/10, fine < 1/200
```

All three values must decrease monotonically and the finest adjacent pair
must have order at least `3/2`. The `4096 eps64` enclosure may classify a
component as numerically zero for order bookkeeping only. Raw norms still
decide both magnitude guards, and the enclosure is never subtracted.

Spatial spectra use the six HLT2 Minkowski-reference deviations, regular-
centre limits, proper radial distance, deterministic resampling, and a
two-edge compact window. Temporal spectra use only causal-past samples along
predeclared normal-flow tracers. Each tracer's generally nonuniform proper
time is resampled onto a uniform FFT grid, and the reverse-interpolation
mismatch is retained as a nonnegative public error component. At least 64
past samples are required.

The certificate reconstructs both eligible amplitudes (`5/2`, then `3`) for
both methods and all three grids. It verifies all twelve continuum state
hashes against ID2, preserves `u` and `p` bitwise, initializes only
`q=D_h u`, hashes every projected state and expanded run configuration, and
passes all four amplitude/method initial common events. These are initial
premises, not collapse outcomes.

## Trapped-sphere selection rule

The calibration observable keeps the fixed physical-null orientation. On
nodes shared by all three grids it forms

```text
T_h = max_r min(-theta_plus(r), -theta_minus(r)).
```

Thus `T_h>0` means that both expansions are negative at the same sphere;
compactness alone is never substituted for the two signs. The conservative
sign error adds, without cancellation:

```text
primary finest-pair Richardson estimate
+ comparator finest-pair Richardson estimate
+ fine-grid cross-method disagreement
+ binary64 roundoff enclosure.
```

Constraint and spectral admission are hard vetoes. They are not multiplied by
an arbitrary unit-conversion factor and called an expansion error. A common
event counts only when both methods pass their full event admission and the
fine sign margin exceeds four times the additive error. Eight consecutive
common events are required. The first amplitude satisfying this frozen rule
ends calibration; otherwise the next declared amplitude is tried.

This is a calibration sign contract. It does not supply the still-missing
constraint-to-Raychaudhuri stability map and cannot be reused as DEF1 error
authorization.

## Frozen numerical campaign

`FGC-1-CAL2-RUN1-PLAN` fixes the two amplitudes, method order, three grids,
physical pulse, cutoff, domain, measurement window, CFL ceiling, dissipation,
common-event and checkpoint cadence, proper-time tracers, stop thresholds,
selection rule, and provenance behavior. Its twelve expanded input hashes and
one campaign identifier are serialized by HLT3 before any trajectory exists.

The runner may reduce a proposed step to hit a common event exactly or to
satisfy the frozen CFL ceiling. It may not change a physical parameter or a
classification threshold in response to an outcome.

The authorized runner is itself included in HLT3's implementation hash map.
Before it creates the fresh namespace it reproduces HLT3, verifies every
authorized implementation hash, and requires a clean tracked worktree. It
then creates the campaign manifest exclusively; overwrite or an unrelated
resume is refused.

The atomic checkpoint owns the six member states, their monitor and causal
ledgers, all already completed amplitude records, and the exact event-log
bytes. Checkpoints are installed at the frozen `1/4` coordinate-time cadence.
If interruption lands after a line is flushed but before its checkpoint is
installed, resume preserves that uncommitted tail in a hash-named sidecar and
restores the checkpoint-owned history. A divergent history is rejected rather
than guessed through. An interrupt while the resolution members are between
common events never writes a mixed-time checkpoint.

## Reproduction

```bash
python3 scripts/reproduce_fgc_hlt3_mon3.py \
  --output results/fgc-1-hlt3-mon3.json
```

This command must leave both fresh run namespaces untouched.

After the authorized calibration has launched, use `--check`; the verifier
then preserves the historical fresh-namespace observation while requiring the
live manifest to bind the exact HLT3 result, campaign identifier, input
manifest, implementation hashes, and an ancestor Git checkpoint. It still
rejects any nonempty holdout namespace.

## Nonclaims

HLT3 does not show that collapse occurs, that a trapped sphere forms, that a
regulator activates, or that metric-null defocusing occurs. It supplies no
SGB-L source or health owner, no resolved holdout manifest, no FGC-QR
trajectory, no retained-EFT authorization, and no physical transition. It
does not establish singularity resolution, a child domain, a dark-sector
mechanism, varying locally measured light speed, or a rejection of the
general gradient programme.
