# FGC-1-HLT9-MON9: PROTO11 runtime preflight and GR-0 authorization

`FGC-1-HLT9-MON9` is the pre-trajectory runtime gate for the frozen PROTO11
reference-state differentiation map. Its canonical machine record is
`results/fgc-1-hlt9-mon9.json`.

The two positive decisions are deliberately narrow:

```text
PROTO11_successor_runtime_implemented = true
PROTO11_fresh_GR0_dynamic_calibration_authorized = true
```

They mean that the reference-balanced numerical map has been implemented,
bound to the inherited campaign, attacked with exact and perturbed controls,
and found fit to launch one fresh GR-0 calibration. They do not say that a
trajectory has run, that a trapped sphere forms, that FGC-QR is healthy, or
that collapse defocuses.

## Immutable evidence boundary

HLT9 reads PROTO11, PRO11-FRZ1, HLT8, CAL7/PREF9, ID2, and the CAL7 run plan
from immutable commit `b53c0f1a38823b84590e3946fb31146be99ed7b5`. It verifies
their configured hashes, their canonical records, and the implementation
hashes carried by those records against the current checkout. The current
CAL8 plan is then normalized back to immutable CAL7. Exact equality after
removing only the PROTO11 map and namespace delta proves that no physical
input, source tolerance, solver, CFL rule, retry rule, constraint threshold,
spectral rule, trapped-sign rule, boundary rule, or stop was silently changed.

The PROTO10 campaign and its CAL7 diagnosis remain disclosed history. CAL7
found a binary64 spherical-coordinate source floor before the first evolved
common event; it did not observe a candidate or holdout outcome. HLT9 does not
reinterpret that calibration-contract failure as mechanism evidence.

## Runtime binding

For the exact spherical Minkowski reference,

```text
u_ref = (1, 0, 1, r, 0, 0)
q_ref = (0, 0, 0, 1, 0, 0),
```

the operative map is

```text
q(0) = D_h(u-u_ref) + q_ref
q_r  = D_h(q-q_ref)
C_q  = q - [D_h(u-u_ref) + q_ref].
```

`p_r`, `du/dt`, and `dq/dt` remain `D_h(p)`, `p`, and `D_h(p)`. The `q`
array remains an independently evolved state variable. It is initialized once
and is never reprojected inside a Runge-Kutta stage or after an accepted step.
The reference map is used by the source, reduction diagnostic, and committed
common-event constraints so no second numerical geometry is hidden in the
admission test.

The complete unredefined REF1 source is still evaluated and its raw residual
must remain below `1e-12`. For the bitwise exact fixed reference only, after
that raw evaluation passes, the analytically known zero acceleration is
selected. This is an equality-only well-balanced representation, not an
epsilon floor. A one-bit perturbation is explicitly shown to bypass that
branch and take the ordinary affine-source path.

## Rebuilt inputs and initial events

HLT9 reconstructs the two frozen amplitudes, two methods, and three grids:

```text
2 amplitudes x 2 methods x (2049, 4097, 8193) = 12 inputs.
```

For every input, `u` and `p` match immutable HLT8 bitwise, `q` matches the
PROTO11 formula bitwise, the complete raw source passes the unchanged gate,
the operator leaves its input unchanged, and no trajectory advances. The four
amplitude/method ladders then pass the reference-balanced physical, gauge, and
reduction constraints together with the complete inherited PROTO10 spectral
contract. Raw PROTO9 evidence, the direct coarse-to-medium ratio veto, the
strict `<1/4` ceiling, the conditioning witnesses, and the statement that
round-trip interpolation is not a continuum-error bound all remain public.

## Adversarial controls

The machine record requires:

- bitwise zero `du/dt`, `dp/dt`, and `dq/dt` for exact Minkowski on both
  methods and all three grids;
- ordinary-source routing for a one-bit perturbation;
- visible reduction and source-derivative defects for a generic perturbed
  `q` field;
- bitwise agreement between committed reduction rows and the frozen map;
- immutable input arrays across operator calls;
- a demonstrated distinction from the old direct fourth-order derivative map;
- no new tolerance, damping, cleaning, or interior reprojection.

The second-order operator happens to differentiate the linear reference
exactly in binary64. That coincidence does not weaken the contract: the
fourth-order boundary stencil exposes the old-map discrepancy, and both
methods are bound to the same explicit PROTO11 algebra.

## Decision and nonclaims

HLT9 authorizes only a fresh, outcome-neutral GR-0 calibration in
`runs/fgc-2-sf1/proto11/calibration`. The holdout namespace remains empty and
closed. In particular:

```text
PROTO11_resolved_holdout_manifest_authorized = false
classical_spherical_diagnostic_authorized = false
FGCQR_holdout_execution_authorized = false
retained_EFT_evolution_authorized = false
physical_transition_claim_authorized = false
```

Therefore HLT9 is a numerical-premise result. It is not evidence for
defocusing, singularity resolution, a child domain, a dark-sector mechanism,
or varying locally measured light speed. Its value is that the next trajectory
will test the frozen experiment rather than a known spherical-coordinate
roundoff artifact.

## Reproduction

Before the first trajectory:

```bash
python3 scripts/reproduce_fgc_hlt9_mon9.py
python3 scripts/reproduce_fgc_hlt9_mon9.py --check
python3 scripts/run_fgc_gr0_calibration_v11.py --check-authorization-only
```

The authorization command must leave both PROTO11 output roots absent or
empty. Once the historical calibration namespace exists, canonical repository
verification reuses the serialized pre-launch namespace observation rather
than pretending the run never occurred.
