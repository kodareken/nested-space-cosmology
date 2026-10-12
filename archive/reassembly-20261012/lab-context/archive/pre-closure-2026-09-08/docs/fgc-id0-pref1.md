# FGC-1-ID0-PREF1: frozen-protocol initial-data preflight

## Result in one sentence

`FGC-1-ID0-PREF1` rejects **FGC-2-SF1-PROTO1 as an executable study
protocol**, before any FGC-QR holdout is opened: PROTO1 leaves the regulator's
unit-normal momentum unspecified, and under the minimal completion
`Pi_phi=0` the complete declared `chi`-amplitude box has an exact conservative
initial Misner--Sharp compactness bound below `1/64`, while PROTO1 requires at
least `1/10`.

This is a pre-holdout protocol obstruction. It is not a failure of the
FGC-QR action, a collapse result, or evidence for or against a gradient phase
transition.

In compact form: FGC-2-SF1-PROTO1 as an executable study protocol is rejected;
this is not a failure of the FGC-QR action.

## Why this gate exists

RUN1 requires a nonzero-width, finite-mass, constraint-compatible initial-data
family. The frozen PROTO1 protocol already supplied a profile, a candidate
amplitude list, an intended ingoing momentum relation, and an initial
compactness interval. Those declarations must be mutually compatible before
the expensive constraint and evolution machinery is allowed to inspect the
FGC-QR branch.

The preflight found two premise-level problems:

1. PROTO1 freezes `phi(r)` but does not freeze its unit-normal momentum
   `Pi_phi`. Therefore its initial data are under-specified. An arbitrary
   unbounded `Pi_phi` would also make a finite compactness conclusion
   impossible.
2. With the least invasive completion `Pi_phi=0`, the stated amplitudes are
   normalized as a bump of `r chi`, not a bump of `chi`. Dividing by the pulse
   radius suppresses both the field and its energy enough that the original
   candidate box cannot reach the protocol's own compactness floor.

The frozen PROTO1 file is preserved. Its revision policy says a premise-based
repair requires a new protocol version; that was the only permitted next move,
and it has now been performed as **FGC-2-SF1-PROTO2** without inspecting a
candidate outcome.

## Declared GR-0 slice

The preflight uses a regular maximal polar-areal slice solely to test whether
the protocol's prescribed matter scale is plausible before constructing the
full FGC family:

```text
alpha = 1,
shift = 0,
R = r,
K^r_r = -2 k,
K^theta_theta = K^phi_phi = k.
```

The trace is zero. With the repository's extrinsic-curvature convention,

```text
lambda_t = 2 lambda k,
R_t = -r k.
```

Regular vacuum data inside the compact support set `lambda=1` and `k=0` at
the inner support boundary. This removes an otherwise arbitrary Schwarzschild
mass seed and the singular homogeneous `r^-3` momentum mode. The preflight is
therefore testing the independently declared scalar pulse, not a hidden black
hole inserted underneath it.

## Specialized physical constraints

For

```text
rho = 1/2 [Pi_phi^2 + Pi_chi^2
           + lambda^-2 (phi_r^2 + chi_r^2)] + V(phi),
P   = Pi_phi phi_r + Pi_chi chi_r,
V   = mu^2 phi^2/2 + g4 phi^4/4,
```

the unredefined GR-0 Hamiltonian and radial-momentum constraints reduce to

```text
H = M_Pl^2 [(1-lambda^-2)/r^2
            + 2 lambda_r/(r lambda^3)
            - 3 k^2] - rho = 0,

M = 2 M_Pl^2 [k_r + 3 k/r] - P = 0.
```

The resulting first-order radial system is

```text
lambda_r = r lambda^3/2
           [rho/M_Pl^2 - (1-lambda^-2)/r^2 + 3 k^2],

k_r = P/(2 M_Pl^2) - 3 k/r.
```

The certificate maps two nontrivial rational ADM fixtures back into the full
four-dimensional spherical two-jet and evaluates the complete unredefined
ACT1/VAR1 metric residual. At both fixtures the full normal projections equal
the specialized `H` and `M` fractions exactly. The comparison does not use the
REF1 gauge extension.

## Misner--Sharp identity and finite exterior

On the declared slice,

```text
C(r) = 2 m(r)/r = 1 + r^2 k^2 - lambda^-2.
```

Differentiating this identity and using both constraints gives

```text
m_r = [r^2 rho + r^3 k P]/(2 M_Pl^2).
```

Outside the pulse, `rho=P=0`, so `m` is constant and

```text
k(r) = k_out (r_out/r)^3.
```

Consequently `2m/r` decreases through the exterior. A maximum cannot be
hidden beyond the exact outer vacuum buffer.

## Entire-box analytic obstruction

PROTO1 declares

```text
r chi = A B((r-12)/2),
0 < A <= 1/8,
phi = epsilon B((r-12)/2),
epsilon = 1/131072,
support = [10,14],
M_Pl = 2,
Pi_phi = 0                 (minimal completion),
Pi_chi = A B_r/r.
```

For the compact bump

```text
B(x) = exp(1 - 1/(1-x^2)),  |x|<1,
```

write `y=1/(1-x^2)>=1`. Then

```text
|B_x| <= 2 y^2 exp(1-y) <= 8/e < 3,
```

so the half-width-two pulse obeys `|B_r|<3/2`. Across the whole support this
implies the exact rational bounds

```text
|Pi_chi| <= 3/160,
|chi_r|  <= 1/50,
|P|      <= 3/8000.
```

Assume provisionally that `lambda^-2<=4`. The matter and regulator terms then
give

```text
rho <=
1439999732348065153649
----------------------------------------------
1475739525896764129280000

     ~= 9.757817738689482e-4.
```

The momentum constraint and regular interior give

```text
|r^3 k| <= 333/1000.
```

Integrating the exact mass identity across `[10,14]` produces

```text
|m| <=
3927454576397784264080597
--------------------------------------------
55340232221128654848000000

    ~= 7.096924640114358e-2,
```

and therefore, everywhere,

```text
|2m/r| <=
3927454576397784264080597
---------------------------------------------
276701161105643274240000000

       ~= 1.4193849280228715e-2
       < 1/64
       < 1/10.
```

The same estimates improve the provisional metric bound to

```text
lambda^-2 <= 1.0142049381802287 < 4.
```

Thus a first-exit bootstrap argument closes the assumption. The separation is
strict, applies to the complete original amplitude interval under the stated
minimal completion, and is not a finite parameter scan.

## Independent floating preflight

The exact obstruction does not depend on numerics. As a regression and scale
check, the repository also integrates the two constraint ODEs with independent
RK4 and SSPRK3 implementations.

For every original amplitude it records three nested grids and requires an
observed outer-mass convergence order above `5/2`. At 4096 radial steps, the
largest original candidate `A=1/8` reaches only about

```text
max_r(2m/r) = 4.2748453e-4,
```

more than two orders of magnitude below the required `1/10` floor. Both
methods agree well inside the frozen `1e-8` comparison tolerance.

A separate GR-0-only premise-design scan checks

```text
A in {2, 5/2, 3, 7/2, 4, 9/2, 5}.
```

Those values span initial peak compactness from approximately `0.109` through
`0.680`, within the original `[0.1,0.75]` window, while remaining initially
untrapped and finite-mass. This scan neither selects a dynamical calibration
case nor inspects SGB-L or FGC-QR. It is evidence for how a successor protocol
may repair the failed premise without outcome tuning.

## Successor protocol and non-circular provenance

**FGC-2-SF1-PROTO2** uses exactly that allowed repair. It freezes
`Pi_phi=0` in the future-slice unit-normal frame, replaces PROTO1's amplitude
list with the seven GR-0-only values above, and moves all eventual holdout
output under `runs/fgc-2-sf1/proto2/holdout`. Its validator removes those
three amendments from a deep copy and requires the remainder to pass the
complete immutable PROTO1 validator. The amendment itself, including its ID0
paths and the exact compactness bound, remains inside the PROTO2 semantic
contract hash.

Updating the live RUN1, CON4, and CTR1 scope bindings would otherwise make ID0
appear to depend on its own successor. The current ID0 reproducer therefore
checks its original canonical certificate, the earlier RUN1 `2/8` stop, and
the original CON4/CTR1 positive gates directly from immutable Git checkpoint
`d4f0cc8f58408ef4ee12fb32fe3619916e231795`. Every historical source,
predecessor, implementation, and derivation ledger is rehashed from that Git
object. The checkpoint must be an ancestor of `HEAD`; missing or mutated
history fails closed. Active RUN1 then consumes the current ID0 result and
PROTO2 in the forward direction. This makes the chronology executable rather
than relying on a claim that no holdout was viewed.

## Decision and claim boundary

The exact decision is:

```text
PROTO1_initial_data_preflight_obstruction_verified = true
PROTO1_initial_data_calibration_authorized = false
PROTO1_FGCQR_holdout_execution_authorized = false
PROTO2_premise_revision_required = true
classical_spherical_diagnostic_authorized = false
retained_EFT_evolution_authorized = false
physical_transition_claim_authorized = false
```

`PROTO2_premise_revision_required=true` is ID0's historical decision. The
requirement is now satisfied by PROTO2; it is not rewritten to false because
doing so would alter what the preflight concluded at its own evidence boundary.

What has been closed is the route **PROTO1 as written**. What has not been
closed includes:

- the FGC-QR action;
- collapse-induced regulator activation;
- metric-null defocusing;
- the broader gradient mechanism;
- singularity resolution or a child domain;
- a dark-sector mechanism;
- variable locally measured light speed;
- retained-EFT validity or applicability to nature.

This failure is valuable because it prevents a non-event caused by bad input
normalization from being mistaken for evidence against the candidate physics.
It is not, by itself, the publication-grade positive or negative FGC-2-SF1
result.

## Reproduction

From the repository root:

```bash
python3 scripts/reproduce_fgc_id0_pref1.py \
  --output results/fgc-1-id0-pref1.json
python3 -m unittest \
  tests/test_fgc_initial_data_preflight.py \
  tests/test_fgc_id0_pref1_reproduction.py -v
```

The current result binds the frozen PROTO1 config, this derivation, and the
current implementation files by SHA-256. It separately validates the original
ID0, RUN1, CON4, CTR1, and every ledger they consumed from immutable checkpoint
`d4f0cc8f58408ef4ee12fb32fe3619916e231795`. Reproduction therefore requires
that Git object; a shallow archive without it fails closed. Any weakening of
the proof contract, mutation of the source assumptions, historical blob drift,
or attempted promotion fails closed.
