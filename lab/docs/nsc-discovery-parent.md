# Global parent preparation from the measured responsive collar

This constructor prepares one periodic Cauchy slice on the unchanged leading
Einstein/Dirac action. It authenticates the sealed `.0013` finite collar, keeps
its initial geometry through distance one from the carrier center `L/2=4`,
blends through `1..1.2` to magnetic radius and `Q=2`, and continues the real
Dirac spinor through the actual blended `Q`. A flat initial envelope reaches
zero at distance `2.8`, before the parent boundary at three. The second real
source column is an annular packet at distances `1.2..3`. These are initial
inputs; no wall, pulse, new force, density replacement or evolution is added.

Both spinor columns are projected into the AP carrier before normalization.
They remain generally nonorthogonal. For the projected raw child column with
norm squared `ZC`, the normalization identity `c=w0 ZC` preserves its raw
sampled covariance. Physical population interventions subsequently change
weights and are reported as changes to that covariance. Actual finite-CAR
admissibility uses the eigenvalues of `sqrt(c) Gram sqrt(c)`.

Let `a_i` and `b_i` be the actual normalized-column content in the child and
parent annulus. The fixed trace is `T=2 w0 ZC`. Define
`delta(c)=sum c_i(a_i-b_i)/sum c_i(a_i+b_i)`.
The two pure-column endpoints must straddle zero. The common symmetric
`delta_star=min(.5,max endpoint,-min endpoint)` determines child-heavy,
balanced and parent-heavy target imbalances `+delta_star,0,-delta_star`.
A two-equation solve finds each weight pair at fixed trace. A balanced case
is not assumed to have equal weights, and the parent-heavy case is not
produced by reversing the child-heavy weights. Failure to straddle is a
preparation obstruction. The actual columns and weights supply the full
source kernel, including its current; neither density nor current is deleted.

An independent orthonormal local observer uses two child spinor component
probes followed by QR. It is distinct from the physical nonorthogonal source.
Its finite-band leakage is part of this carrier, rather than exact spatial
support. All geometry modes remain retained; the geometry map is identity.

With `g=8pi A`, compute the zero-momentum actual finite constraint `C_s` and
its periodic primitive `J` of `4g r_x C_s/(r²Q)`. Removed mean/Nyquist modes
and the primitive reproduction gap are reported. A common positive offset
above the largest population `kmin` seeds

```text
P=±sqrt(r³(k+J))
p_r=3QP/(2r)+2gr C_s/P.
```

This continuum relation is only a seed. Both retained momenta are corrected
against the actual finite even `h=QC`, odd `D`, and the fixed kinetic anchor
`mean(P²/r³)` measured on that seed. The anchor is not forced to equal `k`.
Analytic finite momentum Jacobians are solved in the physical metric
`dxq Bfine.T diag(SP^-2,SR^-2) Bfine`, where
`SP=2g r² omega*/Q`, `SR=2g r omega*` and
`omega*=max(kgeometry,kfield,1)`. Cholesky whitening precedes SVD, with
physical gravity/source/kinetic row scales. Actual unscaled residuals,
excluded parity components and projection tails remain visible.
The scalar kinetic anchor is a hard equality: its analytic derivative defines
a nullspace for the C/D SVD solve and is never subjected to the relative
singular-value cutoff. Each trial retracts retained `P` to the exact quadratic
anchor level without changing fields or weights. Convergence requires actual
anchor relative error at most `5e-11`, and trials retain a strictly positive
declared momentum-sign margin (at least one tenth the initial minimum and
`1e-10` of the momentum scale). These requirements correct the discarded-anchor
defect in the original immutable pilot; that evidence remains unchanged.

Correction permits at most twelve accepted momentum updates in total, eight
backtracking candidates per update, and one fallback with four even geometry
parameters (means and first cosines of `r,Q`) under the same CPU allowance.
The fallback preserves the fields, weights and kinetic anchor and recomputes
the actual source when geometry changes. Positive chart and momentum sign
checks bound trial updates. A small correction is a numerical stationarity
indicator, not a continuous constraint certificate or a universal `1e-8` gate.
The two momentum signs preserve field arrays byte for byte; they do not
complex-conjugate the source or claim complete dynamical time reversal.

The stable API is
`prepare_parent(nf,population=0,sign=1,k_override=None,profile=None,cpu_limit=30)`
returning `(ParentPair, leading.State, report)`. The pair exposes its full
FFT grid, geometry map, weights, physical source columns and independent
observer columns. Diagnostics use the rank-general parent timestep owner,
with twelve local real variables for a rank-two source.
The expansion diagnostic is the normal angular rate
`Hangular=projected rdot/(r²Q)`, using the actual decoded/prolonged metric
velocity. Coordinate `projected rdot/r` and normal radial
`(projected rdot/r+projected Qdot/Q)/(rQ)` rates are reported separately.
Unprojected pointwise rates are labeled raw diagnostics and do not supply
expansion or acceptance decisions.

Default CLI operation only authenticates and previews the input, without
solving a spatial ODE or calling the constructor. Root
freezes the four producer files before creating a new record, for example:

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_parent.py --prepare --nf 64 --population 0 --sign 1 --cpu-limit 30 --output results/development/nsc-discovery-parent-v1/nf64-child-plus
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_parent.py --check --output results/development/nsc-discovery-parent-v1/nf64-child-plus
```

Use absolute `--output` paths when invoking the root launcher, whose working
directory is `lab`; relative examples should therefore omit the leading
`lab/`. Creation-only `parent.json/parent.npz` records bind both the sealed
input chain and producer commit, retain finite failed preparations, and stay
under 64 MiB. Read-only checks authenticate historical producers; numerical
replay requires matching current bytes. The thirty-CPU-second allowance is
for one constructor call, including population-offset comparison, correction
and diagnostic work. Root owns the aggregate feasibility pilot and chooses a
shared `--k` when comparing separately prepared cases/resolutions. No
trajectory, global stationary solution or autonomous renewal is claimed.

For the single bounded NF64/128 feasibility pilot, root may instead run:

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_parent.py --pilot --population 1 --sign 1 --cpu-limit 30 --output results/development/nsc-discovery-parent-v1/feasibility-nf64-128
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_parent.py --check --output results/development/nsc-discovery-parent-v1/feasibility-nf64-128
```

The pilot authenticates the frozen producer before computation, compares
both resolution seed families, and uses the largest suggested `k` for both
cases. It divides the measured remaining CPU allowance across the two
preparations and writes immutable `nf64/parent.*`, `nf128/parent.*` plus
`pilot.json`. The aggregate ledger records the shared `k`, each `kmin`, and
each actual projected kinetic anchor. Equal declared offsets do not replace
inspection of differing anchors/source/projection readouts. Deferred
preparations remain named numerical obstructions rather than static-force
vetoes or evidence of a completed global solution.

A successor must use a new output prefix. The thin `--transition-end` option
forwards an existing profile parameter; for example root may repeat the bounded
pilot with `--transition-end 1.4` and a new prefix. The protected inner collar,
actual source continuation, finite constraints and hard anchor are unchanged
in their definition. This option does not launch an additional solver or fit
the source to the constraints.
