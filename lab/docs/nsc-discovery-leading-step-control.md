# One numerics-only leading-Einstein successor

This successor addresses the raw-coordinate reaction-norm bottleneck of the
frozen e6e784 leading-Einstein batch. It changes only the timestep policy.
The action, RHS, source, canonical maps and field generator remain the exact
same committed files. Removing a primitive numerical bottleneck is progress;
it is not evidence of regeneration or global stability.

Owners are the [step-control module](../src/recursive_horizons/nsc_discovery_leading_step_control.py),
[driver](../scripts/derive_nsc_discovery_leading_step_control.py), and
[focused tests](../tests/test_nsc_discovery_leading_step_control.py).
The e6 leading files and v1 results remain preserved.

## Full physical local norm

At each actual prolonged state, let
\(\omega_*=\max(k_{\rm geometry\ band},k_{\rm field},1)\), and choose
\[
S_{\rm geo}=
\left(Q,\ r,\ \frac{2|a|r^2\omega_*}{Q},\ 2|a|r\omega_*\right),
\qquad \psi=\Phi/\sqrt{\Delta x_q}.
\]
All 24 real components of the six two-component complex spinors have unit
scale in \(\psi\). The local reaction Jacobian \(J_0\) is the existing new
leading analytic Jv with only the spatial D/P operators set to zero. No
geometry/source/field links are discarded.

For each of the 28 input components, the implementation applies its local
scale as a tangent, divides the output by its scale, and sums absolute
row entries. This computes \(\|S^{-1}J_0S\|_\infty\) pointwise.
The moving-coordinate scale rates use the actual projected/lifted rates:
\[
g=\left(\dot Q/Q,\ \dot r/r,\
2\dot r/r-\dot Q/Q,\ \dot r/r\right).
\]
Field scales are fixed. The reaction majorant is
\[
b=\max_x\|S^{-1}J_0S\|_\infty+\max_x|g|.
\]
The absolute triangle term avoids cancelling genuine growth through the
moving-coordinate sign. The admission frequency is
\[
\omega_{\rm admit}=\max(k_{\rm quadrature},k_{\rm field})+
b+2\max(\|DQ/Q\|_\infty,\|Dr/r\|_\infty),\qquad
dt=\min(dt_{\rm cap},1.4/\omega_{\rm admit}).
\]
There is no implicit damping, reaction removal, or changed evolution
coordinate. The scaling applies to the numerical indicator only.

The independently reviewed endpoint local scaled norms were about228–438,
compared with raw norms of \(10^5\)–\(10^6\); full local eigen magnitudes were
below123 and real growth about29–38. That growth remains present.
Whole finite weighted-norm reference values3435/7799 are retained as review
context. They are not computed as a global Jacobian every step, and this
cheap local admission indicator is not a global stability certificate.
The bounded paired-step check compares one new step with two half steps at
each saved endpoint and verifies bit-identical old/new RHS evaluation.

## Exact handoffs and bounded continuation

Every case copies its exact latest v1 NPZ bytes, complete per-array hashes,
geometry, momenta, fields, W, source weights, normal clocks and work ledger.
No constraint solve, interpolation, source selection, eigenframe regeneration
or state reset occurs. Two already-completed frozen arms stay terminal and
are never evolved again.

The old observation file is copied by bytes and its prefix hash/length are
recorded. Its numerical domain ends at the old checkpoint time and is never
relabelled as the new policy. New rows append afterward. Carried work
quadrature can bridge this numerics boundary; it remains a sampled ledger,
not a closure or propagated-error certificate.

A distinct successor schema rejects accidental use of the old raw-policy
CLI. A scoped adapter replaces only the leading schema, timestep function
and pool initializer hook. The custom spawn-safe initializer installs the
new schema/function before calling the existing leading worker initializer.
The same single episode pool, admission ledger, event retention, source
operators and immutable64-MiB chunks are reused.

The fresh successor budget is at most165 CPU seconds, including preparation,
with parent plus successor at most600. Parent reservations must be zero.
Source hashes at e6 and the new frozen producer are authenticated. Reader
and process overhead are reported in the successor accounting. If genuine
growth, a chart event or the remaining budget stops progress, the state is
retained and the outcome reported; the budget is not escalated.

Default invocation is pure preflight. After freezing the new four files:

~~~sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_leading_step_control.py --paired-check
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_leading_step_control.py --prepare --producer-commit HEAD --cpu-budget 165
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_leading_step_control.py --run --workers 4 --cpu-budget 165
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_leading_step_control.py --check
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_leading_step_control.py -q
~~~

The default successor owner is
lab/results/development/nsc-discovery-leading-einstein-v2.
Preparation creates a new directory exclusively. Checks read bindings,
prefixes, handoffs and constraints without evolving or writing. A paired
endpoint check takes only bounded local steps and never publishes a
trajectory or checkpoint.
