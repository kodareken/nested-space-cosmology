# Coupled-history retained finite-source evaluation

This control evaluates the committed coupled `w`/`U` history on every original
nonzero-angular source family. It restores the frozen upstream columns with
`RetainedUpstreamArchive` and does not repeat source preparation. The geometric
seed in the [group14 throughput pilot](nsc-ks-coupled-history-pilot.md) remains
solver control only. The physical state is the original `C_up` evolved through
that history.

Each of the sixty positive-angular families is paired with its explicit
negative-energy opposite-angular partner, giving one hundred twenty distinct
signed contributions. A per-family `KSHistoryEvaluator` uses the archived
positive interpolation interval `[0,160]` or `[0,320]`, degree32, the64-point
spatial cell, and the common union of eight solve and seventeen verification
nodes. Every operator carries all sixteen retarded directions. Family receipts
retain only the signed matter corrections and retarded matter tangents.
Baseline and local/reference geometry are added once in the aggregate.

The unit-angular source-cutoff phase and its sixteen tangents are computed once
on that complete target union. Original signed angular-square weights, one
hundred seven thousand nine hundred fifty-two for each energy sign, contract
the derived edge. Groups 0, 13 and 23 keep their exact zero radius responses
and their baseline allocations. Both matching group14 operators are reused from
the committed pilot payload; at most fifty-eight new operators are evolved.

Per-family operator arrays and bound JSON are written immediately and resume
only for the same inputs, hashes and solver configuration. Each new operator
is limited to forty CPU seconds; the full-source computation is limited to six
hundred CPU seconds. A budget stop or timeout leaves the record explicitly
OPEN and never issues a NON-EXISTENCE claim. Missing source, field, finite
tail, phase-quadrature and between-node errors remain `None`. This evaluation
does not close those bounds and does not claim a physical root.

Replay reloads the saved operators and original sources, reconstructs the
matter assembly and the derived edge, and does not call Dirac evolution.

```sh
python scripts/derive_nsc_ks_coupled_retained_control.py --run
python scripts/derive_nsc_ks_coupled_retained_control.py --check
```
