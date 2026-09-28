# Incoming evaluator with the derived source-cutoff edge

`KSCutoffBridgeEvaluator` wraps the actual
[`KSHistoryEvaluator`](nsc-ks-history-evaluator.md). It keeps the source-fixed
Dirac columns and retarded Jacobian and adds the
[same-action cutoff edge](nsc-source-cutoff-bridge.md) to N,beta and their
history derivatives. A covariance or matter dictionary cannot replace the
wrapped state owner.

The pure-radius phase coefficient is proportional to ell squared. One
unit-angular geometry integral is therefore contracted with the original
sum of `multiplicity*ell²`, separately for positive and negative source
labels. Physical source evolution keeps every original angular eigenvalue,
source covariance and weight. Each represented signed family enters once,
irrespective of how many spectral batches it contains; no energy-sign
folding or extra spin factor is applied.

The phase and state share the analytic profile, all amplitude directions,
upstream coordinate and complete target-node union. Exact node subsets
reuse both calculations. The geometric Jacobian remains available as a
preconditioner, while the full Jacobian includes the state and edge
derivatives. The returned record retains the raw arrays separately so the
addition can be replayed and inspected.

The finite input-energy tail, phase quadrature and other source/field
errors remain explicitly unbounded. A partial source selection remains a
partial source selection; the edge does not fill its missing energy panels
or angular families. Neither this wrapper nor a small nodal residual gives
physical EXISTENCE. No local/reference coefficient, baseline action or
`Gamma_rest` is changed.
