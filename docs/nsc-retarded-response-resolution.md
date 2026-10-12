# Targeted resolution of the retarded incoming response

## Reuse, gap, decision and stop

Reuse the exact pulse, source covariance, four signed frequencies and weights
of the [retarded pilot](nsc-retarded-compatible-response.md), and its frozen
geometry, SBP/SAT and prepared-field owners. The failed spatial comparison
is the reason for this extension: common-node forcing agrees, whereas its
discrete spatial transport differs substantially. This record changes only
spatial sampling, followed by the declared time check. It selects neither
physical initial data nor a new source law.

Use the authenticated group14_1 at_zero columns with the existing radial
reference sampler once at3201 points on[-2,1.2], retaining its original
ODE tolerance. Within-spin subsampling supplies1601 and801 nodes. Compare
all new801 columns against the archived801 columns before propagation;
maximum error must be below3e-11. This evaluates the same reference modes
on missing grid nodes, without horizon, scattering or source-state solves.

Run1601/128 and3201/128. All three spatial max-entry relative indicators
(field tangent, actual rho1 trace tangent, coherent covariance tangent)
must be below25%; otherwise stop OPEN. Only then run3201/256, requiring
the same three time indicators below10%. At most three field solves and
120 CPU seconds including reference resampling are allowed. Checkpoint
completed cases. No automatic refinement follows an OPEN result.

The physical pulse remains lambda0, U=0, normal radii0.007/0.03, the same
compact axial bump and PG interval[0,0.3]. Preparation, right incident and
source-covariance tangents vanish by the unchanged past/upstream geometry.
The fixed-surface constructor uses(N,q,beta,r), without a second stationary
phase. All source coherence and the single sqrt(weight/(2*pi)) column
weight are retained. The full sampled projector tangent is reported and
never forced to zero. Report reference drift separately; do not subtract it.

The slow characteristic is the first characteristic row. Its peak magnitude
and time are reported both on the full trace and after the known global
pulse stop; this late window is a diagnostic, not a physical duration.
Comparison with the old801/128 result is context only and never reruns it.
The independent causal diagnosis also supplies the trace-support formula
`[.09-Xplus,.21+Xplus]`, where `Xplus=integral_1^rho_plus d rho/(beta^2-1)`.
Its numerical endpoint estimates are[.0540479548,.2459520452]. Maximum
entries and trapezoidal sampled energy fractions outside that interval are
reported without cropping the trace or changing the acceptance criteria.
The endpoints are diagnostic estimates, not directed interval bounds.

Three likely failure modes are an archived/new-node mismatch, unresolved
spatial transport despite finer sampling, and exhaustion of the CPU budget.
Each yields OPEN at its declared gate, rather than a wider computation.

## Record and scope

The [JSON record](../results/development/nsc-retarded-response-resolution.json)
and its authenticated arrays contain the source binding, resampling residual,
completed cases, comparisons, late slow trace, stationary drift, algebra,
Hermiticity, runtime and stopping reason. Replay only contracts saved arrays;
it does not sample modes or propagate fields. PASS means these sampled
numerical indicators pass. It is not a continuum error bound, fullC0/CAR
verification, stationarity, constraint solution or physical metric history.
The fixed rho0 seam and all Gamma terms remain unchanged.

```sh
python3 scripts/derive_nsc_retarded_response_resolution.py --run
python3 scripts/derive_nsc_retarded_response_resolution.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_retarded_response_resolution.py
```
