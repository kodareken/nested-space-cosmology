# Local incoming-gate certificate v2

Certificate v2 is a decision owner, not a blocked-status result. It consumes
the declared local class, the unchanged source identity, a continuous residual
bound, all nine componentwise errors, the real mutation results and a complete
provenance graph.

EXISTENCE is emitted only when both exact componentwise sums satisfy

\[
\sup_I |\widehat{\mathcal E}_B|+\epsilon_B\le3\times10^{-11}.
\]

NON_EXISTENCE is emitted only from a necessary relation covering the entire
declared class with strict error-controlled separation. An optimizer stall is
explicitly rejected. Missing bounds or mutation controls produce OPEN, and the
recording CLI refuses to write an OPEN certificate.

The current audit remains OPEN because the corrected best residual is nodal,
six budget components are missing and the final mutation campaign has not run.

```sh
python3 scripts/derive_nsc_local_incoming_gate_certificate_v2.py --audit-open
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_local_gate_certificate_v2.py
```
