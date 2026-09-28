# Local incoming gate after the n=64 truncated steps

Three evolved steps from iterate6 stay in the declared `(w, U)` class.
None improves both components. Rank 64 matches beta and leaves an N
offset that does not shrink when the step is damped by 32. The gate is
OPEN. Iterate6 remains the best measured residual.

```sh
python3 scripts/derive_nsc_local_incoming_gate_status_after_n64_truncated.py --check
```
