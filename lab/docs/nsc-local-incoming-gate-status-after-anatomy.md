# Local incoming-gate status after leftover anatomy

The leftover is now classified. It is still not a certificate.

- Solve-node leftover after the unclipped jet-constrained step is about
  `6.26e-4` and equals the orthogonal leftover. Jet rows are `3e-12`.
- All-node leftover of that step is about `6.34e-4`. All-node collocation
  makes it worse.
- The leftover direction is not sign-stable on the five accepted histories,
  so there is no scoped NON-EXISTENCE.
- Dropping the solver jet kills the solve-node leftover and leaves all-node
  maxima about `(1.52e-4, 7.77e-6)`. That is still far above `3e-11`. The
  unit clip stays on the current residual scale. No family was evolved.

State law, class and interval are unchanged. OPEN is not Done.

```sh
python scripts/derive_nsc_local_incoming_gate_status_after_anatomy.py --check
```
