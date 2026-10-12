# Reference value accuracy on the approved local interval

The completed boundary-buffer record has fine-grid sampled errors
(7.94880e-11,8.17611e-11) on I=S(1)+[.12,.18], above3e-11. Its coarse/fine
ratio is near256 for a fourfold spatial step change. This specific observed
spatial defect motivates one h=.0005 control, with the same operator and source.

The owned radial sampler evaluates the same at_zero source columns on the
new6401-node prefix. The same canonical continuation supplies new nodes above
rho1.2. Common prefix columns must agree within3e-11 before any propagation,
and every existing3801-grid column and coordinate is then copied bitwise into
the nested7601 grid. No normalization, occupation, phase or source refit occurs.

One exact-harmonic reference evolution supplies61 nodes on[0,.18], including
the21 local nodes[.12,.18]. Its raw matter insertion uses actual PDE F_z and
the original coherent source/multiplicity. The budget is90 CPU seconds; no
extra case follows a failure. No existing field case or certificate is rerun.
The shorter sampled coordinate interval evaluates the already approved local
I and does not select a physical duration.

The numerical target is3e-11 on these local samples. A PASS is strictly a
four-energy reference-value control. Full source, changed-history, continuous
space/time and between-sample errors remain OPEN. Drift is never subtracted.

```sh
python3 scripts/derive_nsc_local_reference_accuracy.py --run
python3 scripts/derive_nsc_local_reference_accuracy.py --check
```

The saved record and payload support array-only replay and future nonzero
history evaluation with the same refined initial columns.
