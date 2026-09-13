# Extended action completion: OPEN after the first-differential calculation

The [F→G→H record](../results/development/nsc-extended-action-completion-gate.json)
adds an executable Gaussian action differential and endpoint chain rule to the
v0.19.0 binding gate. The previous component and non-existence certificates
remain unchanged.

| Step | Computed contribution | Physical gate |
|---|---|---|
| F | First differential of the existing CTP action on all 33 NSC channels, checked against that action | OPEN: admissible physical link/embedding variations still required to select $V_c$ |
| G | Conditional contraction into the existing four-field, two-endpoint covector | OPEN: physical branch jets and same-action remainder not supplied |
| H | Composition of these computed and missing residual channels | **OPEN** |

The common-branch action is flat for all channel unitaries. The
difference-history derivative supplies the finite Gaussian source response;
it cannot be treated as an unrestricted equation for an independent unitary
without defining that variational problem. This distinction identifies the
specific physical data the remaining completion must supply.

The necessary input is the transmitting action's functional dependence on
link, embedding and metric, including its normal-jet boundary data. The
existing owners evaluate supplied histories and vertices. Naming their
uncomputed dependence does not evaluate it, and no optimizer is called to
choose it implicitly.

No full-history metric variation, finite stress, null sign or updated
constraint is emitted. The homogeneous frequency-diagonal no-interface
non-existence certificate remains valid within that old domain and is not
extended to the tilted class. All locked coefficients and scales are
unchanged; metric stepping and coupled evolution remain closed.

```sh
python3 scripts/derive_nsc_transmitting_action_completion.py --check
python3 -m pytest -q tests/test_nsc_transmitting_ctp_variation.py
```

Only this new first-differential calculation and its binding are executed.
The old covariance-freedom probe, source integrals and history controls are
read as authenticated inputs.
