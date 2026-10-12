# Conformal retained-region reduction feasibility

The finalized initial conformal episode supplies an unchanged-source
interface check. Its regional throughflow is not yet resolved. This note
records only the two-step reduction probe on \([0.025,0.035]\); the full
\([0.025,0.05]\) comparison was forecast and was not run. Production
reduction remains held for the finalized same-source continuation with a
resolved regional exchange. This is not a programme-completion verdict.

The source is `nsc-spherical-conformal-episode-v1`, primary
`nf512_dt_0.0005`. Its positive-chart and sampled constraint-forcing
assessment passes. The time integral and continuum constraint are still
uncertified in that episode. No geometry is regenerated here. The
[module](../src/recursive_horizons/nsc_conformal_local_response.py) is
separate from the sealed prescribed-gauge local-response consumer.

The field operator is the representative block

$$
H(t)=\sigma_2P+\kappa Q(t)\sigma_1.
$$

Coarse stored \(Q\) is prolonged to its recorded quadrature and interpolated
linearly between frames. Canonical column multiplication is
\(U_f^\dagger Q U_f\). The initial observer is still the original
\(T=0\) region-0 pair, columns 0 and 1, without QR or a phase reset.
The initial state is the actual transported \(\Phi(0.025)\), with the six
original Gaussian weights. The isotropic copy factor stays in the source
energy and force; it does not multiply this generator or the occupation.

The actual initial cross block has Frobenius norm
\(0.176498626520663\). The two retained eigenvalues are
\(0.7278969653465316\) and \(0.7281824992528759\). The full weighted
covariance remains admissible, with nonzero eigenvalues reproducing the
original weights to about \(10^{-14}\). No synthetic cross block is used.

| Probe quantity | Value |
|---|---:|
| Autonomous occupation change | \(0.02051180779146944\) |
| Conditional full occupation versus saved autonomous frames | \(2.858824288409778\times10^{-12}\) |
| Streamed occupation versus conditional full occupation | \(5.915941789003121\times10^{-6}\) |
| Memory omission occupation movement | \(0.003654370068278623\) |
| Initial exterior drive omission movement | \(0.016803736908439615\) |
| Drive-off error versus independent full projected-initial evolution | \(1.7964660790070752\times10^{-6}\) |
| Actual initial cross omission movement | \(0.017170317079074918\) |
| Cross error versus independent full retained/exterior column split | \(4.25038609645767\times10^{-6}\) |

The headline reduction error is about \(0.02884\%\) of the occupation
change. Drive-off error is about \(0.01069\%\) of its effect; cross error
is about \(0.02475\%\) of its effect. The cross superposition residual
\(1.11\times10^{-15}\) is an algebra check, not that reduction error.
The memory row is a probe separation; its own temporal refinement is left
to production controls. These controls preserve the source preparation and
geometry except for their named conditional omission, and do not evolve
a new coupled geometry.

An independent quadrature Fourier action matches the matrix on the actual
initial columns with relative error \(4.21\times10^{-13}\). Tests also
compare it directly with the owned fine-grid Dirac action, check the
representative normalization, and reject partial, changed-source,
failed-chart and failed-sampled-budget inputs. The production entry point
rejects this initial window as feasibility only.

The streamed history is \((3,1022,6)\), \(294336\) bytes. Time-indexed
exterior propagator storage is zero. The Hamiltonian cache holds at most
two matrices. The comparison measured \(2.456365\) CPU seconds and
forecast a five-step comparison at \(15.352281\) seconds before scaling.
That full comparison was not admitted for this scientific scope. The
probe and its controls together used \(13.461717\) CPU seconds under the
\(60\)-second feasibility limit; saved JSON plus NPZ is \(7545\) bytes.

The record is
[JSON](../results/development/nsc-conformal-local-response-feasibility-v1.json)
and the series are
[NPZ](../results/development/nsc-conformal-local-response-feasibility-v1.npz).
The bound conformal-episode NPZ hash is
`5949bc30763d7cc9855c0e1d6398b725d402013a9ee21c5c424aedd2eeec1d5b`.
The old local-response v1/v2 payloads and their consumer code are unchanged.

```sh
python scripts/lab.py scripts/derive_nsc_conformal_local_response.py --budget-s 60
python scripts/lab.py -m pytest tests/test_nsc_conformal_local_response.py -q
```

The future state-dependent Duhamel map remains conditional as stated in
[the local-response note](nsc-coupled-local-response.md#conditional-geometry-to-occupation-error-relation).
A controlled \(\delta H\Phi W^{1/2}\) integral would imply an occupation
error on the same fixed observer. The current probe does not supply that
integral for an unknown geometric path, claim a stress, or establish
maintained throughflow.
