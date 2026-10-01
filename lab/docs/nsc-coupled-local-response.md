# Conditional local response on stored geometry

This note records one streamed retained-region comparison on geometry that
the spherical episode already generated. The episode driver is not run
again. Nothing here regenerates \(Q\), fits the old six-mode \(H\) and \(B\)
fractions, or certifies a stress.

The incoming-gate campaign stays paused. The episode JSON verdict is
recorded and is not an acceptance test. The physical arrays are the input.

## Owners

Paths are relative to `lab/`.

| Role | Owner |
|---|---|
| Consumer | [nsc_coupled_local_response.py](../src/recursive_horizons/nsc_coupled_local_response.py) |
| Driver | [derive_nsc_coupled_local_response.py](../scripts/derive_nsc_coupled_local_response.py) |
| Tests | [test_nsc_coupled_local_response.py](../tests/test_nsc_coupled_local_response.py) |
| Record | [nsc-coupled-local-response-v1.json](../results/development/nsc-coupled-local-response-v1.json) |
| Series | [nsc-coupled-local-response-v1.npz](../results/development/nsc-coupled-local-response-v1.npz) |
| Geometry | [nsc-spherical-feedback-episode-v1.npz](../results/development/nsc-spherical-feedback-episode-v1.npz) |
| Reducer | [nsc_evolving_reduction.py](../src/recursive_horizons/nsc_evolving_reduction.py), called with `backend="streamed"` and not edited |

Module sha256 `3e4942a3462962a6e308f14a2eb71e9629d14909824e444d78ad3eece5b34a15`.
Test sha256 `d349a4fe488d7f77c666c431a494c41bcd1f0af99b9ecbb548b9c30e9b62a5c4`.
Driver sha256 `29c0053260223acd0e27602037af28a198f1dbe23563d09ef39bc6d042deedb7`.
JSON sha256 `33fe043b30b037432801f12905b3545a20ca2070b9edc46bd3a11aa0d76201a1`.
NPZ sha256 `6f4bfce40a4e34a958f0ac1e3b6f763484505a9800b3e34c7362125cbe160ad1`.

## Operator

\(H(g(t))\) is the Fourier–Galerkin form of `apply_dirac` on the
antiperiodic fermion band. The prolongation sits in that band, and the
fine antiperiodic momentum preserves it, so the Galerkin matrix is the
band Fourier multiplication by the nodal profiles \(L/Q\), \(\beta\), and
\(\kappa L\), rotated back to the coarse nodal spinor basis
\([\phi_0;\phi_1]\). No fine identity image is built. For the saved
\(n_f=512\) grid that image would have been \(2048\times 2048\) or larger.
The matrix that is formed is \(1024\times 1024\).

Only the stored fine conformal factor \(Q\) enters the evolved part of
\(H\). Lapse density and shift stay the calibration samples. Stored \(r\),
\(\chi\), and the geometric momenta are not matrix entries of this Dirac
operator. A geometry-rate series for \(Q\) is not in the episode payload.
The final coarse \(\dot Q\) is kept as a norm and is not used to redraw
\(Q\). Between stored frames, \(Q(t)\) is piecewise linear. That
interpolation scale is reported and is not removed.

On an independent \(n_f=24\) control the Fourier action matched
`apply_dirac` to \(4.531090817357182\times 10^{-15}\) on random columns
and \(3.753764274373907\times 10^{-15}\) on the rank-6 packet. The band
leak was \(1.9134213540303943\times 10^{-14}\). Distinct times did not
commute: the Frobenius commutator was \(91.21078530359848\). The
occupation-weighted source covariance did not commute with \(H(0)\); that
commutator was \(24.114781824373615\).

On the saved \(n_f=512\) frames the same gate passed before the consumer
ran. The initial action residual was \(2.845628578658056\times 10^{-13}\),
the final action residual was \(2.5876300098573036\times 10^{-13}\), and
the stored final spinor rate matched \(-iH\phi\) with relative residual
\(0\). The initial fine \(Q\) matched the prolonged coarse \(Q\) with gap
\(0\).

## Preparation

The source columns and the six occupation weights are the saved initial
arrays, not a refit. The fixed observer \(V\) is the first packet's two
initial modes, columns 0 and 1, copied without a QR rephase. Their
one-body block is diagonal with entries \(0.75\) and \(0.75\). The
initial cross block between \(V\) and its complement has norm
\(7.106664191995381\times 10^{-16}\). Dropping it on the small control
changed the retained covariance by less than \(10^{-8}\). That cross
block is the source block. It is not a hand-set zero used to hide a term.

The column Gram defect on the saved \(n_f=512\) preparation is
\(2.220446049250313\times 10^{-16}\).

The stored \(dt=0.0005\) and \(dt=0.00025\) fine-\(Q\) frames differ by at
most \(4.8995807411245096\times 10^{-12}\). The consumer uses the primary
\(dt=0.0005\) frames and keeps that gap. The uniform-frame second
difference of \(Q\) has maximum \(0.0007112125979157613\), so the
piecewise-linear midpoint scale is \(8.890157473947016\times 10^{-05}\).

## Finite domain that was executed

| Item | Value |
|---|---|
| Case | `nf512_dt_0_0005` |
| Band | \(n_f=512\), quadrature \(2048\), length \(8\) |
| Spinor dimension | \(1024\) |
| Observer | initial mode columns 0 and 1, weights \(0.75\) and \(0.75\), no phase QR |
| Window | stored frames from \(t=0\) through \(t=0.035\) (7 steps of \(0.005\)) |
| Stored frames available | 11, through \(t=0.05\) |
| Reference | independent midpoint evolution of the same six columns under the same \(H(Q(t))\) |
| Reduction | streamed Volterra, memory and exterior drive on |

The finite-domain label names those initial mode columns and the stored
weights. The v1 JSON file retains the observer string from the run that
wrote it.

This is not the coupled episode's final spinor, not an RK4-node geometry,
and not a new forced-regeneration model. Output nodes land on stored
frames, so nodal \(Q\) there is the stored sample. Midpoints use the
linear interpolant.

## Observed occupation

The resolved signal is the region-0 occupation. Over \([0,0.035]\) the
full-band occupations move from \((0.75,0.75)\) to about
\((0.69758119,0.56234705)\). The maximum absolute change is
\(0.1876529543525184\). The streamed occupations differ from that
full-band series by at most \(0.0004995332964643495\), which is
\(0.0006660443952857995\) of the occupation level and
\(0.002662006032294826\) of the change. The trapezoid residual was
\(2.8118846054444253\times 10^{-16}\).

Coherence between the two local modes stays below
\(0.0007659885446204019\). Its absolute error is
\(2.560339698158102\times 10^{-6}\), and the phase error on the samples
where that coherence is resolved is \(0.0012901833396748129\) radians.
Occupation is the reported effect because its change is the larger signal.

On the first \(0.01\), halving the output step reduced the
occupation discrepancy against the finer full evolution from
\(5.788866508071866\times 10^{-5}\) to \(1.4486698050486524\times 10^{-5}\).
The ratio is about four. That is the observed step pattern on this
window, not a certified bound. The same pattern was required on the
independent \(n_f=24\) pulse against a DOP853 sample of the full field.

## Exterior-to-local omissions

These omissions are on the probe window \([0,0.01]\), where the full
occupation change is \(0.01764381842805629\) and the streamed discrepancy
is \(5.787699418158265\times 10^{-5}\).

| Control | Occupation separation from the full series | History norm |
|---|---:|---:|
| Memory omitted | \(0.019196365211694877\) | \(0\) |
| Exterior initial drive omitted | \(0.0014796716168330448\) | \(0.20002579315843946\) |
| Coupling block set to zero | \(0.017643878295767523\) | \(0\) |

Each separation is larger than the reduction discrepancy on that window.
The drive omission is the exterior-to-local number: removing the initial
exterior amplitude moves the local occupation by \(0.0014796716168330448\)
while the recorded initial exterior norm stays \(2\). Memory is larger
than the net occupation change on this short window. The headline run
keeps memory and the exterior drive on. Its history norm is
\(1.020114910796508\).

## Cost and allocation

The two-step probe took \(1.9496839999999998\) seconds of process time.
The selector then ran the 7-step comparison, the short refinement, and
the three omissions. The whole pilot used \(18.504784\) seconds. The JSON
plus NPZ payload is \(18124\) bytes. Eleven tests passed in \(0.99\)
seconds.

The 7-step allocation does not store \(W(t)\):

| Record | Value |
|---|---:|
| Time-indexed exterior propagator | \(0\) bytes |
| History shape | \((8,1022,6)\) |
| History bytes | \(784896\) |
| Dense frame, once | \(16777216\) bytes |
| Largest exterior block formed | \(16711744\) bytes |
| Peak transported columns | \(228928\) bytes |
| Dense \(W\) that was not stored | \(133693952\) bytes |

The same probe extrapolates later runs as
\(\mathrm{probe}\times(\mathrm{steps}/2)^2\times\mathrm{substeps}\).
That attributes the whole probe to the quadratic replay, so the fixed
frame cost is not removed. Under that rule:

| Target | Estimate |
|---|---:|
| All 11 stored frames, substeps 1 | \(48.7421\) s |
| All 11 stored frames, substeps 2 | \(97.4842\) s |
| 100 output nodes, substeps 1 | \(4874.21\) s |

The 100-node figure assumes a \(Q\) sample at every node. Those samples
are not in the episode payload. The replay inside the streamed backend
still grows quadratically with output nodes, and each sample still forms
the dense exterior block from the dense \(H(t)\) callback.

## Full stored window, appended

The pilot phase above is unchanged: it still ends at \(t=0.035\), and its
process time remains \(18.504784\) seconds. The stored-frame phase is a
separate record, `phases.full_window`, on the same source and observer.
Its declared geometry is the same piecewise-linear fine \(Q\). The exterior
midpoint uses two substeps. Output nodes are the 11 stored frames from
\(t=0\) through \(t=0.05\). The 100-node case was not opened.

The conditional full-band occupation changes by at most
\(0.3305919037025379\). The streamed occupation differs from that series by
at most \(0.0007353273811995242\), which is \(0.0009804365082660327\) of
the occupation level and \(0.0022242752256303342\) of the change. Coherence
stays below \(0.001504387091847762\); its absolute error is
\(5.5090520756831094\times 10^{-6}\), and the phase error is
\(0.0014932598604091177\) radians. The trapezoid residual is
\(4.616179710183972\times 10^{-16}\). History shape is \((11,1022,6)\),
\(1079232\) bytes, and the time-indexed exterior propagator is \(0\) bytes.
The dense frame is \(16777216\) bytes and the largest exterior block is
\(16711744\) bytes. The dense \(W\) not stored would have been
\(183829184\) bytes. The streamed comparison took \(41.583679\) seconds.
The whole appended phase, including the Hermite indicator and the two
omissions, took \(93.133264\) seconds. Payload is \(30058\) bytes.

Saved autonomous spinors exist only at the initial time and at \(T=0.05\).
Projected onto the same initial observer, with the saved weights, the final
occupations are \(0.647437990920533\) and \(0.41940775506601946\). The
conditional full-band endpoint differs from that saved projection by
\(3.4123144243558556\times 10^{-7}\), which is
\(1.0321822355625109\times 10^{-6}\) of the saved occupation change. The
streamed endpoint differs from the saved projection by
\(0.0007356686126419598\). That larger gap matches the Volterra discrepancy
against the conditional full evolution, so the geometry-path discrepancy
and the reduction error are separated.

A nodal \(Q\) rate series is not stored. The only saved rate is the final
coarse \(\dot Q\). Its maximum absolute difference from the last
stored-frame secant, pulled back to the coarse grid, is
\(0.07335046359043451\), or \(0.060582014306918446\) of the saved rate's
maximum. A cubic Hermite curve using finite-difference slopes of the
stored frames, not those saved rates, differs from the declared linear
midpoints by at most \(8.693988543928555\times 10^{-5}\) in \(Q\). An
independent full-band propagation on that Hermite schedule differs from
the declared full-band endpoint occupation by
\(2.915936925806939\times 10^{-7}\). Both indicators are smaller than the
Volterra error and smaller than the memory and drive effects below.

On this same declared geometry, the full-window omissions are:

| Control | Occupation separation from the conditional full series | Fraction of the occupation change |
|---|---:|---:|
| Memory omitted | \(0.36709975607817125\) | \(1.1104317799884253\) |
| Exterior initial drive omitted | \(0.028764901996227665\) | \(0.08701030386427713\) |

Both exceed the reduction error \(0.0007353273811995242\). The drive
separation is the exterior-to-local consequence over \([0,0.05]\): the
initial exterior amplitude, whose norm remains \(2\), moves the local
occupation by \(0.028764901996227665\). Memory remains larger than the net
occupation change. The field cross block is still
\(7.106664191995381\times 10^{-16}\) and is not an active omission. The
active cross control is the prescribed six-mode covariance: its cross norm
is \(0.0023549907515443206\) and dropping it moves the retained covariance
by \(0.0014197589797623497\).

## What this does not claim

No stress is claimed. No effective-action variation was derived, so
occupied and empty kernels were not turned into a force. A force
comparison would need a kernel metric variation; that variation is not
computed. The run does not close the incoming gate and does not regenerate
the spherical trajectory. The pilot window remains \([0,0.035]\). The
\(T=0.05\) comparison uses the 11 stored frames, not the unstored RK4
nodes. A passing \(n_f=24\) control is only the operator and omission test.

## Interface gaps

The reducer was not modified. Three limits remain:

- `evolve_retained_region` accepts a dense Hermitian callback. A
  `LinearOperator` is rejected, so every sample still builds the dense
  exterior block. The rejection was observed as `TypeError`.
- `backend="streamed"` still completes one dense `null_space` frame and
  replays the interaction history from the initial time at each new
  output node.
- A full covariance argument evolves one column per mode and stores
  history shaped `(time, N_outside, N)`. The field comparison uses the
  six source weights, so the stored history is `(time, N_outside, 6)`.

## Reproduction

From the repository root:

```sh
python scripts/lab.py -m pytest tests/test_nsc_coupled_local_response.py -q
python scripts/lab.py scripts/derive_nsc_coupled_local_response.py --domain pilot --case nf512_dt_0_0005 --budget-s 60
python scripts/lab.py scripts/derive_nsc_coupled_local_response.py --domain stored-frames --case nf512_dt_0_0005 --substeps 2 --budget-s 600
```

`--resume` returns a finished record only when domain, case, substeps,
source, and declared settings match. A different request is rejected and
does not rewrite the record. Sixteen tests passed in \(1.22\) seconds
after the appended window.

## Local-boundary successor

The JSON and NPZ hashes above remain the immutable v1 payloads. The v1
review is not rewritten. `results/development/nsc-local-boundary-review-v2.json`
binds those payloads and the current bytes of this note and of the consumer
module. The live finite-domain label remains initial mode columns 0 and 1,
weights \(0.75\) and \(0.75\), with no phase QR. The stored JSON keeps its
original observer string. This note adds no stress and no
\(1024\)-dimensional reintegration.
