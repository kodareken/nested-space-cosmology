# Regeneration realization audit

Independent reading of the completed finite chain against the stored
episode, the conditional local response, and the initial-constraint
certificate. No trajectory was integrated. No source array was edited.
Global eternity, a spherical continuum limit, and a ΛCDM fit are outside
this chain and are not scored.

The episode record's `goal_complete: true` is its own maintained-episode
flag. It is not this audit's verdict.

## Verdict

The maintained chain is measured from \(T=0\) through the arc crossing at
\(T=0.085\). Episode 1 supplies the bitwise Cauchy data of episode 2.
Window 0 stays the leader, the packet flux is cancelled by the reservoir,
and the field exchange requires the geometric velocity. The fixed \(T=0\)
observer then moves, and memory, the exterior drive, and the actual cross
block each move it by more than the reduction error. Renewal is false.
That does not fail the maintained branch.

An evolved error estimate at the scale of the claimed observables remains
open: `evolution_error_bound` is null. The initial radius distance and the
Hamilton residual have different units. Neither a comparison between them
nor demanding that the evolving trajectory stay within its initial radius
uncertainty would supply that estimate.

## Requirement to evidence

| Requirement | Status | Evidence |
|---|---|---|
| Bitwise handoff into episode 2 | PROVED | All eight Cauchy arrays of each regeneration run at \(T=0.05\) equal the predecessor `dt=0.0005` final. Hashes `cd6b1f88a12635b137aaddc47444c06cc681c7568271ef989312f48f43a0649e` (\(n_f=512\)) and `f24f5dbfb5591e2943efc4d60eef0b09f573643de5ea868cea2872d46d222441` (\(n_f=256\)). The half-step finals are not the handoff. |
| Two connected episodes, no reset | PROVED | Frames start at the handoff. The run stops at \(T=0.085\) on `BRIDGE_ARC_ENTERED_PACKET`, short of the requested \(T=0.15\). Joined reversal is 0. |
| Positive chart on the stored samples | MEASURED | On all four runs, \(r\ge 4.21235\), \(Q\ge 0.16260\), static \(L\ge 0.19114\), and \(G\ge 0.02644\). This is the sample chart, not a continuum enclosure. |
| Admissible covariance \(C=\Phi W\Phi^\dagger\) | MEASURED | Nonzero eigenvalues reproduce the weights \((0.75,0.75,0.5,0.5,0.25,0.25)\) within \(2\times 10^{-9}\). Gram gap \(\le 1.14\times 10^{-9}\). Trace is \(3\) within \(10^{-8}\). Eigenvalues stay in \((0.25,0.75)\). |
| Timestep safety | MEASURED | Largest \(dt\,\omega\) is \(1.20101\), below the cap \(1.4\) and below \(2\sqrt{2}\). The stored product matches \(dt\times\omega\) with gap 0. No step was rejected. |
| Declared source and \(T=0\) observer | PROVED | Occupations are those six weights. The observer is predecessor columns 0 and 1. Its Gram gap is \(2.2\times 10^{-16}\). It is not the \(\Phi(0.05)\) pair. |
| Window ledger, work, and coordinate split | MEASURED | \(\Delta E-\int(F_{\mathrm{boundary}}+W_{\mathrm{pressure}}+W_{\mathrm{lapse}})\) is at most \(1.54\times 10^{-7}\), under \(6\times 10^{-7}\) of the field exchange \(-0.28608\). Matter plus gravity equals the coordinate window energy within \(5\times 10^{-15}\). Packet plus reservoir flux ends below \(5\times 10^{-13}\). Coordinate work tracks the field exchange within \(10^{-4}\) of that exchange. Pressure work is \(-5.91\times 10^{-4}\), present and smaller than the exchange. Field plus gravity cancels to the same relative size. |
| Space and time agreement of those effects | MEASURED | Across the four runs the relative movement is at most \(2.6\times 10^{-5}\) for the proper-velocity change and at most \(1.5\times 10^{-6}\) for \(\chi\), field energy, leader content, and end flux. |
| Generated geometry versus frozen geometry | MEASURED | The frozen step keeps the column rates and sets geometry velocity and fieldwork to zero. The coupled field change \(-0.41230\) falls to \(-3.24\times 10^{-8}\). \(\chi\) and proper velocity do not move. Shell content still rises by \(0.17849\). The final frozen spinor differs from the initial spinor by \(0.217\). The frozen Hamilton maximum reaches \(24.230\). The control isolates the geometric exchange. It is not a constraint solution. |
| Uniform and reversed sources, same chart | MEASURED | At \(n_f=512\) the initial strong residuals are \(8.87\times 10^{-6}\) and \(8.86\times 10^{-6}\), with \(G>0.048\) and field changes \(-0.54691\) and \(-0.68145\). At \(n_f=256\) those initial residuals are \(1.37\times 10^{-3}\) and \(1.61\times 10^{-3}\). The coarse grid does not meet the fine initial level. |
| Maintained structure and throughflow | MEASURED | Leader stays window 0. Share moves \(0.608704\to 0.616530\). Leader content rises by \(0.081175\). End packet flux is \(-4.987559\). Renewal is false and the joined reversal is 0. |
| Arc crossing ends the episode | MEASURED | Recomputed on `nf512_dtmax_0_00025`: at \(T=0.05\) the both-positive arc starts after \(x=4\); at \(T=0.085\) it runs from \(x=3.96875\) to past \(x=5.9\). The same edge is stored on all four runs. The stop is that entry into \([0,4)\), not an imposed time cut. |
| Realized rate jet and its projection | MEASURED | Stored indicator equals the one-step increment divided by its \(dt\), gap 0. The \(0.005\) secant of \(\dot Q\) differs from that indicator by more than \(1.5\). Pullback of the prolonged coarse rate returns it within \(1.5\times 10^{-13}\). The indicator is not a curvature bound. |
| Weyl factor on the stored scalar | MEASURED | `weyl_C2` equals \(\chi^2/(3r^4)\) within \(4\times 10^{-16}\). The \(3\) is the owned product \((R_h-2)^2/(3r^4)\) with \(R_h\) replaced by \(\chi+2\). |
| Independent metric curvature | INCOMPLETE | The episode records `metric_curvature_computed_here: false`. The stored jets are not inserted into that scalar. |
| Local occupation, full versus reduced | MEASURED | On the stored linear \(Q(t)\), the autonomous occupation changes by \(0.267489\). Conditional full differs from it by \(1.16662\times 10^{-7}\). Streamed retained differs from conditional full by \(2.57030\times 10^{-4}\). No stress and no force are claimed. Renewal is false. |
| Memory and exterior drive | MEASURED | Recomputed from the stored occupation series: memory \(0.069528\), drive \(0.155260\). Both exceed the reduction error. |
| Actual cross omission | MEASURED | Recomputed \(C_{AE}=0.442238710589459\). \(C_{AA}\) eigenvalues are \(0.419397830429754\) and \(0.647447915556798\). The NPZ full-window series is `covariance_retained - covariance_without_cross`: its diagonal reproduces occupation movement \(0.204342\) and its Frobenius norm reaches \(0.235132\). The explicit probe diagonal has maximum \(0.087971\) on \([0.05,0.06]\). The [independent full split](nsc-coupled-local-response.md#independent-control-replay) gives movement \(0.204373\), with streamed occupation error \(1.196\times10^{-4}\). The sealed \(4.863\times10^{-16}\) is superposition/assembly roundoff, not reduction error. |
| Unique positive critical point of the initial lapse | PROVED | Conditional on the sealed enclosure \(G\ge 0.053733032411767\) and \(Q^2=0.063164222082737\). The direct-method theorem returns a unique smooth positive critical point. The sealed field `minimizer_existence_proved` is still null. This audit did not recompute the numerator enclosure. |
| Initial radius distance | MEASURED | Sealed certificate \(4.060849694930819\times 10^{-6}\). It is an initial-chart bound. It was not re-enclosed here and it is not propagated. |
| Initial shift mean | MEASURED | Recomputed mean lies in \((4.238643756404872\times 10^{-17},\,4.238643756404882\times 10^{-17})\) and does not contain 0. Symmetry defect is 0. The source array is unchanged. |
| Same-action \(p_r=\lambda r'\) profile | MEASURED | \(\lambda\in[-1.203880937969956\times 10^{-15},\,-1.203880937969953\times 10^{-15}]\). Critical-point leftover \(\le 3.857333633418152\times 10^{-21}\). The profile is not installed: initial \(p_r\) is 0, and by the handoff \(\lvert p_r\rvert\) is already \(6.135\). |
| Evolved error at the observable scale | INCOMPLETE | `evolution_error_bound` is null. On `nf512` the nodal Hamilton maximum grows from \(2.826\times 10^{-5}\) at \(T=0.05\) to \(1.396\times 10^{-4}\). The geometry-band projection reaches \(1.234\times 10^{-4}\) and the held-out complement \(3.034\times 10^{-5}\). On `nf256` the nodal maximum is \(9.241\times 10^{-3}\). Endpoint window smearing stays about \(10^{-8}\) to \(10^{-7}\). These are constraint diagnostics. Converting them into geometry or occupation error requires a derived estimate with the appropriate norm and units. |

The recomputed relations agree with the stored arrays. The Gram-identity flag in the controls report is a column-Gram statement. The covariance spectrum was recomputed separately and matches the weights.

## What the maintained chain uses

Episode 1, \(T=0\) to \(T=0.05\), already carries a field exchange of about
\(-0.412\). Episode 2 continues that state. Over the next \(0.035\) of
coordinate time the primary run exchanges a further \(-0.28608\) between
field and gravitational energy, keeps window 0, and stops because the
both-positive arc has crossed \(x=4\) into the packet. The four resolution
pairs agree on that exchange, the content change, and the flux.

Freezing the geometry removes the exchange and the growth of \(\chi\) while
the shell still rearranges, so the exchange is not radial bookkeeping of a
fixed chart. Changing the initial occupations to uniform or reversed, and
keeping the same chart and source formula, still produces a positive chart
and a field exchange. The local response then reads the generated \(Q(t)\)
as a prescribed schedule. Its observer stays the \(T=0\) pair. The
occupation movement and the three omissions are larger than the difference
between the conditional full series and the saved autonomous spinor.

## Constraint scope

The initial lapse statement is a proved existence result sitting on a
sealed distance. It says the saved positive \(y=\sqrt{r}\) is within
\(4.061\times 10^{-6}\) in \(r\) of the unique positive critical point of
\(J\), if the sealed lower bounds on \(G\) and \(Q^2\) are granted. The
shift statement is separate. A periodic \(p_Q\) cannot cancel a nonzero
mean current. The recomputed mean is about \(4.239\times 10^{-17}\). The
existing action already contains the profile \(p_r=\lambda r'\) that
removes that mean at the saved radius, with a critical-point leftover
below \(4\times 10^{-21}\). The evolved model does not use that profile.

A high-Fourier cutoff would not close the evolved constraint. At \(n_f=512\)
the geometry-band projection is the same order as the nodal maximum, and
the held-out complement is smaller. The unresolved piece is in the band
that the evolution keeps. The window smearing is a different functional:
it is small at the stored endpoints, and it is not `evolution_error_bound`.

## Remaining primitive

The open object is the field the records already leave null,
`evolution_error_bound`: an evolved estimate at the scale of the claimed
geometry, transfer, and local-response changes through \(T=0.085\).
The coarse and fine nodal Hamilton maxima differ by a factor of about
\(66\). Each keeps its own numerical domain; the coarse case need not have
the fine case's initial error. At \(n_f=512\) the held-out complement is
the smaller piece, so a high-Fourier-tail estimate alone would leave the
in-band diagnostic unaccounted for. The helper now also provides physical
source normalizations and proper-volume/coordinate smearings. Those are
constraint budgets, not radius or occupation error bounds.

The shift profile \(\lambda r'\) is quantified and is not installed.
Installing it would not by itself propagate the initial error. Until an
observable estimate exists, the maintained chain remains a measured finite
episode with refinement indicators and reported constraint budgets.

## Hashes and what this audit did not do

| Object | sha256 |
|---|---|
| `nsc-regeneration-episode-v1.npz` | `17364df24efead5814c3ab2a2ab5f82d9e09bf412c8c82e68d57a3b38498d7a1` |
| `nsc-spherical-feedback-episode-v1.npz` | `2864d3a8d3ba413e840c395961f14a3d71eabee702c27fd9f47c7d333eb24b77` |
| `nsc-coupled-local-response-v2.json` | `5330ccef1ba2fc4d80e0e0bdc0a98b1117e851f3873a46a0e057d63aae7457bb` |
| `nsc-coupled-local-response-v2.npz` | `84c09302413d48c2f5dfaba32090b8169d0a4c1dbbe9dfc7a971debaaa525d14` |
| `nsc-regeneration-controls-v1.npz` | `8509a1d0a8151c189ce3346d7d47f7c4a6453065d39ec01fda35515af037c96b` |
| `nsc-spherical-cauchy-error-v1.json` | `e0a37e1177b2a483eca0eb3e44243893bbdc04667496b732dd7e63467a97bdd2` |
| `nsc_spherical_episode_assessment.py` at this reading | `b1839b42d2a4841df0fa2231523eb787d1f45e581cfbd3b0f25c450915b331a5` |

At that helper hash, `inspect_saved_episode` still reports the predecessor
payload as missing a realized rate. This audit did not call `assess` on
the four regeneration runs. Feeding the one-step indicator in as
\(\partial_t\dot Q\), or feeding \(\chi+2\) in as \(R_h\), would manufacture
a metric completion the episode does not claim. A later edit of that helper
is outside this hash.

The \(65536\)-node lapse enclosure was not rerun. The shift enclosure was
rerun and left the sealed records unchanged. No RK4 step, no new source
solve, and no omission evolution was started.

`lab/tests/test_nsc_regeneration_realization_independent.py` passed, nine
tests in \(1.98\) s.
