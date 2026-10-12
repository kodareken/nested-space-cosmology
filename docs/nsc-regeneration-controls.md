# Regeneration controls on the saved spherical episode

Thin consumer of the owned Galerkin step and the saved \(T=0.05\) episode.
The column source, radius Newton, shift momentum, and `rk4_step` stay with
their owners. No added radius force, no mean subtraction, and no edit of the
episode driver. Coarse CPU was 302.92 s. The fine budget was not opened.
The v1 JSON and NPZ are immutable. The v2 file is an assessment of those
bytes, not a new evolution.

## Active note

\(G > 0\) is the initial convex-radius hypothesis for constant \(Q\),
\(\chi = p_r = p_\chi = 0\) and \(r = y^2\), with \(G\) independent of \(y\).
On that domain \(J(y)\) is strictly convex. The evolving chart is the owned
`chart_failure`: a finite state with \(r > 0\), \(Q > 0\) and \(L > 0\).
Lapse \(N = r L\) is positive when \(r\) and \(L\) are. \(\min G\) remains a
diagnostic. A dynamic step stops on `chart_failure`, checks the column Gram
as \(Q\) changes, and checkpoints the exit. A failed chart step is not kept.

The Galerkin generator is \(H_G = U_f^\dagger H_{\mathrm{fine}} U_f\) on
\(2 n_f\). On the \(n_f = 14\), \(n_q = 64\) control,
\(\tfrac12 I_{2 n_f}\) commutes with \(H_G\): the Frobenius commutator is 0
and the pullback gap against `apply_dirac` is \(7.11\times 10^{-15}\).
The \(2 n_q\) embedding of that identity has Frobenius commutator 11.08
with the unprojected \(H_{\mathrm{fine}}\). The fine-grid identity remains
an algebraic source control. Population control accepts six weights and
refuses a quadrature-sized occupation.

Window reversal of at least 0.10, a new leader that gains at least 0.10,
and a leader share of at least 0.50 are proxies. `candidate_regime` is
their disjunction with the stable-throughflow proxy. They are not programme
requirements. The programme allows a maintained structure with throughflow.
Splitting one drift is not complete renewal. Uniform six weights stay
localized relative to the empty complement. A four-window share below 0.50
is a separate proxy.

## Criterion, fixed before the evolutions

Normal-shell content is measured on the four fixed regional windows. The
leader is the largest share. Renewal, on the joined series from \(T=0\),
asks for an end state that stays localized and either a reversal of at least
0.10 in the original leader's shell content or a different window becoming
the leader after gaining at least 0.10. Two segments of one drift are not
complete renewal. Content, flux, reversal, and energy accounting are reported
separately.

A changed-source initial state is admitted when \(G\), \(r\), and \(Q\) are
positive, under the convex-radius hypothesis above. The column Gram stays
within \(10^{-8}\) of the identity, and the projected and full Hamilton
residuals are at most \(10^{-4}\) and \(10^{-2}\). A Newton stall above
\(10^{-10}\) is a solver floor, not a physical failure.

The 0.10 shell scale is the saved episode's leader change, about \(+0.108\).
Field energy, \(\chi\), and proper velocity are compared with floors 0.10,
1, and 0.01. Those floors sit on the saved coupled effect, not on a constraint
residual.

## Initial data

Same saved mode columns. Uniform occupations \((1/2)^6\) and reversed
\((1/4,1/4,1/2,1/2,3/4,3/4)\) recompute the source and solve radius and
\(p_Q\). At \(n_f=512\) both solves stall at a projected residual of about
\(6\times 10^{-9}\) with a full residual of about \(9\times 10^{-6}\). That is
the solver floor. \(G_{\min}=0.06316\), \(r\) and \(Q\) stay positive, and the
Gram gap is \(2\times 10^{-16}\). The current mean is about \(10^{-14}\) and
the shift-residual mean matches it to \(10^{-27}\). The current array is not
edited. The source does change: the density gap from the saved preparation is
13.76 (uniform) and 27.52 (reversed).

Uniform shell shares start at 0.474, below the 0.50 window proxy, split
between the first two windows. The six weights themselves remain the whole
occupied support; the quadrature complement carries occupation 0. Reversed
occupations move the leader to the second window, share 0.596.

## Controls against the saved coupled run

The saved baseline was not rerun. Frozen geometry uses that saved initial
state and the owned column rates, with geometry velocities and applied
fieldwork set to zero.

| Control, \(n_f=512\) | What changes by \(T=0.05\) | Against the saved coupled run |
|---|---|---|
| Frozen geometry | Shell leader \(+0.178\), field energy \(-3.2\times 10^{-8}\), \(\chi\) and proper velocity unchanged, applied work 0 | The \(0.412\) field-energy exchange, \(\Delta\chi_{\max}=11.15\), and proper-velocity change \(0.088\) are absent. The shell rearrangement is not. |
| Uniform | Share \(0.474\to 0.487\), field energy \(-0.547\), \(\Delta\chi_{\max}=11.50\), proper-velocity change \(0.094\) | The window-share proxy stays below 0.50. Field-energy exchange is \(0.135\) larger than the baseline. \(\chi\) and proper velocity stay inside the declared floors. |
| Reversed | Leader is window 1, share \(0.596\to 0.584\), field energy \(-0.681\), \(\Delta\chi_{\max}=11.57\) | Field-energy exchange is \(0.269\) larger. Proper velocity is \(0.011\) larger. The shell leader follows the reversed occupations. |

\(n_f=256\) reproduces these effects well inside one percent. Frozen geometry
keeps each column energy at its initial value and drives the full Hamilton
residual to about 24, so that trajectory is not a constraint solution. It is
the control that removes metric work. Coupled runs keep the total-energy
defect below \(10^{-6}\) of the field-energy exchange.

## Second episode

The state at \(T=0.05\) is the saved `nf512_dt_0_0005` final arrays, hash
`cd6b1f88a12635b137aaddc47444c06cc681c7568271ef989312f48f43a0649e`. Nothing is
reset. The same step runs to \(T=0.10\).

The four readings at \(n_f=512\), copied from the immutable series:

| Ledger | \(T=0.10\) |
|---|---|
| Content | Window 0 remains the leader. Share \(0.609\to 0.620\). Leader shell content \(+0.117\). Maintained. |
| Flux | Packet proper flux \(-5.534\), reservoir \(+5.534\), sum \(2.6\times 10^{-14}\). Balanced throughflow. |
| Reversal | 0 on the joined series from \(T=0\). One drift. Splitting it at \(T=0.05\) is not complete renewal. |
| Accounting | Field energy \(-0.396\). Total-energy defect \(5.0\times 10^{-8}\) of that exchange. |

Column energies stay on the original occupations. Proper clock advance on
this slice is 0.153, and the leader-window clock is 0.092, over coordinate
time 0.05.

\(\chi_{\max}\) grows from 11.15 to 89.31, \(Q_{\min}\) falls from 0.220 to
0.1248, and \(G_{\min}\) falls to 0.0156 while remaining positive as a
diagnostic. Proper velocity at the end is \([-0.102, 0.664]\). \(n_f=256\)
matches the leader-content change to \(6\times 10^{-10}\) and
\(\Delta\chi_{\max}\) to about 0.002.

The proxy `candidate_regime` is false: reversal is 0, the leader does not
change, and the shell content moves by more than 0.10, so the
stable-throughflow proxy is false as well. The continuation is still a
maintained structure with balanced throughflow. The \(\pm 5\%\) imbalance
perturbations were not run. The proxy does not authorize them.

## Successor

`nsc-regeneration-controls-v2.json` is hash-bound to the immutable v1 files.

| File | sha256 |
|---|---|
| `nsc-regeneration-controls-v1.json` | `3abe72de1be7f5090569086a3134664727e25849eecada6f12b89040bdddc922` |
| `nsc-regeneration-controls-v1.npz` | `8509a1d0a8151c189ce3346d7d47f7c4a6453065d39ec01fda35515af037c96b` |

No dynamics, action, or force was added. The preserved endpoint is
\(T=0.10\), \(Q_{\min}=0.12477163551088903\), leader shell
\(+0.11695656049787839\), packet flux \(-5.534007341562886\), field-energy
exchange \(-0.395776510604783\).

On the stored final arrays the owned `subspace_frequency` is
\(\omega = 2998.824498749641\). With the recorded \(dt=5\times 10^{-4}\),
\(dt\,\omega = 1.49941\). That is above the owned cap 1.4 and below the
absolute RK4 limit \(2\sqrt{2}\). `stable_timestep` is
\(4.668495940938623\times 10^{-4}\). `chart_failure` is empty. The column
Gram gap is \(1.30694\times 10^{-9}\). The series maximum of the stored
unitarity is \(1.30695\times 10^{-9}\).

The lapse-peak reconstruction \(\omega \approx (L_{\mathrm{peak}}/Q_{\min}
+ |\beta|_{\mathrm{peak}}) k + m_{\mathrm{peak}}\) matches that endpoint
frequency to relative gap \(6.12\times 10^{-4}\). It crosses \(dt\,\omega =
1.4\) at the stored sample \(T=0.0965\), \(Q_{\min}=0.1345054206385445\).
There the leader shell change is \(+0.10855\), the balanced flux magnitude
is 5.408, reversal is 0, and the field-energy change from \(T=0.05\) is
\(-0.371\). The same four readings hold as at \(T=0.10\).

## Next question

The stored endpoint is already past the owned stability cap, and the sample
where the lapse-peak diagnostic crosses that cap does not change the four
ledgers. The smallest probe that can change a ledger is one owned RK4 step
from the stored \(T=0.10\) state at `stable_timestep`
\(4.6685\times 10^{-4}\), checkpointing `chart_failure`, the Gram, and the
four ledgers. The measured cost of one \(n_f=512\) continuation step is
0.483 CPU seconds. Another full \(\Delta T=0.05\) is estimated at 47.3 CPU
seconds. That one step is not run in this record.

## What this does not claim

The shell pattern is one maintained drift that carries balanced throughflow.
It is not a completed renewal. Frozen geometry shows that the field-energy
exchange and the proper motion require the coupled metric step, while the
shell rearrangement does not. This record does not claim a continuum limit,
a cosmological regeneration, or a closed constraint bound. It does not
replace the v1 evolution.

```sh
python scripts/lab.py scripts/derive_nsc_regeneration_controls.py --check
python scripts/lab.py scripts/derive_nsc_regeneration_controls.py --assess
python scripts/lab.py -m pytest tests/test_nsc_regeneration_controls.py -q
```
