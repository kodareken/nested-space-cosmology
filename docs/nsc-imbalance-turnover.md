# Mean circulation from an incoherent imbalance

This note records one executed calculation on the frozen six-mode window.
The preparation has no interregion coherence and no interregion channel
current. Unitary evolution still has a nonzero Cesaro mean current. The
geometry is an input. Nothing here chooses a radius, relaxes a state, or
closes the incoming gate.

The core is
[nsc_imbalance_turnover.py](../src/recursive_horizons/nsc_imbalance_turnover.py),
sha256 `565221169e8d188eb464928f51f980ae4d96974a635e0662b00110c24039ab32`.
Its `compute_result` was executed once, in about 15 seconds. An independent
manual cross-check agrees with the readings below. The saved
[record](../results/development/nsc-imbalance-turnover-v1.json) replays
34/34 checks. Its schema is `NSC-IMBALANCE-TURNOVER-v1` and its verdict is
`PASS_IMBALANCE_MEAN_CIRCULATION`. The physical local gate in that record
remains `OPEN`. An independent review accepted the 34 conditions, the
finite-time bound, and the forced optional-failure semantics. The driver is
[derive_nsc_imbalance_turnover.py](../scripts/derive_nsc_imbalance_turnover.py),
sha256 `8afc016f7299ff270a7cd3745dca7aadd8e68e40bcb5272ad85c0e802a17770f`.
The tests are
[test_nsc_imbalance_turnover.py](../tests/test_nsc_imbalance_turnover.py).

The stationary witness in [nsc-finite-turnover.md](nsc-finite-turnover.md)
remains valid. That state is prepared already carrying current. The
preparation here removes that presupposition. It does not replace the
earlier numbers.

## Window and preparation

The generator is the same frozen window as that witness,

$$
J=\texttt{finite\_window}(H,B,3/2,0,3),
$$

with \(H=\begin{pmatrix}1&i/5\\-i/5&2\end{pmatrix}\),
\(B=\begin{pmatrix}1/4&i/7\\1/9&1/6\end{pmatrix}\), and \(\hbar=1\).
Region labels start at zero. Diagonal blocks are \((3/2)^r H\), nearest
links are \((3/2)^r B\) and \((3/2)^r B^\dagger\), and the corner block is
absent. Those bytes are reused, not recomputed.

Let \(P_0\) and \(P_2\) be the projectors onto the first and last regional
blocks. The recorded family is the incoherent preparation

$$
C_0(\delta)=\frac I2+\delta(P_0-P_2),\qquad \delta\in\left\{+\frac14,-\frac14,0\right\}.
$$

Off-blocks are zero. The executed member is \(\delta=1/4\). Its eigenvalues
are \(3/4,3/4,1/2,1/2,1/4,1/4\), so \(0\le C_0\le I\). Regional populations,
meaning traces of the three diagonal blocks, start at \(3/2\), \(1\), and
\(1/2\). The trace is \(3\). The executed energy is

$$
\operatorname{Tr}(JC_0)=\frac{99}{16}.
$$

Both scalars are constant under \(C(t)=e^{-iJt}C_0 e^{iJt}\), and both agree
between \(C_0\) and the mean below. The three regional populations are not
separately constant. Their Cesaro values differ from \(3/2\), \(1\), and
\(1/2\).

## Commutant projection

On this \(n=6\) generator the core builds \(I,J,\ldots,J^{10}\) and solves
the Gram system

$$
G_{ab}=\operatorname{Tr}(J^{a+b}),\qquad
b_a=\operatorname{Tr}(J^a C_0),\qquad
a,b=0,\ldots,5,
$$

$$
G\alpha=b,\qquad
\bar C=\sum_{k=0}^{5}\alpha_k J^k.
$$

Gram entries are rational. A nonzero determinant means those powers are a
basis of the polynomials in \(J\). Moments may lie in \(\mathbb Q(\sqrt{29})\).
The solved \(\bar C\) matches the moments of \(C_0\) through degree 5 and
commutes with \(J\). For nondegenerate Hermitian \(J\), polynomials in \(J\)
are the commutant, so \(\bar C\) is the Cesaro mean

$$
\bar C=\lim_{T\to\infty}\frac1T\int_0^T e^{-iJt}C_0 e^{iJt}\,dt.
$$

That identification is an average. It is not a statement that \(C(t)\)
approaches \(\bar C\) pointwise.

## Current and sign

Embedded spectral projectors of the onsite \(H\) are labelled plus and minus
on each region. The core's channel current is

$$
j_{ab}=2\operatorname{Im}\operatorname{Tr}(P_a J P_b C).
$$

The core declares that formula as `j_ab=2*Im*Tr(P_a*J*P_b*C)`. It is the
block trace implemented by `_current`. Evolution gives
\(\dot C=-i[J,C]\) and \(\dot n_a=\operatorname{Tr}(P_a\dot C)\). With the
mode resolution of the identity and Hermitian \(J\) and \(C\),

$$
\dot n_a=\sum_b j_{ab},\qquad j_{ba}=-j_{ab}.
$$

A positive value adds occupation to \(a\) and removes the same occupation
from \(b\). It is inflow to \(a\) from \(b\). The finite-turnover core checks
the regional form of this identity as
`sign_derivative_matches_incoming_mode_sum`, using the same `_current`.
The channel token `0:plus->1:plus` only names the argument order \((a,b)\).
It is not a flow arrow from \(0+\) toward \(1+\).

Because \(C_0\) has zero off-blocks, every interregion \(j_{ab}(C_0)\) is
zero. The absent corner also carries no current, and neither does an
intraregion plus–minus pair: the onsite block is a multiple of \(H\). The
executed mean is not current-free. The signed current into \(0+\) from \(1+\) is

$$
j_{(0,+),(1,+)}=0.0012511140479781876.
$$

The arrow is into the region-0 plus mode, from the region-1 plus mode.
Because \(\bar C\) commutes with \(J\), the mean accumulation of every mode
is zero: the mean exchange is balanced. That balance is not a fixed
instantaneous regional population.

Executed readings of this core:

| Quantity | Reading |
|---|---|
| \(\delta\) | \(1/4\) |
| \(\operatorname{Tr}(C_0)=\operatorname{Tr}(\bar C)\) | \(3\) |
| \(\operatorname{Tr}(JC_0)=\operatorname{Tr}(J\bar C)\) | \(99/16\) |
| Initial regional populations, regions \(0,1,2\) | \(3/2\), \(1\), \(1/2\) |
| Cesaro regional populations, regions \(0,1,2\) | \(1.3837895948404104\), \(1.0058995585808204\), \(0.6103108465787691\) |
| \(j_{(0,+),(1,+)}\), inflow to \(0+\) from \(1+\) | \(0.0012511140479781876\) |

The three mean decimals are the executed readings. As written, they sum to
the conserved trace \(3\). Their exact rational strings are not copied here.
A channel current accounts for that mode-occupation
transfer. By itself it is not a net heat flow and not the energy flow. The
energy datum is the separate conserved scalar \(99/16\).

## Controls and the exterior probe

The core compares the same window at the other two values of \(\delta\).
At \(\delta=0\), \(C_0=I/2\) is already stationary and every mode current is
zero. At \(\delta=-1/4\), every adjacent mean current is the negative of the
\(\delta=+1/4\) current. Adjacent mean currents scale exactly with \(\delta\).
A mode-filling map, projecting one embedded onsite projector at a time, is
the core's linear check that those imbalance currents are reproduced. The
saved record marks `filling_map_reproduces_imbalance_currents` true. Forcing
that optional control to fail keeps the primary readings, sets the verdict
to `FAIL_IMBALANCE_TURNOVER`, and does not add an optional-pass flag.

The filtered response reuses `_filtered` at \(z=1+2i\), retained region 0:

$$
F\bar C F^\dagger,\qquad
F=\bigl[S^{-1},\; S^{-1}V(zI-J_{EE})^{-1}\bigr].
$$

This is a resolvent-weighted covariance. It is not the equal-time occupation
and not a stress. The exterior self-energy is not \(\bar C\).

The exterior probe adds \(1/8\) to \(J_{2,2}\) only. The retained block and
the coupling block stay fixed, and the probe keeps the same \(C_0\). The
mean is then recomputed by the same Gram projection. Holding the unperturbed
\(\bar C\) fixed instead is a different shortcut. Executed Frobenius changes
of the filtered response from the unperturbed baseline are

| State used with the probed generator | Frobenius change |
|---|---|
| Same \(C_0\), mean reprojected | \(0.004920298239642384\) |
| Frozen unperturbed mean | \(0.00030421076749148856\) |

## Finite-time bound

The core's bound, once a positive rational floor exists for an adjacent
mean current, is an averaging remainder. It is not a cooling law. The
docstring states

$$
\|C_T-\bar C\|_F\le \frac{2\|Y\|_F}{\gamma T}.
$$

Here \(Y=C_0-\bar C\) is off-diagonal in the eigenbasis, \(\gamma>0\) is a
rational lower bound on the eigenvalue gap, and

$$
C_T=\frac1T\int_0^T e^{-iJt}C_0 e^{iJt}\,dt.
$$

Each off-diagonal oscillation contributes at most \(2/(|\lambda_m-\lambda_n|T)\).
The value returned by the function encloses this with rational bounds.
\(\|Y\|_F^2=\operatorname{Tr}(Y^2)\) is required to be rational, and
`deviation_upper` is the integer-square-root upper bound of that square, hence
an upper bound of \(\|Y\|_F\). `gamma` is `gamma_lower_exact`: the minimum
separation of the six simple real isolating intervals of the monic
characteristic polynomial, after those intervals have been shrunk until they
are disjoint. For the selected adjacent channel,
\(K=-i(P_aJP_b-P_bJP_a)\) and

$$
\|K\|_F^2=\operatorname{Tr}(P_b J_{ba}P_a J_{ab})+\operatorname{Tr}(P_a J_{ab}P_b J_{ba}).
$$

`current_operator_upper` bounds \(\|K\|_F\). The square is first replaced by
a rational upper bound in \(\mathbb Q(\sqrt{29})\), and the same
integer-square-root step then bounds the square root.
`current_abs_lower` is a positive rational floor on \(|j_{ab}(\bar C)|\).
The returned statements are exactly

$$
\texttt{2*deviation\_upper/(gamma*T)}
$$

for the state, and

$$
\texttt{2*deviation\_upper*current\_operator\_upper/(gamma*T)}
$$

for the mean-current error. The sign of that selected mean current stays
separated from zero when

$$
T>\frac{2\cdot\texttt{deviation\_upper}\cdot\texttt{current\_operator\_upper}}{\gamma\cdot\texttt{current\_abs\_lower}}.
$$

The selected channel is the adjacent pair with the greatest positive floor,
not a channel chosen in this note. Those rational enclosures are not copied
here. The saved record carries them in `finite_time_bound`, with `obtained`
true and `finite_time_bound_error` null.

## What this does not claim

Instantaneous regional populations are not fixed. Only the total trace and
\(\operatorname{Tr}(JC)\) are the conserved scalars recorded here. Cesaro
averaging is not relaxation. A channel current is not, by itself, net heat
or energy flow. The filtered covariance is not a stress. The generator is
the frozen window, so a change in the state does not change the geometry.
Autonomous regeneration is therefore still missing. The calculation is not
a self-regulating radius, a thermal mechanism, a matter–antimatter or clock
identity, a cosmological history, or an eternity theorem. It does not
compare with LambdaCDM and does not refute it. It does not close the paused
incoming gate. That gate remains OPEN, and its campaign remains paused.

The research direction remains regeneration: sustained balanced turnover at
stable overall regional size and content. The next scientific blocker is
geometry feedback through the complete action, including the CW mapping. It
is not a new gate.

## Reproduction

From the repository root, replay the saved record with the validation
interpreter. `--record` is creation-only and refuses to overwrite an
existing file. Existing evidence is replayed by `--check`.

```sh
PYTHONPATH=lab/src .venv/validation/bin/python lab/scripts/derive_nsc_imbalance_turnover.py --check
PYTHONPATH=lab/src .venv/validation/bin/python -m pytest -q lab/tests/test_nsc_imbalance_turnover.py lab/tests/test_nsc_finite_turnover.py
```

This replay passed 34/34 source-bound checks. Thirteen focused tests passed:
seven for this calculation and six for the existing finite-turnover witness.
The record sha256 is
`7c86876f61b4af4b58e50da43ba12324b13154e830a1a0d8d9eed37721bc5e8e`.
It binds the core and driver above, the frozen window constructor
[nsc_nested_qualities.py](../src/recursive_horizons/nsc_nested_qualities.py)
at `02e445310e5f44a1cd1b972071cdf4f03c1a84ba12ccf041f9c3f55a621bddc1`, and
the frozen helpers
[nsc_finite_turnover.py](../src/recursive_horizons/nsc_finite_turnover.py)
at `778fd15a46dc421aa625e12aa27d16868f29b7d85430ed4c5555db4f5d744308`.
The finite-turnover record
[nsc-finite-turnover-v1.json](../results/development/nsc-finite-turnover-v1.json)
remains `2caa664b2c46994e503467e8f5703d801bfef17dca9c563174b8bb4915d3e382`.
The physical local gate remains `OPEN`.
