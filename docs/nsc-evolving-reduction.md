# Retained-region reduction with retarded memory

This note records one finite numerical control for a time-dependent
retained-region reduction. The sine schedule below is a prescribed method
control. It is not an autonomous trajectory, not regeneration, and not a
coupled geometry. No actual \(H(g(t))\) was consumed.

The incoming-gate campaign stays paused, and that gate stays open. These
grid discrepancies are observed truncation errors. They are not certified
bounds.

## Owners

Paths are relative to this directory.

| Role | Owner |
|---|---|
| Core | [nsc_evolving_reduction.py](../src/recursive_horizons/nsc_evolving_reduction.py) |
| Dense tests | [test_nsc_evolving_reduction.py](../tests/test_nsc_evolving_reduction.py) |
| Streamed-backend tests | [test_nsc_evolving_reduction_streamed.py](../tests/test_nsc_evolving_reduction_streamed.py) |
| Note | [nsc-evolving-reduction.md](nsc-evolving-reduction.md) |

The frozen window and correlated Gaussian state are reused from
[finite turnover](nsc-finite-turnover.md). The occupied, empty, retarded,
and initial-cross kernels use the same ordering as
[the boundary state](nsc-boundary-state.md). Neither owner is modified here.

## API

`evolve_retained_region` takes a callable Hermitian \(H(t)\), one fixed
orthonormal basis \(V\), a uniform time grid, and either initial columns
or a full one-body covariance. Column weights in \([0,1]\) are optional.
A later rank-6 source is those columns. Building the full map is a
separate path.

The frame is \(T=[V,\mathrm{null\_space}(V^\dagger)]\). QR is not applied
to \(V\), because that rephases the observer. A callable or stacked basis
is rejected: a moving projector is unsupported without a Berry term.

In that fixed frame the blocks are \(A=H_{AA}\), \(B=H_{AE}\), and
\(E=H_{EE}\). The free exterior propagator solves

\[
i\dot W(t)=E(t)W(t),\qquad W(t_0)=I,
\]

by midpoint matrix exponentials. `propagator_substeps` refines each output
interval. The reducer does not integrate the full system. The independent
DOP853 sample lives only in the test.

Retained columns \(X\) and the history accumulator

\[
Z(t)=\int_{t_0}^{t} W(s)^\dagger B(s)^\dagger X(s)\,ds,\qquad Z(t_0)=0,
\]

obey

\[
\dot X=-iAX-BWZ-iBWY_0.
\]

On a uniform step \(h\), with \(K=W^\dagger B^\dagger\), the trapezoid step
solves only the retained matrix

\[
\begin{aligned}
&\Bigl(I+i\frac h2 A_{\mathrm{new}}+\frac{h^2}4 B_{\mathrm{new}}W_{\mathrm{new}}K_{\mathrm{new}}\Bigr)X_{\mathrm{new}}\\
&\quad=X_n+\frac h2\Bigl[f_n-B_{\mathrm{new}}W_{\mathrm{new}}\Bigl(Z_n+\frac h2 K_n X_n\Bigr)-iB_{\mathrm{new}}W_{\mathrm{new}}Y_0\Bigr].
\end{aligned}
\]

\(f_n\) is \(\dot X\) at the left endpoint. \(Z\) is the semiseparable
memory sum. It is not a dense kernel and not a relabelled full Hamiltonian
evolution.

Identity columns, or the covariance path, return the retained map \(R\).
The same-time check is \(\|RR^\dagger-I\|\). For a covariance \(C\),

\[
\begin{aligned}
C_A={}&A_{\mathrm{resp}}C_{AA}A_{\mathrm{resp}}^\dagger
+V_{\mathrm{resp}}C_{EE}V_{\mathrm{resp}}^\dagger\\
&+A_{\mathrm{resp}}C_{AE}V_{\mathrm{resp}}^\dagger
+V_{\mathrm{resp}}C_{EA}A_{\mathrm{resp}}^\dagger.
\end{aligned}
\]

\(A_{\mathrm{resp}}\) and \(V_{\mathrm{resp}}\) are the parent and child
responses in the boundary-state reconstruction. The last two terms are the
initial cross covariance. Dropping them does not change \(X\).

## Three contributions that stay separate

| Piece | Where it is kept | What it is not |
|---|---|---|
| Initial drive | Exterior amplitude \(Y_0\) in \(\dot X\) | Not the history \(Z\) |
| Initial cross covariance | \(C_{AE}\) and \(C_{EA}\) in the four-term assembly | Not exterior noise |
| Retarded memory | Accumulator \(Z\) and \(\Sigma^R=-i\theta B(t)W(t)W(s)^\dagger B(s)^\dagger\) | Not a full-system ODE |
| Occupied exterior kernel | \(F(t)C_{EE}F(s)^\dagger\), \(F=BW\) | Not added again on top of \(RCR^\dagger\) |
| Empty exterior kernel | \(F(t)(I-C_{EE})F(s)^\dagger\) | Not the occupied kernel |

\(\theta(0)=1/2\), matching the boundary-state owner. That equal-time value
does not enter the history integral. `memory=False` leaves \(Z=0\).
`outside_drive=False` omits \(Y_0\) from \(\dot X\) but still records the
supplied initial exterior amplitude. `coupling=False` zeroes \(B\) only.
`drop_cross_covariance=True` zeroes the two cross terms in \(C_A\).

The callable contract for a future geometry is

```python
evolve_retained_region(lambda t: operator(g(t)), local_basis, times, columns)
evolve_retained_region(
    lambda t: operator(g(t)),
    local_basis,
    times,
    columns,
    backend="streamed",
)
```

`backend="dense"` is the default and stores \(W(t)\). `backend="streamed"` is the opt-in column transport below. This module does not build \(g(t)\). None was passed, and no generated \(H(g(t))\) was compared. The active spherical loop already evolves geometry; that comparison belongs to a later trajectory consumer.

## Method control

The control reuses the inherited six-mode window \(J\) and the correlated
Gaussian \(C=I/2-J^3/2048\), which satisfies \(0\le C\le I\). The retained
basis is the first region, two coordinate columns, left unphased. The
prescribed generator is

\[
H(t)=J+0.2\sin(2t)\,W,
\]

with \(W_{22}=1\), \(W_{02}=W_{20}=0.1\), and every other entry of \(W\)
equal to zero. Indices are zero-based. This \(H(t)\) is a fixed input.
It does not regenerate a history, and it is not the output of a coupled
trajectory.

Negative controls on this same prescribed \(H\) omit memory, omit the
exterior initial drive, drop the initial cross covariance, or set the
coupling block to zero. Their separations below are observed distances,
not bounds.

## Observed comparison

The test integrates \(i\dot\psi=H(t)\psi\) with DOP853, outside the
reducer, on the same initial columns. Frobenius discrepancies for
\(t\in[0,1]\) and identity columns:

| Steps | \(\|X-V^\dagger\psi\|\) | \(\|C_A-C_A^{\mathrm{oracle}}\|\) |
|---:|---:|---:|
| 32 | 0.000794511018284106 | 0.00004041021196050106 |
| 64 | 0.0001987203539841231 | 0.000010117337491779154 |
| 128 | 0.00004968587992978299 | 0.000002530259501193949 |

The maximum over each grid equals the final-time amplitude discrepancy.
Successive ratios are about four, which is the observed second-order
behavior of the midpoint exponential plus the trapezoid. That pattern is
not a certified error bound. The trapezoid residual of the reduced step
and the four-term assembly residual against \(RCR^\dagger\) were both near
\(10^{-16}\).

On the 128-step grid the same-time CAR defects \(\|RR^\dagger-I\|\) were
observed at \(8.371778439790898\times 10^{-5}\), \(2.0960285419146274\times
10^{-5}\), and \(5.242001100107288\times 10^{-6}\) for 32, 64, and 128
steps. They track the discretization and are not a roundoff identity.

Ablation distances on the 128-step grid, from the full Volterra solution:

| Control | Observed separation |
|---|---:|
| Memory omitted | 0.051700309986861086 |
| Exterior initial drive omitted | 0.33489209407248566 |
| Coupling block set to zero | 0.33881278554664335 |
| Initial cross covariance dropped | 0.0006662699405864039 |

Each separation is larger than the corresponding 128-step discretization
error. With the coupling removed, the reduced solution still matched a
decoupled DOP853 sample to \(4.332142371574326\times 10^{-5}\). Omitting
memory left the stored history at norm zero. The full history norm was
zero at the initial time and \(0.33889957128313775\) at the final time.
The initial exterior amplitude remained present, with norm \(2\).

Final four-term norms on that grid were parent \(0.664656777695798\),
child \(0.0467113559240982\), each initial cross \(0.0005467079943855282\),
and total \(0.703671592062021\). At the initial time the child response is
zero, so the cross terms do not yet contribute to \(C_A\).

Exterior kernels at the final time and the midpoint were occupied
\(0.05232107508861198\), empty \(0.054026706900863294\), retarded
\(0.10634232610583307\), and initial-cross \(0.0005587354736209597\).
Occupied and empty differ. They match the fixed-frame boundary-state
kernels and are not a second copy of \(C_A\).

## Streamed column backend

`backend="streamed"` solves the same retained history. It does not replace that history by a full-system ODE. With \(U_E(t,s)\) the exterior propagator of \(E(t)\),

\[
\dot X=-iAX-\int_{t_0}^{t}B(t)U_E(t,s)B(s)^\dagger X(s)\,ds-iB(t)U_E(t,t_0)Y_0.
\]

The causal kernel in the integral is \(B(t)U_E(t,s)B(s)^\dagger\). Its source cohort is the set of columns of \(B^\dagger\). Adjoint midpoint transport to \(t_0\) gives \(G(t)=U_E(t,t_0)^\dagger B(t)^\dagger\), and

\[
B(t)U_E(t,s)B(s)^\dagger=G(t)^\dagger G(s).
\]

`causal_memory_kernel` returns that product and the two cohorts. It does not form \(W\). Occupied, empty, and initial-cross kernels are then \(G(t)^\dagger C_{EE} G(s)\), \(G(t)^\dagger(I-C_{EE})G(s)\), and \(G(t)^\dagger C_{EA}\). `exterior_kernels_from_cohorts` keeps \(\theta(0)=1/2\). These are the same blocks as \(F=BW\), because \(F(t)^\dagger=G(t)\).

Each output step advances only the current memory image \(S=UZ\), the drive columns \(UY_0\), and the source cohort, using `expm_multiply` on the same midpoint product as the dense propagator. The generator at one time may be the dense block \(E(t)\) or a `LinearOperator`. The multiply is that product's action, not a commutator-free or averaged surrogate. The returned `history` is still the interaction-picture accumulator \(Z\), rebuilt by the adjoint product on the new cohort. `memory=False`, `outside_drive=False`, `coupling=False`, and `drop_cross_covariance=True` keep the dense meanings, including the recorded initial exterior amplitude when the drive term is omitted.

The output grid stores \(A(t_n)\), \(B(t_n)\), \(X(t_n)\), and \(Z(t_n)\). \(Z\) has shape `(times, exterior_dimension, columns)`. The allocation record sets `time_indexed_exterior_propagator_bytes` to zero. Replaying \(Z\) from \(t_0\) grows quadratically with the number of output nodes. That cost is recomputation, not a stored \(W(t)\), and it is not a large-exterior timing claim.

The constructor still completes the observer with one dense `null_space`. A dense callable \(H(t)\) still forms the current exterior block at a midpoint and drops it. Those two arrays are not a time-indexed propagator. The fixed columns of \(V\) are copied into the frame unchanged.

No coupled trajectory is claimed. The incoming-gate campaign stays paused. The streamed tests are methodology checks on prescribed controls.

## Reproduction

From the repository root, the dense control and the streamed backend are

```sh
python scripts/lab.py -m pytest tests/test_nsc_evolving_reduction.py -q
python scripts/lab.py -m pytest tests/test_nsc_evolving_reduction_streamed.py -q
```

The dense file reported 13 passed. That comparison is the table above.
The streamed file reported 11 passed. It checks dense agreement, the
independent DOP853 sample, the omission controls, and the absence of a
time-indexed exterior propagator. Neither run closes the paused incoming
gate or identifies this reduced covariance with a metric source.
