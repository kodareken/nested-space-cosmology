# Direct-order conformal ADM source

This owner is a raw finite canonical trace on one radial block. It is not
a renormalized stress, not the completed common action, and not a physical
matter branch \(V\). Those pieces are not supplied here. Active action
ownership is being audited independently. The historical
[vacuum-matched CTP prototype](nsc-vacuum-matched-ctp.md) is one recorded
prescription, \(B=\Gamma_{\mathrm{heat}}-\Gamma_{\mathrm{canonical\,sea}}\).
It is not declared mandatory for every future model, and this module does
not evaluate it.

The spherical feedback chart in
[nsc-spherical-feedback-action.md](nsc-spherical-feedback-action.md) is a
different domain: a bulk and boundary primitive with lapse and shift free.
Its matter-source branch is not coupled. The recorded finite imbalance in
[nsc-imbalance-turnover.md](nsc-imbalance-turnover.md) and
[nsc-imbalance-turnover-v1.json](../results/development/nsc-imbalance-turnover-v1.json)
is also a different domain: a six-mode mean circulation on a frozen window,
with geometry as an input. Neither result is this source, and this source
does not replace them.

The historical factorized owner
[nsc_adm_source.py](../src/recursive_horizons/nsc_adm_source.py), described in
[nsc-adm-source-constraints.md](nsc-adm-source-constraints.md), is unchanged.
Its sha256 is `0dba93a1090728c0c26fca2f04ac13ae74201fbb66b1aee09155a7146e5d2d09`
in both the lab copy and the publication copy. That file keeps the
finite-grid product rule \(\sqrt{N}\{P,1/q\}\sqrt{N}\). The rule is not
marked false. On a nonconstant grid it is a different ordering from the
direct anticommutator below.

## Operator and frame

Signature \(+---\), on the existing `CovariantStaticMetric` Fourier grid.
The canonical half-density is \(u=r\sqrt{q}\,\psi\). Physical fields and
conformal densities are

\[
N=rL,\qquad q=rQ.
\]

One positive angular block, with the existing momentum \(P\), is

\[
H=\sigma_2\otimes\{\mathrm{diag}(L/Q),P\}/2
+\sigma_1\otimes\mathrm{diag}(\kappa L)
-I\otimes\{\mathrm{diag}\,\beta,P\}/2.
\]

`direct_hamiltonian` builds that matrix. `conformal_source` returns the
same one-block matrix (`hamiltonian_is_one_block`) and weights only the
trace. The weight is \(M=4\kappa\) once, for that positive sector and its
identical isotropic copies. A separate \(-\kappa\) block is not added.
`conformal_source` rejects \(\kappa\le 0\). `block_source` is the same
trace without \(M\) and allows the massless angular control \(\kappa=0\).

At fixed \(L\) and \(Q\), this \(H\) does not depend on \(r\). The
historical product rule does, because \(\sqrt{N}\) does not commute with
\(P\) when \(N\) and \(q\) vary in space. Constant \(N=q\) still agrees:
on the 12-point cylinder the two Hamiltonians match to \(<10^{-12}\).

## Coupling API

```text
conformal_source(metric, shift, state, kappa,
                 state_interface="covariance", relative_name=None)
```

Schema `NSC-CONFORMAL-ADM-SOURCE-v1`. The state is held fixed. Nodal
values are partial derivatives of \(M\operatorname{Tr}(CH)\), not
coordinate densities. Keys are `L`, `Q`, `beta`, `N`, `q`, and `r`.
`energy` is the multiplied trace. `sector_energy` is the one-block trace.
`spacing` is the grid length per node.

If a constraint or a gravity momentum is a density per coordinate length,
divide the matching nodal partial by that node spacing. Do not insert
another angular factor. \(M=4\kappa\) is already inside `conformal_source`.
`block_source` and `ordering_comparison` omit \(M\); multiplying those
again by \(4\kappa\), or dividing a `conformal_source` force by \(4\kappa\),
would mistreat the copies.

The chain at fixed state, node by node, is

\[
F_N=F_L/r,\qquad F_q=F_Q/r,\qquad
F_r=-(LF_L+QF_Q)/r,
\]

and \(F_\beta\) is the same in both sets of variables. The direct ordering
makes the massless combination \(NF_N+qF_q+rF_r\) zero on the finite grid,
including nonzero angular \(\kappa\), because \(N\kappa/r=\kappa L\). That
is an ordering identity of this raw trace. It is not a trace anomaly and
not a renormalized stress.

`state_interface` is `covariance` or `relative_difference`. A covariance
must pass the existing fermionic check \(0\le C\le I\). A relative
difference must be Hermitian and must carry an explicit `relative_name`.
Nothing is subtracted unless the caller puts it in that matrix.
`vacuum_subtracted` stays false. The naive one-block spin identity
\(C=I/2\) has raw energy and nodal forces at roundoff zero. That does not
mean the renormalized source is zero, and it is not the recorded six-mode
\(C_0\). `frame_mapping` is `missing`. There is no six-mode embedding and
no full-source closure.

`prescribed_work_ledger` evolves \(C\) unitarily on a geometry written as
a function of time. The canonical identity used there is
\(d\operatorname{Tr}(CH)/dt=\operatorname{Tr}(C\,dH/dt)\) when
\(i\dot C=[H,C]\). The path is a control. It is not a solution of
\(F_A=0\), and it does not include a vacuum branch or a physical \(V\).

## Executed control

Cell: `smooth_metric(12, general=True)`, shift \(0.08\cos(2\pi x/L)\),
\(\kappa=1\). Fixed-\((L,Q)\) probe:
\(\exp(0.12\cos(4\pi x/L))\). The Gaussian \(C\) is the fixed unitary
sample in the module, not a recorded 6J state. The sea comparison is the
negative projector of this block's direct \(H\). The same matrix is used
on both orderings. Ten tests passed, in 0.17s and again in 0.22s:

```text
PYTHONPATH=lab/src .venv/validation/bin/python -m pytest -q lab/tests/test_nsc_conformal_adm_source.py
python scripts/lab.py -m pytest tests/test_nsc_conformal_adm_source.py -q
```

Direct order, both states: fixed-\((L,Q)\) Hamiltonian change and Ward
\(\ell_2\) are below \(10^{-12}\). Finite differences at nodes 0, 5, and
11 match \(L,Q,\beta,N,q,r\) to relative \(10^{-6}\) and absolute
\(10^{-7}\). Every one-block nodal norm is greater than \(10^{-3}\).
One-block \(F_r\) and \(F_\beta\) agree with the historical owner to
\(<10^{-12}\). The product-rule split is in \(F_N\) and \(F_q\).

Gaussian state, one block before \(M\):

| Quantity | Value |
|---|---:|
| Product-rule Ward \(\ell_2\) | 0.010872496849216692 |
| Physical \(\\|F_r\\|_2\) | 0.03308026372835154 |
| Product-rule fixed-\((L,Q)\) trace change | 0.00014262357291895732 |
| Product-rule \(\\|H\\|_F\) change | 0.0533295845381964 |
| \(\\|F_N^{\mathrm{direct}}-F_N^{\mathrm{product}}\\|_2\) | 0.021077411660050126 |
| One-block energy | -0.8700053565189911 |
| Coupled energy, \(M=4\) | -3.4800214260759645 |

Canonical sea of the direct block, same probe, still one block:

| Quantity | Value |
|---|---:|
| Product-rule fixed-\((L,Q)\) trace change | -0.0461283850516527 |
| Product-rule Ward \(\ell_2\) | 0.13566530208932837 |
| Physical \(\\|F_r\\|_2\) | 0.6496334743012157 |

These are the executed magnitudes. They are not a recalled pair of digits,
and they do not retract the historical factorization.

Prescribed ledger on the Gaussian state, \(M=4\), times
\((0,0.17,0.41,0.73)\), amplitude \(0.04\). Time 0 is the supplied metric.

| Quantity | Value |
|---|---:|
| Total work | 0.020566298598263436 |
| Energy change | 0.02056629859826309 |
| Ledger closure | \(-3.469446951953614\times 10^{-16}\) |
| Unitary residual max | \(2.6645352591003757\times 10^{-15}\) |
| Work residual max | \(1.824929096727601\times 10^{-15}\) |
| Force-versus-trace linear residual | \(4.1633442832782175\times 10^{-15}\) |
| Alignment of the prescribed step with \(F\) | -0.09978728262390682 |

The alignment is not steepest descent. A named relative matrix
`gaussian_minus_canonical_sea` reproduces source(Gaussian) minus
source(sea) at \(\kappa=2\). The bare Gaussian source is not that
difference.

## Bytes

| Path | sha256 |
|---|---|
| [nsc_conformal_adm_source.py](../src/recursive_horizons/nsc_conformal_adm_source.py) | `7e935b981de6cb4876142819a1a15587b81574ef1c15bf11adc62a7d8cdaeaea` |
| [test_nsc_conformal_adm_source.py](../tests/test_nsc_conformal_adm_source.py) | `27aa469b8bdec245597bd17a44f17d24244c629cbc2f5ad2628e0aeb43423cc9` |
| Historical [nsc_adm_source.py](../src/recursive_horizons/nsc_adm_source.py) | `0dba93a1090728c0c26fca2f04ac13ae74201fbb66b1aee09155a7146e5d2d09` |

No parent or child frame map is included. No common-action branch and no
physical \(V\) are included. The raw trace above is the whole claim of
this note.
