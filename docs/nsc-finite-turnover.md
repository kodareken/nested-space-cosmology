# Finite turnover on one inherited window

This note records a computed finite result for the author's fountain
proposal: an inherited law, balanced ongoing exchange, and a local response.
The values below are the saved record for one frozen window and one chosen
stationary state. That record replays 37/37 checks. The incoming-gate
campaign stays paused, and that gate stays OPEN. See the
[claim ledger](claim-ledger.md).

The geometry and the state are inputs. Nothing here selects a radius, a
temperature, or a continuing cosmological history.

## Owners

Paths are relative to this directory. The first four are the calculation
package; this note is the fifth. The frozen window constructor is reused
and is not one of these five.

| Role | Owner |
|---|---|
| Core | [nsc_finite_turnover.py](../src/recursive_horizons/nsc_finite_turnover.py) |
| Driver | [derive_nsc_finite_turnover.py](../scripts/derive_nsc_finite_turnover.py) |
| Tests | [test_nsc_finite_turnover.py](../tests/test_nsc_finite_turnover.py) |
| Record | [nsc-finite-turnover-v1.json](../results/development/nsc-finite-turnover-v1.json) |
| Note | [nsc-finite-turnover.md](nsc-finite-turnover.md) |

The window is the frozen `finite_window` in
[nsc_nested_qualities.py](../src/recursive_horizons/nsc_nested_qualities.py),
as specified in [nested qualities](nsc-nested-qualities.md).

## Fixed window and stationary state

Let \(H=\begin{pmatrix}1&i/5\\-i/5&2\end{pmatrix}\) and
\(B=\begin{pmatrix}1/4&i/7\\1/9&1/6\end{pmatrix}\). The operator is the
\(6\times 6\) matrix

\[
J=\texttt{finite\_window}(H,B,3/2,0,3).
\]

Region labels start at zero. Diagonal blocks are
\((3/2)^r H\), nearest-neighbor links are \((3/2)^r B\) and
\((3/2)^r B^\dagger\), and the corner block between the first and third
regions is absent. The maximum absolute row sum is \(379/70<8\).

The chosen state is

\[
C=\frac I2-\frac{J^3}{2048}=\frac I2-\frac{J^3}{4\cdot 8^3}.
\]

It is admissible in the same sense as the finite-window covariances,
\(0\le C\le I\), and it is stationary: \([J,C]=0\). Unitary evolution
therefore leaves \(C\) unchanged. The expected total occupation

\[
\operatorname{Tr}(C)=\frac{7246127239}{2477260800}
\]

is constant. No thermal ensemble was used to choose \(C\).

## Balanced exchange

Let \(P_a\) be the natural spectral projectors of the onsite \(H\), embedded
on the regional blocks. The mean channel current is

\[
j_{ab}=2\operatorname{Im}\operatorname{Tr}(P_a J P_b C).
\]

Inter-region currents are nonzero in both directions. Every onsite mode and
every region has net accumulation zero, so the signed currents into each of
those subspaces sum to zero. There is no direct current along a corner that
was not already in \(J\): a return stays on modes carried by the existing
nearest-neighbor links.

Across the first regional cut, \(G\) is the sum of the absolute mean channel
currents,

\[
G=\frac{8986745}{950450651136}.
\]

The one-way positive sum on that cut is \(G/2\). Because \(C\) is stationary,
\(G\) does not grow. The activity integral \(Gt\) still grows with the
evolution parameter. Stored occupation does not. \(Gt\) is not a count of
newly created particles.

## Filtered local response

Fix \(z=1+2i\). Split \(J\) into the first region \(A\) and its complement
\(E\), with coupling \(V\) from \(A\) into \(E\). The exterior self-energy,
Schur factor and retained resolvent row are

\[
\Sigma=V(zI-J_{EE})^{-1}V^\dagger,\qquad
S=zI-J_{AA}-\Sigma,
\]
\[
F=\bigl[S^{-1},\; S^{-1}V(zI-J_{EE})^{-1}\bigr].
\]

The filtered covariance \(FCF^\dagger\) equals the retained block of the full
resolvent covariance \((zI-J)^{-1}C\bigl[(zI-J)^{-1}\bigr]^\dagger\). The
cross blocks \(C_{AE}\) are kept. Dropping the memory kernel, the exterior
drive, or those cross correlations produces a different retained result.
\(FCF^\dagger\) is a filtered covariance at this complex frequency. It is
not the equal-time occupation and not a stress.

With \(J_{AA}\), \(V\) and the initial \(C\) held fixed, the outside-only
update \(J_{22}\mathrel{+}=1/8\) changes that filtered response by Frobenius
norm \(6.879146355717442\times 10^{-5}\). A farther diagonal probe
\(J_{5,5}\mathrel{+}=1/10\), again outside the retained block, changes it by
\(2.8630169396700523\times 10^{-7}\). The same \(J_{AA}\), \(V\) and \(C\)
stay fixed.

## Assumptions and scope

The finite depth, scale \(3/2\), frozen \(H\) and \(B\), and this \(C\) are
chosen data. The row-sum comparison with \(8\) only justifies the cubic
shift inside that choice. Evolution is the linear law of the fixed
Hamiltonian \(J\), in units \(\hbar=1\), on the inherited window. The local
blocks and the coupling are not refitted when the exterior block changes.
The coarse region graph remains a path; no corner link is added.

This is not a self-regulating size, a thermal mechanism, a
matter–antimatter or clock identity, a cosmological regeneration, or an
eternity theorem. It does not close the paused incoming gate and does not
identify the filtered covariance with a metric source.

## Reproduction

From the repository root, replay the saved record.
`--record` is creation-only and refuses to overwrite an existing file.
Existing evidence is replayed by `--check`.

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_finite_turnover.py --check
.venv/validation/bin/python scripts/lab.py -m pytest -q tests/test_nsc_finite_turnover.py
```

The saved record replays 37/37 checks. Six focused tests passed.
Its schema is `NSC-FINITE-CHANNEL-TURNOVER-v1` and its verdict is
`PASS_FINITE_CHANNEL_TURNOVER`. The physical local gate remains `OPEN`.
