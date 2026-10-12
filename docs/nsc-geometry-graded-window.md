# Geometry-graded window of the initial three-packet compression

The finite-window theorem in [docs/nsc-nested-qualities.md](../../docs/nsc-nested-qualities.md)
is applied to blocks measured from the active initial packets. The packets are
`prepare_rank6`: two parity lobes on \([n,n+1]\) and \([n+1,n+2]\), the recorded
carrier, and one minus-column phase, then the owner Löwdin step. The operator is
the owner conformal block

\[
H=\sigma_2\{L/Q,P\}/2+\sigma_1\kappa L-\{β,P\}/2,
\]

with the supplied gauge \(L=b_0\Omega^s\), \(β=β_0\Omega^s\), \(\Omega=3/2\), and
\(s=x\) on the packet arc \([0,4]\). Initial \(Q=b_0/a_0\) is constant up to
quadrature prolongation at \(10^{-14}\). No trajectory is integrated. The
historical rational link is not a target.

Reproduce from the repository root:

```text
python3 scripts/lab.py scripts/derive_nsc_geometry_graded_window.py --check
python3 scripts/lab.py -m pytest tests/test_nsc_geometry_graded_window.py -q
```

The record is `results/development/nsc-geometry-graded-window-v1.json`. The
writing run used under 4 s of numerical CPU. The record stores that time.

## Derived relation

On the nf=512 prolonged columns the measured window satisfies

\[
J_{nn}=\Omega^n H_{\mathrm{derived}},\qquad
J_{n,n+1}=\Omega^n B_{01},
\]

with defects \(5\times 10^{-14}\). The block coefficient in the published
recurrence, including its factor \(1/\Omega\), is this \(B_{01}\). Removing a
factor \(\Omega^{n/2}\) from each regional frame leaves the constant link

\[
B=\Omega^{-1/2} B_{01}.
\]

That normalized \(B\) is not a substitute for the block coefficient: putting it
into \(J_{n,n+1}=\Omega^n B\) moves the top resolvent entry by \(1.17\times 10^{-2}\).
Dropping \(1/\Omega\) moves the normalized Schur complement by \(0.108\).
Reversing the non-Hermitian link on the terminal block moves it by \(1.36\).
With \(B_{01}\) itself, assembly, nested Schur, and the \(1/\Omega\) recurrence
agree with the direct resolvent at \(10^{-14}\) or smaller. Depth 3 is these
three regions only.

nf=512, initial constant \(Q\):

\[
H_{\mathrm{derived}}=\begin{pmatrix}1& i/5\\ -i/5& 2\end{pmatrix}
\]

to \(1.4\times 10^{-12}\), and

\[
B_{01}=\begin{pmatrix}
2.34731355775+1.47597252128\,i&
0.043082867673+0.023020565378\,i\\
0.043082867673-0.023020565378\,i&
4.69462711550-2.95194504256\,i
\end{pmatrix}.
\]

\(\Omega^{-1/2}B_{01}\) is stored as `derived.normalized_B`. The optional
rational example has the same \(H\) and \(\Omega=3/2\), and a different link.
Its Frobenius gap from \(B_{01}\) is \(5.405\). Nothing was adjusted toward it.
Its own exact controls still pass.

## Translation, Löwdin, support

For coefficients proportional to \(\Omega^x\) and constant \(Q\),
\(HT_1=\Omega T_1 H\). A 96-node Gauss quadrature of the closed
anticommutator has relative residual \(1.8\times 10^{-16}\). An independent
FFT product rule, which does not insert \(\tfrac12\log\Omega\) by hand, has
relative residual \(4.8\times 10^{-15}\) and agrees with the closed form at
\(6.6\times 10^{-15}\).

The owner even and odd lobes on a shared unit interval have overlap
\(-1.4\times 10^{-22}\). Next-nearest supports are disjoint, so the measured
\(3\times 3\) Gram is \(I\) and the Löwdin factor does not move the columns
(\(0\) at 512 points). Envelope mass outside \([n,n+2]\) is \(0\), and packet
mass on the bridge \((4,8)\) is \(0\). The circle momentum therefore does not
open a resolved corner: continuum corner \(0\), nf=512 corner Frobenius
\(1.2\times 10^{-14}\). The nf=256 bridge-image fraction is \(1.4\times 10^{-9}\);
nf=512 is \(8.6\times 10^{-15}\).

`direct_hamiltonian` and `apply_dirac` agree at \(1.8\times 10^{-14}\) (nf=256)
and \(2.6\times 10^{-14}\) (nf=512).

## Leakage, corners, refinement

These are different quantities.

| Quantity | Value |
|---|---|
| Projection leakage, six nf=512 columns | 0.931488, 0.931627, 0.903251, 0.903353, 0.965624, 0.965807 |
| Smallest captured fraction | about 0.259 |
| Corner Frobenius, nf=512 | \(1.2\times 10^{-14}\) |
| Link Frobenius, nf=512 | 6.200536 |
| \(\|J_{256}-J_{512}\|_{\max}\) | \(1.339\times 10^{-8}\) |
| Continuum quadrature versus nf=512 | \(7.9\times 10^{-14}\) |
| Quadrature panels 64 versus 96 | \(2.6\times 10^{-14}\) |

Leakage above \(0.1\) sets `closed_six_mode` false. The projected window is not
an invariant subspace of the Dirac operator.

## Saved \(T=0.05\) frame

The episode payload `nsc-spherical-feedback-episode-v1` is read, not rerun.
Its renewal flag is false. Final \(Q\) lies in \([0.219522,0.261283]\), so
\(q_0/Q\) lies in \([0.961886,1.144873]\). \(L\) and \(β\) remain the supplied
gauge. The common bulk density is still
`nsc_spherical_feedback_action.first_order_density`; the sample uses \(r_x=0.2\)
so the \(L/Q\) term is visible, and that value changes between \(q_0\) and one
saved node.

With the initial packets on that final \(Q\), the kinetic block's scale defect
is \(0.012666\) on the diagonal and \(0.040594\) on the link. The mass and shift
blocks, which still see \(L\) and \(β\) proportional to \(\Omega^x\), stay at
\(10^{-14}\). The saved final columns give a full-window diagonal defect
\(0.013034\) and link defect \(0.041877\). The corner stays \(4.0\times 10^{-10}\).
Nonconstant \(Q\) removes the instantaneous \(\Omega^n\) law for \(L/Q\) without
replacing the action functional.

## Assumptions

The compression is one angular block at the calibration \(\kappa\), in the
orthonormal frame of the prolonged rank-6 columns. Occupations do not enter
\(J\). Radius cancels between \(N=rL\) and \(q=rQ\). The checked shift
covariance is the scaling inside these three regions. Geometry grading is the
supplied initial gauge, not a derived fractal law.

## What remains unconnected

The theorem turns \(H_{\mathrm{derived}}\) and \(B_{01}\) into a finite window
and a Schur response. The Dirac image lies mostly outside the six columns, so
that window is not the evolution. The saved \(T=0.05\) frame keeps the same
action and the same \(L,β\) gauge, while evolved \(Q\) breaks the kinetic
grading. Renewal, a maintained six-mode subspace, and a local response that
stays inside the window are not obtained. No infinite nest, \(\Lambda\)CDM
comparison, or incoming-gate change is made or required.
