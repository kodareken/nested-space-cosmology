# FGC-1-HYP1-CON3-CAU1: conditional gauge Cauchy preservation

**FGC-1-HYP1-CON3-CAU1** composes CON2-MPROP1’s exact metric-derived
subsidiary identity with COMP1’s exact activated compatibility witness. It
closes one boundary-free logical implication; it does not construct a
solution, an initial hypersurface, or an IBVP.

For the actual REF1 gauge vector

\[
C^\alpha[g]=-widetilde g^{\rho\sigma}
(\Gamma^\alpha{}_{\rho\sigma}-\bar\Gamma^\alpha{}_{\rho\sigma}),
\]

CON2 proves

\[
L_{\hat g}C^\nu
=\nabla_\mu\mathcal E_{\rm REF1}^{\mu\nu}
+E_\phi\nabla^\nu\phi+E_\chi\nabla^\nu\chi,
\]

with

\[
L_{\hat g}C^\nu
=\frac F2\widehat g^{\mu\beta}\nabla_\mu\nabla_\beta C^\nu
+\text{terms linear in }(C,\nabla C).
\]

The full REF1 metric equations and **both** unmodified scalar equations are
therefore required to obtain the homogeneous system

\[
L_{\hat g}C=0.
\]

If \(F>0\), \(\widehat g\) is smooth Lorentzian, \(\Sigma\) is spacelike for
\(\widehat g\), and

\[
C|_\Sigma=0,\qquad
\widehat n^\mu\nabla_\mu C|_\Sigma=0,
\]

standard uniqueness for a linear normally hyperbolic operator implies

\[
C=0
\]

through the boundary-free domain of dependence of \(\Sigma\). The REF1 gauge
extension then vanishes and its metric equation reduces to the unredefined
ACT1 equation. This is the same conditional gauge-preservation logic used in
the modified-harmonic Cauchy construction of
[Kovács and Reall](https://arxiv.org/html/2003.08398v2#S4), now tied to this
repository’s two-scalar ACT1 Noether sign and its metric-derived reference
gauge vector.

## Exact COMP1 witness

At the activated COMP1 local two-jet, the reproducer independently checks

\[
F>0,\qquad \widehat g^{tt}<0,
\]

and the Lorentzian sign of the complete hat inverse. It recomputes

\[
C^A=0,\quad \nabla_aC^A=0,\quad \mathcal H=0,\quad \mathcal M=0,
\qquad A\in\{t,r\}.
\]

The exact normal gauge-extension map is

\[
\begin{pmatrix}X_{tt}\\X_{tr}\end{pmatrix}
=
\begin{pmatrix}
-2473901162487/137438953472&0\\
0&2473901162487/137438953472
\end{pmatrix}
\begin{pmatrix}\nabla_tC^t\\\nabla_tC^r\end{pmatrix},
\]

whose determinant is nonzero. COMP1’s exact theorem composition therefore
shows that its conditional QIFT1 REF1 root supplies the zero normal gauge data
at this one local point.

## Exact claim boundary

The theorem is conditional: **if** a sufficiently smooth full REF1 plus
two-scalar solution exists with the compatible data on a hat-spacelike
hypersurface, zero gauge data remain zero until an unhandled boundary enters
the domain of dependence. The artifact does not construct that solution or
hypersurface. It also does not propagate the physical Hamiltonian/momentum
constraints, map incoming gauge constraints to the main boundary fields,
prove a quasilinear IBVP, authorize evolution, or derive collapse or a finite
transition surface.

Reproduce with:

```bash
python3 scripts/reproduce_fgc_hyp1_con3_cau1.py
```
