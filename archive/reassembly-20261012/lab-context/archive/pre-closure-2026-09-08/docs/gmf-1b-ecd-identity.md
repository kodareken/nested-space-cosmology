# GMF-1B-ECD-INT1: ECD contact-interaction identity gate

GMF-1B-ECD-INT1 is a bounded algebraic preflight for the Ventrella--Choptuik
massless left-chiral spinor doublet. It verifies an action/projection identity
for the frozen minimal Einstein--Cartan--Dirac contact term. It does not solve
an evolution, a null constraint, a bounce, a horizon crossing, or a global
parent--child spacetime.

There are two explicitly separated convention layers. The physical model is
first defined in a canonical `(+---)` ECD action convention, then the complete
matter equation is translated into the Ventrella--Choptuik `(-+++)` numerical
representation. No formula is allowed to mix the two layers:

\[
(+,-,-,-):\quad \Gamma^a,\ J_+^a,
\qquad
(-,+,+,+):\quad \gamma^a_{\rm VC},\ A_-^a,
\qquad \kappa=8\pi G.
\]

VC writes an imaginary ``(-,+,+,+)`` matrix representation. Its hermitizer is

\[
H=-i\gamma^0_{\rm VC},\qquad \bar\psi=\psi^\dagger H.
\]

The canonical layer uses the action and field-equation signs displayed by
Cabral--Lobo--Rubiera-Garcia; their paper writes `kappa^2=8 pi G`, while this
repository writes `kappa=8 pi G`. The executable bridge is

\[
\Gamma^a=-i\gamma^a_{\rm VC},\qquad
J^a=A^a=-i\,\bar\psi\gamma^a_{\rm VC}\gamma_\star\psi,
\qquad J_+^2=-A_-^2.
\]

The executable separately verifies both Clifford algebras, the common adjoint,
and the gamma-five relation

\[
\{\Gamma^a,\Gamma^b\}=2\eta_+^{ab},\qquad
\{\gamma^a_{\rm VC},\gamma^b_{\rm VC}\}=2\eta_-^{ab},
\]

\[
\bar\psi=\psi^\dagger\Gamma^0
=\psi^\dagger(-i\gamma^0_{\rm VC}),\qquad
\Gamma^5=i\Gamma^0\Gamma^1\Gamma^2\Gamma^3=\gamma_\star.
\]

The raw VC bilinear is therefore imaginary:
\(\bar\psi\gamma^a_{\rm VC}\gamma_\star\psi=iA^a\). This is a global
convention translation, not two simultaneously active metrics. The tetrad,
orientation/epsilon, index lowering, Einstein--Hilbert/curvature sign, and
stress-tensor sign must all be translated with the metric signature. The
Levi--Civita connection and physical geometry are unchanged by the constant
overall metric-sign convention, while the written scalar curvature and action
signs transform coherently. INT1 tests only the displayed matter-sector bridge;
it does not implement or solve the gravitational equations.

Cabral et al., [arXiv:1902.02222](https://arxiv.org/abs/1902.02222), combine
their effective-action decomposition (Eqs. 8--11 and 23), spin-vector relation
(Eq. 18), and substituted Dirac sector (Eq. 20) to give the full
torsion-eliminated interaction after combining the torsion-dependent Dirac and
gravitational-contortion pieces (rather than using the intermediate Dirac term
alone):

\[
L_{\mathrm{int},+}=-\frac{3\kappa}{16}J_aJ^a.
\]

Using \(J_+^2=-A_-^2\) gives the frozen VC-signature interaction

\[
L_{\mathrm{int},-}=+\frac{3\kappa}{16}A_aA^a.
\]

Their Eq. 83, \(i\Gamma^aD_a\psi=(3\kappa/8)J^a\Gamma_a\gamma_\star\psi\),
becomes exactly the VC left-hand equation

\[
\gamma^a_{\rm VC}D_a\psi
-i\frac{3\kappa}{8}A_a\gamma^a_{\rm VC}\gamma_\star\psi=0.
\]

INT1 computes both matrix expressions and requires their componentwise
residual to vanish to binary64 roundoff.

For the VC pair,

\[
\mathcal C=\frac{1}{2\sqrt\pi r\sqrt a},\qquad
A^2=-16\mathcal C^4|F|^2|G|^2,
\]

so sphere reduction of \(\lambda A^2\) gives

\[
L_{4,\mathrm{red}}=-\frac{3\kappa\alpha}{4\pi r^2a}|F|^2|G|^2.
\]

The Wirtinger sources are

\[
\frac{\partial L_{4,\mathrm{red}}}{\partial F^*}
=-\frac{3\kappa\alpha}{4\pi r^2a}|G|^2F,
\qquad
\frac{\partial L_{4,\mathrm{red}}}{\partial G^*}
=-\frac{3\kappa\alpha}{4\pi r^2a}|F|^2G.
\]

The direct gamma-matrix Hehl--Datta operator is reconstructed pointwise in the
declared local opposite-chirality two-vector basis at three nondegenerate test
angles. Its coefficients are

\[
N_F=-\frac{3\kappa}{2}\mathcal C^2|G|^2F,
\qquad
N_G=-\frac{3\kappa}{2}\mathcal C^2|F|^2G.
\]

The executable weak identity at those test points is

\[
\frac{\partial L_{4,\mathrm{red}}}{\partial F^*}=2\alpha N_F,
\qquad
\frac{\partial L_{4,\mathrm{red}}}{\partial G^*}=2\alpha N_G.
\]

This is deliberately **not** a time-evolution normalization claim. It only
compares the covariant interaction operator with the variation of the same
reduced interaction action and checks local basis reconstruction at the
declared angles. It does not integrate against a complete spin-weighted
spherical-harmonic basis or prove that every omitted `(ell,m)` coefficient
vanishes. A nonzero local reconstruction or coefficient/phase residual is a P0
stop for any later ECD initial-data construction; a genuine full harmonic
closure remains part of that later derivation.

Using the repository's `(-+++)` stress convention
\(T_{\mu\nu}=-2(\sqrt{-g})^{-1}\delta S_m/\delta g^{\mu\nu}\), the reduced
contact action has

\[
T^{(4)}_{\mu\nu}=\lambda A^2g_{\mu\nu}.
\]

Here \(A^2=A^I A_I\) is an internal Lorentz scalar and is tetrad independent at
fixed spinor components. KIN1 now verifies the equivalent coordinate-current
cancellation in all 16 coframe directions. Thus this is the complete
contact-sector tetrad stress, and
\(T^{(4)}_{\mu\nu}k^\mu k^\nu=0\) for a metric-null vector. It is still not the
full ECD effective stress because the torsion-free Dirac kinetic and mass
terms remain separate; it gives no null-defocusing conclusion.

For the exact regression \(\mathcal C=1\), \(\kappa=1\), \(F=3/4\),
and \(G=1/2\), the artifact verifies

\[
A^2=-\frac94,
\qquad
\lambda A^2=-\frac{27}{64}.
\]

One-amplitude cases require care: if \(F=0\) or \(G=0\), then \(A^2=0\)
and the cubic contact sources vanish, but the total axial current can still be
nonzero and null. This artifact distinguishes a nonzero torsion current from a
nonzero axial--axial contact invariant.

Reproduce the deterministic result with:

```bash
python3 scripts/reproduce_gmf1b_ecd_identity.py
```

Primary ingredients: [Ventrella and Choptuik (2003)](https://doi.org/10.1103/PhysRevD.68.044020),
[Hehl and Datta (1971)](https://doi.org/10.1063/1.1665738), and
[Kibble (1961)](https://doi.org/10.1063/1.1703702), with the complete
torsion-eliminated convention bridge taken from
[Cabral, Lobo, and Rubiera-Garcia (2019), especially Eqs. 8--11, 18, 20, 23, and 83](https://arxiv.org/abs/1902.02222)
and the reduced action/stress cross-check from
[Choudhury, Maity, and Lahiri (2024), Eqs. 7--10](https://link.springer.com/article/10.1140/epjc/s10052-024-13618-4).
