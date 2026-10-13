# Local-observer correspondence of bound finite identities

This is Stage 1 of the approved NSC arXiv path: one compact, reproducible
correspondence between a local observer in one nested region and the finite
operators already owned by the repository. It is a v1 append-only record. It
does not close the local incoming gate, publish to arXiv, or rewrite
historical bytes.

Reproduce with

```text
python3 scripts/derive_nsc_local_observer_correspondence.py --check
```

The owner is `src/recursive_horizons/nsc_local_observer_correspondence.py`.
The record is `results/development/nsc-local-observer-correspondence.json`.

## Layers of the statement

| Layer | What is being said |
|---|---|
| Imported algebra | Block Schur/Feshbach reduction, retarded elimination with an initial-state remainder, Klich's determinant identity, and Dirac charge conjugation are established mathematics. They are not NSC novelty. |
| Repository-derived identities | The joined resolvent equals the ordered Schur complement; the Laplace identity retains \(B(zI-H_{BB})^{-1}\psi_B(0)\); finite \(\dot E_A=\operatorname{Tr}(CJ_A)+\operatorname{Tr}(C\dot h_A)\); \(\delta\log\det\mathscr K=\operatorname{Tr}(G_c\delta K_c)+\operatorname{Tr}(S^{-1}\delta S)\); the \(1/\Omega\) chain matches the direct multi-room inverse; a local scalar sheet coupling admits an invariant Dirac sector. |
| NSC interpretation | A local observer in one nested region is described by that reduced block, its retarded memory, the retained initial-state/noise data, and the same-action variation of the remaining Hamiltonian. Inside/outside are relational roles of \(\mathcal T_\Theta\), not gauge-charge labels. |
| Nonclaims | No complete \(\Lambda\)CDM background or perturbation match. Sheet exchange is not charge conjugation and is not antimatter. Antimatter is not mechanical pressure. The retarded map does not select an occupation. The local incoming gate remains OPEN. |

## Bound equations

Exact block reduction of a finite resolvent, independently checked on a
noncommuting \(2+2\) split distinct from the published Dirac lattice:

\[
S=K_A-BK_B^{-1}B^\dagger,\qquad
G_{AA}=S^{-1}=\bigl((zI-H)^{-1}\bigr)_{AA},\qquad
\Sigma_A(z)=B(zI-H_{BB})^{-1}B^\dagger.
\]

Reversed multiplication, dropping \(B^\dagger\), replacing the child resolvent
by \(BB^\dagger/z\), and transposing the coupling are negative controls: each
disagrees with the joined inverse. This is imported linear algebra. The
published [boundary Schur record](https://doi.org/10.5281/zenodo.23090015) remains the
Dirac-operator instance; this correspondence does not rerun that generator.

Eliminating the complementary region in time retains memory *and* the
eliminated initial state ([energy transfer](https://doi.org/10.5281/zenodo.23090015)):

\[
i\dot\psi_A=H_{AA}\psi_A
-i\int_0^t B\,e^{-iH_{BB}(t-s)}B^\dagger\psi_A(s)\,ds
+B\,e^{-iH_{BB}t}\psi_B(0).
\]

The covariance of the last term contains \(B e^{-iH_{BB}t}C_{BB}(0)e^{iH_{BB}s}B^\dagger\)
and the initial \(C_{AB}\) block. The retarded self-energy is independent of
occupation; the noise term is not. The same Schur map has vanishing vacuum
current and a nonzero excited current. Omitting the initial child source
changes the resolvent solution.

The same Hamiltonian owns conservation and metric variation
([common-source derivation](nsc-common-source-derivation.md),
[declared action](nsc-declared-action-scope.md)):

\[
\dot E_A=\operatorname{Tr}(CJ_A)+\operatorname{Tr}(C\dot h_A),\qquad
\dot E_{\rm tot}=\operatorname{Tr}(C\dot H),\qquad
\operatorname{Tr}([H,C]H)=0,
\]

\[
\left.\frac{\delta\Gamma_F}{\delta X^A}\right|_{X_+=X_-}
=-\operatorname{Tr}\!\left(C\frac{\delta H}{\delta X^A}\right),
\qquad
\delta\log\det\mathscr K
=\operatorname{Tr}(G_c\delta K_c)+\operatorname{Tr}(S^{-1}\delta S).
\]

The child's determinant is counted once. Omitting \(\operatorname{Tr}(G_c\delta K_c)\)
fails. The existing two-mode CTP owner reproduces equal-history normalization
and the finite energy/work identity, and still lists `full continuum
coefficient match` as unresolved.

Normalized inheritance on a finite chain
([regulated recursion](nsc-regulated-recursion.md)):

\[
\Gamma_n(x)=xI-H-\frac1\Omega\,b\,\Gamma_{n+1}(x/\Omega)^{-1}b^\dagger,
\qquad H_n=\Omega^n H,\quad B_n=\Omega^n b.
\]

The reduced inverse agrees with the direct multi-room matrix. Dropping
\(1/\Omega\) or reversing a noncommuting link fails. The control does not
select a physical \(\Omega\).

The tested two-sheet embedding
([observable bridge](nsc-observable-bridge.md)):

\[
\Pi_-=\frac{I_8-\tau_3\otimes\gamma^5}{2},\qquad
W^\dagger H_8W=\boldsymbol\alpha\cdot\mathbf p+\Phi\beta,
\qquad
\psi^c=i\gamma^2\psi^*.
\]

Sheet exchange is \(\tau_1\otimes I_4\). Charge conjugation in the restricted
sector is \((\tau_1\otimes i\gamma_2)K\). Those operators are distinct on a
generic spinor. An identity sheet link is gapless at \(|p|=|\Phi|\). The
complementary rank-four sector is equally invariant; the actual throat has
not selected a sector.

## Principal nonclaims

- No complete \(\Lambda\)CDM background \(H(z)\) or perturbation
  (growth, lensing, BAO/CMB) match.
- Sheet exchange is not antimatter and is not charge conjugation.
- Antimatter is not mechanical pressure.
- Schur/Feshbach algebra is not claimed as an NSC theorem.
- The local incoming gate under \(C_\Sigma[g]=U_g C_{\rm up}U_g^\dagger\)
  remains OPEN.

## Missing common-action evidence

For the historical vacuum-matched common-action candidate, the following
stronger claims remain open:

1. the absolute covariantly renormalized continuum stress on the varying throat;
2. same-action metric backreaction of an evolving geometry (the vacuum-work
   pulse is a prescribed drive);
3. a full continuum CTP coefficient match;
4. a derived local observable with an appropriate background/perturbation comparison;
5. a closed local incoming gate on \(I=S(1)+[0.12,0.18]\), only for the
   source-fixed application that requires that particular incoming state.

Those gaps are recorded, not filled. The finite common-metric successors
have their own [status](current-result.md); they do not establish the
continuum vacuum-matched claims above. The incoming campaign remains paused.
Its closure and a cosmological fit are not prerequisites for the
[exterior atlas](nsc-cosmological-atlas.md). Each candidate names
its actual dependencies; this historical list does not start a campaign.

## Transmitting memory and the finite Gaussian differential

For the whole-line PG convention $G^R(t,s)=-i\theta(t-s)U(t,s)$, the
transmitting packet resolvent enters as $G_J^R(z)=-\mathcal R(z)$ and
$D_J^R(z)=-\mathcal K(z)$. This energy-dependent inverse carries eliminated
bulk memory and is not a branch unitary or an instantaneous Hermitian link.
Its harmonic vertex obeys
$\delta\mathcal R_A(z_o,z_i)=-J^\dagger R(z_o)W_AR(z_i)J$ and
$\delta D_{J,A}^R=\mathcal K(z_o)\delta\mathcal R_A\mathcal K(z_i)$.
The 33-channel control used $z_i=0.2+0.25i$, $z_o=0.6+0.25i$;
the real transfer $0.4$ is a response argument, not a selected history.
Retarded and advanced blocks preserve their distinct causal directions.
The Keldysh kernel additionally requires the full Cauchy covariance and
packet/complement cross correlations; retarded response alone determines
no occupation.

For full canonical branch unitaries, the retained first differential is

$$
\delta\Gamma_G=-i\operatorname{Tr}
\left[Q^{-1}C_0(\delta V_-^\dagger V_++V_-^\dagger\delta V_+)\right],
\qquad Q=I-C_0+C_0V_-^\dagger V_+.
$$

The 33-channel finite-difference checks authenticate this imported identity.
Equal-history common changes are flat and select no kernel. With
$V_\pm(s)=e^{\pm isX/2}V$, the difference derivative is
$\operatorname{Tr}(VC_0V^\dagger X)$; $X=I/\sqrt{d_c}$ gives
$\operatorname{Tr}(C_0)/\sqrt{d_c}$ independently of $V$.
Admissible metric/embedding derivatives, the same-action boundary remainder
and reference/local allocation are still required for physical selection.
The finite Gaussian complementary determinant is not silently identified
with the separately unresolved $\Gamma_{\rm rest}$.

## Historical vacuum-matched prescription

That candidate defined
$\mathcal B_\Lambda^E=\Gamma_{\rm heat,\Lambda}^E
-\Gamma_{\rm canonical,vac}^E$ in a common normalization, then proposed

$$
Z_{\rm vm}=e^{i(\mathcal B_\Lambda^L[g_+]-\mathcal B_\Lambda^L[g_-])}
\det(I-C_0+C_0U_-^\dagger U_+).
$$

The real, state-independent phase preserves equal-history normalization and
canonical overlap modulus. This is an explicit quantization choice, without
a proof of causal metric dynamics or a unique Euclidean-to-state prescription.
Its bounded proper-time quadratic radius vertex has
$\mathcal B_{2,L}(\omega)=\mathcal B_2(-\omega^2)\in\mathbb R$ at finite
spatial rank. For the 12-point AP circle, $R=2,a=N=q=1$, one reduced
$\kappa=1$ block and cutoff 4, the canonical response `-0.37527623091`
plus conversion `0.33146152518` gave `-0.04381470572`, differing from the
heat response by `-5.35e-15` at imaginary external frequency `0.7`.
The multiplicity four was not applied. Continuum subtraction, other metric
variations, full curved history and physical source closure remain open;
this historical choice does not replace the later declared canonical action.
