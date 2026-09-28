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
published [boundary Schur record](nsc-boundary-response.md) remains the
Dirac-operator instance; this correspondence does not rerun that generator.

Eliminating the complementary region in time retains memory *and* the
eliminated initial state ([energy transfer](nsc-energy-transfer.md)):

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

A stronger, cosmological result is blocked by the absence of:

1. the absolute covariantly renormalized continuum stress on the varying throat;
2. same-action metric backreaction of an evolving geometry (the vacuum-work
   pulse is a prescribed drive);
3. a full continuum CTP coefficient match;
4. a complete \(\Lambda\)CDM background and perturbation comparison;
5. a closed local incoming gate on \(I=S(1)+[0.12,0.18]\).

Those gaps are recorded, not filled. No scientific campaign is started here.
