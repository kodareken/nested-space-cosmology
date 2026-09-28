# NSC road to arXiv: prove the nonlocal-to-local correspondence and one self-consistent realization

## Summary and claim contract

The central insight is that NSC does not need different local physics. It needs a mathematically explicit commuting reduction:

\[
\text{full nested system}
\;\xrightarrow{\text{hide unresolved regions}}\;
\text{ordinary local equations with an effective source}.
\]

For

\[
\mathbb D=
\begin{pmatrix}
D_{\rm local}&B\\
B^\dagger&D_{\rm outside}
\end{pmatrix},
\]

the exact local resolvent is

\[
P(z-\mathbb D)^{-1}P
=
\left[
z-D_{\rm local}
-B(z-D_{\rm outside})^{-1}B^\dagger
\right]^{-1}.
\]

The outside region disappears from the local description but remains as the effective self-energy

\[
\Sigma_{\rm outside}(z)
=
B(z-D_{\rm outside})^{-1}B^\dagger.
\]

The same effective action must produce its stress tensor,

\[
T^{\rm outside}_{\mu\nu}
=
-\frac{2}{\sqrt{-g}}
\frac{\delta\Gamma_{\rm outside}}{\delta g^{\mu\nu}},
\]

its retarded state response and its conservation law. A local observer may then describe that source using the same Einstein equations and effective components used by \(\Lambda\)CDM. NSC proposes an origin for the effective source; it does not replace the local equations.

The arXiv paper will maintain this claim ladder:

| Level | Claim |
|---|---|
| Exact mathematics | Schur elimination, state-memory terms and action variation turn unresolved regions into a local effective response. |
| NSC postulate | \(\mathcal T_\Theta^*\mathbb D_\Theta=\mathbb D_\Theta\) says the same operator law is inherited across nested regions. |
| New demonstrated result | One declared spherical gradient transition either satisfies both source-fixed incoming constraints or is excluded by a necessary relation. |
| Interpretation | The local observer’s “unknown” effective source can be interpreted as unresolved nested-region response. |
| Future work | Actual \(\Lambda\)CDM background and perturbation matching, global recurrence and a geometric derivation of charge conjugation. |

“Same mathematics” means the same reduced local equation form on the domain actually derived. It does not yet mean that NSC has reproduced every \(\Lambda\)CDM observable.

Antimatter will remain a charge-conjugate field excitation. The repository’s two-sheet construction may motivate a future connection between geometric exchange and charge conjugation, but the first paper will not call antimatter a mechanical pressure or claim that sheet exchange alone produces antiparticles.

## Execution roadmap

### 1. Consolidate the exact correspondence before more numerical work

Use the existing [claim ledger](../lab/docs/claim-ledger.md) and [handover](handover.md). Do not create another roadmap.

Create one compact, reproducible correspondence record that binds the existing:

- block-resolvent/Schur identity;
- retarded memory, initial-state and noise terms;
- common-action metric variation and conservation identity;
- inherited recursion relation;
- tested two-sheet Dirac embedding.

The record must distinguish imported mathematics from NSC’s identification. Its principal nonclaim is that the zero- and finite-momentum outside responses have not yet been shown to reproduce the complete dark-energy and cold-dark-matter background and perturbation laws.

The paper’s conceptual language will be field-based and operational:

- “particle” means a stable pole or localized field excitation;
- the vacuum is a field state with correlations, rather than literal nothingness;
- detector events are raw correlations that every interpretation must reproduce;
- the sphere is the first controlled symmetry class, not a proof that all structures are spherical.

### 2. Make the six certification methods feasible

Do not calculate final history-specific bounds for the present \(10^{-3}\) candidate. Validate reusable bound-producing methods first, then instantiate them on the eventual solved history.

Run three bounded workstreams in parallel, with Codex integrating and Grok Fleet handling independent implementation and review:

1. **Prepared source and low/subgap**

   - Enclose the original modes at \(\rho=1\).
   - Transport that enclosure through the unchanged homogeneous continuation to the recorded \(\rho_{\rm up}\).
   - Bound covariance and source arithmetic.
   - Replace or enclose the remaining low/subgap windows using the existing source machinery.
   - Keep the exact zero response of the three \(\ell=0\) groups separate from their nonzero baseline uncertainty.

2. **Continuous field/space/time**

   - Replace the unusably large global Fourier-tail estimate with a local physical-space or localized Fourier enclosure on the backward cone and cutoff ramps.
   - Use the corrected profile-domain helper and form the history-minus-reference residual before taking norms.
   - Validate one dominant family and one control family before scaling to all 60 positive families.
   - Propagate the resulting \(A,D\) enclosures to \(F,F_z\) and both constraint contractions.

3. **Changed-history UV remainder**

   - Continue the existing coefficient recurrence through the first noncancelling paired order.
   - Preserve the proven cancellation of the leading \(E^{-2}\) vacuum coefficient.
   - Derive numerical \(C_4\) and \(C_M\) constants for the actual operator defect.
   - Bound vacuum, thermal and coherent remainders separately beyond the 160/320 cutoffs.
   - Integrate the directed tail bound into the N and beta constraints.

In parallel, implement reusable methods for:

- the full between-node residual enclosure;
- remaining geometry, field and assembly arithmetic;
- componentwise aggregation under the nine-component gate schema.

The present known N-error sum is about \(8.999\times10^{-12}\), leaving approximately \(1.1001\times10^{-11}\) of the \(2\times10^{-11}\) error reserve for all six missing components. Replace the provisional component allocations with a recorded redistribution only after the pilots produce measured bounds. Preserve the total error reserve.

Proceed to the nonlinear search only when every bound producer returns a finite result and the projected componentwise total can plausibly fit the remaining budget. A loose or excessive bound means the proof method must be improved; it is not evidence against the physical class.

### 3. Build one corrected full-retarded solve path

The historical Newton and LM scripts use different grids, solvers and 47-node rulers. Preserve them as historical evidence. Add one production driver for the corrected path.

Extend the existing evaluation binding so that source identity, profile, target nodes, grid, interpolation, integrator and step-control mode are all immutable cache inputs.

The solve will use:

- starting history `0b0e4ced…`;
- 32 coefficients for each of \(w\) and \(U\);
- 129 nodes on \(I=S(1)+[0.12,0.18]\);
- the corrected DOP853, grid-1024, degree-48 settings;
- the full retarded tangent of the state;
- a rectangular \(258\times64\) Jacobian;
- the existing physical \(H^3(w)\oplus H^1(U)\) metric and geometric principal preconditioner.

At every Jacobian refresh:

1. Evolve the unchanged upstream state and all 64 retarded directions.
2. Assemble both geometric and matter derivatives.
3. Verify selected columns against centred finite differences of the actual constraint path.
4. Reject any mismatch of source, history or numerical binding.

Each proposed step must:

- stay inside the physical trust region;
- preserve positive radius and preparation support;
- be accepted only after a fresh 129-node value-only evolution reduces the declared joint merit;
- record predicted versus measured reduction.

Use Broyden updates between full Jacobians. Refresh after five accepted updates or two poor prediction ratios. Increase from 32 to 64 coefficients only if a controlled coefficient-tail or held-out-node comparison identifies representation error. Do not open 128 or 256 coefficients without a new recorded justification.

The search target is a fresh nodal merit at or below \(10^{-11}\). A stalled search remains OPEN. It becomes NON-EXISTENCE only if a necessary relation gives a strict lower separation over the entire declared \((w,U)\) class.

### 4. Produce the final local certificate

For a plausible final history:

1. Re-evolve all 60 positive and 120 signed source contributions from the same upstream preparation.
2. Evaluate at 257 nodes using certificate numerics.
3. Instantiate all history-dependent field, UV, between-node and arithmetic bounds.
4. Combine them with upstream, low/subgap, covered-region, interpolation and phase bounds.
5. Enclose the residual continuously over \(I\).
6. Execute the real mutation controls: dropped state derivative, coherence, weight duplication, \(k=-E\), and missing signed-energy sector.
7. Verify every JSON and NPZ descendant by hash from a pinned lab commit.

Introduce `NSC-LOCAL-INCOMING-GATE-CERTIFICATE-v2` with:

- declared history class and interval;
- state law and unchanged-source digest;
- witness history or obstruction;
- nodal and continuous residual bounds;
- all nine error components for N and beta;
- numerical settings and source coverage;
- mutation results;
- complete dependency closure;
- verdict `EXISTENCE`, `NON_EXISTENCE`, or `OPEN`.

EXISTENCE requires, for both components,

\[
\sup_{z\in I}
\left(
|\widehat{\mathcal E}_B(z)|+\epsilon_B(z)
\right)
\le 3\times10^{-11}.
\]

NON-EXISTENCE requires an error-controlled necessary relation excluding the whole declared class. No `None` bound, optimizer failure or finite-dimensional stall can close the certificate.

### 5. Write, release and submit the result

Prepare a focused LaTeX article titled:

> **Unresolved Regions as Local Sources: A Certified Incoming-Gate Test of Nested-Space Inheritance**

Keep the existing long Markdown manuscript as the public theory notebook. The arXiv article will be the focused scientific result:

1. The local-observer problem and inheritance postulate.
2. Exact nonlocal-to-local reduction.
3. The declared spherical geometry, state and source.
4. The two incoming constraints and numerical method.
5. The closed result and complete error budget.
6. Interpretation, limitations and future tests.

The abstract will state the local result first. It may say that NSC supplies a candidate origin for effective local sources. It will not claim that dark matter, dark energy, antimatter or global recurrence have been derived unless the displayed equations establish that claim.

The discussion may present the bubble/field picture as motivation. The experimentally defined antimatter sector stays distinct. The combined sheet/charge-conjugation map becomes a future conjecture with explicit tests: reversed current, equal mass spectrum, correct magnetic response, spectroscopy and annihilation correlations.

For public version `0.27.0`:

- extend the development-snapshot/import path to accept content-addressed NPZ payloads from the pinned final lab commit;
- preserve frozen historical manifests and bytes;
- bind the paper build to the certificate and its entire dependency graph;
- synchronize manuscript, metadata, citation and release versions;
- replay the locked certificate rather than rerunning old scientific campaigns;
- build the LaTeX PDF twice under pinned TeX tooling and require deterministic output;
- verify references, embedded fonts, searchable equations and text, page layout and PDF integrity;
- publish GitHub `main`, tag, release and PDF;
- verify remote commit, peeled tag, downloaded PDF hash, lab SHA and public SHA.

Use `gr-qc` as the primary arXiv category and `math-ph` as the planned cross-list. Check the account and endorsement status early; if endorsement is required, prepare a targeted request using the finished abstract and draft rather than mass-contacting potential endorsers. arXiv requires author self-submission and may require endorsement for a first submission or new category. [Submission rules](https://info.arxiv.org/help/submit/index.html), [endorsement](https://info.arxiv.org/help/endorsement.html)

Prepare an arXiv-compatible source archive with only the TeX, bibliography and used figures. Metadata must use ASCII-compatible input, include an abstract of no more than 1,920 characters and report page/figure counts. [Metadata rules](https://info.arxiv.org/help/prep.html)

Douglas Ek remains the sole accountable author. Significant AI assistance is disclosed in the paper; AI tools are not authors. [arXiv AI policy](https://info.arxiv.org/help/moderation/index.html#policy-for-authors-use-of-generative-ai-language-tools)

Assume the arXiv perpetual non-exclusive license unless Douglas selects another license before submission. The choice is irrevocable for that version. [License options](https://info.arxiv.org/help/license/index.html)

Douglas performs the final visual/content review and the final “Submit Article” action. After announcement, record the arXiv identifier in the repository and GitHub release.

## Verification and completion criteria

The route is complete when all of the following agree:

- The exact correspondence record shows how unresolved regions become a local effective source without adding a second copy of the action.
- A v2 certificate closes the declared local gate through EXISTENCE or scoped NON-EXISTENCE.
- The manuscript’s theorem, numerical result, scope and nonclaims match the certificate.
- The public repository contains the complete authenticated dependency chain.
- The release PDF and arXiv PDF are the same scientific article.
- The arXiv metadata, source archive and AI disclosure pass a dry run.
- Douglas has reviewed and submitted the article, and the eventual identifier is recorded.

Global parent/child matching, metric evolution, a full \(\Lambda\)CDM likelihood fit and a derivation equating geometric exchange with charge conjugation remain later research. They are not hidden prerequisites for this submission.

## Failure forecast

1. **The proof bounds remain far larger than the numerical movement.**
   This will appear as a field Fourier remainder, low/subgap source error or UV \(C_M\) term exceeding the remaining \(1.10\times10^{-11}\) N budget. Stop production scaling and improve that specific enclosure.

2. **The corrected retarded solve cannot approach the residual reserve.**
   This will appear as repeated poor trust-region ratios, a residual plateau or representation comparisons showing unresolved modes. Record OPEN; derive a class-wide obstruction only if the mathematics supports it.

3. **The interpretation outruns the proved correspondence.**
   This will appear when the abstract calls the outside response “dark matter,” “dark energy” or “antimatter” without the required stress, perturbation or charge-conjugation mapping. The claim ledger and independent review must block publication until the language matches the equations.


## Active execution cursor — 2026-09-28

Douglas explicitly reactivated this full plan on Mac. The narrower nested-quality
plan at commit da6fb7b9 is historical, and its completed finite theorem is reusable
background, not closure of this plan. The only current roadmap is this PLAN.md.
The current handover owns measured progress and next executable work.

Stage 1 correspondence exists and must be replayed with its declared dependencies.
Stage 2 is active: authenticate and reuse later source/field/UV work before adding
missing bound-producing methods. Stages 3–5 remain incomplete. No production solve
until all required bound producers are finite and the projected budget is viable.
Mac remains the integration owner; targeted retrieval of Windows results is allowed
for reuse without moving the project or restarting duplicate campaigns.
