# Nested-Space Cosmology

**One recursive law connecting matter, geometry and the spaces beyond a boundary.**

Nested-Space Cosmology develops Douglas Ek's proposal that a black-hole-like
collapse in one space can continue as an expanding space within it. Parent and
child inherit a common law, expressed by one organizing equation:

$$
\mathcal T_\Theta^*\mathbb D_\Theta=\mathbb D_\Theta.
$$

The idea becomes concrete through a block Dirac operator. Its diagonal blocks
describe locally resolved fields; its off-diagonal blocks carry the interaction
between rooms. The same interaction appears in a particle's spectral gap, the
response of an unresolved room, and the energy exchanged across their boundary.

[Read the paper](paper/nested-space-cosmology.pdf) ·
[Explore the equations and evidence](docs/current-result.md) ·
[Follow the central idea](THEORY.md)

## The strongest connections

| Connection | Mathematical result | Evidence |
|---|---|---|
| One coupling, mass and external response | In the invariant scalar-link sector, $E^2=\lvert\mathbf p\rvert^2+\Phi^2$ and $\Sigma_p(E)=\Phi^2(E+\boldsymbol\alpha\cdot\mathbf p)^{-1}$ | [Dirac reduction](docs/nsc-observable-bridge.md) |
| Geometry determines boundary transport | Evaluated curved spinor maps agree with direct elimination and independent radial integration | [Boundary channels](docs/nsc-chiral-boundary.md) |
| One coherence, stored energy and transfer | $z_B=\mathrm{Tr}(BC_{cp})$ gives $E_{\rm link}=2\Re z_B$ and $\dot N_p=2\Im z_B$ | [Common-source derivation](docs/nsc-common-source-derivation.md) |
| Changing geometry creates Dirac excitations | The prescribed radius pulse's work equals the increase in field energy | [Vacuum-work experiment](docs/nsc-vacuum-work.md) |
| Parent state determines child vacuum stress | The computed canonical massless state has two negative radial null contractions at the smooth neck and outward power | [Parent-matched source](docs/nsc-unruh-state.md) |
| The same operator supplies geometric forces | Independent lapse, shift, radial-metric and radius variations give the stress projections and conservation identities | [Metric source equations](docs/nsc-adm-source-constraints.md) |

The scalar-link reduction is exact for its stated coupling and invariant sector.
The computed massless throat kernel has vector and axial channels. Connecting
that transport kernel to a physical mass interaction is one of the equations
that the common action must determine.

## Six physical questions

The programme connects six questions through the same operator and state:

1. **Constants:** how local propagation, gravitational coupling and resolution emerge.
2. **Black holes:** how collapse continues through a finite geometry into expansion.
3. **Dark response:** how adjacent spaces contribute to background gravity, clustering and lensing.
4. **Particles and waves:** how field excitations interact with detectors and produce records.
5. **Mass and antimatter:** how chiral coupling, gauge charges and charge-conjugate states arise.
6. **Recursive infinity:** how successive rooms inherit states, clocks and observable probabilities.

The [theory overview](THEORY.md) develops the proposal. The
[current result](docs/current-result.md) connects each step to its equations,
computed quantities and next dependency.

## The equation being completed

The physical feedback loop is

$$
\text{geometry}\ \longrightarrow\ \text{Dirac state and boundary response}
\ \longrightarrow\ \text{quantum stress}\ \longrightarrow\ \text{geometry}.
$$

At the recorded neck, the parent-matched massless source gives
$T_{kk}^{+}=-0.0528075900$, $T_{kk}^{-}=-0.0527953956$, and outward
power $P_0=0.0001422206795$, in the specified throat units. The complete
source equation also includes density, both pressures, flux and the geometric
contribution. The [narrow closure gates](docs/development-update-2026-09-10.md)
now bind canonical causal response to the charged MMP normalization. The
recursive zero-tadpole law
$\Gamma_{\rm rel}=(1-\mathcal P_0)\Gamma_{\rm one}$ fixes the homogeneous
unlinked $V_{\rm full}=0$ while preserving Einstein, gauge, curvature,
Casimir and link response.

## Read, inspect and reproduce

The repository contains the **100-record preprint collection** and an
**twelve-record development snapshot**, with original source hashes, generators
and comparison policies. The paper has received an editorial revision;
the [September 10 development update](docs/development-update-2026-09-10.md)
indexes the subsequent technical derivations.

    python3 -m pip install -e '.[paper]'
    make demonstrate

The default demonstration reads fourteen authenticated stored cases.
[Reproduction instructions](docs/reproducing.md) describe the explicitly
requested calculation routes. Publication changes use path, provenance and
PDF checks; completed scientific calculations are reused.

## Research status

NSC is an active theoretical programme and working preprint, not yet peer
reviewed. Its exact identities and numerical results apply to the operators,
states, domains and approximations specified in their records. A complete
self-sourced solution, identified particle and nuclear sectors, and independent
cosmological predictions remain research objectives. The current results do
not establish a theory of everything or an observational identification of
our universe with a black-hole interior. In particular, the quoted
$3.34M_\odot$ stellar endpoint belongs to an imported, nuclear-calibrated
Einstein–BPS benchmark.

Established Dirac theory, spectral geometry, quantum-field methods and regular
black-universe solutions supply the mathematical foundation. Original sources
are credited in the paper and [prior-art notes](docs/prior-art-and-open-claim.md).
The project contribution is the explicit integration of these relationships
through a common recursive operator.

## Authorship

**Douglas Ek** provides the conceptual synthesis, research direction and
scientific responsibility. ChatGPT/OpenAI Codex assisted with derivations,
software, computation and writing; Grok supplied bounded investigations.
The principal collaborating sessions are identified as GPT-5.6 Sol and
GPT-6 Astra. The equations, source records and reproducible artifacts carry
the scientific evidence.

Code: MIT. Original prose and figures: CC BY 4.0.
[Citation](CITATION.cff) · [Contributing](CONTRIBUTING.md) · [Licence](LICENSE)
