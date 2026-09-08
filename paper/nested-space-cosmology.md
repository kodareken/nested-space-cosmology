# Nested-Space Cosmology

## One Recursive Spectral Gradient–Boundary Equation Across Scales

**Douglas Ek**<br>
Nested-Space Cosmology<br>
Version 0.1.0 — 8 September 2026 (unit-closure revision)

## Abstract

Nested-Space Cosmology proposes that particles, waves, local constants, unresolved (“dark”) metric response, and the black-to-child transition are projections of one recursive spectral operator rather than separately fitted sectors. The configuration of a room can change. The law, the dimensionless coefficient vector \(\Theta\), the energy meaning, and the boundary kernel must not. The one equation is the heat-kernel spectral action of a self-similar block Dirac operator together with the requirement that a parent/child dilation leave that operator invariant.

The construction uses established physics as a launch surface: spectral action and scale anomaly, Skyrme/BPS nuclear structure, Randall–Sundrum outside geometry, two-wall shadow matter, Wetterich flow, and the exact regular black-universe four-geometry. Binding those ingredients is not novelty. What this repository derives is a concrete operator, a two-sheet representation of a noninjective areal-radius chart, Gauss–Bonnet identities of that carrier, a scalar recursive outside kernel, a global horizon-penetrating spectral foliation, a proposed throat Neumann-jump operator whose energy-dependent maps remain unevaluated, and a unitary dilation that makes the recursion map executable.

A ZETA1 numerical chain then produced determinant-only inheritance-scale roots. Those historical numbers used the inconsistent mass conversion \((\lambda_j+\mu^2)/\zeta\). The declared additive-gap operator requires \(q_j=\lambda_j/\zeta+\mu^2\). Compensating the pure cutoff anomaly, and then correcting the units, removes those roots: in the finite, fixed-geometry family the subtracted derivative is strictly negative for every \(\zeta>3\pi/2\). Adding an independently weighted geometric action is forbidden double-counting. The missing owner is the child’s own next-child response. With explicit parent normalization the derived, still unsolved, energy-resolved equation is \(\Gamma_p(x)=K_p(x)-(1/\Omega)\,b\,\Gamma_c(x/\Omega)^{-1}b^{\dagger}\). \(\Omega\) is not selected. Dilation alone does not derive \(\zeta=\Omega^2\). Physical \(\zeta\), a particle spectrum, a dark-sector fit, and identification of our universe with a black-hole interior remain open.

> **Status.** This is a working preprint of a conjectural research programme, not a completed theory and not a peer-reviewed result. It does not prove Nested-Space Cosmology, a final inheritance scale \(\zeta\), a particle spectrum, a dark-sector fit, or that our universe is the interior of a black hole. Numerical values below are copied from committed compact JSON under [`results/`](../results/manifest.json). If a JSON `nonclaims` flag is `false`, the corresponding statement is not claimed.

## Authorship and AI contribution disclosure

**Douglas Ek** is the sole accountable author. He supplied the research direction, the organizing thesis that black holes, universes, singularities, infinity, and local constants belong to one inherited recursive architecture, the scientific decisions about what to import, derive, invalidate, or leave open, and publication responsibility.

The public documents, derivations, scripts, and compact records were developed in iterative collaboration with **ChatGPT** and **OpenAI Codex**. The AI collaborator drafted formal language, executable code, literature comparisons, and adversarial checks; it proposed operator identifications, truncation tests, and claim discipline; and it assisted with numerical investigation of the inheritance-scale chain. Every AI contribution is candidate work until it matches explicit mathematics, a committed compact JSON record, cited primary sources, and the JSON nonclaims. AI credit records material assistance with formulation, computation, drafting, comparison, and audit. It does not transfer scientific accountability and does not imply independent authorship.

Evidence labels follow [`THEORY.md`](../THEORY.md): **Postulate**, **Imported result**, **Repository derivation**, **Numerical diagnostic**, **Open prediction**, **Interpretive hypothesis**. The live checkpoint is [`docs/current-result.md`](../docs/current-result.md). Imported launch surfaces and the one joint claim are in [`docs/prior-art-and-open-claim.md`](../docs/prior-art-and-open-claim.md).

## 1. The one-equation claim

**Postulate.** Spaces form inside spaces through collapse, localization, and renewed expansion. Eternity is then an unending process of finite rooms, not one infinite container. Singularity and infinity are two directions of the same local-chart wall. Local constants are formation parameters of one room. Unresolved rooms act locally through the same kernel that makes a boundary black. Parent collapse and child expansion are two sides of one finite transition.

```text
Parent space
└── Local collapse / black hole
    └── Child space
        └── ...
```

**Postulate.** Particles, waves, forces, space, nuclear binding, collapse, and expansion are states, derivatives, or projections of one relational gradient–boundary system. Let \(\mathbb{D}_\Theta\) be one self-similar block Dirac operator whose diagonal blocks are rooms and whose off-diagonal blocks \(\Phi_n\) are finite parent/child maps. The organizing invariant partition equation is

$$
\boxed{
Z_{\Lambda}^{\mathrm{one}}(\mathbb D_\Theta)
=
\int \mathcal D\varphi\,
Z_\Lambda\!\left(
e^{-\varphi/2}\mathbb D_\Theta e^{-\varphi/2}
\right),
\qquad
\mathcal T_\Theta^*\mathbb D_\Theta=\mathbb D_\Theta.
}
$$

The present concrete regulated action inside that partition function is

$$
S_{\mathrm{one}}[\mathbb{D}_\Theta,\Psi]
=
\operatorname{Tr}\,e^{-\mathbb{D}_\Theta^2/\Lambda^2}
+\langle J\Psi,\mathbb{D}_\Theta\Psi\rangle,
\qquad
\mathcal{T}_\Theta^*\mathbb{D}_\Theta
=
\mathbb{D}_\Theta.
$$

The exponential profile is not an extra coefficient. **Repository derivation.** Composition of successive unresolved-resolution steps, with positivity and \(\mathbb{D}^2\) as generator, selects the heat semigroup as the present heat-kernel diagnostic. Remaining \(\Lambda\) names the unit of resolution. Identifying that positivity with the complete physical anomaly or covariance requires a common regulator, field measure, and derived variation.

The operational test is one frozen object \((S_{\mathrm{one}},\Theta,\mathcal{T})\) whose same fields, operator coefficients, boundary kernel, and scale map jointly produce:

1. finite stationary electron, proton, and neutron charge sectors;
2. quantum amplitudes and detector couplings of those same configurations;
3. nuclear masses and reaction \(Q\)-values without inserted nucleon masses;
4. the nuclear-fitted self-gravitating branch;
5. a stable trapped-to-defocusing-to-expanding transition; and
6. a non-tunable dimensionless prediction identifying our cosmology with the child side of that transition.

No individual item is project novelty. The invariant closure is. Until that joint solution exists, Nested-Space Cosmology remains a constructive working hypothesis: one equation, imported launch surfaces, a sequence of honest derivations and invalidated truncations, and an unsolved recursive child tail.

## 2. Imported foundations

These are launch surfaces, not discoveries. Primary sources and the compute-admission rule are collected in [`docs/prior-art-and-open-claim.md`](../docs/prior-art-and-open-claim.md). Binding an imported result to \(S_{\mathrm{one}}\) is bookkeeping unless a previously missing arrow is closed with the same \(\Theta\), boundary kernel, and transition map.

| Subject | Imported result | Project use |
|---|---|---|
| Electron interference | Individual detections build an interference pattern (Tonomura). | Constraint on any particle ontology. |
| Quantum state | QFT treats particles as field excitations. PBR constrains purely epistemic wavefunctions ([arXiv:1111.3328](https://arxiv.org/abs/1111.3328)). | Distinguish localized energy from the amplitude of the same field state. |
| Bell nonlocality | Loophole-free tests exclude local hidden variables ([PRL 115, 250401](https://doi.org/10.1103/PhysRevLett.115.250401), [250402](https://doi.org/10.1103/PhysRevLett.115.250402)). | A gradient medium cannot be an ordinary local three-dimensional fluid. |
| Solitons / Skyrme | Localized finite-energy field configurations; collective quantization can give fermionic baryons ([arXiv:1309.0820](https://arxiv.org/abs/1309.0820), [arXiv:1312.2960](https://arxiv.org/abs/1312.2960)). | Pattern for stationary charge sectors. |
| Nuclear energy | Masses and \(Q\)-values are differences of complete bound-state energies. | Required low-energy projection. |
| Einstein–BPS stars | Nuclear-fixed parameters produce compact-star branches ([arXiv:1503.03095](https://arxiv.org/abs/1503.03095)). | Gravity benchmark, not a discovery target. |
| Black universes | Exact regular solutions contain a trapped interior, a finite radius minimum, and an expanding cosmological interior ([arXiv:gr-qc/0611022](https://arxiv.org/abs/gr-qc/0611022)). | Limiting four-geometry. Section 6. |
| RS2 outside geometry | Projected five-dimensional Weyl term can carry the exotic null response while the local scalar remains positive-kinetic ([arXiv:0910.4930](https://arxiv.org/abs/0910.4930)). | Structural realization of \(\Gamma_{\mathrm{outside}}\). |
| Two-wall shadow matter | Other-wall matter gravitates locally and deflects light \(25\%\) more weakly at equal Newtonian mass ([arXiv:hep-th/9911055](https://arxiv.org/abs/hep-th/9911055)). | Outside particles acting on local fabric, plus a fixed discriminator. |
| Spectral action | One Dirac operator yields the Standard Model coupled to Einstein plus Weyl gravity ([arXiv:hep-th/9606001](https://arxiv.org/abs/hep-th/9606001)). | Local-room scaffold. Its finite Dirac/Yukawa spectrum remains input until the recursive fixed point determines it. |
| Local Higgs–dilaton anomaly | Common local Weyl mode of the regulated fermionic determinant ([arXiv:1210.2663](https://arxiv.org/abs/1210.2663)). | Import the common local scale; derive only the relative parent/child sheet. |
| Doubled spectral geometry | Two metrics and one off-diagonal field yield a canonical relative-metric interaction and a leading FLRW action ([arXiv:2012.06401](https://arxiv.org/abs/2012.06401), [arXiv:2201.03839](https://arxiv.org/abs/2201.03839)). | Leading two-geometry interaction. |
| Spectral propagation wall | Scalar, gauge, and graviton propagation shuts off near the spectral cutoff ([arXiv:1312.2235](https://arxiv.org/abs/1312.2235)). | High-resolution boundary of local boson propagation. |
| Wetterich flow | Exact effective-action flow with changing resolution ([DOI 10.1016/0370-2693(93)90726-X](https://doi.org/10.1016/0370-2693(93)90726-X)). | Resolution mathematics; the new task is the recursive boundary condition. |
| Baby-universe literature | Black-hole offspring and inherited/mutated constants, including cosmological natural selection ([arXiv:gr-qc/9404011](https://arxiv.org/abs/gr-qc/9404011)). | Prior hypothesis class, not an observation. |

**Imported result.** Observed interference proves wave behavior and localized detections. It does not select a unique ontology in which a wavefunction is a second material substance. Exact black-universe solutions show that a black-hole exterior and an expanding interior can belong to one regular four-geometry. They do not show that our observed universe is such an interior.

A scientific calculation belongs in this repository only when its compact record contains a coefficient shared by two previously separate scales, a stationary configuration of the frozen action, a prediction whose datum was not used to determine \(\Theta\), a stable parent/child solution sourced by that same action, or a dimensionless discriminator against ordinary cosmology. Standalone repetitions of interference, a known soliton, a published nuclear table, or a black-universe metric are not project results.

## 3. The concrete recursive operator

**Repository derivation.** All preceding projections are calculations of one operator:

$$
\mathbb{D}_\Theta
=
\begin{pmatrix}
\ddots & \ddots & & \\
\ddots & D_{n-1} & \Phi_{n-1} & \\
& \Phi_{n-1}^{\dagger} & D_n & \Phi_n \\
& & \Phi_n^{\dagger} & D_{n+1} & \ddots \\
& & & \ddots & \ddots
\end{pmatrix}.
$$

\(D_n\) contains the geometry, causal cone, internal algebra, and locally resolved fields of room \(n\). \(\Phi_n\) is the finite gradient/boundary map. Eigenvalues and poles are masses and propagation scales; eigenvectors and residues are wave amplitudes of those same states; heat-kernel invariants give local gravitational, gauge, vacuum, and boundary coefficients; modes concentrated outside one diagonal block but overlapping it through \(\Phi\) act locally as dark response; a nonlinear spectral geometry in which a diagonal causal branch ends while the off-diagonal map remains finite is a black-to-child transition.

**Imported result.** Wetterich’s exact flow supplies the resolution equation. Nested-Space completion is the additional recursive boundary condition \(\Gamma_{n+1,\Omega k}[\mathcal{T}_\Theta\Phi]=\Gamma_{n,k}[\Phi]\), not a new proof of the flow itself.

**Open prediction.** The joint flow and recursion should leave only one relevant direction: an overall dimensional unit. Compact record [`nsc-1-s-one-relevant-direction-count.json`](../results/nsc-1-s-one-relevant-direction-count.json) closes the \(L0\)–\(L6\) Skyrme subspace to that rank-three, nullity-one condition. It does not close the remaining gauge, spin, mixing, and boundary coefficients, and it does not predict mass ratios.

**Repository derivation.** Nested \(L0\)–\(L6\) self-equality reduces the four Skyrme scale coefficients to one overall normalization ([`nsc-1-s-one-nested-pair-closure.json`](../results/nsc-1-s-one-nested-pair-closure.json)). The constant dictionary then *names* local constants as Hessian and spectral outputs of one stationary \(S_{\mathrm{one}}\): characteristic-cone speed, metric stiffness as \(G^{-1}\), localized-sector masses, representation/topological spin and charge, vacuum-plus-boundary energy as \(\Lambda_n\), Planck length as a derived resolution scale, and absolute zero as the lower spectral bound ([`nsc-1-s-one-constant-dictionary.json`](../results/nsc-1-s-one-constant-dictionary.json)). Those names are not numerical values of nature’s constants. The compact record states that the electron/proton/neutron spectrum is not solved.

**Imported result.** The bosonic spectral action arises, with specified coefficient modifications, as the term required to cancel the scale anomaly of the spectrally regularized fermionic action ([arXiv:1001.2036](https://arxiv.org/abs/1001.2036)). The quantum form is therefore one invariant partition function, not independently weighted fermionic and bosonic actions ([arXiv:1106.3263](https://arxiv.org/abs/1106.3263)). Geometry is the compensating anomaly of the same spectrum. Adding five-dimensional Einstein–Gauss–Bonnet or heat geometry again with a free weight double-counts that geometry.

**Repository derivation.** Splitting local resolution into a mean field and a relative field shows that the off-diagonal \(\Phi\) depends only on the mean resolution, while the relative field changes the two diagonal propagators. At equal sheet resolution, the same \(\Phi\) is the spectral mass gap and, after eliminating the child sheet, the visible outside self-energy. With the exponential profile and the imported doubled Dirac operator, that same \(\Phi\) is also the relative-metric interaction and simultaneously shifts the effective vacuum curvature ([`nsc-1-s-one-local-two-sheet-anomaly.json`](../results/nsc-1-s-one-local-two-sheet-anomaly.json)). Setting \(\Phi=0\) removes the gap and the outside self-energy at once. This identity does not predict an electron mass. Raw finite-determinant invariance of the relative two-sheet symbol does not set the regulated relative Hessian to zero ([`nsc-2-zeta1-unit-closure-check.json`](../results/nsc-2-zeta1-unit-closure-check.json)).

## 4. Particle/wave and local-room restrictions

**Imported result.** Individual electron detections build an interference pattern. Quantum field theory describes particles as field excitations. Loophole-free Bell tests exclude local hidden variables. PBR constrains purely epistemic wavefunctions under preparation independence. A gradient medium therefore cannot be an ordinary local three-dimensional fluid.

**Imported result.** Solitons and Skyrmions are localized finite-energy field configurations; collective quantization can give fermionic baryons. Nuclear masses and fusion/fission \(Q\)-values are differences of complete bound-state energies.

**Repository derivation**, classified as an imported collective-coordinate reduction rather than project novelty. On the BPS core, one charge-localized configuration and its collective wave amplitude share the same field state and the same rest energy. No separate point particle or pilot wave is inserted ([`nsc-1-s-one-particle-wave-identity.json`](../results/nsc-1-s-one-particle-wave-identity.json)). The compact record states that the electron solution, spin-half, electric charge, and Bell reduction are not obtained, and that complete \(\Theta\) is not calibrated.

**Postulate.** Detection is a local transfer of the same conserved charge and energy at a detector boundary, not conversion of a wave substance into a particle substance.

**Imported result.** The Chamseddine–Connes spectral action produces the Standard Model coupled to Einstein plus Weyl gravity from one Dirac operator. This repository imports that local-room reduction ([`nsc-1-s-one-spectral-room-bind.json`](../results/nsc-1-s-one-spectral-room-bind.json)). It does not reconstruct each known Standard-Model term as a discovery. The finite Dirac/Yukawa matrix still contains measured fermion masses and mixings as input. Nested-Space completion is stronger: the combined functional flow and recursive fixed point must determine that finite spectrum with the same \(\Theta\). Until it does, the spectral action is an imported unifying scaffold, not the derivation of particle masses from the room.

**Imported result.** Full Einstein–BPS calculations produce compact-star mass-radius branches from nuclear-fixed parameters. The compact record [`nsc-1-gravitating-bps-observation-link.json`](../results/nsc-1-gravitating-bps-observation-link.json) notes that the nuclear-fitted full-field maximum mass \(3.34\,M_\odot\) lands in the observed compact-object transition region. That is an imported gravity benchmark. The same record forbids reading \(3.34\,M_\odot\) as a precision observational fit, treating a static star as a collapse trajectory, or claiming that BPS stress already defocuses or that the parent/child energy ledger is closed.

**Open prediction.** The same Hessian of \(S_{\mathrm{one}}\) must later supply finite rest energy, propagation, interference, detector coupling, spin/statistics, charge, and no-signalling Bell correlations without a second parameter vector. That joint particle/wave result is not available before a physical \(\zeta\) and a solved \(\Theta\).

## 5. Outside response and two-sheet geometry

**Postulate.** Degrees of freedom outside the gradient that defines one local room are not locally visible, but they are not absent. “Dark” names unresolved response. “Black” names a causal/boundary continuation into another room. Neither name is a license for a second substance.

**Repository derivation.** Recursive self-equality replaces an arbitrary outside source by a self-consistent spectrum. At quadratic scalar order

$$
\Gamma
=
K-\frac{b^2}{\Gamma},
\qquad
\Gamma
=
\frac{K+\sqrt{K^2-4b^2}}{2},
$$

on the branch continuous with \(\Gamma\to K\) as \(b\to 0\) ([`nsc-1-s-one-recursive-outside-kernel.json`](../results/nsc-1-s-one-recursive-outside-kernel.json)). One link coefficient then generates the complete infinite tail. The scalar equation proves only that recursive integration can remove an arbitrary outside function. It does not predict particle masses, lensing, or a black-universe bulk.

**Imported result.** An explicit RS2 construction already produces black-universe solutions with a positive-kinetic local scalar, the exotic null response being carried by the projected five-dimensional Weyl term ([arXiv:0910.4930](https://arxiv.org/abs/0910.4930); [`nsc-1-s-one-outside-black-bind.json`](../results/nsc-1-s-one-outside-black-bind.json)). Garriga and Tanaka’s two-wall result makes matter on the other wall gravitate locally as shadow matter and predicts \(25\%\) weaker light deflection at equal Newtonian mass, i.e. the frozen ratio \(3/4\) ([arXiv:hep-th/9911055](https://arxiv.org/abs/hep-th/9911055); [`nsc-1-s-one-shadow-matter-bind.json`](../results/nsc-1-s-one-shadow-matter-bind.json)). That discriminator is imported, not discovered here.

**Numerical diagnostic.** The frozen constant two-wall tensor ratio predicts \(\gamma=1/2\) for pure shadow dominance. The compared SLACS sample is centered at \(\gamma=0.98\pm 0.07\) ([arXiv:astro-ph/0607657](https://arxiv.org/abs/astro-ph/0607657)). Pure shadow dominance on those kiloparsec scales does not pass ([`nsc-1-s-one-shadow-lensing-observation.json`](../results/nsc-1-s-one-shadow-lensing-observation.json)). That is a rejection of the constant low-energy two-wall projection as the dominant Newtonian mass, not a rejection of recursive outside geometry as a whole.

**Postulate.** A parent-to-child geometry with a positive minimum areal radius cannot be represented globally by one single-valued metric function of the areal radius. The same radius occurs on both sides of the minimum while the parent and child metric states differ.

**Repository derivation.** The complete parent/child geometry is single-valued only after adding a two-sheet branch state ([`nsc-1-s-one-two-sheet-black-geometry.json`](../results/nsc-1-s-one-two-sheet-black-geometry.json)). An explicit compact even five-dimensional carrier then holds both black-universe branches. The transition-fixed curvature \(b_2/F_0=-18/1015\) determines the five-dimensional Gauss–Bonnet coefficient

$$
\frac{\alpha_{\mathrm{GB}}}{L_\star^2}
=
\frac{1015}{144}
$$

([`nsc-1-s-one-gauss-bonnet-null-closure.json`](../results/nsc-1-s-one-gauss-bonnet-null-closure.json)). On the local sheet \(y=0\), the same coefficient supplies a negative effective null curvature while leaving the local matter null source nonnegative, a constant tangential vacuum-form component \(108/1015\), and a normal component equal to minus one half of the room’s own geometric action density ([`nsc-1-s-one-gauss-bonnet-component-closure.json`](../results/nsc-1-s-one-gauss-bonnet-component-closure.json)).

**Numerical diagnostic.** \(108/1015\) is not identified with the observed cosmological constant. \(L_\star\) is not derived. Linear stability is not proved.

**Repository derivation.** That coefficient closure sits at the five-dimensional Einstein–Gauss–Bonnet critical/unique-vacuum relation \(\alpha_{\mathrm{GB}}\Lambda_0=-3/4\) ([`nsc-1-s-one-gauss-bonnet-criticality.json`](../results/nsc-1-s-one-gauss-bonnet-criticality.json)). **Numerical diagnostic.** On the actual background, the pure critical bulk time-principal matrix has rank \(5/10\) on \(y=0\) ([`nsc-1-s-one-gauss-bonnet-kinetic-rank.json`](../results/nsc-1-s-one-gauss-bonnet-kinetic-rank.json)). The induced local spectral Einstein block restores rank \(10/10\) ([`nsc-1-s-one-induced-interface-rank.json`](../results/nsc-1-s-one-induced-interface-rank.json)). A scan of \(401\) positive relative normalizations still retains complex physical characteristic roots in tangent directions ([`nsc-1-s-one-interface-characteristic-scan.json`](../results/nsc-1-s-one-interface-characteristic-scan.json)). The finite local truncation \(R_5+\alpha_{\mathrm{GB}}\mathcal{G}_5+R_4|_{y=0}\) is therefore not the complete healthy action. The omitted off-diagonal \(\Phi\) modes must complete the characteristic system. Adding or rescaling another diagonal coefficient is not the repair. Nested-Space theory is not rejected by that scan.

**Open prediction.** A scale-dependent tensor projection of the same recursive kernel must recover \(\gamma_{\mathrm{eff}}\approx 1\) on solar and galactic dynamical/lensing scales while retaining nonlocal outside response and the nonlinear black-to-child saddle. That prediction is not available from the constant two-wall truncation.

## 6. Exact black-universe benchmark, labeled prior art

**Imported result.** The exact regular black-universe solution of Bronnikov, Dehnen, and Melnikov already exhibits a trapped interior, a finite positive areal-radius minimum, and an expanding cosmological interior in one four-geometry ([arXiv:gr-qc/0611022](https://arxiv.org/abs/gr-qc/0611022); DOI [10.1007/s10714-007-0430-6](https://doi.org/10.1007/s10714-007-0430-6)). This repository reproduces that geometry as a limiting trapped-to-expanding benchmark. It does not present the reproduction as Nested-Space novelty.

Compact record [`nsc-1-exact-black-universe-defocusing.json`](../results/nsc-1-exact-black-universe-defocusing.json) reproduces, on the imported metric:

- a trapped interval with positive complete-null-Raychaudhuri \(Q\) and a later expanding anti-trapped interval on the same affine congruence;
- a finite areal-radius minimum equal to \(1\) in the source normalization;
- exact interior limits \(V_{-\infty}=9\pi\) and \(B_{-\infty}=-3\pi\);
- field-equation residuals consistent with the published equations at working precision.

The compact `nonclaims` are all `false` for: our observed universe being this specific solution; microscopic stability of the phantom effective description; coupling of the nuclear Skyrme sector to this solution.

A curvature singularity or an infinite limit alone also does not encode a child domain. The transition map, its conserved charges, and the recursive tail must be supplied by \(S_{\mathrm{one}}\). Section 8 records the exact child vacuum \(\Lambda_{e,\mathrm{child}}L_\star^2=18\pi\) implied by that same imported asymptotic; that identity constrains \(\zeta\), it does not identify our cosmology with the solution.

## 7. Global foliation and self-adjoint throat

**Repository derivation.** The original constant-\(t\) slices cannot define the inheritance spectrum because \(t\) and \(\rho\) exchange causal roles across the horizon. The exact metric nevertheless admits a global horizon-penetrating time

$$
d\tau=dt+\frac{\sqrt{1-A}}{A}\,d\rho
$$

and a positive spatial metric on which a self-adjoint Dirac Hamiltonian can live ([`nsc-2-zeta1-foliation.json`](../results/nsc-2-zeta1-foliation.json)). The unwarped \((\tau,\rho)\) block has determinant exactly \(-1\). The identity \(1-A=3(x-\sin x\cos x)/\sin^2 x\) with \(x=\pi/2-\arctan\rho\) shows \(A<1\) on the complete real-\(\rho\) domain, so the Painlevé–Gullstrand form remains global. Compact nonclaims: the spatial foliation alone does not fix \(\zeta\), and the relative determinant is not computed.

**Repository derivation.** Cutting at the throat produces a proposed Neumann jump operator rather than an inserted function ([`nsc-2-zeta1-self-adjoint-domain.json`](../results/nsc-2-zeta1-self-adjoint-domain.json)):

$$
\Phi_{\mathrm{throat}}(E)=\mathcal{N}_p(E)+\mathcal{N}_c(E).
$$

The compact record establishes radial partner identities and a two-component flux-cancellation example. The maps \(N_{\mathrm{parent}}(E)\), \(N_{\mathrm{child}}(E)\), the APS projectors, lapse/shift additions, and \(\Phi_{\mathrm{throat}}=R(E)\) remain unevaluated declarations. Flux cancellation and a proposed gluing identity do not establish their physical identification. The cited BFK theorem assumes compactness, a product collar, and spectral conditions; its direct applicability to this noncompact warped geometry is unproved ([arXiv:math/0304347](https://arxiv.org/abs/math/0304347)). The compact record also states that the formal domain is not a computed spectrum, that the lowest partial wave is not the complete Dirac operator, and that the physical gap is not fixed.

## 8. ZETA1 numerical investigation

The inheritance scale is

$$
\zeta=(\Lambda L_\star)^2,
\qquad
\mu^2=1-\frac{3\pi}{2\zeta}.
$$

**Repository derivation.** Matching the local on-room vacuum-form component \(108/1015\) to the doubled gap, then using the exact expanding child’s Ricci sign, selected a *conditional* one-scale gap \(|\Phi|^2/\Lambda^2=1006/1015\) and a boundary cross-projection \(\Xi_{\mathrm{boundary}}=54/503\) ([`nsc-1-s-one-child-orientation.json`](../results/nsc-1-s-one-child-orientation.json)). **Repository derivation.** The exact child asymptotic rejects that magnitude. The actual child vacuum is \(\Lambda_{e,\mathrm{child}}L_\star^2=18\pi\), so a real positive gap requires \(\zeta>3\pi/2\). Direct \(\Lambda L_\star=1\) gives a negative \(\mu^2\) and is not the physical child-matched branch ([`nsc-1-s-one-child-scale-correction.json`](../results/nsc-1-s-one-child-scale-correction.json)). The identities for \(\Phi\) as gap, Schur self-energy, and metric link remain. The numbers \(1006/1015\) and \(54/503\) are **superseded** as child-matched predictions.

The ZETA1 chain then searched for a child-allowed stationary scale of the joined-minus-disconnected Dirac spectrum on the global slice. Every root in the following table is a **numerical diagnostic**. None is a physical \(\zeta\). The later anomaly compensation (Section 9) supersedes all of them as stationary points of the one equation.

| Compact record | What was computed | Diagnostic value | Why it is not physical \(\zeta\) |
|---|---|---|---|
| [`nsc-2-zeta1-lowest-mode.json`](../results/nsc-2-zeta1-lowest-mode.json) | Isolated lowest-mode heat minimum | \(\zeta=5.096657218448003\), \(\mu=0.27458356021684416\) | Lowest sector only; heat trace, not the determinant; **superseded candidate**. |
| [`nsc-2-zeta1-angular-tower.json`](../results/nsc-2-zeta1-angular-tower.json) | Angular tower through \(n=12\) | Diagnostic \(n=1\) heat root \(5.096622100098284\); \(n=2\) removes that heat-only root; \(n\le 12\) minimum derivative \(7.350807275212887\) | Heat-only truncation; **superseded candidate**; absence of a heat root is not a theorem against stationarity. |
| [`nsc-2-zeta1-regulated-determinant.json`](../results/nsc-2-zeta1-regulated-determinant.json) | Exponentially regulated joined-minus-disconnected determinant | \(\zeta_{\det}=4.99072409354754\), \(\mu_{\det}=0.2361577587050819\) | Uncompensated determinant; warp, lapse, shift, and anomaly omitted; **superseded candidate**. |
| [`nsc-2-zeta1-y-boundary-sensitivity.json`](../results/nsc-2-zeta1-y-boundary-sensitivity.json) | Factorized compact-\(y\) controls | Periodic root \(4.8765555133228515\); Neumann root \(4.79786026064662\); Dirichlet and MIT bag: no root | Factorized unwarped control; a zero mode was not selected in order to force a root; **superseded candidate**. |
| [`nsc-2-zeta1-orbifold-parity.json`](../results/nsc-2-zeta1-orbifold-parity.json) | Joint fermion/metric \(S^1/\mathbb{Z}_2\) parity | \(\zeta_{\mathrm{orb}}=4.747863009733466\), \(\mu_{\mathrm{orb}}=0.08643829078229333\) | Unwarped factorization; **superseded candidate**. |
| [`nsc-2-zeta1-warped-y.json`](../results/nsc-2-zeta1-warped-y.json) | Exact conformal compact warp | \(\zeta=4.748389947082489\), \(\mu=0.08707307719677547\) | Factorized radial-plus-\(y\) sum, not the full Dirac operator with lapse and shift; **superseded candidate**. |

Every corresponding JSON record states that the candidate is not promoted, that lapse and shift or the local matrix anomaly remain omitted, and that the electron mass is not predicted.

## 9. Why determinant-only roots are superseded by anomaly compensation

**Repository derivation.** The invariant partition function requires cancellation of pure cutoff rescaling while retaining the physical change of the child-constrained gap. For one eigenvalue the physical child-link derivative after compensation is negative whenever the gap constraint \(\mu^2=1-3\pi/(2\zeta)\) is imposed ([`nsc-2-zeta1-anomaly-decomposition.json`](../results/nsc-2-zeta1-anomaly-decomposition.json)).

At the warped determinant candidate \(\zeta=4.748389947082489\):

| Piece | Value |
|---|---:|
| Determinant derivative | \(0.01201522080638507\) |
| Pure cutoff anomaly | \(12.336590431250533\) |
| Physical child-link derivative | \(-12.324575210444092\) |

The compensated child-link derivative stays negative on the \(257\)-point scan \(3\pi/2<\zeta\le 200\), from \(-19.422735759304963\) to \(-0.5875504421475721\).

Those historical numbers used the inconsistent conversion \((\lambda_j+\mu^2)/\zeta\). They remain numerical diagnostics of that proxy.

**Therefore those determinant-only \(\zeta\) values are invalidated as stationary points of the one equation.** \(\zeta_{\det}\) and \(\zeta_{\mathrm{warped}\,y}\) are uncompensated determinant diagnostics, not inheritance-scale predictions.

**Repository derivation.** The missing owner is the recursive child tail, not an independently weighted Einstein–Gauss–Bonnet or heat action ([`nsc-2-zeta1-anomaly-owner-correction.json`](../results/nsc-2-zeta1-anomaly-owner-correction.json)). Geometry is already the compensating anomaly of the same fermionic determinant. Adding that geometry again is forbidden double-counting. The calculations through the scan joined one parent to one child treated as a terminal asymptotic domain. Recursion requires that the child’s boundary response already contain its own child tail. The compact `nonclaims` include `recursive_tail_is_known_to_restore_a_scale_root: false`.

Related routes already rejected on their own JSON nonclaims, and not repaired by retuning \(\Phi\) or adding a second geometric weight:

- \(\Lambda L_\star=1\) as the physical child-matched gap;
- the lowest-mode heat-only scale minimum;
- independently weighted Einstein–Gauss–Bonnet or heat geometry as the missing scale owner;
- pure critical five-dimensional Einstein–Gauss–Bonnet as the complete interface law;
- the standalone bosonic heat-trace graviton covariance as a positive Stieltjes propagator (\(g''(0)=228/175>0\); [`nsc-1-s-one-reflection-positivity-obstruction.json`](../results/nsc-1-s-one-reflection-positivity-obstruction.json));
- the inverse relative heat Hessian as Osterwalder–Schrader positive metric covariance (the full \(31\)-time reflection matrix has minimum eigenvalue \(-1.3679466754324912\); the source record counts fourteen negative modes while portable Linux counts thirteen because one near-zero mode crosses the numerical counting tolerance; [`nsc-1-s-one-relative-reflection-test.json`](../results/nsc-1-s-one-relative-reflection-test.json));
- constant two-wall shadow matter as the dominant kiloparsec gravitating mass.

The fermionic sheet-resolution observable is Osterwalder–Schrader positive in the controlled flat sector ([`nsc-1-s-one-fermionic-relative-observable.json`](../results/nsc-1-s-one-fermionic-relative-observable.json)). That does not make \(\sigma_3\) the complete graviton, and it does not restore a scale root. Full-exponential positivity in that sector is a heat-kernel diagnostic, not the complete physical anomaly.

## 10. Unit-consistent mass term and finite-family sign theorem

**Repository derivation.** The ZETA1 definitions are \(\lambda_j=L_\star^2\operatorname{eig}_j(D_{\mathrm{spatial}}^2)\), \(\zeta=\Lambda^2 L_\star^2\), and \(\mu^2=|\Phi|^2/\Lambda^2=1-a/\zeta\) with \(a=3\pi/2\). The same declared additive-gap operator therefore requires

$$
\boxed{
q_j=\frac{\lambda_j}{\zeta}+\mu^2=1+\frac{\lambda_j-a}{\zeta}.
}
$$

The historical runners implemented \((\lambda_j+\mu^2)/\zeta\), a physical mass squared smaller by \(1/\zeta\) than the declared mass squared. Compact record [`nsc-2-zeta1-unit-closure-check.json`](../results/nsc-2-zeta1-unit-closure-check.json) reconstructs the old derivative as a control, then changes only that conversion. It retains radius 30, 750 radial half-intervals, angular sectors 1–12, 32 compact intervals, and the historical factorization and angular weights.

At the old warped candidate \(\zeta=4.748389947082489\):

| Quantity | Historical argument | Argument with consistent units |
|---|---:|---:|
| Determinant logarithmic derivative | \(0.01201522080638507\) | \(-24.319952223294422\) |
| Derivative after inherited cutoff subtraction | \(-12.324575210444092\) | \(-36.58292832809548\) |

The independent old-argument reconstruction differs from the stored derivative by \(1.53\times 10^{-8}\). With radial spacings \(0.08\), \(0.04\), and \(0.02\), the corrected subtracted derivative is \(-36.58813190\), \(-36.58292833\), and \(-36.58162543\). Its sign is resolved.

**Repository derivation.** For each mode the determinant piece is one half of \(E_1(q_j)\). At fixed spatial geometry, subtracting the same pure-cutoff term used historically gives

$$
\boxed{
\frac{d\Gamma_{\mathrm{sub},j}}{d\log\zeta}=-\frac{e^{-q_j}}{2q_j}.
}
$$

The disconnected radial matrix is the principal submatrix obtained by deleting the throat node. Hermitian eigenvalue interlacing and the positive decreasing weight \(f(\lambda)=e^{-q(\lambda)}/q(\lambda)\) then imply that the total subtracted derivative is strictly negative for every \(\zeta>3\pi/2\) in this finite, frozen-geometry, additive-gap family. The \(257\)-point scan corroborates the sign; it is not the basis for extending the statement between sample points. A diagnostic root of the unsubtracted finite determinant near \(\zeta=6.09675392014\) still has subtracted derivative about \(-14.15469\), so it does not solve the adopted subtracted scale equation. No physical scale is inferred.

**Consequence.** Refining this family or selecting another determinant-only minimum cannot close its scale equation. A further attempt needs a derived change in the physical operator, its geometry/link dependence on scale, or the anomaly prescription. Merely naming the omitted recursive tail does not prove that it supplies that change.

**Repository derivation.** For the same relative-sheet symbol

$$
D_\delta
=
\begin{pmatrix}
p e^{-\delta} & \Phi \\
\Phi & -p e^{\delta}
\end{pmatrix},
$$

the unregulated finite determinant is independent of \(\delta\). The proper-time functional used later in ZETA1 has a strictly positive regulated relative Hessian at \(\delta=0\) for \(p\neq 0\). Raw determinant invariance therefore cannot justify setting the regulated relative Hessian to zero. The primary anomaly paper distinguishes its normalization scale from the spectral cutoff ([arXiv:1106.3263](https://arxiv.org/abs/1106.3263)).

## 11. The current recursive-tail equation

**The current frontier is the unsolved energy-resolved recursive child tail.**

**Repository derivation.** Norm preservation and first-order Dirac scaling fix a minimal parent/child dilation ([`nsc-2-zeta1-recursion-map.json`](../results/nsc-2-zeta1-recursion-map.json)):

$$
(U_\Omega\psi)(x)=\Omega^{d/2}\psi(\Omega x),
\qquad
\mathcal{T}_\Omega=\sigma_1 U_\Omega.
$$

The factor \(\Omega^{d/2}\) cancels the coordinate Jacobian. The maps compose as a semigroup. First-order spectral quantities scale together, so \(E/\Lambda\) and \(\Phi/\Lambda\) are inherited on a transformed mode. That dilation does not by itself derive \(\zeta=\Omega^2\). The inheritance scale \(\zeta=(\Lambda L_\star)^2\) and the cutoff ratio \(\Omega=\Lambda_{\mathrm{child}}/\Lambda_{\mathrm{parent}}\) remain distinct until the geometry establishes their relationship.

A parent dimensionless energy \(x\) at one common dimensional energy is seen by the finer child as \(x/\Omega\). With \(b=B_{\mathrm{dim}}/\Lambda_{\mathrm{parent}}\), the recursive outside equation must be solved as

$$
\boxed{
\Gamma_p(x)
=
K_p(x)
-
\frac{1}{\Omega}\,
b\,\Gamma_c(x/\Omega)^{-1}b^{\dagger}.
}
$$

The unweighted form \(\Gamma(x)=K(x)-B^{\dagger}\Gamma(x/\Omega)^{-1}B\) is recovered if \(b_{\mathrm{sym}}=B_{\mathrm{dim}}/\sqrt{\Lambda_{\mathrm{parent}}\Lambda_{\mathrm{child}}}\). Both conventions are legitimate; they cannot silently share the same numerical link ([`nsc-2-zeta1-unit-closure-check.json`](../results/nsc-2-zeta1-unit-closure-check.json)). At vanishing argument this recovers the earlier quadratic fixed point on the branch continuous with \(\Gamma\to K\) as \(B\to 0\).

This map removes an undefined symbol from the one equation. It does not compute the mode-resolved tail. The compact `nonclaims` are all `false`:

- `Omega_value_selected`
- `functional_tail_solution_computed`
- `physical_zeta_promoted`
- `scale_root_restored`

The gate records `mode_resolved_tail_solved: false` and `zeta_derived: false`. An implementation fixture with \(\Omega=2\) then \(\Omega=3\) composes to \(\Omega=6\) (and the fixture identification \(\zeta=\Omega^2=4\)); that fixture does not select nature’s \(\Omega\) and does not derive \(\zeta=\Omega^2\) from dilation.

The recorded next result is to solve the functional tail on the warped parent/child mode spectrum and recompute the anomaly-compensated scale derivative. Do not restore a root by adding a second geometric weight, retuning \(\Phi\) after seeing the target, or inserting a dark function.

## 12. Predictions available only after \(\zeta\)

Until the energy-resolved tail is solved and a physical \(\zeta\) is derived without using the later compared datum, the following remain **open predictions**. They are not present results.

- Dimensionless mass ratios of electron, proton, and neutron stationary sectors from one \(\Theta\).
- Detector couplings, spin/statistics, charge, and no-signalling Bell correlations of those same configurations.
- Nuclear masses and reaction \(Q\)-values without inserted nucleon masses, then the nuclear-fitted self-gravitating branch as a holdout.
- A scale-dependent tensor projection of the recursive kernel with \(\gamma_{\mathrm{eff}}\approx 1\) on tested dynamical/lensing scales and a retained nonlinear black-to-child saddle.
- A non-tunable dimensionless identification of our cosmology with the child side of the transition, including any reading of \(108/1015\), \(18\pi\), \(1006/1015\), or \(54/503\) as an observed cosmological constant, electron mass, or \(\Xi_{\mathrm{NSC}}\).

**Open prediction**, restated as a restriction. Solving \(\Gamma_p(x)=K_p(x)-(1/\Omega)\,b\,\Gamma_c(x/\Omega)^{-1}b^{\dagger}\) on the warped transition spectrum could restore a physical scale root. The compact record states that the recursive tail is **not** known to restore a root. Failure of one truncation, coordinate chart, or numerical method is not rejection of Nested-Space Cosmology; neither is survival of a diagnostic root promotion of the theory. Refining the finite additive-gap family cannot close its scale equation.

## 13. Reproducibility

This repository publishes the original \(58\) v0.1.0 compact JSON records, preserved byte-for-byte, plus the 59th unit-closure follow-up and the scripts that generated them. Regeneration checks those records. It does not create a final \(\zeta\), a particle spectrum, a dark-sector fit, or a proof that our universe is inside a black hole. Inspection and regeneration are described in [`docs/reproducing.md`](../docs/reproducing.md). The follow-up itself is reproduced by

```text
python3 scripts/check_nsc_scale_closure.py --check
```

which compares every exact and numeric field. Isolated public reproduction copies the authenticated v0.1.0 inputs as auxiliary provenance files; it does not add them as extra steps of the historical 58-record chain.

Pinned dependencies are `mpmath==1.3.0`, `numpy==2.5.1`, `scipy==1.17.1`, and `sympy==1.14.0`. From the repository root:

```text
make check
make test
make reproduce
make paper
make verify
```

The live frontier record is [`results/nsc-2-zeta1-recursion-map.json`](../results/nsc-2-zeta1-recursion-map.json) (`NSC-2-ZETA1-RECURSION-MAP`). The unit-closure follow-up is [`results/nsc-2-zeta1-unit-closure-check.json`](../results/nsc-2-zeta1-unit-closure-check.json). Committed JSON is the public scientific object. Typical runners refuse to overwrite an existing file. Portable reproduction requires exact generator hashes, graph, schemas, formulas, classifications, gates, and nonclaims, plus the declared headline numerical observables within their manifest tolerances. The follow-up compares every recorded field. Platform-dependent coordinates of an otherwise identical diagnostic argmin are not promoted observables. Same-environment exact mode compares every output byte. Passing either mode does not promote a numerical diagnostic or convert an imported black-universe, Skyrme, spectral-action, or shadow-matter result into Nested-Space novelty.

## 14. Limitations and open claims

The programme is one equation plus a closure test. The present limitations are structural, not cosmetic.

- \(\Omega\) is not selected. The energy-resolved tail is not solved on the warped transition spectrum. Physical \(\zeta\) is not derived. Dilation alone does not derive \(\zeta=\Omega^2\).
- Throat maps \(N_{\mathrm{parent}}(E)\), \(N_{\mathrm{child}}(E)\), and \(\Phi_{\mathrm{throat}}=R(E)\) are unevaluated declarations. Flux cancellation is not a computed boundary response.
- Raw finite-determinant invariance does not remove the regulated relative Hessian. Full-exponential positivity is a heat diagnostic, not the complete anomaly.
- Electron, proton, and neutron stationary sectors are not solved. Standard-Model mass ratios are not predicted. The finite Dirac/Yukawa spectrum remains imported input.
- The exact black-universe four-geometry is prior art. This repository does not identify our universe with that interior, prove phantom microscopic stability, or couple Skyrme matter to the transition.
- \(108/1015\) is an on-room vacuum-form component, not the observed cosmological constant. \(1006/1015\) and \(54/503\) are superseded one-scale diagnostics.
- Constant two-wall shadow matter does not dominate the tested kiloparsec gravitating mass. Recursive outside geometry as a whole is not thereby rejected, and it is not thereby confirmed.
- Pure critical Einstein–Gauss–Bonnet is not the complete interface law. The omitted object is the off-diagonal \(\Phi\) block.
- Standalone bosonic heat-trace covariance and the inverse relative heat Hessian fail necessary positivity tests. The fermionic sheet observable that passes is not the graviton.
- Lorentzian finite-momentum stability of the nested operator, linear stability of the two-sheet carrier, and a complete parent/child energy ledger remain open.

The only project-level open claim is the invariant closure of \((S_{\mathrm{one}},\Theta,\mathcal{T})\) listed in Section 1. Until that joint solution exists, Nested-Space Cosmology is a constructive working hypothesis with a derived unsolved tail, not a proof.

## References

Primary sources support only the established ingredient named in the citing sentence. They do not transfer authority to Nested-Space postulates.

1. A. Tonomura, J. Endo, T. Matsuda, T. Kawasaki, and H. Ezawa, “Demonstration of single-electron buildup of an interference pattern,” *Am. J. Phys.* **57**, 117 (1989). [DOI 10.1119/1.16104](https://doi.org/10.1119/1.16104).
2. M. F. Pusey, J. Barrett, and T. Rudolph, “On the reality of the quantum state,” *Nature Phys.* **8**, 475 (2012), [arXiv:1111.3328](https://arxiv.org/abs/1111.3328).
3. B. Hensen *et al.*, “Loophole-free Bell inequality violation using electron spins separated by 1.3 kilometres,” *Nature* **526**, 682 (2015); L. K. Shalm *et al.*, *Phys. Rev. Lett.* **115**, 250401 (2015), [DOI 10.1103/PhysRevLett.115.250401](https://doi.org/10.1103/PhysRevLett.115.250401); M. Giustina *et al.*, *Phys. Rev. Lett.* **115**, 250402 (2015), [DOI 10.1103/PhysRevLett.115.250402](https://doi.org/10.1103/PhysRevLett.115.250402).
4. C. Adam, C. Naya, J. Sanchez-Guillen, and A. Wereszczynski, “Bogomol’nyi-Prasad-Sommerfield Skyrme model and nuclear binding energies,” *Phys. Rev. Lett.* **111**, 232501 (2013), [arXiv:1309.0820](https://arxiv.org/abs/1309.0820); [arXiv:1312.2960](https://arxiv.org/abs/1312.2960).
5. C. Adam, C. Naya, J. Sanchez-Guillen, R. Vazquez, and A. Wereszczynski, “Neutron stars in the BPS Skyrme model: mean-field limit vs. full field theory,” *Phys. Rev. C* **92**, 025802 (2015), [arXiv:1503.03095](https://arxiv.org/abs/1503.03095), [DOI 10.1103/PhysRevC.92.025802](https://doi.org/10.1103/PhysRevC.92.025802).
6. K. A. Bronnikov, H. Dehnen, and V. N. Melnikov, “Regular black holes and black universes,” *Gen. Relativ. Gravit.* **39**, 973 (2007), [arXiv:gr-qc/0611022](https://arxiv.org/abs/gr-qc/0611022), [DOI 10.1007/s10714-007-0430-6](https://doi.org/10.1007/s10714-007-0430-6).
7. K. A. Bronnikov and E. V. Donskoy, “Possible black universes in a brane world,” *Grav. Cosmol.* **16**, 42 (2010), [arXiv:0910.4930](https://arxiv.org/abs/0910.4930), [DOI 10.1134/S0202289310010068](https://doi.org/10.1134/S0202289310010068).
8. J. Garriga and T. Tanaka, “Gravity in the Randall–Sundrum brane world,” *Phys. Rev. Lett.* **84**, 2778 (2000), [arXiv:hep-th/9911055](https://arxiv.org/abs/hep-th/9911055), [DOI 10.1103/PhysRevLett.84.2778](https://doi.org/10.1103/PhysRevLett.84.2778).
9. A. H. Chamseddine and A. Connes, “The spectral action principle,” *Commun. Math. Phys.* **186**, 731 (1997), [arXiv:hep-th/9606001](https://arxiv.org/abs/hep-th/9606001), [DOI 10.1007/s002200050126](https://doi.org/10.1007/s002200050126).
10. A. A. Andrianov and F. Lizzi, “Bosonic spectral action induced from anomaly cancelation,” *JHEP* **05**, 057 (2010), [arXiv:1001.2036](https://arxiv.org/abs/1001.2036), [DOI 10.1007/JHEP05(2010)057](https://doi.org/10.1007/JHEP05(2010)057).
11. A. A. Andrianov, M. A. Kurkov, and F. Lizzi, “Spectral action, Weyl anomaly and the Higgs-Dilaton potential,” *JHEP* **10**, 001 (2011), [arXiv:1106.3263](https://arxiv.org/abs/1106.3263), [DOI 10.1007/JHEP10(2011)001](https://doi.org/10.1007/JHEP10(2011)001).
12. M. A. Kurkov and F. Lizzi, “Higgs-Dilaton Lagrangian from spectral regularization,” *Mod. Phys. Lett. A* **27**, 1250203 (2012), [arXiv:1210.2663](https://arxiv.org/abs/1210.2663).
13. M. A. Kurkov, F. Lizzi, and D. Vassilevich, “High energy bosons do not propagate,” *Phys. Lett. B* **731**, 81 (2014), [arXiv:1312.2235](https://arxiv.org/abs/1312.2235), [DOI 10.1016/j.physletb.2014.02.053](https://doi.org/10.1016/j.physletb.2014.02.053).
14. A. Bochniak and A. Sitarz, “On stability of Friedmann–Lemaître–Robertson–Walker solutions in doubled geometries,” *Phys. Rev. D* **103**, 044041 (2021), [arXiv:2012.06401](https://arxiv.org/abs/2012.06401); “Spectral interaction between universes,” *JCAP* **04**, 055 (2022), [arXiv:2201.03839](https://arxiv.org/abs/2201.03839), [DOI 10.1088/1475-7516/2022/04/055](https://doi.org/10.1088/1475-7516/2022/04/055).
15. C. Wetterich, “Exact evolution equation for the effective potential,” *Phys. Lett. B* **301**, 90 (1993), [DOI 10.1016/0370-2693(93)90726-X](https://doi.org/10.1016/0370-2693(93)90726-X).
16. L. Smolin, “The fate of black hole singularities and the parameters of the standard models of particle physics and cosmology,” [arXiv:gr-qc/9404011](https://arxiv.org/abs/gr-qc/9404011).
17. A. S. Bolton, S. Rappaport, and S. Burles, “Constraint on the post-Newtonian parameter \(\gamma\) on galactic size scales,” *Phys. Rev. D* **74**, 061501 (2006), [arXiv:astro-ph/0607657](https://arxiv.org/abs/astro-ph/0607657).
18. Z.-Y. Fan, B. Chen, and H. Lü, “Criticality in Einstein-Gauss-Bonnet gravity: gravity without graviton,” *Eur. Phys. J. C* **76**, 512 (2016), [arXiv:1606.02728](https://arxiv.org/abs/1606.02728).
19. Y. Lee, “Burghelea–Friedlander–Kappeler’s gluing formula and the adiabatic limit of the spectral \(\zeta\)-determinant of a Dirac Laplacian,” [arXiv:math/0304347](https://arxiv.org/abs/math/0304347).
20. J. Nemec, D. Tománek, and G. Cuniberti, surface Green-function recursion, Appendix A.3 of [arXiv:0711.1088](https://arxiv.org/abs/0711.1088).

Repository records cited in the text live under [`results/`](../results/manifest.json). The labeled theory is [`THEORY.md`](../THEORY.md). The current compact checkpoint is [`docs/current-result.md`](../docs/current-result.md).
