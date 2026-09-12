# Nested-Space Cosmology

## One Inherited Frequency Spectrum Across Matter, Geometry, and Nested Space

**Douglas Ek**<br>
Version 0.16.0 · joint history/state BVP gate, 12 September 2026<br>
Working preprint

## Abstract

Nested-Space Cosmology asks whether one inherited frequency spectrum is enough
to describe what we separately call particles, waves, matter, antimatter,
dark gravity, black holes, and expanding space. A room is a locally complete
spacetime with its own clocks, rulers, and resolved spectrum. A black-hole
boundary is proposed as the handoff at which a compressed parent gradient
continues in a child room with a new absolute register and the same
dimensionless law.

The mathematics is organized by
$\mathcal T_\Theta^*\mathbb D_\Theta=\mathbb D_\Theta$. The eigenvalues of
$\mathbb D_n$ are the room's natural frequencies. Its diagonal blocks
describe locally resolved fields, while off-diagonal blocks transmit state and
energy between rooms. Eliminating an unresolved block produces a visible
Schur self-energy. Its near-zero-momentum projection is the proposed source of
smooth dark expansion; its finite-wavelength projection is the proposed source
of dark clustering and lensing.

Three exact connections anchor this picture. A scalar sheet link has an
invariant four-component sector with the ordinary Dirac mass shell, and the
same coupling fixes the visible self-energy. Curved full-spinor boundary maps
connect direct propagation to boundary elimination. Canonical metric variation
and state evolution connect geometric work, field excitation, and
cross-boundary energy transfer.

On the smooth expanding-child benchmark, affine-horizon data and the actual
complex exterior Dirac reflection determine a canonical massless source with
mean radial null stress approximately $-0.05280$, in units
$\hbar=c=L_{\mathrm{throat}}=1$ and normalization $\mu=1$. Both null
contractions are negative, and the state carries outward parent Killing power
$0.0001422206795$. A prescribed geometry pulse independently produces Dirac
pairs with a verified work balance. The calculation specifies density,
pressure, flux and their remaining componentwise metric residuals.

The canonical closed-time-path functional now owns causal response. A
recursive zero-tadpole law fixes homogeneous unlinked
$V_{\rm full}=0$ while preserving Einstein, gauge, curvature, Casimir, and
link response. Proper-volume and child-clock projection of the completed
tensor gives the conserved local density derivative $0.38776013531$ and
identifies the parent Killing power as child spatial momentum rather than a
second energy source. The initial backreaction gate then finds that this source
and the locked neck miss the Hamiltonian and momentum constraints already at
$T=0$. The permitted geometry-from-source route then uses the tensor's Landau
normal and density to select a replacement neck that closes both constraints.
The first coupled step then identifies the missing mode-resolved Gaussian
Cauchy state needed to recalculate pressure on that changed geometry. That
seed state is now serialized channel by channel and reconstructs the completed
unit-radius tensor. Conditional blockwise propagation is unitary for every
supplied history, while endpoint-identical histories give different maps.
The joint boundary-value audit therefore selects a measured nonuniqueness stop
and identifies the absent executable same-action history variation.
The paper follows the spectrum from intuitive picture to exact operator equations, using
100 preprint records and twenty-three post-preprint
development records. The remaining integration is the physical recursive
return state, joint geometry/scale solution, particle pole, and cosmological
projection of the same state.

## 1. What if one spectrum is enough?

A musical note is a stable pattern inside an oscillation. It is not a second
substance added to the wave. NSC applies that simple relation to fields and
geometry. For room $n$,

$$
\mathbb D_n u_{n,k}=\omega_{n,k}u_{n,k},
\qquad
\omega_{n,k}=\Lambda_n\widehat\omega_k(\Theta).
$$

The eigenvector is a possible pattern; the eigenvalue is its natural
frequency. Energy, frequency, and wavelength are related by

$$
E=\hbar\omega,
\qquad
\nu=\frac{\omega}{2\pi},
\qquad
\lambda=\frac{2\pi c}{\omega}.
$$

A wider frequency range resolves shorter lengths and supports more stable
combinations. Complexity is the growth of distinguishable patterns within the
spectrum, not the introduction of unrelated fundamental substances.

If $\Lambda_{n+1}=\Omega\Lambda_n$, parent and child can occupy different
absolute registers while preserving the same dimensionless intervals:

$$
\frac{\omega_{n+1,k}}{\Lambda_{n+1}}
=\frac{\omega_{n,k}}{\Lambda_n}
=\widehat\omega_k(\Theta).
$$

The central proposal is that spaces form inside spaces through collapse,
localization, and renewed expansion. A room changes configuration and scale
while inheriting the common law

$$
\mathcal T_\Theta^*\mathbb D_\Theta=\mathbb D_\Theta.
$$

The Perseus cluster provides a concrete low-frequency image of black-hole
feedback. Repeated outbursts from its central black hole drive pressure waves
through the hot gas at roughly $3.3\times10^{-15}$ Hz—one cycle per ten
million years and about 57 octaves below the B-flat above middle C. The waves
carry energy into their environment. [Chandra gives the measured scale and
octave calculation](https://chandra.harvard.edu/chronicle/0303/perseus/index.html).
NSC extends the same frequency language from cosmic pressure waves to field
resonances and geometric boundaries.

The diagonal blocks of $\mathbb D_\Theta$ describe room operators;
off-diagonal blocks describe their boundary interactions. The proposed
invariant partition and spectral prescription are

$$
Z_{\Lambda}^{\mathrm{one}}(\mathbb D_\Theta)=\int\mathcal D\varphi\,Z_\Lambda(e^{-\varphi/2}\mathbb D_\Theta e^{-\varphi/2}),
$$

$$
S_{\mathrm{one}}=\mathrm{Tr}e^{-\mathbb D_\Theta^2/\Lambda^2}+\langle J\Psi,\mathbb D_\Theta\Psi\rangle.
$$

A physical realization specifies the regulator, normalization, measure,
determinant phase, zero modes, state, and boundary domain together. Induced
geometric terms enter once. These ingredients form one chain:

```text
common action, quantum state and recursive domain
       -> boundary response and quantum stress
       -> coupled geometry and state evolution
       -> matter, energy transfer and physical modes
       -> expansion, clustering and lensing
```

The equations below move from this picture to the exact algebra and computed
evidence. General research status is stated once in the
[repository README](../README.md#research-status).

Douglas Ek supplies the conceptual synthesis and scientific responsibility.
ChatGPT/OpenAI Codex assisted with derivations, computation and writing;
Grok supplied bounded investigations. The principal collaborating sessions
are identified as GPT-5.6 Sol and GPT-6 Astra. Equations, source hashes and
reproducible artifacts provide the evidence for each result.

## 2. How one spectrum appears as six parts of reality

| Familiar name | Meaning in the common spectrum | Calculated connection |
|---|---|---|
| Constants and resolution | The local conversion between frequency, energy, length, and gravitational response | Characteristic currents, compact/light matching, and thermal conversion |
| Black hole and child space | The causal/resolution handoff between two spectral rooms | Smooth geometry, horizon transport, and parent-matched stress |
| Dark energy and dark matter | Zero- and finite-wavelength ranges of one unresolved boundary response | Schur self-energy, source projections, and energy-transfer identities |
| Wave and particle | Extended amplitude and stable pole of the same field | Pair creation, response poles, and Gaussian state functional |
| Matter and antimatter | Conjugate positive- and negative-frequency Dirac sectors | Invariant Dirac mass sector and charge-conjugation algebra |
| Recursive infinity | The same law and normalized spectrum continuing through descendants | Endpoint criteria, state transport, and normalized recursion |

<!-- nsc-claim:target-constants -->
### 2.1 Constants and physical resolution

The causal cone follows from the physical principal symbol. A dimensional
speed is expressed using the room's clock and ruler conventions. Gravitational
response supplies $G$, from which

$$
\ell_P=\sqrt{\hbar G/c^3}
$$

can be compared with the spectral resolution scale. Establishing an operational
minimum distance additionally requires localized probes and their backreaction.
The [constant dictionary](../results/nsc-1-s-one-constant-dictionary.json)
records the intended observables.

Temperature is defined by the physical many-body state. Its zero-temperature
limit must respect energy-origin invariance, ground-state degeneracy and zero
modes. Residual spin and vacuum correlations belong to that state's structure;
thermal excitation is determined by its occupations.

<!-- nsc-claim:target-continuation -->
### 2.2 Collapse and an expanding interior

Regular black-universe solutions already combine a trapped region, positive
minimum radius and expanding interior [6]. They provide the geometric
benchmark for the present Dirac calculation.

<!-- nsc-claim:imported-charged-throat -->
Maldacena, Milekhin and Popov supply a complementary source benchmark: a
semiclassical Einstein–Maxwell throat supported by charged massless Dirac
vacuum energy [32]. Their sections 2 and 5 establish the light magnetic
channels, Casimir stress and Einstein matching; sections 5.5–6 discuss longer
separations and stabilization. Applying this construction requires its magnetic
return geometry, state and coefficient normalization. Its optical length and
proper length enter different parts of the matching.

At a lossless field boundary, reflection and transmission satisfy

$$
|R(\omega)|^2+|T(\omega)|^2=1.
$$

The parent reconstructs what returns through $R$; the transmitted component
continues beyond its locally accessible chart. In this sense the shadow,
throat, and expanding interior are three views of one boundary problem. The
calculated complex Dirac reflection below determines the inherited child
covariance rather than serving as a visual analogy alone.

The project's horizon-penetrating Hamiltonian retains shift and spin connection.
Its declared Cauchy evolution conserves Dirac norm while transporting it into
the child ([tetrad](../results/nsc-4-dirac-tetrad.json),
[transport](../results/nsc-4-lorentzian-transport.json)). Section 4.8 evaluates
the parent-matched massless stress. Its negative radial null components,
pressure anisotropy and outward flux define the input to gravitational
backreaction. Complete continuation is assessed through the resulting
constraints, tidal response and finite-affine endpoints
([clock and horizon](../results/nsc-5-clock-horizon.json)). Physical bosonic
propagation requires the Lorentzian response associated with the Euclidean
spectral form factors [13].

<!-- nsc-claim:target-dark-sector -->
### 2.3 Adjacent rooms and gravitational response

Eliminating an unresolved room produces an exact Schur self-energy. The dark
sector proposal identifies background gravity and clustering as physical
projections of this boundary response. Their calculation requires the same
state's density, pressure, anisotropic stress and transfer currents.
A predicted vacuum curvature and physical scale connect the spectral
construction to the cosmological-constant question.

In momentum language,

$$
\Gamma_n^{(2)}(k)=K_{\rm local}(k)+\Pi_{\rm outside,\Theta}(k),
$$

with the proposed observational dictionary

$$
\Pi_{\rm outside}(0)\longrightarrow\text{background expansion},
\qquad
\Pi_{\rm outside}(k>0)\longrightarrow\text{clustering and lensing}.
$$

<!-- nsc-claim:target-measurement -->
### 2.4 Localized excitations, waves and detection

Quantum field theory describes particles as field excitations, and individual
detections build interference patterns [1]. Solitons and collective coordinates
supply additional construction patterns [4]; the project's
[particle/wave record](../results/nsc-1-s-one-particle-wave-identity.json)
imports that reduction. A detector model couples the identified excitation,
probe and environment through their interaction while retaining energy and
charge accounting. Its outputs are localized records, normalized outcome
probabilities, interference and two-detector correlations. Bell and PBR
results [2,3] constrain that probability construction.

<!-- nsc-claim:target-antimatter -->
### 2.5 Chiral mass and charge-conjugate states

The scalar-sheet benchmark in Section 3 gives an exact Dirac mass interaction
on an invariant subspace. The actual boundary problem determines its Clifford
content and physical sector. Chirality, sheet label and gauge charge have
separate operators. With the stated gamma conventions,
$C=i\gamma^2\gamma^0$ enters $\psi^c=C\bar\psi^T$, while the corresponding
column-spinor map is $\psi^c=i\gamma^2\psi^*$. Charge conjugation acts on the
gauge representation. Both complementary scalar-link sectors are retained
until the action or domain selects one. Spin-one fields and the matter–antimatter
asymmetry require their respective gauge and state dynamics.

The frequency picture is explicit in the quantized field:

$$
\Psi(x)=\sum_s\int d^3p\,
\left[a_su_se^{-ip\cdot x}+b_s^\dagger v_se^{+ip\cdot x}\right].
$$

The two phase orientations become particle and antiparticle sectors after
quantization. A real sine combines both,

$$
\sin(\omega t)=\frac{e^{i\omega t}-e^{-i\omega t}}{2i},
$$

while charge conjugation supplies the physical reversal of gauge
representation.

<!-- nsc-claim:target-recursion -->
### 2.6 Recursion, clocks and probability

The nested interpretation treats infinity as successive generation of rooms.
On the specified unweighted spatial chain, the endpoint is uniquely determined
for $\Omega q\leq1$; $\Omega q>1$ requires data at infinite depth
([tail result](../results/nsc-5-tail-limit.json)). Critical recursions can
converge algebraically. The physical construction couples this endpoint
analysis to a state measure, clock-transfer law and limiting observables.
Those quantities determine what an infinite chain means for duration and
probability within a room.

## 3. A precise shared mass and self-energy benchmark

In sheet-first Weyl order, consider

$$
H_8=I_2\otimes\boldsymbol{\alpha}\cdot\mathbf p+\Phi\tau_1\otimes\beta,
\qquad
\Pi_-=(I_8-\tau_3\otimes\gamma^5)/2.
$$

The projector selects parent-left plus child-right. Exact calculation gives

$$
[H_8,\Pi_-]=0,
\qquad
W^\dagger H_8W=\boldsymbol{\alpha}\cdot\mathbf p+\Phi\beta,
$$

where W embeds those four components. Thus one invariant sector is an ordinary massive Dirac Hamiltonian, with

$$
E^2=|\mathbf p|^2+\Phi^2,
\qquad
\Sigma_p(E)=\Phi^2(E+\boldsymbol{\alpha}\cdot\mathbf p)^{-1}.
$$

The complementary sector is equally invariant and spectrally equivalent. The common action or domain must supply physical sector selection. A spinor-identity sheet link instead has energies ±p±Phi and is gapless at p=Phi. A Hermitian pseudoscalar link has the same free gap and self-energy as the scalar benchmark. The spectrum therefore does not determine the complete interaction. The [observable bridge](../results/nsc-7-observable-bridge.json) checks these distinctions and the corresponding charge algebra.

Both chiral components carry the same Dirac charge. A massive positive-energy eigenstate has a constant chirality expectation. The familiar left-to-right oscillation belongs to a specified initially chiral state containing both energy signs. In the full paired angular continuum problem, the massless radial geometry and its radius pulse preserve physical chirality. Their radial gap is consequently not already a four-dimensional chirality-breaking rest mass.

The boundary calculation below determines the Clifford structure, orientation and invariant domain needed to connect this Dirac reduction to the actual geometry.

**evaluated boundary channels.** The full paired-angular calculation now evaluates an energy-dependent four-component spatial link in an explicit current-conserving domain. Direct joined, transfer and Schur calculations agree below 2e-15 relative; an independent continuum integration gives 0.00947%–0.03645% errors at the reported resolution. The link preserves physical chirality. Its action-form Clifford decomposition has vector and axial content; scalar and pseudoscalar terms vanish to numerical precision ([chiral boundary](../results/nsc-8-chiral-boundary.json)).

This finite spatial point-port response fixes the transport structure of the massless realization in its stated units and regulator. Chirality is transported with the spin frame. The scalar mass reduction therefore requires a physical interaction or domain mechanism, followed by its causal response and pole calculation.

<!-- nsc-claim:compact-mode-overlaps -->
**Imported result / repository derivation: compact mass identification.** Applying the established conformal Dirac and Kaluza–Klein reductions [33,34] to the declared free carrier on an interval of coordinate length 2 L⋆, with the right-handed component vanishing at both endpoints, yields one Weyl zero mode and the massive Dirac tower

$$
m_n=\frac{n\pi}{2L_{\star}},\quad n\geq1.
$$

The full spacetime warp cancels from canonically normalized free evolution. The declared compact size fixes this free spectrum; species, occupations and the interacting throat coupling enter the common source problem ([compact mass map](../docs/nsc-compact-mass-map.md)).

**torsion interaction overlaps.** Substituting those NSC profiles into the published higher-dimensional Einstein–Cartan contact [35] produces a mode-matrix vertex. After fermionic antisymmetrization the 3111 coupling remains nonzero, so a single massive level is not an exact interacting sector. The geometric warp changes the overlap coefficients even though it cancelled from free evolution. The five-dimensional stiffness $s_T$, compact boundary action, physical state and absolute stress remain inputs ([compact interaction](../results/development/compact-interaction.json)).

<!-- nsc-claim:charged-sector-parity -->
**Imported result / repository derivation: charged parity and anomaly compatibility.** To use the imported charged throat as a source test, a local realization must contain an actual gauge connection. The present candidate is a U(1) extension with two equally charged bulk copies; it does not derive electromagnetism or the internal algebra. Within this field content, opposite compact parities are required both for an even neutral scalar link and for cancellation of the zero-mode cubic and mixed gravitational–gauge anomalies. The two allowed projectors are the invariant scalar-link sectors identified above:

$$
P_{\mathrm{even}}=\frac{I_8\mp\tau_3\otimes\gamma^5}{2}.
$$

At the zero-mode level captured by the diagnostic, the potential for Callan–Harvey anomaly inflow at the spatial boundary is explicitly neutralized by the field content. The two equally charged bulk copies with opposite compact parities, selected by $P_{\mathrm{even}}=(I_8\mp\tau_3\otimes\gamma^5)/2$, give one left- and one right-handed zero mode and make both the cubic $U(1)^3$ and mixed gravitational–$U(1)$ boundary-anomaly coefficients vanish exactly. This exact cancellation strictly preserves $U(1)$ gauge symmetry across the throat boundary in the [`charged-sector.json`](../results/development/charged-sector.json) diagnostic without an additional localized Wess–Zumino boundary term; the full five-dimensional determinant phase and boundary-domain matching remain the separate completion item recorded there.

Both remain. Sheet swap is not charge conjugation: a right-handed charge-one field is counted as a left-handed charge-minus-one field only in anomaly bookkeeping ([charged sector](../results/development/charged-sector.json)). The imported Einstein and Casimir solutions are not recomputed. Link mass, flux geometry, renormalized coefficients and child expansion remain unresolved.

The two five-dimensional Dirac copies yield one four-dimensional Dirac zero field before a link mass is added. Their ultraviolet spinor trace has rank eight; applying the bulk coefficients to this candidate requires that two-copy multiplicity, not the one-copy normalization of the earlier UV control.


## 4. Geometry, vacuum stress and an explicit energy ledger

The isolated radial throat is gapless. Periodic confinement can create a band gap, and a smooth periodic profile removes the earlier geometric seams:

$$
r(x)=\sqrt{a^2+[\sin(kx)/k]^2},\qquad k=\pi/(2R).
$$

At a=1 and R=2,4,8, its continuum first band edges are approximately 0.76423204, 0.51425546 and 0.28454085 ([smooth geometry](../results/nsc-4-smooth-geometry.json)). R remains an input. These are spatial spectral values, not identified particle masses.

The same profile requires a negative effective Einstein null source,

$$
8\pi G(\rho+p_x)=-2r''/r.
$$

The [covariant vacuum controls](../results/nsc-4-covariant-measure.json) establish negative axial-null vacuum response on a specified compact constant cylinder. The [varying-neck calculation](../results/nsc-4-shape-response.json) differentiates all metric functions and computes a periodic-minus-antiperiodic stress difference. Its local finite counterterms cancel in that difference; its absolute stress is not thereby fixed. Choosing a spin structure because its stress has the desired sign would not derive it.

A finite term proportional to the squared Weyl curvature leaves the anomaly unchanged but changes the varying-neck shape equation. The current flat common-scale integral also diverges in its finite realization. The ultraviolet functional, complete compensator and state must resolve these questions before a stationary root represents self-sourcing.

**finite metric response.** The general-metric calculation now isolates four independent bulk channels in the local basis $M^4,M^2R,C^2,R^2,E_4,\Box R$. On the closed smooth cell, the Euler and divergence terms have exact zero bulk variation. Both R=2 and R=4 profiles distinguish the remaining four coefficients; the constant-cylinder control has rank two. Independent lapse, radial metric and radius variations agree with direct energy changes ([finite terms](../results/nsc-8-finite-terms.json)).

For the declared finite-action sign convention, the exact neck null response is

$$
\Delta(\rho+p_x)_{\mathrm{fin}}=\frac{4M^2}{a^2}c_R+\frac{32(1+a^2k^2)}{3a^4}c_C-\frac{16(1+4a^2k^2)}{a^4}c_{R^2}.
$$

This gives a concrete ultraviolet matching problem: determine these coefficients from the same functional, then solve the independent metric equations. They have not been fitted to make the neck stationary. The result concerns a closed torsionless metric sector through four derivatives; boundary contributions and other fields require their corresponding terms.

<!-- nsc-claim:torsion-uv-finite-matching -->
**five-dimensional ultraviolet matching.** The leading proper-time bulk coefficients of the candidate torsion contact reuse the dimension-independent heat formulae [29,36]. In five dimensions the induced Einstein and $K^2$ terms have stiffness ratio 9. A matching window keeps this induced piece separate from the complementary one-loop determinant. The volume term is retained; it is not set to zero to manufacture an asymptotically flat solution. Leading stiffness 9 does not fix the finite datum $c_T-9c_R$. The displayed derivative scale is a coefficient comparison, not a physical torsion mass ([UV map](../results/development/torsion-uv-map.json)).

<!-- nsc-claim:published-flow-compatibility -->
**Imported result / repository derivation: quantum-flow limits.** A published generalized Palatini flow [37] has a larger connection field space. Its metric contribution does not preserve the constrained Dirac three-form subspace in $d=5$; its fermion contribution does. Copying that paper's four-dimensional fixed-point numbers does not determine the present five-dimensional functional. The next flow must be defined on the actual metric, skew three-form, physical fermions, links and measure ([flow compatibility](../results/development/flow-compatibility.json)).

**covariant determinant source.** The static Euclidean operator on the same smooth cell is

$$
D_\omega=\beta\omega/N+i\beta H_0(q,r).
$$

Physical four-dimensional chirality anticommutes with $D_\omega$. Its square retains the lapse gradient $\alpha\,\omega N'/(qN^2)$. The cutoff modulus per coordinate time is

$$
E_\Lambda=\frac{1}{2}\sum_{\kappa\geq1}2\kappa\int_{\mathbb R}\frac{d\omega}{2\pi}\mathrm{Tr}_{x,4}E_1(D_\omega^2/\Lambda^2).
$$

Frequency integration reproduces the known cylinder vacuum response. Independent metric variations conserve the static stress and recover the established Dirac heat coefficients $a_0=4$, $a_2=-R_E/3$ through $a_4$ [29]. The specified ultraviolet subtraction is

$$
E_{\mathrm{sub}}(M)=E_\Lambda-\frac{A_0\Lambda^4/2+A_2\Lambda^2+A_4\log(\Lambda^2/M^2)}{32\pi^2}.
$$

At $a=1$, $R=2$ and $M=1$, the remainder energy is about $-0.09050$ with neck null stress about $-0.00911$ at $\Lambda=4$. This is the negative null contribution in the stated subtraction prescription ([covariant source](../results/nsc-9-covariant-source.json)). On the same lapse/radius/shape probes the determinant now supplies actual projections $b_A^{\mathrm{determinant}}$ in

$$
b_A^{\mathrm{determinant}}+J_{Ai}c_i+b_A^{\mathrm{compensator,link,state}}=0.
$$

The coefficient and determinant terms retain their common normalization. Uniform clock scaling multiplies the covariant energy, while a fixed spatial-Hamiltonian cutoff does not.

**normalization profiles.** Write $x=d^2/\Lambda^2$. The calculated modulus per eigenvalue is half the exponential integral E1(x). A smooth weighted logarithm and the heat-covariantized old rank term obey

$$
g_{\mathrm{log}}=-\frac{1}{2}e^{-x}\log(d^2/M^2),\qquad
g_{\mathrm{log}}=g_{\mathrm{PT}}+[\log(M/\Lambda)+\gamma_E/2]e^{-x}+h(x),\qquad
h(x)=-\frac{1}{2}[e^{-x}(\log x+\gamma_E)+E_1(x)].
$$

The remainder is not a multiple of the identity or of the heat kernel, so replacing the old finite-rank term by a heat trace changes the finite action. The identity is a comparison of prescriptions, not an insertion into $\Gamma_{\mathrm{one}}$ ([normalization](../results/nsc-10-measure-normalization.json)). Andrianov’s primary projector $P_N=\Theta(1-D^2/\Lambda^2)$ is a different object [11]; its cylinder cutoff wall has a closed form, while an unrestricted $\mathrm{Tr}1$ integral diverges.

Assign half of the cross-region link energy to each region, with $h_A=(P_AH+HP_A)/2$, and use $C_{ij}=\langle c_j^\dagger c_i\rangle$. Then

$$
\frac{dE_A}{dt}=\mathrm{Tr}(C\,i[H,h_A])+\mathrm{Tr}(C\dot h_A).
$$

For the static vacuum, $[H,C]=0$, so continuing net regional injection is zero. A prepared excitation transports conserved energy through the same geometry-derived link ([energy transfer](../results/nsc-6-energy-transfer.json)). Its current depends on occupations and correlations as well as the retarded map.

A prescribed smooth radius pulse gives exactly $H(t)=H_0+\epsilon g(t)V$. It creates Dirac particle–antiparticle excitations. Independent evolution and spectral transition calculations agree, and supplied geometric work accounts for the excitation energy ([vacuum work](../results/nsc-6-vacuum-work.json)). At epsilon=0.02 and tau=0.5, one representative kappa=1 channel produces approximately 1.40181e-5 pairs and energy 3.46640e-5 in throat units. These values describe one channel at the stated pulse amplitude and duration. A uniform clock rescaling produces no pairs.

<!-- nsc-figure:vacuum-work -->

The pulse demonstrates conversion of geometric work into Dirac excitations. The coupled metric equation determines the deformation when its source is evolved dynamically.

**normalized histories and geometric noise.** The same covariant operator supplies the canonical Hamiltonian $H=N^{1/2}H_0 N^{1/2}$. For a Gaussian fermionic covariance $0\leq C\leq I$ and two history unitaries, the normalized closed-time-path amplitude is

$$
Z[+,-;C]=\det(I-C+C U_-^\dagger U_+),\qquad \Gamma_{\mathrm{IF}}=-i\log Z.
$$

The Fock-space trace identity is due to Klich [30]; the present record implements it for this geometry vertex and checks mixed and pure states against an explicit eight-dimensional Fock calculation. Equal histories give $Z=1$. Causal contour derivatives reproduce the retarded commutator response, and the symmetrized noise matrix is positive. For the Gaussian time window already used in the pulse calculation, the smeared vertex variance equals the leading pair-production coefficient:

$$
\mathrm{Var}(V_f)=N_{\mathrm{pairs}}^{(2)}/\epsilon^2.
$$

At $\tau=0.5$ this common weight is $0.03533459211$ on the refined spatial grids ([influence](../results/nsc-10-influence.json)). Bare mean forces and equal-time fluctuations do not converge in the same way; they remain finite-regulator quantities, not renormalized absolute stress. Multiplying by a real local-action phase preserves equal-history normalization while changing the mean source, so $Z=1$ does not select the physical finite terms. The influence approach to gravitational response has prior work by Martín and Verdaguer [31]. The contribution here is the executable connection to this smooth Dirac geometry and its previous energy/pair calculation.

**matching static and causal response.** The radius channel now connects the Euclidean frequency operator to the Lorentzian Hamiltonian response. Independent matrix integration of the ordinary finite Dirac loop equals the spectral vacuum susceptibility at imaginary external frequency, to about $10^{-15}$ on the tested finite grids. The finite proper-time Hessian is calculated separately, with its regulator retained.

The periodic-minus-antiperiodic difference removes state-independent local metric terms on identical geometry. It permits the continuum matching control

$$
K_E^{P-AP}(\nu)\longrightarrow\chi_R^P(i\nu)-\chi_R^{AP}(i\nu).
$$

The independent retarded sum uses a positive Abel transition weight and a stated cutoff extrapolation. These finite regulators are distinct. Through angular label $\kappa=8$, the two methods give:

| Imaginary frequency | Euclidean relative radius response | Retarded extrapolation minus Euclidean |
|---|---|---|
| 0 | 0.185771550 | -1.49e-9 |
| 0.6 | 0.172525449 | -1.41e-9 |
| 1.4 | 0.130625412 | -1.14e-9 |

Spatial, angular, internal-frequency and cutoff refinements are recorded separately ([response matching](../results/nsc-11-response-matching.json)). The omitted angular tower has no proved remainder bound. A change from inverse-radius to log-radius coordinates adds the required contact $\langle H''\rangle$; its relative value is nonzero and is retained. With that term, independent metric differentiation agrees with the earlier covariant source Hessian. This is a verified matter-response connection in the smooth ultrastatic radius channel. Absolute finite coefficients, finite-cutoff causal completion and constrained metric backreaction remain to be determined.

### 4.1 One source must satisfy the geometric boundary equation

<!-- nsc-claim:compact-boundary-adjoint -->
**Repository derivation using established boundary calculus.** The compact chiral domain fixes the leading Einstein boundary term of the determinant magnitude. It is isotropic for the Lorentzian Hamiltonian current, but the Euclidean Dirac operator requires complementary adjoint data. Using $D_E^\dagger D_E$, the known mixed heat coefficients [29] yield the Einstein/Gibbons–Hawking–York ratio 2, no leading boundary-volume term, and a specified geometric $a_3$ contribution ([compact boundary action](../results/development/compact-boundary-action.json)).

The Gaussian compact endpoints carry nonzero geometric momentum. The common metric variation therefore requires

$$
\pi_p^{ab}+\pi_c^{ab}+\frac{\delta\Gamma_{\mathrm{rest}}}{\delta h_{ab}}=0.
$$

Here $\Gamma_{\mathrm{rest}}$ excludes the geometric terms already included in $\pi_p,\pi_c$. The variation uses established gravitational boundary calculus [41]. Reflecting compact endpoints and a transmitting parent–child throat are distinct domains.

### 4.2 A finite Dirac interaction on the curved cell

<!-- nsc-claim:curved-compact-source -->
The published massless interval determinant [38] can be applied to the existing covariant four-dimensional operator. For two Dirac copies and the untwisted transverse interval of conformal length $\ell=2L_\star$, its finite endpoint interaction is

$$
\Gamma_{\mathrm{int}}=-\frac{1}{2}\mathrm{Tr}\log\left(1-e^{-2\ell\sqrt{D_4^2}}\right).
$$

The half-rank zero level and full massive-level multiplicities are retained. Differentiating this function of the same $D_4$ supplies independent lapse, radial-metric, radius and separation responses. These effective four-dimensional sources follow directly from the Dirac determinant, integrated over the transverse interval.

For the unit-neck ultrastatic cell with axial circumference 4 and transverse length 2, the compact interaction gives

$$
(\rho+p_r)_{\mathrm{int}}=-8.9162\times10^{-4}.
$$

The units belong to the recorded $g_4$ geometry. Matrix variation agrees with independent energy differences; an independent Bessel-frequency integral and separate spatial/angular/frequency refinements agree. The local conservation residual decreases with resolution ([source and state record](../results/development/compact-casimir.json)).

The first two local curvature terms of this finite interaction are derived and kept in the same energy account. A decompactified bulk determinant is generally nonlocal in curved $g_4$; it cannot be discarded using a scaleless flat-space subtraction. Sections 4.5–4.6 now provide the full free compact cutoff weight and its light-field matching for the declared reflecting domain. The common completion, state and transmitting phase remain necessary for the physical metric equation.

### 4.3 A conditional state-selection mechanism

<!-- nsc-claim:conditional-holonomy-saddle -->
The null-source sign depends on global geometry and state. With the same transverse length and unit neck, enlarging the axial circle to 8 changes the interaction's neck null source to $+1.6402\times10^{-4}$. Periodic axial data can also reverse the sign. These values resolve the dependence on the prescribed cell geometry and state.

Where the charged ultrastatic geometry contains a closed axial return path, a flat U(1) connection has a physical Wilson line. Its effective spatial phase combines the spin structure and gauge holonomy:

$$
\alpha=\eta_{\mathrm{spin}}-\frac{1}{2\pi}\oint A_x\,dx.
$$

The phase is defined modulo one. Local gauge-invariant terms do not distinguish these flat connections. The general boundary determinant [39,40], applied to the positive radial partner operators, makes each free-fermion contribution minimal at effective AP phase, $\alpha_\star=1/2$. The full free-KK determinant independently gives an AP-minus-P energy of about -0.16995535 and positive second variation 3.19343937 in the recorded cell. This is the one-loop minimum of the real potential for a freely varying holonomy on the stated closed cell.

<!-- nsc-figure:compact-source -->

The Casimir and Wilson-line mechanisms are established physics [38–40]. The project calculation matches them to this operator and identifies their source, normalization and domain dependence. The recursive geometry must supply the actual return path; an unwrapped chain acquires no Wilson loop merely because the development cell has one.

### 4.4 Apply the source to the unwrapped horizon

<!-- nsc-claim:derived-horizon-source -->
**Repository application of established Dirac and conformal-state results.** The unwrapped transition patch has spatial form R times S2 times the transverse interval. It has no non-contractible radial cycle supplying the finite cell's flat Wilson line. A declared magnetic flux instead supplies known angular zero modes: one four-dimensional Dirac field gives $|q_{\mathrm{mag}}|$ complex two-dimensional fields [32,43]. The integer flux is not selected by this calculation.

For this free massless sector, Unruh boundary data and the stored horizon surface gravity give the parent Killing power [44]

$$
P_\infty=\frac{|q_{\mathrm{mag}}|\kappa_h^2}{48\pi}.
$$

At $\kappa_h=0.2383257996$, this is $3.76661306\times10^{-4}$ per unit $|q_{\mathrm{mag}}|$, in the recorded natural units. The local horizon limit is regular under the specified state conditions. Converting this power to a child density rate requires proper volume, time and deposition; interactions, a light-field mass and global state enter the corresponding source problem ([horizon/source record](../results/development/horizon-source.json)).

The result tests an actual source requirement. The sector's horizon null stress is negative, but both null contractions at the old imposed neck are positive, where that geometry requires negative ones. Static magnetic and potential terms have zero radial null contraction. The minimal sector therefore does not source that prescribed neck by itself. The reduced metric equations keep the missing source U(r) explicit and propagate their constraint; they do not choose U from the desired geometry. A nonzero outgoing flux also requires evolution or compensating flux before it can belong to an exactly static solution.

### 4.5 The compact quantum warp supplies a spectral weight

<!-- nsc-claim:derived-warped-source -->
The classical conformal reduction preserves the canonical compact mass tower. Its quantum determinant has a different normalization. For the conformal metric

$$
g_{5,s}=e^{2s\sigma(Y)}(g_4+dY^2),
$$

transport to fixed L2 measure gives

$$
D_s=e^{-s\sigma/2}D_0e^{-s\sigma/2},\qquad L_s=D_s^\dagger D_s.
$$

The compact chiral values and complementary adjoint domain are retained. In the canonical field, the eigenproblem uses both weights, $\int e^{-s\sigma}|D_0\chi|^2=\epsilon\int e^{s\sigma}|\chi|^2$. Its compact eigenvalues define

$$
h_{\Lambda,s}(y)=\sum_j E_1\!\left(\epsilon_j(y;s)/\Lambda^2\right),
\qquad
\Gamma_{5,\Lambda}=\frac{1}{2}\mathrm{Tr}_4 h_{\Lambda,s}(D_4^2).
$$

The two complementary five-dimensional copies are already included. This is a finite-cutoff free Euclidean determinant magnitude, with a separate prescription required for actual spacetime zero modes. It composes with the existing D4 rather than replacing it. The physical compact warp is calculated at s=1; the numerical homogeneous application uses $g_4=\mathbb R^2\times S^2(r)$.

For $q_{\mathrm{mag}}=1$, r=1, transverse length 2 and cutoff 2, the potential per reference two-dimensional area changes from 2.218882505 to 2.154980560. This is a 2.88% change in that regulated contribution. Independent first variations and a dual-heat conformal identity agree within the stated Galerkin convergence. The classical mass tower is preserved while the regulated geometric potential changes ([warped-source record](../results/development/warped-source.json)).

### 4.6 One coefficient account for the retained light field

<!-- nsc-claim:derived-compact-matching -->
Separate one canonical light Dirac field at a matching cutoff nu by defining

$$
H_\nu(y)=h_\Lambda(y)-E_1(y/\nu^2),\qquad
\Gamma_{5,\Lambda}=\Gamma_{\mathrm{light},\nu}+\Gamma_{H,\nu}.
$$

The compact proper length fixes the finite limit $H_\nu(0)$. Two remaining spectral moments determine the vacuum and Einstein terms. Applying the established four-dimensional heat/Mellin formulas [29,45] gives the complementary Euclidean density

$$
\mathcal L_{H,E}=V_D-A_D R_E+C_{F,D}F^2+\cdots.
$$

At cutoff 2, transverse length 2, s=1 and matching cutoff $\nu=1$:

| Dirac contribution | Coefficient | Units |
|---|---:|---|
| Vacuum $V_D$ | 0.1725408633 | inverse length to the fourth power |
| Einstein $A_D$ | 0.00550228445 | inverse length squared |
| Gauge $C_{F,D}$ | 0.00483725070 | dimensionless |

The curvature-squared terms are fixed by the same $H_\nu(0)$. The Dirac a4 combination contains no independent R-squared term; that does not eliminate an undetermined finite term of the completed theory. The source map retains the established static curvature-sign convention. Its local truncation needs small external curvature, derivatives and gauge field; that approximation has not been established at the imposed neck.

Changing nu transfers terms between the light field and its complement:

$$
\partial_\nu\Gamma_{H,\nu}
=-\partial_\nu\Gamma_{\mathrm{light},\nu}.
$$

The numerical vacuum, Einstein and gauge changes cancel. The matching identity fixes how these contributions combine. The complete Newton constant, gauge coupling and U(r) follow after the common functional and causal state are matched ([compact/light matching](../results/development/compact-matching.json)).

<!-- nsc-figure:source-matching -->

The source equation now has a defined interface: combine the computed free determinant with the contribution required by the common measure, compensator, transmitting interaction and state. These records test distinct controlled settings; their assembly into one Lorentzian solution remains part of physical closure. In particular, the spherical Einstein term already present inside the determinant must be counted once. No additional independently weighted gravitational action is introduced.

The next subsections follow that interface rather than nine separate programmes. A definite interaction among already retained light fields, a first-order spherical geometry owner, and a local curvature correction fix the low-energy source language. At the stored child curvature the local truncation is not controlled, so the original spectral weight is used. Homogeneous free-state differences then dilute, a canonical massless reference is distinguished from a parent-matched candidate, and the finite proper-time factor is shown to be state-dependent. The following calculations apply the cited methods with the shared coefficient account.

### 4.7 A retained interaction, a spherical owner and a curvature contact

<!-- nsc-claim:derived-gauge-source -->
The retained Dirac and spherical-photon fields already supply a projected interaction in the magnetic lowest-Landau-level sector. Established two-dimensional bosonization [56] reduces the charged current to one collective mode with screening scale $\mu^2(r)=|q| g_4^2/(4\pi^2 r^2)$. That mode is not an inserted fundamental scalar, a resolution compensator or a predicted four-dimensional rest mass. The physical matched $g_4$ is not selected here, and the Dirac complement coefficient $C_{F,D}$ is not its complete inverse coupling.

Outside the recorded horizon the collective potential is nonnegative and vanishes at both scattering ends, so a constant-mass Boltzmann factor is the wrong flux calculation. With the specified Unruh data and no incoming parent bath, the projected Killing power is the free central-charge term plus the charged Bose integral of the transmission. Reference [44] supplies the free flux normalization; the charged transmission integral is evaluated in the gauge-source record. At the stored $\kappa_h=0.238325799634$, for $q=1$ and unselected test couplings, $g_4=0.5$ gives parent Killing power $0.0003631182$, or $0.964044$ of one free LLL channel. These are numerical controls at the stated test couplings. Stationary current conservation and future-regular Painlevé–Gullstrand components convert that power into a horizon null source $T(K,K)_h=-P_\infty/(4\pi r_h^2)$. A nonzero net power still requires evolution or compensating flux before it can belong to an exactly static geometry ([gauge-source record](../results/development/gauge-source.json)).

<!-- nsc-claim:canonical-spherical-action -->
The local two-derivative spherical sector of the same action now has an explicit first-order owner. Spherical reduction of $S_4=\int\sqrt{|g_4|}[-A_4 R_L-V_4-C_F F^2]$ with dimensionless area $X=8\pi A_4 r^2$ recovers the Grumiller–Kummer–Vassilevich first-order theory [46], including the area-dependent charge potential. The auxiliary torsion of that first-order patch is not the independent five-dimensional skew torsion of the earlier Einstein–Cartan candidate. Coupling the existing charge mode closes the classical constraint algebra, including the $1/X$ area force that is not the factorized matter potential printed in that reference. This identifies an action-specific canonical measure and ghost operator for the reduced theory. It replaces a guessed one-dimensional scale weight; it does not quantize the full nonlocal spectral action, select physical couplings, or evolve a new geometry ([spherical-action record](../results/development/spherical-action.json)).

<!-- nsc-claim:curvature-eft-source -->
The known local curvature-squared contribution can be retained in that same spherical system as a first-order effective-field-theory contact. A checked metric redefinition [47,48] converts the bulk $C^2$ and $R^2$ terms into

$$
\Delta\mathcal L_m=-\frac{c_W}{2A^2}\left(T_{\mu\nu}T^{\mu\nu}-\frac{T^2}{3}\right)-\frac{c_R}{4A^2}T^2
$$

without assigning extra initial data to a fourth-order truncation. The charge/Maxwell source, including its angular pressure, supplies a definite spherical contact. The first-order constraints still close at this EFT order. A change of variables does not cancel physical vacuum curvature: pulling the map back to the original metric removes the apparent cosmological-constant shift. The construction requires $|c_i R|/A$ and $|c_i T|/A^2$ small. That regime is not established at the imposed neck, and the nonlocal remainder $h(D_4^2)$ still belongs to the theory. This is not a replacement of the spectral action or a proof of physical-mode health ([curvature-EFT record](../results/development/curvature-eft.json)).

<!-- nsc-claim:full-spectral-endpoint -->
**Imported spectrum / repository composition.** The stored development geometry reaches curvatures at which that local expansion is not controlled. At the recorded horizon, neck and asymptotic child, the largest sectional magnitudes divided by $\nu^2=1$ are $0.2168$, $4.7124$ and $9.4248$. Composing the existing compact weight $h_\Lambda$ with the four-sphere Dirac spectrum [49] therefore supplies the usable Euclidean endpoint response. At the child-radius control $a=0.325735$, the full logarithmic-radius derivative is $0.001239600$, while the local truncation through $a_4$ gives $-0.164273768$. A local stationary radius near $0.43742$ has full derivative $0.084326$, so it is not stationary in the same determinant. Near agreement at larger $a$ supports the expected low-curvature use of the local coefficients. The round sphere is not the global continuation of the black-universe chart, whose interior has Kantowski–Sachs spatial topology. Use the original spectral functional at these curvatures; do not accept the local root as a solution of $\Gamma_5$ ([spectral-endpoint record](../results/development/spectral-endpoint.json)).

### 4.8 Canonical reference, parent matching and the neck budget

Three objects must be kept distinct: a renormalized massless reference, a parent-matched state on the actual domain, and the finite-$\Lambda$ full action. Asymptotic tail formulae are not substituted at the neck.

<!-- nsc-claim:child-state-decay -->
**Repository derivation using established Dirac stress differences.** For the fixed free compact tower on the actual expanding interior, a defined class of homogeneous state differences has vanishing stress density. Trace-norm conservation and Hölder's inequality give

$$
|\Delta\rho|\leq\frac{M+K/r}{4\pi s_{\min}r^3},
$$

with $O(r^{-4})$ for massless modes ($M=0$) and an $O(r^{-3})$ envelope when finite massive occupation moments are present. Expansion reduces the local stress through redshift and proper-volume growth; it does not erase the state. The free LLL Unruh-minus-Hartle–Hawking difference decays as $r^{-4}$ in the interior, with asymptotic $r^4$ density $-\kappa_h^2/(576\pi^3)$. The bound applies to the specified free-state differences; the absolute reference stress is retained separately ([child-state record](../results/development/child-state.json)).

<!-- nsc-claim:absolute-massless-reference -->
**Imported Hadamard/conformal construction / repository application.** A massless reference can be specified on the child tail by smoothly continuing the conformal metric $\bar g=g/r^2$ to an auxiliary ultrastatic region and propagating its ground covariance back [26,52,53]. In the inherited proper-time convention the leading absolute stress on that tail is

$$
\lim\langle T^\mu{}_\nu\rangle_{\mathrm{ref}}=\frac{11 H_{\mathrm{child}}^4}{960\pi^2}\delta^\mu{}_\nu,
\qquad H_{\mathrm{child}}=\sqrt{3\pi},
$$

Expanding that same conformal source gives the leading radial null combination $-u^2/16+O(u^3)$. That expansion is not valid at the neck, where $u$ is not small. Converting the reference into the adopted finite cutoff requires

$$
\Gamma_5=\Gamma_{\mathrm{light,ren}}+\left(\Gamma_{H,\nu}+\mathcal C_{\nu,\mu}\right).
$$

The pieces depend on normalization; their sum preserves the original functional. Adding the renormalized reference to an unconverted complement would change $\Gamma_5$. The conversion preserves the specified full functional ([massless-reference record](../results/development/massless-reference.json)).

<!-- nsc-claim:finite-angular-reference -->
The same reference is then evaluated at finite child radius, retaining the full angular tower, covariant subtraction, conformal anomaly transport and the energy/work identity. Units remain $\hbar=c=L_{\mathrm{throat}}=1$ and $\mu=1$. At the neck the auxiliary reference is $\rho=0.7101684$, $p_{\parallel}=0.5026775$, $p_{\perp}=0.0376558$, hence $\rho+p_{\parallel}=1.2128459$. After extracting one healthy Einstein term $a_{\mathrm{EH}}=A_{\mathrm{EH}}L_{\mathrm{throat}}^2>0$, the imposed neck requires $L_{\mathrm{throat}}^4 T_{\mathrm{total},kk}=4a_{\mathrm{EH}}(1-3\pi/2)=-14.8495559\,a_{\mathrm{EH}}<0$. The calculated reference is positive and cannot supply that null requirement by itself. A pure vacuum term proportional to $g_{\mu\nu}$ cannot change the null component. No reference width was varied to obtain a stress sign ([angular-stress record](../results/development/angular-stress.json)).

<!-- nsc-claim:parent-matched-source -->
**Imported characteristic-state prescription / repository calculation.** Combining the standard fermionic Unruh boundary data [50] with the complex reflection amplitude of the actual exterior Dirac operator produces a definite massless source candidate on the expanding child. Affine-horizon occupation and empty incoming parent modes are compressed with the unitary scattering matrix; the factor $i$ in the real interior spin frame is required by the known zero-frequency Rindler projector. The state was specified by those horizon and parent data before the neck sign was computed. In the same canonical heat convention the $N=32$ tail-corrected neck values are:

| Component | Value in the declared throat units |
|---|---:|
| $\rho$ | $-0.0033068962$ |
| $p_\parallel$ | $-0.0494945966$ |
| $p_\perp$ | $-0.0429957713$ |
| radial null $+$ | $-0.0528075900$ |
| radial null $-$ | $-0.0527953956$ |
| outward parent Killing power | $0.0001422206795$ |

The stress entries have dimensions $L^{-4}$; Killing power has dimensions $L^{-2}$ with $\hbar=c=1$. Angular refinement from $N=16$ to $32$, after the sixth-order tail correction, changes the tensor by less than $4\times10^{-6}$ across the recorded radii. Frequency, quadrature and time refinements are assessed separately. These are measured sensitivities with an asymptotic angular-tail approximation, not a rigorous error bound on the infinite tower.

 Both radial null contractions are negative after including the flux. The two-derivative geometry with $a_{\mathrm{EH}}>0$ instead requires

$$
(\rho,p_{\parallel},p_{\perp},T_{\hat T\hat z})
=\left(2a_{\mathrm{EH}},\;2a_{\mathrm{EH}}(1-3\pi),\;2a_{\mathrm{EH}}(1-3\pi),\;0\right).
$$

The computed density has the opposite sign to that positive-Einstein requirement, the pressures are unequal, and the state carries outward parent Killing power $0.0001422206795$. Matching only the null equation by choosing $a_{\mathrm{EH}}$ would conceal those residuals. Leave $a_{\mathrm{EH}}$ symbolic. The remaining common-action terms must supply the componentwise difference, or the geometry must evolve under the same action. The source uses the stated canonical massless normalization and domain. A cosmological density rate additionally requires deposition, proper volume, clock matching and the evolved geometry ([parent-matched source](../results/development/unruh-state.json)).

<!-- nsc-figure:parent-source -->


### 4.9 Finite-cutoff conversion is state dependent

The continuation of a regular black-hole geometry through a positive-radius throat into an expanding child region is used here as an established classical benchmark. Its quantum source cannot be obtained by simply carrying the Euclidean spatial determinant through the horizon: the horizon-penetrating spatial Dirac Hamiltonian loses ellipticity there, as described by Finster and Röken [21] and witnessed explicitly for this carrier in Appendix B.1. Consequently, a proper-time Euclidean heat trace cannot be naively analytically continued into the trapped Lorentzian region to define the finite stress. This obstruction is precisely why the construction abandons the static Euclidean trace here and matches the causal Lorentzian closed-time-path functional to the state-dependent reference.

The retarded causal matching constraint is

$$
K_E^{P-AP}(\nu)\longrightarrow\chi_R^P(i\nu)-\chi_R^{AP}(i\nu).
$$

It carries the finite terms into the causal assembly

$$
\Gamma_{\mathrm{one}}^{\mathrm{CTP}}
=\Gamma_{\mathrm{ref,ren}}^{\mathrm{CTP}}
 +\left(\Gamma_{H,\nu}^{\mathrm{CTP}}+\mathcal C_{\nu,\mu}^{\mathrm{CTP}}\right)
 +\Delta\Gamma_{\mathrm{state}}^{\mathrm{CTP}}
 +\Gamma_{\mathrm{remaining}}^{\mathrm{CTP}}.
$$

At the recorded finite endpoints, this transition agrees with the response values in [`nsc-11-response-matching.json`](../results/nsc-11-response-matching.json) and the state-dependent conversion in [`state-regulator.json`](../results/development/state-regulator.json), without introducing another finite term.

<!-- nsc-claim:state-regulator-conversion -->
**Imported thermal determinant / repository compatibility check.** The adopted finite Euclidean proper-time factor cannot be identified, without a state-dependent conversion, with the thermal source of the canonical Dirac field used above. On flat $\mathbb R^3$ times an antiperiodic Euclidean thermal circle, with one ordinary massless Dirac field, $T/\nu=0.5$ gives canonical $\rho=0.0719658654\,\nu^4$ versus raw finite $\rho=-0.00122667356\,\nu^4$ after subtracting that prescription's own zero-temperature energy [51]. Local vacuum coefficients cancel in the thermal-minus-vacuum difference and cannot repair it. A pole-residue argument is insufficient: the regulated inverse has residue one on shell, yet the finite Matsubara determinant differs from the canonical CAR excitation energy. This flat equilibrium calculation fixes the conversion requirement for that state and prescription.

The subsequent [canonical causal functional](../docs/nsc-causal-common-functional.md)
assigns real-time response to the Lorentzian Dirac state while retaining
state-independent induced spectral coefficients once. Equal-history
normalization, covariance trace, unitarity, and finite energy/work balance
agree at approximately $10^{-16}$. The raw proper-time kernel remains a
distinct Euclidean prescription rather than being relabeled as a retarded
correlator. The [state-regulator record](../results/development/state-regulator.json)
provides the equilibrium conversion control.

## 5. Recursion and cosmological predictions

At a common dimensional energy, parent-normalized first-order elimination gives

$$
\Gamma_p(x)=K_p(x)-\frac{1}{\Omega}b\,\Gamma_c(x/\Omega)^{-1}b^\dagger.
$$

Here $b=B_{\mathrm{dim}}/\Lambda_p$. The inheritance scale $\zeta=(\Lambda L_\star)^2$ and cutoff ratio $\Omega=\Lambda_c/\Lambda_p$ remain distinct. Neither is selected by unit dilation alone. Finite chains and spatial endpoint estimates are verified; the following charged-radius calculation binds them conditionally in a fixed flux sector, while physical flux selection and the full recursive domain remain open.

**charged radius/scale binding.** Keep the massless charged lowest-Landau field explicit and evaluate the compact Wilsonian coefficients above MMP's physical magnetic scale $\mu_B=\sqrt{|q|}/r_e$. With the relational value $V_{\rm full}=0$, set $r_e=L_\star=1$ and solve

$$
\frac{q^2 C(\Omega,\mu_B)}{4A(\Omega,\mu_B)}=1,
\qquad \zeta=\Omega^2.
$$

For the adjacent fixed sectors $q=2,3,4$, the first radius root that also places the stored child curvature below the cutoff is the $q=4$ development branch,

$$
\Omega=3.973074368754331,
\qquad \zeta=15.785319939652627.
$$

The radius residual is $-1.55\times10^{-15}$ and the compact-resolution change in $\Omega$ is $8.64\times10^{-9}$ ([scale binding](../results/development/scale-binding.json)). The complete action has not selected $q$, and this small-flux branch does not invoke MMP's parametrically large-$q$ control.

**recursive source binding.** Applying the simplest inherited transparent-LLL state law, $\kappa_c=\Omega\kappa_p$, fixes $t_v=\Omega^2t_u$ and gives parent Killing power $-0.02227623165818482$. On the actual black-universe neck its projected radial null components are

$$
T^{(4)}_{++}=0.007886485863354125,
\qquad T^{(4)}_{--}=0.05267065468742516.
$$

Both have the opposite sign from the required neck source. The imported MMP LLL therefore cannot be substituted for the full charged black-universe CTP tensor ([source binding](../results/development/recursive-source-binding.json)). This fixes the work cursor: the charged angular/compact covariance and its metric variation must be evaluated on the horizon-penetrating domain; another scale scan cannot repair this source mismatch.

**Charged CTP neck result.** The new magnetic angular calculation performs that transport at the locked scale and combines it with the compact Wilsonian local response. In the common child frame its retained tensor gives radial null components $-0.05681746054$ and $-0.06171294116$, so the retained sign gate passes ([charged CTP source](../results/development/charged-ctp-neck-source.json)). The hard full gate remains open because the first positive compact levels still lack their nonlocal horizon-domain CTP covariance and four ADM variations; this is a fail/incomplete result for the requested full owner, and no background projection is started.

**Positive compact CTP completion.** The first two positive compact levels below the locked cutoff are now included with their massive exterior reflection and fourth-order superadiabatic covariance subtraction. The massless limit reproduces the existing $E_2+E_4$ and $P_2+P_4$ subtractions. Reading the preceding angular tensor rather than rerunning it gives the completed free-Gaussian child-frame source

$$
(\rho,T_{01},p_\parallel,p_\perp)
=(0.10074836289,\ 0.00122387016,\ -0.34978809527,\ 0.12958962232),
$$

$$
\boxed{T_{++}=-0.24659199207,\qquad
T_{--}=-0.25148747269.}
$$

Both null signs remain negative under angular, frequency-extent and frequency/time changes, so the hard neck gate passes in the declared free Gaussian realization ([compact CTP completion](../results/development/charged-compact-ctp-completion.json)). Its homogeneous projection has $p_{\rm iso}=-0.03020295021$ and $w_{\rm iso,neck}=-0.29978601484$ with nonzero anisotropic stress. These are local neck data; the following projection determines how they enter the child ledger before any $H(z)$ calculation.

**Proper-volume and child-clock deposition.** The same completed tensor is now projected through the stored Bronnikov child ledger, without another source calculation. At the neck, $V=4\pi a_\parallel r^2=24.21233094747$, the proper-volume rate is $12\pi$, $H_\perp=0$, and the standard anisotropic conservation equation gives

$$
\boxed{\dot\rho_0
=-\frac{\dot V_0}{V_0}\rho_0-H_\parallel p_\parallel
=0.38776013531.}
$$

Volume dilution contributes $-0.15686733379$ and directional pressure work contributes $+0.54462746910$. The proper energy per unit coordinate $z$ is $E_0=2.43935270463$; its child-time rate is $13.18670052493$, entirely accounted for by that pressure work. The stored parent power obeys

$$
P_{\rm parent}=-P_z=-a_\parallel V T_{01}
=-0.057095079694875,
$$

with residual $4.4\times10^{-16}$. The parent Killing direction is spatial in the child, so this quantity is the opposite of a conserved child momentum charge rather than a continuing child energy source. The explicit free compact blocks are decoupled and give regular inter-block $Q=0$; no visible/dark assignment is made. This result fixes the local density Cauchy jet, while later pressure functions require metric backreaction and evolution of the same state ([background projection](../results/development/nsc-background-projection.json)).

**Initial child-backreaction constraint gate.** Before the same covariance can be evolved, its Cauchy data must satisfy the imported ADM constraints. At the unit-radius neck, the locked Einstein coefficient $A=0.04501936182826115$ requires $\rho=2A=0.0900387236565223$ and zero homogeneous momentum density. The completed tensor gives $\rho=0.100748362886786$ and $T_{01}=0.00122387015580509$, hence

$$
\boxed{\mathcal C_H=-0.118944813912720,
\qquad \mathcal C_M=-0.0135927088490713.}
$$

Both null components remain negative, but the metric and source are not on the same constraint surface. The evolution therefore stops at $T=0$ before replaying the CTP covariance. Closing this gate requires the same action to supply $\Delta\rho=-0.0107096392302639$ and $\Delta T_{01}=-0.00122387015580509$, or to derive a different initial geometry; neither is inserted here ([backreaction gate](../results/development/nsc-child-metric-backreaction.json)).

**Constraint-complete geometry from the source.** The second permitted route keeps the completed tensor and every action/state parameter fixed. Its radial block is type-I and uniquely selects the subluminal Landau normal

$$
v=0.00491447570673408,
\qquad T_{01,L}=0.
$$

At a temporal minimum with $H_\perp=0$, the imported Hamiltonian constraint then selects

$$
\boxed{r_\star=\sqrt{\frac{2A}{\rho_L}}
=0.945328394434129.}
$$

The resulting initial metric remains Kantowski--Sachs, retaining $a_{\parallel,0}=1.92675607703328$ and $H_{\parallel,0}=1.55702116929053$ from the seed. It is not the stored Bronnikov or MMP unit-radius profile. The maximum constraint and frame-invariance residual is $2.2\times10^{-16}$, and both null contractions are $-0.249027703022043$. This passes the initial assembly gate and authorizes covariance--metric evolution from the source-selected geometry ([constraint-complete neck](../results/development/nsc-constraint-complete-neck.json)).

**First coupled CTP/metric step.** Re-evaluating the same source on the selected geometry is stronger than carrying over its integrated tensor. The radius change modifies the angular Hamiltonian by

$$
\frac{1}{r_\star}-1=0.0578334533139646,
$$

and the existing metric vertex is

$$
\boxed{
\left.\frac{\partial\rho}{\partial\log r}\right|_C
=-2(\rho_L+p_{\perp,L})
=-0.460687999769215.
}
$$

Thus the unit-radius stress cannot be copied unchanged onto $r_\star$. The source records retain integrated stress rows and covariance-eigenvalue diagnostics, but not the complex matrix $C_j(k,k')$ for every retained mode. Those four moments do not determine the mode evolution

$$
\dot C_j=-i[H_j,C_j]
$$

or its later pressure. The coupled trajectory therefore records a break at time zero during source re-evaluation, before the first metric step. Its next owner is the explicit mode-resolved Cauchy state in the canonical half-density basis, with the existing fourth-order reference carried on the same history ([coupled evolution gate](../results/development/nsc-coupled-ctp-metric-evolution.json)).

**Mode-resolved Cauchy state.** One authorized base extraction now persists the actual physical Gaussian covariance for every retained channel in the canonical half-density basis

$$
\chi=r\sqrt{a_\parallel}\,\psi.
$$

The content-addressed payload contains 33 channels and 1,904 complex $2\times2$ blocks: 48 LLL nodes, 576 massive-angular blocks and 1,280 positive-compact blocks. Identical copies remain multiplicities. A cold, pickle-free reload gives

$$
\lambda_{\min}(C)=-2.67\times10^{-16},
\qquad
\lambda_{\max}(C)=1+4.15\times10^{-14},
$$

and reconstructs the complete old-surface tensor with maximum residual $9.11\times10^{-13}$. The payload separately retains the massive-angular $E_2+E_4/P_2+P_4$ terms, positive-compact fourth-order reference and one-counted local allocation. This passes serialization, CAR and old-surface reconstruction ([mode-resolved state](../results/development/nsc-mode-resolved-cauchy-state.json)).

The local Landau current factor is not the full Cauchy isometry: endpoint rapidity does not specify the tilted hypersurface or intervening metric/gauge/Dirac history, and a direct nonunitary boost can violate $C\le I$. The full map must satisfy

$$
U_{L0}=J_L\mathrm{Res}_{\Sigma_L}E_g
\mathrm{Res}_{\Sigma_0}^{-1}J_0^{-1},
\qquad U_{L0}^\dagger U_{L0}=I.
$$

The finite $r_\star$ stress also requires the existing fourth-order reference and local induced allocation as general Kantowski--Sachs history providers. Those two owners, rather than occupations or integrated moments, now gate metric evolution.

**Conditional Landau Cauchy propagation.** The serialized state is now connected to a blockwise history API. For a declared homogeneous ADM history,

$$
H_{jn}=N\left[-m_j\sigma_1+\frac{\lambda_n}{r}\sigma_2
+\frac{k}{a_\parallel}\sigma_3\right]-\beta k I,
$$

and exact midpoint exponentials give $U[\mathrm{history}]$ and $C_L=UC_0U^\dagger$. Two smooth diagnostic histories share all recorded seed and Landau endpoint data. Both preserve unitarity below $1.2\times10^{-14}$ and retain the covariance spectrum, but

$$
\boxed{
\max\left|U_{0.5}-U_{1.0}\right|=1.99966699255.
}
$$

This finite witness shows that the endpoint tuple does not select the physical Cauchy map. The controls are not metric solutions, and neither is promoted to the history in $E_g$. Consequently the finite $r_\star$ stress, its null signs and its updated constraints remain undefined ([Landau Cauchy gate](../results/development/nsc-landau-cauchy-isometry.json)). A same-action history-selection equation or joint history/state boundary-value solution must choose the intervening geometry before stress can be assigned.

**Joint history/state boundary-value gate.** The two permitted solve routes are now evaluated against the same records. Dynamic route A cannot start because the serialized unit-radius seed has

$$
(\mathcal C_H,\mathcal C_M)
=(-0.118944813912720,-0.0135927088490713).
$$

Two-point route B has unitary conditional maps, but its endpoint-identical controls satisfy

$$
\boxed{
\max|C_A-C_B|=0.192089815976105,
}
$$

so at least two distinct state histories obey every currently executable endpoint, CAR and unitarity condition. The selected result is route C: nonuniqueness bound and stop ([joint BVP record](../results/development/nsc-joint-history-state-bvp.json)). The exact missing selector is

$$
\frac{\delta\Gamma_{\rm one}^{\rm CTP}[g,C[g]]}
     {\delta g_\Delta^A(\tau)}=0,
\qquad A\in\{N,\beta,q_{\rm ADM},r\},
$$

implemented with the general-KS fourth-order reference, local induced history and boundary/interface variation. This is a missing executable part of the declared action rather than a new law. No finite stress or diagnostic duration is selected.

Let Q_b denote internal transfer into the visible component and J_b,J_d denote external room supply. The complete ledger is

$$
\dot\rho_b+3H(\rho_b+p_b)=Q_b+J_b,
\qquad
\dot\rho_d+3H(\rho_d+p_d)=-Q_b+J_d,
$$

$$
\dot\rho+3H(\rho+p)=J_b+J_d.
$$

Internal Q_b cancels. The Bianchi identity requires the full source of ordinary Einstein gravity to be conserved; external room supply therefore needs its reservoir/boundary gravitational accounting. A throat surface integral is power, and a worldtube integral is energy. A cosmological density rate additionally requires proper time, proper volume and deposition into the identified component. The parent-matched Killing power $0.0001422206795$ is therefore not a cosmological $Q$.

For the closed flat constant-G baseline with visible dust, a constant fraction f below unity requires

$$
\frac{Q_b}{H\rho_b}=1-2q_{\mathrm{dec}}.
$$

Zero Q_b alone does not force deceleration: dust plus positive Lambda has zero exchange and approaches de Sitter expansion. A finite regular plateau is an additional condition. A known scalar Q_b also leaves pressure and perturbations undetermined ([plateau conditions](../results/nsc-5-plateau-conditions.json), [source requirements](../results/nsc-7-observable-bridge.json)).

Predicting H(z), BAO and CMB requires the background and thermal history. Predicting clustering and lensing also requires pressure perturbations, anisotropic stress, momentum transfer and constrained metric response. The resulting background and perturbations are compared jointly with Lambda-CDM and the stated observational likelihoods.

The frozen development targets are documented with their model assumptions in [observational targets](../docs/nsc-observational-targets.md). H0DN reports 73.50±0.81 km/s/Mpc, with significance depending on the comparison [23]. DES Y6 reports S8=0.789±0.012, a 2.6-sigma projected difference from its primary-CMB comparison [24]. KiDS-Legacy reports S8=0.815 with +0.016/-0.021 uncertainty, agreeing with Planck at 0.73 sigma [25]. These development targets specify the observational comparisons for the derived background and perturbations.

## 6. Nuclear physics and the remaining common-action equations

The Skyrme/BPS sector is a possible low-energy construction pattern. Its coefficients must be obtained from the evolved common effective action before its particles, binding energies or stellar limits are predictions of NSC.

The often-quoted 3.34 solar-mass maximum is imported from a nuclear-calibrated Einstein–BPS model [5]. The historical runner records the published value rather than independently recomputing it. Its potential and calibrated nuclear coefficients matter; other potentials yield different maxima. It serves here as the stellar control for that imported model ([stellar benchmark](../results/nsc-1-gravitating-bps-observation-link.json)).


The source equation is

$$
\left.\frac{\delta\Gamma_{\rm one}}{\delta g_\Delta^{\mu\nu}}\right|_{g_\Delta=0}=0.
$$

The [common-source derivation](../docs/nsc-common-source-derivation.md) binds
boundary coherence, energy transfer, and metric variation. The canonical CTP
owner supplies causal state response. The recursive zero-tadpole law

$$
\Gamma_{\rm rel}=(1-\mathcal P_0)\Gamma_{\rm one},
\qquad
(V,A,C)\mapsto(0,A,C),
$$

fixes the homogeneous unlinked $V_{\rm full}=0$ and gives
$\lambda_4=\Xi=0$ in the charged-throat map. Einstein, gauge, curvature,
finite Casimir, and link responses remain. The twenty-three
[development records](../results/development-snapshot.json) provide the
boundary state, measure, ADM projections, parent mass balance, MMP embedding,
causal source, relational normalization, scale candidate and direct source
binding. The positive compact completion closes the free-Gaussian neck tensor,
and its proper-volume/clock projection closes the local child ledger. The
initial ADM gate identifies the density and momentum mismatch; the
geometry-from-source follow-up closes it and makes joint evolution the next
owner. Conditional Cauchy maps now exist, and the BVP gate identifies the
general-KS same-action history variation that remains before joint evolution.

The project contribution being developed is the common-action closure across these sectors. The primary ingredients retain the attribution given in their derivations and references.

## 7. Evidence and reproduction

The preprint collection contains 100 records: 58 frozen historical records and
42 scoped follow-ups. A separate twenty-three-record development snapshot, imported
through laboratory commit $234f03b$, carries the common-source, constraint
closure, mode-resolved state, conditional Cauchy maps and joint BVP gate. The 100-record manifest remains byte-preserved;
the development index authenticates its additional source files and comparison
policies separately.

The new records compare their published fields, formulae, scope statements and source/input hashes under the absolute and relative tolerances stored in those JSON files. Earlier horizon-source, warped-source and compact/light policies are unchanged. For the spectral endpoint's small finite-difference residual, portable mode additionally checks the Richardson and subtraction identities, retains the generator's original derivative-accuracy requirement, and propagates the raw operands' comparison budgets. The raw derivatives retain their original tolerances; exact mode is unchanged. Structure, source/comparator hashes and scope fields remain exact. These reproduction policies are distinct from the convergence studies and physical assumptions in each record.

Use the default demonstration to inspect fourteen authenticated stored cases. It launches no scientific generators. Explicit recomputation and the full release gate are documented in [reproducing](../docs/reproducing.md). One integrated gate authenticates the dependency graph, runs relevant tests, regenerates the collection in isolation, and checks deterministic PDF output. Same-environment byte identity and portable all-field agreement are reported separately.

The appendices preserve the normalization corrections, finite controls and historical scale diagnostics supporting this construction.

## Appendix A. Historical scale tests and corrections

This appendix records the historical scale candidates and the normalization corrections that supersede them.

### A.1 Historical scale candidates

The inheritance scale is

$$
\zeta=(\Lambda L_\star)^2,
\qquad
\mu^2=1-\frac{3\pi}{2\zeta}.
$$

Matching the local on-room vacuum-form component $108/1015$ to the doubled gap, then using the exact expanding child’s Ricci sign, selected a *conditional* one-scale gap $|\Phi|^2/\Lambda^2=1006/1015$ and a boundary cross-projection $\Xi_{\mathrm{boundary}}=54/503$ ([`nsc-1-s-one-child-orientation.json`](../results/nsc-1-s-one-child-orientation.json)). The exact child asymptotic rejects that magnitude. The actual child vacuum is $\Lambda_{e,\mathrm{child}}L_\star^2=18\pi$, so a real positive gap requires $\zeta>3\pi/2$. Direct $\Lambda L_\star=1$ gives a negative $\mu^2$ and is not the physical child-matched branch ([`nsc-1-s-one-child-scale-correction.json`](../results/nsc-1-s-one-child-scale-correction.json)). The identities for $\Phi$ as gap, Schur self-energy, and metric link remain. The numbers $1006/1015$ and $54/503$ are **superseded** as child-matched predictions.

The ZETA1 chain then searched for a child-allowed stationary scale of the joined-minus-disconnected Dirac spectrum on the global slice. Every root in the following table is a **numerical diagnostic**. None is a physical $\zeta$. The later anomaly compensation supersedes all of them as stationary points of the one equation.

| Compact record | What was computed | Diagnostic value | Why it is not physical $\zeta$ |
|---|---|---|---|
| [`nsc-2-zeta1-lowest-mode.json`](../results/nsc-2-zeta1-lowest-mode.json) | Isolated lowest-mode heat minimum | $\zeta=5.096657218448003$, $\mu=0.27458356021684416$ | Lowest sector only; heat trace, not the determinant; **superseded candidate**. |
| [`nsc-2-zeta1-angular-tower.json`](../results/nsc-2-zeta1-angular-tower.json) | Angular tower through $n=12$ | Diagnostic $n=1$ heat root $5.096622100098284$; $n=2$ removes that heat-only root; $n\leq 12$ minimum derivative $7.350807275212887$ | Heat-only truncation; **superseded candidate**; absence of a heat root is not a theorem against stationarity. |
| [`nsc-2-zeta1-regulated-determinant.json`](../results/nsc-2-zeta1-regulated-determinant.json) | Exponentially regulated joined-minus-disconnected determinant | $\zeta_{\det}=4.99072409354754$, $\mu_{\det}=0.2361577587050819$ | Uncompensated determinant; warp, lapse, shift, and anomaly omitted; **superseded candidate**. |
| [`nsc-2-zeta1-y-boundary-sensitivity.json`](../results/nsc-2-zeta1-y-boundary-sensitivity.json) | Factorized compact-$y$ controls | Periodic root $4.8765555133228515$; Neumann root $4.79786026064662$; Dirichlet and MIT bag: no root | Factorized unwarped control; a zero mode was not selected in order to force a root; **superseded candidate**. |
| [`nsc-2-zeta1-orbifold-parity.json`](../results/nsc-2-zeta1-orbifold-parity.json) | Joint fermion/metric $S^1/\mathbb{Z}_2$ parity | $\zeta_{\mathrm{orb}}=4.747863009733466$, $\mu_{\mathrm{orb}}=0.08643829078229333$ | Unwarped factorization; **superseded candidate**. |
| [`nsc-2-zeta1-warped-y.json`](../results/nsc-2-zeta1-warped-y.json) | Exact conformal compact warp | $\zeta=4.748389947082489$, $\mu=0.08707307719677547$ | Factorized radial-plus-$y$ sum, not the full Dirac operator with lapse and shift; **superseded candidate**. |

The table records each truncation; the following sections give the corrected scale equation.

### A.2 Anomaly and inherited subtraction

The invariant partition function requires cancellation of pure cutoff rescaling while retaining the physical change of the child-constrained gap. For one eigenvalue the physical child-link derivative after compensation is negative whenever the gap constraint $\mu^2=1-3\pi/(2\zeta)$ is imposed ([`nsc-2-zeta1-anomaly-decomposition.json`](../results/nsc-2-zeta1-anomaly-decomposition.json)).

At the warped determinant candidate $\zeta=4.748389947082489$:

| Piece | Value |
|---|---:|
| Determinant derivative | $0.01201522080638507$ |
| Pure cutoff anomaly | $12.336590431250533$ |
| Physical child-link derivative | $-12.324575210444092$ |

The compensated child-link derivative stays negative on the $257$-point scan $3\pi/2<\zeta\leq 200$, from $-19.422735759304963$ to $-0.5875504421475721$.

Those historical numbers used the inconsistent conversion $(\lambda_j+\mu^2)/\zeta$. They remain numerical diagnostics of that proxy.

**Therefore those determinant-only $\zeta$ values are invalidated as stationary points of the one equation.** $\zeta_{\det}$ and $\zeta_{\mathrm{warped}\,y}$ are uncompensated determinant diagnostics, not inheritance-scale predictions.

**Historical interpretation.** The next missing contribution was identified as the recursive child tail, rather than an independently weighted Einstein–Gauss–Bonnet or heat action ([`nsc-2-zeta1-anomaly-owner-correction.json`](../results/nsc-2-zeta1-anomaly-owner-correction.json)). Geometry is already the compensating anomaly of the same fermionic determinant. Adding that geometry again is forbidden double-counting. The calculations through the scan joined one parent to one child treated as a terminal asymptotic domain. Recursion requires that the child’s boundary response already contain its own child tail. The compact `nonclaims` include `recursive_tail_is_known_to_restore_a_scale_root: false`.

Related routes already rejected on their own JSON nonclaims, and not repaired by retuning $\Phi$ or adding a second geometric weight:

- $\Lambda L_\star=1$ as the physical child-matched gap;
- the lowest-mode heat-only scale minimum;
- independently weighted Einstein–Gauss–Bonnet or heat geometry as the missing scale owner;
- pure critical five-dimensional Einstein–Gauss–Bonnet as the complete interface law;
- the standalone bosonic heat-trace graviton covariance as a positive Stieltjes propagator ($g''(0)=228/175>0$; [`nsc-1-s-one-reflection-positivity-obstruction.json`](../results/nsc-1-s-one-reflection-positivity-obstruction.json));
- the inverse relative heat Hessian as Osterwalder–Schrader positive metric covariance (the full $31$-time reflection matrix has minimum eigenvalue $-1.3679466754324912$; the source record counts fourteen negative modes while portable Linux counts thirteen because one near-zero mode crosses the numerical counting tolerance; [`nsc-1-s-one-relative-reflection-test.json`](../results/nsc-1-s-one-relative-reflection-test.json));
- constant two-wall shadow matter as the dominant kiloparsec gravitating mass.

The fermionic sheet-resolution observable is Osterwalder–Schrader positive in the controlled flat sector ([`nsc-1-s-one-fermionic-relative-observable.json`](../results/nsc-1-s-one-fermionic-relative-observable.json)). That does not make $\sigma_3$ the complete graviton, and it does not restore a scale root. Full-exponential positivity in that sector is a heat-kernel diagnostic, not the complete physical anomaly.

### A.3 Corrected units and the sign theorem

The ZETA1 definitions are $\lambda_j=L_\star^2\mathrm{eig}_j(D_{\mathrm{spatial}}^2)$, $\zeta=\Lambda^2 L_\star^2$, and $\mu^2=|\Phi|^2/\Lambda^2=1-a/\zeta$ with $a=3\pi/2$. The same declared additive-gap operator therefore requires

$$
\boxed{
q_j=\frac{\lambda_j}{\zeta}+\mu^2=1+\frac{\lambda_j-a}{\zeta}.
}
$$

The historical runners implemented $(\lambda_j+\mu^2)/\zeta$, a physical mass squared smaller by $1/\zeta$ than the declared mass squared. Compact record [`nsc-2-zeta1-unit-closure-check.json`](../results/nsc-2-zeta1-unit-closure-check.json) reconstructs the old derivative as a control, then changes only that conversion. It retains radius 30, 750 radial half-intervals, angular sectors 1–12, 32 compact intervals, and the historical factorization and angular weights.

At the old warped candidate $\zeta=4.748389947082489$:

| Quantity | Historical argument | Argument with consistent units |
|---|---:|---:|
| Determinant logarithmic derivative | $0.01201522080638507$ | $-24.319952223294422$ |
| Derivative after inherited cutoff subtraction | $-12.324575210444092$ | $-36.58292832809548$ |

The independent old-argument reconstruction differs from the stored derivative by $1.53\times 10^{-8}$. With radial spacings $0.08$, $0.04$, and $0.02$, the corrected subtracted derivative is $-36.58813190$, $-36.58292833$, and $-36.58162543$. Its sign is resolved.

For each mode the determinant piece is one half of $E_1(q_j)$. At fixed spatial geometry, subtracting the same pure-cutoff term used historically gives

$$
\boxed{
\frac{d\Gamma_{\mathrm{sub},j}}{d\log\zeta}=-\frac{e^{-q_j}}{2q_j}.
}
$$

The disconnected radial matrix is the principal submatrix obtained by deleting the throat node. Hermitian eigenvalue interlacing and the positive decreasing weight $f(\lambda)=e^{-q(\lambda)}/q(\lambda)$ then imply that the total subtracted derivative is strictly negative for every $\zeta>3\pi/2$ in this finite, frozen-geometry, additive-gap family. The $257$-point scan corroborates the sign; it is not the basis for extending the statement between sample points. A diagnostic root of the unsubtracted finite determinant near $\zeta=6.09675392014$ still has subtracted derivative about $-14.15469$, so it does not solve the adopted subtracted scale equation. No physical scale is inferred.

**Consequence.** Refining this family or selecting another determinant-only minimum cannot close its scale equation. A further attempt needs a derived change in the physical operator, its geometry/link dependence on scale, or the anomaly prescription. Merely naming the omitted recursive tail does not prove that it supplies that change.

For the same relative-sheet symbol

$$
D_\delta
=
\begin{pmatrix}
p e^{-\delta} & \Phi \\
\Phi & -p e^{\delta}
\end{pmatrix},
$$

the unregulated finite determinant is independent of $\delta$. The proper-time functional used later in ZETA1 has a strictly positive regulated relative Hessian at $\delta=0$ for $p\neq 0$. Raw determinant invariance therefore cannot justify setting the regulated relative Hessian to zero. The primary anomaly paper distinguishes its normalization scale from the spectral cutoff ([arXiv:1106.3263](https://arxiv.org/abs/1106.3263)).


## Appendix B. Finite operator and boundary controls

### B.1 Finite regulated variations and a spacetime obstruction

For a Hermitian finite Dirac matrix with no unresolved zeros, one explicit proper-time prescription is ([`nsc-3-regulated-recursion.json`](../results/nsc-3-regulated-recursion.json); [`docs/nsc-regulated-recursion.md`](../docs/nsc-regulated-recursion.md))

$$
\Gamma_{\Lambda,M}(D)
=
\frac{1}{2}\mathrm{Tr}E_1(D^2/\Lambda^2)
+N[\log(M/\Lambda)+\gamma_E/2].
$$

Its first and second variations retain noncommuting perturbations through divided differences. A finite scale cocycle $C(t)=\Gamma_{\mathrm{reg}}(D_t)-\Gamma_{\mathrm{reg}}(D)-t\mathrm{Tr}\Sigma$ is derived rather than discarded; its endpoint matches an independent integral of $C'(t)$. Normalization $M$ and cutoff $\Lambda$ remain distinct. This finite orthonormal-basis calculus is not a derived continuum gravitational measure.

An ultrastatic Euclidean product $D_E^2=-\partial_\tau^2+H^2$ produces a frequency factor $\mathrm{Tr}_{\mathrm{sp}}e^{-tH^2}/\sqrt{4\pi t}$. That factor changes the variational functional; a spatial heat trace cannot silently inherit a four-spacetime anomaly.

**Imported result / repository witness.** The recorded horizon-penetrating coordinate-time Hamiltonian is not elliptic in the trapped region: at the throat $\beta^2=3\pi/2>1$, and a nonzero spatial covector makes one principal eigenvalue exactly zero. Elliptic spatial heat calculus therefore cannot be imported as the full spacetime determinant. The general loss of ellipticity at horizons is documented by Finster and Röken ([arXiv:1512.00761](https://arxiv.org/abs/1512.00761)); the project contribution is the explicit witness for this carrier, not discovery of the phenomenon.

Along $D(t)=e^{-t}D$ at fixed cutoff, normalization, and finite domain, $d\Gamma/dt=\mathrm{Tr}e^{-D(t)^2/\Lambda^2}>0$. No isolated common-scale extremum exists in that finite family. Normalized finite recursion with $1/\Omega$ agrees with direct chain inversion in declared-matrix controls. $H$ and $b$ there are not yet computed throat maps. The compact `nonclaims` include `full_covariant_action_derived: false` and `physical_stationarity_solved: false`.

### B.2 The isolated radial throat is gapless

The unwarped spatial radial Dirac operator already declared by the project,

$$
D=-i\sigma_2\partial_\rho+\sigma_1\frac{\kappa}{\sqrt{1+\rho^2}},
$$

has essential spectrum equal to $\mathbb{R}$ ([`nsc-3-radial-spectrum.json`](../results/nsc-3-radial-spectrum.json); [`docs/nsc-radial-spectrum.md`](../docs/nsc-radial-spectrum.md)). A constructive Weyl sequence of compactly supported $H^1$ bumps escaping to infinity shows that every real energy belongs to the essential spectrum. Smooth radial geometry and current-conserving gluing alone do not generate a nonzero spectral mass gap. Finite-box levels are not physical masses. There is no global $L^2$ zero mode. The result does not cover the compact warped, interacting, recursive, or Lorentzian operator, and it does not reject Nested-Space Cosmology.

### B.3 A controlled geometric gap from repeated throat geometry

Repeating the finite segment $\rho\in[-R,R]$ with $w=(1+\rho^2)^{-1/2}$, instead of two asymptotically widening ends, produces an open spectral gap without an inserted constant mass ([`nsc-3-geometric-chain.json`](../results/nsc-3-geometric-chain.json); [`docs/nsc-geometric-chain.md`](../docs/nsc-geometric-chain.md)). At zero energy the exact motif transfer multipliers are $\exp(\pm 2\mathrm{asinh} R)$; neither has modulus one, so zero is absent from every real Bloch fibre. Continuum first band edges and second-order lattice agreement are

| $R$ | Continuum first band edge | Finest-lattice absolute error |
|---|---:|---:|
| 2 | $0.7034881641$ | $6.88\times 10^{-6}$ |
| 4 | $0.4552655377$ | $6.20\times 10^{-6}$ |
| 8 | $0.2479700183$ | $2.43\times 10^{-6}$ |

<!-- nsc-figure:geometric-gap -->

The next-room link is the staggered stencil entry $B_{\mathrm{last\,edge},\mathrm{first\,node}}=1/h+w_{\mathrm{last}}/2$, not a fitted mass. Normalized finite recursion with that geometric $B$ agrees with direct multi-motif inversion. $R$ and $\Omega$ remain declared inputs. Seams are periodic in $w$ but gravitational junction equations are unproved. This is a positive project-level controlled realization of geometric confinement, not a new-to-world Floquet discovery and not a Lorentzian nested cosmology or an electron.

### B.4 Computed finite boundary maps and the threshold distinction

The finite first-order radial calculation now evaluates parent and child Weyl maps with sparse LU, throat reconstruction, independent ODE integration, and direct Schur controls ([nsc-3-boundary-response.json](../results/nsc-3-boundary-response.json)). The parent lies at positive radial coordinate and its outward throat normal is negative; the child lies at negative radial coordinate and its outward throat normal is positive. The oriented relation is

$$
N_p(E)+N_c(E)=E[m_c(E)-m_p(E)].
$$

A nonzero real-energy probability-current probe tests the local conservation law. The compact warp calculation reports mode-mixing leakage: its average-warp expression is a projected approximation. These results compute finite spatial maps; the full horizon boundary maps, physical compensator, and Lorentzian geometry equations remain open.

The same isolated radial geometry supplies an explicit threshold counterexample. Solving the upper squared-Dirac equation with unit throat value and vanishing exterior value gives, for t=asinh(R),

$$
I_p=e^{3t}/6+e^t/2-2/3,\qquad I_c=2/3-e^{-t}/2-e^{-3t}/6,
$$

$$
N_p+N_c=I_p^{-1}+I_c^{-1}\longrightarrow 3/2.
$$

Independent ODE solutions reproduce these maps ([nsc-3-threshold-response.json](../results/nsc-3-threshold-response.json)). The positive boundary stiffness coexists with the gapless radial bulk spectrum; it is not a particle mass of 3/2. A mass identification requires the full pole problem, residues and physical state. At zero energy m(E) can be singular, so E m(E) must be evaluated by its limit or by the squared-operator boundary problem.

## Appendix C. Compact, UV and charged-throat controls

This appendix records technical controls behind the compact-sector and source matching. It does not add a compact-mass result record, recompute the imported wormhole, or select a physical coupling.

The compact-mass application uses four-dimensional signature +--- and a fifth gamma matrix iγ⁵, giving five-dimensional signature +----. It applies Fischmann's conformal law and the interval reduction of Grossman and Neubert. The zero mode is a single Weyl field; it cannot be passed to a massive four-component source as though it contained both chiralities. Opposite endpoint chirality supplies the mirror zero mode and the same massive tower. Self-adjointness alone does not choose between them.

The compact-interaction overlaps are dimensionless profile integrals after the bulk Einstein normalization

$$
\kappa_{4,\mathrm{bulk}}^2=\frac{\kappa_5^2}{I_3},\quad I_3=\int e^{3\sigma}dY.
$$

For the declared warp, I₃/L⋆=1.96509101243726. The 3111 left-current fermionic component remains of order 6 after antisymmetrization; the 0111 vertex vanishes by the even bulk warp. These signs are not energy-density signs.

For one bulk Dirac field, the ultraviolet map retains

$$
\frac{4(\Lambda^5-\nu^5)}{5(4\pi)^{5/2}}\int\sqrt{g}.
$$

Omitting that volume term in an unconstrained metric variation would change the source problem. Finite coefficients $c_R$ and $c_T$ are left unspecified; leading cutoff compensation leaves $c_T-9c_R$ unchanged.

For the imported throat, the short-exterior-length relations of Maldacena, Milekhin and Popov include

$$
r_e^2=\frac{\pi q^2G_N}{g^2},\quad E_{\mathrm{quantum}}=-\frac{q}{8\ell}.
$$

The approximation uses large flux, weak gauge coupling and specified exterior geometry. The optical/redshift scale ℓ is not the proper throat length. Exterior matching and mouth stabilization must satisfy the source's assumptions; these lengths are not substituted as NSC predictions.

### C.1 Vacuum, gravity and gauge coefficients must match together

<!-- nsc-claim:vacuum-charge-coefficients -->
For the retained local action $A R-CF^2-V$, with positive $A,C$, define $G_N=(16\pi A)^{-1}$, $g^2=(4C)^{-1}$ and $\lambda_4=V/(2A)$. The magnetic charge-radius convention of [32] gives

$$
\Xi=\lambda_4r_Q^2=\frac{q_{\mathrm{mag}}^2VC}{8A^2}.
$$

The integer magnetic flux $q_{\mathrm{mag}}$ is distinct from the radial-metric factor and cosmological deceleration parameter. The common five-dimensional positive spectral coefficients obey $\Xi_{\mathrm{vol}}\ge q_{\mathrm{mag}}^2$; the nonnegative proper-time class gives the stronger $\Xi_{\mathrm{vol}}\ge3q_{\mathrm{mag}}^2$. This follows from their Mellin and warp-moment inequalities. In contrast, the standard positive-vacuum charged extremal seed requires $\Xi\leq1/4$ [42]. The controlled Gaussian warp and compact Einstein boundary correction do not bridge that mismatch ([coefficient condition](../results/development/vacuum-charge-matching.json)).

This excludes using those retained coefficients alone as the near-extremal seed of [32]. It does not exclude the complete NSC functional or all charged geometries. A vacuum subtraction, another positive profile or an additional copy count is not selected to force that seed.

### C.2 Compact source and holonomy scope

For the declared compact chiral projector, the Euclidean adjoint has complementary boundary values. The resulting mixed condition has $S=-K P/2$, with the inward derivative and outward mean-curvature convention stated in the [boundary note](../docs/nsc-compact-boundary-action.md). The leading boundary term supplies the required Einstein/GHY ratio; the actual Gaussian endpoints still require a source balance.

The curved compact interaction is a finite piece of the determinant. It does not include the decompactified bulk term, formally $-\ell\mathrm{Tr}|D_4|/2$, or the conformal boundary cocycle. Its flat-space volume coefficient and induced Einstein coefficient are fixed by the known interval and heat formulas [29,38]. Subtracting those two terms from the full interaction is exact bookkeeping, not a short-curvature approximation at arbitrary $\ell^2|R|$.

For the holonomy result, an ultrastatic proper-radial coordinate gives $L_\pm=-\partial_s^2+W^2\pm W'$, $W=\kappa/r>0$. The possible zero solutions are not periodic. With $t=\omega^2+m_n^2$, the discriminant satisfies $\Delta(t)>2$. The positive determinant has phase dependence $\Delta(t)-2\cos(2\pi\alpha)$ [39,40]. Each fermionic $-c\log[\Delta-2\cos(2\pi\alpha)]$, $c>0$, decreases on $0<\alpha<1/2$ and has positive curvature at the AP minimum. This fixes the conditional one-loop real-potential saddle. It neither selects a temperature nor proves that the physical recursive domain contains this loop.

## Appendix D. The source matching equations

### D.1 State and horizon conditions

In the two-dimensional conformal frame $ds^2=A\,du\,dv$, the free sector obeys

$$
T_{uu}=\frac{c_{\mathrm{eff}}(2AA''-A'^2)}{192\pi}+t_u,
\qquad
T_{vv}=\frac{c_{\mathrm{eff}}(2AA''-A'^2)}{192\pi}+t_v.
$$

The mixed component is $T_{uv}=c_{\mathrm{eff}}AA''/(96\pi)$. Here $c_{\mathrm{eff}}=|q_{\mathrm{mag}}|$ is the free field count. Unruh data use $t_u=c_{\mathrm{eff}}\kappa_h^2/(48\pi)$, $t_v=0$; the Hartle–Hawking control has equal incoming and outgoing terms. The parent Killing current is $t_u-t_v$. Its spherical projection divides by the proper area $4\pi r^2$, and the recorded PG transformation controls the horizon limit. This fixes the source within the specified sector and state; it does not construct a global NSC vacuum [44]. The full [source derivation](../docs/nsc-horizon-source.md) also gives the reduced metric constraint and the original-neck comparison.

### D.2 Keep the complex-spinor and measure contributions

For the compact pair, choose $D_0=\lambda\sigma_1-i\sigma_3\partial_Y$. A real matrix representation uses the phase basis $\chi=(f,ig)$, for which

$$
|D_0\chi|^2=(\lambda g-f')^2+(\lambda f-g')^2.
$$

The cross term integrates with the nonconstant weight $e^{-s\sigma}$; it cannot be removed by restricting both physical components to real values. The direct complex-field norm, matrix quadratic form and scalar quadrature agree in the [warped-source check](../docs/nsc-warped-source.md). The independent conformal variation retains both adjoint domains before using the paired-copy symmetry. Omitted-sector heat bounds are distinct from the observed convergence of retained Galerkin eigenvalues.

### D.3 Light-field normalization and the remaining functional

Define the compact norm and proper length by

$$
J=\frac{1}{\ell}\int e^{s\sigma(Y)}\,dY,
\qquad \ell_{\mathrm{proper}}=\ell J.
$$

The light singular branch has $\epsilon_0(y)=y/J^2+O(y^2)$, and

$$
H_\nu(0)=\log\frac{\Lambda^2J^2}{\nu^2}
 +2\sum_{p\geq1}E_1\!\left[\left(\frac{p\pi}{\ell_{\mathrm{proper}}\Lambda}\right)^2\right].
$$

The canonical light mass remains zero. Define the two moments

$$
Q_1=\int_{0}^{\infty}H(y)\,dy,
\qquad Q_2=\int_{0}^{\infty}yH(y)\,dy.
$$

The established trace expansion [29,45] gives

$$
V_D=\frac{2Q_2}{(4\pi)^2},\qquad
A_D=\frac{Q_1}{6(4\pi)^2},\qquad
C_{F,D}=\frac{H_\nu(0)}{3(4\pi)^2}.
$$

For the old static energy basis, $R_E=-R_L$; squared invariants and the static box-R term agree. The [matching record](../docs/nsc-compact-matching.md) gives that conversion and the finite changes at three matching cutoffs. The normalization mass of a determinant, the matching cutoff nu, and the physical Lambda are different quantities.

The unresolved causal assembly is now written with the computed reference and state pieces explicit:

$$
\Gamma_{\mathrm{one}}^{\mathrm{CTP}}
=\Gamma_{\mathrm{ref,ren}}^{\mathrm{CTP}}
 +\left(\Gamma_{H,\nu}^{\mathrm{CTP}}+\mathcal C_{\nu,\mu}^{\mathrm{CTP}}\right)
 +\Delta\Gamma_{\mathrm{state}}^{\mathrm{CTP}}
 +\Gamma_{\mathrm{remaining}}^{\mathrm{CTP}}.
$$

The last term denotes the contribution still required by the common measure, compensator, remaining interactions and boundary/recursive sectors. The free massive compact determinant is already included in $\Gamma_H$ and is counted once. It is not an independently adjustable gravitational action. Adding the canonical parent-matched stress to an unconverted cutoff complement would double-count terms. The Euclidean matching identity does not determine the initial state, cross-correlations, heavy-sector causal continuation or transmitting phase. Those are the remaining inputs to the same-action metric equations, not invitations to retune the computed coefficients or to fit $G$.

## Appendix E. Interaction, spectral endpoint and state conversion

This appendix records the interaction, spectral and state controls behind Sections 4.7–4.9, using the cited formulations and the same coefficient account.

### E.1 Projected gauge transmission

The charged collective mode is the canonical two-dimensional current/electric-displacement field. There is no four-dimensional minimally coupled $r''/r$ term in its exterior potential. Angular pressure follows by varying $r$ in the reduced action and must not be dropped. For this real stationary potential, scattering reciprocity equates horizon emission with the computed transmission. Independent scalar-field initial-value checks agree at three frequencies; relative current defects are below $6\times10^{-12}$ in the three base coupling calculations. Boonserm–Visser tail bounds [55] enclose the $g_4=0.5$ charged-power ratio between about $0.9640376$ and $0.9641120$. Higher angular photons, Landau levels, compact excitations and magnetic catalysis remain outside the calculation.

### E.2 Spherical constraints and the EFT contact

In the first-order variables of [46], with $\mathcal V=U X_+ X_-+V+W_m$, the three smeared classical brackets close after the area-dependent charge potential is included in the structure functions. The local temporal-gauge ghost operator has formal determinant $(\det\partial_0)^2\det(\partial_0+X_+ U)$. Geometry/ghost cancellation in that patch does not remove the matter measure. At a positive-radius minimum the map $r\to X$ remains invertible even though $X'=0$; taking $X$ itself as a spacetime coordinate fails where its gradient vanishes.

The curvature correction is first order in a bookkeeping parameter $\varepsilon$. Its constraint increments are $\Delta G_1=0$, $\Delta G_2=e^+_1 F$, $\Delta G_3=-e^-_1 F$ with $F=d(X)K_0^2+F_0(X,\chi)$ and $d(X)=-c_W/(2AX)$. Closure is not asserted for a resummed $G_0+\varepsilon G_1$ at arbitrary $\varepsilon$. Extra matter or neutral CFT sectors belong in the total $T_{\mu\nu}$ before it is squared.

### E.3 Full spectral versus local radius derivatives

The two-copy compact weight on $S^4$ is

$$
\Gamma_5(a)=\frac{2}{3}\sum_{m=2}^{\infty}(m^3-m)\,h_\Lambda(m^2/a^2).
$$

There is no additional copy factor. The comparison retains the canonical light field exactly and approximates only the matched complement through $a_4$. On $S^4$ the Weyl tensor and $\Box R$ vanish; the integrated Euler term is constant.

| Sphere radius $a$ | Full $\partial_{\ln a}\Gamma_5$ | Local $\partial_{\ln a}\Gamma_5$ |
|---|---:|---:|
| $0.325735$, matching the stored child $H$ | $0.001239600$ | $-0.164273768$ |
| $0.43742254$, local stationary control | $0.084326048$ | approximately $0$ |
| $1$ | $14.8033344$ | $14.8393135$ |
| $2$ | $284.848997$ | $284.854519$ |
| $4$ | $4754.593055$ | $4754.594317$ |

The isotropic action response at the child-radius control is $\Pi_a\approx 0.001045912$. It already includes induced geometric terms and must not be added as a second matter stress.

### E.4 Reference conversion and angular subtraction

At the stored child curvature radius $a=0.325735$ and $\nu=1$, the invariant density projection splits as renormalized massless reference $0.103125000$, original finite-cutoff light contribution $2.8628\times10^{-16}$, converted complement $-0.102079088$, and unchanged total Dirac modulus $0.00104591235$.

The finite-radius angular calculation subtracts two-dimensional adiabatic orders $0$, $2$ and $4$ before the angular sum [54], then restores the four-dimensional proper-time finite part with harmonic coefficient $L_R=\ln(\mu R)+\gamma_E/2$. Smooth windows act on the already convergent remainder; they are not a physical regulator or a change of state. The auxiliary reference at five areal radii is recorded in the angular-stress JSON; only the neck value enters the source budget above.

### E.5 Parent-matched covariance

In the exterior current basis, $H_{\mathrm{ext}}=-i\sigma_3\partial_x+V_\kappa\sigma_1$ with $V_\kappa=\kappa\sqrt{A}/r$. For positive frequency, $f=(1+e^{2\pi\omega/\kappa_h})^{-1}$ and $s=\sqrt{f(1-f)}$. Compressing occupied affine-horizon data with the scattering matrix yields the interior covariance

$$
C_{\mathrm{in}}=\begin{pmatrix}1-f& i s R^*\\ -i s R& f|R|^2\end{pmatrix}.
$$

Its eigenvalues are $0$ and $1-f|T|^2$. The charge-related momentum channel is $C(+\omega)=I-\sigma_3 C(-\omega)^*\sigma_3$. The two channels give equal diagonal stress and a nonzero momentum flux. A development version that omitted the factor $i$ was corrected before any accepted result. Global extension, spin/current trace maps and microlocal regularity remain explicit conditions of physical use; Gérard–Häfner–Wrochna [50] is a characteristic-state prescription on its stated domain, not an automatic theorem on the NSC black-universe extension.

### E.6 Flat thermal conversion table

The control domain is flat $\mathbb R^3$ times an antiperiodic Euclidean thermal circle. The cutoff $\nu$ is the light-sector matching cutoff, distinct from the carrier $\Lambda$ and from temperature.

| $T/\nu$ | Canonical thermal $\rho/\nu^4$ | Finite-endpoint thermal $\rho/\nu^4$ |
|---:|---:|---:|
| $0.10$ | $0.0001151454$ | $0.0001151454$ |
| $0.25$ | $0.0044978666$ | $0.0031350475$ |
| $0.50$ | $0.0719658654$ | $-0.0012266736$ |
| $1.00$ | $1.1514538468$ | $-0.0063244282$ |

At $T/\nu=0.5$ the record's finite-endpoint density is $-0.00122667356\,\nu^4$. The discrepancy is exponentially small at very small $T/\nu$, but it is not an exact identity. A flat control at the stored parent Hawking temperature also retains the logarithm of its tiny mismatch; that number is not an estimate for the curved child source.

## References

The references identify the established formulations used in each derivation.

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
17. A. S. Bolton, S. Rappaport, and S. Burles, “Constraint on the post-Newtonian parameter $\gamma$ on galactic size scales,” *Phys. Rev. D* **74**, 061501 (2006), [arXiv:astro-ph/0607657](https://arxiv.org/abs/astro-ph/0607657).
18. Z.-Y. Fan, B. Chen, and H. Lü, “Criticality in Einstein-Gauss-Bonnet gravity: gravity without graviton,” *Eur. Phys. J. C* **76**, 512 (2016), [arXiv:1606.02728](https://arxiv.org/abs/1606.02728).
19. Y. Lee, “Burghelea–Friedlander–Kappeler’s gluing formula and the adiabatic limit of the spectral $\zeta$-determinant of a Dirac Laplacian,” [arXiv:math/0304347](https://arxiv.org/abs/math/0304347).
20. J. Nemec, D. Tománek, and G. Cuniberti, surface Green-function recursion, Appendix A.3 of [arXiv:0711.1088](https://arxiv.org/abs/0711.1088).
21. F. Finster and C. Röken, “Self-adjointness of the Dirac Hamiltonian for a class of non-uniformly elliptic boundary value problems,” [arXiv:1512.00761](https://arxiv.org/abs/1512.00761).
22. G. Teschl, *Ordinary Differential Equations and Dynamical Systems*, Graduate Studies in Mathematics **140**, American Mathematical Society (2012), [author PDF](https://www.mat.univie.ac.at/~gerald/ftp/book-ode/ode).

Repository records cited in the text live under [`results/`](../results/manifest.json). The labeled theory is [`THEORY.md`](../THEORY.md). The current compact checkpoint is [`docs/current-result.md`](../docs/current-result.md).


23. H0DN Collaboration, “The Local Distance Network,” A&A 708 A166 (2026), [arXiv:2510.23823](https://arxiv.org/abs/2510.23823).
24. DES Collaboration, “Dark Energy Survey Year 6 Results: Cosmological Constraints from Galaxy Clustering and Weak Lensing,” [arXiv:2601.14559](https://arxiv.org/abs/2601.14559).
25. A. H. Wright et al., “KiDS-Legacy: Cosmological constraints from cosmic shear,” A&A 703 A158 (2025), [arXiv:2503.19441](https://arxiv.org/abs/2503.19441).
26. C. Dappiaggi, T.-P. Hack and N. Pinamonti, covariant Dirac stress and trace anomaly, [arXiv:0904.0612](https://arxiv.org/abs/0904.0612).
27. A.-P. Jauho, N. Wingreen and Y. Meir, state-dependent quantum transport, [arXiv:cond-mat/9404027](https://arxiv.org/abs/cond-mat/9404027).
28. D. Chung et al., gravitational Dirac production in inflation, [arXiv:1109.2524](https://arxiv.org/abs/1109.2524). This differs from the anisotropic radius-pulse benchmark here.
29. D. V. Vassilevich, “Heat kernel expansion: user’s manual,” *Phys. Rept.* **388**, 279 (2003), [arXiv:hep-th/0306138](https://arxiv.org/abs/hep-th/0306138). The Dirac $a_0,a_2,a_4$ coefficients used here are this established formula, independently reconstructed on the present operator.
30. I. Klich, “Full Counting Statistics: An elementary derivation of Levitov’s formula,” [arXiv:cond-mat/0209642](https://arxiv.org/abs/cond-mat/0209642). Source of the Fock-to-determinant identity.
31. R. Martín and E. Verdaguer, “Stochastic semiclassical gravity,” *Phys. Rev. D* **60**, 084008 (1999), [arXiv:gr-qc/9904021](https://arxiv.org/abs/gr-qc/9904021). The influence-functional approach to gravitational response is imported; the present record connects it to this Dirac geometry.
32. J. Maldacena, A. Milekhin, and F. Popov, “Traversable wormholes in four dimensions,” [arXiv:1807.04726v3](https://arxiv.org/abs/1807.04726). Source of the semiclassical Einstein–Maxwell throat supported by charged massless Dirac vacuum energy.
33. M. Fischmann, “On conformal powers of the Dirac operator on spin manifolds,” [arXiv:1311.4182](https://arxiv.org/abs/1311.4182). Used for the conformal Dirac identification of the free compact carrier.
34. Y. Grossman and M. Neubert, “Neutrino masses and mixings in non-factorizable geometry,” *Phys. Lett. B* **474**, 361 (2000), [arXiv:hep-ph/9912408](https://arxiv.org/abs/hep-ph/9912408). Used for the five-dimensional spinor reduction and chiral endpoint choice; their Randall–Sundrum phenomenology is not imported.
35. O. Castillo-Felisola, C. Corral, S. Kovalenko, and I. Schmidt, “Torsion in extra dimensions and one-loop observables,” [arXiv:1405.0397](https://arxiv.org/abs/1405.0397). Equations 13 and 18–20 supply the higher-dimensional Einstein–Cartan contact reused on NSC compact profiles.
36. F. Pfäffle and C. A. Stephan, “On gravity, torsion and the spectral action principle,” [arXiv:1101.1424](https://arxiv.org/abs/1101.1424). Proposition 5.4 supplies the dimension-independent Bochner/$a_2$ formulae; the four-dimensional curved $a_4$ is not reused as a five-dimensional formula.
37. C. Pagani and R. Percacci, “Quantum gravity with torsion and non-metricity,” *Class. Quantum Grav.* **32**, 195019 (2015), [arXiv:1506.02882](https://arxiv.org/abs/1506.02882). The Palatini flow is a compatibility benchmark, not an NSC beta function.
38. A. Flachi, I. G. Moss and D. J. Toms, “Fermion vacuum energies in brane world models,” *Phys. Lett. B* **518**, 153–156 (2001), [arXiv:hep-th/0103138v2](https://arxiv.org/html/hep-th/0103138v2). The massless interval determinant, conformal reduction and Wilson-line mechanism are reused.
39. G. V. Dunne, “Functional Determinants in Quantum Field Theory,” [arXiv:0711.1178](https://arxiv.org/html/0711.1178v1), equations18–20 and the periodic examples.
40. K. Kirsten and A. J. McKane, “Functional determinants for general Sturm-Liouville problems,” [arXiv:math-ph/0403050](https://arxiv.org/html/math-ph/0403050v1), equation49 and the general boundary-condition construction.
41. C. Krishnan and A. Raju, “A Neumann Boundary Term for Gravity,” [arXiv:1605.01603](https://arxiv.org/html/1605.01603v1), equations4–5 for the Dirichlet metric variation reused here.
42. M. Montero, T. Van Riet and G. Venken, “Festina Lente: EFT Constraints from Charged Black Hole Evaporation in de Sitter,” [arXiv:1910.01648v4](https://arxiv.org/html/1910.01648v4), section2 for the charged de Sitter geometry.
43. V. Borokhov, A. Kapustin and X. Wu, “Topological Disorder Operators in Three-Dimensional Conformal Field Theory,” [arXiv:hep-th/0206054v2](https://arxiv.org/html/hep-th/0206054v2), section4.1 and appendix for the magnetic angular spectrum and degeneracies.
44. S. Iso, H. Umetsu and F. Wilczek, “Hawking radiation from charged black holes via gauge and gravitational anomalies,” [arXiv:hep-th/0602146v2](https://arxiv.org/html/hep-th/0602146v2), equation28 for the free fermionic flux normalization.
45. A. Codello, R. Percacci and C. Rahmede, “Investigating the Ultraviolet Properties of Gravity with a Wilsonian Renormalization Group Equation,” [arXiv:0805.2909v5](https://arxiv.org/pdf/0805.2909v5), Appendix A, equations A10 and A14–A15 for the trace/Mellin identities used in source matching.
46. D. Grumiller, W. Kummer and D. V. Vassilevich, “Dilaton gravity in two dimensions,” *Phys. Rept.* **369**, 327 (2002), [arXiv:hep-th/0204253](https://arxiv.org/abs/hep-th/0204253). Spherical reduction and first-order canonical gravity are reused; the present checks concern this source's normalization and area-dependent charge coupling.
47. L. Parker and D. J. Simon, “Einstein equation with quantum corrections reduced to second order,” *Phys. Rev. D* **47**, 1339 (1993), [arXiv:gr-qc/9211002](https://arxiv.org/abs/gr-qc/9211002). Used for perturbative order reduction of local curvature-squared terms.
48. M. Ruhdorfer, J. Serra and A. Weiler, “Effective Field Theory of Gravity to All Orders,” *JHEP* **05**, 083 (2020), [arXiv:1908.08050](https://arxiv.org/abs/1908.08050). Source for removing redundant curvature operators by metric redefinition.
49. R. Camporesi and A. Higuchi, “On the eigenfunctions of the Dirac operator on spheres and real hyperbolic spaces,” *J. Geom. Phys.* **20**, 1 (1996), [arXiv:gr-qc/9505009](https://arxiv.org/abs/gr-qc/9505009). The four-sphere Dirac spectrum is imported and composed with the project's compact weight.
50. C. Gérard, D. Häfner and M. Wrochna, “The Unruh state for massless fermions on Kerr spacetime and its Hadamard property,” [arXiv:2008.10995](https://arxiv.org/abs/2008.10995). Used as a characteristic-state prescription on its stated domain, not as an automatic theorem on the NSC black-universe extension.
51. Y. V. Gusev and A. I. Zelnikov, “Finite temperature nonlocal effective action for quantum fields in curved space,” *Phys. Rev. D* **59**, 024002 (1999), [arXiv:hep-th/9807038](https://arxiv.org/abs/hep-th/9807038). Fermionic thermal images and the canonical free energy are reused on the adopted finite proper-time endpoint.
52. S. Hollands, “Adiabatic Hadamard States for Dirac Quantum Fields on Curved Space,” [arXiv:gr-qc/9901069](https://arxiv.org/abs/gr-qc/9901069). Used for Hadamard propagation of a Dirac reference through a smooth auxiliary metric.
53. A. O. Barvinsky and W. Wachowski, “Notes on conformal anomaly, nonlocal effective action and the metamorphosis of the running scale,” *Phys. Rev. D* **108**, 045014 (2023), [arXiv:2306.03780](https://arxiv.org/abs/2306.03780). Used for the conformal-stress relation beyond conformally flat metrics.
54. J. F. Barbero G., A. Ferreiro, J. Navarro-Salas and E. J. S. Villaseñor, “Adiabatic expansions for Dirac fields, renormalization, and anomalies,” *Phys. Rev. D* **98**, 025016 (2018), [arXiv:1805.05107](https://arxiv.org/abs/1805.05107). The Dirac adiabatic expansion reused here is a unitary spinor expansion, not a scalar WKB ansatz.
55. P. Boonserm and M. Visser, “Analytic bounds on transmission probabilities,” *Ann. Phys.* **325**, 1328 (2010), [arXiv:0901.0944](https://arxiv.org/abs/0901.0944). Applied to omitted potential tails of the collective-mode scattering problem.
56. M. Alimohammadi and H. Mohseni Sadjadi, “Massive Schwinger model and its confining aspects on curved space-time,” *Phys. Rev. D* **63**, 105018 (2001), [arXiv:hep-th/0011232](https://arxiv.org/abs/hep-th/0011232). Curved-space bosonization is imported; the new work maps that action to the NSC exterior.
