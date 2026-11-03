# Paper guide

**Foundation Edition `v0.27.0` (Published on Zenodo: [`10.5281/zenodo.23890258`](https://zenodo.org/records/23890258))** — [Markdown Source (`paper/nested-space-cosmology.md`)](../paper/nested-space-cosmology.md) · [Illustrated PDF (`paper/nested-space-cosmology.pdf`)](../paper/nested-space-cosmology.pdf):
Establishes **Axioms A0–A4** as the fundamental generating geometric law of space, time, and boundary exchange, with 3+1D field equations and general relativity derived as coarse-grained local dictionary projections (Propositions 1.1–1.2b). Completes the parameter-free deductive chain across **Theorems 1–7** and **Propositions 9.2b & 9.7b**:
- **Propositions 1.1–1.2b:** Curve-bundle reduction to $E = Q_{\mathrm{stretch}} + PV$, **A1** metric generation $g_{\mu\nu}$, and smooth orthonormal bounce invariants $R_{\mathrm{Ricci}}(0) = 6\omega^2$, $K(0) = 12\omega^4$, $\rho_{\mathrm{core}}(0) = 3\omega^2/(8\pi G)$.
- **Theorem 1:** 5-variable least-resistance cycle $\sqrt{\theta^\star} = (\sqrt{16e^{2\chi}+17}-3)/(4e^\chi)$, exact cubic budget identity $9DVB = 4(V+B/2)^3$, and open-valve equation of state $w_{\mathrm{DE}} = -1$ ($dE_{\mathrm{DE}}^\star = (-\frac{1}{3} + \frac{4}{3})\rho_{\mathrm{DE}} dV_- = +\rho_{\mathrm{DE}} dV_-$) vs. closed-box illusion.
- **Theorem 2 & Proposition 9.2b:** Universal exponential-wall bypasses (`1.5` helix plateau, `4.0` child-sphere plateau), `(1,1,1)` ($N=3$) 3D volumetric ground state with short-range nuclear hard-core barrier $\Delta E_{\mathrm{core}} = 8(\sqrt{2}-1) \approx +3.3137$, and finite-facet sphere resolution between the Absolute-Zero floor ($\lambda_{\max} = 2R$) and Planck ceiling ($\ell_P$), fixing $\sqrt{R_\star\ell_P} = \lambda_n$, $F_\star = (\lambda_n/\ell_P)^2 = 1.688\times 10^{38}$ facets, base mass $M_{\star,0} = M_P^3/m_n^2 = 1.8481\ M_\odot$, neutron-star radius $R_{\mathrm{NS}} = 4\lambda_n^2/\ell_P = 10.916\ \mathrm{km}$ ($R_S/R_{\mathrm{NS}} = 1/2 = \theta^\star(\chi_{\mathrm{spaghetti}})$), and TOV black-hole collapse threshold $M_{\mathrm{TOV}} \in [2.0414, 2.1454]\ M_\odot$.
- **Theorems 3–6:** Gauss-divergence $3:1$ exterior-to-cavity Virial theorem ($75$% exterior flow, $25$% cavity, added mass $M_{\mathrm{eff}} = \frac{3}{2}M_{\mathrm{cav}}$), horizon-area conservation in same-sign bubble fission, Cardano quadrupole root $p_+ = \phi^{2/3} + \phi^{-2/3} = 2.1038034027$ ($\phi = \frac{3+\sqrt{5}}{2}$), and Laplace harmonicity $\partial_{\zeta_x}v + \partial_{\zeta_y}v + \partial_{\zeta_z}v = 0$ forcing the $1\to 2\to 4\to 8$ orthogonal cascade.
- **Theorem 7, Lemma 7.1 & Proposition 9.7b:** Cascade freeze-out at $\chi_{\mathrm{freeze}} = 3\chi_1 = 3$ ($\chi_{\mathrm{freeze,dyn}} = 3\ln(1.2+2^{2/3}) = 3.075371$), yielding $\Omega_{\mathrm{DE}} = 70.88\text{--}71.14$%, $\Omega_{\mathrm{DM}} = 23.81\text{--}23.88$%, primordial baryon baseline $\Omega_b = e^{-3} = 4.9787$%, late WHIM excess $\Delta\Omega_{\mathrm{WHIM}} = 0.3252$%, $H_0 = 67.03 \pm 0.22\ \mathrm{km\ s^{-1}\ Mpc^{-1}}$, and $S_8 = 0.776 \pm 0.010$ via the $3.114$% uncondensed horizon-slip fraction, while quarantining post-hoc numerical coincidences ($6\pi^5$, $K=7$, ad-hoc CMB exponents) in Section 9.8.

Revision 2026-11-02.01 adds the complete Axiom-to-Theorem Non-Local Multiverse
to Local Observer Atlas (Subsections 1.1--1.2 and Section 9, with Figure
`nsc-nonlocal-to-local-atlas.png` and summary `nsc-nonlocal-to-local-atlas-v1.json`):
the five openly declared geometric axioms (`A0--A4`), the Circle--Pendulum
smooth zero-crossing theorem ($R(\tau)=R_0\cos(\omega\tau)$ vs. unsigned
$r=|R|\geq 0$ cusp and $(-1)\times(-1)=+1$ across two generations), the unified
5-variable least-resistance cycle (`unified_least_resistance_cycle`), the
zero-parameter two-valve dark-sector fixed point $\chi_\star=3.06697$ at
$\Omega_{\mathrm{vis}}(\chi_\star)=e^{-3}$ ($71.14$% Dark Energy, $23.88$% Dark
Matter, $4.98$% Visible filaments), and the $K=7$ 3D volumetric harmonic
$\eta_7=\tanh(\frac{\pi}{4}e^{-21})=5.955329\times 10^{-10}$.

Revision 2026-10-31.05 puts the outward chain in a compact main argument:
offset geometry, smooth-zero convention, inherited crossing, finite initial
admissibility, P006 vacuum/pressure comparison and child/charge status.
The main argument (excluding cover/contents) changes from 207 pages to7;
the full preserved PDF changes from239 to247 pages. The longer total
retains all prior technical chapters rather than deleting their content.
All634 prior displayed equations, every prior table and all152 earlier
claim markers remain. Two scoped P006 claims are added. Scientific
values and frozen P001--P005 evidence are unchanged.

The moved topics have one-paragraph summaries with commits in Appendix A:
GH engine(A.1), conditioning(A.2), observer readouts(A.3), boson transport(A.4),
dense history(A.5), Schur memory(A.6), kick/shear(A.7), epsilon=.25 controls(A.8)
and matter-work ledgers(A.9). The epsilon=1 nonlinear runs are explicitly
**under-resolved, not used for claims**. Original main parts I--VIII and
the earlier characteristic appendix follow those summaries. The complete
section move list below lets earlier references be traced without changing
mathematical or numerical values. Appendices B and C remain in place.
Author review remains pending and submission_ready remains false.

<details>
<summary>Every moved section: prior identifier to current Appendix A identifier</summary>

| Prior section | Current section | Title |
| --- | --- | --- |
| 1 | A.12 | The central geometric proposal and the route through the atlas |
| 1.1 | A.13 | The proposed generator: contrast, circulation and continuation |
| 1.2 | A.14 | One spatial seed and an accumulated offset generate a space |
| 1.3 | A.15 | The generated space supports a sourced physical readout |
| 2 | A.16 | Why begin with the sphere? Energy, relaxation and inheritance |
| 3 | A.17 | The same argument across operator, spacetime and observer descriptions |
| 4 | A.18 | Four arrows that must not be confused |
| 5 | A.19 | A precise finite inheritance construction |
| 5.1 | A.20 | An infinite inherited law with a controlled local response |
| 5.2 | A.21 | Physical measure, pressure and a measurable nesting ratio |
| 5.3 | A.22 | A response shape that recovers the nesting ratio |
| 6 | A.23 | What an unresolved surrounding region can leave locally |
| 6.1 | A.24 | The local effective Lambda target |
| 7 | A.26 | Start with two holes, without manufacturing two children |
| 8 | A.27 | A failed attachment can sharpen the map |
| 9 | A.28 | A distinct fixed-source wall that does connect |
| 10 | A.29 | A sourced initial slice and a conditional common parent |
| 10.1 | A.30 | Two unequal sourced children in the same theory |
| 10.2 | A.31 | A retained grandchild and an explicit second boundary |
| 10.3 | A.32 | A central proper-clock and frequency observation |
| 10.4 | A.33 | Radar and a wall clock recover a physical junction |
| 10.5 | A.34 | Both walls in one host history |
| 10.6 | A.35 | The shared parent evolves into a calculable effective source |
| 10.7 | A.36 | A compact parent and the limits of its global average |
| 10.8 | A.37 | The nested body in an expanding parent |
| 10.9 | A.38 | Two explicit ways surrounding geometry enters local gravity |
| 10.10 | A.39 | Can the enclosing wall select the cap curvature? |
| 10.11 | A.40 | The source law behind a growing mass |
| 10.12 | A.41 | Sourced rooms from an occupied wall spectrum |
| 10.13 | A.42 | Two physical interiors, one parent, and a central beacon |
| 10.14 | A.43 | Energy, frequencies and a moving radius |
| 11 | A.44 | Complete one reference before assigning its causal arrows |
| 12 | A.45 | Different exteriors can be locally indistinguishable |
| 13 | A.47 | The shared source and its cosmological domain |
| 13.1 | A.48 | A field law that changes the inherited vacuum |
| 13.2 | A.49 | Three vacuum regions in one evolving cosmology |
| 13.3 | A.50 | Preparing the scalar, membrane and geometry together |
| 13.4 | A.51 | Physical boundary reduction and the prepared child |
| 13.5 | A.52 | Real parent sources and received energy |
| 13.6 | A.53 | A finite parent source with its gravitational response |
| 13.7 | A.54 | A conditional inward pulse from the excited parent |
| 13.8 | A.55 | A one-sided geometric observer and the full mass balance |
| 13.9 | A.56 | What an additional inherited region can identify |
| 13.10 | A.57 | Three sourced regions on one compact initial slice |
| 13.11 | A.58 | Cosmology in the protected parent domain |
| 13.12 | A.59 | A constrained parent perturbation that grows |
| 13.13 | A.60 | Inheritance must also supply the structure scale |
| 13.14 | A.61 | What the inherited law leaves to the state |
| 13.15 | A.62 | The inherited scalar gives a controlled matter background |
| 13.16 | A.63 | The controlled background predicts an interior light cone |
| 13.17 | A.64 | An epoch-independent vacuum, scalar and curvature constraint |
| 13.18 | A.65 | Continuing the same inherited parent to the emitters |
| 13.19 | A.66 | The inherited geometry meets the supernova distance shape |
| 13.20 | A.67 | One calibrated geometry, a second distance channel |
| 13.21 | A.68 | What this connection closes and what its source still must supply |
| 13.22 | A.69 | One common ruler and the remaining SI selector |
| 13.23 | A.70 | A physical volume connects normalization and evolving matter |
| 13.24 | A.71 | An extremal member is not yet an inherited selection |
| 13.25 | A.72 | Selection mechanisms and the actual evolving-volume distance test |
| 13.26 | A.73 | Survival ranking and the curvature that must be inherited |
| 13.27 | A.74 | The wall, the exit state and the actual pressure history |
| 13.28 | A.75 | A mass-free nucleation radius and the matter-reference gate |
| 13.29 | A.76 | Three targets from the spectral and regional map |
| 13.30 | A.77 | Constructing the pen: a phase history on the sphere |
| 13.31 | A.78 | The inward drawing: golden turns and inherited clocks |
| 13.32 | A.79 | Through the horizon: smooth geometry, charge and an offset |
| 13.33 | A.80 | Deriving the charged pair from internal geometry |
| 14 | A.82 | The energy supply and formation of structure |
| 14.1 | A.83 | A profile selected by scalar self-gravity |
| 14.2 | A.84 | The coherent core and its inherited environment |
| 14.3 | A.85 | Supplying a local core from the inherited parent |
| 14.4 | A.86 | A perturbed centre and the same evolving cosmology |
| 14.5 | A.87 | Evolving the inherited field and its geometry together |
| 14.6 | A.88 | A diffuse cloud develops a compact, breathing concentration |
| 14.7 | A.89 | Conservation selects a second local scale |
| 14.8 | A.90 | A mass--reservoir--vacuum relation and a dimensional example |
| 14.9 | A.91 | Evolving the core together with its reservoir |
| 14.10 | A.92 | Conserved transport, finite persistence and a relativistic check |
| 14.11 | A.93 | Testing the breathing core beyond spherical symmetry |
| 14.12 | A.94 | A finite response theorem for every orientation in a declared span |
| 14.13 | A.95 | The surrounding field must carry the tidal phase |
| 14.14 | A.96 | A local concentration carries the signature of its surroundings |
| 14.15 | A.97 | The age and population are part of the correlation |
| 14.16 | A.98 | Does the inherited parent leave enough time to form a core? |
| 14.17 | A.99 | A compact earlier seed identifies the nonlinear continuation |
| 14.18 | A.100 | Nonlinear formation without installing the final core |
| 14.19 | A.101 | Physical contraction separates formation from density contrast |
| 15 | A.103 | One geometry, several measured consequences |
| 15.1 | A.104 | Three freely falling clocks and a withheld prediction |
| 15.2 | A.105 | A signal from the concentration's evolving geometry |
| 15.3 | A.106 | The same concentration bends and focuses light |
| 15.4 | A.107 | A proper-clock lens map during contraction and breathing |
| 15.5 | A.108 | Two paths, one receiving clock |
| 15.6 | A.109 | A local image delay and the parent's vacuum must share one scale |
| 15.7 | A.110 | The forming field supplies one connected observer map |
| 15.8 | A.111 | Direction, delay and frequency are derivatives of one arrival map |
| 15.9 | A.112 | Motion calibrates a prediction before the new ray is evaluated |
| 15.10 | A.113 | NSC-P004 passes its frozen synthetic cross-observable check |
| 15.11 | A.114 | Two observer distances, on two metrics |
| 15.12 | A.115 | A stretched region, an orbit and a received frequency |
| 15.13 | A.116 | Compensated circulation, a light offset and smooth-horizon exchange |
| 15.14 | A.117 | Depletion changes area: the inherited source and a dynamical first law |
| 15.15 | A.118 | The field continues: two inherited traces and a receiving state |
| 15.16 | A.119 | The zero line: thermodynamic decay, retained state and a crossover test |
| 15.17 | A.120 | The inherited wave changes geometry, spectra and coupling |
| 15.18 | A.121 | Coupled parent feedback closes matter work and links two clocks |
| 15.19 | A.122 | Finite-amplitude geometry exists for the inherited preparation |
| 15.20 | A.123 | Evolving geometry connects pressure work to local readouts |
| 15.21 | A.124 | The inherited interaction supplies a nonzero canonical response |
| 16 | A.126 | Responses across matter, gravity and quantum descriptions |
| 16.1 | A.127 | Flow geometry and a common boundary stress |
| 16.2 | A.128 | Two stable material interfaces in one static geometry |
| 16.3 | A.129 | The gravitational action now composes through the interfaces |
| 16.4 | A.130 | The same exterior mass permits different gravitational returns |
| 16.5 | A.131 | A published nonlinear vortex mechanism enters the atlas |
| 16.6 | A.132 | Fine motion produces an exact resolved circulation |
| 16.7 | A.133 | The correlations retain their energy and momentum budget |
| 16.7.1 | A.134 | Finite circulation now carries a gravitational readout |
| 16.7.2 | A.135 | Energy depletion determines an averaged geometric work term |
| 16.7.3 | A.136 | A rotating horizon supplies an energy-transfer threshold |
| 16.7.4 | A.137 | The average work coefficient becomes a local Einstein constraint |
| 16.8 | A.138 | A material action now supplies the interaction vertex |
| 16.9 | A.139 | An actual gravitational mode carries the cosmohedron seed |
| 16.10 | A.140 | Interior geometry changes the interaction and fixes its scale law |
| 16.11 | A.141 | The actual cap--wall Dirac operator and the missing nested slip |
| 16.12 | A.142 | Local Dirac trace on the cap and material wall |
| 16.13 | A.143 | An actual wall radar chart and a finite clock prediction |
| 16.14 | A.144 | Stress vertices and the inherited quantum preparation |
| 16.15 | A.145 | Prepared Dirac propagation, mean metric and clock |
| 16.16 | A.146 | A conserved compact source and local gravitational tides |
| 16.17 | A.147 | A correlated continuum parent and local vacuum regularity |
| 16.18 | A.148 | Inherited preparation and the cap's local late state |
| 17 | A.150 | What the eliminated geometry contributes |
| 17.1 | A.151 | The current scalar's classical vacuum projection |
| 17.2 | A.152 | A conditional quantum contribution and its physical scale |
| 17.3 | A.153 | The quantum interface and the remaining exterior reduction |
| 17.4 | A.154 | The exterior scalar response on the quiet wedge |
| 17.5 | A.155 | Parent response, slice fit and a separate contracting geometry |
| 17.6 | A.156 | The interface required by a common action |
| 17.7 | A.157 | One state at two physical resolutions |
| 17.8 | A.158 | A prepared-state change propagated to a photon measurement |
| 17.9 | A.159 | A constrained gravitational exterior and its validity test |
| 17.10 | A.160 | The inherited spectral profile and its gravitational response |
| 17.11 | A.161 | Flat-reference curvature channels and the relation between motion and lensing |
| 17.12 | A.162 | Composition and conservation of a nested geometric response |
| 18 | A.164 | Physical inheritance carries law, scale and state |
| 19 | A.165 | Geometric spectra, quantum observables and causal scales |
| 20 | A.166 | The connected dark-sector mechanism and its remaining physical inputs |
| 21 | A.167 | Charged matter and the asymmetry question |
| 22 | A.168 | Formation, continuation and successive generations |
| 23 | A.169 | The connected result and the route to a tested cosmology |
| A.1 | A.171 | A first field that reaches the interior |
| A.2 | A.172 | Where the early signal can originate |
| A.3 | A.173 | A prepared pulse leaves a finite scattering tail |
| A.4 | A.174 | What the signal determines when preparation is unknown |
| A.5 | A.175 | The transmitted probe and the wall's stress |
| A.6 | A.176 | The prepared pulse requires gravitational initial data |
| A.7 | A.177 | An invariant mass response on the same initial surface |
| A.8 | A.178 | The scalar adds no normal-force jump |
| A.9 | A.179 | Exterior mass closes from the known wall field |
| A.10 | A.180 | Cap mass and a regular wall displacement functional |
| A.11 | A.181 | The prepared field drives a smooth exterior metric response |
| A.12 | A.182 | The same field displaces the wall and changes its exterior clock |
| A.13 | A.183 | The cap geometry changes the central arrival clock |
| A.14 | A.184 | Geometry supplies an outgoing scalar boundary correction |
| A.15 | A.185 | The exterior scalar correction propagates on the computed geometry |
| A.16 | A.186 | Reflected cap transport gives the central scalar correction |
| A.17 | A.187 | A nonlinear coupled window and its unresolved accuracy |
| A.18 | A.188 | Profile-knot phase and a finer temporal comparison |
| A.19 | A.189 | A common exterior derivative and a mixed accuracy result |
| A.20 | A.190 | A finite-amplitude onset solution |
| A.21 | A.191 | Amplitude dependence on the onset cell |
| A.22 | A.192 | Uniform geometric remainders and signed endpoint response |
| A.23 | A.193 | The completed full-label numerical trial |
| A.24 | A.194 | A separate finite model with geometric feedback |
| A.25 | A.195 | Memory is a consequence of elimination, not an added reservoir |
| A.26 | A.196 | The retained local incoming gate |

</details>

Revision 2026-10-28.3 adds Section 13.22: the full shared-ruler DESI
comparison with propagated SN shape posterior, a high-redshift scalar/
Cauchy bound and explicitly conditional Planck geometry diagnostic, and
the computed SI route ladder. Fixed checkpoint source data and P004/P005
are unchanged. The stronger joint score is 1203.22 versus 1191.18 with
equal fitted dimension. A fixed microscopic coupling-ratio length is
distinguished from the unselected parent-radius hierarchy; conditional
quantum-period examples and target-inverted lengths retain those labels.
Submission readiness and author review remain unchanged.

Revision 2026-10-28.2 adds Sections 13.18--13.21 and Figure 18. A finite
backward scalar/causal proof provides the earlier emitters required by
Pantheon+, followed by an identical-covariance LambdaCDM control and a
SN-fixed DESI radial/transverse comparison. Earlier failed domains remain
visible. The source, flux and curvature are unchanged; the observer epoch
is calibrated. The five remaining physical/numerical gaps are explicit.
Later figure labels advance by one. Author review and submission readiness
remain unchanged.

Revision 2026-10-28.1 adds Sections 13.15--13.17 and Figure 17:
an epoch-independent density-fraction relation, a uniform inherited-scalar
matter-background theorem, its rational
luminosity-distance enclosure, and NSC-P005's independent unaveraged
comparison. The original source normalization is retained and the
high-frequency conclusion requires no long carrier evolution. Existing
figure labels 17--38 advance to 18--39. The background proof's domain is
A>=1; ancestral selection, early history and halo/perturbation accuracy
remain separate. Author review and submission readiness are unchanged.

Revision 2026-10-27.4 rewrites the atlas as one connected book. Eight
parts carry spherical geometry and inheritance through physical sources,
formation, observation, gravitational/quantum exchange and effective gravity.
The scalar/four-form parent history is the central worked sequence. The
characteristic-field development moves to Appendix A, with evidence in
Appendix B. The abstract, opening argument, chapter transitions and closing
synthesis are rewritten; source attribution remains beside its use.

All 439 original displayed equations, 102 existing claim markers and 38
figure assets are retained. Sections and figures are renumbered by their
new reader order, with claim bindings updated. The spherical energy argument
is stated with its isotropic fixed-volume assumptions; repeated spherical
inheritance has a separate proposal binding. This is an editorial synthesis,
not a new vacuum-selection, quantum-foundations or observational result.
Earlier section numbers refer to the recorded .3 edition and its Git history.

Revision 2026-10-27.3 adds Sections 10.64--10.66 and Figure 34:
an explicit zero-background surface probe, actual normalized gravitational
mode, action-derived scalar-current vertices and selected cosmohedron
exchange kernel. A same-exterior child-radius contrast and the physical
G/M^2 coupling flow are calculated. The result is a finite projected
tree contribution; full constrained four-point and backreaction work
are not implied.

Revision 2026-10-27.2 adds Sections 10.61--10.63 and Figure 33:
the official OpenAI Navier--Stokes paper, an independently derived finite
nonlinear mean-circulation law and its momentum/viscous budget. Published
blowup, the prepared unforced initial-rate result and future gravitational
identification have distinct domains and attribution.

Revision 2026-10-27.1 adds Sections 10.57--10.60 and Figure 32:
the boundary-fluid dictionary, two sourced stiff-scalar interfaces, a
physical axial gravitational Schur response and a same-exterior return
contrast. The radial and axial stability domains, causal front and physical
unit example are explicit. No rotating cosmology or observed vacuum value
is inferred from the static benchmark.

Revision2026-10-26.1 adds Sections10.55--10.56 and Figure31: the
stationarity precheck, finite timelike calibration and prospectively frozen
NSC-P004 image forecast. The independently evaluated moving-observer
null ray passes the unchanged comparison tolerance. Synthetic validation,
empirical calibration and unique nesting identification remain distinct.

Revision2026-10-25.6 adds Sections10.53--10.54 and Figure30: a fixed-observer
arrival map from the actual formed source. Direction, proper delay and
frequency satisfy two independently tested endpoint identities. A frozen
potential control exposes the frequency's sensitivity to source evolution.
The held-clock and weak-metric assumptions are explicit.

Revision2026-10-25.5 adds Sections10.51--10.52 and Figure29: compatible
finite-amplitude compact data and nonlinear concentration from a weak seed.
Physical mass-shell contraction separates the 10% formation example from
the still-expanding high-contrast 5% comparison. The source mass/radius
conversion is conditional, with no astronomical fit.

Revision2026-10-25.4 adds Sections10.49--10.50 and Figure28: a finite-future
all-wavelength growth bound on the inherited parent and a compact earlier
seed without an installed core. Its controlled linear amplification identifies
where a nonlinear formation calculation becomes useful. Earlier parent
completion and physical seed selection remain open.

Revision2026-10-25.3 adds Sections10.47--10.48 and Figure27: the inherited
core/reservoir spatial kernel, a uniform initial envelope-error bound,
vacuum/mass/turnover relation and age dependence from the saved evolution.
Projection and conditional population formulas identify the next route
to observations without claiming an abundance or full sky spectrum.

Revision2026-10-23.3 adds Sections10.36--10.37 and Figure22: nonradial
lensing on the same evolved scalar parent, the direct null-contraction
identity, Weyl shear and independent ray-map/Jacobi checks. The experiment
uses specified proper clocks and screens; a conventional astronomical
source/observer image fit and population inference remain open.

Revision2026-10-23.2 adds Sections10.33--10.35 and Figure21: full finite
scalar/gravity evolution, a same-mass diffuse concentration test, and an
operational received frequency after contraction. The moving outer grid
retains the protected parent and all outgoing radiation. Finite breathing
is distinguished from asymptotic or nonspherical stability.

Revision2026-10-23.1 adds Sections10.31--10.32 and Figure20: a nonlinear
compensated core, inherited reservoir and outward-current pulse, plus the
surrounding parent domain whose evolution remains unchanged. Local
redistribution supplies structure while preserving the inherited vacuum.
Initial compatibility and causal evolution are distinguished from the
still-missing complete formation history.

Revision 2026-10-22.1 adds Sections 10.29--10.30 and Figure 19: physical
profile selection by scalar self-gravity, a controlled local cooling
example, and coherent admission into the shared atlas. It calculates how
the inherited ambient phase changes the geometric mass. The coupling
choice and the missing full parent formation history remain explicit.

Revision 2026-10-21.7 adds Sections 10.27--10.28 and Figure 18: the finite
state-flow obstruction to automatic phase selection, an operational
free-fall clock transfer, and NSC-P003's third-clock prediction after two
state calibrations. Independent metric variation checks the clock formula.

Revision 2026-10-21.6 adds Sections 10.25--10.26 and Figure 17: constrained
same-source perturbations of the protected parent, measured gravitational
growth, independent scalar-phase response, and the connection between the
vacuum and structure scales. The approximation domains remain explicit.

Revision 2026-10-21.5 adds Sections 10.23--10.24 and Figure 16: a compact
three-region scalar initial atlas, a source-independent inherited-law
obstruction, and the actual protected parent's cosmological evolution.
The whole two-wall future remains open; the protected smooth domain and
conditional dimensionless vacuum are calculated under the same action.

The one living manuscript is **Nested-Space Cosmology: The Geometry of a Universe**:
[Markdown source](../paper/nested-space-cosmology.md) and
[illustrated PDF](../paper/nested-space-cosmology.pdf). Its revision is independent
of the software release. Future writing updates these files. The
[catalog](../paper/catalog.json) identifies current, supporting and historical editions.

Revision 2026-10-21.1 adds Sections 10.16--10.17 and Figure 12: a coupled
regular scalar/gravity/membrane preparation, its physical boundary
Schur reduction, and its Lorentzian wall and open-child continuation.
The selected branch's homogeneous excitation is tiny, so it does not
realize the earlier rapid-transfer matter component. The abstract now
summarizes the connected results and their domains rather than listing
every historical benchmark. Quantum transition probability remains open.

Revision 2026-10-20.6 adds Sections 10.14--10.15 and Figure 11: an
attributed four-form/charged-membrane action supplies vacuum inheritance,
two compatible moving boundaries and an expanding cosmological join.
Free-fall tidal traces and the equal-charge second difference define
conditional observer tests. Absolute flux/scale selection and formation
remain open; this separate candidate does not alter the earlier Schur zeros.

Revision 2026-10-20.5 adds Sections 10.12--10.13: a declared occupied
surface-Dirac material, its conserved pressure, a radially stable
horizonless reference and nearby states under one fixed action.
Two unequal sourced interiors fit one compensated dust parent, and
their central beacons reach cosmological observers with distinct
calculated redshifts. The new plate connects spectrum, material,
geometry and signal. It keeps the supplied occupation and cap vacuum,
the unstable wormhole orientation, and the unresolved shape,
fragmentation and formation questions explicit.

Revision `2026-10-20.4` develops the illustrated-atlas direction with a
reader's construction map and Sections 10.8--10.11. Einstein--Straus
matching gives the retained exterior an expanding dust parent and a
clock prediction. Attributed brane examples show induced curvature and
shadow-matter lensing. The existing constant-tension wall has no stable
equilibrium selecting Lambda; its exact critical-mass bound and the
remaining inheritance invariants are displayed. Gravastar material
laws and established throat/interior-universe work are credited as
building blocks. The new four-panel plate compares these calculations
without identifying their different actions or backgrounds.

Revision `2026-10-20.3` adds Section 10.7's explicit compact S3 parent,
its intrinsic geometry figure, and the matched remote-nesting comparison.
Global effective traces and physical volume fractions change while the
retained local cap signal remains identical. A quiet-throat lifetime and
necessary bounce-stress condition identify the next continuation problem.
This establishes a finite atlas connection, not a generated cosmological
Lambda, a regular bounce, or observational validation.

Revision `2026-10-20.2` adds Section 10.6's transported normal-clock parent.
The whole quiet nested slice evolves with material-wall Jacobi impulses;
its effective trace, explicit source exchange and delayed cap-clock signal
come from the same history. Conservation and independent ray derivatives
are checked. The effective source is non-vacuum, and the isolated parent's
large-domain trace vanishes at fixed finite time. Cosmological parent
construction and the origin of the regional vacuum values remain open.

Revision `2026-10-20.1` adds Section 10.5's shared two-wall host history.
Both original wall intersections retain their common chart normalization,
and incoming characteristics locate source collars on the initial slice.
The finite segments have no causal path to one another within the
constructed host patch. Equal elapsed proper times do not define
simultaneity. This supplies the geometry for later field transfer;
physical foliation, field state and cosmological averaging remain separate.

Revision `2026-10-19.3` adds Section 10.4's physical radar/clock inverse:
the quiet cap lapse and signed junction combinations are reconstructed
from calibrated clocks, then predict a later received frequency. Both
existing wall states recover the same membrane tension. Conditional input
bounds, an alternative material law and independent proper-time/high-precision
controls keep the scope explicit; no cosmological Lambda is generated.

Revision `2026-10-19.2` adds Section5.3's scale-free inherited inverse and
held-out response coefficient, and Section6.1.21's associative geometric
composition, wall-volume transport and time-dependent effective-source
balance. The finite application uses the existing nested slice; its
nonzero mixing term is not a generated cosmological constant. A compact
source-hashed record and independent controls accompany the two proofs.

Revision `2026-10-19.1` adds Section6.1.20: the actual quiet wall's causal
isolation and a local late-state theorem from inherited Hadamard data.
The full four-spinor kernel, stress memory bound and scalar-density
detector spectrum are calculated. The previous trace now fixes a full
asymptotic one-loop vacuum source at fixed classical couplings. Its
nesting-specific offset against the same-H smooth cap is zero; the
finite parent state, compact-source slip and origin of H remain open.

Revision `2026-10-18.14` adds Section5.2: an explicit volume-weighted
action, its lapse/scale pressure variation and a calibrated dynamic root
observable. The positive field sector fails vacuum form at every inherited
depth, while the dimensionless inertial response distinguishes different
nesting ratios with identical static compliance. The full gravitational
correlation source remains separate.

Revision `2026-10-18.13` proves the infinite finite-energy inherited-law
construction in Section5.1: normalized subtree scaling, an exact branching
Schur equation, a uniform depth error bound and a spectral counting
exponent. The proof states its domain and positivity assumptions; it does
not identify the operator response with a cosmological vacuum source.

Revision `2026-10-18.12` adds the correlated continuum preparation and
its local Hadamard proposition in Section6.1.19. The local cap vacuum
stress exists under those assumptions; its controlled amplitude and full
source-induced Psi/Phi remain uncomputed.

Revision `2026-10-18.11` adds the conserved compact-source family and
its calculated prepared clock mass derivative in Section 6.1.18. The
local tidal ratio changes by about five parts per million between the
nested geometry and smooth-cap control. It states the physical source,
observer, numerical diagnostics and scalar-potential zero-mode distinction;
the full vacuum Psi/Phi remains open. Exact triangular wall quadrature
resolves the clock precision floor without another field resolution.

Revision `2026-10-18.10` adds the actual Dirac transfer and prepared mean
metric/clock response in Section 6.1.17. The nested logarithmic coefficient
is -.41780 versus -.39929 for the smooth-cap control with common observer
motion. It records the physical scale, numerical spread, contracting
mass-path discrepancies, state assumptions and remaining full-slip inputs.
The result is conditional first-order mean backreaction, not a selected
baseline vacuum or a completed Psi/Phi calculation.

Revision `2026-10-18.9` adds the actual Dirac bulk and Euler-wall trace,
the exact wall radar chart, and action-derived stress/metric vertices in
Sections 6.1.14--6.1.16. New Figure 3 shows the finite causal diamond and
the change in radial Dirac potential. A fixed radar receiver gives a
background frequency ratio of 1.153626 versus 1.080686 for the smooth cap
continuation. Quantum preparation and the full curved slip remain open.
Section 13.22 now includes the archived PDE/boundary-lift diagnosis and
the distinct manufactured shared-action control; P001 remains unconverged.
The ten figures are numbered in reading order.

Revision `2026-10-18.8` identifies Section 6.1.12 as an imported flat Dirac
reference, then derives the actual quiet cap's radial probe operator and
the material-interface domain in Section 6.1.13. No nested slip is claimed:
its quantum state and metric-source/readout are still missing. Section 13.22
now locates P001's active retained field inside one unbisected outer cell
on all three levels. Its growing integrated mass defect rules out denser
maximum sampling as the sole explanation. The trial remains unconverged.

Revision `2026-10-18.7` adds Section 6.1.12: both covariant curvature channels
of the same induced profile and a conditional static motion--lensing relation.
The source-action variation gives a gravitational slip of 0.9989169187 at
the existing reference scale. The scalar zero mode adds no volume term.
Physical length calibration, astronomical source inference and the actual
causal state remain separate; no cosmological Lambda is inferred.

Revision `2026-10-18.6` adds Section 6.1.11: the flat quadratic gravitational
response of the existing warped compact spectral profile. At the fixed
observer scale it differs from the local approximation by 1.59 parts per
million of the Einstein response. The local polynomial's extra root is
above the saved ultraviolet cutoff. The manuscript retains the inherited
spectral-error boundary and the missing curved causal/state connection;
no full-kernel stability or Lambda is inferred.

Revision `2026-10-18.5` adds Section 6.1.10: a constraint-compatible vacuum
gravitational normal map and a fixed-coupling stability diagnostic. The
local auxiliary action's flat exterior has a negative squared scale and a
negative quadratic action weight. The manuscript states its exact causal
domain and distinguishes this nonperturbative reference from the incomplete
spectral/state theory and the finite periodic episode. The independent
comparison ruler is not identified with an ultraviolet cutoff.

Revision `2026-10-18.4` adds Section 6.1.9: one initial paired-phase direction
propagated with the complete original conformal geometry to two received
photon frequencies. Both original time steps finish and reproduce the saved
baseline histories. The table reports the received-redshift derivatives and
their observed temporal differences. The nonzero inherited constraint,
physical preparation gap and remaining effective-action connection stay
explicit. The first stopped derivative check is preserved and the unchanged
radial force is checked through its independent initial Fourier action.

Revision `2026-10-18.3` adds Section 13.22: all three original P001 endpoint
runs have completed. The local wall-radius protocol now has a supplied
physical ruler and a 22.85-micrometre observed resolution envelope. The
mass gap grows while the global wall-flux maximum shrinks. The numerical
extension stays unconverged, and the nine-cut split does not identify an
unstored maximum state. The finite onset theorem is not extended by this
trial. No new evolution was run for its comparison or registration.

Revision `2026-10-18.2` adds Section 6.1.8: a fixed physical resolution map
on one saved conformal state, the action-derived vacuum-volume generator,
and a whole-force comparison. The calculated energy coefficient fails
the vacuum-force test. The same action and full fermion columns are held
fixed, while the filtered constraints and additional response remain
explicit. This supplies no cosmological Lambda or fractal dimension.

Revision `2026-10-18.1` adds Sections 10.1--10.3: a second unequal sourced
body, a two-level initial geometry with a retained grandchild, regular bulk
charts and a central wall-emitted flash. Figure 5 shows the actual radius
profile. Section 6.1.5 retains the new hierarchy's negative trace and failed
vacuum-form test. The flash has an exact clock/frequency relation under
declared inputs; no cosmological scale, global evolution or unlimited
inheritance is claimed. All nine figures are numbered in reading order.

Revision `2026-10-17.5` connects the exterior Schur response to whole-body
gravity, geometric averaging and causal distance measurements in Sections
6.1.4-6.1.7. It includes the quiet-cap supernova comparison and the separate
evolving conformal candidate's light and radar distances. The two actions,
their interface requirements and each calculation's finite domain remain
explicit. The original three-grid endpoint comparison remains unfinished;
this working revision does not declare a converged cap-wall spacetime or a
derived cosmological Lambda.

Revision `2026-10-17.4` added Sections 6.1.1-6.1.3: the current candidate's
zero added classical vacuum trace, conditional quantum curvature terms at
the fixed P001 scale, a necessary cap-curvature condition, and the formal
interface heat coefficient. The linked compact record regenerates these
exact calculations. The manuscript distinguishes them from the remaining
full renormalized stress, joint quantum state and cosmological tests.

Revision `2026-10-17.3` added Section 13.21: explicit uniform geometric
remainders and strict fixed-label endpoint signs across the entire admitted
positive-amplitude range. It retains the exact canonical-to-Jordan frame
conversion, material radius and quiet-centre clock convention. The
[complete comparison](nsc-atlas-finite-amplitude.md#certified-finite-amplitude-endpoint-signs)
and [material-radius result](nsc-atlas-finite-amplitude.md#certified-jordan-radius-endpoint-comparison)
use the same original cell as the existence and parameter theorems in
Sections 13.19-13.20. Seven directly linked compact certificates carry
the inputs and quantitative conclusions. The older numerical histories,
later central arrival, physical preparation selection and observational
tests retain their own domains and unresolved requirements.

## The common research thread

The current working manuscript carries the prepared field through the
[smooth exterior response](nsc-atlas-bulk-response.md),
[moving wall](nsc-atlas-wall-evolution.md) and
[cap metric and central null clock](nsc-atlas-cap-response.md).
Independent radial and timelike mass paths check this quadratic geometry.
The central arrival map now includes both travel through the changed cap
and the moved wall event at a fixed exterior label. The higher-order scalar
calculation uses that same source and clock convention.
Section 13.13 now derives and evaluates the cap's outgoing scalar backscatter
and its exterior wall derivative, with five numerical comparisons. This
boundary input is used by Section 13.14, which now computes the exterior
scalar correction with separate forcing-cache, spatial and ODE comparisons.
The [central reconstruction](nsc-atlas-central-scalar-response.md) now
evaluates the formal cubic coefficient at five interior arrival cuts. Its
[radial](../results/development/nsc-atlas-central-scalar-v1/central.json) and
[ray](../results/development/nsc-atlas-central-ray-v1/ray.json) comparisons
retain separate incoming-quadrature and exterior cache/grid/ODE sensitivities.
The [exact nonlinear model](nsc-atlas-finite-amplitude.md) specifies the
source, smooth equations and material interface. A full finite-amplitude
solution and uniform remainder are not established by these coefficients.

The [characteristic-background coupled window](current-result-history-20261102.md#characteristic-background-coupled-evolution-16-october-2026)
now supplies six finite-amplitude numerical histories with a common moving
wall and cap. Carrying the known incoming pulse analytically removes the
dominant cap transport error. Paired wall geometry and clock responses
agree closely under the declared spatial/time refinements, while the
scalar time sensitivity, cap outgoing response and exterior equation
residuals remain explicit. The [comparison figure](figures/nsc-atlas-background-evolution.png)
and saved fields connect this numerical advance to the manuscript without
claiming a uniform error enclosure or central-arrival solution.

Section 13.17 now adds the [fixed-input SSPRK3 comparison](current-result-history-20261102.md#midpoint-complete-ssprk3-comparison-16-october-2026)
and its [saved-history figure](figures/nsc-atlas-ssprk3-refinement.png).
The paired scalar's 256/512-step difference is 1.68% of its fine peak.
An exact local profile-knot identity explains why crossing phase matters
to quadrature, without assigning the full coupled error. The observed
difference ratio, small geometry changes and earlier extrapolant agreement
do not establish convergence order or resolve the remaining spatial,
constraint and computational-support defects.

Section 13.18 adds the [same-grid exterior-operator comparison](current-result-history-20261102.md#fixed-grid-exterior-operator-comparison-16-october-2026).
One quadratic departure derivative drives both radial field reconstruction
and moving-grid transport. Final kinematic defects decrease substantially
and paired wall responses remain close, while the independent near-wall
curvature remains unresolved. The comparison keeps this distinction visible
in the numerical table and links the matching source and evidence records.

Regions, boundaries and fields define possible nested geometries. Attachment and
initial state then determine what an interior observer can access and measure.
Spherical symmetry supplies the first tractable charts; it does not assert global
multi-hole topology or a literal shape for every particle.

The paper moves from motivation and mathematical language to the map, causal
access, calculated responses and missing physical connections. The scalar pulse
and finite transfer/memory results retain their different domains. A mode-space
complement is not automatically an exterior spacetime. Descent, attachment,
influence and correlation remain distinct.

The [observer-inverse proof](nsc-atlas-observer-inference.md) connects the
calculated signal back to a preparation within the fixed geometry. Its
positive compact-source exclusions explain why the next comparison must
state what preparations are physically allowed before inferring an exterior.

The [conformal stress check](nsc-atlas-scalar-wall-consistency.md) connects
the same prepared probe to the wall's quadratic stress conditions. Its
surface coefficient and tangential flux jump vanish; the smooth bulk
gravitational evolution requires compatible initial data. The
[characteristic initialization](nsc-atlas-characteristic-initialization.md)
now computes the required source-leg radius correction with the original
compact preparation, corner geometry and null normalization fixed. This
connects the field to one gravitational constraint; transverse evolution,
mass change and the coupled wall response require additional equations.
The [source-leg mass successor](nsc-atlas-source-mass.md) now supplies
the mixed-radius connection on that same initial surface. The geometric
mass has a negative constant quadratic coefficient after the compact
trace, checked through full frame conversion and improved-stress transport.
The paper distinguishes this local initial-data result from asymptotic
energy and from the subsequent bulk and wall-response calculations.
The [normal-force check](nsc-atlas-scalar-normal-force.md) now adds zero
regular normal-stress jump for the same conformal probe. It supplies a
necessary mean-curvature condition for a later sourced wall, with bulk
metric and wall displacement both included in that variation.
The [wall-mass closure](nsc-atlas-wall-mass-closure.md) now shows that the
fixed leading trace supplies an analytic exterior-side mass functional,
using both null fluxes and the full frame correction. Its numerical
evaluation and compatible cap/wall response are now recorded in the
successor calculations above. Those formal coefficients retain their
perturbative domain; formation, stability and a physical preparation
mechanism remain later questions.

## Supporting and historical editions

- [Foundation v0.26.0 PDF](../paper/archive/foundation-v0.26.0/paper/nested-space-cosmology.pdf),
  [source](../paper/archive/foundation-v0.26.0/paper/nested-space-cosmology.md) and
  [archive receipt](../paper/archive/foundation-v0.26.0/archive-manifest.json),
  retaining metadata, claim bindings and the renderer.
- [Regional Transfer, Geometric Feedback, and Local Memory](../paper/finite-regeneration.pdf):
  the finite realization article, [LaTeX](../paper/finite-regeneration/main.tex)
  and [evidence manifest](../paper/finite-regeneration-manifest.json).
- [Local incoming-gate companion](../paper/local-incoming-gate-draft.pdf):
  the preserved OPEN application checkpoint; its campaign remains paused.
- Atlas [v1 evidence](https://doi.org/10.5281/zenodo.23090015),
  [v2 manuscript](https://doi.org/10.5281/zenodo.23090015),
  [v3](https://doi.org/10.5281/zenodo.23090015), [v4](https://doi.org/10.5281/zenodo.23090015),
  [v5](https://doi.org/10.5281/zenodo.23090015), [v6](https://doi.org/10.5281/zenodo.23090015):
  source-bound historical stages, not competing current drafts.

Historical text retains original logical paths. Materialize its original layout
with `python scripts/replay_history.py prepare` when following old `lab/` links.
The living paper uses root paths. The [reproduction guide](reproducing.md)
distinguishes authentication, numerical replay and publication builds.
No arXiv submission occurred in reassembly.
