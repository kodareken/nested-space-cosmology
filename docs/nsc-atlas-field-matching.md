# A conformal scalar signal through the fixed atlas wall

This note owns the field model and matching assumptions for the first
conditional response at the central interior observer. The geometry stays the
isolated spherical reference of
[the causal classification](nsc-atlas-causal-classification.md), with
\(H\ell=1\), \(k\ell=2\), \(M/\ell=0.4\), and the expanding wall. The
[observer construction](nsc-atlas-causal-observers.md) supplies the existing
wall and central clocks. No parameter is changed to obtain a signal.

The calculation answers: **given a specified incoming scalar characteristic
trace on the early Schwarzschild side of the wall, what classical field does
the central de Sitter observer receive?** It is an exact transfer inside this
test-field model, evaluated on a finite interval. The incident trace is an
explicit input. Its generation in the exterior remains a separate question.
The result is neither a quantum vacuum selection nor a signal from the
designated right-hand parent observer, whose access is excluded by the
isolated reference causal map. The nonspherical glued parent's full future
is still unknown.

## Borrowed equations and reading boundary

The focused reading is sufficient to import the listed components; it is not
a survey establishing priority for their combination.

| Primary source and fixed version | Read for this connection | Imported component and limit |
| --- | --- | --- |
| [Mars and Senovilla, gr-qc/0201054v1](https://arxiv.org/html/gr-qc/0201054v1) | Distributional curvature and junction discussion, especially Eqs. (53)-(58), (69)-(71), (74) | Ricci-scalar delta coefficient at a non-null shell. This is distributional geometry, not a global evolution or stability theorem for our wall. |
| [Santiago et al., 1512.02120v2](https://arxiv.org/html/1512.02120v2) | Sections II-III through the junction conditions, especially Eqs. (5), (8), (17), (29), (31)-(33) | Independent spherical-shell statement of scalar continuity and the normal-derivative jump, with \(\mathcal R=-8\pi S\delta\) and coefficient \(-2[K]\). Their background is a static charged shell with flat interior, not our moving de Sitter wall. Its quantum vacuum and instability analysis are not imported. |
| [Tomasini, Rubin de Celis and Simeone, 1901.07524v1](https://arxiv.org/html/1901.07524v1) | Sections 2-3 through Eq. (18), especially action (1), field equation (4), curvature discussion after (7), and integrated jump (15) | Nonminimal scalar action and the normal-derivative jump induced by shell curvature. Their examples are static conical/cylindrical spacetimes. Their Green functions, stability classifications and cylindrical reduced-field cancellation at coupling 1/4 are not imported. |
| [Huguet, Queva and Renaud, gr-qc/0603031v1](https://arxiv.org/html/gr-qc/0603031v1) | Sections 4-5, especially Eqs. (21), (24), (30), (36) | Exact conformal weight \(\Omega^{-(d-2)/2}\) of the massless conformal scalar. In four dimensions \(a\phi\) obeys the flat wave equation in the cap's conformal chart. Their parameter named \(\xi\) is a curvature deformation, not our curvature coupling; their metric and curvature sign conventions must also be translated. No global quantum state on the joined spacetime is imported. |
| [Dolan and Ottewill, 1106.4318v1](https://arxiv.org/html/1106.4318v1) | Sections II-III, especially Eqs. (4), (6), (10)-(12) and the discussion separating prompt, ringing and tail contributions | The Schwarzschild scalar radial operator retains its smooth curvature potential. Their usual exterior black-hole Green-function boundary conditions are not automatically the conditions on our white-hole/opposite-exterior wall history. |
| [Chu, 1305.6933v2](https://arxiv.org/html/1305.6933v2) | Section IV result and four-dimensional flat-slicing specialization, Eqs. (80), (82), (85), (87)-(89) | A minimally coupled massless scalar in de Sitter has a direct response and a nonzero interior-light-cone tail. Its exact Green function is a reusable future comparison; the moving-wall boundary problem has not been implemented here. |

The earlier signed-wall equations and global patches remain borrowed from
[Aguirre and Johnson, gr-qc/0508093v2](https://arxiv.org/html/gr-qc/0508093v2),
as recorded in the classification owner. The conformal transformation,
one-dimensional d'Alembert solution, centre regularity, and shell jump are
established mathematics. Applying them with consistent signs, clocks and
input data closes a finite atlas connection; it is not a discovery of those
formulas.

Searches covered nonminimal test scalars on thin shells, spherical scalar
propagation, de Sitter conformal fields and Schwarzschild propagation.
Static shell self-force results supplied a usable distributional matching
equation, rather than a dynamical solution for this wall. Finite-thickness
wall microphysics, self-consistent scalar backreaction, quantum states on
the entire joined history, and a complete characteristic exterior problem
remain unread or unverified boundaries. No exhaustive predecessor/successor
search or worldwide priority claim is made.

## Field action and the shell term

Use signature \((-+++)\), geometrized units, and one unit spacelike normal
\(n\) pointing from the de Sitter side \((- )\) to the Schwarzschild side
\((+)\). A bracket means plus trace minus minus trace using this same normal.
The timelike wall has induced metric \(\gamma_{ab}\), surface stress
\(S_{ab}=-\sigma\gamma_{ab}\), and \(k=4\pi\sigma>0\). With
\(K_{ab}=\tfrac12\mathcal L_n\gamma_{ab}\), the Israel condition gives

\[
[K_{ab}]=-k\gamma_{ab},\qquad [K]=-3k,\qquad
\mathcal R_{\rm sing}=-2[K]\delta_\Sigma=6k\delta_\Sigma.
\]

Here \(\delta_\Sigma\) uses proper Gaussian normal distance. The scalar
curvature \(\mathcal R\) is distinct from areal radius \(R\). As an
independent sign check, the shell stress trace is \(-3\sigma\), so the traced
Einstein equation gives the same scalar coefficient \(24\pi\sigma=6k\).

Choose a real, massless classical test scalar with action, up to an irrelevant
positive normalization,

\[
S_\phi=-\frac12\int\sqrt{-g}\,
  \bigl(g^{\mu\nu}\partial_\mu\phi\partial_\nu\phi
          +\xi\mathcal R\phi^2\bigr)\,d^4x,
\qquad (\Box-\xi\mathcal R)\phi=0.
\]

The selected candidate sets \(\xi=1/6\). This is an explicit field-theory
choice, not a coupling inferred from the atlas. The scalar is a tracer,
distinct from the matter supporting the background wall. Its amplitude is
treated at linear order; a finite-amplitude backreaction bound is not being
claimed. No independent surface scalar source, surface mass, kinetic term or
other direct scalar-wall interaction is added. Such an interaction would
define a different candidate and generally change the transmission law.

A continuous scalar has no delta-prime term. Integrating its wave equation
through the wall gives

\[
[\phi]=0,\qquad [n\phi]=\xi\,6k\phi.
\]

This is the local timelike-shell version of the integrated jump in
Tomasini et al. Eq. (15), with their extrinsic-curvature trace replaced by
our \([K]=-3k\). It does not require importing their static cylindrical
solution. Santiago et al. Eqs. (31)-(33) state the same continuity and
proper-normal-derivative rule for a spherical shell. Applying that local
integration in Gaussian normal coordinates around our moving timelike wall
does not import their static spectral boundary conditions. Smooth mass or
angular terms cannot contribute a delta in this integration.

The background angular junction condition is already known:

\[
[nR]=\beta_+-\beta_-=-kR.
\]

For the spherical reduced field \(\psi=R\phi\), the product rule therefore
gives

\[
[\psi]=0,\qquad
[n\psi]=R[n\phi]+\phi[nR]=(6\xi-1)k\psi.
\]

At \(\xi=1/6\), \([n\psi]=0\). Tangential derivatives of the shared trace
are also continuous, so both radial characteristic derivatives of \(\psi\)
match in the wall's orthonormal frame. The cancellation removes the
distributional radial barrier. It does not remove the Schwarzschild Weyl
curvature, its smooth scalar scattering potential, or scattering from other
parts of the joined geometry. In particular, it is not a statement that
every incoming wave from an unspecified exterior has unit global
transmission.

## Exact interior solution and incoming characteristic

Use the expanding flat de Sitter chart already fixed by the atlas,

\[
ds_-^2=a^2(\eta)(-d\eta^2+dr^2+r^2d\Omega^2),\quad
a=e^{Ht},\quad \eta=-\frac{e^{-Ht}}H,\quad R=ar.
\]

The central observer has proper time \(t_C=t\), with \(t_C=0\) on the v4
cap slice. The conformal field transformation and spherical symmetry imply

\[
(\partial_\eta^2-\partial_r^2)\psi=0,\qquad
u=\eta-r,\quad v=\eta+r.
\]

Regularity at the centre requires \(\psi(\eta,0)=0\), giving

\[
\psi(\eta,r)=G(v)-G(u),\qquad
\phi_C(t)=\lim_{r\to0}\frac{\psi}{ar}
                  =\frac{2G'(\eta)}a.
\]

The outgoing term is necessary to cancel the apparent \(1/R\) singularity
at the centre. It is fixed by regularity, not an independently emitted
central signal.

Let \(U_w\) be the unit future wall tangent and \(\tau\) its proper time.
In the cap chart,

\[
U_w=\dot\eta_w\partial_\eta+\dot r_w\partial_r,
\qquad n_- =\dot r_w\partial_\eta+\dot\eta_w\partial_r,
\qquad \dot v_w=\dot\eta_w+\dot r_w>0.
\]

The specified incoming datum is

\[
I(\tau)=\bigl[(U_w+n)\psi\bigr]_+
       =\bigl[(U_w+n)\psi\bigr]_-
       =2\dot v_wG'(v_w).
\]

The **plus** sign in the derivative is deliberate: the inward-travelling
field is \(G(v)\); differentiating it along \(U_w+n\) isolates that component.
Its future propagation direction is inward, \(U_w-n\). These two statements
must not be conflated. Only this one incoming characteristic is prescribed;
the outgoing characteristic is obtained from the interior solution.

An ingoing radial signal reaches the centre when \(\eta_C=v_w\). A finite
arrival exists only when \(v_w<0\), equivalently \(HR_w<1\). The arrival map is

\[
T(\tau)=t_w(\tau)-H^{-1}\log(1-HR_w(\tau)),\qquad
T'(\tau)=a_C\dot v_w
     =\frac{\gamma(1+v_{\rm pec})}{1-HR_w}>0,
\]

where \(v_{\rm pec}\) is the wall speed in the comoving cap frame and
\(\gamma=dt_w/d\tau\). Here \(T\) denotes central proper arrival time, not
the exterior Schwarzschild time. Consequently the exact conditional transfer
is

\[
\boxed{\quad\phi_C(T(\tau))=\frac{I(\tau)}{T'(\tau)}.\quad}
\]

This is also an independent normalization check:
\(\phi_C\,dT=I\,d\tau\). A late wall event at or beyond \(HR_w=1\) is
outside this arrival map, not a divergent finite-time detector measurement.
No geometric-optics frequency assumption was used.

## Explicit finite classical input

Let \(\tau_*\) be the existing expanding wall event with \(R_w=0.5\ell\),
and freeze \(\delta=0.03\ell\). Set

\[
z=\frac{\tau-\tau_*}{\delta},\qquad
B(z)=\begin{cases}(1-z^2)^4,&|z|<1,\\0,&|z|\ge1,\end{cases}
\qquad I(\tau)=\frac{\varepsilon}{\delta}B'(z).
\]

The small amplitude \(\varepsilon\) has the dimensions of \(R\phi\).
Normalized transfer results are independent of it. The extended \(B\) is
\(C^3\), and \(I\) is \(C^2\), with compact support and zero integral in
wall proper time. It is not described as a \(C^\infty\) bump. Its duration,
centre, profile and sign are declared input data, not inferred cosmological
observations or adjustable matching constants.

Choose the initial cap slice at the retained wall event \(R_w=0.05\ell\),
before the first support event and after the reference wall's past endpoint,
and set \(\phi=0\) and its initial normal derivative to zero throughout
that cap. Equivalently, take \(G\) constant
before the first incoming support label; the constant cancels from
\(G(v)-G(u)\). Prescribe the above incident trace and zero incident trace
outside its support over the finite boundary history needed to reach
\(t_C=0\). This is a classical initial-boundary specification. It neither
sets fluctuations to zero nor selects a global quantum vacuum.

Integrating the characteristic datum gives

\[
2G(v_w(\tau))=\varepsilon B(z)+\text{constant},\qquad
\phi_C(T(\tau))=
 \frac{\varepsilon}{\delta\,T'(\tau)}B'(z).
\]

The pulse support must be checked to lie entirely in the recorded regular
wall history and in \(HR_w<1\); its received endpoint must precede the
chosen final central event. The producer records these inequalities and the
numerical accuracy of the existing wall quadrature. If they hold, the direct
central response is exactly zero before and after its received support,
including at \(t_C=0\) when that event is later than the support. There is
no interior curvature tail for this conformal candidate. This statement
assumes the stipulated incident trace; it does not exclude a later return
generated by an as-yet-unsolved exterior field problem.

## What remains identifiable

The first record can report the arrival interval, the clock-dependent
stretch/compression and the normalized classical central waveform. It
cannot identify an exterior from a waveform when the incoming trace is
freely selectable. On any monotone arrival interval, a compatible target
waveform \(Y(T)\) can be produced by choosing

\[
I(\tau)=T'(\tau)Y(T(\tau)).
\]

If the incident family is restricted to zero integrated trace, the target
must obey the corresponding zero-integral condition. This is a concrete
input-versus-geometry degeneracy, not an inference that all exterior
geometries have the same response. To discriminate geometries, constrain
the incident source independently or add measurements, for example a known
source clock/profile, a second interior observer, or a measured wall-side
trace. The selected polynomial pulse is a calibration state and does not
remove this inverse-problem ambiguity by itself.

The next exterior calculation must propagate independently specified data
on an admissible exterior initial/characteristic surface to this incoming
trace. In Schwarzschild vacuum the scalar curvature vanishes for both
minimal and conformal couplings, while the reduced equation still contains

\[
V_l(R)=\left(1-\frac{2M}{R}\right)
       \left(\frac{l(l+1)}{R^2}+\frac{2M}{R^3}\right).
\]

Borrow the operator before solving any such problem, and specify the
retained patches and their initial boundaries. The ordinary right-exterior
ingoing-at-horizon scattering prescription cannot silently provide our
early white-hole-side incident state.

Minimal coupling is a distinct comparison, with
\([n\psi]=-k\psi\) and cap equation
\((\partial_\eta^2-\partial_r^2-2/\eta^2)\psi=0\).
Chu's exact de Sitter Green function already supplies its bulk tail. The
surface matching and moving-boundary response would still need integration;
the present conformal answer is not transferred to that theory. Geometric
optics remains another option only with a stated wavelength/curvature
regime and controlled crossing errors. It offers no accuracy advantage over
the exact conformal interior calculation used here.

## Compact admission and stopping conditions

| Item | Bounded decision |
| --- | --- |
| Borrowed | The fixed signed-wall trajectory, distributional curvature, nonminimal scalar jump, four-dimensional conformal transformation and regular spherical d'Alembert solution, with the versions above. |
| Missing connection | Turn an allowed early wall-to-centre causal ray into a sourced field response with an explicit transmission law and initial data. Exterior generation of the incident trace is intentionally still open. |
| Decision and unlock | Establish the conditional transfer, its finite received support and normalization, then expose what source information is required to distinguish an exterior. A vanishing signal at the later central event is an outcome, not a reason to retune the pulse. |
| Finite domain | One isolated spherical reference, one linear massless scalar with \(\xi=1/6\), no additional wall coupling, one declared compact pulse, the regular cap, and the response through \(t_C=0\). No many-body evolution or global singular-past initial state is computed. |
| Budget | Evaluate the existing wall trajectory and this explicit transfer at at most 513 distinct wall radii and within 60 CPU seconds per record; paired 32/64 quadrature and fixed support/peak bracketing, with targeted integration checks. No new PDE campaign, parameter sweep or quantum determinant. |
| Stop or redirect | Stop if support exits the verified wall history, touches the arrival-map boundary, or violates the fixed matching conditions. Stop a broader inference if it requires an unprovided exterior source/state. A need to change \(\xi\), wall stress or incident trace starts a new candidate. |
| Verification | Independently check the shell trace and radius jumps, the incoming-characteristic sign, centre regularity, direct \(\eta_C=v_w\) arrivals, \(\phi_C dT=I d\tau\), zero-input/linearity/unit-rescaling controls, compact support and numerical refinement. Passing these checks validates the adapter, not a new universal transmission theorem. |

The fixed-background transfer and the demonstrated input degeneracy satisfy
the [existing admission rule](../../docs/prior-art-and-open-claim.md) as
finite compatibility and observer-map work. The source components keep
their published attribution, and no claim about dark matter, dark energy,
antimatter, branching formation or observed cosmic signals follows from this
test alone.
