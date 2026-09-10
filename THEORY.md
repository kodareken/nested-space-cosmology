# The central idea: one law across nested spaces

Nested-Space Cosmology begins with a geometric proposal: collapse in a parent
space and expansion inside a child space can be two descriptions of a connected
process. A locally complete room inherits the law of the surrounding structure
while expressing it at its own scale.

The organizing condition is

$$
\boxed{\mathcal T_\Theta^*\mathbb D_\Theta=\mathbb D_\Theta.}
$$

Here $\mathbb D_\Theta$ describes the room and its links, $\mathcal T_\Theta$
transports the configuration between scales, and $\Theta$ is the shared set
of dimensionless coefficients. The proposal is that matter, gravitational
response and inherited boundary effects are projections of this same object.
The [README](README.md#research-status) states the programme's research status.

## 1. A room and the spaces it connects to

Write the parent–child Hamiltonian in a common spin frame:

$$
H=\begin{pmatrix}H_p&B\\B^\dagger&H_c\end{pmatrix}.
$$

The diagonal blocks govern the two rooms. The link $B$ couples their fields.
Eliminating the child gives the exact parent response

$$
G_{pp}(E)^{-1}=E-H_p-B(E-H_c)^{-1}B^\dagger.
$$

This Schur complement makes an unresolved room dynamically accessible through
its boundary. Its response depends on energy, geometry and the link. The
calculated [curved boundary maps](docs/nsc-boundary-response.md) and
[full spinor channels](docs/nsc-chiral-boundary.md) evaluate this construction
on the project's smooth geometry.

## 2. A shared coupling and the Dirac mass shell

For a real scalar coupling in sheet-first Weyl order,

$$
H_8=I_2\otimes\boldsymbol\alpha\cdot\mathbf p+\Phi\tau_1\otimes\beta,
\qquad \Pi_-=(I_8-\tau_3\otimes\gamma^5)/2,
$$

$$
[H_8,\Pi_-]=0,\qquad
H_8|_{\Pi_-}\simeq\boldsymbol\alpha\cdot\mathbf p+\Phi\beta.
$$

The invariant sector therefore has

$$
E^2=\lvert\mathbf p\rvert^2+\Phi^2,\qquad
\Sigma_p(E)=\Phi^2(E+\boldsymbol\alpha\cdot\mathbf p)^{-1}.
$$

One coupling sets both the mass shell and the visible self-energy. The
complementary sector has the same spectrum. Physical sector selection and the
scalar interaction belong to the action/domain problem; the evaluated
massless throat's action kernel has vector and axial channels.
[The observable bridge](docs/nsc-observable-bridge.md) fixes the conventions.

Chirality, sheet exchange and charge conjugation have distinct mathematical
roles. In the stated gamma conventions,
$\psi^c=i\gamma^2\psi^*$; its gauge representation is conjugated.
The phase gradient of a wave gives its local frequency and momentum, while
the charge representation determines the particle–antiparticle relation.

## 3. Coherence stores energy and carries transfer

A Gaussian fermion state is specified by
$C_{ij}=\langle c_j^\dagger c_i\rangle$, with $0\le C\le I$ and
$i\dot C=[H,C]$. Define

$$
z_B=\mathrm{Tr}_p(B C_{cp}).
$$

Then

$$
\boxed{E_{\rm link}=2\Re z_B,\qquad
\dot N_p=2\Im z_B=-\dot N_c.}
$$

The real and imaginary parts of one frame-invariant coherence determine
stored link energy and occupation transfer. Energy currents include the room
Hamiltonians; the complete account satisfies
$\dot E_{\rm tot}=\mathrm{Tr}(C\dot H)$.
The [common-source derivation](docs/nsc-common-source-derivation.md) carries
these identities into boundary elimination and metric variation.

## 4. Geometry, state and source form the feedback loop

Regular black-universe geometries supply an established realization of a
trapped region, positive minimum radius and expanding interior. The project
uses their geometry with its own Dirac transport and state construction.
Its prescribed radius-pulse experiment converts geometric work into pairs;
its parent-matched vacuum calculation gives negative radial null stress on
the smooth neck. Both are linked in the [evidence map](docs/current-result.md).

The common action must satisfy all independent metric equations. In a
closed-time-path description this takes the form

$$
\left.\frac{\delta\Gamma_{\rm one}}{\delta g_\Delta^{\mu\nu}}\right|_{g_\Delta=0}=0,
\qquad g_\Delta=g_+-g_-.
$$

The completed source maps identify lapse, shift, radial-metric and radius
forces. The [next defining equation](docs/nsc-causal-source-definition-gap.md)
is the causal geometric contribution with its state or admissible history
prescription. Solving it would couple the current stress calculation to the
geometry's own evolution.

## 5. Recursion carries the same response across scales

At common dimensional energy, parent-normalized elimination gives

$$
\Gamma_p(x)=K_p(x)-\Omega^{-1}b\,\Gamma_c(x/\Omega)^{-1}b^\dagger,
\qquad b=B_{\rm dim}/\Lambda_p.
$$

The child response contains its own descendants. The cutoff ratio
$\Omega=\Lambda_c/\Lambda_p$ and the room parameter
$\zeta=(\Lambda L_\star)^2$ have separate definitions; their relation is a
geometric output to determine. The recursive state must carry occupations and
cross-boundary correlations alongside the response kernel.

In this picture, infinity describes repeated generation of rooms. A physical
clock-transfer law and convergent observable marginals give that interpretation
its mathematical content. The [endpoint](docs/nsc-tail-limit.md) and
[clock calculations](docs/nsc-clock-horizon.md) provide the existing controls.

## 6. From the source to observations

The same source has several observable projections: a homogeneous density and
pressure determine expansion; anisotropic stress and momentum transfer enter
clustering and lensing; poles, residues, spin and charges identify matter.
Internal conversion $Q_b$ cancels when regional energy equations are summed:

$$
\dot\rho+3H(\rho+p)=J_b+J_d.
$$

Here $J_b,J_d$ denote external room supply, including its reservoir and
boundary account. Proper time, volume and deposition relate a boundary power
to a cosmological density rate. The programme seeks an unfitted dimensionless
relation across these sectors from one solved state and parameter set.

## Reading the evidence

The repository uses six labels where useful: **Postulate** for a defining
requirement; **Imported result** for an attributed established formulation;
**Repository derivation** for an equation obtained here; **Numerical diagnostic**
for a computed case; **Open prediction** for an observable still to derive;
and **Interpretive hypothesis** for the proposed physical reading.

[Current equations and evidence](docs/current-result.md) ·
[Source provenance](results/development-snapshot.json) ·
[Prior mathematics](docs/prior-art-and-open-claim.md) ·
[Working paper](paper/nested-space-cosmology.pdf)
