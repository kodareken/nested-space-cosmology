# Matching a published torsion action to a published bounce source

The reusable construction is the spin-fluid collapse/bounce system. The
missing NSC interface is its matter source. This comparison uses published
coefficients directly; it does not rederive heat kernels or simulate collapse.
It tests a candidate torsion extension, not the current torsionless operator
or a completed quantum functional.

## Normalize the same connection

Write the connection difference as K, with
\(\nabla_XY=\nabla_X^gY+K(X,Y,\cdot)^\sharp\). Physical torsion is 2K.
Use the component norm \(\|K\|^2=\sum K_{ijk}^2\).

[Pfäffle and Stephan](https://arxiv.org/pdf/1101.1424), Lemma2.5 and
Propositions5.4–5.5, give for one complex four-component Dirac field

\[
R(\nabla)=R_g-\|K\|^2,\qquad
a_2=\frac{1}{4\pi^2}\int\left(-\frac{R_g}{12}+\frac34\|K\|^2\right),
\qquad
\Delta a_4=-\frac{3}{32\pi^2}\int\|\delta K\|^2.
\]

Their equation45 defines \(F_2=\int_0^\infty F(s)ds\), \(F_0=F(0)\).
With \(A=\Lambda^2F_2/(48\pi^2)\), comparison at the same Einstein
coefficient gives

\[
S_{\rm spectral}^{(2)}=-A\int R_g+9A\int\|K\|^2,
\qquad
S_{\rm EC}^{(2)}=-A\int R_g+A\int\|K\|^2.
\]

**Matching consequence:** the torsion stiffness differs by nine. Keep the
same canonically normalized fermion source and the same connection when
making this comparison. A field rescaling changes the fermion coupling too;
it cannot remove the physical ratio.

For an algebraic component \(C K^2+J K\), eliminating K gives
\(-J^2/(4C)\). Thus the induced current-current term is one ninth of its
Einstein–Cartan value in this leading-order, same-source comparison. This
is a conditional matching statement, not a derived NSC fluid equation.

## Cross-check the other paper's conventions

[Iochum, Levy and Vassilevich](https://arxiv.org/pdf/1008.3630) use connection
difference \(4T_{\rm ILV}\), so \(K=4T_{\rm ILV}\). Their axial covector
is \(B=-(3/2)*K\), up to orientation; their Clifford convention also
differs by a factor of i. Equations21 and25 give

\[
\Delta a_2=\frac8{16\pi^2}\int B^2,\qquad
\Delta a_4=-\frac{2/3}{16\pi^2}\int F_B^{\mu\nu}F^B_{\mu\nu}.
\]

Using \(\|K\|^2=6|*K|^2\) recovers the same coefficients. Only the
coefficient conversion was checked; the published calculation is reused.

## What must match before using the collapse result

The spectral coefficients are
\(\gamma_1=3\Lambda^2F_2/(16\pi^2)\) and
\(\gamma_2=3F_0/(32\pi^2)\). The published geometric torsion equation
has operator \(3\gamma_1-\gamma_2d\delta\). Its coefficient scale is

\[
\mathcal M_T^2=\frac{3\gamma_1}{\gamma_2}
=6\Lambda^2\frac{F_2}{F_0}.
\]

Neglecting the derivative term requires the relevant eigenvalues of dδ to
be small compared with this scale. It is not yet a physical Lorentzian mass
or a positivity proof. For the exponential profile F2=F0=1, higher heat
terms still exist; the compact-support cutoff remainder used in the paper
cannot be carried over unchanged.

[Popławski's collapse model](https://arxiv.org/pdf/2008.02136) uses
\(\tilde\epsilon=\epsilon-\alpha n_f^2\),
\(\tilde p=p-\alpha n_f^2\), with \(\alpha=\kappa(\hbar c)^2/32\).
Under the same averaging and a justified algebraic limit, the coefficient
comparison would replace α by α/9. The required quantum current correlation
has not been derived for NSC. The paper's equations29 and34 can be reused
once that source is matched; equation33's phenomenological production rate
must also be supplied by the selected state. A vacuum control is not that
thermal spin fluid.

Further interface conditions are explicit:

- The spectral critical-point calculation fixes volume. Its discarded a0
  volume term cannot be omitted from an unconstrained NSC metric variation.
- The sources above use Riemannian compact domains. Lorentzian reality,
  constraints and source signs need a compatible prescription.
- Chiral bag boundaries are not automatically the transmitting parent-child
  domain. Keep their boundary terms only where the domain matches.
- The current NSC common functional is not a separately weighted bare heat
  action. These coefficient relations concern a candidate spectral extension;
  its place in the same regulated functional must be established first.

## Reuse the microscopic quantum source as well

[Lucat and Prokopec, arXiv:1512.06074v1](https://arxiv.org/html/1512.06074v1)
provide a closer state-level starting point: a renormalized in-in 2PI Dirac
functional and semiclassical backreaction. Reuse equations17,21,24,26 and
33 rather than constructing that general machinery again.

Their local interaction coefficient is
\(\alpha_5=3\pi G_N\xi^2/2\), with minimal Einstein–Cartan at \(\xi=1\).
The conditional one-ninth matching above therefore maps to

\[
|\xi_{\rm eff}|=1/3,\qquad
\alpha_5^{\rm match}=\pi G_N/6.
\]

These are coefficient conversions for the candidate's algebraic limit,
not an adopted value for the completed NSC action. The source-sign and
Lorentzian/domain conditions above still apply.

With \(X=iS^{aa}(x,x)\) and \(A^\mu=\gamma^5\gamma^\mu\), the two
contractions in their equation17 have structure

\[
\mathcal W[X]=\operatorname{tr}(XA^\mu)\operatorname{tr}(XA_\mu)
-\operatorname{tr}(XA^\mu XA_\mu).
\]

Use this inside the renormalized functional, including its counterterms.
A vanishing mean axial current does not justify omitting the second term.

| Interface | NSC ownership still required |
|---|---|
| G_N and interaction coefficient | Same action normalization; no second gravitational weight |
| Scalar/pseudoscalar mass | Physical pole/sector of the full operator; not the radial gap or subtraction mass |
| Contour propagators and initial state | Full spinor, geometric measure and parent-child state correlations; a retarded map alone is insufficient |
| Finite renormalization terms | Match their dimensional prescription to the NSC regulator, including the interaction's curvature–fermion term R psi-bar psi, before using absolute stress |
| Room boundary and recursion | Supply the NSC transmission domain and scale map; neither is provided by a homogeneous cosmological example |

The paper's sectionVI is a required reuse condition: its near-equilibrium
relativistic bounce reaches the effective scattering cutoff; the authors
argue for a controlled nonrelativistic nonequilibrium regime. Its illustrated
large-coupling examples do not establish a bounce at our conditional matched
coupling. Mass, state and scale separation must be fixed before a run.

The next computation must resolve one of those NSC inputs or evaluate the
matched source where the imported approximation applies. Rebuilding either
the collapse geometry or the general 2PI derivation would not do that.
