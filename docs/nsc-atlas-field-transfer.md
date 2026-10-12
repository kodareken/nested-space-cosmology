# A prescribed early scalar pulse and the central response

This successor fixes one classical test field and one incoming state on the
previous [causal reference](nsc-atlas-causal-observers.md). It calculates the
central response conditional on an exterior-side wall datum. It does not
derive that datum from the right parent, a remote source or a quantum state.

## Admission before calculation

- **Borrowed:** the v5 fixed de Sitter/Schwarzschild wall and proper clocks;
  the nonminimal scalar jump, regular-centre condition and exterior radial
  potential in [1512.02120v2](https://arxiv.org/html/1512.02120v2), Section III,
  equations (8), (26)--(34); the conformal reduction in
  [gr-qc/0603031v1](https://arxiv.org/html/gr-qc/0603031v1), equations
  (21), (24), (30) and (36), independently reviewed for this successor.
  The static charged-shell model in the first source is not substituted for
  our moving pure-tension geometry: its distributional matching law is.
- **Missing:** geometric access alone did not specify a propagating field,
  incoming data, transmission law, or central waveform.
- **Decision and unlock:** determine the finite central scalar response to
  one declared early wall proper-time pulse, and whether this particular
  state leaves a central field after the pulse. This supplies a conditional
  response operator for later comparison with physically supplied inputs.
- **Fixed theory:** four-dimensional real massless conformal scalar,
  angular mode $l=0$, $\xi=1/6$, no scalar self-interaction or additional
  scalar term in the wall action. The source remains $H\ell=1$,
  $k\ell=2$, $M/\ell=0.4$. Geometry is held fixed at linear probe order;
  scalar stress and backreaction start at quadratic amplitude order.
- **State:** zero field and its time derivative in the retained initial cap
  at the v5 wall event $R/\ell=0.05$. Prescribe only the exterior-side
  incoming characteristic $I=(\partial_{\tau_W}+n)(R\phi)$ at the wall.
  Its compact pulse is centered on the wall proper time at $R/\ell=0.5$
  and has half-width $\delta=0.03\ell$, as defined below. Its origin is
  uncomputed. The incoming datum is zero outside this support throughout
  the recorded response window; this includes any otherwise unspecified
  exterior incoming tails.
- **Finite budget:** at most 513 displayed or checked wall radius samples,
  fixed endpoint/peak bracketing and paired 32/64 quadrature, reusing the
  existing v5 clock quadrature. Sixty CPU seconds; no PDE evolution, mass
  sweep or global exterior scattering calculation. The central window ends
  at $t_C=0$. Stop on loss of support coverage, a nonpositive clock map,
  an unsupported return to the wall, or failed scalar matching.
- **Independent checks:** scalar jump from the surface stress trace and
  from the extrinsic-curvature trace; regular centre and radial wave
  identity; proper-time input normalization; an analytic input norm;
  paired quadrature; coordinate and length rescaling; finite-domain exits;
  signed area and the postpulse field. Squared amplitude integrals are
  diagnostics, not stress-energy integrals. No priority claim is made.

## Field law and the wall

With metric signature $(-+++)$ choose the action

$$
S_\phi=-\frac12\int\sqrt{-g}\left((\nabla\phi)^2+
\frac16\mathcal R\phi^2\right)d^4x,
\qquad (\Box-\mathcal R/6)\phi=0.
$$

Here $\mathcal R$ is the four-dimensional Ricci scalar, while $R$ is
areal radius. The normal $n$ points from the cap toward the retained
exterior. The existing pure-tension jump is $[K_{ab}]=-k h_{ab}$,
$[K]=-3k$ and $[nR]=-kR$. Its curvature distribution is
$\mathcal R_{\rm sing}=-2[K]\delta_\Sigma=6k\delta_\Sigma$.
The borrowed scalar conditions therefore give

$$
[\phi]=0,\qquad [n\phi]=6\xi k\phi,
\qquad [\chi]=0,\quad [n\chi]=(6\xi-1)k\chi,
\qquad \chi=R\phi.
$$

At the fixed conformal value, $[n\chi]=0$. This cancellation concerns
the distributional interface potential. It does not remove the smooth
Schwarzschild radial potential $f(R)2M/R^3$ or solve exterior scattering.
An additional scalar surface interaction would change these conditions.

## One incoming characteristic, one regular centre

In the expanding flat de Sitter chart, $a=e^{Ht}$ and
$\eta=-e^{-Ht}/H$. The reduced spherical conformal equation is
$\chi_{,\eta\eta}-\chi_{,rr}=0$. Put $u=\eta-r$ and $v=\eta+r$.
The regular solution for the chosen incoming data is

$$
\chi(\eta,r)=F(v)-F(u),\qquad
\phi_C(t_C)=\frac{2F'(\eta_C)}{a_C}.
$$

The subtraction is required: an ingoing spherical pulse alone is singular
at the centre. The expansion is
$\chi=2rF'(\eta)+r^3F'''(\eta)/3+O(r^5)$ where the derivatives exist,
so $\phi$ has a finite even central limit.

Write $v_{\rm p}$ for wall peculiar velocity, to distinguish it from the
incoming null label $v$. The incoming and outgoing wall characteristics are

$$
I=(\partial_{\tau_W}+n)\chi
 =\frac{2\gamma(1+v_{\rm p})}{a_W}F'(v_W),
\qquad
O=(\partial_{\tau_W}-n)\chi
 =-\frac{2\gamma(1-v_{\rm p})}{a_W}F'(u_W).
$$

Only $I$ is prescribed. The centre and initial data determine $O$; setting
both independently would overprescribe this interior problem. Continuity
of the two reduced field derivatives identifies the incoming datum on the
exterior and interior sides of the wall.

## Fixed proper-time pulse and exact transfer

Let $\tau_c$ be the wall proper time at $R/\ell=0.5$, and define

$$
z=\frac{\tau_W-\tau_c}{\delta},\qquad
B(z)=\begin{cases}(1-z^2)^4,&|z|<1,\\0,&|z|\geq1,\end{cases}
\qquad I(\tau_W)=\frac{\epsilon}{\delta}B'(z),
\qquad \delta=0.03\ell.
$$

The primitive is $C^3$ and the input is $C^2$, sufficient for the stated
classical field equation. The amplitude $\epsilon$ is formal and small;
reported waveforms are divided by it. Since $v_W$ increases strictly,
$F(v_W(\tau_W))=\epsilon B(z)/2$ defines a compact primitive with no step.
The input is a wall proper-time pulse, not a pulse chosen in centre time.

For $R<\ell$, its central arrival map and derivative are

$$
T(\tau_W)=t_W-\ell\log(1-R/\ell),\qquad
T'=\frac{\gamma(1+v_{\rm p})}{1-R/\ell}>0,
\qquad
\phi_C[T(\tau_W)]=\frac{I(\tau_W)}{T'}.
$$

Thus geometric clock stretching changes the waveform and its pointwise
amplitude. The signed area obeys
$\int\phi_C\,dt_C=\int I\,d\tau_W=0$. This is not an energy balance.
The separate quantities $\int I^2d\tau_W$ and
$\int\phi_C^2dt_C$ measure squared amplitude, not energy or flux.

The finite certificate must check this pulse's actual support. The earliest
central arrival must occur after the wall passes $R=0.8\ell$. Thereafter,
until $t_C=0$, $r_W=R/a_W\geq0.8\ell$. If the earliest incoming null label
is $v_-$, the reflected wave occupies at most
$r=-\ell-v_-$ by this endpoint. Checking $-\ell-v_-<0.8\ell$ excludes
a centre-reflected return to the wall within the window. It does not
eliminate exterior tails; their incoming trace remains part of the
prescribed state. The conditional central field vanishes after the last
arrival through $t_C=0$.

The [adapter](../src/recursive_horizons/nsc_atlas_field_transfer.py) and
[focused tests](../tests/test_nsc_atlas_field_transfer.py) record the finite
support, waveform and checks. A quantum detector response, a remote source,
generic-field tails, persistent memory and stress-energy backreaction are
separate missing calculations.

## Recorded finite response

The proper-time pulse has $\tau_c/\ell=-1.1138558046$ and reaches the
wall over $R/\ell\in[0.4749900536,\,0.5228504193]$. Its central arrival
support is $t_C/\ell\in[-0.5059339731,\,-0.3479958418]$. The support
width is stretched by a factor $2.6323021882$ relative to the $0.06\ell$
wall proper-time input width. The normalized field $\ell\phi_C/\epsilon$
has extrema $+23.8414479985$ at $t_C/\ell=-0.4553839554$ and
$-24.4653347800$ at $t_C/\ell=-0.3958775671$. The negative lobe's magnitude
is $1.0261681581$ times the positive lobe's; the symmetric source was held
fixed in the wall clock.

The dimensionless squared-amplitude integrals are $96.9933769934$ for
$\ell\int I^2d\tau_W/\epsilon^2$ and $36.9002136729$ for
$\ell\int\phi_C^2dt_C/\epsilon^2$. Their ratio is $0.3804405498$.
These are waveform diagnostics, with no energy interpretation. The signed
output area divided by $\epsilon$ is numerically $-8.8\times10^{-16}$,
consistent with its exact zero. The central field vanishes after the final
arrival through $t_C=0$ for the stated initial-boundary data.

The initial null-label clearance is $2.7903396486\ell$. The first centre
arrival follows the wall's $R=0.8\ell$ event in conformal time by
$0.2446111895\ell$. Until the chosen endpoint, the reflected front has
comoving radius at most $0.6585338232\ell$, below the conservative wall
bound $0.8\ell$ by $0.1414661768\ell$. This checks the actual proper-time
pulse, rather than reusing a bound for a different profile.

The record uses 440 distinct wall radius samples and 620,000 evaluations
of the reused clock quadrature. The paired integral difference is below
$10^{-13}$; the independently integrated input differs from its exact
rational polynomial norm by $3.9\times10^{-13}$. A direct symmetric
radial limit of the composed primitive agrees with the central transfer
to relative error $2.8\times10^{-8}$ at $r/\ell=10^{-5}$. These are
implementation diagnostics, not rigorous quadrature error enclosures.
The 27 focused tests pass in 5.80 seconds, including repeated unit-scaled
records. The adapter also enforces its own finite sample and CPU budgets.

An independent mathematical review of the adapter and note found no
scientific blocker; its companion derivation is
[the scalar matching owner](nsc-atlas-field-matching.md). The present
connection closes an explicitly supplied-state response question. It leaves
the exterior origin of that state open and does not establish lasting
memory, energy transfer, a dark component or a generic property of other
fields or attachments.
