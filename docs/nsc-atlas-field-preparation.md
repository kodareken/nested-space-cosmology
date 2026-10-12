# A finite exterior preparation and its interior transfer operator

This connection replaces an arbitrary wall trace by independently specified
classical data on a finite null surface in the **retained Schwarzschild
white-hole sector IV**. Its future reaches the early moving wall and the
regular interior centre. A bounded integral operator connects the data to
the wall trace; this record does not numerically solve a new waveform.

The source surface is outside the wall. It is not the designated right
parent exterior. The existing reference excludes a causal signal from that
right parent to the wall or centre at every parent emission time.

## Admission and borrowed mathematics

- **Borrowed:** the fixed v5 signed wall and regular Kruskal charts; the v6
  conformal scalar matching law and regular-centre solution; the Schwarzschild
  scalar radial potential from [Dolan and Ottewill, 1106.4318v1, Sections II–III,
  especially equation (11)](https://arxiv.org/html/1106.4318v1); ordinary
  characteristic integration and the contraction mapping argument.
- **Reading boundary:** [Bär and Tagne Wafo, 1408.4995v2, Section 4.2,
  Theorem 23](https://arxiv.org/html/1408.4995v2), supplies characteristic
  initial-value context. Its smooth wave-operator and past-compactness
  hypotheses matter. It is not applied as a theorem for our thin wall.
  The explicit finite integral argument below supplies this certificate.
- **Missing:** v6 specified an incoming characteristic at the wall without
  an exterior preparation that determines it.
- **Decision and unlock:** certify one admissible exterior data surface,
  its retained domain and the operator carrying its data to the interior.
  This makes a later source-fixed response comparison possible without
  assuming that a source in the ordinary right parent can drive the signal.
- **Fixed model:** four-dimensional massless classical scalar, spherical
  mode zero, coupling xi=1/6, linear probe on the expanding reference with
  H ell=1, k ell=2 and M/ell=0.4. No additional surface scalar interaction.
  The preparation changes the declared input; it does not retune the theory.
- **Finite domain and budget:** wall radii 0.45<=R/ell<=0.55, exterior
  R/ell>=0.4, at most 513 boundary samples and 60 CPU seconds. Use existing
  32/64-order geometric quadratures. No PDE evolution, operator iteration,
  parameter scan or new detector calculation.
- **Stop or redirect:** reject an initial surface on the excised side,
  incompatible corner/cap data, a lost no-return condition, a changed
  coupling, or a domain for which this contraction certificate fails.
  No extension to a later observer time is inferred from truncating the record.

## A preparation inside the retained vacuum

Put x=R/ell and m=M/ell=0.4. In the recorded Kruskal normalization,

$$
UV=F(x),\qquad F(x)=(1-x/0.8)e^{x/0.8}.
$$

The future timelike wall has increasing U and V. Throughout the chosen
window it lies in IV, where both are negative. With the normal directed
from the cap into the retained vacuum, its null components are
(-Udot,Vdot). Therefore the retained side at a fixed V is U<=U_W(V).
Being in quadrant IV alone would not establish this membership.

Choose x_a=0.45, x_b=0.55 and the radial floor x_f=0.4. Define

$$
V_a=V_W(x_a),\qquad V_b=V_W(x_b),\qquad
U_0=F(x_f)/V_a.
$$

The finite domain is

$$
D=\{(U,V):V_a\le V\le V_b,\quad U_0\le U\le U_W(V)\}.
$$

It is a characteristic strip truncated by a timelike wall. A rectangular
past diamond continued past the wall would include the cap and would use
the wrong exterior equation. In D, negativity and wall monotonicity give

$$
F(0.55)\le UV\le F(0.4).
$$

Since F is strictly decreasing for positive x, every point has
0.4<=R/ell<=0.55. The strip avoids both R=0 and the Schwarzschild horizon
R/ell=0.8. Its initial source leg U=U_0 lies on the retained side, and each
constant-V ray reaches the wall by increasing U entirely within D.

The recorded coordinates are

| Boundary | Dimensionless coordinate |
| --- | ---: |
| U_0 | -0.32130105902052203 |
| U_a | -0.299270297350313 |
| U_b | -0.274746803092629 |
| V_a | -2.565695357068247 |
| V_b | -2.262011613052105 |

## Independent data and the initial cap

For chi=R phi, prescribe chi=0 on V=V_a from U_0 to U_a. On U=U_0 use

$$
g(V)=\epsilon B\left(\frac{4(V-V_{\rm mid})}{V_b-V_a}\right),
\qquad V_{\rm mid}=\frac{V_a+V_b}{2}.
$$

B(z)=(1-z^2)^4 when |z|<1, and B=0 otherwise. This is a compact **C3
polynomial**, with a C2 first derivative; it is not infinitely smooth.
Its support is the middle half of the V interval. The source corner and
an open initial interval are zero. Epsilon is a positive formal small
amplitude used for normalization, not a fitted density or energy fraction.
The profile is specified on the exterior null surface independently of
the desired wall trace or central waveform.

Supply zero field and zero flat-time derivative on the cap at
t_a=t_W(0.45), with t_a/ell=-1.179595952350965. This is earlier than the
v4 slice at t=0; the latter cannot initialize this early pulse. These are
explicit data for the new finite problem. They are consistent with keeping
v6's original zero cap at R/ell=0.05 **if** no incoming drive is additionally
assumed before this preparation. Zero data at 0.05 alone would not imply
zero data at t_a in the presence of an unspecified earlier signal.
Neither this conditional extension nor the null input derives a global
past state or remote physical source.

## Why the outgoing wall characteristic is zero

In the conformal cap put u=eta-r and v=eta+r. The zero cap data make the
generating function constant over the initial cap's full null-label
interval. Subsequent incoming data have v>=v_a. The wall labels satisfy

$$
v_a/\ell=-1.789182750439567,\qquad
u_b/\ell=-4.423774662616959.
$$

Both wall labels increase, while u_b<v_a. Thus no incoming signal can
reflect at the regular centre and return to this wall window. The outgoing
cap characteristic is zero over the entire window, rather than just at
sampled events. The numerical label margin is 2.634591912177392 ell.

There is also a clock-independent analytic inequality. On the selected
window beta_- decreases and exceeds one. The inherited wall proof gives
Rdot^2>1/8, while Rdot^2<beta_-^2. Consequently

$$
\Delta t/\ell\le(b-a)\sqrt8[1+\beta_-(a)^2]
<\log\frac{1+b}{1-a},\qquad a=0.45,\ b=0.55.
$$

The first inequality follows from
gamma=(1+Rdot^2)/(beta_-+x Rdot), and the second is a direct bound using
the fixed endpoint a. This proves u_b<v_a without treating numerical
quadrature agreement as a rigorous error enclosure.

Conformal matching has continuous chi and continuous normal derivative
of chi. The outward normal is n=(-Udot,Vdot), so

$$
(\partial_\tau-n)\chi=2\dot U\,\partial_U\chi=0
$$

on the exterior wall limit. This condition is inherited from the prepared
cap and absence of a return; it is not an independently adjustable wall law.

## The exterior-to-wall operator

The Schwarzschild scalar curvature vanishes. Its spherical scalar potential
f(R)2M/R^3 transforms to the regular null equation

$$
\partial_U\partial_V\chi+K(R)\chi=0,\qquad
K(R)=\frac{16M^4}{R^4}e^{-R/(2M)}=-\frac{R_{UV}}{R}.
$$

Implicit differentiation of UV=F(R/ell) supplies the last identity
independently of transforming the tortoise-coordinate wave equation.
K is positive and decreases with R on the certified interval.

For an integration variable p, let L(p)=V_a when p<=U_a, and let L(p)
be the inverse wall coordinate V_W(p) otherwise. Integrating from the
zero bottom derivative or zero outgoing wall derivative gives

$$
\chi(U,V)=g(V)-\int_{U_0}^{U}dp\int_{L(p)}^{V}dq\,
K(p,q)\chi(p,q).
$$

The lower limit follows the actual domain. It supplies the bottom data
below U_a and the wall condition above U_a. Differentiation recovers
chi_UV=-K chi, the source trace and chi_U=0 on the wall.

The integral operator has sup-norm bound K_max DeltaU DeltaV, where
K_max=16 exp(-1/2). For F_f=F(0.4), F_b=F(0.55) and r=V_b/V_a>0,

$$
\Delta U\Delta V=F_f+F_b-F_fr-F_b/r
\le(\sqrt{F_f}-\sqrt{F_b})^2.
$$

The inequality is AM-GM. Hence an analytic bound, independent of the
integrated wall coordinates, is

$$
q_* =16e^{-1/2}(\sqrt{F(0.4)}-\sqrt{F(0.55)})^2
=0.13882099234357653<1.
$$

The standard contraction argument establishes a unique continuous solution
on this finite domain. The integral identities give the required first and
mixed derivatives for the compatible data. This is an application of known
linear-wave and fixed-point mathematics, not a new general existence theorem.
No numerical operator iteration has been performed.

It also establishes the reduced-field bound

$$
\|\chi-g\|_\infty\le\frac{q_*}{1-q_*}\|g\|_\infty.
$$

For the stated positive normalized pulse, the wall value at V_mid is at
least about 0.8388 epsilon. This is a bound on chi, not a central phi peak,
an energy flux or a persistent detector record. The new central waveform
has not been evaluated.

## Interface to the existing interior response

The trace needed by v6 is fixed by the admitted operator:

$$
I(\tau)=2\dot V\left[g'(V)-\int_{U_0}^{U_W(V)}K(p,V)\chi(p,V)dp\right].
$$

This is exactly I=(partial_tau+n)chi in the existing convention. In the
no-return window chi_U=0, so I=2 dchi_W/dtau. The cap generating function
satisfies F(v_W)=chi_W, giving the same regular-centre normalization as v6:

$$
\phi_C[T(\tau)]=\frac{I(\tau)}{T'(\tau)},\qquad
T=t_{\rm flat}-\ell\log(1-R/\ell).
$$

The admitted central proper-time window is
-0.581758951595344<=t_C/ell<=-0.250230700560700. Its endpoints are
the images of the wall window, not a claim that the received pulse has
exactly that support. Smooth exterior scattering can add a tail. Data and
propagation after this window remain uncomputed, including the interval
from its upper endpoint to t_C=0. The previous compact wall-input result
cannot set this new tail to zero.

## Evidence and implementation

The [adapter](../src/recursive_horizons/nsc_atlas_field_preparation.py)
exports preparation_record(length_scale=1), a deterministic JSON-compatible
record with 65 wall samples, 65 source samples and 33 bottom samples.
The [tests](../tests/test_nsc_atlas_field_preparation.py) check full retained
membership, independent mixed radius derivatives and the scalar potential,
boundary compatibility, theory changes, unit rescaling, clock exits,
no-return and the analytic operator bound.

The 32/64-order geometric comparison differs by at most about 2.87e-11 in
the checked dimensionless null coordinates and wall clocks. This is an
integration diagnostic. The analytic contraction and no-return proofs
above do not rely on calling this difference a rigorous error bound.

The right-parent control remains excluded: P has V>0 and the wall V<0,
whereas a future causal curve cannot decrease V. This is a statement about
driven signalling in the isolated reference, not initial correlations,
regular formation, the glued many-body future or arbitrary atlas branches.
The source of the admitted signal remains the declared IV preparation.
