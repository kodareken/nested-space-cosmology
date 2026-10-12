# What an interior signal determines about its preparation

The [prepared response](nsc-atlas-prepared-response.md) connects independent
exterior data to a measured finite interior waveform. The next connection is
an inverse: **if this geometry is given, what exterior data are required by a
complete central field history?** That question must precede interpreting a
waveform as evidence for the geometry itself.

On the same retained wedge, compatible signed field histories determine a
unique preparation. The inverse is explicit as a convergent integral
resolvent, and it preserves an initially quiet interval. It does not preserve
the original pulse's late cutoff, positivity or polynomial shape. Those
restrictions give a genuine test of whether a proposed history belongs to
the allowed response family.

This is a conditional application of established contraction and resolvent
mathematics. No new waveform is evolved and no second geometry is shown to
have the same observations. The [admission](prior-art-and-open-claim.md#observer-preparation-inverse-admission-12-october-2026)
limits this cycle to proof, source comparison and evidence packaging.

## Fixed model and data class

Retain H ell=1, kappa ell=2, M/ell=0.4, the classical massless conformal
spherical scalar, no extra surface interaction, zero initial cap, and the
no-return wall interval 0.45<=R/ell<=0.55. Use the exact smooth wall and
clock, not a numerical interpolant as the definition of the geometry.
The [preparation certificate](nsc-atlas-field-preparation.md) supplies

$$
D=\{(U,V):V_a\le V\le V_b,\ U_0\le U\le w(V)\},
\qquad w(V)=U_W(V),\qquad U_a=w(V_a).
$$

Both wall null coordinates increase, R/ell>=0.4 everywhere in D, and

$$
K(U,V)=16M^4R^{-4}e^{-R/(2M)}>0,\qquad f_{UV}+Kf=0.
$$

Here f=chi/epsilon, g(V)=f(U_0,V) and h(V)=f(w(V),V). All three are
signed classical reduced fields. Set L(p)=V_a for p<=U_a and
L(p)=w^{-1}(p) otherwise. The inherited forward equation is

$$
f=Eg-Af,\qquad
(Af)(U,V)=\int_{U_0}^{U}dp\int_{L(p)}^Vdq\,K(p,q)f(p,q),
$$

where E extends a function of V constantly in U. The inherited analytic
bound is ||A||<=q_*<1 in the uniform norm, with

$$
q_*=16e^{-1/2}[\sqrt{F(0.4)}-\sqrt{F(0.55)}]^2,
\qquad F(x)=(1-x/0.8)e^{x/0.8}.
$$

The recorded q_* is approximately 0.138821; this cycle reuses that bound
without a new geometric or kernel evaluation.

Fix V_s=V_a+(V_b-V_a)/4, the original preparation's onset. Define the
enlarged preparation and wall-trace class

$$
\mathcal S=\{a\in C^2([V_a,V_b]):a=0\text{ on }[V_a,V_s]\}.
$$

Members need not stop at V_e=V_a+3(V_b-V_a)/4, remain positive, or keep the
original profile. The original C3 pulse belongs to this class, but the
class is larger. The formal small physical amplitude epsilon must still
be compatible with the linear-probe approximation.

## Inverse from a wall trace

Subtract the forward equation at w(V) from its value at U. This gives

$$
f=Eh+\mathsf Bf,\qquad
(\mathsf Bf)(U,V)=\int_U^{w(V)}dp\int_{L(p)}^Vdq\,K(p,q)f(p,q).
$$

Every integration point lies in D. The same rectangular majorant used for
A bounds ||B||<=q_*<1. The borrowed Neumann/contraction argument yields

$$
f=(I-\mathsf B)^{-1}Eh,\qquad
g=S(I-\mathsf B)^{-1}Eh,\qquad
\|g\|_\infty\le\frac{\|h\|_\infty}{1-q_*},
$$

where S takes the source trace at U_0. This is a convergent analytic
representation; no inverse iteration has been run.

The following checks establish that it is the inverse of the physical
finite problem, not just an abstract integral equation. Direct differentiation
gives

$$
f_U=-\int_{L(U)}^V K(U,q)f(U,q)dq,
\qquad
f_V=h'(V)+\int_U^{w(V)}K(p,V)f(p,V)dp.
$$

The moving upper-limit term in the second identity is zero because
L(w(V))=V. Thus f_UV=-Kf, f_U=0 on the wall, and f(w(V),V)=h(V).
At V_a the integral vanishes and h(V_a)=0 supplies the zero bottom.
Integrating f_U from U_0 recovers f=Eg-Af with g=Sf. Conversely every
forward solution gives the B equation by subtraction. Uniqueness proves
that the constructions are mutual inverses. Source and wall data have
not been independently imposed in the same forward problem.

Both operators depend only on V-prefixes. Uniqueness on each prefix
implies that an initially zero h gives f=g=0 on that entire initial
strip; the analogous statement holds in the forward direction.

### Regularity at the initial characteristic and wall corner

The lower limit L has a derivative switch at U_a. Its smoothness must not
be silently assumed. Put H=Kf. Away from the switch,

$$
f_{UU}=-\int_{L(U)}^V H_U(U,q)dq+L'(U)H(U,L(U)).
$$

The extra term has zero one-sided limit at U_a because
f(U_a,V_a)=h(V_a)=0. In the chosen initial-zero class, f is identically
zero near this corner, so the term is zero in a neighborhood of the
switch. The integral of H_U is continuous across it. The first-derivative
identities and smooth K give f_UV=-Kf and

$$
f_{VV}=h''+w'K(w(V),V)h(V)
+\int_U^{w(V)}(K_Vf+Kf_V)dp.
$$

Hence f has continuous second derivatives up to the retained boundary.
The same argument with the forward identities proves C2 regularity from
g in S. Its wall restriction is C2. This proves a bijection S<->S between
source and wall traces. Smooth joining to the zero cap is compatible
because the full initial open interval and its derivatives vanish.
No theorem for arbitrary incompatible corner data is claimed.
This regularity statement concerns the reduced exterior field and the C1
central trace below. It does not assert C2 regularity of the full
four-dimensional scalar at the centre for every C2 input: recovering that
scalar from the odd reduced radial field can lose a derivative. The
original C3 pulse and the smooth exclusion example have greater regularity.

## Inverse from the central observer

Use the known strictly increasing clock theta(V)=T(V)/ell and the signed
observable Y(theta)=ell phi_C/epsilon. The inherited central response is

$$
Y(\theta(V))=\frac{2h'(V)}{\theta'(V)},\qquad
h(V)=\frac12\int_{\theta_a}^{\theta(V)}Y(s)ds.
$$

The zero initial cap fixes the integration constant. Consequently the
preparation map is a bijection from S to

$$
\mathcal Y=\{Y\in C^1([\theta_a,\theta_b]):
Y=0\text{ on }[\theta_a,\theta(V_s)]\}.
$$

For two such histories within this same fixed model, linearity gives

$$
\|g_1-g_2\|_\infty\le
\frac{\theta_b-\theta_a}{2(1-q_*)}\|Y_1-Y_2\|_\infty.
$$

This bounds reconstructed source amplitude in response to uniform errors
in complete compatible field histories. It does not bound source derivatives,
frequency content, energy, backreaction or uncertainty in the geometry.
The forward map differentiates h, so this is not a claim of uniform-norm
stability in both directions. Intensity-only detections, finitely sampled
traces, missing early history and unknown clocks are different inverse
problems. No inference beyond the recorded finite clock interval follows.

## The preparation restriction is observable

Let P be the source class allowed by a physical preparation mechanism.
A proposed complete trace is allowed by this fixed model precisely when

$$
S(I-\mathsf B)^{-1}E
\left[\frac12\int_{\theta_a}^{\theta(V)}Y(s)ds\right]\in P.
$$

For example, keeping the original source cutoff requires this reconstructed
function to vanish throughout [V_e,V_b]. A single zero-area condition does
not enforce that functional restriction. Positivity and a fixed profile
are further restrictions. The original single normalized polynomial
preparation admits only its own response (or its amplitude multiples if
that amplitude is left free).

There is also a direct exclusion for the larger class of nonzero
nonnegative sources supported in [V_s,V_e]. The positive-kernel argument
already given in the [response owner](nsc-atlas-prepared-response.md#an-independent-sign-and-magnitude-check-after-the-input)
uses positivity and support, not the particular polynomial. It therefore
gives, after the source ends,

$$
h(V)<0,\qquad h'(V)>0,\qquad
Y(\theta(V))>0,\qquad
\int_{\theta_a}^{\theta(V)}Y(s)ds=2h(V)<0.
$$

These are joint constraints on the finite signed trace, not on energy.
The strict signs follow because the continuous nonzero source is positive
on an interval and the kernel and retained transverse interval are
positive. They do not extend past V_b.

To see that the enlarged response class is genuinely larger, choose a
smooth, nonzero, nonnegative Y supported strictly between theta(V_s) and
theta(V_e). Its integrated h is nonnegative and becomes a positive
constant. Positivity of B makes f=(I-B)^(-1)Eh>=Eh and therefore g>=h>0
after V_e. This trace has a unique initially quiet preparation in S,
but that preparation cannot switch off at the prescribed cutoff. This is
an exact within-model counterexample to unrestricted fitting with a
compact source, not a calculated alternative universe or a new wave run.

## What this unlocks for the atlas

Three questions remain distinct:

1. **Fixed geometry, unknown preparation:** the complete compatible signed
   field trace has the unique preparation just derived.
2. **Different geometries, different permitted preparations:** a trace is
   degenerate only if it lies in both restricted response ranges. No second
   compatible geometry is supplied by this result.
3. **Different geometries, the same independently identified preparation:**
   comparing response operators is a stronger experiment. Identifying the
   physical input across regions is part of that experiment.

The next useful connection is an independently justified preparation
restriction and its response-range test, or an additional field channel
whose preparation is physically tied to the first. A second detector in
the unchanged cap need not help: a fixed full wall history and zero cap
already determine the interior solution. An extra observable must have
an independent dependence before it earns another calculation.

This inverse is a reconstruction by the analyst; it supplies no remote
source mechanism and creates no causal route from the excluded right
parent. Formation, quantum-state selection, descendants, connected mouths
and exterior-response alternatives retain their separate open connections.

## Borrowing, independent review and evidence

[Brunner's *Volterra Integral Equations*, Section 1.2.1](https://assets.cambridge.org/97811070/98725/excerpt/9781107098725_excerpt.pdf)
supplies the standard Neumann/resolvent method (equations 1.2.4-1.2.10,
Theorem 1.2.3 and Remark 1.2.4). Its one-variable triangular-kernel theorem
is not imported as a ready-made thin-wall theorem: the retained-domain
norm and differentiated identities above supply this application.

An independent analytic review recovered the inverse, PDE and wall
conditions, source/bottom compatibility, the corner correction and central
normalization. The forward-map cross-check is h=(I-Q)g with
Q=RA(I+A)^(-1)E, where R is wall restriction. Here ||Q||<=q_*/(1-q_*)<1,
so its resolvent gives the same inverse by uniqueness. Its sufficient
q_*<1/2 condition is stronger than the B proof's q_*<1 and is not a
necessary physical bound. The review checks the C2 statement used here;
no computation is substituted for the proof.

The [literature library](arxiv-studies.md#observer-inverse-library)
contrasts arbitrary bulk-forcing nullspaces, controlled nonlinear
source-response maps, and recorded stochastic boundary inputs with this
one passive finite trace. Those results preclude a universal claim that
geometry can never be inferred. No full proof audit of those external
inverse theorems or exhaustive priority search has been performed.

The creation-only [evidence record](../results/development/nsc-atlas-observer-inference-v1/inference.json)
binds this note to the original v7/v8 inputs and explicitly distinguishes
an analytic result from numerical reconstruction. Its
[manifest](../results/development/nsc-atlas-observer-inference-v1/manifest.json)
and [checker](../scripts/record_nsc_atlas_observer_inference.py) authenticate
the declarations and prerequisite bytes. Passing that checker authenticates
the record; the mathematical argument above supports the theorem.
