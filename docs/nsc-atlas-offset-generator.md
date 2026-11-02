# One spatial seed and a repeated offset: a constructive geometric map

## Active start: parent compression, not an isolated offset

Douglas's correction, 2 November 2026: the child's proposed generating
seed arises through collapse inside an evolving parent, itself inside an
enclosing region. The offset must be connected to that inherited drive.
The isolated-loop transfer and fixed-alpha calculations below are retained
controls. Their open/closed comparisons are not a task to select a preferred
percentage before constructing the parent-child relationship.

### An additive measure from the enclosing geometry

Let F_p(u), F_c(u) describe corresponding local material-patch deformations
in a common comparison frame. Their orientation-preserving volume factors
are J_p=det F_p and J_c=det F_c. The label u indexes the supplied history;
it is not an independently postulated physical clock. For spatial dimension
d define normalized volume factors j_a(u)=J_a(u)/J_a(u_0), and

$$
\chi_{p\to c}(u)=\frac{1}{d}\log\frac{j_p(u)}{j_c(u)}.
$$

For the constructed three-dimensional patches d=3. This is logarithmic
relative LENGTH compression. If the child contracts in every direction
by half while its parent patch is unchanged, chi=ln 2. Its normalization
comes from taking the d-th root of a volume factor, not fitting a target.
It requires nonzero-volume comparison patches; it is not itself a birth
solution at a degenerate one-dimensional seed.

The logarithm is forced if a continuous scalar measure of a positive
relative scale must add under composition: C(ab)=C(a)+C(b) gives C(a)=k ln a.
Taking the logarithm of the relative length ratio fixes its geometric unit.
This does not establish a conversion to angular phase, seconds or energy.
Identifying inherited compression with the proposed generating offset
remains the explicit physical proposal being developed.

The increment delta chi is calculated from the deformation history.
It is not an independently chosen alpha. Where derivatives exist,

$$
\frac{d\chi}{du}=\frac{1}{d}
\left[\operatorname{tr}(F_p^{-1}F'_p)
-\operatorname{tr}(F_c^{-1}F'_c)\right].
$$

A common invertible frame change F_a -> B(u)F_a cancels in the relative
ratio, including its initial normalization. For three nested regions,

$$
\chi_{g\to c}=\chi_{g\to p}+\chi_{p\to c}.
$$

Thus zero ADDITIONAL child compression does not erase the inherited parent
contribution. This is an exact finite-chain log-ratio telescoping identity
(bookkeeping, not dynamical evidence for the mechanism); it assumes neither
a first universe nor a common increment for every level. It is not a
derivation of an infinite physical chain.

### Executed controls and the next connection

Analytic controls F_g=exp(-u)I, F_p=exp(-2u)I, F_c=exp(-4u)I give
chi_gp=u, chi_pc=2u and chi_gc=3u, with zero sampled composition residual.
These are test deformations, not selected collapse histories. Independent
controls check common-frame cancellation, reverse sign under exchanging
regions, ancestor inheritance under common compression, and log-Jacobian
evaluation at compression 400 without determinant underflow. Pure rotation
has unit volume factor: rotation and compression cannot be identified by
volume data alone.

```text
python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --inheritance-only
```

**Next:** obtain the actual relative deformation, including turning/shear,
from a shared parent-child compression construction and connect its inherited
increment to the drawing operation. A pressure-to-deformation law, angular
phase relation and receiving geometry have not been derived by this
kinematic identity. No free increment is used to stand in for those links.

## Intrinsic deformation response and a conservative neck control

The next candidate uses the existing angular phase with a unit target
metric and the actual spatial metric:

$$
Q_g=\frac12\int g^{ij}\partial_i\psi\partial_j\psi\,dV_g.
$$

Choosing this functional and target metric is explicit. It is not the old
fixed-drawing-coordinate functional with a new name, and has no assigned
joule conversion. On a circle of radius R with phase m theta it gives
Q=pi m^2/R. Replacing theta by any smooth increasing periodic parameter
leaves this value unchanged when both the metric and volume factor are
transformed. The calculation checks this directly.

### How inherited compression changes the cost

For an affine deformation A of a reference patch, let
M=integral grad(psi) grad(psi)^T dV_reference. Then

$$
Q(A)=\tfrac12\det(A)\operatorname{tr}(A^{-T}MA^{-1}).
$$

The code represents M=C^T C, keeping cost positivity explicit. It does not
introduce another field. For uniform A=lambda I in dimension d,
Q(lambda)=lambda^(d-2) Q(1). Thus half-size costs relative to the original
are 2 in d=1, 1 in d=2, and 1/2 in d=3. These are scale controls, not
a selection of spatial dimension.

For the SAME carried reference phase pattern, this connects directly to
the enclosing-deformation result above:

$$
\frac{Q_c}{Q_p}=e^{-(d-2)\chi_{p\to c}}.
$$

This is a derived response to inherited compression, with no isolated alpha
or independently chosen damping rate. It is the phase-gradient portion of
a geometric budget, not a claim that total energy disappears.

### Work and directional response

Let Q_i be the contributions of the squared current phase gradient in
orthogonal directions. Changing a current principal length gives
partial Q/partial ln(lambda_i)=Q-2Q_i. These derivatives define the
generalized response; they are not inserted pressure labels. Under a
constant isotropic enclosing compressive load P>0, the load potential is
P V (equivalently, input work on compression is -P dV). In d=3,

$$
\frac{d}{d\lambda}[Q_0\lambda+P V_0\lambda^3]
=Q_0+3P V_0\lambda^2>0.
$$

Consequently this candidate alone has no finite isotropic holding radius
under that load. A one-dimensional circle's resistance cannot simply be
carried over to a three-dimensional region. The constant load is a symbolic
control; a physical enclosing load has not been obtained from this identity.

At FIXED volume, an isotropic M instead favours equal principal lengths
by the arithmetic/geometric-mean inequality. This is shape stability,
not stability of overall size. For A(s)=diag(exp(-s/2),exp(-s/2),exp(s)):

| Log axial stretch s | Q for normalized axial gradient | Q for M=I |
| --- | --- | --- |
| 0 | 0.5 | 1.5 |
| 1/2 | 0.1839397206 | 1.8326609913 |
| 1 | 0.0676676416 | 2.7859494701 |
| 2 | 0.0091578194 | 7.3982139184 |

An inherited directional phase pattern can therefore favour elongation and
transverse narrowing. These two reference patterns are analytic controls,
not selected preparations of a universe. A rigid rotation changes neither
cost; this functional alone supplies no torque that creates a vortex.

### Exact motion and its budget

For the normalized axial pattern Q(s)=exp(-2s)/2. Declare the geometric
inertial metric sum_i(d ln lambda_i)^2=(3/2)ds^2. From rest at s=0 its
conservative motion obeys

$$
\ddot s=\frac23e^{-2s},\qquad
s(\tau)=\log\cosh(\sqrt{2/3}\,\tau),\qquad
Q+\frac34\dot s^2=\frac12.
$$

This adds an explicit candidate inertia for deformation; it is not silently
deduced from the earlier phase kinetic rule. Offset tau is normalized and
has no calibration to the generating offset or seconds. Fixed volume makes
the isotropic load do zero work on this shape mode. The decreasing Q is
accounted for by increasing geometric motion:

| tau | Log axial stretch | Q | Geometric kinetic | Total |
| --- | --- | --- | --- | --- |
| 0 | 0 | 0.5 | 0 | 0.5 |
| 1 | 0.3017834121 | 0.2734288042 | 0.2265711958 | 0.5 |
| 2 | 0.9772952024 | 0.0708112364 | 0.4291887636 | 0.5 |
| 4 | 2.5742942146 | 0.0029037984 | 0.4970962016 | 0.5 |

### A0 least-stretch justification, enclosing compression minimization, and break continuation

Adopting the 3D volume-weighted functional $Q_g=\frac12\int|\nabla_g\psi|^2\,dV_g$
without deriving it from the generated geometry dilutes the 1D phase gradient
by the 3D volume factor $\det(A)\propto\lambda^3$, so $Q_g(\lambda)=\lambda Q_0=e^{-\chi}Q_0$
decreases to zero under 3D compression and offers no holding resistance.

In the offset generator, the geometry is drawn by 1D phase curves (spokes/windings).
By Cauchy-Schwarz on a 1D curve of physical length $L_i=\lambda_i L_{i,0}$ carrying
phase change $\Delta\psi_i$, the least-stretch cost with respect to induced arc
length $ds_g$ is uniquely
$$
\frac12\int_0^{L_i}\left(\frac{d\psi}{ds_g}\right)^2 ds_g
\ge \frac{(\Delta\psi_i)^2}{2\lambda_i L_{i,0}}=\frac{m_i}{\lambda_i},
\qquad
Q_{\rm stretch}(A)=\tfrac12\operatorname{tr}\!\left[(A^TA)^{-1/2}C^TC\right].
$$
Under uniform compression $A=\lambda I=e^{-\chi}I$, $Q_{\rm stretch}(\lambda I)=\lambda^{-1}Q_{\rm stretch}(I)=e^{+\chi}Q_{\rm stretch}(I)$
in every dimension $d\ge1$, resisting compression rather than vanishing.

Minimizing total cost $E(\lambda_1,\lambda_2,\lambda_3)=\sum_{i=1}^3\frac{m_i}{\lambda_i}+PV_0\lambda_1\lambda_2\lambda_3$
under enclosing compression $P>0$ gives the exact unique equilibrium for $m_i>0$:
$$
\lambda_i^*=\frac{m_i}{(m_1m_2m_3)^{1/4}(PV_0)^{1/4}},
\qquad
Q_{\rm stretch}^*=3PV_0\lambda_1^*\lambda_2^*\lambda_3^*=3W_{\rm ext}^*.
$$

With $PV_0=1$, the three generator regimes give:

1. **Isotropic 3D winding ($m=(1/2,1/2,1/2)$) $\to$ uniform contraction:**
   $\lambda_1^*=\lambda_2^*=\lambda_3^*=2^{-1/4}\approx0.8408964153$,
   $Q_{\rm stretch}^*\approx1.7838106765$, $W_{\rm ext}^*\approx0.5946035588$,
   $E_{\rm total}^*\approx2.3784142353$.
2. **Raw two-accumulator moment ($m=(1/4,1/4,1/2)$) $\to$ localized prolate neck equilibrium:**
   $\lambda_\perp^*=2^{-5/12}\approx0.7491535384$, $\lambda_3^*=2^{7/12}\approx1.4983070769$,
   aspect ratio $\lambda_3^*/\lambda_\perp^*=2$, $s_*=\tfrac13\log2\approx0.2310490602$.
3. **1D axial seed ($m=(0,0,1/2)$) $\to$ transverse neck-pinch break:**
   $\lambda_\perp^*\to0$ ($s\to+\infty$), because zero transverse winding cannot
   hold the transverse radius against enclosing compression.

**Break continuation into the child frame.** In parent coordinates, the axial
neck motion $s(\tau)=\log\cosh(\sqrt{2/3}\,\tau)$ drives the transverse throat
radius $a(\tau)=e^{-s(\tau)/2}\to0$ while the normalized velocity fraction
$\dot s/\sqrt{2/3}=\tanh(\sqrt{2/3}\,\tau)$ hits the $0.999\dots$ wall
($0.9998888764$ at $\tau=6$, where $a(6)\approx0.1220964547$). Applying the
break rule via throat inversion $R_{\rm child}(\tau)=a(\tau)^{-2}=e^{s(\tau)}=\cosh(\sqrt{2/3}\,\tau)$
continues the trajectory into an expanding child frame obeying the exact
Friedmann-like balance
$$
\left(\frac{\dot R_{\rm child}}{R_{\rm child}}\right)^2+\frac{2}{3R_{\rm child}^2}=\frac23,
\qquad
H_{\rm child}(\tau)\to\sqrt{\frac23}\approx0.8164965809,
\qquad
\Lambda_{\rm child}=3H_\infty^2=2.
$$
At $\tau=6$, $R_{\rm child}(6)\approx67.0801295585$ and $99.9778\%$ of the initial
parent phase cost has transferred into child expansion kinetic budget.

**Turning and spiral comparison.** Both $Q_{\rm vol}(A)$ and $Q_{\rm stretch}(A)$
are $SO(3)$-invariant functions of symmetric stretch and produce zero rotational
torque on their own. When unit logarithmic compression $\lambda(\chi)=e^{-\chi}$
is paired with a supplied or inherited angular turn rate $\omega=2\pi\alpha$,
$z(\chi)=e^{(-1+2\pi i\alpha)\chi}$ has decay-per-radian $\gamma/\omega=1/(2\pi\alpha)$
and opposite-lobe asymmetry $\eta=(A_+-|A_-|)/(A_++|A_-|)=\tanh(1/(4\alpha))$.
For $\alpha_{\rm golden}=(3-\sqrt5)/2\approx0.3819660113$, $\gamma/\omega\approx0.4166730505$
and $\eta\approx0.5746972340$ ($|A_-|/|A_+|\approx0.2700854214$), whereas a golden
logarithmic spiral has $\gamma/\omega=2\ln\varphi/\pi\approx0.3063489625$ (which
would require $\alpha=1/(4\ln\varphi)\approx0.5195217303$).

- **Coupling $c$:** Still supplied along the 1D drawing coordinate; the 3D neck shape wave speed $\sqrt{2/3}$ follows only after declaring the log-strain kinetic metric.
- **Increment $\alpha$:** Still supplied (or inherited from external rotational/chiral parent history) because symmetric compression costs are $SO(3)$-invariant and generate zero torque.

```text
python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --intrinsic-only
```

**Next connection:** derive a chiral/rotational coupling from the two-sided
(`+`/`0`/`-`) contact so turning and $\alpha$ emerge from the compression
response rather than remaining supplied.

## Retained foundation control: opposite windings, 2 November 2026

The geometry-only calculation is:

```text
python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --winding-only
```

It does not evaluate the supplied parent action, select an observational
target or write another result file. The older parent comparison below is
retained as a conditional construction, outside these premises.

### Closed orientations and their contact

Take a circle of radius R and a normalized drawing coordinate x in [0,1].
The two orientations are z_+(x)=R exp(2 pi i m x) and
z_-(x)=R exp(-2 pi i m x). Closure requires integer m. Their degrees are
+m and -m. Define the geometric squared variation, explicitly in this
drawing coordinate, by Q[z]=(1/2) integral |dz/dx|^2 dx. Each loop has
Q=2 pi^2 R^2 m^2: reversing orientation preserves its positive value.
This functional has no supplied conversion to joules or particle mass.

Join the loops at one contact and use a combined coordinate x in [0,2].
Let f(x)=x for x<=1 and f(x)=2-x for x>=1. The joined path

$$
Z_h(x)=R\exp[2\pi i m(1-h)f(x)],\qquad 0\leq h\leq1,
$$

is closed, continuous and never passes through zero. It is piecewise
smooth with finite squared variation at the join. At h=0 it follows the
two loops; at h=1 it is constant. During contraction the two halves need
not remain separately closed: allowing their common contact to reconnect
is an explicit operation, not a consequence of the offset update alone.

Its variation is exactly

$$
Q[Z_h]=4\pi^2R^2m^2(1-h)^2.
$$

Cauchy-Schwarz proves this is the least variation at each h for a fixed
circle radius, fixed endpoint lift and midpoint lift: on each length-one
half, integral (phase derivative)^2 >= (phase change)^2, with equality
for the constant phase slopes used above. This is a sharp constrained
minimum, not a claim that nature selects the constraints.

For R=m=1:

| Contraction h | Q |
| --- | --- |
| 0 | 39.4784176044 |
| 1/4 | 22.2066099025 |
| 1/2 | 9.8696044011 |
| 3/4 | 2.4674011003 |
| 1 | 0 |

This supplies a regular geometric continuation from opposite winding to
zero variation. The parameter h is not yet the offset's selected history.
No received energy, new region or outgoing signal is inferred from the
decrease without a transfer law.

### What forces the 1:1 and 180-degree pair?

Suppose exactly two outputs have positive scalar amounts a,b, unit
directions u,v, total a+b=Q_out, and balanced weighted directions
a u+b v=0, with no remaining receiver. Taking norms gives a=b, hence
u=-v. Therefore each amount is Q_out/2: **ratio 1, angle 180 degrees**.
This is an exact geometric consequence of those explicit output rules.

Balance alone does not select two outputs. Three equal outputs separated
by 120 degrees, or four separated by 90 degrees in a plane, also balance.
The calculation constructs all three cases. Nor does opposite orientation
select a contact rule: multiplying conjugate unit-circle maps gives a
constant, adding them gives 2 cos(phase), and retaining both keeps two
loops. The offset recurrence alone specifies none of these interactions.

**Next missing connection:** the contact composition/evolution law that
selects output channels and transfers the lost geometric variation.
Deriving that law, its coupling to the offset and an instrument readout is
what can distinguish this construction from a change of names. No 511 keV
scale, electron charge or matter excess has been derived here. Integer
winding is conditional on the circle target; it is not a charge assignment.

## Conservative contact evolution, 2 November 2026

The preceding contraction supplies a path but not its motion or a receiver.
Three distinct choices were considered: prescribing h directly, descending
Q by a geometric gradient, and giving the same circle phase conservative
inertial motion. The first prescribes the history; the second dissipates Q
without specifying a receiver. The third permits a closed transfer ledger
and an explicit open/closed comparison, so it is the worked candidate here.

### The added rule and its exact solution

Write Z=R exp(i psi). Its squared spatial variation is
Q=(R^2/2) integral psi_x^2 dx. **Add** geometric inertia using the same
circle metric: K=(R^2/2) integral psi_tau^2 dx. Stationarity of the
integral of K-Q, with fixed variation at the endpoints, gives

$$
\psi_{\tau\tau}=\psi_{xx}.
$$

This equation is derived from that newly declared local quadratic
inertial rule, not from the offset recurrence or from boundary minimization
alone. A ratio of inertia and spatial coefficients can be absorbed into
the dimensionless tau unit; no physical propagation speed or duration is
thereby predicted. No new field species, potential or observed scale is used.

Let g(x)=2 pi m f(x), with the triangular f above. Start with psi(x,0)=g(x)
and psi_tau(x,0)=0. Compare two drawing-coordinate continuations: extend g
by zero outside [0,2], or identify that interval periodically. These are
different domain assumptions for the SAME circle-valued target.
For either continuation the exact weak solution is

$$
\psi(x,\tau)=\tfrac12[g(x-\tau)+g(x+\tau)].
$$

The moving corners have finite energy and need no singularity cutoff.
Derivatives exist almost everywhere; integrals are evaluated separately
on their exact intervals. The two characteristics are derived from the
second-order equation. Equal amounts follow from zero initial velocity;
they are not selected by opposite winding alone.

### Where the variation goes

The local conserved density and flux are

$$
e=\frac{R^2}{2}(\psi_\tau^2+\psi_x^2),\qquad
j=-R^2\psi_\tau\psi_x,\qquad
\partial_\tau e+\partial_x j=0.
$$

Equivalently, e=e_right+e_left and j=e_right-e_left, with
e_right=(R^2/4)g'(x-tau)^2 and e_left=(R^2/4)g'(x+tau)^2.
Each integrates to 2 pi^2 R^2 m^2: exactly half the initial budget.
These are geometric conserved amounts, not calibrated joules.

For R=m=1, the initial total is 39.4784176044:

| Domain | Offset | Spatial Q | Inertial K | Total |
| --- | --- | --- | --- | --- |
| Open | 0 | 39.4784176044 | 0 | 39.4784176044 |
| Open | 1/2 | 9.8696044011 | 29.6088132033 | 39.4784176044 |
| Open | 1 | 19.7392088022 | 19.7392088022 | 39.4784176044 |
| Periodic | 1/2 | 0 | 39.4784176044 | 39.4784176044 |
| Periodic | 2 | 39.4784176044 | 0 | 39.4784176044 |

On the open line, the supports [tau,2+tau] and [-tau,2-tau] separate
at tau=1. Each outgoing packet carries 19.7392088022, with j=+e or j=-e.
Across each end of the original drawing interval [0,2], the outward
flux is pi^2 R^2 m^2 for 0<tau<2 and zero afterward, apart from corner
representatives that do not affect integrals. Thus the explicitly receiving
left and right intervals each acquire

$$
B_{\rm left}=B_{\rm right}=\pi^2R^2m^2\min(\tau,2),\qquad
B_{[0,2]}=4\pi^2R^2m^2-2B_{\rm left}.
$$

At tau=2 the original interval is empty of this geometric budget and each
receiving side holds exactly half. Direct boundary-flux and spatial
integration independently check these identities. These drawing-domain
boundaries are not being identified with the black-hole horizon (0).
Every fixed bounded interval eventually clears. The opposite propagation
directions have angle 180 degrees if read along a straight spatial axis;
the solution alone does not identify that axis or a three-dimensional
detector embedding.

On the periodic coordinate, psi is constant at tau=1/2 but its velocity
is not. Zero spatial variation is not zero total budget. The full initial
state recurs at tau=2: no irreversible drain or child has been produced.
This applies the continuation question concretely: changing the domain
changes where the conserved amount can go, rather than erasing it at a
flat drawing. Neither continuation has been derived to be an event horizon.

Reproduce without files or the conditional parent:

```text
python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --transport-only
```

**Next connection:** select the inertial contact rule, initial velocity and
receiving geometry from the offset construction. The present result is an
exact conditional transfer mechanism. Its 1:1 ratio and opposite directions
do not establish electron/positron identification or the 511 keV scale.

## Retained offset-history construction

Douglas's starting proposal is one spatial seed plus an offset operation,
applied repeatedly to generate a circle, a sphere and continued structure.
Here the offset is the primitive update; no independent proper-time
coordinate is used to generate the spatial map. The construction below
realizes part of that proposal with an explicit metric and a proof.

The [calculation](../results/summaries/nsc-atlas-offset-generator-v1.json)
uses one continuous seed coordinate and one irrational offset. A history
rule supplies the two phase values rather than choosing two independent
initial phases. Retaining the earlier draws, choosing a cross-pass distance
law and taking its metric completion are explicit assumptions.

## Preparation and choice

We compared supplying three independent spatial coordinates, claiming a
single smooth instantaneous trajectory fills a volume, and retaining an
offset history whose distance completion has higher dimension. The third
formulates the proposed drawing without making the false instantaneous-rank
claim. The earlier [suspension and golden-map controls](nsc-atlas-golden-vortex.md)
remain valid in their scope; this entry replaces independently supplied
phase coordinates with a single deterministic history.

The most likely mistakes are hiding a second freely supplied phase,
confusing a dense history with its completed space, and using an embedding's
component count as the intrinsic spatial dimension. The formula, completion
proof and induced-metric calculation check those separately. Finite coverage
uses fixed counts 100,200,1000,4000 and at most 120 CPU seconds; no physical
target selects a phase, count or distance rule.

## 1. The same accumulation operation is applied twice

Let alpha be irrational, theta_0=phi_0=0, and set, modulo one,

$$
\theta_{n+1}=\theta_n+\alpha,\qquad
\phi_{n+1}=\phi_n+\theta_n.
$$

Equivalently, Delta theta=alpha and Delta phi=theta. One fixed offset feeds
the first accumulator, whose history feeds the second. Thus

$$
\theta_n=n\alpha\quad\bmod1,\qquad
\phi_n=\frac{n(n-1)}{2}\alpha\quad\bmod1.
$$

The index n counts applications. It is not a separately supplied clock in
seconds. Both phases are determined by that count and alpha; there is no
second independent phase input. The inverse update exists, so the map
alone does not select a thermodynamic arrow or explain why updates occur.

**Density proof.** For any nonzero integer pair (k,l), the corresponding
Fourier character is the polynomial phase

$$
\alpha\left[kn+\frac{l n(n-1)}{2}\right].
$$

If l=0, its average tends to zero by the geometric-sum formula. If l is
nonzero, subtracting the value at n from the value at n+h gives a linear
phase with irrational slope l h alpha for every positive integer h.
Its average again tends to zero. The van der Corput bound then gives an
upper limit of order 1/H for the squared original average for every H;
letting H grow makes it zero. Weyl's criterion proves joint equidistribution
on the two-phase torus. The general polynomial theorem is attributed to
Weyl; [Tao's notes](https://terrytao.wordpress.com/2010/03/28/254b-notes-1-equidistribution-of-polynomial-sequences-in-torii/)
give the character/differencing framework, and
[Yifrach](https://arxiv.org/abs/2201.07138) records the classical statement
before developing a further generalization. No generalization beyond this
quadratic special case is needed here.

A rational offset need not work: alpha=1/4 returns to the initial state
after eight updates and supplies only finitely many directions. A zero
offset leaves the original line. Nonresonant accumulation is an essential
condition, not a word added after the calculation.

## 2. Line, circular boundary, spherical boundary

Put Theta_n=2pi theta_n and Phi_n=2pi phi_n. A first rotation of the retained
seed interval r in [-R,R] gives lines (r cos Theta_n,r sin Theta_n). Their
completion is a disk with a circular boundary.

For the second stage define the direction correlation

$$
C_{nm}=\cos\Theta_n\cos\Theta_m
+\sin\Theta_n\sin\Theta_m\cos(\Phi_n-\Phi_m),
$$

and the distance rule

$$
d^2[(r,n),(s,m)]=r^2+s^2-2rsC_{nm}.
$$

This is a declared metric kernel. Its Gram factorization is

$$
v_n=(\sin\Theta_n\cos\Phi_n,\sin\Theta_n\sin\Phi_n,
\cos\Theta_n),\qquad C_{nm}=v_n\cdot v_m.
$$

The distance is therefore Euclidean distance between r v_n and s v_m;
coincident labels, including every representation of the origin, are
identified. Density of the phase pairs implies density of directions on
the sphere. Every point of the ball is a limit of retained seed points.
Conversely every seed point lies in the ball, which is complete and closed.
**The metric completion is exactly the three-dimensional closed ball.**

This is a conditional constructive theorem: one spatial interval and one
offset history suffice for the stated completion. Finite histories and
the countable union of lines have zero three-dimensional volume before
completion. If motion were restricted to walking along the spokes through
their shared origin, a different distance would result. The distance law,
memory and completion must not be hidden inside the phrase "one plus time."

## 3. The same history supplies the closed parent geometry

For a radial seed r in [0,pi R], choose the positive-curvature distance law

$$
d=R\arccos\left[
\cos(r/R)\cos(s/R)+\sin(r/R)\sin(s/R)C_{nm}\right].
$$

Its completion is the round closed three-dimensional space used by the
parent model. A four-component representation is
(R cos(r/R), R sin(r/R) v_n), with one fixed-norm constraint. The intrinsic
metric, calculated by differentiation, is

$$
d\ell^2=dr^2+R^2\sin^2(r/R)
\left(d\Theta^2+\sin^2\Theta\,d\Phi^2\right).
$$

There are three intrinsic spatial directions. Four embedding components
do not make four independent spatial dimensions. Likewise a separate
five-coordinate auxiliary fluid/gravity model is not the count of primitive
variables in this generator.

The radius is still a parameter of this geometric construction. The
existing constant-field parent at its retained minimum supplies a compatible
choice through its Hamiltonian constraint. For zero extrinsic curvature,

$$
{}^{(3)}R=\frac{6}{R^2}=2\Lambda_*,\qquad
R=\sqrt{\frac{3}{\Lambda_*}}=3.684522144189
$$

in the current model units, with Lambda_star=.220983024213. Momentum
constraints vanish for those constant fields. This joins the generated
spatial geometry to existing admissible parent initial data. It is not a
primordial formation calculation or an SI scale selection.
The earlier action's phase ratio and radius are compatibility inputs;
this construction has not derived that action or its auxiliary
internal-dimensional interpretation from the offset alone.

The scalar Laplacian then has eigenvalues l(l+2)/R^2. The first two positive
values are .220983024213 and .589288064567. Under radius rescaling,
R^2(-Delta) preserves its dimensionless spectrum. This is an explicit
normalized operator inheritance law across separately specified regional
metrics; it does not yet attach those regions physically.

## 4. Finite resolution and measure

For N complete radial lines, let h be the spherical covering radius of
their two endpoint directions. The largest distance from the unit ball
to that line family is exactly sin(h). For N oriented north-to-south
meridians in the closed parent, it is min(h_oriented,pi/2) times R.
Spherical Voronoi vertices and a nearest-neighbor check evaluate these
finite errors; all phases are calculated at 80-digit precision before
converting the direction coordinates.

| Offset | Passes | Ball distance / R | Closed-parent distance / R |
| --- | --- | --- | --- |
| Retained action ratio | 100 | .35147195 | .50145491 |
| Retained action ratio | 1000 | .16691160 | .19978947 |
| Retained action ratio | 4000 | .08015006 | .12788791 |
| Golden-angle fraction | 100 | .39247396 | .57017792 |
| Golden-angle fraction | 1000 | .15218010 | .22441514 |
| Golden-angle fraction | 4000 | .08859494 | .13519849 |

The retained action offset is alpha=.0573586470504441, from the previously
derived ratio of linear mode frequencies. The golden comparison uses
alpha=(3-sqrt(5))/2=.381966011250105. Both are irrational; neither was fitted
to these errors. Their physical identification as the fundamental offset
is still a hypothesis. This construction does not select a cosmic pass
count, the 10^60 hierarchy, or a golden branching average.

Uniform phase passes are also not uniform physical volume. Their direction
second moments tend to (1/4,1/4,1/2). The metric area weight |sin Theta|
instead gives (1/3,1/3,1/3). The closed-parent volume element is
R^2 sin^2(r/R)|sin Theta| dr dTheta dPhi; the ball has r^2 in place of
R^2 sin^2(r/R). This is the metric measure used by the covariant action,
not an identification of update frequency with energy or particle density.

## What this adds to the atlas

The generator is now explicit and does not require an independent second
starting phase. Its completed spatial geometry and a compatible parent
slice are proved under the stated construction rules. Existing spacetime
calculations can be read as tests on a generated geometry, rather than
being asserted as the primitive origin of it.

What remains is to derive the physical selection of this distance kernel,
the offset law, two accumulation levels followed by nesting, and a duration
or irreversible arrow. Further accumulation levels can be constructed, so
exactly three is not selected by irrationality alone. A coherent physical
theory must also supply regional stress, formation and observable predictions.

## A conditional dimension-selection lemma

A more restrictive candidate can select three local directions. Let a
nontrivial real inner-product space V, dim(V)>=2, have a bilinear
antisymmetric bracket satisfying Jacobi and ad-invariance, with one c>0:

$$
\langle[x,u],v\rangle=-\langle u,[x,v]\rangle,
\qquad |[x,y]|^2=c^2(|x|^2|y|^2-\langle x,y\rangle^2).
$$

Then dim(V)=3. Ad-invariance and polarization give
ad_x²=-c²(|x|² Id-x tensor x). For orthonormal x,y define z=[x,y]/c.
The three vectors are orthonormal and satisfy [x,y]=cz,[y,z]=cx,[z,x]=cy.
Their orthogonal complement W is invariant under these adjoint maps.
On W put I=ad_x/c,J=ad_y/c,K=ad_z/c. The squared identity implies
I²=J²=K²=-Id and IJ+JI=0. Jacobi gives [I,J]=K, so K=2IJ and
K²=-4Id, contradicting K²=-Id unless W=0. This proves the claim.
Dimensions0 and1 are excluded because their zero brackets satisfy the
norm equation vacuously.

The quaternionic reference frame realizes this algebra with c=2r.
[Baez, Section4.1](https://arxiv.org/html/math/0105155v4) supplies the
related quaternion/octonion cross-product context; the seven-dimensional
norm-compatible cross product fails Jacobi. The short argument above
does not require a classification theorem.

These are candidate composition/isotropy postulates. The existing single
inheritance identity and offset generator do not derive them, and the
physical sourced metric need not be isotropic. The lemma selects neither
global topology nor radius, internal d=12 or an SI scale. The physical
dimension-selection gate therefore remains open.
