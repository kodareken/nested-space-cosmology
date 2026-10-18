# Unlimited inherited levels and a finite local response

This extends the [finite inherited-law construction](nsc-nested-qualities.md)
to an actual infinite branching operator. The result establishes existence,
exact composition, a controlled response at the root, and a spectral scaling
dimension. Its assumptions and proof are explicit below. The construction
uses the same local matrix and link at every level, rather than adding a
new correction separately for each eliminated region.

## Theorem

Let the hierarchy be a rooted tree with b>=2 children per vertex. Fix
Omega>1, a positive Hermitian d-by-d matrix H and a d-by-d link B. Put

\[
h=\lambda_{\min}(H),\quad c=\sqrt{b/\Omega}\,\|B\|,
\qquad h>2c,
\]
\[
\delta=h-2c>0,\quad M=\lambda_{\max}(H)+2c,\quad
\gamma=\frac{h+\sqrt{h^2-4c^2}}2,\quad r=\frac{c^2}{\gamma^2}<1.
\]

The Hilbert space is the square-summable regional fields
`l2(vertices; C^d)`. For a vertex v at depth n, the diagonal block is
Omega^n H. Its block to each child is Omega^n B; the reverse block is
the Hermitian adjoint. These are exactly the existing finite-window
scalings, now with branching. At infinite depth select the finite-energy
form domain

\[
\mathcal Q=\left\{u:\sum_v\Omega^{|v|}\|u_v\|^2<\infty\right\}.
\]

Then the following hold.

1. The corresponding closed positive form defines one specified
   self-adjoint operator D with form domain Q, compact inverse and
   D>=delta I. This is the finite-energy/Friedrichs realization; no
   uniqueness among other possible boundary conditions at infinity is asserted.
2. Every isolated descendant subtree has the same normalized law:
   `U* D_subtree U = Omega D`. Thus its normalized pullback is exactly
   `T* D_subtree = D`, for arbitrarily many generations.
3. For lambda>=0, the root response
   `R(lambda)=P_root(D+lambda I)^(-1)P_root*` satisfies

   \[
   R(\lambda)=\left[H+\lambda I-
   \frac{b}{\Omega} B R(\lambda/\Omega)B^\dagger\right]^{-1}.
   \tag{1}
   \]

   This equation has a unique solution in the class
   `0<=R(lambda)<=(gamma+lambda)^(-1) I`.
4. Let R_N retain depths0,...,N, with a Dirichlet cut beyond the last
   generation. Then

   \[
   0\le R_N(\lambda)\le R(\lambda),\qquad
   \|R(\lambda)-R_N(\lambda)\|\le\frac{r^{N+1}}\gamma.
   \tag{2}
   \]

   The bound is uniform on lambda>=0. If B is invertible, adding each
   generation strictly increases the response as a positive matrix.
5. If N_D(E) counts the eigenvalues of D not exceeding E, including
   multiplicity, then

   \[
   \lim_{E\to\infty}\frac{\log N_D(E)}{\log E}
   =\frac{\log b}{\log\Omega}.
   \tag{3}
   \]

   For a second-order spectral interpretation, set d_s to twice this
   exponent. Then the regional length ratio is Omega^(-1/2) and
   `d_s=2*log(b)/log(Omega)=log(b)/log(L_parent/L_child)`.
   A first-order operator instead has dimension equal to the counting
   exponent. Neither identification is a claim about the local
   dimension of spacetime or the Hausdorff dimension of an unspecified embedding.

## Proof

**Existence and the infinite domain.** Define the diagonal operator
`(Su)_v=Omega^(|v|/2)u_v`. Let C have diagonal H and block
`Omega^(-1/2) B` from a vertex to each child. The one-way child-sum
operator has norm at most c: Cauchy--Schwarz gives
`||sum_children u_w||^2<=b sum_children ||u_w||^2`, and every nonroot
vertex has one parent. Its adjoint has the same norm. Therefore

\[
\delta I\le C\le M I.
\]

For finite-support u, the desired form is exactly
`q[u]=<Su,C Su>`. This expression extends to Q=Dom(S), with
`delta ||Su||^2<=q[u]<=M ||Su||^2`. Since S is closed, the form norm is
equivalent to its complete graph norm. Hence q is closed and positive.
Its self-adjoint operator has the explicit inverse

\[
D^{-1}=S^{-1}C^{-1}S^{-1}.
\tag{4}
\]

Precisely, `Dom(D)={u in Q: C S u in Dom(S)}` and `D u=S C S u` there.
The operator domain is not identified with the larger form domain Q.

Indeed, for u equal to the right-hand side applied to f,
`<Sv,C Su>=<v,f>` for every v in Q. The finite-support vectors are a
form core by truncating the convergent weighted sum. Each depth contains
only d*b^n coordinates, while the diagonal entries of S^(-1) tend to
zero. Thus S^(-1), and consequently (4), are compact. This also gives
the discrete spectrum and the stated lower bound.

**Inheritance and composition.** Removing the parent link and relabeling
one child as the new root multiplies every block by Omega, including
every link. It maps the stated form domain onto itself under the same
constant normalization. Equality therefore holds for the closed forms
and their operators, not merely for finite matrices.

Schur elimination at the root sees b identical child subtrees. Each has
operator Omega D, whose root resolvent is
`Omega^(-1) R(lambda/Omega)`. Summing the b self-energies gives (1).
Equivalently apply finite Schur elimination and take the form-core limit.
For that limit one can use the variational identity

\[
\langle f,(D+\lambda I)^{-1}f\rangle
=\sup_{u\in\mathcal Q}\{2\operatorname{Re}\langle f,u\rangle
-q[u]-\lambda\|u\|^2\}.
\]

Restricting the supremum to successive finite depths increases it to the
same supremum, because finite-support vectors are a form core. In the
finite-dimensional root block this proves Loewner and norm convergence.

**Finite response and explicit depth error.** Start at R_(-1)=0. Finite
Schur elimination is repeated application of the map in (1). The identity
`h-c^2/gamma=gamma` implies that if the incoming response is at most
`gamma^(-1)I`, then the next denominator is at least
`(gamma+lambda)I`. All finite responses, and their limit, have the
claimed bound.

For two incoming responses X,Y in this bounded positive class, the
resolvent identity gives

\[
\|F_\lambda(X)-F_\lambda(Y)\|
\le\frac{c^2}{(\gamma+\lambda)^2}\|X-Y\|
\le r\|X-Y\|.
\]

The same estimate applies to bounded response functions with their
supremum norm on lambda>=0, because argument scaling preserves this
half-line. It proves uniqueness. Apply it N+1 times to the infinite
response and the terminal value0 to obtain (2). Order preservation of
matrix inversion and positive congruence gives monotonicity. With B
invertible, each initially positive response increment remains positive
definite under these operations, proving strict influence of every depth.

**Spectral scaling.** The form inequalities are
`delta S^2<=D<=M S^2`. The min--max characterization of ordered
eigenvalues therefore bounds the counting function by

\[
N_{S^2}(E/M)\le N_D(E)\le N_{S^2}(E/\delta).
\]

The eigenvalue Omega^n of S^2 has multiplicity d*b^n, so, for x>=1,

\[
N_{S^2}(x)=d\sum_{n=0}^{\lfloor\log x/\log\Omega\rfloor}b^n.
\]

The geometric sum is bounded above and below by positive constants times
`x^(log(b)/log(Omega))`. Taking logarithms and the limit proves (3).
Log-periodic variations in a prefactor are allowed; no nonexistent
constant prefactor limit is assumed. This completes the proof.

## Corollary: an explicit nested geometry with the same dimension

The hierarchy also has a specified metric realization. Let X be the set
of infinite addresses with symbols1,...,b. If two addresses agree for
exactly k initial symbols, set

\[
d(x,y)=q^k,\qquad q=\Omega^{-1/2},\qquad d(x,x)=0.
\]

The common-prefix inequality proves the stronger triangle inequality
`d(x,z)<=max(d(x,y),d(y,z))`, so this is an ultrametric. A depth-n cylinder
is one fixed prefix; it is a ball of diameter q^n containing b subballs
of diameter q^(n+1). Prefixing any symbol scales every distance by q.
The normalized geometric law is therefore identical at every depth.

This space is compact: at each depth it has a finite cover by those
balls, and a Cauchy sequence fixes every successive prefix and hence
has a limiting infinite address. Its Hausdorff dimension is

\[
d_H=\frac{\log b}{\log(1/q)}=\frac{2\log b}{\log\Omega}=d_s
\]

for the second-order operator interpretation above.

**Proof.** The depth-n cover has b^n sets of diameter q^n. For s>d_H,
its total s-power diameter is `(b*q^s)^n`, tending to zero, giving the
upper bound. Give each cylinder the uniform measure b^(-n). A set of
diameter r with `q^(n+1)<r<=q^n` lies in a depth-n cylinder and has
measure at most b^(-n)<=q^(-d_H)*r^(d_H). Consequently every cover of
X has sum of d_H-power diameters at least q^(d_H), proving the lower
bound. Singleton sets have zero measure and do not affect the argument.
This proves the equality and the dimension.

Thus the noninteger dimension can belong to an actual defined nested
metric, as well as to its inherited operator's spectrum. This ultrametric
is a mathematical geometry of addresses and resolution; it is not assumed
to be a Lorentzian metric for the cap--wall spacetime.

## Corollary: deriving the matrices from edge geometry

The inherited H and B can themselves be derived from one specified
geometric action. Give every parent--child edge at depth n length
`L*q^n`, q in(0,1). Give a vertex half the length of each incident edge
as its kinetic weight. The root has a fixed parent port with edge
length L/q. Thus `m_n=m0*q^n`, where `m0=L*(1/q+b)/2`.

For real vertex fields phi, define the weighted network action

\[
\mathcal L=\tfrac12\sum_v m_{|v|}\dot\phi_v^2
-\tfrac12\sum_{v\to w}\frac{(\phi_v-\phi_w)^2}{Lq^{|v|}}
-\frac{q}{2L}\phi_{\rm root}^2.
\]

The edge potential is the integral of the squared gradient of the
linear field along that edge. The kinetic rule is the stated half-edge
(lumped) rule; it is not claimed to be the exact continuum edge mass
matrix. Close the action on finite-support fields as in the theorem.

In canonical coordinates `u_v=sqrt(m_n)*phi_v`, direct expansion gives

\[
\Omega=q^{-2},\qquad H=\frac{b+q}{Lm_0},\qquad
B=-\frac1{Lm_0\sqrt q}.
\]

Consequently `c=sqrt(b*q)/(L*m0)`, and the required positive margin is

\[
h-2c=\frac{(\sqrt b-\sqrt q)^2}{Lm_0}>0.
\]

It follows from geometry for every admitted b and q. No independent
H or B is selected to make the theorem apply. The exact static result is

\[
\gamma=\frac b{Lm_0},\qquad r=\frac qb,\qquad
\frac{\phi_{\rm root}}{f_{\rm root}}\bigg|_{\rm static}=\frac Lb.
\]

For the last equality a force is coupled as `f_root*phi_root` in the
action: its compliance is `R(0)/m0=1/(gamma*m0)`. This is a local
measurement predicted by this explicit geometry and action, including
all descendant levels and the declared parent port. The boundary choice
at infinity is the specified form closure. Other boundary choices are
not silently identified with it. This is a geometric network realization,
not an Einstein cap--wall solution or a generated cosmological Lambda.

## Metric variation: a spatial volume and its pressure

The preceding length-weighted network has one-dimensional edge measure.
To obtain a three-dimensional volume response, its edges must carry a
section as well as a length. Give an edge at depth n length `L*q^n` and
section `A*q^(2*n)`, with A>0. Its conductance is `(A/L)*q^n` and its
volume is `A*L*q^(3*n)`. The same half-edge rule gives

\[
m_n=m_0q^{3n},\quad m_0=\frac{AL}{2}(q^{-3}+b),\qquad
H=\frac{A(b+q^{-1})}{Lm_0},\quad B=-\frac{A}{Lm_0q^{3/2}}.
\]

The operator still inherits `Omega=q^(-2)`. Its positive form margin is
`A*(sqrt(b)-1/sqrt(q))^2/(L*m0)`, so the theorem applies when b*q!=1.
This calculation does not decide the critical case. The bounded branch gives

\[
\gamma=\frac{A}{Lm_0}\max(b,q^{-1}),\qquad
\chi_0=\frac{R(0)}{m_0}=\frac{L}{A\max(b,q^{-1})}.
\tag{5}
\]

These identities follow by expanding the geometric action and solving its
scalar Schur quadratic. They apply on both sides of b*q=1. If all edge
volumes are regarded as disjoint physical pieces, their sum is finite only
when b*q^3<1; that requirement is separate from the operator's positive
margin. This construction has not embedded those pieces in the actual
Lorentzian atlas or supplied the parent gravitational action.

**Dilation proposition.** Fix a positive kinetic matrix M0 and a closed
positive spatial form K0, using the inherited operator for their canonical
combination. Under a common physical spatial dilation by a>0,

\[
M(a)=a^3M_0,\qquad K(a)=aK_0,\qquad
D(a)=M(a)^{-1/2}K(a)M(a)^{-1/2}=a^{-2}D_0.
\]

The physical static root compliance therefore scales as
`chi(a)=chi0/a`. It scales oppositely to the original length-only action's
`a*chi0`. The same normalized operator scaling alone cannot distinguish
these different physical measures. For the already-owned spatial metric
`h=r^2(Q^2 dx^2+dOmega^2)`, direct variation r -> a*r at fixed Q similarly
gives `M_ij=integral r^3 Q phi_i phi_j dx` proportional to a^3 and
`K_ij=integral (r/Q) phi_i' phi_j' dx` proportional to a. Thus this spatial
measure is the one needed for a kinetic-plus-gradient lift of that radial
operator. It does not replace the atlas's improved conformal-scalar action.

Include a homogeneous lapse N and use real fields. The declared lifted action is

\[
S_\phi=\frac12\int dt\left[
\frac{a^3}{N}\dot\phi^TM_0\dot\phi-Na\phi^TK_0\phi\right].
\tag{6}
\]

Let `T0=pi^T M0^(-1) pi/2` and `G0=phi^T K0 phi/2`, with canonical
momentum pi. At fixed canonical data its proper-time energy is
`E=a^(-3)*T0+a*G0`. In a declared homogeneous isotropic volume
`V=a^3*Vc`, lapse and spatial-scale variations yield

\[
E_{\rm kin}=a^{-3}T_0,\quad E_{\rm grad}=aG_0,\qquad
\rho=\frac{E_{\rm kin}+E_{\rm grad}}{V},\quad
p=\frac{E_{\rm kin}-E_{\rm grad}/3}{V}.
\tag{7}
\]

**Proof.** Legendre transformation of (6) gives N*E. Its lapse derivative
therefore gives E, and its scale derivative at fixed phi,pi gives
`p=-partial_a E/(3*a^2*Vc)`, which is (7). Equivalently differentiate
(6) before substituting `pi=a^3*M0*phi_dot/N`; both routes agree. The
finite-energy infinite hierarchy has finite T0 and G0. Multiplication
by fixed positive powers of a preserves its form domain, so the same
identity passes to that hierarchy without an additional depth term.

For source-free evolution, Hamilton's equations give
`dE/dtau=-3*(a_dot/a)*p*V`, where tau is homogeneous proper time.
Consequently `rho_dot+3*H*(rho+p)=0`: the pressure and energy belong to
the same action and obey the corresponding energy balance.

For every nonzero state in this positive-action class,

\[
\rho+p=\frac{2E_{\rm kin}+2E_{\rm grad}/3}{V}>0,\qquad
\rho+3p=\frac{4E_{\rm kin}}{V}\ge0.
\tag{8}
\]

Thus even infinitely many inherited levels cannot turn this classical
kinetic-plus-gradient sector into a nonzero vacuum-form source. Its trace
projection `Lambda_tr=kappa*(rho-3*p)/4` can be nonzero while (8) fails
the vacuum test. A flat Einstein cosmology containing only this lifted
source and no supplied volume term would have the dimensionless prediction

\[
q_{\rm dec}=\frac{\rho+3p}{2\rho}
=\frac{2E_{\rm kin}}{E_{\rm kin}+E_{\rm grad}}\in[0,2].
\tag{9}
\]

It cannot accelerate. This is a conditional homogeneous consequence of
the declared lift, not an identification of the atlas with that cosmology.
The hierarchy's branch count changes allowed modes and state energies;
it cannot change the signs in (8). No observed expansion rate selects them.

**Root work and state restrictions.** With a fixed root value xi, static
elimination gives `E_stat=a*xi^2/(2*chi0)` and root force `f=a*xi/chi0`.
Its pressure is the gradient result `p=-rho/3`. Holding f fixed instead
selects a different family, `xi=chi0*f/a`, with stored energy
`E_stat=chi0*f^2/(2*a)`. Differentiating that latter energy without its
root work would produce the wrong pressure. Direct differentiation gives

\[
dE_{\rm stat}=\frac{E_{\rm stat}}a\,da+f\,d\xi.
\tag{10}
\]

For fixed f, the second term is `-2*E_stat*da/a`, explaining the changed
total derivative. The force-controlled Legendre potential is
`E_stat-f*xi`; it is not the stored positive energy. Equations (7)--(10)
refer to metric variation or quasistatic virtual work, respectively, not
an assumed finite-rate forced trajectory.

This settles the metric source of this positive field lift. A cosmological
Lambda now needs the actual gravitational correlation action, or another
specified sector with its own metric dependence and state. Relabeling the
positive root response as a vacuum coefficient would omit the calculated
pressure. The nonminimal and distributional terms of the cap--wall action,
and a renormalized quantum determinant if used, retain their own equations.

**Bounded check.** The inputs are the proved inherited operator, the
explicit edge action and the existing homogeneous energy-variation API.
The missing connection is the physical measure and pressure. At most30
CPU seconds are allowed for independent three-dimensional edge assembly,
lapse/scale differentiation and the root-work identity, followed by one
illustration. No physical grids, observational fits or field evolution are
run. A failed identity stops the calculation. The positive-action result
redirects the Lambda calculation toward the metric correlation sector;
it does not reopen this regional proof or select a vacuum amplitude.
No novelty claim is made for the underlying variational identities.

## A calibrated observable that distinguishes the nested geometries

The same volume-weighted action supplies an observer relation that retains
the length ratio q. Start with zero initial fields and apply a known smooth
compact root force f(t). Its amplitude and time history are calibrated, so
they cannot be adjusted separately for each geometry. The root readout is
the physical field xi(t) in the same proper clock. For real Laplace s>0,

\[
\widetilde\xi(s)=\chi(s)\widetilde f(s),\qquad
\chi(s)=\frac{R(s^2)}{m_0},\qquad Z(s)=\chi(s)^{-1}.
\]

**Proposition.** Write `chi0=1/(m0*gamma)` for (5). Then

\[
Z(s)=\chi_0^{-1}+\mu_*s^2+O(s^4),\qquad
\mu_*=\frac{m_0}{1-r/\Omega},\qquad \tau_*^2=\chi_0\mu_*.
\tag{11}
\]

Here mu_* is the effective root inertia. The scale tau_* is a response
coefficient, not a first-arrival time. It can be determined from the static
and low-frequency force-to-field response. In c=1 units it has length units;
its physical time scale is tau_*/c.

**Proof.** The positive lower bound on D makes its resolvent differentiable
at lambda=0. Differentiate the scalar equation(1):

\[
R'(0)=-\frac{1}{\gamma^2(1-r/\Omega)}.
\]

Since Z(s)=m0/R(s^2), this gives (11). Independently, the static profile
with unit root value has `phi_n=(b*q)^(-n)` when b*q>1, and `phi_n=1`
when b*q<1. Substitution in the conductance equation gives these profiles;
the stated form closure selects the respective branch. Summing the full
profile's kinetic weight gives

\[
\mu_*=\sum_{n\ge0}b^nm_0q^{3n}\phi_n^2
=\begin{cases}
m_0/(1-q/b),&bq>1,\\
m_0/(1-bq^3),&bq<1.
\end{cases}
\]

Both sums converge and agree with (11). This directly includes all
descendants rather than assigning a root mass by hand.

For an exact-arithmetic expansion error, a second differentiation gives
the scalar moment

\[
A_2=\langle e_0,D^{-3}e_0\rangle
=\frac{1}{\gamma^3(1-r/\Omega)^2(1-r/\Omega^2)}.
\]

The resolvent identity gives
`0<=R(lambda)-R(0)-lambda*R'(0)<=lambda^2*A2` for lambda>=0.
This uses `0<(D+lambda)^(-1)<=D^(-1)`; it retains the actual root moment
instead of replacing it by a worst-case norm. Substituting
lambda=s^2/Omega in the Schur equation gives

\[
0\le\chi_0^{-1}+\mu_*s^2-Z(s)
\le\frac{m_0c^2 A_2}{\Omega^2}s^4.
\tag{12}
\]

This controls the Laplace expansion error. It does not include measurement
error or an approximation to a physical cap--wall field.

The dimensionless relation is explicitly

\[
F_b(q)=\frac{\tau_*^2}{L^2}=
\begin{cases}
\displaystyle\frac{q^{-3}+b}{2(b-q)},&bq>1,\\
\displaystyle\frac{q^{-2}+bq}{2(1-bq^3)},&bq<1.
\end{cases}
\tag{13}
\]

For b=3, q=1/2 and q=2/5 have the same static compliance `L/(3*A)`.
Their dimensionless inertial coefficients are respectively **2.2** and
**3.581730769230769**. Their calibrated dynamic readouts therefore differ
even though a static readout cannot separate them.

With b and L fixed independently, (13) determines the admitted q uniquely.
For q>1/b, the derivative has the sign of `b*q^4+4*q-3*b`, which is
negative for b>=2 and0<q<1. For q<1/b, its sign is that of
`x^2+3*x-1`, where `x=b*q^3<1/b^2<=1/4`; this is also negative. The two
expressions meet at q=1/b, although that case is outside the present strict
operator theorem. Thus F_b is strictly decreasing across the admitted
branches. Equation (11)'s independent physical scale and calibrated source
are essential: an unknown L or an adjustable incident force restores a
degeneracy. This is an inverse relation for the specified geometry and
preparation, not a determination of nature's q from absent observations.

A finite-band readout avoids estimating a derivative at arbitrarily small
frequency. Put `x=s*L`, `mathcal_Z=Z(s)/Z(0)` and
`G_b=r*F_b(q)^2/(Omega^2-r)`. Equations (11)--(12) imply

\[
1+F_b(q)x^2-G_bx^4\le\mathcal Z\le1+F_b(q)x^2.
\tag{14}
\]

At x=.1, b=3, the q=1/2 interval lies inside
`[1.0219789,1.0220001]`; the q=2/5 interval lies inside
`[1.0357893,1.0358174]`, with outward rounding. They are disjoint.
This is an analytic response-expansion uncertainty for the declared action,
not an observational uncertainty. A finite-depth calculator comparison
retains its separate proved depth error. The time/frequency conversion uses
the independently supplied proper length L and c; no cosmological ruler is
inferred from these example inputs.

**Finite preparation and readout.** One admissible preparation is an
externally calibrated smooth nonnegative force, supported in [0,Tf], with
`F1=integral f(t)dt>0`, acting on zero initial data. Its origin is that
root actuator; it is not a claimed primordial preparation. The same wave
energy satisfies `E'=f*xi_dot` and
`|xi_dot|<=sqrt(2*E/m0)`, hence `E(t)<=F1^2/(2*m0)`. This gives a finite
injected-energy budget from the same action.

If observation ends at T>=Tf and additive root-field error is bounded by
epsilon, define chi_T from the observed Laplace numerator over [0,T]
divided by the known complete force transform. Then

\[
|\chi_T(s)-\chi(s)|\le
\frac{e^{-s(T-T_f)}}{m_0\sqrt\delta\,s}
+\frac{\epsilon(1-e^{-sT})}{s\widetilde f(s)}.
\tag{15}
\]

Indeed the retarded kernel is bounded by `1/(m0*sqrt(delta))`, and
`f_tilde(s)>=exp(-s*Tf)*F1>0`. Integrating its tail and the bounded
measurement error proves (15). Positive compliance intervals can then
be inverted to intervals for Z. This separates finite coverage and
measurement error from the expansion and depth errors.

A common constant clock-rate error `t_measured=k*t_proper` rescales the
response argument to k*s: the factor k in the two Laplace integrals
cancels. Thus a declared rate interval gives the corresponding interval
of x in (14); Z is increasing for s>=0 by positive-resolvent order.
This changes the readout calibration, not the proper-time energy or its
equation of state. No observational noise level or clock precision is
assumed in the example intervals above.

The existing30 CPU-second allowance also covers one direct static-profile
sum, an independently differentiated finite Schur recurrence and this paired
readout. No evolution, parameter scan or observational dataset is needed.

## Corollary: global evolution, energy and a local measurement

For u0 in Q and v0 in the Hilbert space, the inherited wave equation

\[
\ddot u+Du=0
\]

has the unique global finite-energy solution

\[
u(t)=\cos(t\sqrt D)u_0+D^{-1/2}\sin(t\sqrt D)v_0,
\qquad
E=\tfrac12\bigl(\|\dot u(t)\|^2+q[u(t)]\bigr)=E(0).
\]

It is continuous in the form norm and continuously differentiable in the
Hilbert norm. Strong twice differentiable solutions require the stronger
initial domain `u0 in Dom(D), v0 in Q`.

**Proof.** Compact inverse gives an orthonormal eigenbasis with positive
eigenvalues mu_j. In each mode the pair `(sqrt(mu_j)*u_j, u_j_dot)` is
rotated by angle `t*sqrt(mu_j)`. Its squared norm is conserved. Summing
over modes proves energy conservation and convergence of the energy-data
series. Dominated convergence gives the stated time continuity. A second
solution with the same data has zero in every mode, proving uniqueness.

The root field `P_root u(t)` is a defined local observable of this
mathematical realization. A causal root force f(t), with zero initial
data and compact time support, has the retarded kernel

\[
K(t)=\mathbf1_{t\ge0}\,P_{\rm root}D^{-1/2}
\sin(t\sqrt D)P_{\rm root}^*.
\]

Its Laplace transform at real s>0 is exactly `R(s^2)`, since each mode
has transform `1/(s^2+mu_j)`. Consequently equation(2) also controls the
finite-depth error in this measured transfer function. The same argument
with a force gives `dE/dt=Re<P_root u_dot,f>`: energy entering the root
is accounted for by the full hierarchy. No probability law or particle
interpretation is required for this deterministic wave result.

Nonzero data can occupy every generation with finite total energy.
For a unit regional vector e and real a with b*Omega*a^2<1, take
`(u0)_v=a^(|v|)*e` and v0=0. Summing the diagonal and link contributions
as geometric series gives

\[
\|u_0\|^2=\frac1{1-ba^2},\qquad
q[u_0]=\frac{\langle e,He\rangle+2ba\operatorname{Re}\langle e,Be\rangle}
{1-b\Omega a^2}.
\]

Thus the global solution need not start with all but finitely many
generations empty. The inequality is an admissible initial-state condition,
not a value selected by the inherited law.

## What this establishes for the map

There is now a constructed infinite hierarchy with one inherited law,
nonzero coupling, an associative local reduction and a finite controlled
observer response. With an invertible link, the response changes at
every added level and converges with an explicit error bound. The
noninteger spectral dimension follows from the same construction.
The corollary supplies global finite-energy evolution and a defined root
input/output measurement for this operator realization.

This is a mathematical operator realization. A cap--wall spacetime,
its metric-dependent action and a vacuum-form stress still have to be
connected to these H,B,Omega,b data before this response is identified
with a cosmological Lambda. The proof does not use a flat-space Dirac
slip or any P001 grid result.

## Reproduction and independent controls

The [calculator](../src/recursive_horizons/nsc_infinite_inheritance.py)
implements the same recurrence without building the exponentially large
tree. A dense tree is limited to256 coordinates for independent checks.
The [nine controls](../tests/test_nsc_infinite_inheritance.py) compare a
noncommuting matrix link with the complete finite tree, the actual subtree
scaling, the scalar closed form and depth bound, the spectral bounds and
the zero-link limit. They also reject missing branch/scale factors and a
reversed noncommuting link. A separately assembled edge action checks the
geometric matrices, their automatic positive margin and the root compliance
L/b. The volume-weighted extension checks both branches of its static
response, independent lapse/scale variations and the root-work identity.
The dynamic readout is checked against the full static-profile mass,
an independently differentiated Schur recursion and a finite-band response
with its separate depth bound.
The proof above establishes the infinite result;
finite tests check the calculator against independent identities.

Independent native Grok review accepted the operator, response, spectral,
ultrametric and wave arguments under their stated assumptions. The edge-action
corollary was added afterward and checked against the separate action assembly.
The manuscript's Section5.1 contains the main theorem and proof; revision
2026-10-18.13 builds to66 pages and the new proof pages pass visual inspection.

A subsequent native Grok review accepted the three-dimensional action,
its pressure, both static profile branches, root inertia and monotonic
inference. It also checked the exact second moment, disjoint finite-band
intervals and finite preparation/readout bounds. Its terminal verdict is
valid under the stated assumptions; it retains the separation from the
actual constrained gravitational and nonminimal cap--wall sectors.
Revision2026-10-18.14 includes the volume action, pressure and dynamic
observable in Section5.2. Its matching67-page PDF builds; pages9--11 pass
visual inspection after correcting inline inequality and Roman-label
rendering. The existing two-build illustration check also passes.

The allowance is30 process-CPU seconds for these small controls and one
illustration, with no physical grid, evolution or large data output. The
first four tests use .375 CPU seconds. The extended nine controls and
illustration pass together in .9375 process-CPU seconds; pytest takes .77
wall seconds. The volume example b=3,q=.5,L=A=1 gives compliance1/3,
falling to1/6 when all spatial lengths double. Equal kinetic and gradient
energies give p/rho=1/3 and the conditional flat q_dec=1. These examples
use declared inputs and no observational calibration.
Neither floating-point arithmetic
nor measured/model uncertainty is included in the analytic depth bound.

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_infinite_inheritance.py -q
python scripts/run_python.py scripts/check_nsc_infinite_inheritance.py
```

For H=3, B=1, b=3, Omega=4, the static response follows exactly from (1):

\[
R(0)=\frac{2}{3+\sqrt{6}},\qquad
D_{\mathrm{eff}}(0)=\frac{3+\sqrt{6}}{2}.
\]

The positive bounded branch selects this quadratic root. The second-order
spectral dimension is `log(3)/log(2)`, which is noninteger: rationality
would require `3^n=2^m` for positive integers. These are derived numbers
for the specified operator data; no observation was used to choose them.
