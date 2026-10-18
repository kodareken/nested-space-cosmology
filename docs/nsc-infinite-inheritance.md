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
The [four controls](../tests/test_nsc_infinite_inheritance.py) compare a
noncommuting matrix link with the complete finite tree, the actual subtree
scaling, the scalar closed form and depth bound, the spectral bounds and
the zero-link limit. They also reject missing branch/scale factors and a
reversed noncommuting link. The proof above establishes the infinite result;
finite tests check the calculator against independent identities.

The allowance is30 process-CPU seconds for these small controls and one
illustration, with no physical grid, evolution or large data output. The
first four tests use .375 CPU seconds. Neither floating-point arithmetic
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
