# Connected geometric influence vertices on the saved finite state

This calculation uses the authenticated `nf128_baseline_dt0.0005` state at
$T=0.3$, its full ambient fermion band, and its frozen canonical geometry
basis. It takes no evolution step and performs no new source preparation.
The state is the actual mixed rank-six Gaussian, with occupation weights
$(0.75,0.75,0.5,0.5,0.25,0.25)$ and covariance
$C=\Phi\operatorname{diag}(c)\Phi^\dagger$. Gaussian elementary fields obey
Wick's theorem. Their quadratic force observables need not be Gaussian.

Owners are
[the calculation](../src/recursive_horizons/nsc_discovery_geometric_vertices.py),
[CLI](../scripts/derive_nsc_discovery_geometric_vertices.py), and
[independent tests](../tests/test_nsc_discovery_geometric_vertices.py).
Inputs are authenticated through
[the coupled-memory handoff](nsc-discovery-coupled-memory.md), the preserved
pair record, and its portable basis. The source energy and forces remain
those of the [same conformal action](nsc-spherical-conformal-gauge.md).

## The actual mixed vertex

For canonical geometry $Q=W a_Q$ and the quadrature map $A_g$, let
$q_j=A_gW_{:j}$. The conformal representative generator is
$H_G=U_f^\dagger(\sigma_2P+\kappa Q\sigma_1)U_f$ and therefore

\[
G_j=\frac{\partial H_G}{\partial a_{Qj}}
=\begin{pmatrix}0&g_j\\g_j&0\end{pmatrix},\qquad
g_j=\kappa U_f^\dagger\operatorname{diag}(q_j)U_f.
\]

The geometric interaction is already present: the Hamiltonian variation
is $\delta a_{Qj}\,\varphi^\dagger G_j\varphi$; its sign in the canonical
Lagrangian is negative. No new cubic term, mass, or force is added.
Indices $j=5,6$ are the first two child-detail geometry modes; $j=0$ is the
broad constant parent mode. They are geometry-coordinate indices, not source
columns or quantum-observer indices. The recorded profile hashes, extrema,
and sampled power outside the child interval identify these modes. Finite
Fourier profiles are not asserted to have exact compact support.

The angular factor is absent from $G_j$. The mean energy gradient
$M\operatorname{Tr}(CG_j)$, with $M=4\kappa$, agrees with the existing
`Q_energy_gradient`; the corresponding source term in $\dot\pi_Q$ has the
opposite sign. All covariance and ambient paths enter this trace.

## A third connected bilinear vertex

For one representative-channel observable $\hat f_j=\hat c^\dagger G_j\hat c$,
the characteristic function follows from the existing normalized Gaussian
influence determinant:

\[
Z_j(\lambda)=\det[I-C+C e^{-i\lambda G_j}],\qquad
\log Z_j(\lambda)=\sum_{n\ge1}\frac{(-i\lambda)^n}{n!}\kappa_n.
\]

In particular,

\[
\kappa_3=\operatorname{Tr}(CG_j^3)
-3\operatorname{Tr}(CG_jCG_j^2)
+2\operatorname{Tr}[(CG_j)^3].
\]

For the saved state, child mode $6$ gives
$\kappa_3=1.14989407922813\times10^{-5}$; mode $5$ also gives a nonzero
third cumulant, about $1.11915902\times10^{-6}$. These are connected
three-bilinear influence vertices, not non-Gaussian elementary interactions.
For the sign convention $\Gamma=-i\log Z$, its real cubic term is
$+\lambda^3\kappa_3/6$. This is a counting-field/difference-source influence
coefficient; it is not identified with a complete second-order retarded
metric response or a local term installed in the action.

The independent counting-field calculation applies the full-band exponential
to the six actual columns, then uses Sylvester's exact thin determinant.
It does not replace $G$ by $\Phi^\dagger G\Phi$. Derivative steps are
$h=0.04,0.02$; the held-out counting field is $\lambda=0.03$. The odd phase
after subtraction of the linear mean is compared with the predicted cubic
term, retaining the $O(\lambda^5)$ remainder. Counting-field amplitudes are
not the exterior occupation perturbations used in the trajectory prediction.
Small independent Fock-space tests check the cumulant and contraction formulas.

## Occupied/empty cut and an instantaneous causal slope

The existing Wick contraction is

\[
S^>_{AB}=\operatorname{Tr}[C A(I-C)B],\quad
S^<_{AB}=\operatorname{Tr}[C B(I-C)A],\quad
S^>_{AB}-S^<_{AB}=\operatorname{Tr}(C[A,B]).
\]

Take $A=i[H_G,G_5]$ and $B=G_0$. The first is the instantaneous Heisenberg
derivative of the actual child force vertex; it is not another interaction.
With $H_{\mathrm{int}}=+JB$, the retarded convention is
$\chi=-i\theta\langle[A(t),B(s)]\rangle$. On this saved slice the occupied/
empty cut yields the instantaneous child/parent response slope
$-0.0150914068335309$, with a cut identity residual about $2\times10^{-17}$.
The parent profile is broad and includes the child: this coefficient is a
response between specified canonical geometry modes, not a solely spatial
exterior contribution.

The actual covariance is nonstationary. No frequency fit, physical pole, or
off-axis resolvent sample is reported. A stationary frozen-slice spectrum
would have its own condition $[H,C]=0$ and would remain distinct from the
time-dependent coupled kernel. The result supplies a finite occupied/empty
cut and factorization relationship; no cosmohedron equivalence is derived.

## Angular multiplicity and scope

The existing source verifies the mean $M\operatorname{Tr}(CG)$, and the
influence function above derives a single representative-channel cumulant.
These two possible extensions are different:

\[
\kappa_3(M\hat f)=M^3\kappa_3(\hat f),\qquad
\kappa_3\!\left(\sum_{s=1}^{M}\hat f_s\right)=M\kappa_3(\hat f)
\quad\text{for an explicit product of identical channel states}.
\]

The implemented mean copy count does not construct or authenticate that
angular product state. Total stress noise therefore remains unmapped.
Neither extension is silently substituted into the geometric force.

## Commands and preservation

The default command prints a bounded calculation, with no evidence writes:

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_geometric_vertices.py
```

An explicit record is exclusive and requires every computational source and
read input to match one full frozen Git commit. Replace the placeholder below
after the source checkpoint:

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_geometric_vertices.py --write --science-commit FULL_FROZEN_COMMIT
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_geometric_vertices.py --check
```

`--check` authenticates those bytes and replays the finite saved-slice
calculation. It takes no evolution step. The numerical CPU ceiling is ten
seconds. This adds a derived geometric influence connection without changing
the dynamics, earlier records, paper, or incoming-gate domain.
