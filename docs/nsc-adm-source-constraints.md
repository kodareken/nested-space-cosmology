# Coupled spherical sources and their constraint equations

The source construction now retains lapse, shift, radial metric and sphere
radius as independent variables. The canonical Dirac operator, Euclidean
measure factors, first/second geometric vertices and stress projections
are explicit. Their continuum coordinate-invariance identities give the
energy/momentum balance and homogeneous constraint equations of the same
variational problem.

This joins the [vacuum-matched prescription](nsc-vacuum-matched-ctp.md)
to the geometric equations. It is not an evaluation of the complete
renormalized source on the expanding child or a solution of those equations.
The [record](../results/development/adm-source-constraints.json) is produced
by [derive_nsc_adm_source.py](../scripts/derive_nsc_adm_source.py); the
[operator owner](../src/recursive_horizons/nsc_adm_source.py) reuses the
existing spatial Hamiltonian without changing its bytes.

## One spherical frame before gauge fixing

Use signature +--- and

\[
ds^2=N^2dt^2-q^2(dx+\beta dt)^2-r^2d\Omega^2,
\qquad n=N^{-1}(\partial_t-\beta\partial_x),\quad
e_{\hat1}=q^{-1}\partial_x.
\]

The Cauchy half-density is u=r sqrt(q) psi. Reuse the general Hermitian
Dirac Hamiltonian of [Obukhov–Silenko–Teryaev, Eqs.2.14–2.16](https://arxiv.org/pdf/1308.4552),
with their radial shift K^x=-beta. In the nonrotating diagonal spherical
triad, the antisymmetric triad-rotation and spatial pseudoscalar terms
vanish; the sphere connection remains in its Dirac eigenvalue kappa.
The density transformation includes the time derivatives of q and r.
It gives

\[
\boxed{H_L=\sigma_2\tfrac12\{N/q,P_x\}
 +\sigma_1 N\kappa/r-\tfrac12\{\beta,P_x\}I,
\qquad P_x=-i\partial_x.}
\]

This specializes to the previously derived PG operator at N=q=1.
The code retains the previous factorized Fourier representation
H_perp=sqrt(N)[sigma2 {q^-1,P_x}/2+sigma1 kappa/r]sqrt(N),
then adds the shift. It converges to the displayed continuum operator;
the finite Fourier derivative does not obey a continuum product rule
on arbitrary grid vectors. No new convergence claim is made here.

For the canonical state, i dot C=[H_L,C]. At fixed C the continuum
Hamiltonian vertices are

\[
\begin{aligned}
\delta_NH&=\sigma_2\tfrac12\{\delta N/q,P_x\}
                 +\sigma_1\kappa\delta N/r,\\
\delta_\beta H&=-\tfrac12\{\delta\beta,P_x\}I,\\
\delta_qH&=-\sigma_2\tfrac12\{N\delta q/q^2,P_x\},\\
\delta_rH&=-\sigma_1N\kappa\delta r/r^2.
\end{aligned}
\]

The implementation differentiates its actual finite factorization and
retains the mixed and second-variation contacts. It does not differentiate
a re-prepared ground state in place of the prescribed evolving C.

## Keep the lapse in the covariant operator

The Euclidean density is chi=sqrt(N)u. For a real Euclidean shift b_E,
write B_E={b_E,P_x}/2, W=N^-1/2, beta_D=sigma3 and P_tau=-i partial_tau.
Then a Hermitian Euclidean representation is

\[
\boxed{D_E=W\beta_D[P_\tau-B_E-iH_\perp]W.}
\]

The factors W act on both sides of the time derivative, including when
N varies with time. For t=-i tau, the metric continuation is b_E=-i beta_L.
A physical real shift consequently does not become a real Euclidean
shift by declaration. The continued operator is not passed to a
Hermitian eigensolver without the required analytic prescription.

At zero shift, this Euclidean representation is related to the old
static fiber by constant beta_D conjugation. Its modulus is unchanged.
For a log-lapse variation a=delta log N,

\[
\delta D_E=-\tfrac12\{a,D_E\}
 +W\beta_D[-\delta B_E-i\delta H_\perp]W.
\]

The first term is a required measure contribution. For a spatially
uniform lapse depending only on time, it cancels the spurious change of
the spatial term from H_perp=N(t)H0: only the properly weighted time
operator changes. The second variation carries the corresponding contact
terms. Rescaling H while holding this measure fixed would be a different
source prescription.

The history operator is executable with arbitrary supplied time-momentum
matrices and time slices of the existing spatial geometry. Its finite
Euclidean example is an operator-construction control, not a physical
thermal state or a continuum vacuum calculation.

## Four independent source projections

For a contribution to the common functional define
F_A=-delta Gamma/delta A_Delta, per dt dx after angular integration.
The canonical fermion contribution is Tr(C delta H/delta A); the geometric
branch contributes -delta B_L/delta A. Sum the contributions once.
The total variational equations are

\[
\boxed{F_N=F_\beta=F_q=F_r=0.}
\]

These include gravitational Euler derivatives. They do not mean that
the matter stress is zero. For an individual stress contribution, its
spherically averaged orthonormal components are

\[
\boxed{
\rho=\frac{F_N}{4\pi q r^2},\quad
T_{\hat0\hat1}=\frac{F_\beta}{4\pi q^2r^2},\quad
p_r=-\frac{F_q}{4\pi N r^2},\quad
p_\perp=-\frac{F_r}{8\pi N q r}.}
\]

These follow from the inverse-metric variation, including the mixed
time/radial entry. The outward energy flux in an orthonormal frame is
T^hat0hat1=-T_hat0hat1 with this signature. Nodal numerical gradients
must first be divided by the spatial spacing, and the complete angular
multiplicities must be restored. A single radial block is not the full
four-dimensional source.

## Coordinate invariance gives the balance equations

For an infinitesimal vector xi=(xi^t,xi^x), the induced variations are

\[
\begin{aligned}
\delta N&=\xi\cdot\partial N+N(\partial_t\xi^t-\beta\partial_x\xi^t),\\
\delta q&=\xi\cdot\partial q+q(\partial_x\xi^x+\beta\partial_x\xi^t),\\
\delta r&=\xi\cdot\partial r,\\
\delta\beta&=\xi\cdot\partial\beta+\partial_t\xi^x
-\beta\partial_x\xi^x+\beta\partial_t\xi^t
-(\beta^2+N^2/q^2)\partial_x\xi^t.
\end{aligned}
\]

Substitute these into the same action variation and integrate by parts.
For interior transformations, with matter/link equations imposed and
the state, regulator, domain and boundary data transported consistently,
define

\[
\mathcal E=N F_N+\beta F_\beta,\qquad
\mathcal J=-N\beta F_N-(\beta^2+N^2/q^2)F_\beta+q\beta F_q.
\]

The two Ward identities are

\[
\boxed{\partial_t\mathcal E+\partial_x\mathcal J
=\sum_A F_A\partial_t A,}
\qquad
\boxed{\partial_tF_\beta+\partial_x(qF_q-\beta F_\beta)
=\sum_A F_A\partial_x A.}
\]

For the canonical matter sector, the integral of E is Tr(C H_L);
the first identity is its geometric-work and flux balance. For geometric
Euler derivatives these are variational balance quantities, not an
independent definition of ADM or quasi-local energy. Actual interface and
asymptotic energy accounts retain their boundary terms.

On a static zero-shift background the radial identity reduces to the
previously used equation p_r'+(N'/N)(rho+p_r)+2(r'/r)(p_r-p_perp)=0.
The old conservation calculation is reused rather than rerun.

## Constraint evolution is homogeneous

When the two spatial metric equations F_q=F_r=0 hold, the identities give

\[
(\partial_t-\beta\partial_x)F_\beta
=N_x F_N+2\beta_xF_\beta,
\]
\[
(\partial_t-\beta\partial_x)F_N
=\beta_xF_N+\frac{N}{q^2}\partial_xF_\beta
 +2\left(\frac{N_x}{q^2}-\frac{Nq_x}{q^3}\right)F_\beta.
\]

No unrelated energy injection or pressure appears. The constraints form
a homogeneous system, so zero constraints are preserved in an evolution
with suitable uniqueness and compatible boundary data. This is not a
proof that the nonlocal metric evolution exists, is hyperbolic, or admits
arbitrary Cauchy data. It is also not the quantum constraint algebra.

## New binding evidence and remaining computation

The derivation checks the ADM Lie variation, stress projections and the
two constraint identities symbolically. A four-time-slice, twelve-spatial-
point Euclidean example then compares the new first/second operator
vertices and regulated action derivatives with independent matrix
differences for each field and their mixed direction. It also matches
the old static fiber by its exact unitary map and evaluates the canonical
four-component force budget at one fixed state.

No old source integral, spectrum campaign or clock-production experiment
is repeated. These checks establish the new combined operator/variation
interface. They do not establish finite-grid diffeomorphism invariance.

```sh
python -B scripts/derive_nsc_adm_source.py --check
```

The next physical computation must evaluate the vacuum-matched branch and
canonical state contribution with these independent vertices, the common
continuum subtraction and actual parent–child domain. The full source
values and a solution of the four coupled equations remain required.
The radius-only real-time prototype is not promoted to that result.
