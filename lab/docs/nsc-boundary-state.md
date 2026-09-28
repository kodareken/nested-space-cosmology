# Explicit real-time boundary state and preparation

The finite quadratic parent–child operator now has an executable boundary
state interface: full two-time occupied/empty kernels, mixed preparation
terms, reconstructed parent covariance and normalized Gaussian endpoint
action. This fills in the state dependence retained abstractly in the
[common-source derivation](nsc-common-source-derivation.md).

The owner is [nsc_boundary_state.py](../src/recursive_horizons/nsc_boundary_state.py).
It extends the accepted [retarded elimination](nsc-energy-transfer.md)
and [Klich state functional](nsc-influence.md). No mass, scalar source,
Markov approximation or uncorrelated-vacuum assumption is introduced.

## The boundary kernels are explicit

Use hbar=1, C_ij=<c_j^dagger c_i>, and one fixed canonical spin frame.
U_c(t,t0) is the isolated child's unitary. Define F(t)=B(t)U_c(t,t0)
and eta(t)=F(t)c_0. The existing parent equation is

\[
(i\partial_t-H_p(t))p(t)-\int_{t_0}^t ds\,\Sigma^R(t,s)p(s)=\eta(t),
\qquad\Sigma^R(t,s)=-i\theta(t-s)F(t)F(s)^\dagger.
\]

The occupied and empty source correlations are respectively

\[
L(t,s)=F(t)C_{cc}F(s)^\dagger,\qquad
G(t,s)=F(t)(I-C_{cc})F(s)^\dagger.
\]

Thus the lesser, greater and Keldysh kernels are

\[
\boxed{\Sigma^<=iL,\quad\Sigma^>=-iG,\quad
\Sigma^K=-i(G-L),\quad\Sigma^R=\theta(t-s)(\Sigma^>-\Sigma^<).}
\]

The implementation chooses theta(0)=1/2; changing that isolated value
does not change the ordinary memory integral. These fermionic correlations
have energy-squared units. They are not classical random-force or metric
stress-noise covariances.

Initial parent–child correlations supply a second explicit kernel:

\[
\boxed{A_0(t)=\langle p_0^\dagger\eta(t)\rangle=F(t)C_{cp},
\qquad\langle\eta(t)p_0^\dagger\rangle=-A_0(t).}
\]

Here A_0 has energy units. The link, state and history must transform
together under a change of spin frame. Neither C_pc nor C_cp is set to zero.

## Reconstruct the state

Write p(t)=A(t)p_0+V(t)c_0. With A(t,s) the parent response to an impulse,

\[
V(t)=-i\int_{t_0}^t ds\,A(t,s)B(s)U_c(s,t_0).
\]

The complete two-time parent covariance is then

\[
\boxed{\begin{aligned}
C_{pp}(t,s)={}&A(t)C_{pp}A(s)^\dagger+V(t)C_{cc}V(s)^\dagger\\
&+A(t)C_{pc}V(s)^\dagger+V(t)C_{cp}A(s)^\dagger.
\end{aligned}}
\]

The owner returns all four terms. The last line is part of the prepared
state, not an adjustable energy injection. Its force and current follow
from the same Hamiltonian as in the common-source derivation.

## Normalized boundary action and pure preparations

For full forward/backward histories, set P=U_minus^dagger U_plus and
Q=I-C_0+C_0P. The existing exact normalized functional Z=det Q has the
finite Grassmann representation

\[
Z=\int d\bar\xi\,d\xi\,e^{-\bar\xi Q\xi},\qquad
\int d\bar\xi\,d\xi\,e^{-\bar\xi\xi}=1.
\]

Integrating the child endpoint variables, when Q_cc is invertible, gives

\[
Q_{p,\mathrm{eff}}=Q_{pp}-Q_{pc}Q_{cc}^{-1}Q_{cp},\qquad
\boxed{Z=\det Q_{cc}\,\det Q_{p,\mathrm{eff}}.}
\]

This retains the child's determinant once and all correlated preparation
data. Both factors use the complete history matrix Q. Neither is an
independently normalized influence or an uncoupled-child determinant.
A singular Q_cc invalidates this pivot, not necessarily the full
functional; use the unfactored determinant or another admissible pivot.

Pure covariances are supported directly because C_0 and I-C_0 are never
inverted. Identical histories give Q=I and Z=1. This is an endpoint
representation of the canonical contour functional, not an identification
with the finite Euclidean heat determinant. It does not assert positivity
of a ratio of influences.
The returned -i log Z uses the principal branch near identical histories;
it does not determine a global continuum determinant phase.

## One new integration calculation

Use the existing general-lapse/general-radial-metric smooth AP cell,
12 spatial points and one representative angular block. A relative phase
of 0.37 is applied on the parent subspace of the ground covariance,
preserving its purity and cross-correlations. This is a declared
preparation, not the expanding child's Unruh state. At t=0.63, s=0.27,
V is assembled independently by a 32-node causal convolution.

| Comparison | Observed value |
|---|---:|
| Maximum covariance error against full unitary evolution | 1.0201e-15 |
| Frobenius error if mixed preparation terms are omitted | 0.67734 |
| Endpoint factorization error against the existing Klich functional | 6.4737e-16 |

The integration passed in 0.007 seconds. The two equality tolerances are
2e-12; omitting the mixed terms must change the result by more than 0.01.
This is one finite implementation comparison, not a convergence study.
No previous scientific generator or full suite was run.

```sh
PYTHONPATH=src python -B -m unittest discover -s tests -p test_nsc_boundary_state.py -v
```

## Remaining absolute-source construction

For a specified quadratic Hamiltonian, state and partition, these kernels
are determined. They should no longer be described as an unspecified
Gaussian state influence. The formulas accept time-dependent links and
child Hamiltonians through their supplied unitary histories; the numerical
example connects them on the existing static cell.

The room partition is distinct from the compact/light allocation
Gamma5=Gamma_light,ren+Gamma_H+C_nu,mu. Applying this construction to that
allocation requires its canonical Hilbert-space decomposition and joint
state; a spectral weight alone is not that decomposition. Absolute
continuum vacuum and measure/phase/compensator matching, interacting
sectors and self-consistent geometry remain inputs to the full source
equation. No coefficient is chosen here to obtain a desired neck stress.
