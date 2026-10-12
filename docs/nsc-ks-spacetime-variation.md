# The transmitting CTP variation in KS spacetime coordinates

The new owner evaluates a raw KS **spacetime** variation through the existing
whole-field PG Dirac operator. It retains the spatial phase that is lost if
KS time is identified with PG time. This is a first variation at the owned
reference geometry, not a selected metric history or an endpoint interpolation.

## Same field, fixed PG Cauchy representation

Reuse the [owned coordinate map](nsc-pg-ks-metric-pullback.md):

$$
T=T(\rho),\quad z=\tau+S(\rho),\quad
T'=-a_0^{-1},\quad S'=\beta_0/a_0^2.
$$

Let $\xi^B(T,z)$ be a smooth raw-field direction for
$(N_K,\beta_K,q_{\rm ADM},r_K)$. The same metric direction in PG variables is

$$
\delta x_P^A=J^A{}_B(\rho)\xi^B(T(\rho),\tau+S(\rho)),
$$

$$
\partial_\rho\delta x_P
=J'\xi+J\left(-\frac{\partial_T\xi}{a_0}
+\frac{\beta_0}{a_0^2}\partial_z\xi\right),
\qquad x_P=(\log N_P,\beta_P,\log q_P,\log r_P).
$$

The PG Cauchy representation and canonical $C_0$ are held fixed. This field
variation does not require a change between two Cauchy surfaces. Variations
of a physical initial-state preparation or of the Cauchy surfaces themselves
would be additional chain-rule contributions and are not inferred here.

The [ADM source owner](nsc-adm-source-constraints.md) already includes the
time-dependent half-density transformation $\chi=r\sqrt q\,\psi$. Differentiate
its canonical operator once:

$$
\delta H=-i\left(\delta v\,\partial_\rho
+\tfrac12\partial_\rho\delta v\right)+\delta M,
\quad v=(N/q)\sigma_2-\beta I,
\quad M=-Nm\sigma_1-N\lambda\sigma_3/r.
$$

No separate measure force, instantaneous hopping matrix or new action term
is inserted. The variation remains smooth through the fixed transmitting
cut. Compact spatial support removes the exterior integration boundary
term, and the common seam traces cancel its two internal contributions.

## Fourier component and conditional CTP derivative

The evaluated direction uses the existing compact collar profile,
$f(T(\rho))=s(\rho)$, with $\xi=f(T)e^{-i\omega z}$. Consequently

$$
\xi(\tau,\rho)=e^{-i\omega\tau}s(\rho)e^{-i\omega S(\rho)},
$$

$$
\partial_\rho\xi=e^{-i\omega\tau-i\omega S}
\left(s'-i\omega S's\right).
$$

The source-mode vertex $M_\omega$ is computed from the archived normalized
fields before any energy integration. The full operator still acts with
$dE/(2\pi)$; a finite sample does not close the unobserved bulk. The paired
vertex obeys $M_{-\omega}=M_\omega^\dagger$. A real cosine direction therefore
has $H_1(\tau)=(e^{-i\omega\tau}M_\omega+
e^{i\omega\tau}M_\omega^\dagger)/2$.

The owner integrates the reference Duhamel kernel on caller-supplied PG time
coordinates, $\delta U=-iU_0\int U_0^\dagger H_1U_0\,d\tau$, and contracts
it through the existing Gaussian `ctp_first_variation`, with the same
$C_H\oplus n_{\rm in}$ source fibers and the existing spectral weights.
Time coordinates and $\omega$ are response arguments; neither is a physical
duration or a chosen spacetime history.

## Verification and limits

The focused record compares the weak vertex with the strong canonical
variation and with a separately perturbed existing Hamiltonian. It also
checks the zero-transfer limit against the archived raw-KS vertex, conjugate
frequency exchange, and the conditional CTP tangent/trace identities.
No horizon/scattering/state generator is rerun.

The record covers all 33 retained groups and 63 signed families on the
already recorded derivative-control frequency sample. It does not assert
convergence of a new full spectral trace or compute an absolute stress.

The [author's action inventory](nsc-declared-action-scope.md) applies:
the full Gaussian determinant, reference allocation and local action are
counted once. The smooth geometric seam match does not set their bulk
variation to zero. The number $93.54264532195464$ remains a diagnostic
time-node derivative, not a target for an interface counterforce.

This owner supplies the missing spacetime chain rule for an explicit
direction. **Complete physical `EndpointBranchJets` remain OPEN:** the old
two diagnostic end-node directions have not been identified with the
physical history variation, and the nonlinear history/state evolution has
not been evaluated. No matrix is reshaped to fill that slot. B2, kernel
selection and extended stationarity remain unevaluated; metric stepping and
moving-surface Z3 stay closed.
