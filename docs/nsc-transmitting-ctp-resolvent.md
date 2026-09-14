# Transmitting resolvent in the causal CTP memory kernel

This calculation extends the [static transmitting response](nsc-transmitting-cross-resolvent.md)
to a time-dependent frequency-transfer vertex. It uses the same whole-line
PG Dirac operator, transparent seam and 33 retained channels. The packet
projection, scales and source state are not changed.

## 1. Which object enters the CTP functional?

Use the existing canonical convention $G^R(t,s)=-i\theta(t-s)U(t,s)$.
For the stationary reference operator this fixes the sign conversion

$$
G_J^R(z)=-\mathcal R(z),\qquad
D_J^R(z)=-\mathcal K(z),\qquad
G_J^A(\bar z)=G_J^R(z)^\dagger.
$$

Here $D_J^R$ is the inverse retarded kernel. It contains the memory of the
eliminated bulk, including its energy-dependent cross block. It is not a
time-independent Hamiltonian and cannot replace a full branch unitary.
The published upper-half-plane zeros in the parent-from-child block are
preserved; the advanced block has the adjoint direction. Neither branch is
symmetrized into a Hermitian hopping matrix.

In the triangular retarded/advanced Keldysh convention, the projected contour
propagator and its inverse have the imported form

$$
\widehat G_J=
\begin{pmatrix}G_J^R&G_J^K\\0&G_J^A\end{pmatrix},\qquad
\widehat D_J=
\begin{pmatrix}D_J^R&-D_J^R\circ G_J^K\circ D_J^A\\0&D_J^A\end{pmatrix}.
$$

Products are two-time convolutions, with the preparation at the initial
surface retained. These expressions specify where the computed response
enters the existing Gaussian theory. They do not determine $G_J^K$ from
$G_J^R$. The numerical evaluations are at nonreal Laplace energies; no
real-axis state spectrum or fluctuation relation is imposed there.

## 2. State and preparation stay in the same full Gaussian action

On one common PG Cauchy slice, a specified full covariance $C_0$ would give

$$
G_J^<(t,s)=iJ^\dagger U(t,t_0)C_0U(s,t_0)^\dagger J,
$$

$$
G_J^>(t,s)=-iJ^\dagger U(t,t_0)(I-C_0)U(s,t_0)^\dagger J,
\qquad G_J^K=G_J^>+G_J^<.
$$

The packet/complement split must retain its two preparation cross blocks as
well as its two diagonal covariance blocks. These are exactly the four
terms already owned by `GaussianBoundaryState.reconstruct`; the complementary
bulk is part of the existing field, not an added reservoir. The packet
projector is distinct from the physical parent/child support projectors.

The normalized closed-time-path action remains the existing full-space one:

$$
Q=I-C_0+C_0V_-^\dagger V_+,\qquad \Gamma_G=-i\log\det Q.
$$

In a finite regulated packet/complement basis, `overlap_schur` preserves
both determinant factors and the initial correlations. Dropping the
complementary factor or applying `influence` directly to $\mathcal R$ would
change this action. The Gaussian complementary factor is not silently
identified with the separately unresolved $\Gamma_{\rm rest}$.

The serialized NSC covariance is specified on the old spacelike seed surface,
in its transported frequency basis. Z1 explicitly did not map that state
onto a global constant-PG-time spatial slice. Thus its occupations cannot be
placed in the four packet slots, or copied into the unobserved complement.
The full PG preparation kernel, including its cross correlations, remains
unevaluated. Retarded response alone cannot select it.

## 3. New temporal metric variation

Apply the existing compact spatial metric vertex $W_A=\delta H_D/\delta g_A$
with harmonic time dependence $e^{-i\omega t}$, where the four coordinates
are $(\log N,\beta,\log q_{\mathrm{PG}},\log r)$. This is a linear-response
source direction about the frozen reference, not a metric trajectory or a
chosen evolution duration. The real perturbation includes its conjugate
harmonic. With $z_o-z_i=\omega$ real and equal nonzero damping,

$$
\delta\mathcal R_A(z_o,z_i)
=-J^\dagger R(z_o)W_A R(z_i)J,
$$

$$
\delta D_{J,A}^R(z_o,z_i)
=\mathcal K(z_o)\,\delta\mathcal R_A(z_o,z_i)\,\mathcal K(z_i).
$$

These are the double-frequency forms of the imported retarded Dyson
insertion. Unlike the diagonal static derivative, they include frequency
transfer from the time-dependent source. Their zero-transfer limit is the
locked static metric kernel.

Reuse the archived incoming fields at $z_i=0.2+0.25i$. The outgoing frequency
is $z_o=0.6+0.25i$: its real part is the other existing energy, and equal
damping makes the transfer $\omega=0.4$ a real harmonic. This frequency is
a response probe, not a fitted physical oscillation. Only the new outgoing
response and driven sensitivity are solved.

For $U_i=R(z_i)J$, the new sensitivity solves

$$
(H_D-z_o)\,\delta U_{oi}=-W_AU_i.
$$

Its continuous seam trace cancels the projected distributional pieces,
including the varied principal coefficient. Causal inflow is fixed exactly
as in Z2a; no condition is imposed at the outflow endpoint. An independent
weak-resolvent integral uses $L_o=R(\bar z_o)J$:

$$
\delta\mathcal R_A=\int s(\rho)
\left[\frac{i}{2}(L_o^\dagger V_AU_i'-L_o'^\dagger V_AU_i)
-L_o^\dagger M_AU_i\right]d\rho.
$$

The compact profile and fixed canonical packets remove outer endpoint
variations in this test. This does not supply the physical KS-history
endpoint pullback. The adjoint check reverses both frequency labels and
uses the conjugate harmonic; comparing an unswapped advanced vertex would
test the wrong identity.

## 4. Interface inventory and stopping gate

`CTPRetardedMemory` now supplies the computed retarded/advanced inverse
kernels and their off-diagonal-frequency metric vertices. Its conversion
to `EndpointBranchJets` is rejected until full unitary history derivatives,
the physical preparation and the KS endpoint-coordinate map are supplied.
In particular, a Laplace response derivative is not $dV_\pm/dg_\Delta$.

The existing `endpoint_ctp_pullback` and `match_endpoint_variations` receive
explicit missing transmitting jets and boundary derivative. The eight
stored Weyl coefficients, their exact extraction residual zero and the
diagnostic value $93.54264532195464$ are preserved. No numerical two-sided
boundary mismatch is fabricated.

**The retarded-memory part of C1 passes; full state-dependent C1 and C2--C4
remain OPEN.** The next input is the same state's full PG Cauchy preparation
and its packet/complement correlations, together with the physical
history-to-KS-endpoint pullback. More static link matrices cannot fill these
slots. No new stress, selected history or metric timestep is produced.
The fixed cut keeps moving-neck Z3 outside scope.

The [record](../results/development/nsc-transmitting-ctp-resolvent.json)
contains the new temporal-vertex residuals, full missing-field inventory and
authenticated numerical data. Earlier generators and certificates are reused.

```sh
python3 scripts/derive_nsc_transmitting_ctp_resolvent.py --check
python3 -m pytest -q tests/test_nsc_transmitting_ctp_resolvent.py
```
