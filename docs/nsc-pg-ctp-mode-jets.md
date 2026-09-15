# Transmitting CTP variation in the physical PG mode basis

The prepared common-PG state now has a metric-variation kernel in its own
source-mode representation. This connects the same Dirac field that generated
the covariance to the first-order time-ordered CTP insertion. It preserves the
whole transmitting operator and its frequency mixing, with no assigned
instantaneous link matrix.

## Physical insertion

Reuse the globally normalized fields $\Phi_E(\rho)$, their source covariance
$C_H\oplus n_{\mathrm{in}}$, and the existing smooth compact metric direction
$s(\rho)$ on $(-1,1)$. For the four PG coordinates
$(\log N,\beta,\log q_{\mathrm{PG}},\log r)$, the evaluated kernel is

$$
M_A(E_o,E_i)=\int d\rho\,
\Phi_{E_o}(\rho)^\dagger\,delta H_A\,\Phi_{E_i}(\rho).
$$

Its action on source amplitudes uses $dE_i/(2\pi)$. The sign is fixed by
the existing resolvent derivative: the archived weak-resolvent integrand is
minus this Hamiltonian insertion. The new calculation compares that weak
form with direct differentiation of the Z1 Hamiltonian acting on the same
fixed canonical half-density columns.

Both sides use the continuous transmitting field at $\rho=0$. The compact
profile vanishes at the outer integration endpoints, while the two seam
terms cancel. Neither a sharp-projector distribution nor its cancellation
is assigned a stress tensor or a Hamiltonian hopping block.

## Time-ordered and state coupling

Import the first-order propagator identity on the source-mode space, where
the reference Hamiltonian is multiplication by $E$. The computed insertion
therefore determines the kernel

$$
\frac{\delta U(t_f,t_i)}{\delta g^A(\tau)}(E_o,E_i)
=-i\,e^{-iE_o(t_f-\tau)}M_A(E_o,E_i)e^{-iE_i(\tau-t_i)}.
$$

The time coordinates are arguments of this conditional derivative; the
calculation does not select a duration or a geometry history. Hermitian
exchange of the incoming and outgoing energies gives the unitary tangent
condition for every insertion time.

The CTP control contracts this insertion with the physical source-fiber
covariances through the existing `ctp_first_variation`. Spectral quadrature
weights multiply the sampled operator as $\sqrt{w_o w_i}/(2\pi)$, while
$C(E)$ remains the multiplication operator on its source fiber. Closed
infinity columns are removed before the finite matrix control. The horizon
correlations are retained; negative frequencies use the previously derived
paired map and the complement on the open source fiber.

These sampled CTP differentials check the actual action contraction. They
are not an absolute stress, a completed infinite trace or a closed evolution
of the eight packet observables. The unsampled bulk remains in the full
source-mode kernel.

## Evaluated scope and result

All 33 retained groups are covered, with 63 signed angular families. For
each family two existing positive-frequency nodes near 0.2 and 0.6 and their
negative-frequency partners are used. Massive columns are read from the
authenticated PG-mode artifact and continued only through the compact
variation region. The LLL uses its separate exact characteristic map.
No horizon/scattering restart or covariance integration is repeated.

| Check | Maximum | Tolerance |
|---|---:|---:|
| Weak form versus direct Hamiltonian variation | $5.63\times10^{-10}$ | $3\times10^{-8}$ |
| Spatial quadrature refinement | $1.86\times10^{-11}$ | $3\times10^{-8}$ |
| Hermitian energy exchange | $3.56\times10^{-15}$ | $3\times10^{-11}$ |
| CTP tangent condition | $2.78\times10^{-17}$ | $3\times10^{-11}$ |
| CTP differential versus source trace | $2.39\times10^{-18}$ | $3\times10^{-11}$ |
| Seam continuity | $0$ | $3\times10^{-11}$ |

**PG mode-space insertion PASS; full B1 OPEN.** The actual PG first-history
derivative is now represented on the prepared field's source fibers.
`EndpointBranchJets` uses raw KS endpoint coordinates and a Cauchy-slice
pullback, which this compact PG direction does not supply. Its transmitting
fields are therefore not filled with these matrices under a new label.
The remaining boundary action, two-sided Weyl mismatch and stationary
history follow that pullback. The Weyl diagnostic and prior certificates
are preserved; moving-neck Z3 and metric stepping remain outside scope.

The [record](../results/development/nsc-pg-ctp-mode-jets.json) and authenticated
field/vertex artifact retain each family's residuals and source data.

```sh
python3 scripts/derive_nsc_pg_ctp_mode_jets.py --check
python3 -m pytest -q tests/test_nsc_pg_ctp_mode_jets.py
```
