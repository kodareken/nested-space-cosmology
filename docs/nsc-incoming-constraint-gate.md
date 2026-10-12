# Incoming-data incompatibility of the compact response family

The current compact PG pulse is a valid conditional response benchmark.
With its fixed incoming data it cannot be a smooth solution of the complete
local equations of the **retained 33-group free-Dirac realization**. This
necessary-condition result does not exclude the larger transmitting theory
or replace the earlier homogeneous no-interface certificate.

## Bound from the already owned source

Reuse the [compact source's Killing-power normalization](../src/recursive_horizons/nsc_compact_ctp_neck.py),
the [horizon occupation law](../src/recursive_horizons/nsc_paired_horizon_preparation.py)
and the [current-normalized mode sewing](../src/recursive_horizons/nsc_pg_massive_modes.py):

$$
P_{m,\lambda}=\frac{n_{\rm copy}d_\lambda}{\pi}
\int_m^\infty dE\,E\,\mathcal T(E)[f_H(E)-n_{\rm in}(E)],
\qquad 0\le\mathcal T\le1.
$$

The factor $1/\pi$ already includes the signed particle/antiparticle part;
the retained degeneracy includes angular signs. Neither is doubled again.
The normalized interior sewing matrix obeys

$$
S=\begin{pmatrix}0&1&0\\R&0&\sqrt{\mathcal T}\end{pmatrix}P,
\qquad \mathrm{Tr}(SC_{\rm source}S^\dagger)
=1+\mathcal T(n_{\rm in}-f_H).
$$

Its first two columns are orthogonal. Horizon coherence is retained in the
state but cancels from this current; it is not set to zero in the covariance.
Below the parent mass threshold, the incoming infinity channel is closed
and this Killing flux vanishes.

For massless angular modes, $\Omega>1$ makes $f_H-n_{\rm in}<0$. The LLL has
unit transmission and its already known contribution is

$$
P_{\rm LLL}=\frac{|q|\kappa^2}{48\pi}(1-\Omega^2).
$$

For every retained massive mode, irrespective of its detailed transmission,

$$
P_{m,\lambda}\le\frac{n_{\rm copy}d_\lambda}{\pi}
e^{-\beta_Hm}\left(\frac m{\beta_H}+\frac1{\beta_H^2}\right),
\qquad \beta_H=\frac{2\pi}{\kappa}.
$$

Applying this to the locked parameters gives
$P_{\rm LLL}\simeq-0.02227623165818447$. The massive multiplicity is512 at
each of $m=\pi/2$ and $m=\pi$, with positive-power bounds about
$1.03\times10^{-17}$ and $2.11\times10^{-35}$. Thus $P_K$ is strictly negative.
No scattering curve or stress integral needs recomputation for this sign.

The [record](../results/development/nsc-incoming-constraint-gate.json) also
contains an exact rational sign witness. It uses outward arithmetic bounds
$0.238<\kappa<0.239$, $3.97<\Omega<3.98$, $m>1.5$ and
$3.14159<\pi<3.14160$. These are enclosures of the unchanged inputs, not new
parameters. Together with $\beta_H>26$ and $e^{-39}<10^{-16}$ they give a
strictly positive rational margin without relying on a near-zero floating
point decision.

## Apply the bound at the actual incoming slice

The [common KS preparation](nsc-common-ks-trace.md) fixes its incoming Cauchy
slice at $\rho=1$. Here $r_1^2=2$ and $a_1^2=3\pi/2-4$. Reuse the existing
proper-frame projection of the conserved Killing power:

$$
T_{\hat0\hat1}(\rho=1)=-\frac{P_K}{4\pi r_1^2a_1^2}
\ge0.001244184174475\ldots>0.
$$

The more conservative exact rational witness exceeds $0.00123871$. This is
a **bound**, not a newly evaluated full stress tensor. The old neck's
numerical shift residual is not copied to this different surface.

`SuppliedKSHarmonicMetric` fixes all incoming metric jets to the homogeneous
Bronnikov KS data: its compact profile is flat as $\rho\to1^-$. The incoming
state is likewise unchanged. On these jets the geometric, magnetic,
parity-completed reference and local contributions have zero normal-axial
momentum. The new LLL allocation has $\Delta T_{01}\propto\partial_zH_a=0$.
The nonzero state momentum therefore violates the necessary shift constraint.

For a well-defined continuous renormalized source and a smooth metric,
continuity leaves a nonzero residual in an adjacent interior collar. A
surface term on the downstream smooth $\rho=0$ cut cannot cancel that bulk
residual. Consequently no member of this compact family can satisfy **all
pointwise history equations**, whatever its four amplitudes. A zero of four
projected action derivatives would not establish full stationarity.

## Correct physical continuation

Keep the pulse, kernels and reference calculations as reusable response
controls. Before a physical stationarity search, allow **constraint-compatible
incoming spherical Cauchy data**. In particular, the normal sphere expansion
$k_\perp(z)=n(r)/r$ may need spatial dependence, together with the normal jets
required by the retained higher-derivative action. A homogeneous shift change
alone cannot supply this physical geometric momentum.

The incoming quantum preparation must be transported or differentiated
consistently with that geometry while retaining its occupation law. Old
numerical covariance slots and stress tensors cannot be copied into a new
metric by assumption. Test the incoming constraints before launching another
complete source calculation aimed at stationarity in a fixed-endpoint pulse.

This result is conditional on the stated retained inventory and source law;
it is not a bound for an omitted infinite tower. The larger transmitting
EXISTENCE/NON-EXISTENCE gate remains OPEN. No coefficient, state, new boundary
term, metric evolution, PDF or public release is introduced.

```sh
python3 scripts/derive_nsc_incoming_constraint_gate.py --check
```
