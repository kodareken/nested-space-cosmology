# Incoming charged source: finite spectral panels and matched local allocation

On the unchanged incoming KS surface at $\rho=1$, the 32 non-LLL channel
groups now supply all four state-minus-reference moments from their actual
global source columns. The result uses the existing horizon covariance,
incoming occupation, mode maps and fourth-order reference. It does not use a
seed covariance as a spatial state and does not evolve the metric.

The channel contribution is

$$
\frac{n_{\mathrm{copies}}d_j}{4\pi^2r_1^2a_1n_{\mathrm{signs}}}
\sum_{\eta}\int_0^{E_{\max,j}}dE\,
\mathrm{Tr}\big[(C_j-P_j^{(0\ldots4)})\,V_A\big].
$$

Here $r_1^2=2$, $a_1^2=3\pi/2-4$, and the incoming canonical momentum is
$k=-E$ on this homogeneous intrinsic surface. The degeneracy already includes
both angular signs. Negative-frequency folding supplies the factor two in
this prefactor; it is not added again. The vertices in the record's order are
the normal Hamiltonian, $(k/a_1)\sigma_3$, $-(k/a_1)I$, and
$\lambda\sigma_2/(2r_1)$.

The new finite sum selects existing refined low/middle panels where available.
Those replace the corresponding base panel. The old retained **packet** record
used generic base panels and treated their refinements as controls; this new
observable's choice does not rewrite that record. Group 14's low32 grid
replaces only $(0.25,1)$ and $(m,8)$ in low16. Group 13 uses its original low
grid and fast24 middle grid. All 38 compact signed families retain their eight
physical real subgap nodes on $(1,m)$. Integrating the selected weights gives
exactly $E_{\max,j}$ for every signed family: no panel is counted twice.

## Stable source arithmetic

The incoming source columns obey the existing physical sewing relations:
the partner column has unit norm, is orthogonal to the outer and incoming
columns, and the latter two squared norms sum to one. With horizon-pair trace
one and no incoming cross correlation, the current can be evaluated as

$$
T_{01,j}(E)=\frac{E}{a_1}\|F_{\mathrm{in}}\|^2(n_{\mathrm{in}}-f_H).
$$

The implementation verifies those relations before using the identity. A
generic coisometry is insufficient and is rejected. Source columns, covariance
and reference are never normalized, clipped or replaced. Raw matrix currents
remain in the payload as roundoff controls. This matters for compact modes:
their physical thermal current is much smaller than the error in subtracting
two unit traces. Energy and pressure kernels retain the direct contraction.

For the middle panels the existing order-16 Riccati boundary mode is exposed
directly. The vacuum Bloch subtraction is rationalized before subtracting
higher adiabatic orders; occupations are read without computing `1-(1-f)`.
This improves arithmetic, not the underlying mode approximation. Its inherited
equation/order diagnostics remain in the record and are not relabeled as
newly passed source-tail bounds.

## Matched allocation, counted once

The [charged light restoration](nsc-magnetic-light-restoration.md) uses the
actual magnetic angular spectrum. Add its **nonzero-angular** restoration,
the existing physical-metric LLL geometry, and the LLL occupation contribution
once each. Equivalently the new general light-restoration action already
includes the geometric LLL term, so that action receives only the additional
LLL state term. These are alternative accounting routes, never cumulative.
The compact induced Einstein, magnetic and curvature terms remain outside
this quantum-source sum, on the local-action side of the equations.

At the current finite resolution the quantum-source approximant in the order
$(\rho,p_\parallel,T_{01},p_\perp)$ is approximately

$$
(0.01225949982,\;-0.00338653087,\;0.00319035235,\;0.00188597873).
$$

These are finite-source values on the incoming surface, not a new neck tensor,
constraint solution or physical history. The JSON retains every channel and
the separate allocation contributions.

## Remaining numerical owner

The measured summed absolute base/refined changes are about $1.95\times10^{-10}$
in density and $1.85\times10^{-10}$ in parallel pressure. The current changes
by less than $3\times10^{-18}$. These are measured quadrature differences,
not rigorous error bounds. The eight-node subgap contribution has not yet
received an independent source refinement. The omitted energy tail and the
retained angular/compact truncation must also be assessed for this observable;
old packet-covariance accuracy does not bound a stress integral. The arithmetic
indicators likewise do not bound upstream modal error.

The next computation therefore targets the subgap source and asymptotic source
remainder, followed by the same incoming source/reference/local balance on
the [permitted normal Cauchy data](nsc-incoming-cauchy-data.md). No optimizer
or metric timestep is started from the known incompatible frozen data.

Reproduce the record from the authenticated spectral payload:

```sh
python3 scripts/derive_nsc_incoming_spectral_source.py --check
```

For an intentional change to the numerical source producer, its explicit
`--prepare-panels` path regenerates only the new observable from archived mode
inputs. It never runs the old horizon, scattering or packet generators.

The magnetic flux, $A$, $\Omega$, $\zeta$, $V_{\mathrm{full}}$, original seeds,
and all preceding certificates remain unchanged. GitHub and PDF are unchanged.
