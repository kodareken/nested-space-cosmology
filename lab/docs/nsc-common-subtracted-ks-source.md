# Common equal-KS-time state/reference pairing

The [kernel owner](../src/recursive_horizons/nsc_common_subtracted_ks_source.py)
puts the physical Gaussian state and the fourth-order reference on the same
spatial kernel before metric contraction. It imports the existing canonical
operator, reference band action and actual raw-KS vertices. It neither
identifies incoming frequency with instantaneous momentum nor introduces a
new physical cutoff.

For fixed KS time let $z_\pm=Z\pm\eta/2$. In the fixed reference chart these
are events with the **same** $\rho$ and PG times $\tau\pm\eta/2$. Resolved
fields, already including their source phases, give

$$
C_\eta=\int\frac{dE}{2\pi}\,
\Psi_E(z_+)C_{\rm source}(E)\Psi_E(z_-)^\dagger,
$$

$$
\partial_\eta C_\eta=\frac12\int\frac{dE}{2\pi}
\left[\Psi_{E,z}(z_+)C_{\rm source}\Psi_E(z_-)^\dagger
-\Psi_E(z_+)C_{\rm source}\Psi_{E,z}(z_-)^\dagger\right].
$$

The same separation is used for the reference Weyl kernel,

$$
P_\eta=\int\frac{dk}{2\pi}e^{ik\eta}P_{\rm ad4}(T,Z,k),
\qquad D_\eta=C_\eta-P_\eta.
$$

Incoming $E$ and local $k$ are separate integration variables. Their finite
numerical extents do not become a common cutoff by assigning them equal
values. The inherited policy is subtraction followed by controlled spectral
and coincidence limits, as in the existing compact completion. No claim of
equality with the entire raw Euclidean heat-vacuum functional is added.

## Metric pairing and reference remainder

For the already owned first-order vertex
$\delta_BH=\mathrm{Op}^{W}(M_B+V_Bk)$, the integrated symmetric pairing
has bulk density

$$
J_{B,\eta}(D)=\mathrm{Re}\int dT\,dZ\,
\mathrm{Tr}\left(\overline M_BD_\eta
-i\overline V_B\partial_\eta D_\eta\right),
$$

where bars average coefficients at the two endpoints. The API requires a
constant spatial split and compact boundary-identical vertices; otherwise
the spatial endpoint term must be supplied and the call is rejected.
The measure remains $dT\,dZ=d\tau\,d\rho/a_0$.

The subtraction action is not reduced to its projector vertex alone. The
existing band-action identity supplies

$$
\mathcal B_{B,\eta}=\mathrm{Re}\int dT\,dZ\,\frac{dk}{2\pi}
e^{ik\eta}\left[
\delta_B\mathrm{Tr}(\Pi h)
-\frac12\mathrm{Tr}(P\star\delta_BH+\delta_BH\star P)
\right]_{\le4}.
$$

The code retains its temporal, star-exchange and antisymmetric-exchange
pieces explicitly. Compact spacetime support does not discard finite
momentum-boundary terms. The combined action derivative is

$$
\delta_B\Gamma=\lim_{\eta\to0}\left[-J_{B,\eta}(D)
+\mathcal B_{B,\eta}\right]
+\delta_BS_{\rm LLL,geo}+\delta_BS_{\rm local,locked}.
$$

The [LLL allocation](nsc-lll-geometric-history.md) and
[induced local terms](nsc-spherical-local-history.md) enter once. Force is
minus this action derivative; no sphere-area or channel multiplicity is
implicitly added by the kernel API.

## Implementation scope and next input

Three focused tests verify the coincidence sign against `weak_KS_source`,
finite-separation adjoints and endpoint coefficients, and the exact retained
band-action remainder. These are implementation checks, not a full stress
calculation. `SpectrumDeclaration` makes input coverage explicit;
`require_complete_trace` rejects control or partial spectra. A declaration
is an evidence claim, not proof supplied by a finite array.

The [spectral field recovery](nsc-pg-spectral-mode-recovery.md) now supplies
many previously omitted real-frequency bulk inputs. They still need
conditional evolution and the common pairing across the varied support.
The earlier four-frequency response sample cannot provide the vacuum/CAR
singularity needed for a complete ultraviolet subtraction. Integrated packet
subgap and tail matrices likewise cannot replace their differentiated source
observables. The full source, physical endpoint equations and stationarity
remain OPEN; no finite total stress or metric step is asserted.

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_common_subtracted_ks_source.py
```
