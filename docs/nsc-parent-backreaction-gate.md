# Energy balance selects the evolving parent branch

The recorded outward Dirac power imposes a source-closure decision. Under
the state-independent, symmetry-preserving analytic vacuum branch used in
the working prescription, geometric terms cannot provide the compensating
radial Killing flux of an exactly static solution. This applies to its
nonlocal remainder as well as to the local terms already evaluated.

With the currently computed empty-incoming state retained, the next
constructive branch is an evolving geometry. A stationary alternative
would need a separately derived directed state, interaction or boundary
contribution. It cannot be obtained by fitting an unspecified geometric
pressure to the outgoing power.

This is a scoped flux obstruction, not a no-go theorem for all NSC
realizations. The [record](../results/development/parent-backreaction-gate.json)
uses the [existing quantum power](nsc-adm-neck-source-map.md),
[local source contribution](nsc-warp-local-neck-source.md), and
[mass/energy normalization](nsc-spherical-action.md). No radiation,
angular spectrum or collapse calculation is rerun.

## The symmetry includes the operator domain

For stationary real Euclidean fields on the noncompact radial line,

\[
ds_E^2=N^2d\tau^2+q^2(d\rho+b_Ed\tau)^2+r^2d\Omega^2,
\quad F_E=N^2+q^2b_E^2>0,
\]

the smooth shear

\[
dt_s=d\tau+\frac{q^2b_E}{F_E}d\rho
\]

gives ds_E²=F_E dt_s²+N²q² d rho²/F_E+r²dOmega². There is no radial
periodic time twist in this stated domain. The diagonal metric is
stationary and invariant under t_s to -t_s.

The compact projector must be preserved too. For the paired neutral
Dirac copies, the combined spin reflection and copy exchange
R_spin=tau1 tensor gamma^hat0 preserves
P=(I-eta tau3 tensor gamma5)/2 at both compact ends. The Dirac expression
changes by an overall sign and unitary transport under the reflection;
its D^dagger D modulus is unchanged. The same statement requires the
domain, reference subtraction and remaining background data to respect
that transformation.

For any differentiable covariant geometric functional with this symmetry,
its variational tensor is stationary and transforms as a tensor. Its
mixed component therefore satisfies

\[
T_{t_s\rho}=-T_{t_s\rho}=0,
\qquad\boxed{\mathcal P_{\rm geometric}=0.}
\]

Locality was not used. The identity continues to a symmetry-preserving
analytic Lorentzian branch wherever that continuation exists, so it does
not require calculating every term of the nonlocal remainder.

There is an important domain boundary to this argument. The ingoing PG
patch with its future horizon is not assumed to possess a global
Lorentzian time reflection: the corresponding static-coordinate map is
singular at the horizon. The argument instead uses the specified
Euclidean vacuum domain and its analytic continuation. An extension with
extra directed boundary/state data, or an unresolved determinant phase,
needs its own calculation. It is not silently included in the zero above.
Existence of the full nonperturbative continuation is not proved here.

## Apply the actual recorded flux

The canonical massless source gives

\[
\mathcal P_0=0.00014222067954246644>0.
\]

The local warp correction has zero Killing power, and the static
Einstein tensor requires zero total radial power. Thus an exactly static
match would require

\[
\mathcal P_{\rm other}=-\mathcal P_0.
\]

The symmetry-preserving geometric branch cannot supply that value.
Additional independent free channels with the same empty-incoming
condition cannot cancel the unchanged recorded channel either. Their
stationary outward powers have the form

\[
\mathcal P_\alpha=\int_{\omega_{\min,\alpha}}^\infty
\frac{d\omega}{2\pi}\;\omega\,\mathcal T_\alpha(\omega)
\frac1{e^{\omega/T_H}+1}\ge0,
\]

where the transmission probability is nonnegative and the channel label
includes the necessary degeneracies and particle/antiparticle counting.
This statement keeps the original channel/operator fixed; it is not a
bound on an interacting realization that changes that operator or its state.

The earlier incoming thermal control cancels net power only after changing
the incoming state. A recursive incoming source could be a legitimate
different construction, but it must be derived. No such flux is inserted
here. This closes the route of repairing the current static free-source
calculation with geometric coefficients alone.

## The required initial mass balance

Reuse the spherical mass normalization and energy balance already derived
from the common local action. In its asymptotically flat Einstein regime,
with fixed matching coefficient and reference charges,

\[
G=\frac1{16\pi A_{\rm EH}},\qquad
M_B=\frac{m_B}{G},\qquad
\frac{dM_B}{du}=-\mathcal P_\infty.
\]

Here u is the asymptotic retarded time with the inherited parent
normalization. For a radiating completion, M_B is the Bondi or
remaining-system mass; the energy carried by the outgoing radiation
remains in the full energy account. It is not a decrease of total ADM
energy including that radiation. The asymptotic distinction is standard
in the [Bondi–Sachs formulation](https://arxiv.org/abs/1609.01731).
The matching of the complete NSC theory to this asymptotic regime remains
an assumption of this application.

With a_EH=A_EH L_throat² and dimensionless mhat=m_B/L_throat,
uhat=u/L_throat, the recorded channel fixes the initial contribution

\[
\boxed{a_{\rm EH}\frac{d\hat m}{d\hat u}
=-\frac{\mathcal P_0}{16\pi}
=-2.829390519884\times10^{-6}.}
\]

No value of a_EH is fitted. For the leading outgoing asymptotic metric
coefficient g_uu=1-2m_B(u)/r+..., this gives

\[
a_{\rm EH}\hat r\,\partial_{\hat u}g_{uu}
=5.658781039768\times10^{-6}
\]

for the same known-channel contribution. Additional outgoing channels,
incoming energy or boundary work change the complete rate and must be
included when calculated.

This is an initial balance condition, not a full trajectory. The old power
is not extrapolated as a constant for all time, and no mass function is
inserted into the old geometry and presented as a solved field equation.
The diagonal source constraints and the remaining forces still have to
be satisfied on admissible initial data and throughout the evolution.

## Consequence for the work sequence

Keep the source calculations and geometric vertices already completed.
Use them to construct the evolving state/geometry problem, with its
initial constraints, work and boundary flux. Do not continue static
coefficient searches in the hope that a symmetry-preserving geometric
term supplies the missing power. Generic collapse and mass-loss laws
are imported tools; the NSC work is matching the same action and state
to their conditions and obtaining the actual finite evolution.

```sh
python -B scripts/derive_nsc_parent_backreaction.py --check
```

The command authenticates the existing records, checks the new
symmetry/domain map and applies their energy normalization. It compares
all fields and leaves the complete evolving parent-to-child solution,
physical stability and cosmological predictions open.
