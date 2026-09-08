# Dirac energy transfer and the cosmological source

The current smooth static Dirac vacuum supplies **zero net continuing energy
transfer** into either spatial half-domain. The same geometry-derived link
transports a prepared finite-energy excitation. This separates a working
transport channel from a state that actually supplies energy. It does not
establish the cosmological transfer Q or rule out transfer on an evolving
parent/child geometry.

The [new record](../results/nsc-6-energy-transfer.json) binds the existing
[smooth geometry](nsc-smooth-geometry.md),
[quantum stress difference](nsc-shape-response.md), and
[plateau requirements](nsc-plateau-conditions.md) without changing their bytes.
The new owner is `src/recursive_horizons/nsc_energy_transfer.py`;
`python3 -B scripts/check_nsc_energy_transfer.py --check` compares every field.

## State and energy accounting

Use the previously derived static Hamiltonian

\[
ds^2=N^2dt^2-q^2dx^2-r^2d\Omega_2^2,\qquad
H_\kappa=-i\sigma_2[c\partial_x+c'/2]+v\sigma_1,
\quad c=N/q,\quad v=N\kappa/r.
\]

Here x is the compact axial coordinate; N=q=1 and
r=sqrt(1+[sin(kx)/k]^2), k=pi/(2R), R=2 in the recorded experiment.
This is the compact static benchmark, not the trapped PG metric or a solved
expanding cosmology. The two radial spinor components are not two rooms.
Parent/child-like elimination below uses an actual spatial partition of H.
Its cross-region block B is taken from the original stencil, not inserted as
an adjustable Phi.

Let C_ij=<c_j† c_i> and P project onto the sites at x>=0. Assign half of every
crossing interaction link to each side:

\[
h_A=\tfrac12\{P,H\},\quad h_B=H-h_A,\quad
E_A=\operatorname{Tr}(C h_A),\quad \dot C=-i[H,C].
\]

Then, in units hbar=c=a_throat=1,

\[
\boxed{\dot E_A=\operatorname{Tr}(C J_A)+
                  \operatorname{Tr}(C\dot h_A),\qquad
J_A=i[H,h_A]=\tfrac i2[H^2,P].}
\]

Positive J_A is energy entering A. For static H, J_A+J_B=0. Alternatively,
use bare energies PHP and (1-P)H(1-P) and keep the complete link energy
separately. Omitting link storage would give an inconsistent balance. The
half-link prescription fixes one regional allocation; total conservation
does not uniquely identify a physical matter/dark split.

For the global static vacuum C=theta(-H), [H,C]=0, so cyclicity gives

\[
\boxed{\langle J_A\rangle
=i\operatorname{Tr}([C,H]h_A)=0.}
\]

This is exact at every finite resolution, for either periodic or antiperiodic
spin structure. Any stationary spectral occupation f(H) has the same net
identity. The numerical spectral weighting in the record is a stationarity
control, not a new covariant regulator or an absolute stress prescription.
All angular multiplicities multiply zero; no angular truncation can turn
this particular equilibrium result into an injection rate.

More generally, stationarity allows circulating currents with zero regional
net gain. In this time-reversal-invariant static P/AP vacuum, both recorded
cut currents vanish separately. The continuum static vacuum and a
symmetry-preserving covariant stress prescription likewise cannot generate
a time-odd radial flux from this background. Absolute energy density and
pressure on the varying throat remain uncomputed. In particular, the
previous P−AP stress difference is stored stress, not power.

## Evaluated response and transport

Six static cases (P/AP, kappa=1,2,4) have nonzero Schur self-energy and
cross-region vacuum correlations, but zero energy current. Their joined and
reduced inverses agree to floating-point accuracy at z=0.3+0.4i.

For a transport control, prepare one normalized positive-energy mode above
the AP vacuum by projecting a smooth bump centered at x=-0.6, width 0.35,
momentum 4, and positive sigma2 seed. Projection is nonlocal; the resulting
state is not claimed to have compact support. The excitation is an explicit
state input, not vacuum particle creation or a predicted particle identity.

At t=1.2, in units 1/a_throat:

| Axial grid | Total excitation energy | Energy gained by x>=0 | Throat contribution | Compact-seam contribution |
|---|---:|---:|---:|---:|
| 64 | 5.76336922 | 5.65332623 | 5.65230009 | 0.00102614 |
| 128 | 5.82712077 | 5.75723400 | 5.75609973 | 0.00113427 |
| 256 | 5.84334389 | 5.77937518 | 5.77814481 | 0.00123037 |
| 512 | 5.84744708 | 5.78342676 | 5.78218434 | 0.00124242 |

The compact region has two interfaces, both included. The last change in
integrated transferred energy is about 0.00405; the table demonstrates
refinement, not an interval-certified continuum error bar. The sharp lattice
partition approaches the throat as spacing vanishes. Local Dirac energy is
not a probability and can be negative even in a positive-energy state.

Total-energy drift is below 4e-14 in this run. The independent sparse matrix
exponential agrees with spectral evolution within 5e-15; time quadrature
agrees with regional energy change within 6e-11. Tests also compare the
current with finite time differentiation, check Pauli occupancy, retain
both interaction-energy allocations, and reproduce an analytic two-level
exchange/recurrence control. These numerical residuals describe finite
algebra and time integration, not the spatial continuum uncertainty.

The one-particle Hamiltonian and retarded Schur map are identical for the
vacuum and the excited state. Yet at grid 64 and t=0.6 the excited state
has inward regional power about 11.1937 in units 1/a_throat^2, while the
vacuum has zero. Thus even a completely evaluated retarded boundary map
does not, by itself, determine an expectation value of energy current.

## The missing part of boundary elimination is explicit

Write the spatial blocks as H_AA, H_BB and B. Eliminating psi_B in time gives

\[
\begin{aligned}
i\dot\psi_A(t)={}&H_{AA}\psi_A(t)
-i\int_0^t B e^{-iH_{BB}(t-s)} B^\dagger\psi_A(s)\,ds\\
&+B e^{-iH_{BB}t}\psi_B(0).
\end{aligned}
\]

The first integral is retarded memory. The last term carries the eliminated
initial state. Its covariance contains

\[
B e^{-iH_{BB}t} C_{BB}(0)e^{iH_{BB}s}B^\dagger,
\]

and initial C_AB correlations must also be retained. Dropping this data is
not a justified vacuum choice. The existing common-energy Schur recursion
specifies the one-particle retarded part; an in-in/occupation construction
must accompany it to determine transport. Tests directly verify the
initial-source term through a block resolvent and show that omitting it
changes the solution.

This distinction is established nonequilibrium Green-function physics;
see [Jauho, Wingreen and Meir](https://arxiv.org/abs/cond-mat/9404027).
The initial density matrix also appears explicitly in the gravitational
influence-functional construction of
[Martín and Verdaguer](https://arxiv.org/abs/gr-qc/9904021), which studies scalar
fields and is methodological prior art, not a Dirac solution for this throat.
Dirac stress must obey its covariant renormalization requirements; see
[Dappiaggi, Hack and Pinamonti](https://arxiv.org/abs/0904.0612).

## Conversion into a cosmological Q

Use signature +---, Q positive into visible dust, and include every
interaction contribution in T_b+T_D:

\[
\nabla_\mu T_b^{\mu\nu}=I^\nu,\quad
\nabla_\mu T_D^{\mu\nu}=-I^\nu,\quad
I^\nu=Q u^\nu+F^\nu,\quad u_\nu F^\nu=0.
\]

For any vector xi, the exact work identity is

\[
\nabla_\mu(T^{\mu\nu}\xi_\nu)
=I^\nu\xi_\nu+\tfrac12T^{\mu\nu}\mathcal L_\xi g_{\mu\nu}.
\]

On the static metric xi=partial_t=N u, with local outward flux
F_out=-T^{mu nu}u_mu s_nu and outward unit spacelike s, it gives

\[
E_{\xi,A}=\int_A N\rho\,dV,\quad dV=4\pi q r^2 dx,\qquad
\frac{dE_{\xi,A}}{dt}=-\oint_{\partial A}N^2F_{\rm out}\,dA
                         +\int_A N^2Q\,dV.
\]

One lapse redshifts energy and the other converts proper time. For a
stationary receiver with uniform lapse N_c, inward proper power is inward
Killing power divided by N_c^2. Only a derived matching to the cosmological
receiver and its proper volume V_c could then identify Q_b=P_deposited/V_c.
Moving boundaries or nonstationary geometry require their work terms.
A spatial-window Ward source and its equivalent surface flux are two
representations of one contribution; do not count them twice.

For the flat constant-G late-time dust/dark baseline, a constant 0<f*<1 gives

\[
\boxed{\frac{Q_*}{H\rho_b}=1-2q_{\rm dec};
\qquad q_{\rm dec}<0\Longrightarrow\frac{Q_*}{H\rho_b}>1.}
\]

If the illustrative f=0.95, q_dec=-0.55 were already a plateau, this would
require epsilon=0.105 and Q/(H rho_b)=2.1. Those numbers use specified
benchmark inputs; neither is an NSC prediction or a proposed source term.

This exchange preserves total conservation. Moreover, for fixed baryonic
mass and conserved comoving baryon number, rho_b=m_b n_b forces Q_b=0.
Massless Dirac energy arriving through a throat could supply radiation or
heat instead. A mechanism producing the specified ordinary dust component
is therefore a necessary physical dependency. A spatial partition of the
free Dirac field has not identified that component.

## Result and next controlled realization

The tested static vacuum supplies no persistent injection, so it cannot
provide the active-transfer branch of an accelerating finite dark plateau.
This does not prevent a stationary geometry from being self-sourced by
vacuum stress, which is a different tensor-matching question.

The next calculation must evolve a specified initial state with the
time-dependent geometry, carry retarded memory and occupation data together,
and evaluate stress, flux and geometric work from the same functional.
A maintained incoming state, evolving throat or collapse may supply energy;
its source cannot be inferred from a nonzero Phi or inserted to satisfy the
plateau identity. In a fixed finite closed system, bounded regional energy
also implies time-averaged net power tends to zero. Infinite-volume and
long-time limits need a separate analysis before invoking steady transport.

This result adds an executable energy ledger and identifies the missing
state and cosmological conversion equations. It is new within this project;
the equilibrium and transport identities are established physics. Absolute
metric-source matching, physical scale selection and full cosmological
closure remain open.

The next controlled step has now been computed in the
[geometry-pulse vacuum experiment](nsc-vacuum-work.md): changing the radius
produces pairs and its supplied work accounts for their energy. That result
does not change the static-equilibrium conclusion or supply a self-consistent
cosmological source.
