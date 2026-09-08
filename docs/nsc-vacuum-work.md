# Geometry does work on the Dirac vacuum

A smooth pulse of the existing sphere-radius geometry produces Dirac
particle–antiparticle excitations from its initial vacuum. The effect is
computed from the same Dirac operator and requires no added scalar matter
field. Its final excitation energy equals work supplied by the prescribed
geometry. A uniform change of clock produces no excitations.

This is a controlled positive response, not yet a self-sourcing solution:
the pulse is an input. The common action must still supply its motion and
backreaction, and relate the resulting state to cosmological matter and Q.

The [record](../results/nsc-6-vacuum-work.json) and
`scripts/check_nsc_vacuum_work.py --check` preserve every scientific field.
The calculation follows the [static energy-transfer result](nsc-energy-transfer.md).

## A time-dependent geometry with a fixed canonical domain

Take the same smooth static axial profile r0(x), now deforming only the
sphere radius:

\[
ds^2=dt^2-dx^2-r(t,x)^2d\Omega_2^2,\qquad
r(t,x)=\frac{r_0(x)}{1+\epsilon s(x)g(t)}.
\]

The periodic spatial profile s(x)=exp(2[cos(2pi(x+0.6)/L)-1]) is smooth,
0<s<=1, and g(t)=exp[-t²/(2tau²)]. The recorded epsilon values keep the
denominator positive. The axial spin structure remains AP. There are no
metric seams or new matter fields.

The tetrad Dirac equation contains both radius connections:

\[
i\gamma^0(\partial_t+r_t/r)\psi+
i\gamma^1(\partial_x+r_x/r)\psi+
\frac{i}{r}\slashed D_{S^2}\psi=0.
\]

The canonical field u=r psi removes both connections and the r² Hilbert
measure. Angular separation therefore gives exactly

\[
\boxed{H_\kappa(t)=-i\sigma_2\partial_x+
       \frac{\kappa}{r(t,x)}\sigma_1
       =H_{0,\kappa}+\epsilon g(t)V_\kappa,\qquad
V_\kappa=\frac{\kappa s(x)}{r_0(x)}\sigma_1.}
\]

This affine relation is also exact in the existing staggered stencil.
Finite geometric interventions check it directly. The potential is an
angular Dirac connection inherited from geometry; it is not an independently
chosen scalar field. The metric pulse changes the sphere relative to the
axial direction, so it is not a conformal rescaling of the entire spacetime.

## Vacuum transitions and supplied energy

Start with the static filled negative spectrum of H0. In leading perturbation
theory, a negative state m can mix into a positive state n:

\[
P_{nm}^{(2)}=2\pi\epsilon^2\tau^2|V_{nm}|^2
             e^{-(E_n-E_m)^2\tau^2},\qquad E_n>0>E_m.
\]

Consequently,

\[
\boxed{N_{\rm pairs}^{(2)}=\sum_{n+,m-}P_{nm}^{(2)},\qquad
\Delta E^{(2)}=\sum_{n+,m-}(E_n-E_m)P_{nm}^{(2)}\ge0.}
\]

The gap counts both the particle and its negative-sea hole. Do not add
another factor two to the energy. Full four-component angular counting
multiplies each representative radial channel by 4 kappa. The generic
omitted order is epsilon³; averaging opposite pulse signs removes odd terms.

At tau=0.5, in units hbar=c=a_throat=1, the coefficient of epsilon² in the
full excitation energy through kappa=8 is:

| Axial grid | Energy coefficient |
|---|---:|
| 32 | 0.61744877 |
| 64 | 0.61299047 |
| 128 | 0.61188656 |
| 256 | 0.61161123 |

The last spatial change is about 0.0002753. At grid 128, adding kappa=7 and 8
after kappa=6 changes the coefficient by about 3.27e-10. These separately
demonstrate numerical convergence; an infinite angular remainder and a
continuum interval enclosure have not been proved.

One-channel kappa=1 coefficients at grid 128 decrease from about 0.194633
for tau=0.25 to 0.0875495 for tau=0.5, 0.00662354 for tau=1, and 1.21e-8
for tau=2. Slow deformation suppresses transitions in this gapped compact
benchmark. This is not a general statement about the gapless isolated throat.

An independent finite real-time calculation evolves all occupied modes in
the interaction picture and integrates external work. For kappa=1, grid 32,
epsilon=0.02, tau=0.5, tightened integration gives:

- Pair number: approximately 1.4018114e-5.
- Final excitation energy: approximately 3.4664005e-5.
- Supplied work: approximately 3.4664005e-5.
- Work/energy residual: below 6e-13; occupied-mode orthonormality residual
  below 6e-15.

These are one representative channel, not the full angular sum in the
preceding table. Opposite-sign pulse averages at epsilon=0.04,0.02,0.01
approach the independently calculated leading coefficient with the expected
even-power scaling. A constant-cylinder momentum calculation independently
checks both transition matrix elements and Gaussian normalization.

The same evolved coherent vacuum also gives a nonzero regional transfer.
For epsilon=0.02, the region x>=0 receives net integrated inflow about
1.05945e-6 and direct geometric work about 4.32787e-6 by t=3.5. Their sum
matches its energy change, about 5.38732e-6, with residual below 1e-11 in
the tightened integration. The net inflow includes both compact-region
boundaries. It is not a separately measured throat-only flux or cosmological
Q. For a smaller positive pulse, the recorded regional energy change can
even be negative because coherent polarization redistributes local energy;
total excitation energy remains positive. There is no universal positive
injection law inferred from the sign of a radius perturbation.

The real-time calculation truncates the pulse at t=±7tau. Its small endpoint
Hamiltonian contribution is included. A finite-matrix Duhamel bound on the
interaction-picture propagators is

\[
\|U_{\rm full}-U_{\rm truncated}\|
\le |\epsilon|\|V\|\sqrt{2\pi}\tau\,
   \operatorname{erfc}(7/\sqrt2).
\]

For the cited epsilon=0.02 case the bound is about 5.20e-14. This is a
finite-operator pulse-tail estimate, separate from numerical integration and
spatial/angular errors.

## Clock, work and physical interpretation

The uniform lapse control has H(t)=N(t)H0 and
U=exp[-iH0 integral N dt]. It changes phases but never mixes positive and
negative modes. The computed pair coefficient vanishes to numerical
precision, as the exact identity requires. The radius pulse passes this
control because it changes the geometry, not just the clock.

The energy balance is

\[
\boxed{\frac{d\langle H\rangle}{dt}=\operatorname{Tr}(C\dot H)
       =-2\int \frac{r_t}{r}\,p_\perp\,dV,
\qquad dV=4\pi r^2dx.}
\]

The second equality is the covariant Ward identity on this closed geometry,
with the corresponding regulated/renormalized stress. The numerical test
checks the finite canonical work identity; it does not calculate an absolute
renormalized p_perp during the pulse. At identical asymptotic metrics, common
local geometric counterterms cancel from final excited-minus-vacuum energy.
This does not fix the absolute vacuum stress needed by the metric equations.

Local flux can include coherent vacuum-polarization terms at order epsilon.
The diagonal final pair probabilities alone therefore cannot reconstruct a
time-dependent local current; the evolved covariance is required. Likewise,
generic produced Dirac excitations do not establish baryonic dust or its Q.

Particle creation by evolving geometry is established prior physics, for
example [Chung et al. on massive Dirac fields in inflation](https://arxiv.org/abs/1109.2524).
Their FLRW problem differs from this anisotropic radius pulse. The new
project contribution is the explicit geometry vertex, spectrum, production
and work check on our smooth operator, not a new-to-world claim for the
general mechanism.

The remaining feedback equation must make this geometry dynamical, with
the same common-action metric source and quantum state. Until it does,
the result is conversion of supplied geometric work into Dirac excitations,
not a sustained self-sourced dark-sector plateau.

## Integrated verification

The laboratory gate passes 131 tests and reproduces all 16 active result
records. The archive verifier preserves all 1,731 historical tracked-file
hashes. Both new records compare every scientific field with the tolerances
stated by their reproduction owners.

An independent derivation checked the temporal spin connection, angular
counting, Gaussian transitions and covariant work identity. A separate
Schrödinger-picture calculation at 16 axial points agreed with the
interaction-picture integrated regional inflow within 1.2e-13 and regional
work within 5e-16. These review checks are additional finite controls, not
continuum bounds. A requested Grok worker returned malformed output after
its timeout; it supplied no accepted scientific evidence. The completed
independent derivation and numerical review were performed by a separate
Codex agent.
