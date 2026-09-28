# Progress toward a self-sourcing geometry

This continuation follows the user's requirement: derive the supporting stress
from the same Dirac/spectral system, without inserting an ad-hoc scalar field.
The verified 64-record preprint was pushed to GitHub; its preserved-history
checkout issue was fixed and both GitHub jobs passed at `aed05a5`.

## The precise feedback condition

The fundamental condition is the metric variation of the common action:

\[
\frac{\delta\Gamma_{\rm one}}{\delta g^{\mu\nu}}=0.
\]

After isolating the Einstein term already contained in that same action, the
condition can be written as an Einstein tensor matched to the remaining
renormalized quantum/boundary stress. The Einstein term must not be added a
second time with a new weight. All independent metric equations, the state,
and stress conservation must agree. A null projection alone tests the
defocusing mechanism but is not the entire self-sourcing condition.

## What is now calculated

| Part of the loop | New evidence | Scope |
|---|---|---|
| Propagation on the original horizon | [Tetrad derivation](nsc-dirac-tetrad.md) and [time evolution](nsc-lorentzian-transport.md) include the shift and spin connection, correct incoming/outgoing data and a conserved Dirac-norm ledger | Massless fermionic probe on the fixed four-dimensional benchmark; no vacuum or metric backreaction |
| Smooth geometry with a gap | [Smooth profile](nsc-smooth-geometry.md) removes the old seams while retaining a spatial band gap | Imposed static geometry, not a stationary solution |
| Stress required by that geometry | Independent curvature and lapse variations determine density and both pressures | Effective Einstein requirement; not assumed equal to quantum stress |
| Quantum source with the required sign | [Covariant-measure calculation](nsc-covariant-measure.md) computes compactification stress, all angular sectors and a finite proper-time cutoff | Constant-radius compact axial product; not a substitute for a varying neck or transverse compactification |
| Quantum response on the actual smooth neck | [Shape response](nsc-shape-response.md) differentiates the Dirac vacuum energy with respect to lapse, radial metric and sphere radius before evaluating the static gauge | Counterterm-independent periodic-minus-antiperiodic stress difference, not absolute vacuum stress |

No new scalar or negative-kinetic-energy field enters these calculations.
Quantum vacuum null stress and positive-norm fermionic propagation are
compatible; one must still calculate their gravitational response.

## Geometry and stress now share an explicit profile

Use the smooth profile

\[
r(x)=\sqrt{a^2+[\sin(kx)/k]^2},\qquad k=\pi/(2R),\quad a=1.
\]

Its throat has `r''(0)=1` and its periodic extension has no seams. At R=2,4,8,
the independently computed continuum band edges are approximately
0.76423204, 0.51425546 and 0.28454085.

The required radial-null Einstein source is

\[
8\pi G(\rho+p_x)_{\rm required}=-\frac{2r''}{r},\qquad
8\pi G\int(\rho+p_x)_{\rm required}\,dx
=-2\int(r'/r)^2dx<0.
\]

A cosmological-constant term has zero null projection and cannot provide this
part of the source. The remaining components and lapse constraint are also
retained explicitly in the result.

For the same varying radius, the quantum calculation uses

\[
H_\kappa=-i\sigma_2(c\partial_x+c'/2)+v\sigma_1,
\qquad c=N/q,\quad v=N\kappa/r.
\]

Full four-component angular counting and polar-factor Hellmann--Feynman
derivatives give the local stress change. On 512 axial points:

| R | E_P − E_AP | Neck Δ(ρ+p_x) | RMS conservation residual |
|---|---:|---:|---:|
| 2 | 0.15865448 | 0.01860824 | 4.44e-7 |
| 4 | 0.01945979 | 0.00179889 | 1.53e-7 |

The local trace/Weyl identity is near floating-point zero. Conservation
improves about fourfold per grid doubling, with angular truncations checked
independently. These figures are numerical approximations, not exact continuum
values or rigorous error bounds.

Switching P to AP lowers the null source at the neck by the displayed amount.
That is the required direction of change. However, neither G nor a baseline
stress has been fitted to equate it with `-2/a²`: Delta(P−AP) is not the absolute
AP expectation value. The energy has inverse-length units; the stress has
inverse-length-to-the-fourth units. The geometry/source comparison includes G.

## Why absolute closure is still open

Anomaly matching leaves a finite Weyl-invariant part undetermined. For example,
`alpha int sqrt(g) C²` leaves the four-dimensional Weyl anomaly unchanged but
changes a spacing variation. This is an uncertainty that the same ultraviolet
operator/measure must fix, not a proposed tuning parameter. On the constant
cylinder its null projection vanishes; on a varying neck this cancellation
cannot be assumed for absolute stress. It does cancel in the same-geometry
spin-structure difference above.

The written flat average over the common Weyl scale also diverges in the
current fixed-finite-matrix prescription. A physical dilaton measure or an
identified gauge volume is needed; arbitrary integration bounds would add
new inputs. These findings do not refute a fully specified ultraviolet
functional or the nested architecture.

The next absolute calculation must fix the covariant determinant's finite
part and complete compensator from the common operator, retain N,q,r and the
physical link through variation, and solve their full residuals together.
The boundary state and topology must be part of that solution. The compact
axial circle used in the vacuum calculation is not automatically the
unwrapped throat chain or the carrier's transverse compact direction.

The current result closes a previously missing *response calculation* in the
feedback loop. It does not yet close the self-consistent metric solution,
Lorentzian gravitational stability, particle identification or observations.
