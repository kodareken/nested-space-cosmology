# The evolved scalar concentration as a lens

The [finite scalar/gravity evolution](nsc-atlas-core-evolution.md) now
supplies nonradial light bending and beam distortion. Photons cross the
actual time-dependent concentration, reservoir and ambient field. No
fitted lens profile, static metric, extra matter species or new coupling
is introduced. The streamed final fields are required to match the
previous evolution's saved hashes exactly at each retained resolution.

This is a finite-screen optical experiment. Its source, emitting angle,
proper age and receiving sphere are defined physically. It supplies a
conditional lens map and a new way to check the same source law; it is
not an astronomical catalogue fit or a prediction of the number of lenses.

## Why this connection comes next

The rethinking compared a population/preparation law, nonspherical
persistence and lensing. Preparation still needs a physical state-selection
input. Nonspherical perturbations can be defined on the evolving background,
but require a new coupled perturbation system. Lensing uses quantities
already present and tests an additional observable immediately.

The exact lens-map approach is attributed to
[Frittelli and Newman](https://arxiv.org/html/gr-qc/9810017) and the
[Frittelli--Kling--Newman example](https://arxiv.org/html/gr-qc/0001037).
Their static worked example is not imposed on this evolving metric. The
ray, affine normalization, screen projection and curvature equations are
derived for the existing NSC field and checked independently.

The optical campaign uses a finite fan at two proper ages, three spatial/time
resolutions, a separate optical-step comparison, a larger outer chart and
an independent homogeneous control. The admitted budget is4500 CPU seconds
including six workers. Compact completed checkpoints are retained before
validation or budget reporting. Full fields remain regenerable local data.

## Physical source and receiving screen

A source and a spherical family of receivers are held at areal radius
Rs=.6. Their clocks are zero on the initial slice; spherical symmetry makes
their proper-time rates identical. The source emits at proper ages2.8 and3.4,
bracketing different phases of the already computed contraction and breathing.
Its inward rest-frame sky angle is psi, with unit emitted energy fixing the
affine normalization. The primary declared angles are .02, .04 and .08 radians.
A wider finite fan examines the geometry between them and nearer the axis.

In the retained metric,

\[
ds^2=-\alpha^2dt^2+
\frac{(dR-\alpha v\,dt)^2}{B}+R^2d\Omega^2,
\qquad v=H(t)R,\qquad \beta=\sqrt{B},
\]

a held observer has velocity w=-v/beta relative to the slice normal,
Lorentz factor gamma=(1-w^2)^(-1/2), and clock rate alpha/gamma.
For mu_e=-cos(psi), the emitted normal-frame quantities are

\[
E_e=\gamma_e(1+w_e\mu_e),\qquad
q_e=\gamma_e(\mu_e+w_e),\qquad
L=R_s\sin\psi,\qquad p_{R,e}=q_e/\beta_e.
\]

L is conserved angular momentum. The source's measured energy remains
one after the boost. At the first later outward recrossing of R=Rs,

\[
\nu_o=\gamma_o(E_o+v_op_{R,o}),\qquad
\mu_o=\frac{q_o+v_oE_o/\beta_o}{E_o+v_op_{R,o}},
\qquad q_o=\beta_op_{R,o}.
\]

The arrival event is solved separately for every ray. Holding coordinate
time fixed while taking angular derivatives would give a different,
incorrect screen map. The receiving worldline is part of the output:
comparisons below match source proper age, launch angle and screen radius,
not an assumed common receiving worldline in the two spacetimes.

## Nonradial null evolution

For an equatorial null ray, define

\[
E=\sqrt{Bp_R^2+L^2/R^2},\qquad
\mathcal H=-p_t=\alpha(E+vp_R).
\]

Hamilton's equations with coordinate time as parameter are

\[
\dot R=\alpha\left(v+\frac{Bp_R}{E}\right),\qquad
\dot\varphi=\frac{\alpha L}{ER^2},
\]

\[
\dot p_R=-\alpha_R(E+vp_R)
-\alpha\left(\frac{B_Rp_R^2-2L^2/R^3}{2E}+Hp_R\right).
\]

All metric coefficients evolve. A second frequency equation transports
normal energy independently, with mu=q/E:

\[
\frac{d\log E}{dt}
=\alpha\left[K_r\mu^2-H(1-\mu^2)\right]
-\mu\beta\alpha_R.
\]

It is compared with E reconstructed from the canonical momentum. In the
outward radial limit, this returns the previous clock equation. The affine
rate is d lambda/dt=alpha/E.

## Vacuum, scalar excitation and Weyl focusing

The canonical Einstein equation in the computed parent is

\[
R_{ab}=\kappa\left(\nabla_a\phi\nabla_b\phi+Ug_{ab}\right).
\]

For any null ray, therefore,

\[
R_{ab}k^ak^b=\kappa(k\phi)^2
=\kappa(E\Pi+q\beta F)^2.
\]

This exact contraction removes the direct potential term, including the
constant inherited vacuum. It does not remove the potential's effects on
the field, expansion, light path or redshift. It is not a statement that
changing Lambda leaves a lens experiment unchanged.

For the spherical scalar, define the Weyl curvature scalar

\[
W=\frac{m}{R^3}-\frac{\kappa}{6}
\left(U+\frac{\Pi^2-BF^2}{2}\right).
\]

The ray's orbital plane fixes two principal screen directions. Their
optical tidal eigenvalues are

\[
T_\perp=\frac{\kappa}{2}(E\Pi+q\beta F)^2
+\frac{3WL^2}{R^2},\qquad
T_\parallel=\frac{\kappa}{2}(E\Pi+q\beta F)^2
-\frac{3WL^2}{R^2}.
\]

The perpendicular direction is normal to the orbital plane. The Ricci
term focuses both directions; Weyl curvature stretches them differently.
The factors and signs follow from the spherical warped-product curvature.
An exact symbolic check substitutes the Einstein areal-radius Hessians
and the null relation E^2=q^2+L^2/R^2. Static vacuum and Minkowski
calculations independently verify the implemented signs.

The optical/Jacobi framework is attributed to
[Seitz, Schneider and Ehlers](https://arxiv.org/abs/astro-ph/9403056)
and its Sachs-equation foundation. The result used here is the explicit
reduction on the already evolved scalar source, not a claim that the
general focusing law or the absence of direct vacuum null focusing is new.

In the homogeneous parent, F=0 and m=kappa rho R^3/6, hence W=0.
The two eigenvalues coincide. At a pure vacuum extremum, Pi=F=0, both
vanish. In the inhomogeneous solution the scalar gradient, momentum and
Weyl curvature are calculated from the same sourced metric used by the clock.
No separate lens mass is adjusted to reproduce the optical result.

## Two independent computations of the beam

The signed principal Jacobi scales satisfy

\[
\frac{d^2D_i}{d\lambda^2}=-T_iD_i,
\qquad D_i(0)=0,\qquad
\frac{dD_i}{d\lambda}(0)=1.
\]

They are integrated through the evolving curvature. Separately, the
landing ray map supplies the exact screen identities

\[
D_\perp=\frac{R_s\sin\varphi_o}{\sin\psi},\qquad
D_\parallel=-R_s\mu_o\frac{d\varphi_o}{d\psi}.
\]

The derivative is of the full landing map, including the angular dependence
of arrival time. Projection into the receiver's ray-orthogonal screen
removes the temporal component and supplies mu_o. The sign convention
gives both D_i=2Rs cos(psi) for a Minkowski chord. Image parity is compared
relative to the homogeneous determinant, avoiding a convention-dependent
sign from coordinates on the antipodal receiving sphere.

The angular derivative uses two symmetric angle steps and fourth-order
Richardson combination at each primary angle. The curvature result does
not use that derivative. These are independent calculations of the same
optical observable, not two independent physical detections.

The signed beam area per emitted solid angle is D_perp D_parallel.
For a fixed isotropic emitted luminosity, photon conservation gives

\[
F_o=\frac{L_e}{4\pi}
\frac{\nu_o^2}{|D_\perp D_\parallel|}.
\]

One frequency factor is photon energy and the other is reception rate;
the source-to-observer area distance already has its own observer
normalization. A relative beam flux compares this expression with its
homogeneous counterpart at the same emission protocol. It is not yet a
fixed-source/fixed-receiving-worldline astronomical image magnification.

Zeros of a Jacobi scale are optical caustics. A geometric-optics zero does
not predict infinite physical brightness: finite source size and wavelength
would need a further calculation. A grazing screen crossing is a separate
failure of the landing-time coordinate, not automatically a caustic.

## Independent homogeneous light cones

In the closed reference geometry, write R=a sin(theta) and use conformal
time eta. The photon follows a great circle of the unit three-sphere.
Let sigma=eta-eta_e, mu_e^(n)=q_e/E_e and s_e^(n)=sin(psi)/E_e. Then

\[
\cos\theta=\cos\theta_e\cos\sigma
-\sin\theta_e\mu_e^{(n)}\sin\sigma,
\]

\[
\sin\theta\cos\varphi=\sin\theta_e\cos\sigma
+\cos\theta_e\mu_e^{(n)}\sin\sigma,\qquad
\sin\theta\sin\varphi=s_e^{(n)}\sin\sigma.
\]

The outward root a sin(theta)=Rs determines reception. Normal energy
redshifts as E_o=E_e a_e/a_o, followed by the same held-observer boost.
The independent source-area distance is

\[
D_\perp=D_\parallel=\frac{a_o\sin\sigma_o}{E_e}.
\]

A separate unaveraged DOP853 parent integration supplies the clocks,
scale factor and conformal time. It checks the streamed homogeneous rays
and curvature equations without treating the control as a flat straight line.

## Reuse and numerical boundaries

The field stepper and its source files remain unchanged. It advances a
prefix and then finite blocks. Auxiliary optical states do not feed back
into the scalar; final-state hashes authenticate the complete source replay
against the previous fine, medium, coarse and outer-domain records.
The moving outer boundary remains in the protected parent annulus.

Between field states, the scalar's exact frequency rotation constructs
the midpoint before the geometry is evaluated. The optical RK4 step uses
those metric states; emission and reception are solved within the relevant
block. Refinement tests the remaining interpolation and event errors.
Completed per-run checkpoints survive a validation or compute-budget
failure, so a diagnosis need not discard finished scientific work.

The [compact result](../results/summaries/nsc-atlas-core-optics-v1.json)
records the declared probes, angular caustic brackets, independent optical
comparisons, hashes and finite numerical error envelopes. Those envelopes
are not rigorous interval bounds or uncertainties in the physical preparation.

## Reproduce and continue

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_core_optics.py --output results/summaries/nsc-atlas-core-optics-v1.json --figure
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_core_optics.py
```

The calculation adds a lensing test to the same evolving scalar whose
concentration and radial clock signal are already retained. Its initial
cloud, coupling choice and scale remain explicit inputs. Population
selection, nonspherical persistence, physical calibration and a comparison
with astronomical lensing data remain further connections. A useful next
observer problem is inversion of this map at fixed source and receiving
worldlines, rather than silently identifying different landing points.
