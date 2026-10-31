# Gravitational completion of the inherited scalar kicks

The [canonical response](nsc-atlas-boson-transport.md) held the metric
fixed. The same source can now be completed to gravitational initial
data. This supplies the missing constraint connection before a full
coupled response is claimed. The action, mean scalar profile, vacuum and
observer are unchanged.

## Exact momentum completion

On the existing round-reference slice, write

$$
h=\psi^4\bar h,\qquad K_{ij}=\lambda\psi^{-2}\bar A_{ij},\qquad
\Pi_\chi=\lambda f\psi^{-6},\qquad \operatorname{tr}K=0.
$$

The source f is the original Y1 or Y2 shape normalized in the physical
reference volume. Scalar displacement is unchanged. Both conformal
weights cancel from the momentum equation:

$$
\bar D_j\bar A^{ij}=-\kappa f\bar D^i\chi=:j^i.
$$

The six round-sphere Killing and four gradient-conformal-Killing charges
vanish exactly for both f profiles against chi=aX1+bX2X4. Thus the source
is orthogonal to the York operator's kernel. This is an exact polynomial
moment calculation, not a fitted removal of source components.

For Q=div(L)/r^2, these inherited sources obey
(Q+5)(Q+16)(Q+20/3)(Q+28)j=0. Therefore a closed vector potential is

$$
W=-\frac{3Q^3+167Q^2+2984Q+20080}{44800r^2}\,j,\qquad \bar A=LW.
$$

The finite harmonic inverse and the ambient polynomial inverse are
independent checks. The freely specified conformal TT tensor is zero.
This minimizes the round-reference longitudinal norm; it does not select
a physical radiation state or minimize every physical weighted norm.

## Hamiltonian and GH data

The remaining conformal equation is

$$
-8\bar\Delta\psi+\bar R\psi-\kappa|\bar D\chi|^2\psi-2\kappa W\psi^5
-\lambda^2\left(|\bar A|^2+\kappa f^2\right)\psi^{-7}=0.
$$

It is continued numerically from the existing positive solution. The
additional source is even in lambda. That symmetry alone is not an
existence, uniqueness or differentiability theorem for a finite branch.
At linear order delta h=0 satisfies the Hamiltonian constraint, since
the added normal kinetic and shear terms are quadratic.

[Choquet-Bruhat--Isenberg--Pollack](https://arxiv.org/html/gr-qc/0610045v2),
section2 and4.1, supplies the conformal weights and Fredholm framework.
The potential coefficient above is derived directly from this action's
ADM equation; it is not transferred by notation from another convention.
The positive-potential case does not supply a blanket existence theorem.

For lapse1, shift0 and F=psi^4, GH-compatible initial derivatives are
v_ij=-2K_ij, v_00=2trK=0, v_0i=-E_iF/(2F). Separate finite projections
of these fields need independent ADM and GH diagnostics. Exact conformal
identities do not erase projection errors.

## A full gravitational local relation

With initially constant s, normal(s)=0 and trK=0, the scalar equation gives

$$
\partial_\tau^3\ln g_{\rm KK}=-\kappa\mu_*^2\chi\Pi_\chi.
$$

For K_s=d(d+2)/(2kappa), mu^2=omega0^2 exp(-d s), and
ln g_KK=constant-(d+2)s/2, the internal dimension cancels exactly.
This includes gravitational response: initial expansion is zero, so its
term multiplies normal(s)=0 at the relevant order. A tracefree TT addition
does not alter this scalar identity.

Let sigma be the shear of the initially normal timelike congruence, with
norm sqrt(sigma_ij sigma^ij), without the alternative factor1/2. Then

$$
\frac{H_*^2\|\sigma\|}{|\partial_\tau^3\ln g_{\rm KK}|}
=\frac{F_f(z)}{z},\qquad z=\mu_*^2/H_*^2.
$$

H_star is the vacuum reference curvature scale, not the initial expansion
rate, which is zero. The factors lambda, psi^-6, epsilon and the overall
kick normalization cancel. The matched two-harmonic shape and conformal
TT=0 preparation remain specified. A freely added TT tensor changes the
shear and this ratio.

At the fixed z=8329.074464994 the two ratios are6.92116807623616e-5 and
9.46983923630246e-6. They require a nonzero signal; lambda0 is a0/0
limit, not a measured ratio. They do not identify a Standard Model
coupling or determine an SI ruler, noise state or present cosmology.

## Finite continuation and budget

The next numerical comparison is sourceY1 at lambda=+/-1/32, on the
existing epsilon=.25, degree16 route, to proper endpoint4.364851930252.
The two signs have the same initial scalar energy and intrinsic geometry,
with opposite momentum and extrinsic shear. No parameter is fitted.
Both fields, the metric and the initially normal observer then evolve.

Initial solves have120-second wall budgets and independent denser
quadrature. The two evolutions each have1800 seconds and4 threads,
running together on the8 physical cores. Tolerances are3e-9/3e-11,
maximum step.025. This is one fixed physical grid, not another grid scan.
Failure risks are a positive-branch loss, represented constraint error,
and mistaking finite differences for a certified linear derivative.
Each is recorded before further refinement. A later full-source
measurement must reproduce the same initial shear/coupling relation.
