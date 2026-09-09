# Finite curved Dirac interaction and a conditional holonomy saddle

The existing compact domain now has an evaluated finite quantum interaction
on the smooth four-dimensional cell. Its effective metric source can have
negative radial null stress without a scalar source. The sign depends on
the global geometry and state. In the ultrastatic charged candidate, a
dynamical flat axial Wilson line has a one-loop minimum at effective
antiperiodic phase. That state-selection result requires a closed axial
return path; it does not turn the development circle into a derived part of
the recursive spacetime.

## Import the interval determinant, preserving its domain

[Flachi, Moss and Toms, equations6–14](https://arxiv.org/html/hep-th/0103138v2)
give the massless conformal interval determinant, its adjoint domain and
finite interaction. Their conformal separation beta is ell=2 L_star here.
Two project Dirac copies give one four-dimensional Dirac zero level and
two Dirac fields at every n>=1 level, with m_n=n pi/ell. The full rank is
eight; the zero level has half the rank of a massive level.

For the **untwisted transverse interval**, the finite interaction is

\[
\Gamma_{\rm int}[g_4,\ell]
=-\frac12\operatorname{Tr}_{4,\;{\rm two\ copies}}
 \log\left(1-e^{-2\ell\sqrt{D_4^2}}\right).
\]

This follows by applying the interval determinant to the eigenvalues of
the same covariant D4, including the half-rank zero level. It is the
interaction between the endpoints, normalized to vanish at infinite
separation for each nonzero D4 eigenvalue. The decompactified bulk
determinant and individual-boundary terms are separate. In curved g4 the
bulk term is generally nonlocal: formally it includes
-ell Tr|D4|/2 and needs its ultraviolet prescription. It cannot be discarded
using the scaleless flat-space subtraction in the published flat example.
Zero modes need their own determinant prescription.

The reference four-geometry is Euclidean
N² d tau²+q² dx²+r² d Omega², with the already implemented paired-angular
Dirac operator. The state is the zero-temperature free vacuum with no
chemical potential and continuous Euclidean time frequency.
The static continuation uses the inherited (+---) source
convention. The fifth-dimensional warp is treated by conformal reduction:
the boundary cocycle remains outside this interaction. The formula is not
the complete finite-cutoff determinant on the warped geometry.

On flat four-space the imported result becomes

\[
V_{\rm int}=\frac{3(8)\zeta_R(5)}{128\pi^2\ell^4}>0.
\]

The twisted transverse interval instead gives -15/16 of this value and
has no n=0 mode. It is a different physical domain. This transverse twist
must also be distinguished from the axial-circle phase below.

## Differentiate the actual curved operator

For an eigenvalue lambda of D4, set
f(lambda)=-log(1-exp(-2 ell |lambda|))/2. The calculation uses

\[
\delta\Gamma_{\rm int}
=\operatorname{Tr} f'(D_4)\delta D_4,\qquad
f'(\lambda)=-\frac{\ell\,{\rm sgn}(\lambda)}
 {e^{2\ell|\lambda|}-1}.
\]

It retains independent N, q and r variations and the ell derivative.
Frequency and angular integration reuse the existing covariant operator.
The resulting rho, radial pressure and sphere pressure are effective
four-dimensional sources, integrated over the transverse interval.
They are not pointwise five-dimensional stress or a surface tension
localized at one endpoint. In particular, they are not yet the unrestricted
variation with respect to one transmitting interface's induced metric.

For the recorded smooth unit-neck cell, axial length4, ell2 and two copies,
the ultrastatic effective neck null source is approximately

\[
(\rho+p_x)_{\rm int}=-8.9162\,10^{-4}.
\]

In the separate nonconstant-lapse/radial-metric control it is approximately
-7.8386e-4. These are development units fixed by the recorded geometry,
not an observed density or a fitted source. The interaction supplies negative
null stress in this finite cell. The decompactified bulk source, boundary
cocycle and absolute local completion remain outside this result.

The numerical record compares matrix differentiation with metric and
interval finite differences, and the frequency quadrature with an
independent Bessel-series integral on the ultrastatic cell. Radial, angular
and frequency refinements are separate. The local radial conservation
residual decreases with radial resolution. A frequency-tail bound applies
only to the energy at fixed radial/angular truncation.

## Keep the same induced terms in one energy account

The first two finite heat moments are
Q2[f]=3 zeta(5)/(8 ell4) and Q1[f]=zeta(3)/(4 ell2).
Together with the existing four-dimensional Dirac heat coefficients, these
give the volume coefficient above and the positive induced Einstein
coefficient

\[
A_{\rm int}=\frac{8\zeta_R(3)}{768\pi^2\ell^2}.
\]

In the repository's Lorentzian-curvature energy convention this multiplies
the existing R_L functional. The record evaluates both terms and their
metric sources, then records the remaining full interaction minus those
terms. No coefficient is fitted or set to zero. This is exact accounting;
it does not justify a short-curvature expansion when ell²|R| is not small.
The next heat moment is infrared singular, so the remainder is kept as
the full spectral function rather than assigned a constant local R² weight.

An ordinary radius minimum also needs joint metric variation.
[Ponton and Poppitz, sections2.1–2.2](https://arxiv.org/html/hep-ph/0105021v1)
explain the Einstein-frame radius and the distinct local counterterms.
Their fitted bulk/brane coefficients are not imported here.
For any reduced action A(ell)R-V(ell), a constant-size metric saddle
requires V'-2A'V/A=0, equivalently stationarity of V/A². A minimum of V
at fixed four-metric does not establish this condition.

## The sign is not determined by the local neck alone

The axial circle is the existing finite spatial development domain.
Enlarging it also changes its global smooth geometry; the table is a
domain-family experiment, not a proved infinite-domain error bound.
The neck radius, transverse ell2 and ultrastatic lapse are held fixed.

| Axial length | Effective neck rho+p_x |
|---|---|
| 4 | -8.9162e-4 |
| 6 | -2.5074e-4 |
| 8 | +1.6402e-4 |
| 12 | +4.9198e-4 |

Periodic axial data on the length4 general-metric control give a positive
source, approximately +1.7211e-3. Thus the negative result cannot be
exported to the isolated throat or obtained by silently choosing a spin
structure with the preferred sign.

## Let the charged operator select its flat holonomy where a loop exists

Use the already proposed common unit-charge U(1) field. A flat connection
along the axial circle has F=0 but can have a physical Wilson line. Define

\[
\alpha=\eta_{\rm spin}-\frac{1}{2\pi}\oint A_x\,dx
\quad({\rm mod}\;1).
\]

Local gauge-invariant bulk/boundary counterterms cannot distinguish these
flat connections. The following result concerns the real free-fermion
one-loop potential when this holonomy is allowed to vary. It is not a
derivation of the full quantum gauge measure or determinant phase.
This is a spatial phase; no temperature is selected by the AP minimum.

For N=1, changing to proper axial distance puts the radial square into the
two real scalar operators

\[
L_\pm=-\partial_s^2+W^2\pm W',\qquad W=\kappa/r>0.
\]

For each compact mass and Euclidean frequency, add
t=m_n²+omega². These factorizable operators have positive periodic lowest
eigenvalue: their possible zero solutions exp(plus/minus integral W ds)
are not periodic. At spectral parameter -t the transfer matrix has
determinant1 and discriminant Delta(t)>2.

Reuse the general boundary determinant formula discussed in
[Dunne, equations18–20 and the periodic examples](https://arxiv.org/html/0711.1178v1) and
[Kirsten–McKane, equation49](https://arxiv.org/html/math-ph/0403050v1).
After removing the unit-modulus boundary phase and fixing the common
high-frequency normalization, the positive determinant depends on alpha as
Delta(t)-2 cos(2 pi alpha). Each fermionic contribution is therefore

\[
-c\log[\Delta(t)-2\cos(2\pi\alpha)],\qquad c>0.
\]

It decreases strictly for 0<alpha<1/2 and is symmetric and periodic.
Its global minimum is alpha=1/2 modulo1, with positive curvature there.
This is an application of established determinant machinery to the actual
positive radial partners; the Wilson-line mechanism itself is prior art.
No new search over heat-kernel cutoffs proves this minimum.

An independent numerical check uses the already derived ultrastatic
frequency integral for the **full** free KK determinant, not just
Gamma_int. The half-rank n=0 level is retained. For axial length4 and
ell2, the cutoff-removed energy difference between alpha1/2 and alpha0 is
approximately -0.16995535. The second derivative in alpha is approximately
3.19343937>0. The current vanishes at the saddle; integer phase shifts
give the same energy within the recorded resolution. Radial, compact,
angular and proper-time-cutoff controls are recorded separately.

Thus the axial effective AP phase need not be inserted by hand in this
charged ultrastatic cell. It may be selected by its fermionic quantum
potential. The bare spin structure and gauge holonomy remain individually
distinct. Other charged sectors, interactions, a prescribed external
holonomy, or a different global domain can change the physical problem.

## What this completes and what comes next

This computes a finite curved quantum source and a conditional gauge
holonomy saddle. It does not solve the metric/scale equations, the local
vacuum completion, the boundary cocycle, a transmitting throat, collapse,
or cosmological observations. A static equilibrium source is not a
continuous injection rate Q.

The next geometric match must establish the actual fermion return path
and apply the computed source together with the other terms of the same
functional. An unwrapped chain has no axial Wilson loop merely because
the development cell has one. The sign reversal above prevents carrying
the finite-cell negative stress into that domain without a new calculation.
If the physical geometry supplies the loop, reuse the holonomy result;
do not rederive generic Casimir, KK or boundary-determinant theory.

The [record](../results/development/compact-casimir.json) authenticates the
prior operator and boundary records. Reproduction runs only this new
interaction/source calculation:

    python scripts/check_nsc_compact_casimir.py --check

All fields are compared, with exact structure/hashes and float tolerances
3e-9 absolute and 3e-8 relative, with no exceptions. Convergence tolerances
and physical scope are distinct from these reproduction tolerances.
