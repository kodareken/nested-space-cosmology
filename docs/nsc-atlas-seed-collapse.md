# Nonlinear continuation of the compact inherited-parent seed

## Decision and admitted calculation

The [formation window](nsc-atlas-formation-window.md) supplies a compact
5% density seed, zero current and a declared earlier parent at A=.01.
Its linear extrapolation alone cannot establish formation. Compare another
linear estimate, a carrier-resolved full-GR run at omega=1e9, and nonlinear
continuation of the existing verified envelope. The third resolves the
missing finite-amplitude dynamics without supplying the stationary core.

Complete: the exact compact seed and its gradient correction satisfy the
local initial constraint. The 5% envelope is checked through A=16, and the
separate 10% preparation through its first compression/rebound at A=1.2.
Fixed-mass contraction, conservation and refinement distinguish the two.
Next: use the formed field for a finite observer calculation. No stationary
profile, core mass or observed collapse endpoint was fitted.

Retain omega=1e9, Lambda=9/32, D=3/8, C=.0227765023676067, comoving
length L=.0004815623798469198 and initial peak seed amplitude .05.
An outer-annulus correction is solved from exact initial geometric mass
compensation; its value and the envelope/GR difference remain explicit.
The earlier homogeneous parent is still a conditional local continuation,
not a completed earlier shared-wall atlas or inherited seed probability.

Reuse the existing symmetric radial Poisson and sine kinetic operators.
Logarithmic scale supplies time, and integrated coefficients retain the
strong early expansion. Three resolutions, independent time/domain,
no-gravity and small-amplitude controls fit a1200 CPU-second budget in
six workers, beginning with one bounded coarse endpoint run. Store full
fields only under `.build`; public results contain summaries and hashes.

Likely failures: insufficient nonlinear growth, unresolved central
concentration, or artificial outer influence. Conserved mass transport,
expansion work, refinement, a larger domain and the prior linear limit
distinguish them. A regular finite envelope concentration is not by itself
an asymptotically stable object, a full-GR collapse or a population model.

The first coarse A=4 trial exposed a boundary inconsistency before any
nonlinear interpretation: the infinity-referenced potential drove a phase
at the finite boundary whose field phase was held fixed. Its mass and work
checks failed and its large central value is rejected. The corrected
potential uses the same boundary clock as the fixed field, subtracting
V(L); this preserves bulk spatial forces and the symmetric Green kernel.
Conservation then passes. At A=4 the central contrast is nonlinear but
the sampled mass shells have not turned around, so the declared finite
endpoint is extended to A=16 before the refinement campaign. Amplitude,
state, action and support remain fixed. The extra future distinguishes
concentration from contrast caused only by the falling ambient density.

At A=16 the corrected coarse 5% case still has expanding fixed-mass shells,
despite a large contrast. Retain that non-formation result. Admit one
separately declared 10% amplitude comparison with the same length, phase,
epoch, action and constraint prescription. It tests a different preparation,
not a parameter fit or a continuation of the identical 5% history. Begin
with its A=4 endpoint and require resolved turnaround and contraction
before calling it finite core formation. The total budget remains1200
CPU seconds; no amplitude scan or inferred observational threshold is added.

The 10% coarse A=4 trial passes through a rapid first compression whose
binding phase outgrows a fixed log-time step. Its late work error is
unacceptable and that endpoint is not retained as evidence. The formation
test is defined on A=.01 to1.2, covering turnaround, first compression and
rebound. Potential-phase-limited steps and work integration at every step
resolve that event. Longer longevity is a different calculation, not
implied by this finite formation claim.

## The compact source is an exact real-space construction

Use r=R/(AL), g(r)=(1-r^2)^4_+ and

\[
S(r)=\frac{-\Delta(g*g)(r)}{\int|\nabla g|^2d^3x}.
\]

The preceding spectral preparation is evaluated here as an exact factored
polynomial on 0<=r<=2, zero outside. Angular integration first gives

\[
(g*g)(r)=\frac{\pi}{5r}\int s(1-s^2)^4
\left[(1-|r-s|^2)^5-(1-\min(1,r+s)^2)^5\right]ds,
\]

with s restricted to 0<=s<=1 and |r-s|<=1. Its two radial pieces simplify
to the same polynomial. Differentiating gives the implemented S; the
normalization is262144 pi/255255. Independent symbolic differentiation
and the previous Bessel integral agree. The radial mass primitive obeys
N_S(r)=integral_0^r s^2 S(s)ds>0 for0<r<2, with N_S(2)=0.
The density changes sign at r=.732791287268. Thus the initial seed
contains a positive central contrast and a compensating outer deficit.

At finite amplitude epsilon, use the real envelope

\[
f_i(r)=\sqrt{1+\epsilon S(r)-\lambda W(r)}.
\]

W is the existing nonnegative C4 bump supported on1<r<2. The coefficient
lambda is not a new coupling or empirical fit. It is fixed by the initial
spherical constraint below. The square root is essential: a linear field
lift 1+epsilon S/2 would introduce an additional positive quadratic density
and fail the intended mass preparation. Both epsilon=.05 and .10 are
declared state inputs. L, phase, epoch and the action are identical.

## Exact local initial gravity and its small envelope difference

Take the real scalar at a turning point,
chi_i=sqrt(2 Eb) f_i/(omega A_i^(3/2)), with Pi_i=0 and the homogeneous
slice expansion H_i. The original potential remains U=Uv+omega^2 chi^2/2.
The vacuum and curvature are unchanged. Define

\[
M_* =\frac{D L^3}{2},\quad n=M/M_*,\quad
q=\frac{1}{\omega^2A_i^2L^2},\quad
\gamma=\frac{D L^2}{A_i},\quad c=C L^2.
\]

The full initial mass constraint becomes the affine equation

\[
n_r=r^2\left[\epsilon S-\lambda W+qB(f_i')^2\right],
\qquad B=1-cr^2-\gamma n/r,\qquad n(0)=0.
\]

Choose lambda so n(2)=0. At lambda=0 the omitted gradient mass is positive;
at zero gradient order the derivative of the endpoint with respect to
lambda is -integral r^2 W dr<0. The implicit-function theorem therefore
gives a nearby compensating correction for small q. The retained positive
B margin and positive squared field verify this branch for both finite
amplitudes. This correction changes the earlier linear seed only beyond
its retained order; it is explicitly part of the new finite-amplitude
preparation rather than an unreported change to the old calculation.

For epsilon=.05 the solved correction is9.75930e-12; for .10 it is
3.85046e-11. The original uncorrected .05 gradient mass is about5.18e-9
of the initial positive excess, agreeing with Grok's independent estimate.
Here it is compensated because the declared initial datum requires the
original exterior geometric mass. Retaining the tiny exterior monopole
would be a different admissible preparation, not the one tested here.

An independent proper-distance Hamilton solve evolves r and h=n/r using
dr/dxi=sqrt(B). Its endpoint agrees with n(2)=0. The field and its
derivative match the homogeneous parent at r=2; B then has the unchanged
local background value. This is a local Einstein-scalar Cauchy constraint
construction. It does not create an earlier evolution of the distant
material walls.

The envelope omits the positive gradient term, so its initial monopole
is -lambda integral r^2 W dr. That small nonzero value is retained in the
evolution's mass balance. Matching a real envelope to exact Pi_i=0 also
has the already identified O(H_i/omega) carrier-quadrature difference.
For this epoch H_i/omega=3.53e-7. Initial compatibility and approximation
error are not the same assertion.

## A boundary phase consistent with the same gravitational clock

Normalize the total envelope by its homogeneous amplitude, f=psi/b, and
evolve u=r(f-1). Write N=ln A. With the inherited H(A),

\[
i f_N=-\frac{a(N)}{2}\Delta_r f+b(N)V f,\qquad
\Delta_r V=\frac{|f|^2-1}{2},
\]
\[
a=\frac{1}{\omega L^2 A^2 H},\qquad
b=\frac{\omega D L^2}{A H}.
\]

The finite boundary fixes f=1. Its clock must also have V=0 there. If
V_infinity is the free-space spherical kernel and Lbox the outer radius,
use

\[
V=V_{\infty}-V_{\infty}(L_{\rm box})
=V_{\infty}+\frac{1}{2L_{\rm box}}
\int_0^{L_{\rm box}}r^2(|f|^2-1)dr.
\]

This changes no bulk spatial force. It is the time-dependent phase choice
associated with the boundary clock, together with the stated boundary
preparation. It is not a vacuum term. The added constant is a symmetric
rank-one change of the density Green kernel, so its variational energy
and conservation identities remain available. A direct radial Laplacian
eigenfunction test verifies the kernel against V=-delta/(2k^2) with
delta(Lbox)=0. The incompatible fixed-field/free-space-phase trial is
rejected rather than used as a formation signal.

## Finite evolution and independent balances

Potential half-kicks rotate the full f, preserving density at that step.
The existing sine kinetic map evolves u. Two-node time quadrature integrates
a and b over each log-scale step. A potential-phase cap resolves the rapid
compression; the time-control run independently tightens it. The boundary
has no absorber. Spatial and outer-domain controls remain necessary.

Let T0=one half integral r^2 |f_r|^2 dr and W0=one half integral
r^2 V(|f|^2-1)dr. The physical-time envelope Hamiltonian is

\[
E=T+W,\qquad T=\frac{T_0}{\omega L^2 A^2},\qquad
W=\frac{\omega D L^2}{A}W_0.
\]

Its exact explicit expansion work is

\[
\frac{dE}{dN}=-(2T+W).
\]

The work integral is evaluated every step. For a fixed comoving cut,

\[
\frac{d}{dN}\int_0^{r_c}r^2(|f|^2-1)dr
=-a r_c^2\operatorname{Im}(f^*f_r)|_{r_c}.
\]

A fixed physical cut additionally has -r_c^3(|f|^2-1), because its comoving
radius decreases. The finite outer boundary current is retained; a change
in central mass is not mistaken for total mass creation. These identities
are checked independently of the time-step updates.

To distinguish actual contraction from increasing contrast, track shells
of fixed **oscillatory-source mass**,

\[
M_{\rm shell}=M_*\int_0^{r_m(N)}r^2|f|^2dr,
\qquad R_m(N)=A L r_m(N).
\]

This is an exactly conserved mass label in the envelope continuity law.
It excludes the vacuum's volume contribution and higher-order gradient
corrections to full GR mass. The labels start at comoving radii
.05,.1,.25,.5 and1. An exact-volume cell integral reconstructs them;
uniform density recovers the cubic volume law even inside the first cell.
Their physical turnaround and contraction are measured separately from
the central density contrast and its physical half-height radius.

## Resolved finite formation and a contrasting expanding history

The [compact record](../results/summaries/nsc-atlas-seed-collapse-v1.json)
contains eleven current-source runs, using six workers and about276 CPU
seconds. The earlier rejected pilots and interrupted first campaign were
also below the total1200-second work allocation. Raw histories and fields
remain local, authenticated by hashes and reproducible from the producer.

The 10% case has three grids of4095,8191 and16383 interior points, log-step
limits .002,.001 and.0005 and potential-phase caps .08,.04 and.02. The
medium grid has an independently smaller time step and a larger domain
with the same spatial spacing. The 5% comparison has its own three-grid
continuation to A=16. Separate no-gravity and two small-amplitude controls
complete the campaign. The small-amplitude derivative tends toward the
previous packet's linear solution: halving epsilon from .0005 to .00025
reduces the central gain discrepancy from .20251 to .10088.

At A=1.2 the 10% case gives:

| Quantity | Fine value and numerical comparison range |
| --- | ---: |
| Central density contrast above ambient | 59971 +/-1514 |
| Physical half-height radius | (4.22172 +/- .02216)e-6 |
| Radius of the shell initially at r=.25 | (4.38056 +/- .03788)e-6 |
| That shell's final radius / largest sampled radius | .109898 +/- .000950 |

Ranges are twice the maximum fine/medium, time-control/medium or
outer-control/medium difference. They describe the numerical envelope
calculation, not total approximation or astrophysical uncertainty.
The selected mass shell turns around near A=.5852 (the saved log-scale
spacing is about .008) and contracts by about89% from its maximum radius.
The physical half-height radius also contracts and undergoes a first
rebound. This is finite regular concentration from a weak initial density
seed without an installed stationary core. It does not establish an
asymptotic soliton, nonspherical stability or exact full-GR evolution.

The no-self-gravity control keeps all sampled mass shells expanding.
The 5% case likewise has no sampled shell turnaround through A=16,
although its final central contrast is about2257. A large density contrast
can arise while the ambient density declines and the physical region
continues expanding. This is why it is not used alone as the formation
criterion. These two amplitudes do not determine a universal threshold,
nor do they establish the 5% history's fate after the declared interval.

The fixed-comoving mass-flux discrepancy, relative to initial positive
excess, contracts from8.88e-5 to2.11e-5 to5.48e-6. The expansion-work
discrepancy contracts from.007108 to.001689 to.0004185. The fine global
mass/current discrepancy is6.69e-10. The outer-domain change in final
central contrast is9.25e-6, compared with a medium value of60727; the
remaining central uncertainty is spatial resolution rather than the outer
phase instability of the rejected trial. The maximum retained potential
depth relative to the boundary clock is1.45e-7 in units c squared.
This indicates weak gravitational potentials, not a rigorous GR error bar.

## A conditional physical mass and radius, without a galaxy fit

Using the same previously supplied Lambda=1e-52 m^-2 ruler, the selected
fixed-source shell contains4.251e9 nominal solar masses. Its largest sampled
radius is68.51 kpc; its final radius is7.529 +/- .066 kpc. The final density
half-height radius is7.256 +/- .039 kpc and the scalar source mass inside
that radius is(3.919 +/- .069)e9 nominal solar masses. The latter is an
output aperture mass, not the preset total core mass used by older examples.

The dimensional ruler, seed length, epoch and amplitude are supplied.
No galaxy measurement selected them. Formation of this finite concentration
does not fix the dark-matter abundance, the primordial seed distribution or
the physical value of Lambda. It does provide a same-source mass profile
for the next operational clock or lens measurement.

## Attribution, verification and continuation

The general scalar concentration mechanism is established work, including
[Seidel--Suen,gr-qc/9309015](https://arxiv.org/abs/gr-qc/9309015).
[Banik et al.,1504.05968](https://arxiv.org/abs/1504.05968) discusses the
expanding Schrodinger--Poisson/real-scalar dictionary. Their abstracts and
bibliographic records were checked here; our earlier owners retain the
more detailed reduction references. No new general collapse mechanism is
claimed. The contribution is the compact source, exact local initial
constraint correction, phase-consistent finite boundary and converged
nonlinear formation comparison on this inherited parent.

Grok independently recovered the square-root preparation, expansion-work
and moving-cut identities and the initial positive gradient-mass fraction.
It recommended keeping that tiny monopole; the present calculation instead
solves the explicitly declared initial compensation condition, preserving
the input amplitude and length. This is a stated scientific choice checked
by the independent Hamilton constraint, not an appeal to model agreement.

Tests verify the exact polynomial convolution, its mass primitive, agreement
with the preceding spectral seed, positive initial metric, independent
proper-distance endpoint, the phase-matched Poisson eigenmode, volume-shell
inversion and source identities. The two physics controls and shrinking
balance discrepancies are evaluated in the compact producer record.

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_seed_collapse.py --output results/summaries/nsc-atlas-seed-collapse-v1.json --figure
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_seed_collapse.py
```

The next connection is an observer response from this actual forming field.
The old omega2000 lens and the installed omega1e9 core are different
histories and cannot supply its signal. The finite formation region is
complete within its stated checks; further longevity or shape studies
should be driven by a later observational requirement or a failing test.
