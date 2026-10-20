# From an occupied sphere spectrum to two sourced rooms and their clocks

The [constant-tension test](nsc-atlas-parent-observables.md#d-does-the-actual-wall-select-its-own-cap-scale)
left a specific source question. This construction replaces the material
law by an explicitly declared surface field and state, derives its pressure,
and follows the resulting geometry into the already checked dust parent.
It supplies a linear-radially-stable, horizonless configuration and two
unequal rooms under the same action. It does not stabilize the former
wormhole orientation or select the cosmological vacuum.

![Two sourced rooms and their state-dependent clocks](figures/nsc-atlas-occupied-wall.png)

## The material candidate and what is being imported

Keep Einstein gravity, a de Sitter core and a Schwarzschild collar, with
the scalar background zero. Change the wall action explicitly to

$$
S_\Sigma=-\sigma_0\int_\Sigma\sqrt{-\gamma}\,d^3y
+\sum_{j=1}^{g}\int_\Sigma\sqrt{-\gamma}\,
\overline\chi_j\,i\hbar\slashed D_\gamma\chi_j\,d^3y,\qquad g=2.
\tag{F1}
$$

The two species are a parity-paired massless worldvolume Dirac sector.
They carry conserved fermion number, not an asserted electric charge.
Use the conformal vacuum as reference, its normal-ordered finite
excitation stress, a specified renormalized sigma_0>0, and no additional
surface-curvature terms. These are the effective theory's stated
renormalization/material choices. The bulk vacuum difference remains
supplied independently.

The original bulk Dirac probe was transparent and had no confining
interaction. It does not become F1 through a reinterpretation.
[Macpherson--Campbell](https://arxiv.org/html/hep-ph/9408387), equations
(1)--(6), gives a scalar/Yukawa localization mechanism.
[Ogure--Yoshida--Arafune](https://arxiv.org/html/hep-ph/0212332), equations
(29)--(36), gives the bound radial state and finite spherical spectrum.
These motivate F1; a finite-thickness four-dimensional realization of
this gravitating candidate has not been constructed. Those papers also
make fragmentation an essential question: radial stability alone does
not establish general stability or formation.

## A state-derived pressure, with exact spherical counting

For one two-component spinor the positive unit-sphere eigenvalues are
n=1,2,..., with degeneracy 2n. The negative branch defines the reference
sea and is not counted again as occupied positive particles. This is the
spherical Dirac spectrum in
[Abrikosov, Sections 2.1--2.2](https://arxiv.org/html/hep-th/0212134).
Fill every angular state through level q in each species. Rotational
invariance of the filled projectors makes the mean stress isotropic.

$$
N_q=gq(q+1),\qquad
\mathcal A_q=\frac{gq(q+1)(2q+1)}{3},\qquad
E_q(R)=\frac{\hbar c\,\mathcal A_q}{R}.
\tag{F2}
$$

The finite sum is exact; no continuum density of states is substituted.
Its surface density and pressure are

$$
\sigma=\sigma_0+\frac{\hbar c\,\mathcal A_q}{4\pi R^3},
\qquad
p=-\sigma_0+\frac{\hbar c\,\mathcal A_q}{8\pi R^3}.
\tag{F3}
$$

Pressure follows from -dE/d(4 pi R^2) at fixed occupation. The gas
has p_g=sigma_g/2, while the bare membrane has p_0=-sigma_0.
Thus d(sigma R^2)/dR+2pR=0 and the gas adds no surface trace,
although it changes the full gravitational source. This is a calculated
material response, not a chosen potential that was made stable.

The transport law is exact within this massless surface theory.
For gamma=-d tau^2+R(tau)^2 dOmega^2, use d eta=d tau/R and
chi=R^-1 chi_cyl. The rescaled Hamiltonian is the fixed unit-sphere
operator. Its evolution is exp(-i D_S2 eta); occupied projectors commute
with it. N_q and A_q are conserved, E_q scales as 1/R, and F3 obeys
the surface Ward identity without an adiabatic approximation.
This transports a specified state along a wall. It does not decide which
state a newly formed child receives.

## Signed junctions distinguish two different geometries

Use c=1 and geometric masses. Write k_0=4 pi G sigma_0,
d=G hbar A_q=ell_P^2 A_q and
s(R)=4 pi G sigma R=k_0 R+d/R^2. Let
f_-=1-H^2R^2 and f_+=1-2M/R. The signed rest/dynamic junction is

$$
\sqrt{\dot R^2+f_-}-\epsilon\sqrt{\dot R^2+f_+}=s(R).
\tag{F4}
$$

The ordinary enclosed region has epsilon=+1. The previous
wormhole-side exterior has epsilon=-1. The latter cannot borrow
the former's stability.

At a static point put a=sqrt(f_-)>0, b=sqrt(f_+)>0 and
F(R,z)=sqrt(z+f_-)-epsilon sqrt(z+f_+)-s(R).
If F=F_R=0, the potential curvature is V''=F_RR/F_z.
On the wormhole branch,

$$
F_z=\frac{1}{2a}+\frac{1}{2b}>0,
$$

$$
F_{RR}=-\frac{H^2}{a^3}
-\frac{2M}{bR^3}-\frac{M^2}{b^3R^4}-\frac{6d}{R^4}<0.
\tag{F5}
$$

For M>0, d>=0, every allowed static point remains unstable.
This sign proof needs no radius scan. A nonnegative inner Kottler
mass adds further nonpositive terms to the corresponding concavity.

On the ordinary branch the exterior term has the opposite sign.
It permits the following fixed reference, chosen before any sky-data
comparison:

$$
a=\frac{24}{25},\quad b=\frac45,\quad
H^2R_0^2=\frac{49}{625},\quad \frac{M_0}{R_0}=\frac9{50}.
\tag{F6}
$$

The angular and temporal Israel conditions independently require

$$
4\pi G R_0\sigma=\frac4{25},\quad
4\pi G R_0p=\frac{11}{150},\quad
k_0R_0=\frac1{225},\quad \frac{d}{R_0^2}=\frac7{45}.
\tag{F7}
$$

Both material components are positive, sigma>|p|, and
dp/dsigma=1/2 for the fixed-tension occupied gas. There is a regular
centre, HR_0=.28<1, and 2M_0/R_0=.36<1. This configuration is
horizonless. It is a different admissible attachment, not an evolved
change of the earlier black-hole bridge into this orientation.

The full potential follows from F4 by retaining its signed factors.
Direct differentiation, with M,H,sigma_0 and q fixed, gives

$$
V(R_0)=V'(R_0)=0,\qquad
R_0^2V''(R_0)=\frac{35159}{7200}>0,\qquad
\frac{\omega_{\rm radial}R_0}{c}
=\frac{\sqrt{35159}}{120}=1.562561109916.
\tag{F8}
$$

This proves linear radial stability within the spherical mean-stress
model. It does not prove nonspherical or fragmentation stability.
Finite-amplitude disturbances carry energy: their exterior mass changes
at second order. No finite oscillation at exactly unchanged minimum
mass is asserted.

## Hold the action fixed and change only the occupied state

Take q_0=1000 and g=2 as the reference state. It has
N_0=2,002,000 and A_0=1,335,334,000. F7 determines the reference ruler,

$$
\frac{R_0}{\ell_P}=\sqrt{\frac{45\mathcal A_0}{7}}
=92651.4435937.
\tag{F9}
$$

This sets a physical scale in Planck units; it is not a cosmological
scale prediction. The dimensionless shape F6 and the initial state
are supplied, and determine which action parameters the benchmark uses.
In particular Lambda_cap R_0^2=147/625=.2352 is the specified cap
curvature, not a vacuum energy generated by the wall occupation.

After this declaration keep G,hbar,sigma_0,H and the parent density fixed.
For any admitted nearby q, define the shell's geometric rest mass
m_s=k_0R^2+d/R. The signed static equation yields

$$
M(R,q)=\frac{H^2R^3}{2}
+m_s\sqrt{1-H^2R^2}-\frac{m_s^2}{2R}.
\tag{F10}
$$

Equilibrium is M_R=0. At the reference,
R_0 M_RR=35159/86400>0. The implicit-function theorem gives a
locally unique equilibrium family under a continuous variation of d;
the allowed integer occupations sample that family.
Along it the exact first-law derivative is

$$
\frac{dM_{\rm eq}}{d\mathcal A}
=\frac{\ell_P^2\,\beta_+}{R}.
\tag{F11}
$$

The gravitational mass increase is the added local spectral energy,
redshifted by the exterior lapse. This connects state counting to the
mass already used in the parent junction; no interior energy is added twice.

| q | N_q | R/R_0 | M/R_0 | omega_radial R/c |
| ---: | ---: | ---: | ---: | ---: |
| 999 | 1998000 | .9993406734 | .1796270290 | 1.563237585 |
| 1000 | 2002000 | 1.0000000000 | .1800000000 | 1.562561110 |
| 1001 | 2006004 | 1.0006588232 | .1803733226 | 1.561883618 |

All couplings are shared. Only the finite occupation differs.
The adjacent-level spacing on the wall is c/R; the last column is
therefore also a dimensionless radial-mode/spectral-spacing ratio.
This is a theoretical frequency relation. A detector coupling and its
selection rules are needed to turn that spacing into an observed line.

## Two unequal interiors in one evolving cosmology

Reuse the [exact dust join](nsc-atlas-parent-observables.md#a-replace-the-isolated-end-by-a-cosmological-parent).
Fix the common parent H_p(0)R_0/c=.075 and a_p(0)=1. Each exterior
mass fixes its comoving vacuole radius:

$$
R_{b,q}(0)=\left(\frac{2M_q}{H_p(0)^2}\right)^{1/3},
\qquad
a_p(\tau)=\left(1+\frac32H_p(0)\tau\right)^{2/3}.
\tag{F12}
$$

Place q=999 and q=1001 at parent comoving locations
(-6,0,0)R_0 and (6,0,0)R_0. Their join radii are
3.997235341585R_0 and 4.002763442849R_0. The initial comoving
clearance is 4.000001215565R_0. Both holes and their separation
scale with a_p, so they stay disjoint for tau>=0.

The exterior dust metric is the same on the complement of the two
balls. Each complete interior supplies its own compatible timelike
matching. Disjoint local replacements therefore define one shared
piecewise Einstein solution, with two unequal, physically sourced
interiors and no fictitious second child. The construction gives initial
states and their ideal spherical evolution, not a formation or
parent-to-child branching process.

## The beacon crosses the entire connection

Unlike the previous right-exterior experiment, this beacon starts at
the regular centre and passes through its material wall, Schwarzschild
collar, moving dust join and the common cosmological region.
The null probe is transparent at the static wall and has negligible
stress. Let t_c be central proper time and T exterior Killing time.
Matching the wall clock fixes dT=(a/b)dt_c.
For an emission synchronized to T_b(0), the core transit adds

$$
\Delta T_{\rm core}=\frac{a}{bH}\operatorname{atanh}(HR).
\tag{F13}
$$

With rstar the Schwarzschild tortoise radius and v_b=sqrt(2M/R_b),
the join reception solves

$$
\int_0^{\tau_b}\frac{d\tau}{1+v_b(\tau)}
=\Delta T_{\rm core}+r_*(R_b(0))-r_*(R).
\tag{F14}
$$

Put an observer 8R_0 in the +y comoving direction from each centre.
These rays avoid the other vacuole. The dust conformal time is
eta_p=(2/H_p(0))[(1+3H_p(0)tau/2)^(1/3)-1]. It gives
eta_r=eta_b+8R_0-R_b(0). The total frequency ratio is

$$
1+z=\frac{a}{b}(1+v_b(\tau_b))
\frac{a_p(\tau_r)}{a_p(\tau_b)}.
\tag{F15}
$$

| q | Join reception / (R_0/c) | Observer reception / (R_0/c) | 1+z |
| ---: | ---: | ---: | ---: |
| 999 | 6.177140238 | 12.61376979 | 1.903001657 |
| 1000 | 6.184208815 | 12.61774703 | 1.903495521 |
| 1001 | 6.191281841 | 12.62172978 | 1.903990397 |

A piecewise null ODE independently crosses all three metric regions.
The derivative of reception proper time with respect to emission
proper time agrees with F15. The same parent density and the same
observer comoving offset are retained across the state comparison.
No measured sky number selected these occupations or couplings.

## What inheritance now determines, and what it leaves open

The normalized sphere operator and its unitary evolution determine
F2--F3 and transport the occupation exactly. This is a concrete
state-to-pressure law. It is stronger than an unspecified
sigma(R), but does not select the initial q.

Copying the same occupation to a radius scaled by a factor c_R at fixed
fundamental constants multiplies k_0R by c_R and divides
ell_P^2 A_q/R^2 by c_R^2. Hence operator inheritance alone does not
make the gravitating configurations exact scale copies. The shell's
physical couplings and state must also be related by a transition law.
Charge conservation would constrain a branching transfer, but its
partition and actual causal creation have not been constructed.

The common dust parent still has zero cosmological Lambda, even though
its interiors have positive supplied vacuum curvature. The construction
does not derive Omega_Lambda, Omega_c, Omega_b or baryon asymmetry.
The next physical connection is state preparation/transfer or a collective
source law; radial stability does not settle the known shape and
fragmentation questions of localized fermion walls.

## Reproduction

The [component](../src/recursive_horizons/nsc_atlas_occupied_wall.py) reuses
the existing cosmological matching evaluator. The
[compact record](../results/summaries/nsc-atlas-occupied-wall-v1.json)
stores three states and scalar observer comparisons; no field grids.
Independent checks include exact sums, symbolic second derivatives,
both signed Israel conditions, a 50-digit neighboring root, the mass
first law and a piecewise null ODE.

    python scripts/run_python.py scripts/check_nsc_atlas_occupied_wall.py
    python scripts/run_python.py scripts/plot_nsc_atlas_occupied_wall.py
    python scripts/run_python.py -m pytest tests/test_nsc_atlas_occupied_wall.py -q

The maximum saved junction residual is 2.1e-15, the neighboring
radius differs from the independent root by 2.9e-15 R_0, and the
clock-derivative difference is 1.2e-11. These are numerical indicators,
not estimates of finite-thickness, quantum-gravity or observational error.
