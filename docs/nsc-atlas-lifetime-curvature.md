# Lifetime selection and inherited curvature

## Lifetime coefficients, computed first

The [lifetime record](../results/summaries/nsc-atlas-volume-lifetime-v1.json)
uses the same g=kappa_nat/L_unit(0)^2 for all five candidates.
The [CDL equations](https://doi.org/10.1103/PhysRevD.21.3305) and
[Hawking--Moss saddle](https://doi.org/10.1016/0370-2693%2882%2990946-7)
are reused in the declared Einstein/canonical-volume action. The
[global versus barrier-top criterion](https://arxiv.org/html/gr-qc/0311040)
explains why a top ratio below four alone does not exclude CDL.

| (d,j) | Lambda_star L_ref(0)^2 | B_HM g | Regular compact CDL status |
| --- | ---: | ---: | --- |
| (12,2) | 9.702237009265 | 403.781741036 | No nonconstant solution found in bounded shooting; absence unresolved |
| (6,3) | 9.991573647602 | 369.001286854 | Excluded by the global identity below |
| (5,4) | 10.635176064650 | 326.046593819 | Excluded by the global identity below |
| (4,5) | 10.776179311559 | 314.036386367 | Excluded by the global identity below |
| (3,7) | 10.727932243144 | 305.543559374 | Excluded by the global identity below |

The HM coefficients use 24*pi^2*(1/V_false-1/V_top). All five have one
negative homogeneous scalar harmonic: their barrier ratios are 3.76250,
3.58398, 3.49617, 3.40670 and 3.27186. No thin-wall approximation is
introduced for these broad barriers. No finite CDL action is reported
for an absent or unresolved saddle.

### A global exclusion, including thick profiles

For canonical phi=sqrt(K)*s, the regular O(4) equations are
phi''+3(a'/a)phi'=V_phi and a''=-a*(phi'^2+V)/3. Let y=phi'/a and
h=a'/a. Differentiation and integration between regular poles give

$$
0=\int a^5\left[y'^2+
\left(V_{\phi\phi}+\frac43(\phi'^2+V)\right)y^2\right]d\xi.
$$

Complete the square with Z=h/10+13*phi'/20. The same identity becomes

$$
0=\int a^5\left[(y'+Zy)^2+
\left(G+\frac{39}{100}(h+3\phi'/2)^2\right)y^2\right]d\xi,
$$

$$
G=V_{\phi\phi}+\frac{13}{10}V+\frac{13}{20}V_\phi.
$$

The extra boundary term a^5*Z*y^2 vanishes at each regular pole.
Consequently G>0 for every field value rules out every nonconstant
regular compact O(4) solution, including oscillating and thick profiles.
The constants 1/10 and 13/20 choose a mathematical square completion;
they add no parameter to the physical action.

With x=exp(-2s), G=exp(-ds)*(a-b*x+c*x^(d/2)). The coefficients follow
by differentiating the existing potential; a,b,c are positive. The
bracket is convex and its unique minimum is explicit. Its values for
the five members are -.01863743, .02689417, .05554020, .07775262 and
.10379209. Exact rational bounds on the square roots and critical point
certify positivity in the last four cases. The first case fails this
sufficient exclusion test; that is not evidence that a CDL bounce exists.

The extra original scalar cannot evade this restriction in the retained
truncation: its Euclidean equation is Laplacian(chi)=omega(s)^2*chi with
positive omega(s)^2. Integration against chi on a regular compact manifold
forces chi=0. Other internal modes and singular instantons are outside
this calculation.

For (12,2), regular pole series are integrated independently to the first
equator a'=0. Matching the two field values and opposite normal derivatives
would give a CDL solution; the Friedmann constraint then fixes the same
equatorial radius. Eight bounded starting pairs approach the HM corner.
The search found no admitted nonconstant root and checked the independent
equator constraint. Its finite domain and residuals are retained. It is
not an exhaustive root count, so a reference CDL action remains unresolved.
This unresolved step prevents a complete escape-channel lifetime ranking.

### Conditional survival theorem

Under the proposed survival rule, let independent prepared populations
have rates Gamma_i=A_i(g)*exp(-b_i/g), with positive subexponential
prefactors, positive initial populations, no replenishment, and absorbing
escape. Conditional survivors obey

$$
p_i(t)=\frac{n_i(0)e^{-\Gamma_i t}}
{\sum_k n_k(0)e^{-\Gamma_k t}}.
$$

For sufficiently small common g, the largest leading action is the
smallest rate; at late t its conditional survivor weight tends to one.
Therefore the HM-only survival model selects (12,2), with late coefficient
9.702237009265. The gap to the next HM coefficient is 34.7804541813.
This theorem states the ensemble and rate assumptions rather than adding
a fitted physical parameter. Common g alone does not specify prefactors,
initial populations, production or a volume-weighted measure.

The four competing regular scalar CDL channels are excluded. A faster
unresolved channel in (12,2) could still change the ranking; if its action
remained above 369.001286854/g, the reference would still win at leading
exponential order. No unconditional selection or absolute lifetime is
claimed. The finite shooting allowance is 1600 residual calls or 120
process-CPU seconds; no failed pole search is silently repeated.

## Curvature mechanisms, evaluated second

The [curvature audit](../results/summaries/nsc-atlas-curvature-mechanisms-v1.json)
keeps the inherited source and its C=61720739711/2709842745600 fixed.

| Route | Outcome |
| --- | --- |
| Open CDL interior | It fixes open sign, not the observed magnitude; compatible nucleation, post-bubble expansion and exit history are missing. |
| Inherited oscillating matter | Dust domination increases relative curvature; late vacuum acceleration dilutes it only after the needed matter has diluted. |
| Repeated similarity | a -> r*a and H -> H/r preserve Omega_k; a curvature/state transfer law and observer level are still needed. |

[Gott (1982)](https://doi.org/10.1038/295304a0) and
[Bucher--Goldhaber--Turok (1995)](https://arxiv.org/html/hep-ph/9411206)
provide the open-bubble construction and its subsequent inflationary
phase. Immediately after the open-bubble origin, a(t) is approximately t
and Omega_k is approximately +1. Near-flatness requires additional expansion;
its present value also depends on the exit and matter/radiation history.
Their examples do not fix those inputs for the current parent. The HM
sphere is not automatically a localized CDL child with that history.

The exact geometric derivative is

$$
\frac{d\log|\Omega_k|}{dN}=-2-2\frac{\dot H}{H^2}=2q_{\rm dec}.
$$

Dust domination has H^2 proportional to a^(-3), so |Omega_k| grows
proportionally to a. Acceleration can reverse this, but the present
source supplies no early inflation-and-reheating stage that resets
its dimensionless curvature/matter normalization.

There is a quantitative obstruction independent of a preferred observer
epoch. The inherited density is E=E0*exp(-3N-6s), hence

$$
\frac{|\Omega_k|}{\Omega_m}=\frac{3C}{E_0}e^{N+6s}.
$$

Along the expanding averaged solution total energy decreases from 21/32.
Exact rational exponential bounds give V(-1/12)>21/32; since the potential
decreases toward its unique minimum on that side, s>-1/12. For N>=0,

$$
\frac{|\Omega_k|}{\Omega_m}>
\frac{3C}{E_0}e^{-1/2}>0.11051717.
$$

Thus |Omega_k|<=.002 implies Omega_m<.01809674 in this fixed preparation.
Changing only the observer epoch cannot combine the requested near-flatness
with a matter fraction near the fitted .3. This is a restriction on this
source, not a general obstruction to nesting or inflation.

The actual first future crossing Omega_k=-.002 occurs at N=2.52065449,
where Omega_m=.00109975. At the previously declared N=6 endpoint,
Omega_k=-1.8998373e-6 and Omega_m=3.2008896e-8. These are numbers from
the trajectory. Neither endpoint is an inherited identification of our
observation epoch. The infinite-future limit is flat and empty of matter.

## Fixed-flatness diagnostic, declared before the retained fit

Because no tested inheritance mechanism selected a small observed
curvature, a physical refit with a newly derived value is not available.
We instead test the user's causal diagnosis with an explicitly imposed
initial C=0. It is fixed before fitting and is never labelled inherited.
Keep d=12,j=2,A=7,B=15/2,T=25/32,E0=3/8 and s=v=0. The homogeneous
constraint then determines the changed initial H; all subsequent source
and volume evolution is recalculated consistently. This does not preserve
the original shared-wall initial geometry.

V_star=.2209830242 remains the potential minimum in the same source units.
The earlier coefficient 9.7022370 used the original closed reference ruler;
it is not a finite curvature-radius product for the flat control, whose
intrinsic curvature radius is infinite. No SI scale is introduced.

Use the same authenticated SN and DESI covariances, one epoch, one
magnitude intercept and one shared ruler. The closed source is replayed
as a control and must recover its existing joint score. Four independent
model/order jobs compare orders 48/96 and background tolerances, with
180 process-CPU seconds each. The report includes Omega_Lambda, Omega_s,
Omega_m, Omega_k and the posterior of their requested dark-energy sum.
That sum is evaluated at every posterior node, so its covariance is
included; independent component errors are not added in quadrature.
