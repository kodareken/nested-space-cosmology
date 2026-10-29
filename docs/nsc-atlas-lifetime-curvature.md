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
