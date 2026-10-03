# Constraint-compatible continuation of a static collar

The [local collar calculation](nsc-discovery-parent-traction.md) supplies
source-responsive geometry and the gradient charge required at its ends.
This note derives a criterion for extending that collar into a smooth
dynamical parent on one periodic Cauchy slice. It does not yet construct
that parent. It uses the same continuum leading action as the
[leading owner](../src/recursive_horizons/nsc_discovery_leading_einstein.py),
with $g=8\pi A>0$, $\mathfrak m=2\pi C_F\mathrm{flux}^2$ and zero actual
source current. The Dirac density $\rho$ remains the density of the extended
field; it cannot be independently manufactured to make the constraints pass.

Let $C_s$ denote the lapse constraint with both geometric momenta zero:

$$
C_s=-\frac{3g r_x^2}{Q}-gQ(r^2-\mathfrak m/g)
       +2g\,\partial_x(rr_x/Q)+\rho.
$$

For $P=p_Q$, the full constraints are

$$
C=C_s-\frac{Pp_r}{2gr}+\frac{3QP^2}{4gr^2},\qquad
D=p_r r_x-QP_x.
$$

Eliminating $p_r$ where $P$ is nonzero gives

$$
(P^2)_x-\frac{3r_x}{r}P^2=\frac{4grr_x}{Q}C_s,
\qquad
\frac{P^2(x)}{r^3(x)}=J(x)
=4g\int_{x_0}^{x}\frac{r_xC_s}{r^2Q}\,ds.
$$

An anchor $x_0$ in the static collar fixes the integration constant to zero.
A periodic continuation requires $J$ to return to zero and to remain
nonnegative everywhere. In addition, $r^3J$ must admit a smooth **signed**
square root $P$. A plain nonnegative square root can introduce a cusp.

Regularity across momentum zeros is explicit: there must be a smooth $Z$
with $C_s=PZ$ and

$$
P_x=\frac{3r_x}{2r}P+\frac{2grr_x}{Q}Z,
\qquad p_r=\frac{3QP}{2r}+2grZ.
$$

These expressions give $C=D=0$, including turning points by continuity.
The static collar requires $P=Z=0$. On any other open set with $P=0$, the
remaining shift condition is $r_xZ=0$.

The global condition has a useful independent form:

$$
\frac{r_xC_s}{r^2Q}
=\partial_x\!\left(\frac{g r_x^2}{rQ^2}-gr-\frac{\mathfrak m}{r}\right)
 +\frac{\rho r_x}{r^2Q},\qquad
\oint\frac{\rho r_x}{r^2Q}\,dx=0.
$$

Even geometry and density enforce this cyclic integral, but do not enforce
$J\ge0$. For example, a source-free electrovac continuation from a static
BR core with $\mathfrak m=g$ has, at a radius extremum $r_t\ne1$,
$J=4g[-gr_t-g/r_t+2g]<0$. An arbitrary even radius bump therefore cannot
be attached just by choosing real momenta.

The next construction needs a common extended Dirac field and admissible
state, positive geometry, and a smooth exterior satisfying the criterion.
That supplies initial parent data, rather than a maintained wall or pump.
These are smooth continuum identities. A finite SBP realization must also
check its own sampled constraints and admissibility; continuum product
rules cannot be substituted into that check. No new evolution, boundary
controller or matched parent is claimed here.
