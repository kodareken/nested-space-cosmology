# Replacing the supplied wall radius and mass: the nucleation test

The target is to determine the wall's radius and seed mass before the
cosmological preparation, then evolve to E=3/8 without resetting the
state. The existing action, field couplings and reference units remain
fixed. This note separates a solved vacuum-wall geometry from its
admission as a solution of that action.

## Radius and mass in the stipulated thin-wall model

[Parke (1983)](https://doi.org/10.1016/0370-2693%2883%2991376-X)
extends [Coleman--De Luccia (1980)](https://doi.org/10.1103/PhysRevD.21.3305)
to arbitrary false and true vacuum cosmological constants. Equations
(7)--(12) of the [public preprint](https://www.slac.stanford.edu/pubs/slacpubs/2750/slac-pub-2908.pdf)
were inspected directly; the result is also checked against unsquared
Israel equations and the existing isolated-bubble calculator.

For constant vacuum densities U_i=7/32, U_o=21/32 and constant tension
sigma=1/5, with kappa=1, the ordinary positive-cap branch gives

$$
R_c^{-2}=\left(\frac{U_o-U_i}{3\sigma}\right)^2
+\frac{\kappa(U_i+U_o)}6+\frac{(\kappa\sigma)^2}{16}
=\frac{39169}{57600}.
$$

Thus the [retained geometry result](../results/summaries/nsc-atlas-parke-bubble-v1.json) is

$$
R_c=\frac{240}{\sqrt{39169}}=1.212662645008\ L_{\rm unit}(0),
\qquad M_i=M_o=0.
$$

No supplied R=1.9 or M=.12 enters this radius calculation. Its ratio to
the supplied radius is .6382434974. The old number 1.9 was in model
length units, not in L_ref(0) units. Using the historical conversion
L_ref(0)/L_unit(0)=6.6260770761 only for comparison gives
R_c/L_ref_old=.1830136642 versus R_old/L_ref_old=.2867458344.
That historical ruler depended on the removed preparation and is not a
newly selected scale. An independent vacuum ruler is
ell_o=sqrt(3/(kappa*U_o)); R_c/ell_o=.5671710179.

At the turning point, beta_i=187/sqrt(39169) and
beta_o=163/sqrt(39169), so beta_i-beta_o=(sigma/2)R_c exactly.
Regular homogeneous vacuum regions have zero Schwarzschild seed offsets.
The wall energy identity verifies that they are compatible:

$$
0=-\frac{(U_o-U_i)R_c^3}{6}
+\frac{\sigma R_c^2\beta_i}{2}
-\frac{\sigma^2R_c^3}{8}.
$$

The three terms are -.130030971151, +.138947380602 and -.008916409450.
Zero seed offset is not zero total vacuum energy: the interior
Misner--Sharp mass is U_i R_c^3/6=.065015485576. Nor does zero seed mass
by itself determine whether a homogeneous scalar is excited.

The auxiliary wall has R(tau_wall)=R_c*cosh(tau_wall/R_c). Its vacuum
open interior has a(t)=sinh(H_i t)/H_i, H_i^2=7/96 and
Omega_k(t)=sech(H_i t)^2, approaching +1 at its open origin.
Wall proper time and interior cosmic time differ; the finite wall
turning radius does not identify a matter-defined FRW reference slice.

## Admission into the unchanged (12,2) source

The [source check](../results/summaries/nsc-atlas-parke-reheating-admission-v1.json)
does not promote the preceding control to a scalar/volume instanton.
The stated densities 7/32 and 21/32 came from an excited scalar snapshot.
They are not stationary vacua of the fixed action. Writing y=omega*phi,
the derivatives at the supplied grandchild-minimum field are

$$
\left.U_y\right|_{q=1}=\frac{\sqrt3}{4},\qquad
\left.U_y\right|_{q=3/2}=\frac{\sqrt3}{2}.
$$

They are nonzero. The relaxed scalar vacuum densities are 1/8 and 9/32,
which would be different thin-wall inputs. The volume source introduces
another incompatibility: its well minimum is V_star=.220983024213,
whereas U_i=.21875 is below it. Its actual stationary extrema do not
supply the requested pair of vacua. The wrapped wall also has
sigma(s)=.2*exp(-6s); its field value and volume-gradient jump must be
obtained from a compatible profile. A constant reference tension is not
the solution of that field problem.

The Parke formula fixes the geometric radius once its vacuum and wall
assumptions are met. It does not determine a missing field turning point
or establish that these particular source values meet those assumptions.
No s_CDL has been inferred from R_c.

## Empty-start evolution and the reference-event gate

For the stipulated empty retained matter field, chi=Pi=0, in the
homogeneous post-wall region at fixed flux and without incoming sources,
the exact classical scalar equations preserve that state for any s(t).
The averaged matter equation makes the same fact explicit:

$$
\dot E+(3H+6\dot s)E=0,\qquad
E=E_i e^{-3(N-N_i)-6(s-s_i)}.
$$

Therefore E_i=0 implies E(t)=0. This conclusion does not substitute the
HM top for a CDL turning point; it holds for any such turning point in
the supplied classical equations. An empty matter field is an explicit
branch assumption, not a deduction from M_offset=0 alone. The latter
would still allow homogeneous field energy fixed by an actual bounce.

Interpreting the volume's oscillation energy itself as matter also needs
the correct budget. For any positive, energy-conserving transfer between
the canonical volume and pressureless matter on an expanding homogeneous
branch,

$$
\dot\rho_{\rm total}=-3H(K\dot s^2+E)\leq0.
$$

A zero-velocity exit that can reach the well with E=3/8 consequently
requires initial potential energy at least

$$
V_{\rm exit}\geq V_*+\frac38=0.595983024213.
$$

In the literal thin-wall vacuum interior, even its entire initial density
U_i=7/32 is below the matter target 3/8. Converting all of it to matter
would still be insufficient on a source-free expanding branch; the target
is 12/7 times that initial density. Keeping the required positive vacuum
makes the deficit stronger.

Within the previously computed minimum-to-barrier interval the available
excess is at most .133553653265, even before expansion losses and even
with perfect conversion. It is below 3/8 by .241446346735. This bound
is restricted to that interval; it does not identify an unknown CDL exit
with the barrier. In the fixed potential, reaching the necessary energy
from rest on the other, rising side would require s<=-.07885602698.
That is only an energy threshold, not a derived turning point.

The literal constant-vacuum child starts below V_* by .002233024213,
so it cannot roll to that higher minimum without incoming energy. The
parent difference U_o-V_*=.435266975787 would in principle cover 3/8;
the missing issue is delivery into the child, not a claim that all parent
energy is insufficient. No such delivery is derived by the radius formula.

For ordinary on-shell reheating into the retained massive scalar near
the minimum, m_s=1.477737463 while 2*m_chi=2.476949780e9 in the same
inverse-length units. That one-to-two decay is below threshold. This is
not a universal prohibition on quantum production: no nonperturbative or
curved-background quantum yield is calculated, and no phenomenological
reheating rate is introduced as a new input.

| Route | One-line disposition |
| --- | --- |
| Parke/CDL constant-vacuum wall | R_c and zero seed masses are fixed, but these endpoints are not vacua of the unchanged scalar/volume action. |
| Classical empty matter field | E=0 is invariant, so the E=3/8 reference event never occurs. |
| Volume roll within the computed well | Its maximum excess .13355365 is insufficient, even with perfect conversion. |
| Ordinary pair decay into the retained scalar | The one-modulus-to-two-chi channel is below threshold; another production calculation is absent. |

These conclusions concern the specified routes and do not reject nesting.
There is no admitted reference event for the requested forward cosmology.
C/E0 at that event, today's Omega_k, w_s(z) and the requested SN+DESI
refit are consequently undefined. No reset, imposed matter or fictitious
fit substitutes for the missing state.

## Execution and next connection

The code reuses the existing isolated-bubble and field-action owners.
Independent checks cover Parke's expression, the unsquared junction,
mass balance, the zero-gravity limit, scalar stationarity, empty-state
invariance and the energy bound. This calculation has no field grid,
background integration or new observational fit.

The next compatible construction must supply actual parent/child vacua,
the coupled wall profile and an energy-transfer or reheating mechanism
that can populate the reference state. Those must be derived in one action
before the radius and mass result can be connected to the observation map.
The radius/mass step succeeded within its vacuum-wall assumptions; the
source-admission and matter-production steps did not.
