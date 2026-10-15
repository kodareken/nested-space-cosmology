# Pre-arrival profile moment and local boundary refinement

The [scalar-input calculation](nsc-atlas-scalar-forcing.md) derived a simpler
mass representation before central arrival. Its
[first interpolation assessment](current-result.md#scalar-lift-and-exterior-forcing-16-october-2026)
found that the active pulse needs additional boundary samples. This successor
checks the mass representation numerically, then uses it for that refinement.
The [admission](prior-art-and-open-claim.md#profile-moment-and-local-lift-admission-16-october-2026)
fixes the comparisons before execution.

The same conformal scalar, source gauge, quiet centre, $h=1$, $b=2$,
$m=0.4$ and finite reference wedge remain in use. All results are formal
coefficients; no finite source amplitude or new physical preparation is
introduced.

## A cumulative first-derivative moment

For $T<q_a$, set $q=T+\mathrm{atanh}(r)$, $A=1-r^2$ and

$$
\mathcal C(q)=\frac12\int_{q_a}^{q}[f'(z)]^2\,dz,\qquad
\nu_E(T,r)=\mathcal C(q)-\frac{Af(q)^2}{4r}.
$$

The moment and profile vanish on their certified earlier quiet interval.
The [existing derivation](nsc-atlas-scalar-forcing.md#a-pre-arrival-profile-moment)
checks both the complete radial kinetic-plus-potential constraint and the
mass-flux equation. The regular quiet interior fixes the integration constant.
This identity does not extend the domain to a nonquiet reflected centre.

Construct one prefix integral in $q$. Its cells split at the profile/source
edges, the arrival images of the field evolution's time cuts, and the arrival
spline knots. For a later query, add only the integral over its last partial
cell. Cache repeated labels and batch first-profile-jet calls within the
existing 513-point limit. The original field evaluator supplies $f,f'$;
stored waveforms are not differentiated.

At each ray point the remaining jets are algebraic:

$$
p=f/r,\qquad p_T=f'/r,\qquad p_r=f'/(rA)-f/r^2.
$$

The cap potential and its outgoing derivative can therefore be evaluated as

$$
K_1=\frac{A\mathcal C}{2r^3}-\frac{Aff'}{12r^3}-\frac{Af^2}{12r^4},
$$
$$
K_{1,u}=\frac{A(3-r^2)\mathcal C}{4r^4}
-\frac{A(3-r^2)ff'}{24r^4}
-\frac{A(2-r^2)f^2}{12r^5},\qquad u=T-\mathrm{atanh}(r).
$$

Here $q$ is fixed when differentiating in $u$. The joint outgoing integral
still uses $J_-=-\int K_1f\,dq$ and its moving-endpoint derivative. The
existing material wall pullback then gives $\mathfrak b=r_Uc_{1,r}|_W$
and $\mathfrak b_V$. Only the mass evaluation changes.

The numerical identity needs an independent check: $f$ comes from a dense
field solution and $f'$ from its evolution equation, so they need not be
exact derivatives of the same floating interpolant. Compare the new mass
against the original radial metric at wall and interior-ray points, and
compare joint boundary values and slopes against the original radial
boundary evaluator. Distinguish moment quadrature sensitivity from a
field/time interpolation discrepancy that persists after refinement.

## A fixed polynomial on every training cell

The new lift accepts explicit knot times and a shared boundary-jet provider.
Every knot set includes the source edges $t=0.25,0.75$ and endpoints $0,1$.
For each training cell, choose four or five nearby training nodes within the
same source block, including both cell endpoints. The stencil is fixed for
the entire cell. A locally normalized confluent Hermite polynomial matches
the values and first derivatives at those nodes.

Adjacent cells share their endpoint values and slopes even when their
stencils differ. Thus the assembled boundary is $C^1$ at the mathematical
level, and its exactly integrated polynomial lift is $C^2$. Floating
endpoint defects are reported separately. Choose integration constants to
make the lift continuous, with $B(U)=0$ for $U\le U_a$ in the retained
domain. This construction preserves $B_U=b_H$ and
$dB(U_W(V))/dV=U'_W(V)b_H(U_W(V))$.

Withheld points remain outside every training stencil. For the local
refinement, use 16/32 cells inside the active interval, and retain four cells
in each quiet/post-pulse block. The fine mesh has 41 training knots. Its
32 active midpoints and four midpoints in each outer block give 40 withheld
queries, for 81 unique direct queries across the nested comparison.
Report boundary-value, cell-width-scaled slope, lift, continuity and
transport discrepancies. These sampled controls do not establish a uniform
enclosure between all points.

The numerical target remains the one derived from the saved field
sensitivity. A successful value fit must also have controlled slope and
lift changes before it supplies the subsequent scalar evolution. A failed
fit directs further local sampling or a diagnosed field-resolution change;
it does not authorize changing the source or couplings.

## Implementation and next connection

The [profile-moment adapter](../src/recursive_horizons/nsc_atlas_profile_moment.py)
and its [tests](../tests/test_nsc_atlas_profile_moment.py) own the faster
pre-arrival queries. The [comparison producer](../scripts/record_nsc_atlas_profile_moment.py)
records the independent radial and moment controls. The
[local lift](../src/recursive_horizons/nsc_atlas_adaptive_lift.py) and its
[tests](../tests/test_nsc_atlas_adaptive_lift.py) own the assembled polynomial
and its exact integral. Numerical outcomes belong in the
[current result](current-result.md).

Once the input accuracy is established, the next equation remains the
previously derived conservative scalar evolution. Its exterior forcing
requires a checked time cache; its cap readout must retain the clock
pullback and regular reflection. An accurate boundary fit alone is not the
propagated scalar signal or a finite-amplitude validity estimate.
