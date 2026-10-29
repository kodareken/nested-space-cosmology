# Physical selection routes and the evolving-volume observation test

## What can be calculated from the current actions

The five frontier values of Lambda_star*L_ref(0)^2 are 9.702237009265
for (12,2), 9.991573647602 for (6,3), 10.635176064650 for (5,4),
10.776179311559 for (4,5), and 10.727932243144 for (3,7). The
[family theorem](nsc-atlas-scale-selection.md#selection-audit-inheritance-permits-a-family--29-october-2026)
establishes their ordering, not a common action connecting them.

The [route calculation](../results/summaries/nsc-atlas-volume-selection-routes-v1.json)
evaluates the three proposed mechanisms without choosing an extra source
parameter. The existing actions have different d,A,B,T and no field that
changes d or the ancestor label. Their direct sum would preserve those
sectors. An interpolating action needs physical information beyond its
endpoint energies: adding lambda*z^4*(1-z)^4 leaves the endpoint values
and first two derivatives intact while changing a barrier. This is a
nonidentifiability example, not a proposed extra NSC coupling.

| Route | Computed result and disposition |
| --- | --- |
| Attractor | All five matched starts flow toward their own late minimum and retain their own d,j; there is no common flow selecting (12,2). |
| de Sitter entropy | (12,2) has the largest late horizon entropy under common units; the equations supply no inheritance step that maximizes it. |
| Brown--Teitelboim / Coleman--De Luccia | Each cross-member bounce to (12,2) is undefined because its common action and interface are absent; the same (12,2) action also permits a lower-energy runaway direction. |

These dispositions concern these constructions; none rejects nesting.
No parameter-free physical selection was obtained from the supplied laws.

## Five actual attractor evolutions

The existing averaged equations are evolved from N=s=v=0 to N=6 with
the same rho=9/32,E0=3/8,C and each member's fixed coefficients. This is
the strongest directly specified selection test: the whole flow is given.
The respective endpoint s values are approximately -0.03565100,
-0.06486662, -0.06171628, -0.07189444 and -0.09700417. Each is within
3e-5 of its own stationary minimum, and total-energy ledgers close within
1e-8. The discrete labels are fixed coefficients, so their constancy is
exact; a longer evolution cannot change the sector.

The relative late entropies S_i/S_12 are 1, .97104194, .91227799,
.90034109 and .90439022. Reusing
[Gibbons--Hawking](https://doi.org/10.1103/PhysRevD.15.2738), in hbar=c=1
units S=3*pi/(G*Lambda). With g=kappa_nat/L_unit(0)^2 and the owner's
vacuum normalization, S*g=24*pi^2/V_star. The absolute g is still free.
This formula applies to the stationary de Sitter endpoints, not to an
instantaneous matter-filled history. Maximum entropy ranks admissible
endpoints; it supplies neither their accessibility nor a transition law.

## An actual same-action Hawking--Moss coefficient

[Coleman--De Luccia](https://doi.org/10.1103/PhysRevD.21.3305) computes
vacuum decay from an interpolating scalar action including gravity.
[Brown--Teitelboim as reviewed and evaluated by Hirano](https://arxiv.org/html/1804.09985)
uses a specified charged membrane and flux sectors of one theory. Neither
result identifies an interface between the different-dimensional parent
actions here. The four cross-member exponents are therefore undefined,
not zero or infinity; inserting the old wall tension .2 would assume the
missing identification.

There is nevertheless a calculable saddle within the reference action.
Its false vacuum is V_f=.220983024212607, while its barrier is
V_b=.354536677477457. Beyond the barrier V(s) tends to zero from above
as s tends to positive infinity: the internal volume decompactifies.
Thus absence of a lower stable listed member does not establish a terminal
state. Using the standard
[Hawking--Moss saddle](https://doi.org/10.1016/0370-2693(82)90946-7),

$$
B_{\rm HM}=\frac{b_{\rm HM}}{g},\qquad
b_{\rm HM}=24\pi^2\left(\frac{1}{V_f}-\frac{1}{V_b}\right)
=403.7817410357.
$$

The canonically normalized barrier curvature gives
-V_ss/(K*H_b^2)=3.76250309. Scalar harmonics on the Euclidean four-sphere
have eigenvalues l(l+3)*H_b^2+V_ss/K: l=0 is negative while l=1 gives
.23749691*H_b^2>0. This is the single scalar negative-mode regime of
the homogeneous saddle; see the
[published negative-mode criterion and exponent](https://link.aps.org/accepted/10.1103/PhysRevD.95.123502).
The second, heavy scalar is positive at its constant minimum. Other
internal modes and a localized higher-dimensional construction remain
outside the reduction. No determinant, absolute rate or dominance over
all other saddles is calculated. The unknown g prevents an absolute
decay timescale; b is an action coefficient, not a selected SI scale.

## Observation method fixed before the retained refit

Carry (12,2) forward as an explicitly unselected conditional model. Evolve
the actual volume and scalar-dust background, including its energy exchange.
Replacing only 9/32 by V_star in the old Friedmann formula would be a
different model. The numerical homogeneous continuation covers
-log(1+z_max)<=N<=log(20), with s=v=0 at N=0. It does not inherit the
old three-region protected-domain theorem or its full scalar error bound.

The comparison assumes Einstein-frame metric redshift and fixed luminosity
standards, as in the previous distance test. A microscopic photon/detector
completion of the varying-mass model is not derived. For an observer N_o,
N_e=N_o-log(1+z), J_N=exp(-N)/H and

$$
\frac{H_o D_M}{c}=\frac{H_o e^{N_o}}{\sqrt C}
\sin\left(\sqrt C[J(N_o)-J(N_e)]\right),\qquad
\frac{H_o D_H}{c}=\frac{H_o}{H(N_e)}.
$$

Use the same 1322 Pantheon+ rows and full covariance, all 13 DESI DR2
entries and full covariance, one magnitude intercept and one shared r_d*h.
Each model fits one shape (N_o or flat Omega_m), so each has three fitted
quantities. Priors are uniform N_o on [0,log(20)] or flat Omega_m on [0,1],
and the same positive linear-ruler beta prior. Profile both nuisances for
the best fit; integrate the entire shape likelihood and the positive beta
Gaussian for posterior moments. Orders 48/96, with tighter background
tolerance at the higher order, must agree within 1e-5 in reported checks.

The constant late floor and the remaining volume stress must both appear:

$$
\Omega_\Lambda=\frac{V_*}{3H_o^2},\quad
\Omega_m=\frac{E_o}{3H_o^2},\quad
\Omega_k=-\frac{Ce^{-2N_o}}{H_o^2},\quad
\Omega_s=\frac{Kv_o^2/2+V(s_o)-V_*}{3H_o^2}.
$$

These four fractions sum to one. Also report the instantaneous potential
fraction V(s_o)/(3H_o^2), which differs from the constant Lambda fraction.
The supplied 3*Omega_Lambda benchmark is a comparison diagnostic, not an
additional independent likelihood or a validated measurement of this model.

Calculation budget: five finite homogeneous forward solves, one scalar
barrier evaluation, and four independent model/order inference workers
with one BLAS thread each and 180 process-CPU seconds per worker. No
source retuning, new field grid or change to the frozen P004/P005 records.
Failure forecast: missing common-action interfaces leave cross-member
rates undefined; a proposed entropy rule may supply only an ordering;
the actual evolved background may fail the shared-ruler test. Each is
reported separately and redirects the next physical construction.
