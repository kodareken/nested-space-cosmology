# Schur memory and the local effective Lambda target

For $D=\begin{pmatrix}D_p&\Phi\\\Phi^\dagger&D_c\end{pmatrix}$, eliminating
the exterior gives $\Sigma(z)=\Phi(z-D_c)^{-1}\Phi^\dagger$ and the interior
resolvent $(z-D_p-\Sigma)^{-1}$. The corresponding causal memory carries an
exterior initial-state drive and cross correlations. The intended physical
target is a local effective $\Lambda$ response from this unresolved exterior
structure. The finite identities, spectral controls and geometric candidates
below supply its retained dictionary. The current candidate's causal scalar
boundary map and its local classical metric response are now calculated below.
The finite matrix identity alone does not provide a physical reduction of the
gravitational geometry, a quantum state or a cosmological vacuum term.

## Step 7: averaging the explicit sourced initial slice

The geometric calculation uses the existing
[cap--wall--bridge slice](nsc-atlas-connection-slice.md), from its regular
centre to the retained right-end sphere $R_{\mathrm{cut}}/\ell=4$. This is
the isolated sourced body, not an unevaluated common-parent connector.
The reference length remains the already supplied $\ell$.

Before evaluation, the finite question is whether proper-volume scalar
averaging of that slice supplies a vacuum-form correction beyond its fixed
cap volume term. The inputs remain $H\ell=1$, $k\ell=2$,
$M/\ell=0.4$, $R_w/\ell=1.5$ and the existing smooth vacuum slicing taper.
No observation selects them. The calculation imports the already checked
spatial geometry and wall jump; it derives its scalar averages directly.
[Buchert, gr-qc/0102049v2](https://arxiv.org/html/gr-qc/0102049v2), Sections
2.1--3.2, equations (5), (7), (10b), (12)--(16), supply a comparison for
the averaging identities. His fluid-comoving perfect-fluid assumptions
are not imposed on the moving material wall.

The budget is 30 CPU seconds and 20,000 evaluations of the existing spatial
fields, with no evolution or parameter scan. Fixed Gauss orders 32 and 64
are compared with adaptive quadrature at this one domain. The independent
controls are the integrated Hamiltonian including the shell delta,
direct wall-stress projection, the homogeneous cap limit and unit scaling.
A nonzero density is called a vacuum contribution only if its pressure
also has the required sign and magnitude. A non-vacuum result redirects
the next calculation toward a physical common-parent averaging domain
and its clocks rather than refining this isolated slice.

With the owned convention $K=-\mathcal L_T h/2$, set
$\theta=-\operatorname{tr}K$ and
$\sigma^2=\tfrac12(K^i{}_jK^j{}_i-(\operatorname{tr}K)^2/3)$.
Use $dV=4\pi R^2dl$, retaining both branches of the bridge. Then

$$
Q_D=\tfrac23(\langle\theta^2\rangle_D-\langle\theta\rangle_D^2)
      -2\langle\sigma^2\rangle_D,\qquad
3H_D=\langle\theta\rangle_D.
$$

The intrinsic scalar curvature includes
$-4[R']/R_w\,\delta(l-l_w)$. The wall's integrated slice energy is
$E_S=\sigma\gamma$ and its spatial stress trace is
$S_S=\sigma(\gamma^2-3)/\gamma$. Dropping that curvature delta would
remove a real term from the averaged constraint.

For a Gaussian-normal continuation at the initial instant, with the
domain transported along $T$, the scalar identities are

$$
3H_D^2=\kappa\langle E\rangle_D-
       \tfrac12(\langle{}^{(3)}R\rangle_D+Q_D),\qquad
3\ddot a_D/a_D=-\tfrac\kappa2\langle E+S\rangle_D+Q_D.
$$

This continuation specifies the clock for the instantaneous second
identity; it does not construct the later Einstein--shell solution.
Putting the curvature and variance into a spatially flat scalar fitting
model gives

$$
\kappa\rho_g=-\tfrac12(Q_D+\langle{}^{(3)}R\rangle_D),\qquad
\kappa p_g=-\tfrac12 Q_D+\tfrac16\langle{}^{(3)}R\rangle_D.
$$

Their trace projection is
$\Lambda_{g,\mathrm{tr}}=(Q_D-\langle{}^{(3)}R\rangle_D)/4$;
the vacuum condition is $3Q_D+\langle{}^{(3)}R\rangle_D=0$.
Scalar fitting is not a tensor averaging prescription or an observer
redshift--distance calculation. The supplied cap vacuum and the material
wall remain explicit in $E$ and $S$.

### Finite result and the next parent input

The [compact calculation](../results/summaries/nsc-atlas-geometric-average-v1.json)
gives $V_D/\ell^3=355.4931547895$ and

$$
H_D\ell=0.08352417231,\quad Q_D\ell^2=0.48101493580,\quad
\langle{}^{(3)}R\rangle_D\ell^2=0.35791137287.
$$

Consequently $\kappa\rho_g\ell^2=-0.41946315433$ and
$\kappa p_g\ell^2=-0.18085557242$. The trace projection is
$\Lambda_{g,\mathrm{tr}}\ell^2=0.030775890732$, but the necessary vacuum
condition has residual $(3Q_D+\langle{}^{(3)}R\rangle_D)\ell^2
=1.80095618027$. This correction is therefore not a vacuum source.
The largest 32/64/adaptive difference in the trace projection is
$1.0\times10^{-10}$; this is a quadrature comparison, not a certified or
physical-model uncertainty. The wall-inclusive Hamiltonian residual is
$1.2\times10^{-16}$ in these units. Seven focused controls and the
eighteen existing slice controls pass. No physical trajectory is run.

An outward extension can be assessed without scanning cutoffs. Define the
compact-body moments
$I_B=\int_B\theta\,dV$, $J_B=\int_B(2\theta^2/3-2\sigma^2)\,dV$,
and $\mathcal R_B=\int_B{}^{(3)}R\,dV$, including the wall.
On the owned time-symmetric right end, $K=0$ and ${}^{(3)}R=0$; enlarging
the isolated domain changes its volume but none of these moments. Hence

$$
\Lambda_{g,\mathrm{tr}}=
 \frac{J_B-\mathcal R_B}{4V_D}-\frac{I_B^2}{6V_D^2}.
$$

The isolated correction vanishes as the exact AE domain grows. This
analytic limit follows from the supplied vacuum end, rather than a
numerical extrapolation to cosmological scales.

For $N$ unchanged disjoint bodies and a hypothetical time-symmetric vacuum
connector, the vacuum condition would require volume per body
$v_* = 2I_B^2/(3J_B+\mathcal R_B)=23.1714109008\,\ell^3$.
The retained domain of just one body already occupies
$355.4931547895\,\ell^3$. Adding positive connector volume cannot meet
that condition. Joining copies through a $K=0$ vacuum connector is therefore
insufficient under this scalar averaging prescription. The existing
general gluing theorem does not assert a time-symmetric connector;
this is a conditional calculation, not a property assigned to its output.

For a vacuum connector with nonzero $K$, its constraint instead gives
$\mathcal R_C=-J_C$. The smallest additional scalar data are its proper
volume and the two moments $I_C,J_C$, evaluated with a physically selected
normal clock. For one retained body the vacuum test becomes

$$
V_D(3J_B+\mathcal R_B+2J_C)-2(I_B+I_C)^2=0.
$$

Choosing these moments to satisfy the equation would fit a desired form;
they must come from a compatible parent geometry and observer prescription.
The existing existence family does not supply them. This identifies the
necessary Step-5 input to Step 7 without requiring a new formation or
stability campaign. A common-parent development and tensor/observer
interpretation remain necessary before a cosmological Lambda is predicted.

For the local interior target there is a further, already established
[observer boundary](nsc-atlas-causal-classification.md#observer-access-preserve-direction-and-the-emission-event).
In the isolated reference an inward wall signal has finite central reception
only for $HR_w<1$; this averaging slice uses $HR_w=1.5$. Its parent volume
is not a causal coarse-graining domain for that central observer. Parent
moments would complete the parent scalar average, not identify a new central
vacuum term. That identification requires a physical reduction tied to the
retained observer and its compatible earlier data or quantum state. Global
causal exclusions from the isolated reference are not assigned to the
unevaluated many-body development.

### Unequal sourced inputs for the parent reduction

The [checked input pair](../results/summaries/nsc-atlas-unequal-children-v1.json)
now replaces the second vacuum bridge with its own regular cap and material
wall, in the same theory. The
[slice owner](nsc-atlas-connection-slice.md#unequal-sourced-input-18-october-2026)
specifies its geometry, boost and source checks. The following moments use
the already declared instantaneous proper normal frame and the right-end
retention sphere `R/ell=4`. Write
`I=int theta dV`, `J=int(2 theta^2/3-2 sigma^2)dV` and
`mathcal_R=int R3 dV`, retaining the wall curvature distribution.

| Retained input | $V/\ell^3$ | $I/\ell^2$ | $J/\ell$ | $\mathcal R/\ell$ |
| --- | ---: | ---: | ---: | ---: |
| A, mass `.4*ell` | 355.493155 | 89.076815 | 185.877651 | 127.235043 |
| B, mass `.6*ell` | 397.394632 | 106.535682 | 207.136386 | 106.880419 |
| Sum before connector | 752.887787 | 195.612497 | 393.014037 | 234.115462 |

A's moments are preserved from its closed adaptive calculation. B uses
Gauss 64, compared once with Gauss 32; the largest absolute moment difference
is `1.35e-7` in the respective ell units. Values in the table are rounded.
The comparison indicates quadrature sensitivity rather than a total error
bound or an observational uncertainty.

These are the inputs to a parent, not its complete volume or average. If
`V_B,I_B,J_B,mathcal_R_B` denote the summed body moments, a selected vacuum
connector would give `mathcal_R_C=-J_C`, and the instantaneous normal-clock
vacuum-form condition would be

$$
(V_B+V_C)(3J_B+\mathcal R_B+2J_C)-2(I_B+I_C)^2=0.
$$

`V_C,I_C,J_C` remain absent in the record. They must follow from a chosen
compatible connector, with the clock terms added for a different observer
prescription. Solving this equation for a desired connector moment would
select the answer. The two input bodies establish neither a local
cosmological vacuum source nor the cap centre's access to the parent volume.

### Two successive sourced levels

The [two-level record](../results/summaries/nsc-atlas-nested-slice-v1.json)
supplies a different hierarchical body, `A2`, at the same outer mass
`.4*ell`. It contains an inner membrane and Kottler host of mass `.01*ell`;
both membranes use `k*ell=2`. The cap and host vacuum scales remain supplied.
This extends the geometry under the same declared law; it is not a pure
change of resolution of the previous one-body state.

Through the same right-end cut `R/ell=4`, including both wall distributions,

$$
V/\ell^3=360.7573878173,\quad I/\ell^2=100.2281955340,\quad
J/\ell=204.0438671729,\quad \mathcal R/\ell=207.6312086409.
$$

The initial normal-clock projection gives

$$
\Lambda_{g,\mathrm{tr}}\ell^2=-0.0153506330754,\qquad
(3Q+\langle R^{(3)}\rangle)\ell^2=2.11796240047.
$$

The latter would have to vanish for vacuum form. The geometric density and
pressure are `kappa*rho_g*ell^2=-.544841233192` and
`kappa*p_g*ell^2=-.161146233630`; they do not obey `p=-rho`.
Gauss 64-minus-32 differences are `-1.28e-10` in the trace and
`-3.91e-10` in the vacuum residual, as numerical indicators. The integrated
Hamiltonian residual is `-3.33e-16` in ell units. These indicators are not
rigorous, model or observational uncertainty bounds.

Adding a checked level changes the finite response with the same regional
law and couplings. It does not make the scalar projection a tensor vacuum
source or select the observed cosmological scale. `A2` also has its
[central test-flash observable](nsc-atlas-connection-slice.md#a-central-test-flash-protocol).
Its exact end can replace A in the
[conditional two-child assembly](nsc-atlas-connection-gluing.md#retaining-a-grandchild-in-the-parent-initial-data).
The evaluated connector, common observer prescription and effective local
gravitational reduction remain the outward calculation. This finite probe
requires no further quadrature refinement.

### Physical resolution of one fixed geometry

The inherited relation describes a common law. It does not by itself force
every admissible state to be self-similar. A concrete resolution map is
needed to test Douglas's proposal that successive unresolved structure
produces increasingly complex local behaviour. Adding the grandchild above
changes a state. The following construction instead reduces one fixed
state of the separate Einstein/Weyl/Dirac candidate.

On its reference spatial metric `h=r^2(Q^2 dx^2+dOmega^2)`, radial scalars
have the positive operator

$$
\mathcal L=-\Delta_h=-\frac{1}{r^3Q}\partial_x\left(\frac{r}{Q}\partial_x\right).
$$

Use the owned odd geometry basis and its existing even quadrature to form
`K_ij=integral (r/Q) phi_i' phi_j' dx` and
`M_ij=integral r^3 Q phi_i phi_j dx`. Select generalized eigenvalues
`mu*L_ref^2<=1` on the stored `T=.05` slice, with the previously fixed
observer ruler `L_ref=4.80028044552352`. Orthonormalizing this span in the
canonical `dx` pairing gives a fixed symmetric projector `P`. The physical
metric selects the span; the canonical pairing selects its complement.
The 2048-point quadrature is still a numerical representation of the same
511 geometry variables, not an extra physical environment.

Apply `P` to all geometry pairs in the exact nodal canonical coordinates
`a=log(r), pi_a=r*p_r`, `b=log(Q), pi_b=Q*p_Q`, and `chi,p_chi`.
Keep the full fermion columns and occupations. Nodal reconstruction by
exponentiation must also pass positivity after prolongation. At fixed
reference projector and canonical environment `z_e=(1-P)z`, the difference

$$
H_{\rm corr}(z_c;z_e)=H(z_c+z_e)-H(z_c),\qquad
\delta F_c=P X_H(z)-P X_H(Pz)
$$

is an instantaneous Hamiltonian correlation. Here `P` acts as the identity
on retained fermions, so the environment contains geometry only. The
canonical log-momentum rates include the product rule,
`pi_a_dot=r*p_r_dot+p_r*r_dot`, and likewise for `b`. Both evaluations use
the same owned Hamiltonian, quadrature, gauge, couplings and matter state.
This first comparison uses only `T=.3`; it is not a retarded elimination.
Rates refer to the same conformal coordinate-time generator. Filtering also
changes the metric lapse; a matched full/coarse proper-clock reduction has
not been constructed by this instantaneous comparison.

The vacuum comparison follows independently from the four-dimensional
volume. The owned `(+---)` action has `-A R_L-V4`, with `R_E=-R_L`.
Thus the dictionary `A=1/(16*pi*G)`, `V4=2*A*Lambda` gives

$$
H_\Lambda=8\pi A\Lambda\int r^4LQ\,dx,
\qquad \delta V=-8\pi A\Lambda r^4.
$$

The existing `r^2` term in `V` is sphere curvature and cannot replace this
volume generator. Varying before choosing `L=Q` gives the unit-Lambda lapse
constraint `8*pi*A*Q*r^4` and momentum rates
`p_Q_dot=-8*pi*A*L*r^4`, `p_r_dot=-32*pi*A*L*Q*r^3`.
After the conformal gauge choice the total `p_Q` rate is
`-16*pi*A*Q*r^4`. No coordinate, auxiliary-momentum, shift-constraint or
fermion-rate change comes from this potential. Adjoint pullback to the
owned geometry variables precedes the log-canonical change and projection.

The coefficient `delta_H/H_unitLambda` is an energy projection to compare
with the **entire** retained vector, not a derivation of vacuum form. Each
canonical block has its own residual and units; mixing their unweighted
norms would obscure the test. The original full-state constraint defect,
filtered constraint and exact saved-rate replay remain separate. The bare
vacuum term and missing quantum vacuum branch are not introduced. The
[bounded admission](prior-art-and-open-claim.md#physical-resolution-of-one-fixed-geometry-18-october-2026)
permits one reference eigensolve and the two final-state evaluations. No
fractal dimension, repeated-scale flow or observed cosmological constant
has yet been derived from this map.

**Finite outcome.** The [8,875-byte record](../results/summaries/nsc-conformal-resolution-v1.json)
uses frozen method `cd00b35`, fifteen manufactured controls and 7.703125
process-CPU seconds. Both BLAS libraries use eight threads. The cutoff
retains the constant geometric mode out of 511; the first omitted
eigenvalue is `.4795862384393`, above the fixed `.04339770650775` cutoff.
The relative eigen residual is `1.65e-15`, projector defect `1.84e-16`, and
the retained constant has relative residual `1.80e-12`. Replaying the
archived prolonged rates agrees to `8.72e-9` in the largest momentum-rate
absolute difference and `2.23e-15` in the radius rate. These are numerical
controls, not a propagated physical error range.

At `T=.3`, `H_full=-1.363e-11`, `H_coarse=1.015478298556`, and the unit
volume generator is `231.9284384454` in the corresponding action units.
Thus the instantaneous energy coefficient is

$$
\Lambda_E L_{\rm ref}^2=-0.100890404750047.
$$

The whole vector fails the vacuum-volume comparison. Its `pi_a` residual
RMS is `.2956761085522`, or `58.23%` of the larger correlation/candidate
RMS; `pi_b` gives `.2406867478252`, or `94.81%`. The retained `a_dot`
correlation RMS is `.0009888445527823` per coordinate time and the
`chi_dot` correlation is `.6727649529294`; a vacuum potential produces
neither rate. The fermion-rate difference is also nonzero, despite holding
the input columns fixed. The `b_dot` difference is only roundoff scale
(`5.54e-18` RMS) and is not treated as a detected physical correction.

The original full Hamiltonian-constraint maximum is `9.005e-6`, while the
filtered state's is `25.5817144`; its momentum constraint likewise changes
from `7.99e-7` to `6.16370736`. These are unreset defects in their own
constraint-density units. The filtered metric with the retained fine matter
is not a second constraint-satisfying spacetime. Its missing correlations
must supply more than a volume term. The calculation establishes a
nonzero resolution response in this finite model and closes this
single-snapshot probe. A causal elimination with the corresponding
constraint and kinetic/matter operators is the next physical connection;
the energy coefficient is not a generated cosmological vacuum constant.

### Observer clocks in the geometric average

The selected larger-geometry route must compare resolutions in a declared
clock. For zero shift, let $d\tau_o=N_o(t)dt$ name that observer's time and
$\nu=N(t,x)/N_o(t)>0$. Transport a compact domain with the slice normals,
using its proper spatial measure. Then $\partial_{\tau_o}dV=\nu\theta dV$.
Differentiating the volume and using the normal Raychaudhuri equation gives

$$
3H_o=\langle\nu\theta\rangle,
\qquad Q_\nu=\tfrac23\bigl(\langle\nu^2\theta^2\rangle
-\langle\nu\theta\rangle^2\bigr)-2\langle\nu^2\sigma^2\rangle,
$$
$$
\mathcal R_\nu=\langle\nu^2{}^{(3)}R\rangle,
\qquad P_\nu=\langle(\partial_{\tau_o}\nu)\theta\rangle
+\langle\nu D^2\nu\rangle.
$$

With the same Einstein representation and sources as above, including all
surface terms where present, the scalar identities are

$$
3H_o^2=\kappa\langle\nu^2E\rangle
-\tfrac12(\mathcal R_\nu+Q_\nu),\qquad
3a_D''/a_D=-\tfrac\kappa2\langle\nu^2(E+S)\rangle+Q_\nu+P_\nu.
$$

The corresponding spatially flat scalar fitting coefficients become

$$
\kappa\rho_g=-\tfrac12(Q_\nu+\mathcal R_\nu),\qquad
\kappa p_g=-\tfrac12Q_\nu+\tfrac16\mathcal R_\nu-\tfrac23P_\nu,
$$
$$
\Lambda_{g,\mathrm{tr}}=\tfrac14(Q_\nu-\mathcal R_\nu+2P_\nu),
\qquad 3Q_\nu+\mathcal R_\nu+2P_\nu=0
\quad\hbox{is the vacuum-form condition.}
$$

The coefficient of $\nu D^2\nu$ is **one**. For the acceleration
$a_i=D_i\log N$, its four-divergence is
$\nabla_\mu a^\mu=D^2N/N$, so $\nu^2\nabla_\mu a^\mu=\nu D^2\nu$.
Integration by parts retains boundary flux:
$\int\nu D^2\nu\,dV=\int_{\partial D}\nu D_n\nu\,dA
-\int|D\nu|^2dV$. A mean lapse alone cannot recover this term.

An independent exact control is the static periodic metric
$ds^2=-[1+\epsilon\cos x]^2dt^2+dx^2+dy^2+dz^2$, $|\epsilon|<1$,
with the clock anchored where $N=1$. Direct connection contraction gives
$R_{nn}=N''/N$ and $\kappa S=2N''/N$, while $E=\theta=\sigma=R_3=0$.
Thus $P_\nu=-\epsilon^2/2$ cancels the averaged Ricci source exactly:
the volume is static. Doubling the lapse term would falsely predict
$3a_D''/a_D=-\epsilon^2/2$. A homogeneous changing-clock control independently
recovers $3a_D''/a_D=3[\nu' H+(\nu H)^2]$ for $a=e^{Ht_{\rm proper}}$.

The [existing average calculator](../src/recursive_horizons/nsc_atlas_geometric_average.py)
now evaluates these weighted moments as well as its original geodesic case.
The finite check uses at most 30 CPU seconds for the two analytic controls
and existing tests, with no new physical field solve or grid. This completes
the missing clock identity for the chosen averaging route. It does not select
the parent metric, observer/domain or state, and its scalar fitting terms are
not by themselves a local four-tensor or a cosmological Lambda prediction.

### Averaging an existing evolving connected metric

The next application reuses the separately specified conformal common-action
realization already evolved on periodic $S^1\times S^2$. Its full metric is

$$
ds^2=r^2\{Q^2[-dT^2+dx^2]+d\Omega^2\},\qquad x\sim x+8.
$$

Thus the physical lapse and radial scale are both $rQ$, rather than $Q$.
The original observer is the fixed normal worldline $x=1$ with
$d\tau_o=r(T,1)Q(T,1)dT$. Its Einstein/Weyl-auxiliary/Dirac action differs
from the P001 Einstein/conformal-scalar action; this is an evaluation of the
existing connected candidate, not a construction of P001's absent parent.
It does not establish a fractal metric or unlimited composition.

The archived NPZ has now been recovered under ignored `.build/data`, with
its original SHA256 verified. Use the four stored time/space cases at the
predeclared times $T=0.05,0.2,0.3$. The primary case is
`nf512_dt_0.0005`; $L_{\rm ref}$ is its observer sphere's areal radius at
$T=0.05$. The [admission](prior-art-and-open-claim.md#evolving-geometric-average-17-october-2026)
fixes the domain, clock, ruler, 30-CPU-second allowance and stopping conditions
before the averages are evaluated. No evolution or new field solve is needed.

The actual second rates follow from the same saved projected dynamics.
Write $\Pi_g=A_g(n_g/n_q)A_g^{\mathsf T}$, $F_r=-8\pi A r$,
$F_\chi=\mathrm{constant}$ and $Z=-24\pi A$. Differentiating the owned
velocities gives

$$
Q_{TT}=\Pi_g\frac{Q_Tp_\chi+Qp_{\chi,T}}{2F_\chi},\qquad
r_{TT}=\Pi_g\frac{p_{r,T}-(F_{r,T}p_\chi+F_rp_{\chi,T})/F_\chi}{2Z},
\qquad F_{r,T}=-8\pi A r_T.
$$

All directions on the right are the saved lifted projected rates. This
retains the projection after multiplication as well as the projected state
direction. Uniform exponential and truncated Fourier controls check these
jets independently; time samples are not differentiated.

The scalar average reconstructs a metric-defined effective Einstein source.
Relating it to the actual matter and moved higher-curvature terms inherits
the original finite field-equation defects. Its Hamiltonian identity is
therefore not an independent check of those sources. The saved finite
constraint defects stay separately reported.
The output keeps the geometric trace, vacuum-form defect, direct volume
acceleration and refinement indicators together. A trace alone supplies no
cosmological vacuum identification or observational agreement.

The [saved-geometry result](../results/summaries/nsc-conformal-geometric-average-v1.json)
now evaluates that question. Its reference ruler is
$L_{\rm ref}=4.80028044552352$ in the original action units. The primary
case gives the following dimensionless coefficients:

| Stored time | $H_oL_{\rm ref}$ | $\Lambda_{g,\rm tr}L_{\rm ref}^2$ | $(3Q_\nu+\mathcal R_\nu+2P_\nu)L_{\rm ref}^2$ |
| --- | ---: | ---: | ---: |
| 0.05 | -0.00488136 | -0.58208154 | 2.22336490 |
| 0.2 | -0.01941849 | -0.57363917 | 2.23211408 |
| 0.3 | -0.02901471 | -0.57251629 | 2.17744464 |

At the last time, $\mathcal R_\nu L_{\rm ref}^2=2.24615801$,
$Q_\nu L_{\rm ref}^2=-0.01240311$ and
$P_\nu L_{\rm ref}^2=-0.01575203$. The geometric density and pressure are
$\kappa\rho_g L_{\rm ref}^2=-1.11687745$ and
$\kappa p_g L_{\rm ref}^2=0.39106257$. Positive spatial curvature dominates
this correction, and the domain is contracting. It does not supply a
vacuum-form or accelerating cosmological response on the tested interval.

The largest four-case deviations of the final trace and vacuum defect are
respectively $4.65\times10^{-9}$ and $1.92\times10^{-8}$. These are stored
refinement indicators, not continuum or physical uncertainty bounds.
Metric identities agree within $1.19\times10^{-14}$ dimensionlessly, and
the saved first-rate reconstruction within $3.61\times10^{-15}$.
The primary original-action Hamiltonian defect remains between
$8.35\times10^{-6}$ and $9.01\times10^{-6}$ at these times; it was not
replaced by the metric-defined Hamiltonian identity. The entire application
used 0.46875 CPU seconds and 737,280 charged values, with zero field solves
or evolution steps. Fourteen metric/jet controls pass.

This closes the finite averaging probe for this candidate. It does not
provide the desired Lambda mechanism. The next physical construction needs
the actual common-parent/transition geometry and its observer reduction;
the present periodic spherical realization supplies neither by repetition.
No additional precision campaign is opened here.

### Causal radial distances and the observer's clock

An averaged scalar coefficient is not yet a distance an observer measures.
The next finite calculation uses actual null paths in the recorded conformal
metric, keeping the existing observer $x_o=1$ and reception $T_o=0.3$.
Normal test emitters and reflectors stay at fixed coordinate positions;
their worldlines may be accelerated. Transparent geometric optics and
photon-number conservation define the light protocol. They are probes,
not an assumed cosmological matter congruence.

With $N=rQ$, the past ray is $T=T_o-\chi$, $x=x_o+\eta\chi$, $\eta=\pm1$.
The null Hamiltonian makes the radial covariant momentum $p_T$ constant.
Frequency measured by a normal observer is $|p_T|/N$, hence

$$
1+z=\frac{N_o}{N_e}.
$$

Radial screen shear vanishes. The warped-sphere curvature gives
$R_{kk}=-2r_{\lambda\lambda}/r$, so $r$ itself solves the scalar Jacobi
equation. Normalize the past affine parameter by observed frequency
magnitude one. Then $d\lambda/d\chi=N^2/N_o$, and the vertex conditions
$D_A(0)=0$, $D_{A,\lambda}(0)=1$ give

$$
D_A=\frac{r_or_e}{N_o}\int_0^\chi
Q(T_o-s,x_o+\eta s)^2\,ds,
\qquad D_L=(1+z)^2D_A.
$$

The factor $r_o$ matters: this observer has finite nonzero areal radius.
The observer can also send an ideal pulse, reflect it at $x_o+\eta\chi$
and receive it at $T_o$. Its radar distance follows from its own proper
round-trip clock,

$$
D_{\rm radar}=\frac12\int_{T_o-2\chi}^{T_o}N(T,x_o)\,dT.
$$

The ratio $D_L/D_{\rm radar}$ needs no choice of length units. Fix
$\chi=0.025,0.075,0.125$ in both directions before evaluating this ratio.
The longest radar trip starts at $T=0.05$, inside the archived interval.
Use the four existing metric cases, their real Fourier spatial fields and
cubic Hermite interpolation with the saved time derivatives. No source
or geometry is evolved, and no waveform samples are differentiated.

The [existing distance calculator](../src/recursive_horizons/nsc_atlas_cosmological_distance.py)
now supplies the radial cone. Its controls recover Minkowski annular
distances, the off-centre flat-de-Sitter law, and a separately integrated
Sachs equation. The [saved-geometry reader](../scripts/measure_nsc_conformal_lightcone.py)
also checks the Fourier/Hermite representation and the proper radar clock;
for de Sitter its exact radar distance is
$\log(1+2H\chi)/(2H)$ when the reception scale factor is one.
The [admission](prior-art-and-open-claim.md#causal-distance-and-clock-comparison-17-october-2026)
fixes the finite domain, probes, compute ceiling and 32/64-point consistency
pair. Their differences do not bound geometry interpolation or model error.
This supplies a causal observer map for this separate candidate; it does
not supply the missing common parent, a full sky or a cosmological Lambda.

The [frozen-method evaluation](../results/summaries/nsc-conformal-lightcone-v1.json)
now gives the following primary-case observations. The same radar distance
applies to both radial directions at each fixed ray length:

| $\chi$ | $D_{\rm radar}/L_{\rm ref}$ | $z$, smaller $x$ | $z$, larger $x$ | $D_L/D_{\rm radar}$, smaller $x$ | $D_L/D_{\rm radar}$, larger $x$ |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0.025 | 0.00626688 | 0.00000486563 | -0.00170720526 | 0.99921752 | 0.99750604 |
| 0.075 | 0.01881968 | 0.00064380534 | -0.00468723591 | 0.99913797 | 0.99381120 |
| 0.125 | 0.03138286 | 0.00166774473 | -0.00730610408 | 1.00007398 | 0.99110851 |

Thus the recorded contracting geometry can give a redshift in one direction
and a blueshift in the other for these probes. Clock gradients and optical
focusing enter actual observations; the volume-average expansion alone does
not specify them. This is a prediction of the fixed finite numerical metric
and probe protocol, rather than an observation of nature or a Hubble fit.

The largest 32/64-point observable difference is $3.34\times10^{-16}$.
The largest four-case deviations are $2.94\times10^{-10}$ in redshift and
$4.84\times10^{-10}$ in the luminosity/radar ratio. They do not include an
interpolation or total propagated-field error bound. The primary original
Hamiltonian defect at reception remains $9.01\times10^{-6}$ in its own
source normalization. The method was frozen in `f0f639d` before evaluation.
It used 11.6875 CPU seconds and 64,252,080 charged coefficient values, with
no field solve or evolution step. Eleven optical/interpolation controls pass.

![Radial redshifts and luminosity/radar distance ratios for the recorded normal probes](figures/nsc-conformal-lightcone.png)

Markers show the three declared ray lengths; lines connect those samples.
The [plotter](../scripts/plot_nsc_conformal_lightcone.py) reads the compact
result without replaying the metric. This completes the finite forward
observation calculation. A cosmological observer/state, larger connected
geometry and its derived effective dynamics remain the further requirements.

### The wall connection between the two action choices

The inherited-law programme needs a physical relation between the Einstein
cap/wall candidate and the separately specified Weyl-squared realization.
Here the finite question is whether the quiet material wall can be imported
literally, with its existing sources. This is not the same question as deriving
a coarser effective interface from a smooth geometry.

Use signature $(-+++)$, $R^a{}_{bcd}=\partial_c\Gamma^a{}_{db}
-\partial_d\Gamma^a{}_{cb}+\Gamma^a{}_{ce}\Gamma^e{}_{db}
-\Gamma^a{}_{de}\Gamma^e{}_{cb}$, a spacelike common normal, and
$[K_{ab}]=-k h_{ab}$. The singular curvature coefficients are
$R_{anbn}^{\rm sing}=k h_{ab}$, $R_{ab}^{\rm sing}=k h_{ab}$,
$R_{nn}^{\rm sing}=3k$ and $R^{\rm sing}=6k$, each multiplying the
surface delta. Substitution in the **linear** four-dimensional Weyl
decomposition cancels every singular Weyl component. This computation does
not manipulate undefined products of deltas or establish a full quadratic
thin-shell calculus.

The regular Weyl tensor still steps from the conformally flat cap to the
Schwarzschild exterior. Direct Christoffel/Riemann calculation, with zero
exterior Ricci checked independently, gives in an orthonormal wall frame

$$
[C_{anbn}]=\frac{M}{R_w^3}\operatorname{diag}(-2,-1,-1),
\qquad h^{ab}[C_{anbn}]=0.
$$

A radial Lorentz boost leaves these components unchanged. At the already
owned $M/\ell=0.4$, $R_w/\ell=1.5$ event, the dimensionless jump is
$\ell^2[C_{anbn}]=\operatorname{diag}(-32/135,-16/135,-16/135)$.
It scales as inverse length squared, not inverse length cubed.

The principal normal-derivative part of the Bach operator
$(\nabla^\rho\nabla^\sigma+R^{\rho\sigma}/2)C_{\mu\rho\nu\sigma}$
therefore has formal double-layer strength $[C_{\mu n\nu n}]$.
In a controlled thin-layer limit with bounded metric/connection and bounded
Weyl curvature approaching the two one-sided values, its leading weak action
on a test tensor is

$$
-\int_\Sigma[C_{\mu n\nu n}]\,n^\rho\nabla_\rho Y^{\mu\nu}\,d\Sigma.
$$

This is the derivative-of-delta part, before the action's nonzero Weyl
coefficient is applied. Ordinary Einstein and tension terms have only delta
order. The quiet scalar is zero; a continuous scalar with bounded one-sided
first derivatives also supplies no such double layer through its improved
bulk stress. The existing Euler/Box ledger completes global divergences on
smooth periodic metrics and finite time caps; it supplies no additional
timelike-wall source. Their completed global terms are not adjustable Weyl
dipole sources.

This is a **necessary leading-order compatibility check**, conditional on
the stated thin-layer control. Lower-order products and a complete junction
require a regularization or resolved transition. It does not prove that a
smooth connected solution is impossible. It shows why the current literal
tension wall cannot simply be identified with an interface of the
Weyl-squared candidate using the supplied data. If the wall is instead an
effective boundary between resolutions, the reduction must supply this
curvature response and its physical scale.

[Reina--Senovilla--Vera, 1510.05515v1](https://arxiv.org/html/1510.05515),
Sections 3.1--3.2, discusses the hazards of formal quadratic cancellations
and restricts the generic distributional class. That generic theorem is not
applied to this special trace jump; the tensor calculation above is direct.
The [existing gravitational calculator](../src/recursive_horizons/nsc_atlas_gravitational_reduction.py)
now exposes the regular normal Weyl components. Seven
[checks](../tests/test_nsc_atlas_gravitational_reduction.py) include an
independent coordinate-curvature derivation, the boost, all 256 singular
Weyl components, units and the zero-mass limit. No new field solve, wall
trajectory or full junction theorem is claimed.

## Step 7: the current candidate's first vacuum projection

The current conformal-scalar cap–wall candidate gives the exact smooth-bulk
result

$$
\boxed{\Delta\Lambda_{\rm classical}\,L_{\rm ref}^2=0.}
$$

This is the additional trace projection of its classical scalar stress. The
cap already contains the supplied term $\Lambda_{0,-}\ell^2=3(H\ell)^2=3$;
the exterior contains $\Lambda_{0,+}=0$. Those baseline values were inputs
to the geometry. The zero above holds for any positive reference length,
so it requires no choice of physical ruler or observed cosmological value.
For displaying the baseline, $L_{\rm ref}=\ell$ uses the existing reference
units; the P001 value $\ell=10000\,\mathrm m$ remains a fiducial scale.

Reuse the [owned action and trace identity](nsc-atlas-cap-wall-response.md),
with $\Lambda_0=\kappa_4 V_0$, rather than extending the cap computation:

$$
T^\phi_{\mu\nu}=\nabla_\mu\phi\nabla_\nu\phi
-\tfrac12g_{\mu\nu}(\nabla\phi)^2
+\tfrac16\left[G_{\mu\nu}\phi^2+
(g_{\mu\nu}\Box-\nabla_\mu\nabla_\nu)\phi^2\right],
$$
$$
T^\phi{}^\mu{}_\mu=\phi(\Box-\mathcal R_J/6)\phi=0,
\qquad \mathcal R_J=4\Lambda_0,
\qquad
\Delta\Lambda_{\rm classical}L_{\rm ref}^2
=\tfrac14(\mathcal R_J-4\Lambda_0)L_{\rm ref}^2=0.
$$

This is exact on shell, including nonzero classical pulses. The pulse still
has nonzero stress, mass transport and wall response: a zero trace does not
make its stress tensor zero. In the decomposition
$T^\phi_{\mu\nu}=-(\Delta\Lambda/\kappa_4)g_{\mu\nu}
+T^{\rm TF}_{\mu\nu}$ with $g^{\mu\nu}T^{\rm TF}_{\mu\nu}=0$, the
vacuum projection is zero. This calculation concerns the unaveraged smooth
bulk. It does not assign a four-dimensional vacuum density to the material
wall, perform nonlinear metric averaging, or include quantum stress.

The [executable calculation](../src/recursive_horizons/nsc_atlas_effective_vacuum.py)
contracts the full improved stress and checks the result independently in
the Einstein frame. With $q=\Psi/\sqrt6$, $F=\operatorname{sech}^2q$ and
$\widehat V=\Lambda_0\cosh^4q$, use
$\mathcal R_E=6(\nabla_Eq)^2+4\widehat V$ and
$\Box_Eq=\widehat V_{,q}/6$. Then

$$
\mathcal R_J=F\left[\mathcal R_E+3\Box_E\log F
-\tfrac32(\nabla_E\log F)^2\right]=4\Lambda_0.
$$

The apparent field-dependent Einstein-frame potential therefore supplies
no extra Jordan-frame vacuum curvature. Seven
[checks](../tests/test_nsc_atlas_effective_vacuum.py) also reconstruct this
curvature from the actual nonlinear bulk operator's metric and field jets,
in both regions at zero and nonzero amplitude. The exact result has no
quadrature error; these checks do not estimate errors in an evolved solution.
Reproduce with
`python scripts/run_python.py -m recursive_horizons.nsc_atlas_effective_vacuum`
and `python scripts/run_python.py -m pytest tests/test_nsc_atlas_effective_vacuum.py -q`.
No endpoint replay, new spatial grid or observational fit enters this calculation.

For Step 7 this resolves the classical scalar route to an additional smooth
vacuum term in the current candidate. A nonzero contribution needs an
explicit additional mechanism, such as a specified quantum effective action
or a defined geometric averaging operation. The existing finite Schur
matrices have not been derived from this candidate's exterior. The causal
field map and its conformal boundary metric check are now supplied below;
the full quantum state and general metric response remain additional inputs.
In particular, the zero-frequency
coefficient of a scalar self-energy multiplies a field quadratic term;
it is not by itself a metric volume term or a cosmological constant.
The older two-sheet and Gauss–Bonnet results below retain their own actions.

## Step 7: conditional quantum curvature contribution

The next calculation imports one **quantized real massless conformal scalar**
on the quiet geometry, conditional on a Hadamard state with the usual
short-distance singularity and locally covariant renormalization. It uses the
four-dimensional field, including all angular sectors; the classical `l=0`
pulse alone cannot supply its vacuum fluctuations. The calculation applies
inside each smooth region, away from the material wall. It does not quantize
the metric or select a joint quantum state across the attachment.

For the scalar Laplace operator with $E=-\mathcal R/6$, the imported
[heat-kernel coefficient, equation (4.28)](https://arxiv.org/pdf/hep-th/0306138)
contains $(5\mathcal R^2+60\mathcal R E+180E^2)/360=0$ and
$(12\Box\mathcal R+60\Box E)/360=\Box\mathcal R/180$.
Including the $(4\pi)^{-2}$ factor gives the
[conformal trace coefficient, equations (15)–(16)](https://arxiv.org/pdf/hep-th/9308075):

$$
\langle T^\mu{}_{\mu}\rangle
=\frac{\hbar}{2880\pi^2}
  (\mathcal R_{\mu\nu\rho\sigma}^2-\mathcal R_{\mu\nu}^2)
+\hbar b\Box\mathcal R.
$$

The finite curvature-counterterm choice can change $b$. Both quiet metrics
have constant $\mathcal R$, so this ambiguity vanishes here. The renormalized
volume and Newton terms remain supplied action coefficients. The imported
heat-kernel and quantum-stress rules are assumptions of this extension;
the algebraic cancellation and geometric contractions are calculated locally.

Use the existing reference unit $L_{\rm ref}=\ell$, and set
$\zeta=\kappa_4\hbar/\ell^2$, $h=H\ell$, $m=M/\ell$, $r=R/\ell$.
Direct contraction for
$ds^2/\ell^2=-f\,dt^2+dr^2/f+r^2d\Omega^2$ gives

| Quiet region | $f$ | $\ell^4\mathcal R_{\mu\nu\rho\sigma}^2$ | $\ell^4\mathcal R_{\mu\nu}^2$ |
| --- | --- | --- | --- |
| Cap | $1-h^2r^2$ | $24h^4$ | $36h^4$ |
| Exterior | $1-2m/r$ | $48m^2/r^6$ | $0$ |

The local contractions also apply in the Schwarzschild interior; no global
Euclidean continuation or stationary state is inferred from this chart.
Consequently the trace projection $\Delta\Lambda_{\rm tr}=-\kappa_4\langle T^\mu{}_{\mu}\rangle/4$ is

$$
\frac{\Delta\Lambda_{{\rm tr},-}\ell^2}{\zeta}
=\frac{h^4}{960\pi^2},\qquad
\frac{\Delta\Lambda_{{\rm tr},+}\ell^2}{\zeta}
=-\frac{m^2}{240\pi^2r^6}.
$$

At the already supplied $h=1$, $m=2/5$ and quiet wall radius $r=11/20$,
these coefficients are **0.0001055428996274352** and
**−0.002440241780407079**, respectively. They multiply $\zeta$; neither
is an absolute cosmological Lambda. No cosmological data selected the
reference length or the coefficients.

The cap trace is constant, but that alone does not establish a stress
proportional to the metric: the wall can break the symmetries of the state.
The exterior trace projection varies as $r^{-6}$ and requires the rest of
the stress for conservation. A joint state, its traceless stress and boundary
exchange, the gravitational renormalization conditions, and a controlled
backreaction calculation are still needed for a full effective Lambda.
In particular, this calculation does not derive the supplied cap scale $H$.
It isolates a calculable quantum contribution while the classical result
above remains zero. No additional cap-grid evolution enters either result.

For the already declared P001 fiducial scale $\ell=10000\,\mathrm m$,
restore SI units with $\kappa_{\rm SI}=8\pi G/c^4$ and the quantum trace
normalization $\hbar c$. Thus

$$
\zeta=\frac{8\pi G\hbar}{c^3\ell^2}
=6.56537648967\times10^{-77}.
$$

Using [CODATA 2022 inputs](https://physics.nist.gov/cuu/pdf/wall_2022.pdf),
$G=6.67430(15)\times10^{-11}\,\mathrm{m^3\,kg^{-1}\,s^{-2}}$,
exact $h_{\rm Planck}=6.62607015\times10^{-34}\,\mathrm{J\,s}$ and
exact $c=299792458\,\mathrm{m\,s^{-1}}$, gives

| Quiet-region trace projection | $\Delta\Lambda_{\rm tr}\ell^2$ | Standard uncertainty from $G$ alone |
| --- | ---: | ---: |
| Cap at $h=1$ | $6.92928871865\times10^{-81}$ | $1.5573\times10^{-85}$ |
| Exterior at $m=2/5$, $r=11/20$ | $-1.60211060142\times10^{-79}$ | $3.6006\times10^{-84}$ |

These are scale illustrations conditional on the fixed geometry and the
fiducial length, not cosmological predictions or total error intervals.
Keeping $h,m,r$ fixed while doubling $\ell$ divides the dimensionless
contribution by four and the dimensional $\Delta\Lambda$ by sixteen.
The supplied cap baseline remains $\Lambda_0\ell^2=3$; this quantum
contribution does not explain that input.

The same [calculation module](../src/recursive_horizons/nsc_atlas_effective_vacuum.py)
derives the conformal cancellation and contracts Christoffel curvature for
both metrics. An independent sectional-curvature contraction checks the
result. The [focused tests](../tests/test_nsc_atlas_effective_vacuum.py)
cover these identities, normalization and flat limits as well as the
classical calculation and SI scaling. The commands above reproduce the exact expressions
and their displayed decimal coefficients; there is no quadrature error in
these symbolic results.

### Does this term generate the cap curvature?

The same quantum extension supplies a necessary test without a new evolution.
Let $X=(H\ell)^2$, $\lambda_0=\Lambda_0\ell^2$ and
$a=\zeta/(960\pi^2)$. The constant-curvature trace equation is

$$
3X=\lambda_0+aX^2.
$$

Set the supplied renormalized volume coefficient to $\lambda_0=0$ explicitly.
The roots are

$$
X=0\quad\hbox{or}\quad X=\frac{2880\pi^2}{\zeta}.
$$

Since $\zeta=8\pi\ell_P^2/\ell^2$, the nonzero root obeys
$H^2\ell_P^2=360\pi$. It lies outside the small-curvature regime
$H^2\ell_P^2\ll1$ in which neglecting quantum-gravitational corrections is
controlled. Thus this one-scalar anomaly does not generate a nonzero,
small-curvature de Sitter cap when the supplied volume term is zero.
This is a consequence of the displayed equation, not an assumption that a
preferred cosmological model must hold.

For a supplied $\lambda_0\geq0$ and $9-4a\lambda_0\geq0$, the lower branch
is $X=2\lambda_0/(3+\sqrt{9-4a\lambda_0})$ and tends to $\lambda_0/3$ as
$\zeta\to0$. It corrects a supplied curvature. Neither branch establishes
the full tensor equation, a joint state or a compatible quantum wall. The
calculation concerns a constant-curvature smooth cap with the selected free
scalar and renormalized gravitational coefficients; a nonconstant geometry
or another physical mechanism is a separate calculation.

The executable `constant_curvature_trace_condition()` checks the roots and
classical limit. A twelfth test derives the nonzero root directly from the
dimensional SI trace equation, independently of the dimensionless conversion.
No extra grid, fitted scale or observational input enters this test.

### Quantum interface: the actual transmission law

The material interface supplies a separate operator contribution. For the
quiet candidate, $[\phi]=0$ and $[n\phi]=k\phi$. These are the transmission
conditions of $D=-(\nabla^2+E)+U\delta_\Sigma$, with $E=-\mathcal R/6$ and
$U=k$. Use the [Gilkey–Kirsten–Vassilevich transmittal coefficients](https://arxiv.org/html/hep-th/0101105v1)
as imported local polynomials; do not square distributional curvature.
Their theorem concerns an elliptic operator on Riemannian manifolds. The
substitution below checks its invariant algebra on the Lorentzian quiet-wall
geometry. It does not construct a global Euclidean continuation, a Lorentzian
quantum state or a renormalized interface stress.

With the project's common normal from cap to exterior, the source's inward
forms are $L^+=-K_+$ and $L^-=K_-$. Thus
$S=L^++L^-=kI$, $D=L^+-L^-=-2\bar K$. The second heat coefficient's
surface density is proportional to $\operatorname{tr}S/3-U=0$.
For unit smearing, the fourth coefficient also contains curvature
contractions; keeping only its extrinsic-curvature terms would give an
incomplete answer.

In the source index convention $\rho_{ij}=R_{ikkj}$, put
$Q_{ac}=R_{abcb}$, with the tangential metric used for contractions. The
Gauss equation on the shared induced geometry gives

$$
Q_+-Q_-=\frac{k}{2}\big[(\operatorname{tr}D)I+D\big],\qquad
\operatorname{tr}(Q_++Q_-)
=-\mathcal R_+-\mathcal R_-+2\rho_+(n,n)+2\rho_-(n,n)
=-6H^2.
$$

The first identity cancels the difference-curvature and difference-shape
terms in Theorem 7.1's $\mathcal A_1$. The normal-curvature terms in
$\mathcal A_2$ cancel because $S=kI$. This uses no static Schwarzschild
observer and remains a contraction identity inside its horizon. The
remaining three groups are

$$
\mathcal A_1=0,\qquad
\mathcal A_2=56k^3-12kH^2,\qquad
\mathcal A_3=-60k^3.
$$

After the factor $1/[360(4\pi)^2]$, the local surface polynomial is

$$
a_4^\Sigma(1)=\int_\Sigma\sqrt{|\gamma|}\,q_\Sigma\,d^3y,
\qquad q_\Sigma\ell^3
=-\frac{b(b^2+3h^2)}{1440\pi^2}
=-\frac7{720\pi^2}\quad(h=1,\ b=2).
$$

The Schwarzschild mass and radius cancel in the full expression. This is
an interface coefficient, not a four-dimensional volume term. Its value
on this background does not establish that its off-shell metric variation
is merely a shift of the material tension. The full counterterm operators,
their renormalized coefficients and the joint-state stress remain needed
before obtaining a physical effective Lambda from the attachment.

There is an independent boundary-curvature check. For outward curvature
$K$, define
$J=\operatorname{tr}K\operatorname{tr}(K^2)
-\frac23\operatorname{tr}(K^3)-\frac13(\operatorname{tr}K)^3$.
The [Euler boundary polynomial](https://arxiv.org/html/hep-th/0208205v3),
equations (6)–(7) after conversion from inward normals, is
$4(J-2\widehat G_{ab}K^{ab})$.
Adding the cap and exterior sides gives
$B_\Sigma=4[J(K_-)-J(K_+)+k\widehat{\mathcal R}]
=4k(k^2+3H^2)$, using the intrinsic Gauss scalar.
The scalar's Euler coefficient $-1/[360(4\pi)^2]$ reproduces $q_\Sigma$.
This checks the local polynomial and its orientation; it introduces no
new gravitational coupling or global topological assertion.
The calculation is implemented by `conformal_interface_heat_coefficients()`
in the same module. Its checks and the independent Euler-boundary test bring
the focused module to fourteen passing tests.

### Connecting the actual exterior to the Schur operation

The operator to eliminate in this candidate is the exterior conformal wave
operator with the transmission domain above. Its boundary variable is
$\varphi=\phi|_\Sigma$. Define $\mathcal N_\pm$ to map that boundary value
to the corresponding **outward** normal derivative, with the initial or outer
data stated separately. The exterior outward normal is $-n$ and the cap's
is $n$. Therefore the homogeneous boundary equation is

$$
(\mathcal N_-+\mathcal N_++U)\varphi=0.
$$

It is exactly $[n\phi]=U\varphi$ in the common-normal convention. Supplied
incoming data add a driving term to this equation. In real time the maps
are causal history operators; the finite preparation does not justify a
stationary zero-frequency limit.

For an auxiliary elliptic problem with specified outer boundary data, write
the exterior field as its Dirichlet solution plus a fluctuation vanishing
on the interface. Integration by parts removes their bulk cross term.
Gaussian elimination then separates the boundary quadratic form
$\frac12\langle\varphi,\mathcal N_+\varphi\rangle_\Sigma$ from
$\frac\hbar2\operatorname{Tr}\log D_{+,D}$, where $D_{+,D}$ has homogeneous
Dirichlet interface data. The boundary measure is the induced geometric
measure, not a freely normalized matrix inner product.

The first term changes field transmission; the second and the metric
dependence of the first contribute to the effective stress. A numerical
coefficient of the first term alone is not a vacuum volume density. Real-time
expectation values additionally require the joint state's symmetric
two-point function, along with the retarded response. The existing prepared
classical pulse specifies a mean field, not that covariance. This identifies
the missing Step 3 input needed by Step 7: a physical quantum preparation
and its covariance, with compatible boundary exchange and renormalization.
The old finite Schur matrices remain useful identities, but are not yet
evaluations of these physical operators.

### Exterior response to a prescribed wall field

The [prepared pulse](nsc-atlas-prepared-response.md) used the particular
condition $\chi_U|_W=0$ selected by its quiet cap and no-return interval.
The exterior operator needed above instead prescribes the wall value
$y(V)=\chi(U_W(V),V)$. Keep the same quiet geometry and finite wedge,
with source value $g(V)$ at $U=U_0$ and bottom value $b(U)$ at $V=V_a$.
Here $\chi=r\phi$ and all coordinates use the existing reference unit.

For a spherical harmonic of degree $l$, direct reduction of the wave
operator gives

$$
\chi_{UV}+K_l\chi=0,\quad
K_l=-\frac{r_{UV}}r+\frac{C\,l(l+1)}{4r^2},\quad
C=\frac{32m^3}r e^{-r/(2m)},
$$
$$
K_l=\frac{16m^4}{r^4}e^{-r/(2m)}
 +\frac{8m^3l(l+1)}{r^3}e^{-r/(2m)}.
$$

The exterior is Ricci flat, so its conformal-curvature coupling contributes
no additional bulk potential. Define the entry time
$v_0(U)=V_a$ for $U\le U_a$ and $v_0(U)=U_W^{-1}(U)$ otherwise, and
$a(U)=b(U)$ or $y(v_0(U))$ on those respective intervals. The Dirichlet
integral equation is

$$
\chi(U,V)=a(U)+g(V)-g(v_0(U))
-\int_{v_0(U)}^Vdq\int_{U_0}^Udp\,K_l(p,q)\chi(p,q).
$$

Every integration point is retained because $U_W(q)\ge U$ when
$q\ge v_0(U)$. Differentiating in $V$ and then $U$ recovers the bulk
equation. The empty integrals recover all three prescribed faces when
$b(U_0)=g(V_a)$ and $y(V_a)=b(U_a)$.

Writing $I_l(V)=\int_{U_0}^{U_W(V)}K_l\chi\,dU$, the first wall jets are

$$
\chi_V|_W=g'-I_l,\qquad
\chi_U|_W=\frac{y'-g'+I_l}{U_W'}.
$$

At the wall-bottom corner a common first jet additionally requires
$y'(V_a)=U_W'(V_a)b'(U_a)+g'(V_a)-I_l(V_a)$.
Further smoothness requires the corresponding higher-jet transport
conditions. Boundary derivatives are supplied or evaluated through these
equations, rather than estimated from stored waveform samples.

For the proper-time tangent $u=(\dot U,\dot V)$, with
$C\dot U\dot V=1$, the exterior outward normal is
$\nu_+=(\dot U,-\dot V)$. Consequently

$$
\mathcal N_+^R\varphi+\text{initial-data drive}
=\dot U\left(\frac{\chi_U}{r}-\frac{yr_U}{r^2}\right)
-\dot V\left(\frac{\chi_V}{r}-\frac{yr_V}{r^2}\right),
\qquad \varphi=y/r.
$$

The radius terms follow from differentiating the physical scalar and are
essential. This determines the actual causal field-response map associated
with the exterior and its induced boundary measure. Its quantum covariance
and metric variation are additional ingredients for effective stress.
For fixed geometry and identical supplied traces, a uniform field error
$e_\chi$ gives a conditional wall-normal error bound
$|\delta(\nu_+\phi)|\le 2\dot V K_{l,\max}(U_W-U_0)e_\chi/r$.
This follows from $|\delta I_l|\le K_{l,\max}(U_W-U_0)e_\chi$;
it supplies no bound on the numerical field error itself.

Both terms in $K_l$ decrease over $0.4\le r\le0.55$. The inherited analytic
area bound therefore gives
$q_l\le q_0[1+l(l+1)/2]$, with $q_0=0.13882099234357653$.
The bounds for $l=0,1,2,3$ are respectively
$0.138821,0.277642,0.555284,0.971747$, all below one. This proves uniqueness
for each of these mode problems in the finite continuous wedge and bounds
the inverse by $1/(1-q_l)$. It does not supply an angular-tail estimate or
a four-dimensional quantum determinant.

The admitted implementation checks this operator with at most 120 CPU
seconds and 30 million kernel values, using existing geometry and product
integration. It checks the zero-kernel limit, the normal and radius factors,
corner compatibility and recovery of the already solved prepared pulse as
Dirichlet input. Numerical product-integration and geometry errors remain
distinct from the analytic contraction. No source selection, new physical
trajectory, parameter scan or state calibration enters those checks.
The [implementation](../src/recursive_horizons/nsc_atlas_exterior_dtn.py)
uses a fixed-U grid with analytic entry values at the moving wall. Each
time segment owns its active nodes, including the endpoint Runge-Kutta
stages, so a new node retains its prescribed value until entry.
Six [controls](../tests/test_nsc_atlas_exterior_dtn.py) pass, including a
constant physical scalar on the curved exterior and preservation of entry
data. That one 129-node control gives maximum outward derivative
$5.993\times10^{-6}$ instead of its exact zero; it is a numerical diagnostic,
not an error enclosure. The existing 129-node prepared pulse is recovered
within the stated numerical tolerances, with no new physical input.

### Metric response of the classical exterior reduction

This calculation uses that actual exterior map and the original
unit-amplitude prepared source, zero bottom data and prepared wall value.
It completes the boundary metric identification for the same conformal
action; it does not select a new quantum state. The split action's scalar
boundary term is $-K\phi^2/6$. Its variational boundary operator is therefore
$\mathcal B_+=\mathcal N_++K_+/3$, where here $K_+$ uses the exterior's
outward normal. Retain initial-data driving separately when calling its
retarded response a Schur operator.

The predeclared limits were one fixed 129-node prepared solve and one
129-node Dirichlet solve, with their existing combined ceiling of 120 CPU
seconds and 30 million kernel values. Boundary integrals use 32/64-point
Gauss quadrature on each existing time panel. Independent checks use the
signed wall curvature, directly transformed scalar/metric jets and exact
conformal weights. Stop on disagreement or budget exhaustion; no grid scan
or new evolving trajectory is required. The output is a compact response
summary and hashes, rather than field history.

The signed reference motion gives the two outward mean curvatures
$K_-\ell=9/4$ and $K_+\ell=15/4$. The mass-dependent terms cancel in
$\beta'(r)+2\beta/r$, including at the regular $\beta_+=0$ crossing.
Thus $(K_-+K_+)/3=k$, recovering the actual scalar transmission law.
Under a smooth Weyl change $g'=e^{2\omega}g$, with all supplied scalar
data transformed with weight $-1$,

$$
\mathcal N'_+\phi'=e^{-2\omega}
 [\mathcal N_+\phi-(\nu_+\omega)\phi],\qquad
K'_+=e^{-\omega}(K_++3\nu_+\omega),\qquad
\mathcal B'_+\phi'=e^{-2\omega}\mathcal B_+\phi.
$$

The $K/3$ term cancels the normal derivative of the Weyl factor. The
boundary measure has weight $+3$, so its pairing with two weight-$-1$
fields and this response is invariant. This is a boundary Ward identity;
the trace of a boundary stress at fixed $\phi$ need not vanish. A driven
retarded map is also not in general self-adjoint. Its pairing alone is not
an ordinary single-history action or a quantum determinant.

The specific metric question is the additional **local smooth bulk** vacuum
term induced by this classical elimination. Variations supported inside
the retained cap, keeping the interface/exterior data fixed, have no support
in the exterior boundary functional. Boundary-supported response, the
supplied cap vacuum term, and a quantum determinant remain separate objects.
A numerical boundary coefficient is not identified with a four-dimensional
volume term. The calculation also tests the conformal boundary Ward identity
with nonzero normal derivatives, so a zero bulk contribution cannot be
attributed to accidentally evaluating a zero source.

The [compact result](../results/summaries/nsc-atlas-exterior-schur-metric-v1.json)
gives

$$
\boxed{\Delta\Lambda_{\mathrm{Schur,cl,local}}L_{\mathrm{ref}}^2=0,}
\qquad L_{\mathrm{ref}}=\ell.
$$

This exact zero concerns the specified local retained-bulk variation of the
classical exterior influence. It does not remove the supplied cap term
$\Lambda_{0,-}\ell^2=3$ or identify a new cosmological value. The actual
nonzero response control is
$\int_\Sigma\phi\mathcal B_+\phi\,d\Sigma/\varepsilon^2
=0.3703318455512892$ over the original finite strip. The bare-normal pairing
is $0.0804611703811659$. These dimensionless boundary diagnostics use the
constant-on-sphere monopole; neither is a vacuum density.
The 32/64-point quadrature differences are respectively
$5.0\times10^{-16}$ and $4.25\times10^{-14}$. The sampled Weyl-pairing
scaled residual is $9.16\times10^{-16}$, and the independent outward radial
orientation residual is $2.92\times10^{-9}$ from the interpolated wall.
These diagnostics are not continuum field-error bounds. Eight focused
checks pass, including an independent transformed-connection contraction
and the signed transmission sum. Evaluation used 6.5625 CPU seconds and
2,312,034 charged kernel values. The original evolving jobs are unchanged.

A nonzero quantum value needs the eliminated operator's determinant and
the joint state's symmetric covariance, with a declared finite vacuum
subtraction and metric variation. The causal mean-field response supplies
neither by itself. A bulk-volume term from a constrained change of the whole
atlas would instead need that physical metric/observer reduction specified.
This calculation settles the direct classical local contribution and
identifies the remaining Step-7 inputs rather than fitting a target Lambda.

### Joint quantum preparation: the actual Cauchy data needed

On a spacelike slice crossing the wall let $T$ be the continuous future
unit normal, $s$ the outward spatial unit normal, and
$u=\gamma(T+vs)$, $n=\gamma(vT+s)$. Put $\Pi=T\phi$, $D=s\phi$.
Continuity along the moving wall and the owned transmission law give

$$
[\Pi]+v[D]=0,\qquad
\gamma(v[\Pi]+[D])=k\phi,
\qquad
\boxed{[D]=\gamma k\phi,\quad[\Pi]=-\gamma v k\phi.}
$$

The inverse boost $T=\gamma(u-vn)$, $s=\gamma(n-vu)$ independently gives
the same equations. At the existing initial slice's $R_w/\ell=1.5$ event,
`wall_slice_data()` gives $\gamma=1.0092620299352366$ and
$\gamma v=0.13641790596910147$. Hence
$\ell[D]/\phi=2.0185240598704732$ and
$\ell[\Pi]/\phi=-0.27283581193820294$. A static interface condition would
miss the momentum jump. These values come from the same quiet reference
trajectory as the earlier prepared wedge; the different slice events do
not themselves create a geometric incompatibility.

For any two transmitted solutions the symplectic flux jump is exactly
$[f\,n(g)-g\,n(f)]=kfg-kgf=0$. Thus there is no independent scalar wall
symplectic degree of freedom in this action. With proper slice volume,
the canonical commutator is
$[\phi(x),\Pi(y)]=i\hbar\delta_\Sigma(x,y)$. A proposed connected Gaussian
Cauchy covariance $\Gamma$ must satisfy

$$
\Gamma=\Gamma^{\mathsf T},\qquad
\Gamma+\frac{i\hbar}{2}J\ge0,\qquad
SJS^{\mathsf T}=J,\qquad
\Gamma_t=S\Gamma_0S^{\mathsf T},
$$

where $J$ is the canonical commutator matrix, with proper-volume smearing,
and $S$ is the joined Cauchy evolution. The two-point kernel must obey the
moving-wall conditions in both arguments. For a local renormalized stress,
the subtracted bulk remainder must have the coincidence derivatives used
by the stress operator. Finite-order stress regularity and a strict smooth
Hadamard condition are distinct requirements; neither follows from an
instantaneous principal-symbol prescription. The exact $C^3$ pulse remains
the classical mean, and no stronger quantum regularity is assigned to it.

These equations specify the missing preparation object; they do not select
or construct $\Gamma_0$. The full joined causal propagator fixes its
commutator once that problem is well posed. The existing finite-angular
exterior map supplies only part of that propagator. An initial covariance
on disjoint subregions must be checked for transmission and short-distance
regularity after evolution; it is not the same as imposing vanishing cross
correlations at all times. Symbolic inversion, the numerical boost and an
independent Grok derivation agree on the equations above.

Constructing this covariance would determine the state-dependent quantum
tensor and flux. It cannot change the earlier one-scalar anomaly trace or
derive the supplied $\Lambda_{0,-}$. A nonzero cosmological value still
requires a physical reduction of the gravitational geometry tied to observers,
or a specified source of nonzero trace. Further state work addresses that
tensor and its observer interpretation; it does not supply a new scalar trace.

### Common-parent route check

An exact cosmological reference is available from
[Mars--Mena--Vera, 1307.4371v1](https://arxiv.org/abs/1307.4371), Section 2:
comoving Schwarzschild cavities match a dust FLRW parent when
$M=4\pi\rho R_b^3/3$. The parent's Lambda is supplied on both sides of
that outer matching. Choosing flat dust with zero Lambda therefore fixes
the exterior observers' $q=1/2$; it does not derive a new vacuum response.
This is a reference control, rather than the next effective-Lambda model.

For a right-hand vacuum collar, the owned spherical constraints permit
$R_l=\sqrt{1-2M/R+V^2}$, $K_{ll}=-V_R$, $K^\theta{}_\theta=-V/R$ and
$R_3=-4VV_R/R-2V^2/R^2$. A flat-ended interpolation from $V=0$ to
$V=H_pR$ can join the existing quiet body to the cosmological slice when
$2M=H_p^2R_b^3$. Those initial identities still leave the physical parent
state and observer reduction to be selected. The dust source would be an
explicit addition to the current scalar model. No such addition, second
child evolution or averaged cosmological Lambda is claimed here.

The existing vacuum-end gluing remains reusable for initial-data existence.
A parent calculation for the Lambda target must instead supply a compatible
metric and a physically selected clock/domain, then calculate its observable
response without prescribing the desired expansion beforehand. Implementing
the dust reference would not resolve that missing connection.

### Whole-body gravitational reduction at the parent end

The next finite calculation varies the metric on a static tube in the
already specified right Schwarzschild end. It eliminates the whole quiet
cap--wall--bridge body as seen by a parent observer. This is a different
retained region from the cap-centre scalar reduction above. The same source
has $M/\ell=0.4$ and the existing outer cut is $R/\ell=4$; no physical
population density or parent expansion is supplied. Scalar preparation and
P001's finite-amplitude trajectories are unchanged.

Use the Einstein--Hilbert action with its GHY term and a Minkowski boundary
subtraction, with energy conjugate to **proper boundary time**.
[Brown--York, gr-qc/9209012v1](https://arxiv.org/abs/gr-qc/9209012),
Section III, equations (3.14)--(3.19), and Section VI, equations
(6.2)--(6.20), supply the boundary variational method. Their two-sphere
extrinsic-curvature sign differs from the positive outward convention used
here; the resulting energy and surface pressure agree. The formulas below
are checked directly from the extrinsic curvature and by differentiating
the energy, rather than accepted from the citation alone. The remaining
paper is not audited.

In $G=c=1$ units, with $f=1-2M/R>0$ and normal
$n=\sqrt f\,\partial_R$, the tube's mixed extrinsic curvature is

$$
K^a{}_b=\operatorname{diag}\left(
\frac{M}{R^2\sqrt f},\frac{\sqrt f}{R},\frac{\sqrt f}{R}\right).
$$

Subtract the flat tube with $K_0^a{}_b=\operatorname{diag}(0,1/R,1/R)$.
The boundary response
$\tau^a{}_b=-[K^a{}_b-K\delta^a{}_b-(K_0^a{}_b-K_0\delta^a{}_b)]/(8\pi)$
then has $\tau^a{}_b=\operatorname{diag}(-\varepsilon,P,P)$, where

$$
E(R)=4\pi R^2\varepsilon=R(1-\sqrt f),\qquad
P=\frac1{8\pi R}\left(\frac{1-M/R}{\sqrt f}-1\right).
$$

These obey

$$
dE=-P\,d(4\pi R^2)+\frac{dM}{\sqrt f},\qquad
M=E-\frac{E^2}{2R},\qquad
E=M+\frac{M^2}{2R}+O(M^3/R^2).
$$

At the declared cut, $E/\ell=0.4222912360003365$ and
$P\ell=0.0000619768235605939$. These are a finite-tube gravitational energy
and **surface** pressure. The mass monopole retained by the far parent is
the ADM limit $M/\ell=0.4$. Neither $P$ nor $E(R)-M$ is a four-dimensional
vacuum density. In SI units, $M$ in these equations means $GM_{\rm SI}/c^2$;
the physical tube energy and surface pressure are their displayed geometric
values multiplied by $c^4/G$.

For well separated, slowly moving, conserved quiet bodies, matching this
monopole gives the leading covariant parent action
$S_{\rm bodies}=-\sum_i M_i\int d\tau_i$. This is a long-wavelength
approximation, not a constructed common-parent solution. In a comoving
homogeneous domain with proper volume $V=a^3V_c$ and lapse $N$, it is
$S=-\int dt\,N E(V)$ with $E=\sum_iM_i$. Lapse and spatial-scale variations
give the geometric energy and pressure densities

$$
\rho_g=E/V,\qquad p_g=-\partial E/\partial V=0.
$$

Thus the leading mass-monopole source is pressureless. Its trace projection
is nonzero, $\Lambda_{\rm tr}=2\pi\rho_g$, while its vacuum-form defect is
$8\pi(\rho_g+p_g)=8\pi\rho_g$. For identical reference bodies, writing the
unsupplied physical number density as $n$, these coefficients are
$\Lambda_{\rm tr}\ell^2=2.5132741228718345\,(n\ell^3)$ and
$8\pi(\rho_g+p_g)\ell^2=10.053096491487338\,(n\ell^3)$.
There is no generated volume-Lambda term at this leading conserved-monopole
order. A dust trace is not identified with one.

For a general conserved homogeneous energy functional, vacuum form instead
requires $\partial E/\partial V=E/V$, hence $E=\rho_\Lambda V$ with
constant $\rho_\Lambda$. This derivative must come from varying the physical
metric-dependent action. Setting $\sum_iM_i\propto a^3$ by hand makes the
density constant but does not supply that variation: with zero mechanical
pressure the source is $Q=\dot E/V$ in
$\dot\rho+3H\rho=Q$. The exchanging sector must also be included.

This settles the leading parent mass response without selecting a desired
expansion. It does not settle tidal, mutual-interaction, velocity or
finite-size terms. The small quantities $M/d$, body size over separation,
body size over parent curvature scale, and squared velocity describe the
approximation; no numerical remainder bound or allowed population is proved.
A nonzero cosmological vacuum term needs further physical dependence of the
reduced action on the parent metric, with its state and energy exchange.
The body's expanding interior alone does not specify that functional.

The finite budget is 30 CPU seconds for exact identities and one reference
evaluation, with independent boundary-tensor, energy-derivative, unit and
homogeneous-action checks. No evolution or spatial-resolution scan is part
of this calculation.
The [calculator](../src/recursive_horizons/nsc_atlas_gravitational_reduction.py),
[five checks](../tests/test_nsc_atlas_gravitational_reduction.py) and
[compact result](../results/summaries/nsc-atlas-gravitational-reduction-v1.json)
record that finite calculation. The reference evaluation uses about 0.74 CPU
seconds; all exact identities and the focused tests pass.

### A distance observable for the Lambda target

The quiet reference cap supplies one immediate observer prediction. With
zero scalar source and comoving test emitters in its flat de Sitter chart,
$L_{\mathrm{ref}}=c/H$ and $D_L/L_{\mathrm{ref}}=z(1+z)$. This follows
by integrating the radial null equation and including redshift of energy
and arrival rate. It tests the existing quiet geometry, not the nonzero
P001 pulse or a newly generated Lambda.

Use the central event at cap time zero, where the existing reference wall
has $R_w/\ell=1.5$. For $1\le R_w/\ell\le1.5$, its signed dynamics bound
$(dR_w/d\tau)^2\le4569/1600<4$, and $dt/d\tau=\gamma\ge1$.
It therefore needs at least $\ell/(4c)$ to reach $R_w/\ell=1$ backwards
from that event. The redshifts below have lookback at most
$\log(1.25)\ell/c<\ell/(4c)$ and emission radius at most $0.2\ell$.
Their null paths remain inside the quiet cap. The test-source luminosities
are probes; a matter/radiation cosmology has not been supplied by this limit.

Fix the following comparison before calculating its value. Use the official
[Pantheon+ release](https://github.com/PantheonPlusSH0ES/DataRelease) at commit
`c447f0fea703fcd0fff57de5000947b5ca81286b`, the standardized `m_b_corr`
column and full `STAT+SYS` covariance. Retain non-calibrator light curves
in `0.025 <= zHD < 0.075` or `0.15 <= zHD <= 0.25`, preserving the
released row order and correlated repeated light curves. No Cepheid-anchored
distance modulus or diagonal-only errors enter the calculation.

The [released SN-only likelihood](https://github.com/PantheonPlusSH0ES/DataRelease/blob/c447f0fea703fcd0fff57de5000947b5ca81286b/Pantheon%2B_Data/5_COSMOLOGY/cosmosis_likelihoods/Pantheon%2B_only_cosmosis_likelihood.py)
uses $D_L=(1+z_{\mathrm{HD}})(1+z_{\mathrm{HEL}})D_A(z_{\mathrm{HD}})$.
For the quiet cap this is proportional to
$z_{\mathrm{HD}}(1+z_{\mathrm{HEL}})$. Define

$$
y_i=m_{b,i}^{\mathrm{corr}}-5\log_{10}
 [z_{\mathrm{HD},i}(1+z_{\mathrm{HEL},i})],\quad
B=(X^TC^{-1}X)^{-1}X^TC^{-1},\quad
w^T=(-1,1)B,\quad \Delta=w^Ty,
$$

where $X$ contains the two band indicators and $C$ is their joint covariance.
Then $w^T\mathbf1=0$, eliminating the common magnitude/distance zero point,
and $\sigma_\Delta^2=w^TCw$. The **quiet-cap prediction is exactly
$\Delta=0$**, or $\mathcal R=10^{\Delta/5}=1$. A later forward model must
use these same weights on its individual redshift predictions; the statistic
is not a distance ratio evaluated at band midpoints.

The finite calculation reads one pinned table and covariance, with a
120 CPU-second ceiling. Raw data stay in the ignored cache; the saved product
contains aggregate results, source hashes and reproduction commands.
Independent checks use a constrained covariance solve, synthetic correlated
data and exact zero-point cancellation. Hash, dimension, covariance or
solver failure stops this comparison without changing its bands.
Only antisymmetric roundoff below `1e-12*max(1,max(abs(C)))` is averaged;
larger asymmetry or failed positivity stops the covariance solve. No diagonal
jitter or extra error term is introduced.

This is a retrospective consistency check on published standardized data,
not a blinded or withheld-observation test. It adopts the release's
standardization and covariance, including its fiducial peculiar-velocity
corrections; it fits no LambdaCDM parameters. Its error range is conditional
on that measurement model. A discrepancy constrains this quiet test state
and its emitter assumptions; it does not reject every nested geometry.

The [first evaluation](../results/summaries/nsc-atlas-distance-contrast-v1.json)
uses 596 light curves: 335 in the lower band and 261 in the higher band.
The method and zero prediction were committed as
[`723ee8e`](https://github.com/kodareken/nested-space-cosmology/commit/723ee8e3c24f6b729e10a31c1e99b4af3dbe8833)
before numerical evaluation. The result is

$$
\Delta=-0.1046180474751224\ \mathrm{mag},\qquad
\sigma_\Delta=0.014373986469318333\ \mathrm{mag}.
$$

Its standardized contrast is $-7.27829$ under the published covariance.
The transformed contrast is $\mathcal R=0.9529637714530853$, with
one-standard-deviation interval $[0.9466764929149757,0.9592928064641945]$.
The quiet prediction $\mathcal R=1$ does not match this statistic. Changing
the common Hubble or luminosity scale cannot alter that prediction, because
the weights sum to zero. The full matter, expansion and observer dynamics
must therefore be supplied before this quiet building block is used as a
cosmological forward model. This does not determine the missing induced
Lambda or validate/falsify the nonzero-source P001 configuration.

The source table and full covariance match the pinned publisher Git objects.
The independent constrained solve differs in weights by at most
$8.07\times10^{-16}$; the common-zero-point residual is below
$9\times10^{-16}$. Five synthetic geometry/covariance controls pass.
A separate table parse and Cholesky--QR calculation agrees well within the
$10^{-10}\,\mathrm{mag}$ reproduction tolerance for both the contrast and
its standard uncertainty, with the same selected rows and counts.
Evaluation costs 1.125 CPU seconds. Only the 4,975-byte summary is retained
in Git; the two raw inputs remain in the ignored cache. The released
measurement assumptions remain part of the comparison's uncertainty scope.

## The transition and vacuum-form coefficient

The frozen positive-$F$ DEV3 transition uses radius $32/25$, expansion
magnitude $225/64$, thickness $1/50$ and endpoint $F=4$. Its coefficients obey

$$
c_{\rm transition}=-\frac{18}{1015},\qquad
v_{\rm room}=-6c_{\rm transition}=\frac{108}{1015},\qquad
\alpha_{\rm GB}/L_*^2=-\frac1{8c_{\rm transition}}=\frac{1015}{144},\qquad
\alpha_{\rm GB}v_{\rm room,physical}=\frac34.
$$

The explicit two-sheet carrier is

$$
ds_5^2=e^{-36y^2/1015}
\left[A(\rho)dt^2-\frac{d\rho^2}{A(\rho)}-(\rho^2+1)d\Omega_2^2-dy^2\right].
$$

Both areal branches are needed: one radius can have distinct parent and
child metric values, joining at the throat. Pure Einstein gravity requires
a negative null source for this carrier; a canonical positive-kinetic scalar
alone cannot provide it. The transition-fixed Gauss–Bonnet coefficient
cancels that local null projection geometrically. On the room, its full
component identity is

$$
[G^A{}_B+\alpha_{\rm GB}H^A{}_B]_{y=0}
=\operatorname{diag}\left(\frac{108}{1015},\frac{108}{1015},
\frac{108}{1015},\frac{108}{1015},
-\frac{R_4+(1015/144)\,GB_4}{2}\right).
$$

Thus the four tangential components have vacuum form and the normal component
contains the room curvature action. The same $\alpha_{\rm GB}=1015/144$ lies
at the five-dimensional EGB critical point: with
$\Lambda_0=-108/1015$, $\alpha_{\rm GB}\Lambda_0=-3/4$ and the vacuum
discriminant vanishes. The archived
[kinetic and interface records](https://doi.org/10.5281/zenodo.23090015),
authenticated by the [result manifest](../results/manifest.json), already
test this actual background. `NSC-1-S-ONE-GAUSS-BONNET-KINETIC-RANK`
reports numerical time-principal rank 5 out of 10 on $y=0$ at six radial
samples, against rank 10 for its GR controls and the tested $y=1$ points.
This is a sampled numerical rank result, not a complete characteristic theorem.
`NSC-1-S-ONE-INDUCED-INTERFACE-RANK` restores time-principal rank 10 by adding
the local induced Einstein block; it does not establish kinetic signs or
fix the relative bulk/interface normalization. The subsequent
`NSC-1-S-ONE-INTERFACE-CHARACTERISTIC-SCAN` finds no tested positive
normalization with real roots in all sampled directions on its 401-point
normalization grid. That scan does not exclude every positive normalization.
The coupled off-diagonal field and full curved spectral response remain
untested there. Reuse these records before designing any new kinetic
calculation. A vacuum-form background coefficient alone does not establish
a physical four-dimensional Lambda. Bulk-source microphysics and full
stability remain open.

The expanding-child Ricci sign selects positive $\Lambda_{\rm effective}$
in the retained doubled-FLRW convention $R=-2\Lambda_{\rm effective}$.
It selects $\Phi^2/\Lambda^2=1006/1015$. The retarded boundary pole then has
positive weights and subluminal group velocity, and its propagating band
ends at $p^2+\Phi^2=\Lambda^2$. The squared bandwidth is $9/1015$, half the
transition-curvature magnitude. This boundary excitation is not an electron
identification or the full graviton propagator.

The imported RS2 black-universe equation uses
$G^\nu{}_\mu=-\Lambda_4\delta^\nu{}_\mu-\kappa_4^2T^\nu{}_\mu
-\kappa_5^4\Pi^\nu{}_\mu-E^\nu{}_\mu$ and
$\Lambda_4=(\Lambda_5+\kappa_5^4\lambda^2/6)/2$.
The outside-kernel variation projects to $-E^\nu{}_\mu$ in that convention.
It matches the black-universe exterior response dictionary; the bulk metric
and its healthy response spectrum still have to determine that projection.

## Spectral and boundary controls

Room geometry forces originate in the regulated fermionic spectrum and its
scale anomaly. The compact allocation separates one massless Dirac field,
two fields at each $m_n=n\pi/\ell$ for $n\geq1$, and the warped complement.
It supplies a Euclidean subtraction, not a real-time state. With
$J=\ell_{\rm proper}/\ell$, the compact light limit is

$$
\lim_{u\downarrow0}[h_\Lambda(u)-E_1(u/\nu^2)]
=2\log\frac{\Lambda J}{\nu}
+2\sum_{n\geq1}E_1\!\left[\left(\frac{n\pi}{\ell J\Lambda}\right)^2\right].
$$

The light singular slope is $1/J^2$ and its canonical mass stays zero.
The imported spectral bulk/boundary constraint is $c_K=2c_R$, leaving one
coefficient for $c_R\int R+c_K\int K$. This universal ratio differs from
the particular transition coefficient $-18/1015$.

For the exponential heat profile,
$h(z)=\int_0^1e^{-\alpha(1-\alpha)z}d\alpha$ and
$F_2(z)=(1+z/4)h(z)-2$. Its nearest complex zeros lie beyond the unit cutoff
disk. Canceling the stationary vacuum background gives
$K_2(z)=F_2(z)+1$, $K_2(0)=0$, $K'_2(0)=1/12$ and a positive massless
propagator residue $12$, with no additional pole in that disk. These flat
and Euclidean controls leave the curved retarded positivity test open.
Reciprocal room symmetry requires $B(\chi)=B(-\chi)$: for
$B=b_0+b_2\chi^2+b_4\chi^4$, the stationary linear scalar mixing is zero,
while $B''(0)=2b_2$ retains the nonlinear transition channel.

The first-order curvature EFT change of inverse metric is

$$
\delta g^{\mu\nu}=-\frac{2c_W}{A_4}
\left[R^{\mu\nu}-\frac{R}{6}g^{\mu\nu}
+\frac{T^{\mu\nu}-Tg^{\mu\nu}/3}{2A_4}\right]
+\frac{c_R}{A_4}\left(R-\frac{T}{2A_4}\right)g^{\mu\nu}.
$$

Here old inverse metric equals new inverse metric plus this shift. Boundary
data, observables and the quantum Jacobian must follow the same perturbative
map; it is not a fundamental fourth-order theory.

## Finite state and energy controls

The coordinate generator is gravity plus $M\operatorname{Tr}(CH)$ with
$M=4\kappa$ counted once. Nodal force partials must be divided by node spacing
before conversion to stress densities. With $N=rL$, $q=rQ$ and
$V=4\pi qr^2$, the unit normal is $N^{-1}(\partial_t-\beta\partial_x)$,
whose dual is $Ndt$. Coordinate energy
$e=L F_L+\beta F_\beta=NV\rho-\beta qVj$ differs from normal-observer energy.
The latter obeys

$$
\partial_t(V\rho)+\partial_x[V(Nj/q-\beta\rho)]
=-NV(p_rK_r+2p_\perp K_\perp)-V(j/q)N_x.
$$

Independent metric variations control this identity and the chart-work
ledger. Mode-projection leakage is extra exchange; the controls establish
neither a continuum limit nor renewal.

The frozen six-mode incoherent imbalance witness uses the Cesàro projection
onto the commutant, with finite-time remainder
$\|C_T-\bar C\|_F\leq2\|Y\|_F/(\gamma T)$. Its retained mean populations are
$(1.3837895948404104,1.0058995585808204,0.6103108465787691)$;
the channel current is $0.0012511140479781876$. Reprojection moves the
response by $0.004920298239642384$, whereas the frozen shortcut gives
$0.00030421076749148856$. These finite controls do not close the physical
local gate. The smooth-motif Weyl-disk tail also controls a projected
resolvent only; it is not identified with a dark-energy density or a cosmic
fraction. Original calculations are at
[Zenodo](https://doi.org/10.5281/zenodo.23090015).

The retained coordinate Hamiltonian and normal-observer stress API lives in
[nested.energy_exchange](../src/recursive_horizons/nested/energy_exchange.py).
It keeps the generator gravity energy plus `M Tr(C H)` with `M=4 kappa`
counted once, the finite regional controls, and the separate normal-observer
shell energy. Renewal and a continuum limit remain unclaimed.
