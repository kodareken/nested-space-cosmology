# One finite parent–child pair under the inherited law

This construction uses one global state, metric, and selected action
$\Gamma_{\mathrm{one}}$. The finite child collar $I_C=(1,3)$ is strictly
inside the parent collar $I_P=(0,4)$ of the periodic radial chart. The
angular factor is the same two-sphere. These are nested spatial domains;
they are not separately copied metrics or two adjacent packet labels.
The construction does not assert a regular spherical centre, a horizon
junction, a new spacetime manifold, or an eternal continuation.

The method reuses the
[finite inheritance theorem](../../docs/nsc-nested-qualities.md),
[spherical action](nsc-spherical-feedback-action.md),
[direct source](nsc-conformal-adm-source.md), and
[conformal variation](nsc-spherical-conformal-gauge.md).
The earlier [graded compression](nsc-geometry-graded-window.md) is an
initial operator measurement. Its large complement leakage is retained;
that six-mode compression is not substituted for full field evolution.

## Nested domains, modes, and geometries

The new rank-six successor source retains the middle packet pair on $I_C$
and places two exterior source pairs on $(0,1)$ and $(3,4)$.
It reuses the canonical bump/lobe utilities and phase convention without
editing the old preparation. The observer hierarchy is separate: the
original middle pair defines $V_C$, and the original overlapping six
observer columns define $V_P$, with

$$
P_C=V_CV_C^\dagger,\qquad P_P=V_PV_P^\dagger,
\qquad P_CP_P=P_C.
$$

The orthogonal modal partition is
$P_C$, $P_D=P_P-P_C$, and $P_E=I-P_P$. The last projector retains the
ambient complement, including every leakage channel. The spatial and
modal partitions are distinct. The prepared exterior source profiles are
spatially outside the child collar; their finite-band interpolation tails
and exact child-frame orthogonality are measured separately. The original
parent-detail observer modes overlap the child collar, while ambient modes
can also have support there. An exterior omission names its actual projector
and preparation rather than identifying the entire Hilbert complement with
the spatial annulus.

For the local continuum Dirac expression, smooth disjoint source supports
give zero direct source-frame compression link at the initial slice.
Numerical interpolation tails do not establish a physical throat coupling.
The live nonzero observer-hierarchy link is instead measured on the
overlapping fixed observer modes. Physical parent-to-child and reciprocal
consequences are checked through the separated-source interventions and
regional flux/geometry increments. The modal response controls retain their
own interpretation. A later reduction starts with the actual transported
state; initial cross correlations are measured there rather than inserted.

New fixed geometry basis functions separate a broad parent component,
localized child details, and any required exterior components. All their
canonical coordinates evolve jointly after varying the projected action.
They reconstruct a single full $Q,r,\chi$ and their momenta. The physical
parent and child metrics are restrictions of this full metric, with the
coordinate pullbacks below. A broad coefficient field or a detail field
alone is not a physical metric. Positivity is checked on the reconstructed
quadrature metric.

For a basis matrix $W$ and quadrature weights $w$, the symplectic mass
matrix is $G=W^T\operatorname{diag}(w)W$. If coordinates and momentum
density use the same basis, the canonical coefficient momentum is
$\pi=Gp$. Equivalently, a weighted orthonormal or dual basis may absorb
$G$. Rates and source forces must use this same normalization. Independent
parent and child variations require a full-rank basis; a post-step filter
or a prescribed child trajectory does not supply those variations.

The implementation uses a complete real orthogonal nodal $W$. Its actual
canonical map is $a=W^Tg$, $\pi=\Delta x_g W^Tp$, with inverse
$g=Wa$, $p=W\pi/\Delta x_g$. This is an exact coordinate change of the
full existing Galerkin phase space. The basis alone is not a new physical
observation. The separated source, its newly prepared and evolved metric,
and the reciprocal regional responses provide the new pair evidence.

## Inherited functional law and the two scales

On either collar write

$$
x=a+\ell\xi,\qquad t=\ell\theta,
\qquad \widehat Q=\ell Q,\quad \widehat L=\ell L,
\quad \widehat\beta=\beta,\quad \widehat r=r,
\quad \widehat\chi=\chi.
$$

Here all right-hand fields are composed with the chart map. The parent
has $(a,\ell)=(0,4)$ and the child $(1,2)$, giving coordinate ratio two.
The metric retains the same template

$$
ds^2=\widehat r^{\,2}
\left[\widehat L^{\,2}d\theta^2-
\widehat Q^{\,2}(d\xi+\widehat\beta d\theta)^2-d\Omega_2^2\right].
$$

The canonical spinor pullback is
$\widehat\varphi=\sqrt\ell\,\varphi$; it preserves its radial norm.
With $P_\xi=-i\partial_\xi$, the same Dirac template becomes

$$
\widehat H=\sigma_2\tfrac12\{\widehat L/\widehat Q,P_\xi\}
 +\kappa\widehat L\sigma_1-\tfrac12\{\widehat\beta,P_\xi\}
=\ell\,U H U^\dagger.
$$

This identity uses the pulled-back operator domain and its exterior
couplings. Closing either collar with a new independent boundary condition
would change that operator. The global field and measured links retain those
couplings in this construction.

The first-order geometric density transforms as
$\widehat{\mathcal L}_g=\ell^2\mathcal L_g$, exactly matching
$dt\,dx=\ell^2d\theta\,d\xi$. Its coefficients and local functional law
are unchanged. The momenta transform as
$\widehat p_Q=p_Q$, $\widehat p_r=\ell p_r$,
$\widehat p_\chi=\ell p_\chi$; the lapse and shift constraint densities
transform as $\widehat C=\ell C$, $\widehat D=\ell^2D$.
This is a coordinate and clock pullback of the same action. The normal
proper time also agrees, since
$r\widehat L\,d\theta=rL\,dt$.

Unequal evolving normalized profiles need not have identical instantaneous
spectra. The exact fixed-window spectral theorem applies to homologous
model data, as its hypotheses state. Here the inherited statement concerns
the shared functional law; measured dressed operators and spectra retain
their actual state and geometry dependence.

A physical scale ratio is measured separately, for example from the
declared radial proper lengths

$$
\mathscr L_P(t)=\int_{I_P}rQ\,dx,\qquad
\mathscr L_C(t)=\int_{I_C}rQ\,dx,\qquad
\Omega_{\mathrm{proper}}(t)=\mathscr L_P(t)/\mathscr L_C(t).
$$

Normal clocks belong to specified worldlines, with
$d\tau_a=r(t,x_a)L(t,x_a)dt$. Coordinate ratio two, proper-length ratio,
and lapse/clock ratio answer different questions. An independently imposed
$\Omega(t)$ or lapse path would be an external drive. If $\Omega$ is only
a readout of the full metric, it adds no extra force; if a Hamiltonian is
parameterized through it, all its geometry chain derivatives are included
in the action variation.

## Source-consistent initial preparation

The successor initializer keeps the owned projected lapse residual and
Jacobian, with constant initial $Q$, zero auxiliary field and momenta, and
the actual separated Gaussian source. A positive source homotopy solves
for the radius. The source density $F_L/\Delta x$ is independent of radius
in this canonical chart; it is not guessed from a desired radius.

Both line search and stopping use the Newton correction in radius units.
The stop requires the ordinary BLAS-ordered and independently reordered
residual corrections to be no larger than eight radius ULP plus the
residual-ordering difference mapped through the actual Jacobian. This is
a finite arithmetic indicator, not a certified bound. The report retains
the raw projected, full and held-out constraints, Jacobian norm/condition,
linear backward error, positivity bracket, and the mean shift source.
No current is deleted and no historical absolute residual line is applied
as the physical acceptance criterion.

The canonical $W$ transform can change last bits. The reported physical
constraints and correction are therefore evaluated again on the actual
returned reconstructed state, separately from the solver diagnostics.
Space refinement and the measured response scale remain distinct from
this initial arithmetic stationarity statement.

## Reciprocal feedback and energy

The global generator has the form

$$
H_{\mathrm{tot}}=H_g+M\operatorname{Tr}(C H_f(g)),
\qquad M=4\kappa,
\qquad \dot\Phi=-iH_f(g)\Phi.
$$

Its fixed weights and full covariance supply every geometric source once.
The representative field generator has no extra $M$. Each stage recomputes
the geometry-derived links, such as $B_{CD}=V_C^\dagger H_f(g)V_D$.
Link energy and its geometric derivative are included in the same full
trace. Child-to-ambient links are retained when nonzero; a nearest-neighbor
continued fraction is not imposed on a general measured operator.

The finite autonomous Hamiltonian supplies both directions of feedback:
the metric changes the field operator, while state changes alter the
parent and child canonical momentum forces. In the conformal gauge the
field work contains $F_Q\dot Q+F_L\dot L$, with $\dot L=\dot Q$.
Source-force density conversion and the geometry adjoint pullback retain
their quadrature spacing; the nodal work sum has no further spacing or
copy factor. Independent directional derivatives of the full energy
check these reciprocal variations.

For a fixed projector $P_A$, a consistent modal allocation is
$h_A=\{P_A,H_f\}/2$. Then

$$
E_A=M\operatorname{Tr}(C h_A),\qquad
\dot E_A=M\operatorname{Tr}(C\,i[H_f,h_A])
 +M\operatorname{Tr}(C\dot h_A).
$$

Allocate only the disjoint child, parent-detail, and ambient sectors in a
total sum. The parent aggregate already includes its child; adding that
aggregate and the child again double-counts energy. Bare diagonal block
energies instead require a separate complete link-energy account.

Spatial shell energy has its own normal-observer ledger: both child
surfaces $x=1,3$, pressure work, lapse-gradient work, and the measured finite
residual. Parent outer cuts and the periodic seam are retained in the full
account. Internal cut contributions cancel on assembly because the metric
and field are shared; no new interface action is inserted. The previous
periodic cap completion remains the global boundary domain. A moving
coordinate cut or a time-dependent basis would require its additional
transport or connection terms; this construction uses fixed ones.

## Finite completion checks

The independent review checks the following concrete conditions before
reporting the pair as measured.

- Strict spatial and modal containment, positive reconstructed metric, two
  distinct measured regional scales, and continuous declared clocks.
- Full-rank independent parent and child canonical geometry groups with
  the correct symplectic mass matrix; both groups evolve under the varied
  global action rather than a supplied trajectory.
- The affine action and Dirac pullbacks, with spinor normalization and
  frequency conversion; unequal-state spectral identities are not assumed.
- Live geometry-derived coupling and reciprocal full-trace forces,
  including cross/link terms, multiplicity once, and all ambient channels.
- Admissible covariance and Gram, reported initial/full constraints and
  actual finite forcing, autonomous total-energy closure, and separate
  modal and physical shell accounting.
- A quantitative parent-to-child local response and reciprocal source or
  geometry response. Perturbation controls state exactly what initial
  regional data they preserve and whether their exterior is spatial or
  modal.
- Full and reduced propagation use identical preparation, full initial
  cross correlations, fixed child observer, metric schedule, and clock.
  Memory, exterior drive, and actual cross omission each retain their own
  comparison discrepancy.
- Space, timestep, and applicable frame/operator movements are compared
  with each claimed effect using the accepted one-percent criterion.
  Unresolved tiny effects are reported as unresolved; refinement differences
  are indicators unless an independent enclosure is supplied.

These conditions describe a finite parent–child realization. They do not
require a cosmological fit, the paused incoming gate, exact continuum
constraint propagation of the discrete model, or a universal GR certificate.

The parent-source control changes the two exterior pair weights oppositely,
keeping total occupation three and the child covariance unchanged at the
initial slice. The child-source control changes only the middle pair and
declares its changed initial total occupation; weights remain fixed during
each subsequent evolution. Initial constraints are prepared from each
actual source. Reciprocal parent consequences use the physical annulus
$I_P\setminus I_C$ and parent-detail occupation, rather than the parent
aggregate that already contains the child. Where source-dependent initial
constraint solves change a metric, comparisons of subsequent increments
separate that initial difference from the dynamic response.

The primary state effects are declared before trajectories: canonical
probability in the physical child collar and in the parent annulus.
The primary geometry effects use their proper-length-weighted mean radius.
Fixed child and parent-detail modal occupations are additional observables;
they do not replace those physical regional effects. Source-control
increments are compared at the same global coordinate time, with each
continuous proper clock recorded. Equal proper-time endpoints are not
claimed for these comparisons.

## Recorded pilot v1: measured pair, incomplete precision and reduction

The [v1 record](../results/development/nsc-nested-parent-child-v1.json) and
[payload](../results/development/nsc-nested-parent-child-v1.npz) preserve
twelve trajectories: baseline, exterior-parent intervention, and
child-source intervention at $n_f=128,256$ and step caps $0.001,0.0005$.
All reach $T=0.3$. CPU is $401.36$ s and payload is about $15.85$ MB.
The record status is `INCOMPLETE`; reaching the target does not by itself
complete the finite pair objective.

The primary effects compare source-control increments from each initial
slice against baseline increments at common coordinate time. The fine
half-step measurements and their spatial indicators are:

| Primary physical response | Effect | Spatial movement / effect |
|---|---:|---:|
| Parent source to child probability | $0.01150685$ | $1.53375\%$ — unresolved |
| Child source to parent-annulus probability | $0.00848896$ | $0.02622\%$ |
| Parent source to child mean radius | $0.0001759865$ | $0.22315\%$ |
| Child source to parent-annulus mean radius | $0.001109926$ | $0.01600\%$ |

Timestep movements are below one percent for all four. The additional
parent-to-child modal occupation effect is $6.5150\times10^{-6}$, with
$3.62\%$ spatial movement, and is separately unresolved. These are
observed refinement indicators, not total continuum state-error bounds.

The saved initial parent intervention preserves the child covariance
exactly while redistributing the exterior source weights at total
occupation three. The child intervention changes its initial covariance
by $0.1I_2$ and declares total occupation $3.2$. Fine exterior source
support leakage is about $7.16\times10^{-8}$ in power. Geometry stays
positive and covariance remains admissible; fine half-step column Gram gaps stay
below $9.6\times10^{-12}$. Both canonical geometry groups evolve.

On the fine baseline, the measured proper-length ratio changes from
$1.9842168$ to $1.9827345$. Child mean radius increases while the
parent-annulus mean decreases. The full sampled constraint norm changes
from about $0.008999$ to $0.009003$, within its sampled forcing budget;
the full initial Hamilton maximum is about $0.08014$. These quantities
retain the finite projected/held-out and arithmetic domains of the
initializer. They are not silently replaced by its solved projected
residual or by an exact continuum constraint claim.

The child normal shell increases by $0.94679485$. Boundary inflow
$0.92339858$, pressure work $0.02083885$, and lapse-gradient work
$0.00254640$ explain it to $1.10\times10^{-5}$, or $0.001164\%$;
the frame quadrature indicator is $0.0720\%$ of the shell change.
Parent/child/annulus probability allocation closes to
$1.33\times10^{-15}$ and complete field cross/link accounting to
$2.84\times10^{-14}$. Coordinate field exchange is
$1.67652\times10^{-5}$; total-energy drift is
$-5.34\times10^{-10}$, with full metric work including the dynamic lapse.

The attempted local reduction exceeded its admitted CPU ceiling and
returned no local result. Its intended start state at $T=0.1$ is saved,
with actual $\|C_{AE}\|_F=0.3478745$; the corresponding initial $T=0$
cross is roundoff zero. The remaining evidence is therefore specific:
resolve the parent-to-child physical state effect under finer space
comparison, and complete the same-trajectory local reduction with its
individual control discrepancies. Reusing the exact linear-$Q$ conformal
operator schedule can remove repeated dense operator construction without
changing that physical schedule. The sealed v1 states and source remain
the comparison domain for the confirmation.

## Confirmation v2: the declared finite pair is complete

The [confirmation](../results/development/nsc-nested-parent-child-confirmation-v2.json)
and its [payload](../results/development/nsc-nested-parent-child-confirmation-v2.npz)
record `MEASURED_CONFIRMED_NESTED_PAIR`. They preserve the v1 witness and
add only the two $n_f=512$ baseline/parent-source trajectories needed for
the finer spatial comparison, plus the local response on the saved witness.
The confirmation uses $559.82$ CPU seconds; the two campaigns together use
about $961.19$ CPU seconds. No source reset or imposed geometry is introduced.

The parent-to-child probability response is $0.01150720$, with
$0.01348\%$ spatial movement at $256\to512$. The same comparison resolves
its modal response at $0.02578\%$ and improves the child mean-radius
indicator to $0.00033864\%$. The three other primary v1 physical effects
already meet the one-percent target and retain their original comparisons.
All four physical responses therefore resolve within the declared finite
numerical policy. Comparisons use initial-subtracted increments at common
coordinate time; the recorded proper clocks are not matched endpoints.

On the same generated geometry over $[0.1,0.3]$, the fixed child's
occupation changes by $0.294235$. The streamed/full discrepancy is
$0.04779\%$ of that change. Exterior-drive and exterior-cross omissions
change occupation by $0.136885$ and $0.176829$; parent-detail-drive and
child/detail-cross omissions change it by $0.0331275$ and $0.0509377$.
Each omission's own comparison discrepancy is below $0.203\%$ of its
effect. The cross-pair omission is an algebraic diagnostic, while pinching
all child/exterior cross correlations is an admissible covariance control.
The effects are not additive. The exact live Schur reduction retains the
ambient region and its causal response; no memory-off series is claimed
for this successor.

Independent saved-array checks reproduce the physical curves, local
covariances, control discrepancies, source invariances, energy accounts
and direct Schur responses. The
[portable basis](../results/development/nsc-nested-parent-child-replay-basis-v1.json)
freezes the actual canonical frames and source/observer columns at all
three bands. It reconstructs all 434 saved geometry frames with zero
metric discrepancy on the recorded host, including clocks and endpoint
energies. Its [producer](../scripts/derive_nsc_nested_parent_child_replay_basis.py)
does not re-evolve any trajectory. The
[independent consumer](../scripts/check_nsc_nested_parent_child_confirmation.py)
also rejects mismatched coordinate clocks across comparison bands.

Read-only completion replay from the repository root:

```sh
python scripts/lab.py scripts/check_nsc_nested_parent_child_confirmation.py
python scripts/lab.py scripts/derive_nsc_nested_parent_child_confirmation.py --check
python scripts/lab.py scripts/plot_nsc_nested_parent_child.py --check
```

The [four-panel figure](../results/development/nsc-nested-parent-child-figure.png)
and its [data manifest](../results/development/nsc-nested-parent-child-figure.json)
replay the actual source, reciprocal physical responses, length ratio and
proper clocks. The result is one finite nested pair under one evolving
metric and the common inherited functional law. Its refinement indicators
do not become continuum state-error certificates. The existing manuscripts
and release remain unchanged.
