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
