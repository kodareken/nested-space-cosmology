# Source-fixed evolved incoming state

Douglas's adopted formulation is that the same quantum state is the same
upstream/source preparation, not a covariance pinned to Sigma independently
of the metric. For an allowed history family the incoming kernel is the
evolved restriction of the original source columns:

\[
C_\Sigma[g]=F[g]\,C_{\rm src}\,F[g]^\dagger,
\qquad
\delta C_\Sigma
=\delta F\,C_{\rm src}F^\dagger
+F C_{\rm src}\delta F^\dagger,
\qquad
\delta C_{\rm src}=0.
\]

`F` is the actual rho1 normal/half-density image of those columns. All source
channels, coherences and cross-position kernels are retained. A frozen
incoming `C0` or a legacy matter dictionary is not this state and is rejected
by the new interface. Equality of a genuine evolved kernel with a reference
matrix is not a reason to reject the state.

## Physical entry

[evolve_incoming_state](../src/recursive_horizons/nsc_evolved_incoming_state.py)
calls the owned prepared-field evolution for the supplied
`CompatibleIncomingMetric` family, then restricts the rho1 characteristic
trace with `MODE_TO_CURRENT` and the existing normal/half-density map.
The radius family leaves the intrinsic metric and normal fixed at Sigma, so
the restriction map has zero parameter derivative. An unsupported family that
changes that intrinsic map is rejected; the derivative is not omitted.

`EvolvedIncomingState` exposes `z(nz)`, `columns(nz,2,nsrc)`,
`column_tangents(ndir,nz,2,nsrc)`, `source_covariance(nsrc,nsrc)` and
`column_weights(nsrc)=repeated sqrt(w/(2 pi))`. Properties
`weighted_columns` and `weighted_column_tangents` apply that quadrature once.
`covariance()` and `covariance_tangent()` return the full cross-`z` kernels
with shapes `(nz,2,nz,2)` and `(ndir,nz,2,nz,2)`. Source energies remain
initial labels; they are not substituted as outgoing momenta `k=-E`.
Angular multiplicity is not added. Coverage records sampled-source and
continuum error as OPEN; this is not a boolean proof of a physical history
solution.

Finite arrays are stored read-only and bound to the generating source digest
and metric-family amplitudes, so a changed geometry cannot reuse an arbitrary
old matrix as the physical state.

## Bounded adapter control

The control reuses the authenticated retarded group14 `14_1` 801-node grid,
four signed source energies and compact window, through the existing
`inputs`/`support_guard` helpers. It does not regenerate baseline, scattering,
horizon or source grids. Time coordinates `[0,0.3]` and window `.007/.03` are
the same numerical controls as the retarded pilot, not a physical duration.

The new runs are three 801/64 field solves of the same compact radius family
at amplitude `0.001` and the centered pair `0.001±1e-4`. Initial, incident
and source tangents are explicit zeros by the unchanged past/upstream
support. A chart/radial cache evaluates the actual finite radius at every
time; it is checked against `CompatibleIncomingMetric` at several nodes and
nonzero metric values and is not a zero-radius freeze.

Field, incoming-trace and coherent covariance derivatives are compared with
the centered pair against `3e-8`. Flux algebra is compared against `3e-11`.
Executed residuals live in the
[development record](../results/development/nsc-evolved-incoming-state.json).
Omitting the column tangent or inserting a fake `dC_src` changes the kernel.
An array-only comparison with the saved retarded response checks the same
restriction and `dC_src=0` formula without rerunning that solver.

This is adapter verification, not a physical root. Finite source and
discretization errors remain OPEN. Full CAR, cosmology EXISTENCE, metric
evolution, new stress/Gamma and parameter refits are out of scope.

```sh
python3 scripts/derive_nsc_evolved_incoming_state.py --run
python3 scripts/derive_nsc_evolved_incoming_state.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_evolved_incoming_state.py
```

`--check` replays saved arrays and provenance and does not propagate fields.
The consumer of this state contracts actual `z` derivatives of `F`; it must
not restore a stationary `k=-E` substitution on these columns.


## What an exactly fixed incoming covariance would require

The older fixed-$C_0$ result is a necessary-condition theorem in a different
class. On a finite regular trapped KS slab, take smooth positive $N,a,r$,
compact axial changes, unchanged upstream preparation and canonical frame,
and the same intrinsic $N,\beta,a,r$ on $\Sigma$. Fix the $N/\beta$ normal
jets through order four there. The exact physical state then has

$$C_g(T)-P_g(T):H^s(\mathbb R;\mathbb C^2)
\longrightarrow H^{s+6}(\mathbb R;\mathbb C^2).$$

Here the finite auxiliary projector has transport defect
$R_g=\partial_TP_g+i[H_g,P_g]\in\Psi^{-6}$ for the existing symmetric
first-order operator
$H_g=(A_gD_z+D_zA_g)/2+B_g$, with
$A_g=(N/a)\sigma_3-\beta I$ and
$B_g=-Nm\sigma_1+(N\ell/r)\sigma_2$.
The finite separated-band construction uses the nonzero band difference;
individual shifted characteristic velocities need not be nonzero. The
upstream P16 error is a translation-invariant $O(E^{-16})$ multiplier,
the coherent/thermal correction is regularizing, and the signed relation
$C(-E,\ell)=I-\sigma_3\overline{C(E,-\ell)}\sigma_3$ covers the other
momentum cone. The finite Duhamel identity

$$C_g(T)-P_g(T)=U_g(T,T_u)D_uU_g(T,T_u)^\dagger
-\int_{T_u}^{T}U_g(T,t)R_g(t)U_g(T,t)^\dagger dt$$

and Sobolev propagation give the stated gain. This imports the finite
approximate-projector method of
[Gérard and Stoskopf](https://arxiv.org/html/2108.11630), rather than
assuming a new vacuum or an infinite-series convergence theorem.

Exact operator equality $C_g(\Sigma)=C_0$ forces agreement of classical
endpoint coefficients through degree minus five. For the first unmatched
normal derivative of order $k$, the off-diagonal finite differences are

$$\begin{aligned}
(\Delta p_{k+1})_{01}&=-\frac{a_1^k}{(2i)^{k+1}}
\left[(\ell/r_1-im)\Delta a_{T^k}
-\ell a_1\Delta r_{T^k}/r_1^2\right],\\
(\Delta p_{k+1})_{10}&=\frac{(-1)^ka_1^k}{(2i)^{k+1}}
\left[(\ell/r_1+im)\Delta a_{T^k}
-\ell a_1\Delta r_{T^k}/r_1^2\right].
\end{aligned}$$

For one retained $m\ell\ne0$ block their determinant is
$im\ell a_1^{2k+1}/(2^{2k+1}r_1^2)\ne0$. Induction for $k=1,2,3,4$
therefore fixes all permitted twenty normal/mixed $a/r$ jet differences
to zero. The result is finite-amplitude and necessary, not sufficient;
sampled kernels, finite energy moments or diagonal covariance equality do
not establish its operator hypothesis.

In that exact fixed-$C_0$ class the included shift constraint remains

$$\mathcal E_\beta[g,C_0]=P_K<-M,\qquad
M=\frac{5887634380411363}{265464355000000000}>0.$$

Thus those compact prepared histories cannot solve both included constraints,
even on an open interval. This excludes neither different parent/preparation
data nor the source-fixed evolved state above, and proves no global history
injectivity or extended transmitting stationarity. The existing 33 exact
scalar checks and finite lower-jet controls support the algebra; they are
not numerical continuum-error enclosures.
