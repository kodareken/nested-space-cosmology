# Retarded preparation of a compatible incoming variation

A smooth compact axial variation gives an explicit preparation family for
the same Dirac operator. It need not prescribe a new parent state. Let
`w,U` be smooth real functions supported in `[z_minus,z_plus]`, and let the
normal window have compact radial support K around rho1, strictly inside
the trapped chart. In the fixed chart,

\[
\delta r(\tau,\rho)=\chi\bigl(T(\rho)-T(1)\bigr)
\left[(T(\rho)-T(1))w(\tau+S(\rho))
+\frac{(T(\rho)-T(1))^3}{6}U(\tau+S(\rho))\right].
\]

This is the [implemented compatible geometry](nsc-compatible-prepared-history.md)
with compact axial functions. N, beta and a remain fixed. For amplitudes
in a sufficiently small neighborhood of zero, r stays positive on compact
K. The canonical principal velocities are unchanged. Extend the variation
by zero outside the chart support; choose K away from the physical rho0
cut so both the geometry and variation are unchanged near that cut.

## The preparation condition follows from support

Reuse the [whole-field variation-of-constants formulation](nsc-transmitting-history-modes.md)
and the existing transmitting Dirac domain; do not reconstruct scattering,
horizon modes or a new quantum-state law. The perturbation's support obeys

\[
\rho\in K,\qquad
z_- - S(\rho)\leq\tau\leq z_+ - S(\rho).
\]

Consequently any caller initial slice satisfying

\[
\tau_i<z_- -\max_{\rho\in K}S(\rho)
\]

lies strictly before the entire perturbation. The full preceding geometry
and Dirac operator equal their reference values. The authenticated global
reference columns and original horizon/infinity source covariance therefore
supply the in-state on this slice. Their parameter derivatives vanish
because the entire past family is unchanged, not because an unknown
preparation derivative was set to zero.

The corresponding upper support limit is `z_plus-min_K S`. These bounds
are coordinate inequalities; they do not select a physical elapsed duration.
Changing an initial slice farther into the unchanged past only changes the
reference representation by its known evolution/phase, not the prepared
state. Source-frequency labels remain labels of the reference in-state;
no conserved future energy is presumed on the varying geometry.

Place the right computational edge strictly above `max K` while remaining
in the trapped chart. Both principal velocities `±N/q-beta` point toward
decreasing rho. Thus the perturbation cannot change the upstream incident
field. Supply the same authenticated incident columns, with their actual
reference phase, and explicitly zero incident parameter tangents. Left
outflow is retained. These physical causal statements do not turn a finite
SBP grid into a rigorous continuum approximation; its spatial and temporal
errors need their own numerical evidence.

## The incoming state is a result

The prepared field is the solution of the same Dirac equation with those
reference past/incident data. At rho1 the intrinsic geometry and normal are
fixed, but the propagated field generally changes. Use its actual trace
across PG time, with `z=tau+S(1)`, through the completed fixed-surface map:

\[
F_\Sigma[g]=R_\Sigma\Phi[g]|_{\rho=1},\qquad
\delta F_\Sigma=R_\Sigma\delta\Phi|_{\rho=1}.
\]

For the fixed original source law, `delta C_src=0`, and the incoming kernel
response is

\[
\delta C_\Sigma(z,z')
=\delta F_\Sigma(z)C_{\rm src}F_\Sigma(z')^\dagger
+F_\Sigma(z)C_{\rm src}\delta F_\Sigma(z')^\dagger.
\]

The expression includes the full source measure and retained correlations;
its finite sampled representation must preserve cross-position and
source-frequency coherence. The saved real-frequency control columns can
inspect this response, but do not complete its continuous spectrum. The
existing subgap analytic representation and spectral tails keep their
original meanings; complex contour nodes are not canonical source states.

This construction supplies a causal prepared *off-shell variation family*.
It does not establish `C_Sigma[g]=C0`: that is an additional matching equation
whose residual can now be evaluated using physically justified preparation
inputs. It also does not solve the compatible constraint ODE. In particular,
the fixed-C0 nonzero-flux result precludes treating a compact axial variation
as a global constraint solution merely because it has the six allowed jets.
No new state, source occupation, action term or fitted coupling is introduced.

## Next informative check and stopping condition

Reuse authenticated global reference mode columns that cover both sides of
rho1, their incident traces, real energies, source covariance and measure.
Use an explicitly labeled compact variation with nonzero incoming normal
jets and an initial slice satisfying the strict support bound. Propagate
its field tangent through the existing prepared solver and inspect the
actual rho1 coherent kernel response. Independent resolution or derivative
checks must concern that new response, not rerun the old background source.

A nonzero sampled source-mode or partial spectral response would establish
a preparation response on that sampled component. It would not by itself
prove that the full integrated C0 changes: unsampled contributions and
coherence could matter. Neither a sampled zero nor a sampled nonzero result
settles the full matching equation or excludes other variations. Stop once the direction of this preparation
response is established within the named numerical scope; do not expand
into a full spectrum until the result can change the matching strategy.


The dispatched pilot reuses group14's four signed real source frequencies
and coherent source fibers on the saved801-point grid, with a nested
401-point subsample. Its numerical control uses normal radii0.007/0.03,
axial center `S(1)+0.15`, axial halfwidth0.06, and PG coordinates `[0,0.3]`.
The radial support is approximately `[0.9743860383,1.0250344273]`, and its
full PG-time support approximately `[0.0431750872,0.2562033203]`, strictly
inside that initial/final interval and away from rho0 and the right edge.
These checks select a test function, not a cosmological duration.

The bounded run compares401/801 grids with64 midpoint steps, then801 with128
steps. Its stopping rule is a120-second CPU ceiling or completion of those
three solves. Spatial/time differences below25%/10% of the fine sampled
response are declared diagnostic indicators only. Failure records OPEN and
stops; it does not trigger automatic refinement or a new source spectrum.
