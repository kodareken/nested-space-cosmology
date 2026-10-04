# Coupled PG state for classical RN

This note owns the coupled flat-slice geometry and its numerical state.
The observer, the reference sample and the source preparation are
`nsc_rn_observables`, `nsc_rn_reference` and `nsc_rn_source`.
The state module calls `nsc_rn_source.shell_metric_densities` when that
hook exists. Until then the Dirac variation below is evaluated here, once.

There is no Weyl pole and no cosmological term. Charged lowest-Landau
level and electric Gauss are later. A parallel MMP run does not block
this classical sector. No production trajectory is advanced here.

## Action before the areal gauge

Signature `+---`. The pre-gauge chart is

\[
ds^2=N^2dt^2-q^2(dx+\beta dt)^2-r^2d\Omega^2.
\]

The bulk action, with \(V_4=0\), is the spherical reduction already fixed
by `nsc_spherical_action`:

\[
S=\int\sqrt{|g|}\bigl[-A_4 R_L-C_F F^2\bigr]+S_{\mathrm{Dirac}},
\qquad G=\frac1{16\pi A_4}.
\]

The executable monopole is the integer \(q=1\). It is not solved from
\(r_m\). \(A\) is the audited Einstein coefficient, so

\[
G_N=\frac1{16\pi A}
\]

is not set to 1. \(C_F\) is the audited Maxwell coefficient rescaled by
\((\mathrm{flux}_{\mathrm{audited}}/1)^2\). Classical \(P^2=C_F q^2/(4A)\)
therefore matches the flux-4 geometry. This integer-\(q=1\) chart is a
rescaled coupling. It is not the audited flux-4 quantized charge sector.

\[
r_m^2=\frac{C_F q^2}{4A}
\]

is an output. A caller-supplied \(r_m\) that disagrees is refused.

Neutral Dirac has no electric current, so this monopole solves Maxwell.
The canonical half-density is \(u=r\sqrt{q}\,\psi\), held fixed when the
metric is varied. One \(\kappa=1\) radial block has Hamiltonian

\[
H=\sigma_2\tfrac12\{N/q,P_x\}+\sigma_1 N\kappa/r-\tfrac12\{\beta,P_x\},
\qquad P_x=-i\partial_x,
\]

and the closed shell multiplies that block by \(W=4\) once. Rank 4 with
another factor 4 is refused. \(\kappa/r\) is the angular coupling, not a
4-dimensional mass. In the characteristic frame,
\(\partial_t\chi_+=-a_+D\chi_+-\alpha\chi_-\) and
\(\partial_t\chi_-=-a_-D\chi_-+\alpha\chi_+\), with \(\alpha=\kappa N/r\).
The ingoing angular term is \(+\alpha\chi_+\). A supplied frequency is
metadata; it does not by itself prove an exterior Killing mode.

On the areal chart \(q_{\mathrm{ADM}}=1\), the angular-integrated forces
\(F_N\), \(F_\beta\) and \(F_r/N\) do not depend on \(N\) or \(\beta\). The
state caches them at \(N=1\), \(\beta=0\). That cache is the intrinsic
stress. The Dirac operator uses the reconstructed \(N\) and \(\beta\). The
cache is not iterated.

Physical \(\rho=F_N/(4\pi r^2)\). The angular-integrated density is \(4\pi\)
times that value and is not substituted for \(\rho\). The shell factor 4 is
already inside \(F_N\), \(F_\beta\) and \(F_r\). The Einstein equation does
not multiply by 4 again.

The physical radial-pressure integral is

\[
F_N+\frac rN F_r^{\mathrm{actual}}=F_N+r F_r\big|_{N=1}.
\]

## Gauge and the two constraints

Flat spatial slices are \(q=1\) and \(x=r\). On that slice the shift that
preserves the definition of the charge-corrected mass is

\[
\beta=N s,\qquad s=\sqrt{2m_Q/r-Q^2/r^2},
\qquad m_{\mathrm{MS}}=r\beta^2/(2N^2),\qquad m_Q=m_{\mathrm{MS}}+Q^2/(2r).
\]

With \(v=\beta/N\) the Hamiltonian and momentum equations reduce to

\[
\partial_r m_Q=4\pi G r^2(\rho-vj)=G(F_N+v F_\beta),
\]

\[
\partial_r\log N=4\pi G r j/v=-G F_\beta/(r v).
\]

Integrate \(m_Q\) from the inner stock and \(\log N\) inward from the single
normalization \(N(r_{\mathrm{out}})=1\). Those are the only boundary data.
The inner mass is not also fixed as a Dirichlet value, and \(N\) is not held
at 1 once \(F_\beta\neq 0\). Holding \(N(r_{\mathrm{out}})=1\) fixes
\(\partial_t N(r_{\mathrm{out}})=0\). That is a boundary constraint, not a
reason to freeze \(N\) in the interior, and \(N\) is not free to vanish.

Source-free data \(\varepsilon_N=\varepsilon_\beta=0\) give constant
\(m_Q=M\) and \(N\equiv 1\), hence

\[
\beta=\sqrt{2M/r-Q^2/r^2}.
\]

The same geometric Euler–Lagrange operator vanishes on that
Painlevé–Gullstrand slice and on the Schwarzschild–Droste chart
\(N=\sqrt{f}\), \(q=1/\sqrt{f}\), \(\beta=0\), with
\(f=1-2M/r+Q^2/r^2\). The observer's 4-dimensional Ricci scalar is 0 on
the source-free jets, and \(m_Q=M\).

## Inner-mass flux

The local mass flux, with both copies of \(v F_N\), is

\[
\partial_t m_Q=G N\bigl[2v F_N+(1+v^2)F_\beta+v(r/N)F_r^{\mathrm{actual}}\bigr].
\]

The second \(v F_N\) is required. The unit-lapse cache may replace
\((r/N)F_r^{\mathrm{actual}}\) by \(r F_r|_{N=1}\).

The inner stock uses this value at excision. The profile of \(m_Q(r)\) is
rebuilt from the constraints at every RK4 stage; nothing is restored by a
fitted force. The difference between that rebuilt \(\partial_t m_Q\) and the
pointwise \(q\)-equation is reported as the balance residual.

## Characteristics, domain and clocks

Stored \(\phi\) has shape `(2, points, rank)` in the \(\sigma_2\) eigenbasis.
Component 0 is outgoing and component 1 is ingoing:

\[
\lambda_\pm=-\beta\pm N.
\]

The radial derivative is the second-order diagonal-norm SBP operator.
Excision is the midpoint of \((r_-,r_+)\). Both speeds are outward there,
and no component is set to zero. The outer boundary is nonperiodic. Only
the incoming characteristic receives an absorbing SAT, with penalty speed
\(\lambda_-\). Default \(r_{\mathrm{out}}=32 r_m\) with \(r_m=\sqrt{P^2}\), not \(r_-\). The
child collar and the parent exterior are disjoint areal intervals on this
one chart. Metric forces come from `nsc_rn_source.source_evaluate`, which
applies the multiplicity 4 once. Occupations are the fixed Gaussian
weights in \([0,1)\), not the weighted Gram.

A stage is admitted only when \(N>0\), \(s^2>0\), excision is pure outflow
and the outer end has exactly one incoming family. The step obeys

\[
\Delta t\le \Delta r\big/\max(|N|+|\beta|).
\]

Clocks are coordinate time, the outer static Killing clock
\(\sqrt{N^2-\beta^2}\) where \(\partial_t\) is timelike, and the excision
PG normal clock \(N\). Mass-flux stocks accumulate the excision and outer
\(q\)-equation rates. The SAT debit is a separate stock. Probability
remaining and the boundary probability outflows are a third stock. They are
not mass flux and they are not the SAT debit.

`metric_jets` is spatial only: one-sided polynomial stencils at the ends and
a centered stencil in the interior. Provenance is
`one_sided_eighth_and_centered_fourth_nodal_derivative`. An observer must not
read a dynamical \(R_4\) from those arrays. Time derivatives of \(N\),
\(\beta\) and \(m_Q\) come from the characteristic rays and from the
linearised constraints. A finite central difference of one stage checks that
linearisation. Source-free \(R_4=0\) is the reference curvature on the
stationary chart, not a time jet.

The production source state is one packet column and a filling \(\nu\).
`make_state` accepts that weight vector or a Hermitian matrix of weighted
columns. The stored covariance is the physical CAR
\(G^{1/2}CG^{1/2}\) when the Gram is positive definite. Multiplicity does
not enter it, and it is not divided by rank. A supplied positive frequency
is metadata. It does not prove an exterior Killing mode. The reference
owner does not provide `require_exterior_killing_frequency`.

No production trajectory and no scientific record are written here. The
checks are numerical identities on this chart.
