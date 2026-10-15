# Reconstructing the central scalar response

The [scalar metric response](nsc-atlas-scalar-metric-response.md) derives
the forced exterior equation and its transmission into the cap. The
[exterior evolution](nsc-atlas-scalar-evolution.md) supplies the areal-gauge
wall coefficient and its total material derivative. This note completes
the analytic connection to the central Jordan scalar, using the same
leading field, quadratic metric, source continuation and finite domain.
It introduces no new incident wave or physical amplitude. Its numerical
adapter has manufactured controls; a computed central coefficient requires
separate admitted numerical evidence.

All quantities are dimensionless and use the established expansion
`alpha = kappa_4 epsilon^2 / ell^2`, with `h=1`. The cap time is Jordan
central proper time. Let

\[
u=T-\xi,\quad q=T+\xi,\quad \xi=\operatorname{atanh}r,
\qquad c_0(u,q)=f(q)-f(u).
\]

The central-clock null coefficient obeys

\[
\widetilde c_{1,uq}=S(u,q),\qquad
S=-K_1c_0,\quad K_1=\frac{A_0W_1}{4},\quad A_0=1-r^2,
\]
\[
W_1=\frac{2\nu_E}{r^3}-\frac{p^2}{3}
-\frac{p_T^2}{6A_0}+\frac{A_0p_r^2}{6}.
\]

The [cap metric](nsc-atlas-cap-response.md) supplies the mass and full
reflected scalar jets. The pre-arrival mass moment and quiet-profile
potential used on an outgoing wall ray do not apply along the whole
incoming ray to the centre.

## Incoming transport and regular reflection

Parametrize the reference wall by its incoming label `q`. Write
\(a(q)=q-2\operatorname{atanh}x(q)\), and define

\[
H(q)=\widetilde c_{1,W}
=y_W(V(q))-C(q)f'(q),\qquad C=\Delta C_{\mathrm{fixed}\,V}.
\]

The previously computed outgoing derivative is
\(J(q)=\widetilde c_{1,u}(a(q),q)\). Differentiating the wall trace gives
\(\widetilde c_{1,q}|_W=H'-a'J\). Integrating the forced equation along
the incoming characteristic therefore gives

\[
I(u,q):=\widetilde c_{1,q}(u,q)
=H'(q)-a'(q)J(q)+\int_{a(q)}^u S(z,q)\,dz.
\]

Regular reflection imposes \(\widetilde c_1(T,T)=0\), hence
\(\widetilde c_{1,u}+\widetilde c_{1,q}=0\) on the centre. The central
Jordan coefficient is consequently

\[
\boxed{p_{1,C}(T)=2H'(T)-2a'(T)J(T)
-2\int_{a(T)}^T K_1(u,T)[f(T)-f(u)]\,du.}
\]

No second cap evolution is necessary for this central readout. The factor
two is the regular reflection factor, not a convention for the incoming
amplitude. With \(u=T-2\operatorname{atanh}r\), the equivalent radius rule is

\[
\boxed{p_{1,C}(T)=2H'(T)-2a'(T)J(T)
-\int_0^{x(T)}W_1(T-\operatorname{atanh}r,r)
 [f(T)-f(T-2\operatorname{atanh}r)]\,dr.}
\]

Here \(x(T)\) means the radius of the reference wall event whose incoming
label is `T`; it is not the radius of the wall at cap static time `T`.
The integrand uses \(c_0=rp\) for stable evaluation. Regularity gives
\(K_1=O(r^2)\), \(c_0=O(r)\), so this forcing integrand is \(O(r^3)\).

For a full reduced-field reconstruction, the same first derivative gives

\[
\widetilde c_1(u,q)=\int_{\max(u,q_a)}^q I(u,v)\,dv,
\]

with the certified quiet extension below `q_a`. Monotonicity of the wall
outgoing label keeps each inner segment inside the retained cap. This
formula is an analytic identity; the current adapter evaluates the central
readout only. No additional homogeneous function or centre mass constant
remains once the quiet data, wall trace and regular reflection are fixed.

## Arrival pullback and its derivative

Keep material wall time `s` and put
\(Q=q_{+,s}=(\beta_-+v)/A\),
\(Q_-=q_{-,s}=(\beta_--v)/A\),
\(V_s=V_{0,s}\), \(A=1-x^2\). Then

\[
a'=Q_-/Q,\qquad
H'=\frac{y_{W,s}}Q-\frac{C_s}Qf'-Cf''.
\]

The exterior solver returns \(dy_W/dV\); the material derivative supplied
to this reconstruction is \(y_{W,s}=V_s\,dy_W/dV\). It already includes
the wall transport term. The cap profile's second derivative comes from
the differentiated continuum leading-field equation, not sampled values.

Define \(D=\delta_E-2\nu_E/(rA_0)\). The incoming null shift and arrival
coefficient at fixed material time are

\[
\lambda_+=-\int_0^x\frac{D(q-\operatorname{atanh}r,r)}{A_0(r)}\,dr,
\qquad \Delta=T_1+e/A+\lambda_+.
\]

The derivative needed here is a first-jet integral:

\[
D_T(T,r)=-\frac{p_0p_{0,T}}6+
\int_0^r\frac{(1-3z^2)p_Tp_r-2zpp_T}{1-z^2}\,dz.
\]

It follows by differentiating the optical radial constraint and substituting
\(\nu_{E,T}=r^2A_0p_Tp_r/2\). The central anchor is
\(p_0=2f'(T)\), \(p_{0,T}=2f''(T)\). With
\(I_T=\int_0^x D_T(q-\operatorname{atanh}r,r)/A_0(r)\,dr\),

\[
\Delta_s=T_{1,s}+\frac{e_s}A+\frac{2xve}{A^2}
-\frac{vD_W}A-QI_T,
\]
\[
C=\Delta-\frac{Q}{V_s}V_1,\qquad
C_s=\Delta_s-\left(\frac{Q}{V_s}\right)_sV_1
-\frac{Q}{V_s}V_{1,s}.
\]

All remaining derivatives are reference coefficients:

\[
Q_s=\frac{\beta_{-,x}v+v_s}A+\frac{2xv(\beta_-+v)}{A^2},
\quad d=v-\beta_+,\quad d_s=v_s-\beta_{+,x}v,
\quad V_{ss}=V_s\left(\frac{V_s}V-\frac{d_s}d\right).
\]

These expressions avoid differentiating numerical cap clocks. Their
continuum slopes need not be exact derivatives of the finite clock
interpolants. Cap clock quadrature and the leading profile therefore retain
separate numerical consistency checks, just as in the exterior lift.
The arrival spline derivative is reported separately from `Q/V_s`; their
discrepancy is not silently removed.

The terms \(-2C_sf'/Q\) and \(-2Cf''\) include both arrival-rate and
arrival-shift effects at fixed central `T`. Adding another shifted-leading-
waveform term would count the same pullback twice. At the centre,

\[
\partial_r c_{1,\mathrm{static}}=p_{1,C}+p_0^3/12,
\qquad
p_{1,J}=\partial_r c_{1,\mathrm{static}}-p_0^3/12=p_{1,C}.
\]

The cubic frame term thus cancels the central lapse pullback. The outgoing
geometric shift `lambda_minus` is unnecessary for this central observable.

## Domain, reuse and numerical interface

For each `q` in the supplied arrival interval, the incoming segment is
\((T,r)=(q-\operatorname{atanh}r,r)\), `0 <= r <= x`. Its incoming label
never exceeds `q_b`, and its outgoing label lies in `[qmin,q]`. Before the
last wall event, the inward ray stays inside the expanding wall; afterward
it stays inside the terminal inward cone. Constant-time radial metric
integrals remain in that same cap. Both reflected-profile edges and
central-anchor edges split the ray quadrature.

Fixed-central-clock physical statements require a compact subinterval of
`(q_a,q_b)`. Endpoint coefficients are one-sided formal limits. They are not
a continued-domain or finite-amplitude assertion, and the prepared tail is
not set to zero at `q_b`.

The [central adapter](../src/recursive_horizons/nsc_atlas_central_scalar.py)
provides:

- `optical_time_derivative`: paired first-jet radial integrals and their
  central anchor, with at most 513 query endpoints.
- `incoming_response_integrals`: signed forcing, optical time and null
  travel integrals, with sampled defects and finite-domain margins.
- `reconstruct_central_scalar`: algebra from supplied compatible cap/wall
  rows, profile derivatives, incoming integrals and exterior wall response.
- `CentralScalarResponse`: a live adapter to the existing cap and wall
  solutions, retaining at most 129 geometric samples and their shared budget.

The lower-level reconstruction permits explicit authenticated reuse of
existing cap `T1`, `T1_s`, `e`, `e_s`, wall metric and ray data. Changing only
the exterior clock quadrature leaves these cap quantities unchanged, but
requires recomputing `C` and `C_s` with the current `V1` and `V1_s`. The
caller must authenticate source, settings and coordinate compatibility;
the algebraic helper does not authenticate files or infer that compatibility.
The freshly calculated null-travel integral is compared with a supplied
saved value when available.

The [tests](../tests/test_nsc_atlas_central_scalar.py) use analytic profiles,
an independent optical derivative, a closed static mass, and forced null
polynomials to check the integral sign, reflection factor, moving-wall term,
clock derivative, frame cancellation, zero source and cubic scaling. They
do not propagate the physical prepared field or establish a uniform error
bound. Numerical arrival and forcing refinements remain required before
reporting a central waveform coefficient.
