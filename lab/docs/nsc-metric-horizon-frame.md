# Full-profile horizon frame and representative reflection enclosure

This source-method increment closes the missing **local horizon-frame series**
and maps an already enclosed exterior decaying phase to its asymptotic horizon
reflection. It does not yet certify the stored source at rho=1. No source,
occupation, covariance, history or physical parameter is replaced.

Owner: `src/recursive_horizons/nsc_metric_horizon_frame.py`.

## Exact metric and two charts

Use the existing compact coordinate q=pi/2+atan(rho) and profile

$$W(q)=3(\pi-q)+\tfrac32\sin(2q)-\sin^2q.$$

Directed bisection and a strictly negative derivative enclose its simple root
q_h. Set u=q-q_h and H(u)=-W(q_h+u)/u, analytically extended at zero.
H(0)=2 kappa_geometry. The stored source kappa remains an occupation parameter;
it is not silently substituted for this geometric derivative.

Let P(u)=-m csc(q_h+u) sigma_1+ell sigma_2. In the exterior, t=sqrt(u),

$$\frac{dF}{d\log u}=K(t)F,\quad
K(t)=\frac{iE}{H(t^2)}\sigma_3+
\frac{t}{\sqrt{H(t^2)}}P(t^2).$$

Inside, u=-delta and t=-i sqrt(delta). Continue the regular coefficient series
in t, while using exp(lambda_s log(delta)) for the two oscillatory factors.
This gives precisely the anti-Hermitian interior generator
iE/H(-delta) sigma_3-i sqrt(delta/H(-delta)) P(-delta).
The exterior coordinate conversion is exact: A(rho)=-r(rho)^2 W(q),
and d rho/d q=r(rho)^2. Hence the original radial Dirac equation must be
multiplied by u*r(rho)^2 to obtain d/dlog(u). Using rho-h in place of
that Jacobian produces a spurious chart discrepancy; a regression checks
the actual radial Riccati owner and detects that omitted factor.

Thermal occupations and horizon-pair coherence stay in the original source
covariance. They are not inserted again by analytically continuing log(u).

## Recurrence and a bound on the entire omitted series

For each sign s, lambda_s=i s E/H(0),

$$F_s=u^{\lambda_s}\sum_{n\ge0}c_{s,n}t^n,
\qquad c_{s,0}=e_s,$$

$$\left(\tfrac n2 I+\lambda_s I-K_0\right)c_{s,n}
=\sum_{j=1}^nK_jc_{s,n-j}.$$

The inverse diagonal matrix has induced 1-norm at most 2/n because E and
H(0) are real. This also covers zero energy and either angular sign.

Choose a complex t-disk |t|<=R and U=R^2. For n>=2, the real-center derivatives
of W obey |W^(n)|<=2*2^n. Thus on |u|<=U,

$$|H(u)-H(0)|\le d_H=2\frac{e^{2U}-1-2U}{U}.$$

The areal-radius denominator is bounded below by

$$d_s=\sin q_h-|\sin q_h|(\cosh U-1)-|\cos q_h|\sinh U.$$

Require H(0)-d_H>0 and d_s>0. The analytic square-root branch is then defined
throughout the disk. A bound on the induced 1-norm of K-K0 is

$$M=\frac{|E|d_H}{H(0)(H(0)-d_H)}+
R\frac{|m|/d_s+|\ell|}{\sqrt{H(0)-d_H}}.$$

Cauchy's coefficient estimate gives ||K_j||_1<=M/R^j. The recurrence is
therefore majorized columnwise by

$$\sum b_n x^n=(1-x/R)^{-2M},\qquad
b_nR^n=\frac{(2M)_n}{n!}.$$

With xi=|t|/R<1 and truncation order N, bound the omitted tail by its first
term divided by 1-r, where

$$r=\xi\max\left(1,\frac{N+1+2M}{N+2}\right)<1.$$

This follows by bounding **every subsequent** positive-term ratio. It is an
analytic remainder bound, not an observed-order or refinement estimate.
Arb/Acb encloses the root, recurrence and finite arithmetic. Component balls
include the entire column l1 remainder. The root bracket represents one exact
root; setting its constant W coefficient to zero uses that proved identity.

## Phase-to-horizon sewing

At the last certified exterior node, the decaying spinor is proportional to
v=(1,exp(i theta)). Theta includes the earlier directed radial transport error.
The exact frame F determines coefficients c=F^-1 v. Its determinant and the
unknown common amplitude cancel in

$$R_{\mathrm{direct}}=\frac{-F_{21}+F_{11}e^{i\theta}}
{F_{22}-F_{12}e^{i\theta}}.$$

The denominator must exclude zero. This encloses the asymptotic reflection in
the explicitly defined direct compact-coordinate frame, including the whole
remaining horizon collar. It is not a claim that the finite-collar producer's
R equals this value in a different phase convention.

For the representative group14/low16_1 row0, the prior phase upper is about
9.04e-12. Its last radius maps to compact distance about 2.21514e-5. At order16
and R=0.1 the new frame tail is below 2.50e-22. The reflected phase uncertainty
therefore remains dominated by the existing finite-band phase bound.

The input record and its NPZ trace are authenticated by hashes. The calculation
uses the original energy label, signed ell and mass with their analytic-label
rounding enclosures. It does not reuse an A2 row from a different energy mesh.

## Coherent phase conventions and the remaining source task

For S(R,T)=[[0,1,0],[R,0,sqrt(T)]], a horizon phase change D=diag(e^(ic),e^(-ic))
satisfies

$$D S(e^{2ic}R,T)=S(R,T)\operatorname{diag}(e^{ic},e^{ic},e^{-ic}).$$

The right diagonal commutes with the **full** original three-port covariance,
including its off-diagonal horizon coherence. Thus a common matched phase
change leaves the physical covariance unchanged. Changing only R or only the
frame does not. This freedom does not permit independent physical state changes.

The original low-energy amplitudes are restrictions of recovered PG packet
fields, produced with `solve_jost` and `MassivePGModeResolution.section(jost=...)`.
They are not produced by `PairedHorizonSeedMap.at_radius` and are not slices of
the separate 64-node A2 mesh. The latter preparation's `_massive_reflection`
is a different numerical owner. The next source task must match the actual
producer's conventions, enclose interior transport, and propagate errors
through the authenticated recovery/restriction to the saved columns. The
present record keeps that physical source error null and the local gate OPEN.

## Controls and replay

Tests check conserved interior norm and exterior signed current, direct-metric
integration plus the original radial Riccati/Jacobian, discarded coefficients, the exact zero-channel
solution, phase uncertainty, coherent sewing and domain rejection.
Independent integration is a sign/chart diagnostic; the analytic majorant
supplies the certificate remainder.

```sh
.venv/validation/bin/python scripts/lab.py -m pytest -q tests/test_nsc_metric_horizon_frame.py
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_metric_horizon_frame.py --check
```
