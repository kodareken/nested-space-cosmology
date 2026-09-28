# A flux form of the compatible surface equations

The completed coefficients give an exact alternative to eliminating the
lapse equation by division through A(w):

\[
\boxed{\frac{d}{dz}\bigl[\Delta(w)w'\bigr]
=A(w)\,[\Pi_0-S_\beta z-\mathcal F(w)]+dD(w)},
\]
\[
U=\frac{\Pi_0-S_\beta z-e w''-\mathcal F(w)}d,
\qquad \mathcal F_w=F,\quad\mathcal F(0)=0.
\]

Here `Delta=A*e-B*d`, and all coefficients are the fixed included operators.
`D(w)=S_N+sum(D_j*w^j)` retains the exact source constant. No boundary value
or integration constant is selected by this representation.

## Reuse and the missing connection

Reuse the [complete polynomial coefficients](nsc-incoming-surface-coefficients.md),
[principal matrix](nsc-incoming-surface-principal.md), and
[inherited momentum balance](nsc-incoming-surface-momentum-balance.md).
The missing connection is a surface representation that preserves the
momentum first integral and distinguishes an avoidable lapse pivot from
actual principal rank loss. Combine the existing identities and stop at
the equivalent equations and crossing condition. No profile integration,
source calculation, boundary choice or metric timestep is needed.

## Exact reduction

The coefficient certificates establish

\[
A=A_0+A_1w,\qquad
B=-\frac{2A+(H_r+w/r)d}{a^2},\qquad e=-d/a^2,
\]
\[
C=-\frac{A_1+d/r}{a^2},\qquad
\Delta_w=\frac{d(A_1+d/r)}{a^2}=-dC.
\]

The shift equation integrates to
`Pi=dU+e*w''+mathcal F(w)=Pi0-S_beta*z`. Combining d times the lapse
equation with this identity eliminates U without dividing by A:

\[
\Delta w''-dC(w')^2
=A(\Pi-\mathcal F)+dD.
\]

Since `Delta_w=-dC`, its left side is exactly `(Delta*w')'`.
Conversely, reconstruction of U from the first integral satisfies the
lapse equation by this equality, while differentiation of the first
integral satisfies the shift equation. Thus this is an equivalent
representation of both included constraints wherever the reconstructed
functions are sufficiently smooth. The division by d is permitted by
its already certified strictly negative interval.

On `Delta!=0`, a first-order system is

\[
w'=J/\Delta(w),\qquad
J'=A(w)[\Pi_0-S_\beta z-\mathcal F(w)]+dD(w).
\]

The independent variable remains spatial z. This supplies a useful
representation for later numerical work but makes no claim about a
numerically integrated profile or its source error.

## Two distinct normal-data roots

The [regular-branch certificate](nsc-incoming-surface-regular-branch.md)
places the principal root in
`[0.8828574152795067,0.882857415279511]` and the lapse-pivot root in
`[0.9860722276629823,0.9860722276629867]`.
These are values of `w=delta r_T`, not positions or durations.

At the lapse-pivot root, `Delta!=0` and `d!=0`; the flux representation
therefore stays regular. The pole created by solving the lapse equation
for U is avoidable. Starting from w=0 still encounters the principal root
first, so this observation alone does not connect the two branches.

At the principal root `w_D`, define

\[
\Delta_1=-dC>0,\qquad
R(w,z)=A(w)[\Pi_0-S_\beta z-\mathcal F(w)]+dD(w).
\]

Any finite C2 crossing at z_c must satisfy

\[
R(w_D,z_c)=\Delta_1\,[w'(z_c)]^2.
\]

The [analytic crossing result](nsc-incoming-surface-crossing.md) proves
sufficiency for a nonzero slope when this compatibility condition holds:
the inverse-coordinate Briot–Bouquet problem has state eigenvalues 0 and -2,
so it has a unique analytic germ. This does not show that any particular
trajectory reaches the root. The constants `S_N,S_beta,Pi0,z_c` remain
symbolic; no source value, crossing location or physical initial tuple is
assigned.

### The inherited flux excludes a zero-slope classical contact

The source-specific momentum witness strengthens this conclusion. Suppose
`w` is C3 and `U` is C1 near an interior point with `w=w_D` and `w'=0`.
Differentiating the lapse equation at this point gives

\[
A U'+B w'''=0.
\]

All other terms contain `w'`, and the fixed-C0 source constant has zero
spatial derivative. At the principal root, `Ae=Bd` and `A>0`, hence

\[
E_\beta=dU'+e w'''+S_\beta
=\frac dA(AU'+Bw''')+S_\beta=S_\beta< -M<0.
\]

This contradicts the shift constraint. Thus an admissible classical
interior contact with the principal root must have nonzero slope and
`R>0`. This corollary uses the already certified inherited flux; it does
not extend the analytic theorem's uniqueness claim to all smooth germs,
prove extension through a one-sided endpoint, or exclude other NSC
families.

Both roots occur at finite derivative values with the same positive
intrinsic metric and finite reference insertions. Neither identifies a
spacetime singularity. Principal compatibility and global parent/boundary
matching still need their own conclusions. The full extended-stationarity
gate remains OPEN; no missing `Gamma_rest` is assigned.
