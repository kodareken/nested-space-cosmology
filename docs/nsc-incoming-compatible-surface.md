# From compatible normal functions to surface constraints

The two-direction point plane is not a surface ansatz. Across a surface,
`r_Tz` must equal the derivative of the function `r_T(z)`. Holding the latter
at its baseline everywhere would force the former to vanish everywhere.

A minimal functional continuation of that point plane is

\[
\delta r(T,z)=T\,w(z)+\frac{T^3}{6}U(z),\qquad \delta a=0.
\]

It preserves the same intrinsic metric, lapse, shift and canonical C0
identification at `T=0`. It supplies six mutually compatible derivative
slots in the [existing Cauchy owner](../src/recursive_horizons/nsc_incoming_cauchy_jets.py):

| Derivative change | Value |
|---|---|
| `delta r_T` | `w` |
| `delta r_Tz` | `w'` |
| `delta r_Tzz` | `w''` |
| `delta r_Tzzz` | `w'''` |
| `delta r_TTT` | `U` |
| `delta r_TTTz` | `U'` |

The constructor accepts derivative values; its existing factorial division
must not be applied twice. This family is a candidate for operator analysis,
not selected physical initial data or a history with an assigned duration.

## Reuse and the missing connection

Reuse the [grade and reflection certificate](nsc-incoming-constraint-germ.md),
the same local Euler and spatial reference owners, the fixed physical source,
and the certified `c_u,c_v` responses. The missing connection is the
differential constraint system obeyed by compatible functions, rather than
independent choices of jets at each point. Derive its permitted structure
and identify the principal coefficient test. Stop before choosing a spatial
length, boundary data, free-function profile or numerical root.

The physical source insertion remains constant along the inherited axial
chart in this fixed-C0 experiment: the instantaneous intrinsic metric,
normal frame and state are the same. This does not assert that an altered
normal geometry shares the old global parent preparation.

## Structure fixed by the existing operators

The variable derivative grades are `1,2,3,4,3,4` for
`w,w',w'',w''',U,U'`. Under axial reflection, `w,U,w''` are even and
`w',w''',U'` are odd. The included local and reference operators have total
derivative grade at most4; their denominators depend on the fixed intrinsic
fields and channel gaps, not on these derivative variables. Consequently
their compatible-family pullback has the form

\[
\mathcal E_N=A(w)U+B(w)w''+C(w')^2+D(w),
\]
\[
\mathcal E_\beta=d\,U'+e\,w'''+F(w)w'+S_\beta.
\]

Here `A,B` are at most linear in w, `C,d,e` are constants, `D` is at most
quartic, and `F` is at most quadratic. Their coefficients are fixed action
responses and the retained source constants, not fit parameters. Pure
baseline derivative factors can lower the displayed polynomial degrees.
The asymmetric physical momentum source is the constant `S_beta`; no
reflection symmetry of that state is assumed.

At `w=w''=w'''=U'=0`, this reduces to the established point plane:
`A(0)=c_u`, `C=c_v2`, `D(0)=S_N` and `F(0)=c_v`.
The new quantities `B(0),d,e` cannot be inferred from that old plane.

Changing w also changes the first normal jet, so the old premise
`Delta P1=0` no longer applies. The
[compatible-family reference proof](nsc-incoming-surface-integrability.md)
instead establishes

\[
\Delta P_1=\frac{\ell w}{4r^2\omega^3}
\bigl[(k/a)\sigma_1+m\sigma_3\bigr],\qquad
\operatorname{tr}(V_N\Delta P_1)
=\operatorname{tr}(V_\beta\Delta P_1)=0.
\]

The first-order matrix is generally nonzero; its two constraint contractions
vanish exactly. The actual first spatial Weyl terms vanish with the unchanged
intrinsic spatial jets, while mixed derivatives remain in higher orders.
The inherited integrable powers for orders2–4 and the nonzero finite-momentum
gaps then establish finite included lapse/shift reference changes for this
family. This extends the needed integrability argument without identifying
the new first-order projector with the old one.

## Algebraic elimination and the next principal test

Where `A(w)!=0`, write

\[
p(w)=B(w)/A(w),\qquad q(w)=C/A(w),\qquad s(w)=D(w)/A(w).
\]

The lapse equation fixes

\[
U=-p(w)w''-q(w)(w')^2-s(w).
\]

Substitution into the shift equation gives the compatible differential
equation

\[
[e-dp(w)]w'''
-d[p_w(w)+2q(w)]w'w''
-d q_w(w)(w')^3
+[F(w)-d s_w(w)]w'
+S_\beta=0.
\]

Thus the determinant of the principal matrix

\[
\begin{pmatrix}A(0)&B(0)\\d&e\end{pmatrix}
\]

decides whether eliminating U yields a regular third-order equation near
the baseline. A vanishing determinant requires analysis of the lower-order
compatibility relation; it does not by itself imply non-existence. The
already positive `c_u,c_v` do not settle this new principal test.
An extension of an old point-plane zero need not keep `U'=w'''=0`:
differentiating the lapse equation supplies an additional compatibility
relation, which the surface equations must satisfy.

If the coefficient is nonzero and the remaining coefficients are finite
and analytic, the standard local initial-value existence theorem applies
to this differential equation. That would be a local compatible surface
result, without an assigned interval length. It would still not establish
a global physical Cauchy surface or its parent and endpoint matching.

## Domain information that is and is not owned

The inherited axial coordinate has a continuous Fourier momentum and no
declared finite length or periodic identification. The pointwise incoming
domain does not supply boundary/asymptotic conditions for its new normal
functions. Compact response profiles and numerical PG radial outflow
conditions are not physical axial constraint boundary data. The fixed
rho0 transmission seam does not fill that gap either.

Local operator and integrability analysis can proceed without inventing
such data. A global surface solution requires an admissible function class
and actual boundary/asymptotic or parent-matching information. No remaining
`Gamma_rest` is assigned, and no field evolution or metric timestep follows
from the formal elimination above.
