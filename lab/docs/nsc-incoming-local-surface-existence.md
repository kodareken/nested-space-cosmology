# Local compatible solutions of the retained incoming constraints

The included retained incoming constraints admit unique smooth local
compatible solutions within the stated two-function family for each finite
pair of spatial initial derivatives specified below, with fixed canonical C0.
This follows from the certified
principal matrix and the standard local ODE existence theorem. It is
separate from global initial data, parent preparation, transmitting
endpoints and full extended stationarity.

## Reused inputs and the decisive gate

Reuse the [finite exact source](nsc-incoming-source-finiteness.md),
[compatible-family structure](nsc-incoming-compatible-surface.md), and
[reference integrability](nsc-incoming-surface-integrability.md). The
principal gate is now supplied by the
[compatible-family coefficient certificate](nsc-incoming-surface-principal.md):
`A(0)!=0` and `A(0)*e-B(0)*d!=0`. In particular,

\[
A(0)\in[0.057399277586276304,0.05739927758627644],
\]
\[
\det\begin{pmatrix}A(0)&B(0)\\d&e\end{pmatrix}
\in[-0.0029161631968961416,-0.0029161631968961255].
\]

The eliminated third-derivative coefficient lies in
`[-0.05080487628982585,-0.050804876289825454]`. These intervals establish
the needed nonzero values. The existence argument requires no new source
integration or metric calculation.

The missing connection is from a compatible differential constraint system
to functions satisfying both constraints on a neighborhood, rather than
independent pointwise roots. Stop after checking the principal gate and
the hypotheses of the existing ODE theorem. Do not choose a numerical
interval length, physical boundary data or an initial-data tuple.

## First-order system in the spatial coordinate

The owned family is

\[
\delta r=T w(z)+T^3U(z)/6,\qquad\delta a=0,
\]

with its forced companion mixed derivatives. Its two included constraints
have the exact form

\[
\mathcal E_N=A(w)U+B(w)w''+C(w')^2+D(w),
\quad
\mathcal E_\beta=dU'+e w'''+F(w)w'+S_\beta.
\]

The functions `A,B,D,F` are fixed finite polynomials of the degrees already
established by grade and reflection. `C,d,e` and the exact source constants
are finite. No coefficient is fitted or supplied as an extra force.

Set `p=B/A`, `q=C/A`, `s=D/A`, and `K=e-d*p`. With `Y=(w,v,xi)`, define

\[
\frac{dY}{dz}=
\begin{pmatrix}
v\\
\xi\\
\displaystyle
\frac{d(p_w+2q)v\xi+dq_wv^3-(F-ds_w)v-S_\beta}{K}
\end{pmatrix},
\qquad
U=-p(w)\xi-q(w)v^2-s(w).
\]

All functions on the right are evaluated at w. The independent variable
is the axial spatial coordinate z, not physical elapsed time.

## Existence on a neighborhood

Nonzero `A(0)` and determinant imply nonzero `A(w),K(w)` on some open
interval J containing w=0, by continuity. The vector field above is smooth
on `J x R^2`, and therefore locally Lipschitz. For every finite real pair
`v0,xi0`, the initial value `Y(0)=(0,v0,xi0)` is in that domain.

The Picard–Lindelöf theorem gives a unique local solution on an interval
around z=0. Smoothness follows from the smooth vector field. This imports
[Teschl, Theorem2.2 and Lemma2.3, printed page38](https://www.mat.univie.ac.at/~gerald/ftp/book-ode/ode.pdf#page=49);
the contraction theorem is not rederived here.

Reconstruct U with the displayed formula. The lapse constraint vanishes
identically by that definition. Differentiating it and using the third
component of the ODE makes the shift constraint vanish identically as well.
The functions consequently supply compatible mixed derivatives at every
point of the local interval, including `U'` and `w'''`.

The corresponding local Taylor metric is smooth and has positive a,r near
the initial point because their unchanged intrinsic values are positive.
This is a kinematic realization of the initial jets. It is not a solution
of the spacetime evolution equations or a metric timestep.

The old algebraic values `(u*,v*)` can be recovered by the formal choice
`v0=-S_beta/c_v`, `xi0=0`: then U(0) equals its old lapse solution. A surface
extension generally has nonzero `U'(0),w'''(0)`, determined by compatibility.
This observation does not assert a certified numerical value or select
that tuple as physical initial data.

## What remains outside this result

The theorem supplies an unspecified local spatial neighborhood. It does
not assign a physical axial length, duration or periodic identification,
and it does not prove continuation over the entire inherited axial chart.
No admissible global boundary/asymptotic or parent-matching data have been
provided for these new normal functions.

The exact retained source is used in the mathematical statement. Computing
particular solutions with a numerical residual tolerance still requires
source and coefficient error propagation; the unfinished low/subgap
numerical work is not replaced by this existence theorem. Nor does this
prove completeness of the retained physical inventory.

Full parent preparation, all transmitting endpoint variations and the
remaining declared-action equations are separate requirements. No
`Gamma_rest` is set to zero or assigned a compensating value. Extended
EXISTENCE and observational predictions do not follow from this local
constraint construction.
