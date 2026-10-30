# From averaged boundary work to a local bulk construction

The previous [boundary-work relation](nsc-atlas-flow-gravity.md#the-next-order-averaged-work-relation)
was progress: the same finite fluid energy history selected a positive
average expansion coefficient. It did not yet supply a local Einstein
constraint or a bulk metric. This entry supplies that local connection,
applies it to the retained finite vortex, and checks an independent
fourth-order bulk example directly. The full expanding-child objective
remains open; no narrower result is substituted for its completion.

The [compact record](../results/summaries/nsc-atlas-local-work-v1.json)
contains the constraint fields, independent tensor controls and exact
reading/verification scope. One fixed fine-state replay was used, with
eight FFT workers, a 900 process-CPU-second limit and no new refinement
campaign. It used 170.13 process CPU seconds and 171.13 wall seconds.
Only small slices and summaries are stored.

## The useful published connection

[Compere, McFadden, Skenderis and Taylor, 1103.3022v2](https://arxiv.org/html/1103.3022v2)
extend the Rindler fluid/gravity construction beyond cubic order. We read
the seed, radial recursion, fluid frame, fourth-order incompressibility,
metric coefficients and horizon-regularity conditions directly. Their
construction is a derivative expansion, not a convergence theorem.

Three routes were considered: extrapolate the old average relation,
embed an arbitrary expanding boundary in an unrelated flat bulk, or use
the local constraint and reconstruction equations in the same fluid/gravity
sector. The third retains the source, frame and radial boundary conditions.
An exact analytic flow is used to verify the tensor implementation; it does
not replace the retained vortex as the main numerical state.

The likely failures are inversion of a nonzero Poisson zero mode, confusion
between local and globally periodic data, and calling a coordinate dilation
the same metric as a shift-free expanding boundary. Each is checked below.

## 1. The average expansion is the zero mode of a local equation

Use the published shear convention sigma_ij=partial_i v_j+partial_j v_i.
With slow coordinates y=epsilon x, t=epsilon^2 tau and
v=epsilon u+epsilon^3 w, the fourth-order scalar constraint is

$$
\nabla\cdot w=S,\qquad
S=\frac{u\cdot\nabla p}{\nu}-u\cdot\Delta u
+\frac12\sigma_{ij}(u)\sigma_{ij}(u).
$$

For a periodic divergence-free leading field,

$$
\int u\cdot\nabla p=0,\quad
-\int u\cdot\Delta u=\int|\nabla u|^2,\quad
\frac12\int\sigma_{ij}\sigma_{ij}=\int|\nabla u|^2.
$$

The last equality uses integration by parts; the corresponding integrands
are not generally equal pointwise. Consequently,

$$
\langle S\rangle=\frac{2D}{V\nu}=3b_2'.
$$

This independently connects the earlier averaged work coefficient to the
local scalar channel of the Einstein constraints. The mean alone is the
same enstrophy balance; the nonzero spatial modes below are the new work.

If D>0, no strictly periodic w can solve this divergence equation, since
the integral of a periodic divergence vanishes. A boundary flux or changing
regional description is required. The explicit local completion is

$$
w=\beta y+\nabla\varphi+w_T,\qquad
\beta=\frac{\langle S\rangle}{3},\qquad
\Delta\varphi=S-\langle S\rangle,\quad
\langle\varphi\rangle=0,\quad\nabla\cdot w_T=0.
$$

Only the zero-mean source is inverted. The gradient solution is the
minimum-norm periodic correction: its inner product with a periodic
divergence-free addition is zero. We select that representative for this
calculation; a transverse field remains free until its next-order evolution
and state are specified. Opposite faces differ by beta L in their normal
component, giving total outward coefficient 3 beta V=integral S.

## 2. The retained vortex satisfies the local completion

The fine 128-grid, dt=.0005, t=.05 case was replayed once with unchanged
initial data and source. Final energy and dissipation matched its saved
values at the reported precision. The spatial constraint uses twice-sized
padding: the pressure source is quadratic and u dot grad(p) is cubic, so
the original strict 2/3 velocity band requires this larger source grid.
It is a dealiased evaluation, not a new physical-grid campaign.

| State | Mean S | beta=mean S/3 | Local scalar residual RMS | Periodic correction RMS |
| --- | --- | --- | --- | --- |
| t=0 | 53.5644254361 | 17.8548084787 | 1.59e-13 | 78.35803697 |
| t=.05 | 51.7916146209 | 17.2638715403 | 1.64e-13 | 75.80962030 |

The relative residuals are about 4.45e-16. These measure the constraint
solve for the retained band-limited data, not the continuum evolution error
or the total Einstein residual of the numerical vortex. The source has
nontrivial positive and negative local structure even though its mean is
positive. Its mean agrees with 2D/(Vnu), independently evaluated from the
curl norm, to rounding precision.

The pressure here is the full quadratic Poisson pressure of each retained
band-limited snapshot. The constraint requires no time derivative; it is
therefore evaluated without silently replacing the projected Galerkin
time derivative by an exact continuum history.

## 3. A horizon-regular local bulk example is checked directly

To independently test the reconstruction, use the exact unforced shear
flow u=(A(t) sin(y+phase),0,0), p=0 and A=exp(-nu t). Its scalar source
is S=A^2, so w=A^2(x,y,z)/3 is the required affine correction. Insert this
field into the published seed and fourth-order metric coefficients in
their stated fluid frame.

The implementation builds all metric first and second derivatives with
the hydrodynamic derivative weights, inverts the metric as a formal
epsilon series, forms the Levi-Civita connection and contracts its Ricci
tensor. It does not use the field equation as a substitute for this tensor
calculation. The metric coefficients are polynomial in the bulk radial
coordinate and are finite at the background future Rindler horizon r=0.

Two independent levels of checking are retained:

- Numerical contractions at r/nu=0,.5,1 and nu=1,.2 give maximum Ricci
  coefficients through fourth order below 9e-15. Omitting w leaves
  fourth-order residuals of .4094--.4804.
- An exact symbolic calculation at a nontrivial shear-flow jet, with the
  radial coordinate left symbolic, gives zero for all 125 component/order
  entries through fourth order. The point and parameters are in the record;
  this is not a claim that the entire general-profile tensor was symbolically
  simplified here.

A curved de Sitter-times-line control also verifies that the tensor engine
returns the expected nonzero Ricci tensor. The general radial reconstruction
and its horizon boundary condition remain attributed mathematics; the
analytic witness is our independent application and check. The reconstructed
metric of the full numerical vortex has not yet been independently contracted.

## 4. Local expansion and global geometry must still be connected

The affine correction is valid on an open local patch. It does not descend
to an unchanged periodic vector field on the original torus. That is the
geometric reason a larger region, flux boundary or changing identification
must be supplied before calling the construction global.

There is also a precise distinction between a moving lattice and a new
intrinsic metric. Pull back a flat boundary by X=b(tau)x:

$$
-\nu\,d\tau^2+dX^2
=-(\nu-\dot b^2|x|^2)d\tau^2
+2b\dot b\,x_i\,dx^i d\tau+b^2dx^2.
$$

This is still flat and contains shift and lapse terms. The earlier
shift-free ansatz -nu d tau^2+b^2 dx^2 is a different metric in general.
On the slow spatial domain the shift enters at the relevant expansion
order. Equality of beta and b2' is a coefficient-level connection, not a
proof that these two geometries are identical. This distinction prevents
a coordinate change from being presented as a new cosmic region.

Extending this expansion further requires the fifth-order spatial Ward
identity, which determines the transverse part of w and the corresponding
pressure correction. The programme's next outward connection is its actual
four-dimensional parent and global attachment; higher-order polishing is
useful only when that connection requires it. This work closes the local
scalar constraint and supplies an independently checked local bulk witness.
The global region and finite-amplitude remainder still need their common
geometry and state.
