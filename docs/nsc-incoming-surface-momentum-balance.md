# The inherited current fixes the axial momentum balance

For the compatible two-function family, define the constraint primitive

\[
\Pi(z)=dU(z)+e w''(z)+\mathcal F(w(z)),\qquad
\mathcal F(w)=\int_0^w F(s)\,ds.
\]

The [included shift equation](nsc-incoming-compatible-surface.md) is exactly

\[
\boxed{\Pi'(z)=-S_\beta},\qquad
\Pi(z_2)-\Pi(z_1)=-S_\beta(z_2-z_1).
\]

This is an integrated constraint relation, not an additional source or a
chosen boundary force. Its integration constant is a value of the normal
data; it is not assigned here.

## Reuse the existing exact current bound

The [incoming constraint gate](nsc-incoming-constraint-gate.md) already
bounds the same retained source on the same intrinsic slice. Its exact
rational witness gives

\[
P_K<-M,\qquad
M=\frac{5887634380411363}{265464355000000000}>0.
\]

The margin is approximately `0.02217862499999807` in the owned raw
Killing-power normalization. It combines the analytic unit-transmission
LLL power with an upper bound on every retained massive positive
contribution, using only `0<=transmission<=1`. Horizon coherence cancels
from this current by the exact sewing identity; it is not removed from
the state. No scattering curve, low-band field or quadrature is needed
to establish the sign.

At the unchanged incoming slice,

\[
T_{\hat0\hat1}=-\frac{P_K}{4\pi r^2a^2},\qquad
\mathcal E_\beta^{\rm matter}=-4\pi a^2r^2T_{\hat0\hat1}=P_K.
\]

The second equality is the owned
[raw action conversion](../src/recursive_horizons/nsc_incoming_joint_constraints.py).
The baseline homogeneous geometric, local and reference momentum is zero;
the compatible spatial changes are the terms in Pi'. Consequently the
same fixed-C0 constant obeys

\[
S_\beta=P_K<-M,\qquad \Pi'(z)>M.
\]

The old gate's normalized miss, about `0.013757531`, is
`T01/(2*A_EH)`, not raw `S_beta`; `A_EH` here denotes the locked Einstein
coefficient, not the function A(w) in the compatible lapse equation.
Neither that normalized number nor a neck residual is substituted into
the present raw equation.

## Necessary global conditions, without choosing a boundary class

Every continuation in this included two-function family must have the
stated affine growth of Pi. Therefore:

- Periodic Pi with any positive period is incompatible with the inherited
  nonzero current.
- Uniformly bounded Pi on an unbounded axial interval is incompatible
  with the same constant source.
- On any finite axial interval, the endpoint difference must equal
  `-S_beta` times the coordinate interval length. No length or endpoint
  value is chosen by this statement.

These are conditional restrictions on boundary data. The current inherited
chart supplies neither a periodic identification nor a prescribed bounded
class for the new normal functions. Unbounded Pi or nonperiodic endpoint
differences are not excluded by this argument. Local compatible existence
is also unaffected.

The relation makes the remaining global question concrete: the parent and
endpoint matching must supply admissible normal data carrying this balance,
or the candidate functional family must be reconsidered within the declared
action. A periodic numerical grid or a compact response profile cannot
silently supply that physical choice.

## Decision and scope

Reuse the exact flux witness and the polynomial structure of the compatible
shift equation. The missing connection is the boundary balance forced by
the actual source, before imposing a global function class. One symbolic
integration and the existing raw normalization close that connection; no
old source, mode, reference or scattering calculation is rerun.

This concerns the same retained fixed-C0 experiment and included bulk
operators. It does not assume the altered future metric retains the old
Killing symmetry. It is not an omitted-tower bound, a global NSC
non-existence result, or a prescription for `Gamma_rest`. No source state,
coupling, physical initial tuple, metric timestep or boundary condition
is changed.
