# Responsive local collar and required parent traction

This calculation uses the unchanged leading Einstein/Dirac action on the
finite symmetric collar `|x|<=1`, with the child readout at `|x|<=.5`.
It relates a declared local standing source to radius, lapse, radial length
and the boundary gradient charge that an exterior continuation must supply.
It does not complete a global static parent, global bound spectrum, physical
vacuum matching or autonomous renewal.

The locked common coefficients are `g=8pi A`, `mag=2pi C_F flux²`,
`mag/g=1`, `M=4`, `kappa=1`, with fixed coordinate frequency `epsilon=1.5`.
Set `h=log Q`, `n=u²+v²`, `S=2uv` and `K=epsilon n-kappa Q S`.
The independently reviewed static equations are

```text
r_x=s, h_x=v_h
r_xx=Q²(r²-mag/g)/r - s²/r - Mw kappa Q S/(2gr)
h_xx=Q²-3r_xx/r
u_x=epsilon v-kappa Q u
v_x=kappa Q v-epsilon u.
```

The first integral is
`I=Q²(r²-mag/g)-3s²-2rs h_x+Mw K/g`.
Symmetric center data `Q=1,r_x=h_x=0,u=v=1/sqrt(pi)` enforce
`r0²=mag/g-Mw/(g pi)` and `I=0`. The weight `w` multiplies the local
source amplitude once; multiplicity remains outside the one-particle
operator. It is not assumed to be a globally normalized CAR occupation.
The rank-one covariance eigenvalue on a finite collar is
`w integral n dx`, which is measured and checked separately.
The real standing spinor has zero current.

The baseline BR geometry is `r=1,Q=sec x`, with
`u=cos x sqrt(1+sin x)/sqrt(pi)` and
`v=cos x sqrt(1-sin x)/sqrt(pi)`.
The first-order response `r=1+w a`, `h=log(sec x)+w q` obeys

```text
a=−M/(4pi g) [2+sin²x+3x tan x]
q_xx−2sec²x q=−3a_xx, q(0)=q_x(0)=0.
```

The predictor uses this analytic radius response and a linear equation for
`q`; the exact nonlinear measurement integrates the actual `u,v,r,Q` source in x
without importing the predicted endpoint. No force is fitted or added.
The held source weight is fixed at `.0013` before measurement. Tests use
other small weights `.0003,.0006` and check the expected quadratic prediction remainder.
An earlier `.002` scratch calculation was observed before sealing and is
development evidence; it is not the new held measurement.

At each declared cut, the right endpoint charge is
`B(b)=g[r r_x+r² h_x](b)`. Reflection gives `B(-b)=-B(b)`;
the full-collar bracket is `B(b)-B(-b)=2B(b)`. The exact static identity is

```text
B(b)−B(−b) = M w epsilon integral[-b,b] n dx
             + mag integral[-b,b] Q² dx.
```

Both endpoint charges and the bracket are reported to preserve the factor
two. This is required exterior boundary data, not an imposed material wall
or a completed pressure-supported parent. Exterior matching still needs the
endpoint geometry/gradients, Dirac phase and frequency continuation, global
mode normalization/state, and compatible boundary charge.

The metric proper lapse and radial length density are `N=rQ`. In these
coordinate units `N(0)=r0`, rather than one; the separately reported ratio
`N(x)/N(0)` refers to the center proper clock. The center proper frequency
is `epsilon/N(0)`. Proper radial length is `integral[-b,b] rQ dx`.
Changing source weight keeps the declared coordinate frequency and center
amplitude fixed; it does not retune a global normalized mode.

Default operation only previews authenticated common-action inputs. Root
freezes all four producer files before executing the immutable stages:

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_parent_traction.py --prepare
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_parent_traction.py --predict
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_parent_traction.py --measure
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_parent_traction.py --check
```

The default prefix directory is
`lab/results/development/nsc-discovery-parent-traction-v1`;
`--output` may select a new child directory. JSON/NPZ stages are exclusively
created, payload/array/input/source bound and limited to 64 MiB. The prediction
is locked before either held nonlinear IVP. Two DOP853 tolerances and nested
quadrature grids provide numerical indicators within an aggregate 30 CPU
second budget; they are not continuous error certificates. Read-only checks
replay stored endpoints and predictions and authenticate original producer
bytes through their recorded Git commit. Historical evidence is never healed.
Explicit stage creation requires a non-None immutable producer commit whose
bytes match the current producers. If current producers differ from an
authenticated historical stage, `--check` reports authentication only and
does not claim numerical replay using current equations. Measurement verifies
the prepared-stage binding before either solve and requires the captured
prediction hash to remain unchanged before writing its result. Both declared
cuts must occur on every solver grid. Replay derives the endpoint gradients,
charges and center clock frequency directly from the saved numerical arrays.
PREPv2 supplies the common coefficient authentication, not the local collar
spinor or a reset of an existing trajectory.
