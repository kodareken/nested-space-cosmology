# Local lapse and shift forces on the incoming Cauchy data

The existing light and locked compact actions now supply local Euler
coefficients on the same incoming KS surface:

$$
\delta S_{\rm local}=\int dT\,dz\,
  (\mathcal E_N\,\delta N+\mathcal E_\beta\,\delta\beta),
\qquad F_B^{\rm local}=-\mathcal E_B.
$$

On the inherited incoming geometry the evaluated local action gradients are
approximately $(\mathcal E_N,\mathcal E_\beta)=(-0.563010381,0)$.
The zero shift value is resolved numerically; the record retains its raw
value. A mixed normal/spatial derivative control produces a nonzero shift
force and passes an independent weak-variation comparison.

The [owner](../src/recursive_horizons/nsc_incoming_local_constraints.py)
applies the Euler operator to the existing densities, without introducing
metric equations from another action:

$$
\mathcal E_B=\sum_{|\alpha|\le2}(-\partial)^\alpha
\frac{\partial\mathcal L}{\partial B_\alpha}.
$$

`IncomingCauchyJets` supplies derivatives through total order four.
Lapse and shift are varied independently before imposing their incoming
values. The intrinsic geometry, normal, canonical frame and state remain
unchanged. The finite Taylor germ and complex Cauchy contours are numerical
differentiation devices; neither selects a history or a physical duration.
The mixed derivative $B_{Tz}$ is one independent variable and therefore has
no extra factor of two in the Euler sum.

The light action is evaluated in its owned KS foliation. An identity
integration adapter extracts the exact densities returned by
`LightRestorationAction.actions` and `LockedSphericalLocalAction.actions`.
No curvature or restoration formula is copied into the new owner. The light
geometric channels and locked compact action are each included once. Euler
bulk cancellation is checked, while constant Euler and BoxR endpoint forces
remain separate from these local bulk equations.

The [record](../results/development/nsc-incoming-local-constraints.json)
keeps every channel gradient, its opposite force, resolution indicators and
the raw Euler cancellation. Independent controls compare:

- All four light metric variations with the inherited homogeneous tensor
  on the actual incoming surface.
- All four compact Weyl variations with the committed charged-neck tensor.
- The two pointwise constraint coefficients with a direct compact weak
  variation of the same densities on a nonhomogeneous Taylor control.
- Two numerical Cauchy radii as well as two contour resolutions.

The numerical verification tolerance is $3\times10^{-11}$. This does not
change the separate constraint stationarity tolerance. The physical
state-minus-reference source and reference-band remainder are absent from
this owner, so these numbers are the local contribution to two constraints,
not a complete constraint residual or a solution of all four metric
equations. No normal data, root, state or history is selected.

```sh
python3 scripts/derive_nsc_incoming_local_constraints.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_local_constraints.py
```
