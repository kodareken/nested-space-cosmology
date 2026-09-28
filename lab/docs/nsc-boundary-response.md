# First finite curved Dirac boundary response

This note records the first numerical parent/child boundary maps of the
first-order spatial Dirac operator on the nested-space throat. It is the
finite calculation demanded by remaining-closure items 2–3 of
[`nsc-closure-verification-2026-09-07.md`](nsc-closure-verification-2026-09-07.md)
and by Stages 3/4 of the accepted spectral plan, with the covariant-action
and Lorentzian-Hamiltonian limitations left explicit.

Geometry-source checkpoint: `5f38712ca01ddd71e715fd265088925a73369aba`.
Schema: `NSC-3-BOUNDARY-RESPONSE-v2`.
Machine record: [`results/nsc-3-boundary-response.json`](../results/nsc-3-boundary-response.json).
Code provenance is the SHA-256 of the four owned source files, stored in the
record at generation time.

The calculation does not select a stationary scale, an electron mass, a
complete covariant action, a Lorentzian graviton, or a nonlinear transition.

Independent root results that this note **cites but does not edit**:

- [`nsc-radial-spectrum.md`](nsc-radial-spectrum.md): a Weyl sequence proves
  `spec_ess = R` for the unwarped complete-line spatial Dirac.
- [`nsc-geometric-chain.md`](nsc-geometric-chain.md): repeating the finite
  motif produces a spatial gap. That is a different geometry than the isolated
  box used here.

## Claim

The spatial throat Dirac operator

```text
r = sqrt(1 + rho^2),   w = kappa / r,
D = [[ 0,   -d + w ],
     [ d + w,   0  ]]
```

is a self-adjoint first-order operator in the flat radial measure obtained
after the spinor rescaling that produces `A = d + w` and `A^dagger = -d + w`.

In the black-universe chart the **parent** is `rho > 0` (asymptotically flat)
and the **child** is `rho < 0` (cosmological). On a finite box with the
separated condition `u(+-R) = 0`, a staggered mimetic discretization with
sparse LU yields:

- energy-resolved parent and child Weyl maps at nonreal energy, with
  `u(0) = 1` so the data have no SVD phase;
- an oriented DtN jump `N_parent + N_child = E (m_child - m_parent)` in which
  the `w(0)` terms cancel;
- a joined resolvent that agrees with its Schur complement;
- second-order agreement with an independent `solve_ivp` Cauchy map, both at
  the staggered edge and at a reconstructed throat value;
- a regulator-independent finite response whose imaginary part has the sign
  of `Im E`.

The decaying superpotential does not produce an asymptotic mass gap. Finite-box
positive eigenvalues scale as `1/R` by a product/ratio criterion, not merely
by `E_1 R > 2`. They are not massive bound poles. The compact even mode of
`Q_y` does not generate a radial mass; multiplication by `exp(-sigma)` mixes
`y` modes, so the average-warp factor is a **projected approximation** with a
measured leakage, not an invariant-subspace theorem.

The spectral Dirichlet-to-Neumann map of `D^2` is `E m(E)`, not a constant
matching field `Phi`. The mixed edge trace `v(+-h/2)/u(0)` is not the
continuum Calderon map.

## Orientation and outward normals

| Side | Interval | Asymptotics | Outward normal at `rho = 0` |
|---|---|---|---|
| Parent | `[0, R]` | `R_4 -> 0` as `rho -> +inf` | `-d/d rho` |
| Child | `[-R, 0]` | `R_4 = -36 pi` as `rho -> -inf` | `+d/d rho` |

Ordinary logarithmic derivatives obey `u'/u = -w + E m`. Therefore

```text
N_parent = - (u'/u)_parent =  w(0) - E m_parent
N_child  = + (u'/u)_child  = -w(0) + E m_child
N_parent + N_child         = E (m_child - m_parent)
```

The `w(0)` terms cancel exactly. This identity is algebraic once `m = v/u` and
the ODE for `u'` are used; it is not a fit. It is not an identification of
that jump with `Phi`.

## Operator, measure, and discretization

After the four-dimensional conformal reduction of the global `tau` slice
(`AGENTS.md` §20.39–20.40), the radial partners on the areal radius
`r = sqrt(1+rho^2)` are

```text
A = d/d rho + w,     A^dagger = -d/d rho + w,     w = kappa / r.
```

The two-component operator is exactly `D` above. Then

```text
D^2 = diag(-d^2 + V_minus, -d^2 + V_plus),
V_plus  = w^2 + w',
V_minus = w^2 - w'.
```

The inner product is the flat radial `L^2` measure in which `A` and
`A^dagger` are adjoints.

Separated self-adjoint data at a regular endpoint are Lagrangian lines for
the Green's current

```text
J(psi, chi) = v_psi^* u_chi - u_psi^* v_chi,
alpha = [[0, -1], [1, 0]].
```

The condition `u = 0` is one such line. Joined calculations use `U = I`.

The local probability-current control uses complex Cauchy data `(1,i)/sqrt(2)`
at real energy, giving nonzero flux. It is a local conservation probe, not a
reflecting-wall eigenstate. Using real standing-wave data alone would give
identically zero current and would not exercise this numerical check.

Naive centred first derivatives on a collocated grid are a negative control.
They are not the continuum proof.

The physical scheme is the one-dimensional Whitney pair: `u` on nodes, `v` on
edges, `d` the incidence matrix, `w` on edges. After a square-root mass
rescaling the Hamiltonian is Euclidean-Hermitian. Same-chirality Dirichlet
walls make the node and edge counts differ by one, so the joined massless
operator has a single `E = 0` mode supported on `v`. That mode is a
finite-interval index, present at every box radius. It is not a mass gap.

Primary half-domain solves are **sparse LU** with the throat unknown fixed at
`u(0) = 1`. Dense SVD nullspaces are retained only as a small independent
control. Dense inversion of the joined resolvent is likewise a control, not
the primary path.

## Edge trace versus throat Cauchy data

The staggered unknown `v` lives at `rho = +- h/2`, not at the throat. The
ratio `m_edge = v(edge)/u(0)` is a discrete edge trace. It is compared to
the continuum solution at the same edge. It is **not** the continuum Calderon
projector.

The throat value is reconstructed by the oriented implicit Euler step of the
ODE at `rho = 0`:

```text
parent: v(0) = (v(+h/2) + (h/2) E u(0)) / (1 + (h/2) w(0))
child:  v(0) = (v(-h/2) - (h/2) E u(0)) / (1 - (h/2) w(0))
```

with `u(0) = 1`. That reconstruction is compared to the continuum
`m(0) = v(0)/u(0)`. Both residuals are recorded, and both are required to
converge under mesh refinement. A rank-one projector built from
`(1, m_throat)` is a discrete reconstructed Cauchy projector, not a claim
that the mixed edge data equal the continuum Calderon map.

The supersymmetric DtN of `D^2` remains `N = E m`, with `m` the throat ratio.
Dimensions:

| Symbol | Dimension | Owner |
|---|---|---|
| `m(E)` | 1 | first-order Weyl function at `rho = 0` |
| `N(E^2)` | 1/length | DtN of the squared partners |
| `Phi` | energy | off-diagonal two-sheet matching field |

## Conservation identities

Two different conserved quantities must not be conflated.

- The **Abel Wronskian** `u_1 v_2 - v_1 u_2` of two solutions of the same
  linear ODE is constant for any complex energy. That is linear ODE algebra.
- The **stationary probability current** `2 Im(u^* v)` is conserved along a
  real-energy solution. It is not the Abel Wronskian.
- The **Green boundary form** of self-adjoint `D` vanishes for discrete
  domain functions because the joined Hamiltonian is exactly Hermitian. That
  residual is the Hermitian residual of `H`, not a complex-energy Wronskian.

## Independent continuum comparison and sparse controls

The continuum Cauchy problem is integrated with `solve_ivp` (`DOP853`) from
the exterior wall of the correct half-line (`+R` for parent, `-R` for child).

At the primary family `R = 8`, `kappa = 1`, `h = 0.1`, parent and child maps
at four probe energies have edge and throat relative residuals below `5e-3`.
Halving the parent mesh at the first energy gives edge-trace orders near 2
and a converging throat reconstruction. Sparse LU edge ratios agree with the
dense SVD control at `N = 20` to about `1e-14`. Sparse joined Schur agrees
with dense inversion on that small control.

The quadratic form `f^dagger (H - E)^{-1} f` at the throat probe has `Im` of
the same sign as `Im E`. No heat-kernel regulator enters this resolvent.

## Spectra and the gap statement

Eigenvalues are stored as algebraically sorted signed bands: positive
ascending, negative descending toward zero, and an explicit index-kernel
count. Near-degenerate `+-` pairs are not ordered by `argmin |E|`.

Massless joined spectra at `R = 4, 8, 12` have one index eigenvalue at
`E = 0`. The first **positive** band satisfies the inverse-radius criterion:

- `E_1` decreases as `R` increases;
- `E_1 R` increases, with last relative change below 5%;
- the finest `E_1` ratio matches the inverse-radius ratio to better than 10%.

A massive control `m sigma_3 = 0.75` holds the first positive value at `0.75`.
That is what a genuine mass edge looks like in the same code.

As `|rho| -> infinity`, `w -> 0`. Together with the independent Weyl-sequence
proof that `spec_ess = R` on the complete line, the isolated throat has no
asymptotic mass gap. The finite-box index kernel and `1/R` levels are not
bound poles. The geometric chain of repeated motifs is a different operator
and is allowed to have a gap without contradicting this box.

The compact warped operator
`Q_y = exp(-sigma_edge/2) d_y exp(-sigma_node/2)` has a rectangular kernel.
Multiplication by `exp(-sigma(y))` does **not** preserve that kernel as an
invariant subspace. The record stores the projected scale
`s = <psi_0, exp(-sigma) psi_0>` and the relative leakage
`||(exp(-sigma) - s) psi_0|| / ||exp(-sigma) psi_0||`. A finite Kronecker
2D matrix is a diagnostic of that projection, not a compact gap theorem.

## Stage 3/4 status and limitations

Stage 3 asked for one curved operator including warp, lapse, and shift. This
record supplies the first-order radial throat operator and a controlled
spatial warp coupling. It does **not** assemble the Lorentzian Hamiltonian.

Stage 4 asked for the matrix anomaly of one regulated measure. The maps and
resolvent here are regulator-independent finite responses of `D`. They do not
replace that variation.

Also not claimed: recursive Calderon/DtN tail, BFK gluing, electron mass,
`Xi`, child vacuum matching, or a healthy graviton.

## Reproduce

```bash
python3 scripts/check_nsc_boundary_response.py --check
python3 -m unittest tests.test_nsc_boundary -v
```

`--check` recomputes the record and compares every numeric field, exact flag,
convention string, and source hash. It is the full numerical reproduction.
The unit tests stay on small grids and do not rebuild the published record.
`--output PATH` writes a new file and refuses if `PATH` already exists.
Historical published JSON is never opened for writing. The previous
uncommitted v1 candidate was copied to `/tmp/nsc-3-boundary-response-previous.json`
before regenerating v2.

## Next equation

The same first-order maps must be inserted into the energy-resolved recursive
tail with the explicit `1/Omega` convention, and the matrix anomaly of one
regulated measure must be derived, before any stationary `zeta` is discussed.
