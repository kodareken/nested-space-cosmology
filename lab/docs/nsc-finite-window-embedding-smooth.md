# Seam-regular compression of the same finite window

The target, gauge and operator are unchanged from the v1 compression. `H`, `B`
and `Ω = 3/2` still come from `finite_window`, and the slice is still
`Q = b0/a0`, `L = b0 Ω^x`, `β = β0 Ω^x`, `κ = 1` on the packet arc `[0, 4]`
inside the circle of length 8. The restrictive real-lobe and conjugate-column
obstructions stay in the v1 record. They are not repeated here.

## Representative

Each region `n = 0, 1, 2` carries a plus mode and a minus mode, supported on
`(n, n+2)`. A lobe on a unit interval is `[u(1-u)]^4` times six Legendre
polynomials `P_k(2u-1)`. The plus right lobe has one constant phase
`0.676570503790893`. The minus right-lobe phase is `0`. The column phase on
the minus carrier is `0.6147015926934642`. Carriers remain
`exp(+i(x-n))` and `exp(iφ) exp(-i(x-n))` in the `σ2` frame
`[1, i]/√2` and `[1, -i]/√2`.

The flat factor removes the join derivative jumps of the degree-4 sine lobes.
The coefficients were obtained by enforcing the Gram, onsite and link targets
on this basis. They are not a smoothing of the saved sine columns. Counts of
4 and 5 Legendre terms, continued from a solved flat-power-3 lobe and checked
with random starts, stalled near `1.7e-2` and `4.3e-3`. A 1200-evaluation
continuation at count 5 stayed at `4.349e-3`. Count 6 reaches a continuum
residual below `1e-14`. That is the smallest count this search found. It is
not a certificate that every other 5-coefficient chart is empty.

The Legendre coefficients are large, of order `10^3` to `10^4`, because they
multiply a flat factor whose maximum is `1/256`. The scalar profiles are
normalized in the continuum `L2` inner product. `half_density_columns` multiplies
by `√dx` and the same spinors as `apply_dirac`.

## Seams

For flat power 4 the lobe and its first three derivatives vanish at each
endpoint. The assembled carrier profile therefore matches the zero extension
through third order at `x = n`, `n+1` and `n+2`. A Cauchy jet of each lobe,
radius `0.05`, has maximum size `9.984757302269112e-12` through order 3. The
sampled profile and its classical derivative are `0` at the seams. At distance
`10^{-3}` the largest derivative is `0.0001412652295117061`, and shrinking
that distance by 10 shrinks the derivative by at least `842.6045862320477`.
A fourth-order finite difference inside the lobes agrees with the returned
derivative to `6.102929128632061e-08`.

## Continuum match

Direct quadrature of `V† H V`, with no inserted factors of `Ω`, gives:

| Quadrature per unit panel | Matrix error |
|---|---:|
| 64 | `7.721884926240916e-14` |
| 128 | `2.7541733819492414e-14` |
| 256 | `5.559950713426307e-14` |

At 256-point panels the Gram defect is `8.881785873244357e-15`, the Hermiticity
defect is `2.2221408060595165e-14`, and the corner is `0`. The gap between the
128-point and 256-point matrices is `6.1114139796286e-14`. The minimum
complement fraction is `0.9934510425652145`. The projected sign-control current
of this compression, for the coefficient vector `[1, 0, i, 0, 0, 0]/√2`, is
`0.2500000000000039`. That current is the regional trace of `V† H V` at that
instant. It is not the current of the full Dirac evolution.

The second derivative of subspace population for the same vector is
`-387.10353461425734`, equal to `-2‖(1-P)Hψ‖²`. About `99.3%` of `‖Hψ‖` lies
in the complement. That complement is retained outside memory. It was not
minimized, and these six modes are not an invariant subspace: they still
vanish together on `(4, 8)`. The closed evolution `e^{-iJt}` is not the Dirac
solution.

## Fourier check

`discrete_compression` samples the same columns and calls the existing
`apply_dirac`. It is not the fitting objective. The v1 sine lobes, whose
derivatives jump, had compression errors `0.21315073261549564` at `N = 128`
and `0.053042488066192287` at `N = 256`. This representative measures:

| Grid | Matrix error | Gram defect | Bridge image max |
|---:|---:|---:|---:|
| 128 | `0.010838594026637356` | `9.336700162476674e-05` | `0.32869759365962326` |
| 256 | `7.746956440221127e-05` | `4.0662163769411563e-07` | `0.04186908812033639` |
| 512 | `3.24320182936768e-07` | `5.122049451244948e-10` | `0.003623808036605499` |

The `N = 256` error is `238.8675404062514` times the `N = 512` error. Sampled
mass outside `[0, 4]` is `0`. Complement fractions stay near `0.993`.
`N = 1024` was not required.

`column_profile` and `half_density_columns` are the calls a later coupling can
use. No coupled evolution is computed here. Nothing in this successor changes
the action, the Dirac operator, or the v1 files.
