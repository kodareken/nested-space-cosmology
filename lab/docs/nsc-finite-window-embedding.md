# Finite window as a compression of the existing Dirac operator

The target is the frozen depth-3 window from `finite_window`, with the recorded
`H = [[1, i/5], [-i/5, 2]]`, `B = [[1/4, i/7], [1/9, 1/6]]` and `Ω = 3/2`.
The operator is the continuum form of `apply_dirac` on the calibrated slice
`Q = b0/a0`, `L = b0 Ω^x`, `β = β0 Ω^x`, `κ = 1`, on the packet arc `[0, 4]`
inside the circle of length 8. No entry of `B` is copied in, and the
inherited-law check is the existing `exact_controls` identity.

## Mode choice and support

Two modes sit in each region `n = 0, 1, 2`. Each mode is supported on
`(n, n+2)` and is assembled from degree-4 sine lobes on `(n, n+1)` and
`(n+1, n+2)`. The plus carrier is `exp(+i(x-n))` with spinor `[1, i]/√2`.
The minus carrier is `exp(iφ) exp(-i(x-n))` with spinor `[1, -i]/√2`,
`φ = 0.5298187870937372`, and an envelope that is not the conjugate of the
plus envelope. The plus right lobe carries one constant phase
`0.6332154732379276`. The minus right-lobe phase is `0`. The six profiles
vanish on `(4, 8)`.

The sine series is zero at each lobe endpoint, so the profile is continuous
and the right-lobe phase does not open a jump in the value. The derivative
may jump at `x = n+1` and at the outer endpoints, where a nonzero sine slope
meets the zero extension. Continuum matrix elements use that piecewise
classical derivative. Agreement of those integrals with `H` and `B` does not
mean the antiperiodic Fourier derivative at `N = 512` evolves the same
operator. The Fourier symbol sees the kinks. On these modes its compression
error is `0.21315073261549564` at `N = 128` and `0.053042488066192287` at
`N = 256`.

## Restrictive packet, then one phase

`prepare_rank6` at `N = 128` matches the onsite blocks
(maximum Frobenius error `2.5346486206910493e-6`) and misses the links
(`5.9845119231619615` and `8.97676788474295`). Its link obeys
`B_{-+} = conjugate(B_{+-})`. The frozen link does not:
`|1/9 + i/7| = 0.18098022620621237`.

Any real left and right lobes at the recorded carrier `k = 1`, in this spin
frame, with neighbor `L2` orthogonality, need weighted overlap
`0.2071801534525706` to make `B++ = 1/4`. Popoviciu's bound on
`ρ = Ω^{1+t} ∈ (Ω, Ω²)` caps that overlap by `(Ω² − Ω)/4 = 0.1875`.
The gap is `0.019680153452570598`. This obstruction stops at real lobes,
this carrier, this spin frame and this one entry. The constant plus-lobe
phase is the freedom that leaves the obstruction.

## Compression, leakage, and what is kept

Direct quadrature on 256-point panels, with no hand-inserted `Ω` factors,
gives `V† H V` against the frozen window. The saved record holds these values:

| Quantity | Value |
|---|---:|
| Matrix error | `4.450601090179064e-14` |
| Gram defect | `1.7763578458222943e-15` |
| Hermiticity defect | `4.8827941617840446e-14` |
| Corner Frobenius error | `0.0` |
| Minimum complement fraction | `0.9901480649567598` |

The sign-control current of that compression is `0.2499999999999967`. For the
same initial vector, the second derivative of subspace population is
`-306.29864965687835`. A state prepared in these six modes therefore has the
stage-1 block currents at that instant, including the regional sign control.
Its first population derivative in the subspace vanishes because `H` is
Hermitian. The state does not stay in the window: about `99%` of `‖Hψ‖`
lies in the complement, and the modes vanish together on an open arc, so
they are not an invariant subspace. The closed evolution `e^{-iJt}` is not
the Dirac solution. Matching the blocks is stronger than matching the
spectrum, and it still does not retain the finite-time stage-1 state.

The saved record repeats this quadrature at 256-point panels. Nothing here
changes the action, the Dirac operator, or any older result file.
