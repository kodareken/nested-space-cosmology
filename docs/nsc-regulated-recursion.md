# Regulated finite operator, recursion and the spacetime boundary

The [calculation](../scripts/check_nsc_regulated_recursion.py) and
[record](../results/nsc-3-regulated-recursion.json) implement one finite
proper-time prescription, its full matrix variations and a normalized nested
chain. They expose an additional obstacle to treating the old spatial scale
search as the stationary point of a covariant spacetime action.

## Explicit finite prescription

For a Hermitian finite Dirac matrix in an orthonormal basis, with no zero modes,

\[
\Gamma_{\Lambda,M}(D)=\tfrac12\operatorname{Tr}E_1(D^2/\Lambda^2)
+N\left[\log(M/\Lambda)+\gamma_E/2\right].
\]

Here Lambda is the proper-time cutoff, M the determinant normalization and N
the finite matrix dimension. This declared finite normalization has the limit
`-log|det(D/M)|` as Lambda increases with D and M fixed. It is not a derived
continuum gravitational measure. Zero modes require an explicit reduced
determinant and their own measure; the implementation rejects unresolved zeros.
The unwrapped eigenvalue-log phase is retained separately and can change only
through a zero crossing in this Hermitian family.

For `f(d)=E1(d²/Lambda²)/2`,

\[
f'(d)=-e^{-d^2/\Lambda^2}/d,\qquad
f''(d)=e^{-d^2/\Lambda^2}(2/\Lambda^2+1/d^2).
\]

The gradient is `Tr f'(D) D1`. In the eigenbasis of D, the Hessian includes
the full divided difference `(f'(d_i)-f'(d_j))/(d_i-d_j)` times
`|D1_ij|²`, with the coincident-eigenvalue limit `f''`, plus `Tr f'(D) D2`.
This retains noncommuting perturbations. Centered differences and unitary
basis changes independently check the implementation, including degeneracy.

## Derive the finite regulator defect rather than discarding it

Set `D_t=exp(-t Sigma/2) D exp(-t Sigma/2)` with Hermitian Sigma. The exact
finite raw determinant changes by `t Tr Sigma` in the negative logarithm.
The regulated transformation differs by

\[
C(t)=\Gamma_{\rm reg}(D_t)-\Gamma_{\rm reg}(D)-t\operatorname{Tr}\Sigma,
\qquad C'(t)=\operatorname{Tr}\Sigma(e^{-D_t^2/\Lambda^2}-I).
\]

The endpoint expression agrees with independent integration of its derivative.
For the traceless two-sheet rescaling the raw determinant is invariant, while
the regulated Hessian is `4p² exp[-(p²+Phi²)/Lambda²]/Lambda²`, as in the
September 7 correction. A compensating finite defect would cancel that change
along this orbit. This does not decide whether a physical relative-sheet field
is gauge, or derive its continuum compensator dynamics. The full metric and
field measure remains open. The cutoff derivative at fixed dimensional D and
M is `Tr exp(-D²/Lambda²)-N`; normalization must remain explicit.

The distinction between normalization and cutoff is consistent with the
primary spectral-anomaly setup, whose continuum assumptions must still be
checked for this geometry. [Andrianov, Kurkov and Lizzi](https://arxiv.org/html/1106.3263v1)

## A spatial heat trace omits a frequency factor

For the explicitly restricted ultrastatic Euclidean product
`D_E²=-partial_tau²+H²`, integrating frequency gives

\[
K_E(t)/T=\frac{\operatorname{Tr}_{\rm spatial}e^{-tH^2}}{\sqrt{4\pi t}},
\]

\[
\Gamma_E/T=\sum_j\left[
\frac{\Lambda e^{-e_j^2/\Lambda^2}}{2\sqrt\pi}
-\frac{|e_j|}{2}\operatorname{erfc}(|e_j|/\Lambda)
\right].
\]

Direct proper-time quadrature agrees with the closed expression. This is the
unnormalized modulus per time under the stated product assumption. The
additional frequency factor changes the variational functional; it is not
optional. A four-spatial-dimensional slice produces a five-dimensional product
heat trace here. It cannot silently inherit a four-spacetime-dimensional anomaly.

## The actual coordinate-time Hamiltonian is not elliptic in the trapped region

The recorded horizon-penetrating metric has
`A(rho)=1+3rho+3(1+rho²)(atan(rho)-pi/2)` and `beta=sqrt(1-A)`.
Up to an overall sign/orientation, the Hamiltonian's principal symbol has
eigenvalues

\[
-\beta k_\rho\pm\sqrt{k_\rho^2+k_\perp^2}.
\]

The common conformal warp cancels between lapse and inverse spatial tetrad at
this order. At the throat, `beta²=3pi/2>1`. The nonzero spatial covector
`k_rho=1`, `k_perp=sqrt(3pi/2-1)` makes one eigenvalue exactly zero.
The determinant identity `(beta²-1)k_rho²-k_perp²=0` is checked symbolically.

Consequently the standard elliptic spatial heat-operator argument does not
apply to this coordinate-time Hamiltonian in the trapped region. This is not
an instability or ill-posedness proof for covariant Lorentzian Dirac evolution,
nor a proof that every regulated trace diverges. It establishes the failure of
that proposed route to a spacetime determinant. The next physical functional
needs a state-defined Lorentzian construction or a justified continuation;
simply adding lapse/shift to the old positive spatial operator is insufficient.

This phenomenon has prior literature: Finster and Röken construct
self-adjoint Dirac Hamiltonians in settings where horizons destroy
ellipticity. Their boundary and asymptotic assumptions are not automatically
satisfied by this carrier. Our contribution here is the explicit witness and
its implication for this project's determinant route, not discovery of the
general phenomenon. [Finster and Röken, arXiv:1512.00761](https://arxiv.org/abs/1512.00761)

## Normalized finite recursion

For `H_n=Omega^n H` and `B_n=Omega^n b` in one common dimensional energy frame,

\[
\Gamma_n(x)=xI-H-\Omega^{-1}b\Gamma_{n+1}(x/\Omega)^{-1}b^\dagger.
\]

Depths 2, 4, 8, 16 and 32 at three Omega values agree with direct full-chain
inversion. Noncommuting complex links test ordering, and the inverse response
has positive imaginary part in the upper half-plane. Depth changes are
reported as observations, not rigorous infinite-tail error bounds. H and b
in this control are declared matrices, not yet the computed throat maps.

Changing D, Lambda and M together while changing L_star inversely leaves the
action and zeta invariant. No unit convention can select Omega, and the code
does not identify zeta with Omega². A stationary geometry, physical boundary
field, internal Dirac spectrum and retarded metric state are still required.

There is also an exact finite stationarity obstruction. Along the physical
matrix dilation `D(t)=exp(-t)D` at fixed cutoff, normalization and domain,
`dGamma/dt=Tr exp(-D(t)²/Lambda²)>0`. No choice of finite Hermitian links can
create an isolated common-scale extremum in that family. The runner compares
the variation with the independent heat trace. Scaling D, cutoff and
normalization together instead gives a flat change of units. Neither case
selects an inherited physical scale. This does not cover the unknown full
continuum field measure, constrained geometry equations or nonuniform
physical deformation. Those equations, rather than a different scalar root
finder, are the next required layer.

## Reproduction and current conclusion

`python3 scripts/check_nsc_regulated_recursion.py --check` compares every exact
and scientific numerical field. This calculation completes the finite
variation and normalized-chain controls. It does not solve the common curved
field equations, identify particles, validate a metric Hessian, or predict
child expansion. The new project-level finding is the explicit spacetime
normalization/domain obstruction, not a new law of nature or priority claim.
