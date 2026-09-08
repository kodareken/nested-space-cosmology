# The physical radial Dirac operator and its throat current

The [exact runner](../scripts/check_nsc_dirac_tetrad.py) and
[record](../results/nsc-4-dirac-tetrad.json) derive the massless Lorentzian
Dirac equation on the existing **fixed, unwarped four-dimensional**
black-universe background. The spin connection gives

\[
\boxed{H_\kappa=-i(\sigma_2-\beta I)\partial_\rho
+\frac{i}{2}\beta'I+\frac{\kappa}{r}\sigma_1,\qquad u=r\psi.}
\]

This equation supports a well-defined whole-line fermion evolution. Its
physical current also establishes that the earlier reflecting spatial
half-domain construction is not the Lorentzian throat boundary problem.
The result does not derive the background, its gravitational action, a quantum
state or a physical metric Hessian.

## Tetrad and spin connection

Use signature `+---`, units `c=hbar=L_throat=1`, and

\[
ds^2=d\tau^2-(d\rho+\beta d\tau)^2-r^2d\Omega_2^2,
\quad r=\sqrt{1+\rho^2},
\]
\[
A=1+3\rho+3(1+\rho^2)(\arctan\rho-\pi/2),\qquad
\beta=\sqrt{1-A}.
\]

The orthonormal coframe is

\[
e^0=d\tau,\quad e^1=d\rho+\beta d\tau,\quad
e^2=r\,d\theta,\quad e^3=r\sin\theta\,d\phi.
\]

Its dual is `E_0=partial_tau-beta partial_rho`,
`E_1=partial_rho`, `E_2=partial_theta/r`, and
`E_3=partial_phi/(r sin(theta))`. With `q=r'/r`, the nonredundant
**lowered-index** connection forms are

\[
\omega_{01}=-\beta'e^1,\quad
\omega_{02}=-\beta q e^2,\quad
\omega_{03}=-\beta q e^3,
\]
\[
\omega_{12}=q e^2,\quad\omega_{13}=q e^3,\quad
\omega_{23}=\frac{\cot\theta}{r}e^3.
\]

The runner computes exterior derivatives from the coordinate coframe and
checks every component of `de^a+omega^a_b wedge e^b=0`, raising the first
connection index with the Lorentz metric. It also checks metric
compatibility, reconstruction of the declared metric and all Clifford
anticommutators using explicit four-by-four gamma matrices.

For `nabla=d+omega_ab gamma^a gamma^b/4`, direct matrix contraction gives

\[
\gamma^a\Omega_a
=-\gamma^0\left(\frac{\beta'}2+\beta\frac{r'}r\right)
+\gamma^1\frac{r'}r+\gamma^2\frac{\cot\theta}{2r}.
\]

The covariant equation is consequently

\[
i\left[
\gamma^0\left(\partial_\tau-\beta\partial_\rho
-\frac{\beta'}2-\beta\frac{r'}r\right)
+\gamma^1\left(\partial_\rho+\frac{r'}r\right)
+\frac{\gamma^2}{r}\left(\partial_\theta+\frac{\cot\theta}2\right)
+\frac{\gamma^3}{r\sin\theta}\partial_\phi
\right]\psi=0.
\]

The normal-slice inner product is
`integral r² psi_dagger psi d_rho d_Omega`. Rescaling `u=r psi` gives flat
radial measure and cancels both radial measure terms, including
`beta r'/r`. The angular spin connection remains part of the unit-sphere
Dirac operator. On its `kappa=plus_or_minus(ell+1)` sector, choose the
constant radial basis in which the zero-shift operator is
`-i sigma2 partial_rho+(kappa/r)sigma1`. This yields the boxed Hamiltonian.
Its shift term is derived from the tetrad, not guessed from Hermiticity.

An independent known control is `r=rho>0`, `beta=sqrt(2M/rho)`.
Before rescaling, the shift becomes `i beta(partial_rho+3/(4rho))`, as in
[Dolan, Doran and Lasenby, section IV.2](https://arxiv.org/html/gr-qc/0605031).
The angular sign convention may change the label `kappa`.

## Current, the throat, and finite boundaries

Set `C=sigma2-beta I` and `V=(kappa/r)sigma1`. The operator is
`H=-i(C partial_rho+C'/2)+V`. Its conserved density and radial flux are

\[
n=u^\dagger u,\qquad
j=u^\dagger C u
=2\operatorname{Im}(u_1^*u_2)-\beta n,
\qquad\partial_\tau n+\partial_\rho j=0.
\]

The exact check uses generic complex amplitudes and derivatives. Removing
the `i beta'/2` term gives the nonzero false source `-beta' n`; this is an
explicit negative control. The Green identity is

\[
\langle u,Hv\rangle-\langle Hu,v\rangle
=-i[u^\dagger Cv]_{\rho_L}^{\rho_R}.
\]

The characteristic speeds in the `sigma2` eigenbasis are
`a_plus=1-beta` and `a_minus=-1-beta`. The actual horizon is
`rho_h=1.90069160547`. At the throat
`beta=sqrt(3pi/2)=2.17080376367`, so both speeds are negative:
`-1.17080376367` and `-3.17080376367`.

For an interval with a trapped left cut and an exterior right cut:

| Endpoint | Incoming fields for forward evolution | Outgoing fields |
|---|---:|---:|
| Left, `beta>1` | 0 | 2 |
| Right, `0<beta<1` | 1 | 1 |

Thus a forward evolution prescribes only the incoming right characteristic.
The left cut is an outflow surface. Setting a radial component to zero
there changes or overconstrains the physical evolution.

At an exterior endpoint a reflecting line can be written in characteristic
variables as

\[
v_-=e^{i\vartheta}\sqrt{\frac{1-\beta}{1+\beta}}\,v_+.
\]

The old condition `u1=0` has current `-beta |u2|²`, so it is not reflecting
for nonzero shift. At a trapped cut `C` is negative definite and has no
nonzero flux-null line. The endpoint form
`C_right direct_sum (-C_left)` has inertia `(3 positive, 1 negative)`.
No endpoint-only self-adjoint finite-interval realization can have those
unbalanced signs. The exact runner checks the throat matrix and an exterior
representative `beta_right=1/2`; all exterior values in `(0,1)` have the
same signs. That representative is not substituted for the actual numeric
geometry probes.

At a smooth throat, continuity `u_child=u_parent` cancels the opposite-normal
fluxes. A generalized transmission matrix must obey
`U_dagger C(0) U=C(0)`, rather than the intrinsic spatial relation involving
`sigma2` alone. A standalone parent half-domain loses flux into the child;
the child half-domain requires both incoming throat components. The
[spatial boundary maps](nsc-boundary-response.md) remain valid within their
stated scope and do not become Lorentzian retarded maps by renaming them.

## Whole-line evolution and its precise scope

There is a direct self-adjointness argument for each fixed finite `kappa`
on the initial core `C_c^infinity(R; C²)`. Write
`x=pi/2-atan(rho)`, which lies in `(0,pi)`. Then

\[
\beta^2=3(1+\rho^2)x-3\rho
\le(3\pi+3/2)(1+\rho^2),
\]

using `2|rho|<=1+rho²`. Thus both smooth characteristic vector fields
`a_plus partial_rho` and `a_minus partial_rho` have at most linear growth
and complete flows. The horizon zero of `a_plus` is a stationary point;
it does not introduce a boundary reachable in finite flow time.

After a constant unitary diagonalization of `sigma2`, each uncoupled
transport generator is `T=-i(a partial_rho+a'/2)`. If `Phi(t,rho)` is its
complete flow, its unitary group is

\[
(U(t)f)(\rho)=
\sqrt{\partial_\rho\Phi(-t,\rho)}\,f(\Phi(-t,\rho)).
\]

The closure of the smooth compact-support transport generator is
self-adjoint. One may see the core property by changing to flow coordinates
on each interval where `a` is nonzero; the horizon is at infinite flow
distance. Compactly supported functions away from the zero are already a
core, and enlarging to all smooth compact-support functions preserves the
same closure. The potential is a bounded Hermitian perturbation,
`||V||<=|kappa|` because `r>=1`; consequently the full radial Hamiltonian
is essentially self-adjoint. This is a proof for the boundaryless whole
line, not for a reflecting finite box.

The argument is consistent with
[Chernoff's hyperbolic-generator framework](https://www.sciencedirect.com/science/article/pii/0022123673900037).
The horizon literature also makes clear that loss of Hamiltonian ellipticity
does not by itself rule out self-adjointness;
[Finster and Roken](https://arxiv.org/abs/1512.00761) provide a different
boundary theorem whose assumptions are not imported here.

Causal curves obey `|d_rho/d_tau+beta|<=1`. The same linear-growth bound
prevents radial escape at finite `tau`; angular sections are compact. On
the declared smooth boundaryless product
`R_tau x R_rho x S²`, each `tau` slice is therefore Cauchy. Norm-preserving
causal fermion evolution on this fixed background is justified. It neither
selects a vacuum nor establishes positivity of a physical metric sector.
The five-dimensional compact warp, constrained gravitational variation,
stationary geometry and observed physical scales remain separate work.

## Reproduction

Run `python3 scripts/check_nsc_dirac_tetrad.py --check`. The runner asserts
the exact matrix and polynomial identities, recomputes the horizon and
characteristic probes, and compares every recorded field including both
source hashes. It reuses the active comparator: floating fields have
`rtol=1e-8, atol=1e-8`; symbolic expressions, labels, flags and provenance
must agree exactly. All numerical fields must be finite.
`--output PATH` uses exclusive creation and cannot overwrite a result.

The complete-flow and bounded-perturbation reasoning is a mathematical
argument documented above, not a claim that symbolic algebra proves an
abstract domain theorem. The executable checks its concrete coefficients,
the bound's square identity and the boundary signs used in that argument.
