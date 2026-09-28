# Finite fixed-C0 preparation fixes the incoming geometric germ

Within the smooth finite-history class stated below, exact prepared
`C_g(Sigma)=C0` requires

\[
\boxed{\partial_T^k a_g(z)=\partial_T^k a_{\rm ref}(z),\qquad
\partial_T^k r_g(z)=\partial_T^k r_{\rm ref}(z),\quad k=1,2,3,4.}
\]

Their spatial derivatives agree too, so all twenty permitted normal/mixed
jet differences vanish. This is a finite-amplitude necessary condition;
it neither fixes higher jets nor makes those conditions sufficient for C0
matching. The physical-state remainder and finite endpoint induction were
reviewed independently through Grok and integrated with the existing
operator and source owners.

## Physical connection

For a smooth finite metric history in the class below, let `C_g(T)` be the
covariance obtained from the original upstream state by the exact canonical
Dirac evolution. A finite local projector approximation `P_g(T)` obeys

\[
C_g(T)-P_g(T):H^s(\mathbb R;\mathbb C^2)
\longrightarrow H^{s+6}(\mathbb R;\mathbb C^2).
\]

This operator regularity statement is different from a pointwise numerical
bound on a truncated symbol. It permits extracting the local coefficients
of degrees minus two through minus five from the physical covariance.
Neither the incoming state nor the action subtraction is replaced by P_g.

The unchanged source and fixed canonical frame are essential. Equality of
two covariance operators is stronger than matching finitely many sampled
energies, diagonal pairs, integrated moments or symbol traces.

## History class and actual operator

Use the existing fixed KS chart on a finite slab between an unchanged
upstream slice and incoming Sigma. The slab stays inside the regular
trapped chart and away from the horizon and physical rho0 seam. Metric
functions are real and smooth, agree with the reference near the upstream
slice, and differ from it only in a compact axial region. N,a,r remain
positive. On a compact time interval these assumptions give uniform lower
bounds and bounded derivatives; amplitudes need not be infinitesimal.
No numerical interval length or cosmological duration is selected.

The incoming slice is spacelike. In the owner's +--- convention,
`ds2=N^2 dT^2-a^2(dz+beta dT)^2` and `g^(TT)=1/N^2>0`.
At the reference incoming slice `T=T(rho)`, its radial normal has squared
norm `a0^2=beta_PG^2-1>0`. A fixed-radius surface outside the horizon has a
different causal character; that exterior description is not used here.
Where the old PG representation is also used, retain its existing connected
chart conditions `F>0,s_K>0` from the metric-pullback owner.

For one retained massive nonzero-angular block the existing symbol is

\[
h_g(T,z,k)=A_g(T,z)k+B_g(T,z),\quad
A_g=\frac Na\sigma_3-\beta I,\quad
B_g=-Nm\sigma_1+\frac{N\ell}{r}\sigma_2.
\]

Its Weyl quantization is the already owned symmetric first-order operator
`H_g=(A_g D_z+D_z A_g)/2+B_g`, with `D_z=-i partial_z`. The half-density
ordering is included once. Scalar shift may make an individual principal
velocity zero; the two bands still differ by `2(N/a)|k|` at nonzero k.
Thus inversion of the off-diagonal commutator needs the band difference,
not ellipticity or a chosen sign of the full shifted Hamiltonian.

The endpoint comparison keeps intrinsic N,beta,a,r fixed as functions
along Sigma, and also fixes the N/beta normal jets through order four
there. The free incoming slots are exactly the twenty a/r normal/mixed
slots already owned by `IncomingNormalJetChange`.

The existing signed-real-energy spectrum is unbounded. The necessary
conditions use arbitrarily large |k|; no finite physical energy cutoff is
introduced or covered by the theorem. The bounded-frequency partitions
used below are mathematical devices for estimating the same complete state.

## Imported method and finite-order use

[Gérard and Stoskopf, Lemma6.4 and Proposition6.8](https://arxiv.org/html/2108.11630)
provide Sobolev propagation and the time-dependent approximate-projector
method. Their full Hadamard conclusion is not an assumption about this
state. For the present shifted operator, use the separated-band Weyl
recursion rather than their spectral sign selection.

Smooth bounded Hermitian A_g,B_g give a unitary evolution U_g on L2 and
bounded evolution on each Hs on the finite slab. The Sobolev estimate uses
the scalar multiplier `Lambda^s=(1+D_z^2)^(s/2)` and bounded order-zero
commutator `[Lambda^s,H_g]Lambda^(-s)`. This does not require a nonzero
individual propagation velocity.

The existing projector construction solves idempotence for the diagonal
coefficient and transport for the off-diagonal coefficient. In the
classical large-|k| expansion it can be retained through degree minus six,
using the nonzero band difference. Smoothly joining the two large-momentum
cones across bounded k affects only a regularizing operator. The
finite approximation has

\[
R_g=\partial_T P_g+i[H_g,P_g]\in\Psi^{-6}.
\]

The old derivative-order-four numerical SymbolJet is not executed at order
six or treated as supplying uncomputed derivatives. This construction is
the finite classical-symbol extension of its imported recursion. Its
coefficients and remainder use the smooth history, not a bare point germ.

For the finite extension, put `G=p#p-p` and
`F=i partial_T p-[h,p]_#`. Associativity and the Leibniz rule give
`[p,G]_#=0` and
`i partial_T G-[h,G]_#=F#p+p#F-F`. At each next classical order these
identities supply the same diagonal/off-diagonal compatibility used by the
owned construction. The next diagonal correction cancels G; the next
off-diagonal correction cancels F by division by the nonzero band difference.
Continue only to degree minus six. Then idempotence has remainder order
minus seven and transport has remainder order minus six. This is a finite
symbol existence argument, not an evaluation of higher stored SymbolJets
or a convergence claim for an infinite series.

## Actual upstream state and both momentum cones

The saved physical affine-vacuum bound controls P16 by a positive
full-collar defect integral. At any unchanged upstream slice inside that
collar it yields an O(E^-16) projector error for positive source energy.
The upstream geometry and covariance are translation invariant, so the
error is a Fourier multiplier. Its norm decay alone implies an Hs to
Hs+16 bound; differentiability estimates in energy are not required for
this particular remainder.

The `per_sign` rows in that receipt are **angular signs**. The other energy
cone follows from the existing source and mode relation

\[
C(-E,\ell)=I-\sigma_3\overline{C(E,-\ell)}\sigma_3.
\]

The corresponding complementary conjugation carries the projector error
with the same norm. Both angular receipts are therefore used to cover
both energy cones of a fixed angular block. The physical convention is
`k=-E` in the incoming axial Fourier representation.

For positive energy, the exact source remainder relative to its affine
vacuum is bounded by `max(sqrt(f_H),n_in)` after canonical coisometry.
Its exponential decay, and the signed relation above, retain the complete
thermal and coherent correction as a regularizing multiplier. On bounded
frequency sets the covariance is bounded by the canonical fermionic
positivity bounds, which is sufficient for every Sobolev gain after a
compact frequency cutoff. No low-energy quadrature accuracy is inferred.

The background P16 projector and the finite local construction solve the
same idempotence/transport recursion and have the same band. Uniqueness of
each successive diagonal/off-diagonal coefficient gives agreement through
degree minus six near the unchanged upstream slice. Their difference is
then at least order minus seven, sufficient for the required gain of six.

## Propagating the remainder without assuming the answer

Let `T_u` be the unchanged upstream slice and
`D_u=C_up-P_g(T_u)`. Exact evolution and the defect identity give

\[
\begin{aligned}
C_g(T)-P_g(T)
={}&U_g(T,T_u)D_uU_g(T,T_u)^\dagger\\
&-\int_{T_u}^{T}U_g(T,t)R_g(t)U_g(T,t)^\dagger\,dt.
\end{aligned}
\]

The upstream difference and R_g both map Hs to Hs+6. Sobolev bounds for
U_g and its inverse preserve that mapping property on both sides of each
term, proving the stated finite regularity gain. This uses the physical
source; it does not promote the formal action reference into a vacuum.

If the exactly prepared endpoint covariance equals the reference C0,
subtracting the two remainder relations implies that
`P_g(Sigma)-P_ref(Sigma)` also improves six Sobolev orders. A nonzero local
classical coefficient of degree minus j, for j<=5, contradicts that gain:
apply the operator to a localized high-frequency spinor with nonzero
principal coefficient. Therefore those endpoint coefficients agree.

## Finite endpoint induction

On the cone `k=-E`, `E>0`, write the classical symbol as
`p=p_0+E^(-1)p_1+...`, with `p_0=diag(1,0)` and
`h=E A_-+B`, `A_-=-(N/a)sigma3+beta I`. Its first off-diagonal coefficient
is

\[
(p_1)_{01}=\tfrac a2(m+i\ell/r),\qquad
(p_1)_{10}=\tfrac a2(m-i\ell/r).
\]

The triangular structure is as follows. Coefficient p_n uses
normal derivatives through order n-1. Its highest normal derivative comes
only from the transport term `i partial_T p_(n-1)` followed by inversion
of `[A_-, .]`. Idempotence determines the diagonal from lower coefficients.
Spatial Weyl terms differentiate in z and k, and therefore supply no
additional time derivative. Nonlinear products and differentiated inverse
gaps that do not contain this highest derivative involve lower normal jets.

Consequently, once the two histories agree in lower normal jets as functions
along Sigma, their first differing degree-(k+1) coefficient is

\[
\begin{aligned}
(\Delta p_{k+1})_{01}
&=-\frac{a_1^k}{(2i)^{k+1}}
\left[(\ell/r_1-im)\Delta a_{T^k}
-\ell a_1\Delta r_{T^k}/r_1^2\right],\\
(\Delta p_{k+1})_{10}
&=\frac{(-1)^ka_1^k}{(2i)^{k+1}}
\left[(\ell/r_1+im)\Delta a_{T^k}
-\ell a_1\Delta r_{T^k}/r_1^2\right].
\end{aligned}
\]

These are finite differences, not derivatives with respect to an amplitude
parameter. The displayed values use the owned endpoint N1=1 and beta1=0.
The linear appearance of the first unmatched derivative follows from
highest-derivative structure, not from assuming that the whole nonlinear
projector is linear. A nonlinear term containing lower normal derivatives
can be nonzero in either history while canceling in their difference.

The corresponding determinant is
`i*m*ell*a1^(2*k+1)/(2^(2*k+1)*r1^2)`, the same nonzero map already checked
for the tangent coefficients. With the physical remainder bridge established,
equality of physical covariances supplies equality of these coefficients;
induction k=1..4 then fixes the whole permitted order-four germ. Its mixed
spatial derivatives agree only after equality of the normal functions is
established along Sigma.

Grok endpoint review69ac44c5-8b55-41bf-ba3c-200822a36b06 independently
confirmed the finite-difference formulas, derivative count, common lower-jet
cancellations and physical sign convention. The committed module checks33
exact scalar identities and its three focused tests pass. In particular,
the test with common nonzero first derivatives keeps the nonlinear lower-jet
terms before taking a finite second-jet difference.

The corrected Grok physical reviewa4589cc3-11b0-4c92-9dea-a7b98d37649b
accepted the finite symbol construction, upstream P16 identification,
Sobolev estimate and exact prepared-state Duhamel connection. Its initial
timelike-slice objection was withdrawn after checking the actual trapped
metric. The physical and formal arguments together establish the boxed
necessary condition.

## Consequence for the included incoming constraints

The [included assembly](nsc-incoming-joint-constraints.md) depends on the
order-four metric germ for its local and reference contributions, and on
the same instantaneous intrinsic geometry and C0 for its physical matter
insertion. Every one of those arguments now equals its reference value.
This is equality of the finite coefficients, not only zero first variation.
Consequently

\[
\mathcal E_\beta[g,C0]=\mathcal E_\beta[g_{\rm ref},C0]
=P_K<-M,\qquad
M=\frac{5887634380411363}{265464355000000000}>0.
\]

The exact raw-current witness is reused from the
[incoming momentum balance](nsc-incoming-surface-momentum-balance.md).
It is independent of the unresolved numerical accuracy of the lapse
source. Thus no member of this finite compact-prepared, fixed-C0 class
satisfies both included incoming constraints, even on an open interval.
The conditional local surface solutions remain valid as geometric data;
they are not realizable by these histories with this exact fixed C0.

This class has unchanged upstream preparation and source law, the same
intrinsic slice and canonical identification, and fixed N/beta endpoint
jets through four. The result does not exclude different parent data,
other preparation classes, an evolved incoming covariance, or the full
nested-space architecture. It gives no global injectivity theorem for
histories and no extended transmitting stationarity verdict. No such
alternative is adopted or fitted here.

## Verification and reuse

The module's33 scalar identities have exact residual zero; three focused
tests pass, including a finite second-jet change with common nonzero lower
jets and the actual signed source on open/closed fibers. These checks
verify finite algebra. The Sobolev, symbol-order and derivative-count
arguments above are analytic, not numerical error enclosures. There is no
claimed finite-energy threshold or computed norm constant.

The P16, source, tangent, constraint and operator evidence is reused with
its original hashes and scope. No field propagation, reference/source
integration, horizon calculation, new action term, physical initial tuple,
metric timestep or publication is part of this result.
