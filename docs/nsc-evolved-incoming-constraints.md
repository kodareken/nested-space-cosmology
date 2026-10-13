# Both incoming constraints with the evolved source state

The new assembly uses the declared source-fixed state law:

\[
\mathcal E_B[g]=\mathcal E_B[g_{\rm ref}]
+\Delta\mathcal E_B^{\rm local+reference}[j_\Sigma g]
+\bigl(G_B[C_\Sigma[g]]-G_B[C_{\rm ref}]\bigr),\quad B=N,\beta.
\]

The reference covariance is a fixed computational reference for this
difference. It is never imposed as the evolved physical state. Consequently
the derivative includes `delta G_B[C_Sigma[g]]` and does not subtract a
derivative of the reference history. Bare incoming C0 matrices and legacy
matter dictionaries are rejected by the new evolved-state interface.

## Actual momentum and complete source derivative

Let F denote evolved canonical source columns weighted once with the
original square roots of `dE/(2*pi)`, and let Csrc be the fixed full source
matrix. For the actual axial derivative F_z, define

\[
D=FC_{\rm src}F^\dagger,\qquad
J=\frac{F_zC_{\rm src}F^\dagger-FC_{\rm src}F_z^\dagger}{2i}.
\]

With the inherited signed-family multiplicity mu, the raw action
coefficients are

\[
G_N=\mu\left[m\operatorname{tr}(\sigma_1D)
-\frac\ell r\operatorname{tr}(\sigma_2D)
-\frac1a\operatorname{tr}(\sigma_3J)\right],\qquad
G_\beta=\mu\operatorname{tr}J.
\]

These are the existing
[near-diagonal KS pairing and vertices](../src/recursive_horizons/nsc_common_subtracted_ks_source.py),
with its Gaussian minus sign and symmetric momentum ordering. The new
implementation calls those raw vertex coefficients. No additional
`4*pi*r^2` factor or folded-energy factor is inserted. Here
`mu=copy_count*degeneracy/number_of_actual_angular_signs`; signed source
energies remain explicit.

The density derivative reuses the owned
[coherent restriction derivative](../src/recursive_horizons/nsc_prepared_history_jets.py).
The momentum derivative retains all four terms involving delta F and
delta F_z. The fixed source law gives delta Csrc=0 explicitly. The fixed
intrinsic identification makes the instantaneous vertex variation zero for
this family; normal-jet variations still enter the local/reference terms.

Only in the stationary reference does `F_z=-iE F`. That substitution is
wrong for evolved columns and their tangents. The complete source matrix
is retained, including its horizon coherences and closed source columns.

## Reuse the geometric response without a frozen matter insertion

For slots `(w,w_z,w_zz,w_zzz,U,U_z)`, reuse the certified coefficient functions

\[
\Delta\mathcal E_N^{\rm local+reference}
=A(w)U+B(w)w_{zz}+Cw_z^2+\sum_{j=1}^4D_jw^j,
\]
\[
\Delta\mathcal E_\beta^{\rm local+reference}
=dU_z+ew_{zzz}+F(w)w_z.
\]

The old constants S_N and S_beta are already in the baseline assembly and
are not inserted a second time. `D[0]` is explicitly zero. In contrast,
`F[0]=c_v` is a geometric response coefficient and is retained. The new
state correction is additional and evaluated through the actual history.
No fixed-C0 root is promoted by this decomposition.

The reference count is exact:
`-(Cref-Pref) + (Pg-Pref) - (Cg-Cref) = -(Cg-Pg)` in the common vertex
trace. Existing local allocations and compact-test reference-band scope
remain unchanged. Physical endpoint terms are not inferred from this
incoming bulk assembly.

## Verification and present accuracy

The `evaluate_source_fixed_history` entry point invokes the actual owned
evolution afresh for each supplied history, constructs fixed incident phases
and zero source/preparation tangents, then calls the new constraint assembly.
It accepts a fixed source preparation rather than an incoming covariance.
Channel and initial-reference consistency are checked, and a binding includes
the fixed source, initial columns, initial slice and channel. Cached states
are not substituted by this evaluation path.

Seven focused tests pass: stationary raw-vertex signs, the complete
coherent tangent, a finite-difference derivative, rejection of stationary
source-energy substitution on general columns, the existing physical
source-kernel contraction, and the geometric Jacobian with its constant
source terms removed. They also check fresh evolution delegation against
the actual recorded preparation inputs and the integrated receipt. These
tests contain no new source quadrature or field propagation. Grok
independently reviewed the signs and subtraction split.

The integrated numerical check consumes the new nonzero-history state
control and its centered changes. It uses the same original group14_1
source fibers and source weights. Axial derivatives estimated from saved
traces carry an explicit continuum-error hole; a discrete identity check
does not certify those derivatives or the full spectrum.

The [integrated record](../results/development/nsc-evolved-incoming-constraints.json)
has whole discrete constraint-derivative residual `1.352287829051212e-9`
against `3e-8`. Dropping delta C changes the derivative by
`0.0017487730205880414`; outgoing k=-E changes it by
`0.0017507437685358041`; dropping source coherence changes it by
`7.177503074526619e-6`. Each is separately resolved by this control.

The value error remains substantial. Reusing the saved zero-amplitude
801/64 trace gives the exact diagnostic decomposition
`Gdisc(g)-Gref = [Gdisc(g)-Gdisc(0)] + [Gdisc(0)-Gref]`.
The history-change maxima are `(1.6414751165e-6,1.7487615238e-6)` for N,beta,
while the numerical reference-field/derivative drift indicators are
`(0.00523815833,0.00441748961)`. Their sum residual is zero. The main
diagnostic is unchanged: this drift is reported, not deleted or used as a
physical compensating force. It identifies a concrete numerical value
error that must be resolved before a physical root is meaningful.

The returned quantity is therefore a *partial constraint diagnostic*, with
the evolved matter correction and its tangent exposed separately. Missing
changed-history energy ranges and signed families, baseline low/subgap
accuracy, and field/derivative errors keep physical constraint closure
OPEN. The baseline's certified high-energy regions are reused only for
the baseline; they are not relabeled as bounds on the changed state.

For the current pure-radius family, the old Dirac generator is independent
of r when ell=0. The source-fixed state correction of groups0,13,23 is
therefore exactly zero. Their baseline contributions and local allocations
remain present; only their response computation can be omitted.

No new force, source law, parameter fit, physical duration, metric timestep
or publication is introduced. The frozen-C0 theorem remains a regression.


## Conditional covariance and incoming-energy error bounds

For the same source-fixed state, auxiliary $P$ only decomposes
$C_g=P_g+D_g$. The physical difference retains both $P_g-P_{\rm ref}$
and $D_g-D_{\rm ref}$; the auxiliary projector supplies no new action.
On a computational Bloch cell of enclosing length $L>0$, use
$k_n=q+2\pi n/L$, fiber inner product $L^{-1}\int dz$, Weyl matrix
$A_{mn}(q)=\hat a_{m-n}((k_m+k_n)/2)$ and reconstruction measure
$dq/(2\pi)$ without another $1/L$. A physical local application requires
a separately proved causal embedding of that padded cell.

The uniform evaluation constants are

$$C_0^2=S^2\le1+L^2/4,\qquad C_1^2\le1/4+L^2/4.$$

For the two-spinor Hilbert--Schmidt norm
$H=\|\Lambda^2A\Lambda^2\|_{\rm HS}$, Schatten Hölder gives
$d_{11},d_{21}\le\sqrt2 HS$; point evaluation adds no further spin factor.
Parseval supplies

$$\int_{\rm BZ}\frac{dq}{2\pi}H^2\le\frac8L I[a],\qquad
I[a]=\int dz\int\frac{dk}{2\pi}
\left(\langle k\rangle^8\|a\|_F^2+
\frac1{256}\|\partial_z^4a\|_F^2\right),$$

so $\int_{\rm BZ}H\le\sqrt{8I}/L$. Keep $L,a,r$ as positive
enclosing intervals in denominators. Replacing them with upper endpoints
before inversion can underbound the error. With
$K_j=\int\|\partial_z^jB\|$, let
$M_1\le1+K_1$ and
$M_2\le1+\sqrt{4K_1^2+(K_2+K_1^2)^2}$. The actual defect
$R=\dot P-[G,P]$ then propagates as

$$d_{11}^{\rm ev}\le M_1^2\left(d_{11}(D_{\rm up})+\int d_{11}(R)\right),
\qquad
d_{21}^{\rm ev}\le M_2M_1\left(d_{21}(D_{\rm up})+\int d_{21}(R)\right).$$

Writing $\|\rho\|_1\le C_0^2d_{11}^{\rm ev}$ and
$\|j\|_1\le C_0C_1d_{21}^{\rm ev}$ gives

$$\epsilon_N\le\mu\left[
\sqrt{m^2+\ell^2/r^2}\,\|\rho\|_1+\|j\|_1/a\right],
\qquad \epsilon_\beta\le\mu\|j\|_1.$$

These are conditional propagators of supplied certified norms. Fourth axial
derivatives of the true fourth-order defect require mixed geometry order
9; spatial-only bounds and formal projector order do not supply them.
Instantaneous $I[a]$ also does not replace the time-integrated defect norms.
Missing genuine norm inputs remain `None`, and the gate stays OPEN.

For a normalized vacuum projector $P=(I+n\cdot\sigma)/2$, exact $P_*$,
$d=\|P_*-P\|_{\rm op}$ and $H=h\cdot\sigma$, rank-one geometry improves
the incoming energy estimate to

$$|\operatorname{tr}[H(P_*-P)]|
\le2|h\times n|d+2|h\cdot n|d^2,$$

because $n\cdot(n_*-n)=-2d^2$. Vacuum current error is exactly zero
when both traces are one. At the incoming surface
$h=(-m,\lambda/r,-E/a)$; the inherited $c_1=\lambda/r+im$ cancels
the constant cross-product coefficient. Directed bounds $B(E)$ for the
cross product and $D(E)$ for the projector distance yield
$2B(E)D(E)+2(E/a+M)D(E)^2$, with
$M=\sqrt{m^2+(\lambda/r)^2}$. Integrate the positive polynomial products
on the original finite bands and count signed multiplicity once; the lapse
action projection remains $4\pi ar^2$. Pressure bounds remain separate.

The historical eight-cell middle lapse enclosure improved from `4.18017e-5`
to `1.15107e-7`; group14's 128-cell values were `3.24091e-10` at order16
and `2.51409e-10` at order24. The latter exceeded `3e-11`. Its conditional
aggregate near `4.998e-9` required the separately recorded order24 lapse
correction `-3.30349e-12`; it is not the current missing low/subgap bound.
Centered depth4 Taylor enclosure at sixteen cells moved `4.37e-6` only
to `4.35e-6`. These finite numerical outcomes leave modal, quadrature,
field and complete-spectrum errors separate from the physical constraints.

The discrete search ledger declares its merit before any step: `m(g) = max over nodes of max(|E_N|, |E_beta|)`. Rows remain append-only, and a predicted merit remains a prediction. This finite search diagnostic does not close the physical gate.
