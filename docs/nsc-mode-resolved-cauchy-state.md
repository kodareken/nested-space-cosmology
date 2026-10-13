# The retained Gaussian mode state is now serialized

The previous source records retained integrated stresses but discarded the
mode matrices needed by unitary evolution.  This owner performs one authorized
base extraction from the unchanged kernels and stores every retained physical
Gaussian covariance in a deterministic, content-addressed NPZ.

The state contains 33 channels and 1,904 physical $2\times2$ blocks:

| Sector | Channels | Nodes per channel | Covariance blocks |
|---|---:|---:|---:|
| LLL, $(j,n)=(0,0)$ | 1 | 48 | 48 |
| Massive angular, $j=0$, $n=1,\ldots,12$ | 12 | 48 | 576 |
| Positive compact, $j=1,2$, $n=0,\ldots,9$ | 20 | 64 | 1,280 |

Identical compact copies and angular degeneracies remain multiplicities; their
matrices are not duplicated.  The canonical basis is

$$
\chi=r\sqrt{a_\parallel}\,\psi,
$$

with channel, frequency and then upper/lower spinor ordering.  In this basis

$$
H_{jn}(k;T)
=-m_j\sigma_1+\frac{\lambda_n}{r(T)}\sigma_2
 +\frac{k}{a_\parallel(T)}\sigma_3.
$$

The LLL covariance is reconstructed on the same 48-node quadrature from its
existing transparent occupation law.  Its two analytic stress constants are
reproduced within $2.8\times10^{-13}$.  The state-independent conformal term
remains a separate, one-counted local allocation.

For massive angular channels, the payload stores the physical covariance and
the existing $E_2+E_4/P_2+P_4$ subtraction arrays.  For positive compact
channels it stores both the physical covariance and the fourth-order
superadiabatic seed reference.  The locked Wilsonian compact tensor remains a
separate local contribution added exactly once.

Reloading the payload with `allow_pickle=False` verifies Hermiticity, the CAR
eigenvalue interval, channel order, offsets, dtypes, array hashes and the full
payload digest.  Reducing the stored blocks and the declared local allocation
reconstructs the old unit-radius tensor within its inherited tolerance.

## The Landau endpoint is not yet a Cauchy map

The recorded stress rapidity fixes the local current metric

$$
B_L=e^{-\eta\sigma_2},
\qquad B_L^{1/2}=e^{-\eta\sigma_2/2},
$$

but $B_L^{1/2}$ is not unitary in the old Cauchy norm.  Applying it directly as
$C\mapsto B_L^{1/2}CB_L^{1/2}$ can move a covariance eigenvalue above one.
The CAR-preserving map must instead be the full Dirac Cauchy isometry

$$
\boxed{
U_{L0}=J_L\,\mathrm{Res}_{\Sigma_L}E_g
       \mathrm{Res}_{\Sigma_0}^{-1}J_0^{-1},
\qquad U_{L0}^\dagger U_{L0}=I.
}
$$

It needs the Landau hypersurface embedding and the intervening metric, gauge,
Dirac and boundary history.  Endpoint $(v,\eta,r_\star)$ does not determine
that history or its frequency-mixing kernel.

The source-selected finite stress also needs the existing fourth-order
reference and local induced allocation generalized from the Bronnikov profile
to the same Kantowski--Sachs history.  The serialized seed reference alone is
not an instantaneous replacement on $r_\star$.

Consequently this record passes state serialization, CAR, and old-surface
moment reconstruction.  It retains the already verified nonzero direction

$$
\left.\partial_{\log r}\rho\right|_C=-0.460687999769215,
$$

but does not label the linearized displacement as the finite $r_\star$ stress.
The full user gate remains open for the physical Cauchy isometry and the
general-history reference/local providers.  Metric evolution remains stopped.

The fast verifier is

```sh
python3 scripts/derive_nsc_mode_resolved_cauchy_state.py --check
```

It authenticates and reduces the stored payload without rerunning a source.
The original extraction was one base run that produced the state artifact and
legacy reductions from the same in-memory vectors; no standalone historical
generator was launched.

## Common-history light allocation

The state-independent LLL term uses $c=|q_{\rm mag}|=4$ once on the actual
KS history $ds_2^2=N^2dT^2-a^2(dz+\beta dT)^2$. Write $L=\log a$ and
$H_a=[a_T-(a\beta)_z]/(Na)$. Subtracting the consistently normal-ordered
chiral state moments from the covariant conformal tensor gives

$$
\Delta\rho=\frac{c}{24\pi}
\left[\frac{2L_{zz}-L_z^2}{a^2}-H_a^2\right],\qquad
\Delta T_{\hat0\hat1}=\frac{c}{12\pi a}\partial_zH_a,\qquad
\Delta p_\parallel=\Delta\rho-\frac{cR_2}{24\pi},\qquad
\Delta p_\perp=0.
$$

This reference-dependent difference becomes covariant with the state source.
Its radial four-dimensional allocation is $\Delta t_{ab}/(4\pi r^2)$.
Up to the retained boundary terms, the action representative is

$$
S_{\rm LLL,geo}=\frac{c}{24\pi}\int dT\,dz
\left[-NaH_a^2+\frac{2N_zL_z}{a}-\frac{NL_z^2}{a}\right],
$$
$$
(\delta_NS,\delta_\beta S,\delta_aS,\delta_rS)
=(-a\Delta\rho,-a^2\Delta T_{\hat0\hat1},N\Delta p_\parallel,0).
$$

It replaces the old homogeneous geometric evaluation. Adding that old
tensor afterward would count it twice. The thermal/state piece remains
separate; $93.54264532195464$ remains the historical seam time-node
diagnostic, not a new matching criterion.

The same charged angular continuation supplies

$$
Z_q(-1)=\frac1{30}-\frac{q^2}{6},\qquad
Z_{q,\ne0}(0)=-|q|-\frac13,\qquad
\mathrm{FP}_{s=1}Z_q(s)=4\gamma_E-2H_{|q|}.
$$

Here $Z_q(s)=4\sum_{\ell=b+1}^\infty\ell(\ell^2-b^2)^{-s}$,
$b=|q|/2$, the physical sector is $|q|=4$, $H_4=25/12$, and $\mu=1$.
The finite cylinder density and fourth-order harmonic are

$$
C_q(R;\mu)=-\frac{Z_q'(-1)+[\log(\mu^2R^2)+1-\gamma_E]Z_q(-1)}
{16\pi^2R^4},\qquad
h_q=\log\mu+\frac18\mathrm{FP}Z_q(1).
$$

For $a=b+1$, $\varrho=(b/a)^2$ and first omitted $J\geq3$, the differentiated
series tail obeys
$|\Delta Z_q'(-1)|\leq4a^3\varrho^J[1+a/(2J-4)]/[J(J-1)(1-\varrho)]$.
The charged trace adds $q^2/(48\pi^2r^4)$ to $\mathcal A_{4,0}$, matching
the imported $(2/3)F^2$ heat coefficient. Nonzero-angular restoration
subtracts the LLL cocycle; adding the LLL geometric term restores exactly
$\mathcal A_{4,q}$.

In the canonical KS foliation the complete light-restoration action is
$S_{\rm cyl}+S_{\bar R_2^2}+S_{\rm LLL,geo}[\bar g_2]+W_{4,q}$, with
$\bar g_2=g_2/r^2$, $\sigma=\log r$,
$S_{\rm cyl}=-4\pi C_q\int\sqrt{|\bar g_2|}$ and
$S_{\bar R_2^2}=h_q\int\sqrt{|\bar g_2|}\bar R_2^2/(240\pi)$.
Its Weyl coefficients are $(-1/20,11/360,-1/30)$ and the charged addition is
$(2/3)\sigma\bar F^2$, $\bar F^2=q^2/2$. All four metric directions retain
this same allocation. A silent substitution of PG slices changes the
canonical prescription and is rejected.

## The reference symbol and its band action

On transmitting histories the formal PST subtraction uses

$$
H_T=-\beta_KkI+N_K[-m\sigma_1+\lambda\sigma_2/r_K+k\sigma_3/a_K],\qquad
P\star P=P,\qquad i\varepsilon\partial_TP=[H_T,P]_\star,
\qquad P=\sum_{j=0}^4\varepsilon^jP_j+O(\varepsilon^5).
$$

The band is the negative normal-energy band, rather than the sign of shifted
coordinate energy. Its gap-zero point is excluded without a denominator
floor. Coordinate time and spatial Moyal derivatives stay explicit.
The scalar first correction is
$\mathrm{Tr}P_1=-m\lambda\partial_zr_K/
[2a_Kr_K^2(m^2+\lambda^2/r_K^2+k^2/a_K^2)^{3/2}]$;
pointwise trace one and CAR bounds therefore belong to the physical
covariance, not each truncated symbol. The homogeneous restriction matches
all 1,904 stored reference blocks and the spatial check covers 63 signed
families. It creates no new vacuum.

For regular high-momentum jets, set $\mu=1/|k|$, $e=\mathrm{sgn}k$,
$P_\infty=(I-e\sigma_3)/2$, $b=a[-m\sigma_1+(\lambda/r)\sigma_2]$,
$b_2=a^2(m^2+\lambda^2/r^2)$ and $s=\sqrt{1+\mu^2b_2}$. Then

$$
P_0=P_\infty+\mu B_0,\qquad
B_0=-\frac b{2s}+\frac{e\sigma_3\mu b_2}{2s(s+1)},\qquad
P_j=\mu^{j+1}B_j\ (j\geq1),\qquad \partial_k=-e\mu^2\partial_\mu.
$$

The finite $B_j$ representation preserves decay at $\mu=0$ without
subtracting nearly equal projectors. Scaled defects are
$G_n=\mu^{n+2}\widetilde G_n$ and $F_n=\mu^n\widetilde F_n$.
The retained tables use physical order 9 and $\mu$ order 5, leaving orders
$9-j$ and $5-j$ in $B_j$. These auxiliary jets do not replace
$C_{\rm up}$, and power counting alone does not close the physical UV
remainder or the local incoming gate.

For a fixed rank-one chart $\Pi$, the imported effective-band construction
uses

$$
S=\Pi P+(I-\Pi)(I-P),\quad
U=S\star(S^\dagger\star S)^{-1/2}_\star,\quad
h=U\star H\star U^\dagger+i\varepsilon(\partial_TU)\star U^\dagger.
$$

The inverse square root is a star inverse. With $D=\delta U\star U^\dagger$,
$\delta h=U\star\delta H\star U^\dagger+[D,h]_\star
+i\varepsilon\partial_TD$. Through order four, the paired, complete
**unweighted** momentum trace gives

$$
\left[\delta S_{\rm band}
-\tfrac12\mathrm{Tr}(P\star\delta H+\delta H\star P)\right]_{\rm bulk}=0.
$$

The only nondecaying primitive is odd in $\lambda$ and cancels in the 30
equal angular pairs. Three zero-angular families give the full 63 signed
families from 33 groups; each sign has weight
`copy_count * degeneracy / n_signs`. Higher momentum-boundary orders decay
as $O(1),O(|k|^{-1}),O(|k|^{-2}),O(|k|^{-3})$. Compact time/spatial variations
are required. Finite momentum boundaries and nonconstant regulator weights
retain their exchange terms and finite heat-scheme matching remains open.

## Tilted transmitting interface

The weighted endpoint metrics are
$G_0=\mathrm{diag}(w_i)\otimes I_2$ and
$G_L=\mathrm{diag}(w_i)\otimes B_L$. Every solution of

$$
R_I[K]=K^\dagger G_LK-G_0=0
$$

has $K_c=G_{L,c}^{-1/2}V_cG_{0,c}^{1/2}$ with $V_c\in U(2n_c)$; its
canonical map is exactly $U_c=V_c$. Dense frequency mixing within each
channel is admissible, while the angular and compact labels stay fixed.
The payload has 13 blocks of dimension 96 and 20 of dimension 128, for
total dimension/rank 3,808. The interface has 894,976 real kernel variables,
Jacobian rank 447,488 and linearized nullity 447,488. Identity and dense
Fourier witnesses preserve CAR and the full covariance spectrum, but the
endpoint equation does not select the physical $V_c$. The intervening
Cauchy embedding/scattering operator remains open. Original merged records
are at [Zenodo](https://doi.org/10.5281/zenodo.23090015).
