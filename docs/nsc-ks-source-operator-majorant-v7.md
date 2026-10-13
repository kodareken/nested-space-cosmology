# Family14_1 source insertion v7

v7 joins the accepted v6 aggregate with the 108 disjoint subgap rows. It does not add the v5 total to the v6 total. v6 already contains the v4 rows plus the 18 threshold rows. v5 contains the v4 rows plus these same 108 subgap rows. The two added sets meet only on the shared 1042 rows.

Family `(14, 1)` now has 1168 of 1168 signed rows. Missing keys: none. No new source capture and no new solver. The saved subgap proof replay was not repeated. Its tree hash is unchanged, `c310e9023e945fad5985651cc0811e0f737bdff6d8eab352e091b52ac05ee161`. Original weights are applied once. The local gate stays OPEN.

## Partial N and beta

The continuous insertion of these 1168 rows is

$$
N\le 3646674054457366029024452496769420151764188105426118876609\times 2^{-229},
$$

$$
\beta\le 5793343250846686183515013928067141454707331929528967339969\times 2^{-230}.
$$

Read as decimals, those bounds are about \(4.226958072607968\times 10^{-12}\) and \(3.357610065481873\times 10^{-12}\). They are larger than the accepted v6 bounds because the 108 subgap rows are now inside the sum. They are not a bound on every source family, on quadrature, on numerical field error, or on value-assembly arithmetic.

The stored weighted moments are mantissa `7786760115494097` at exponent `-86` and mantissa `5541269878732561` at exponent `-83`.

## Signed-zero identity

NumPy 2.5.1 drops the negative zero at `C[0, 1].real` of each negative fiber while assembling the negative partner. The numeric arrays remain exactly equal: nonzero differences 0. Restoring that one sign bit on each fiber recovers the negative source and preparation digests stored in v6. The low32 batch's two raw source digests are still `641b187a3c37aa47a5841b559163eba5f78ec3ca2c79b535dc21dae8db657bdf` (signed-zero spelling) and `63e52bb2986c13bd72d91f3965a6f48b3d11380a6d29ae8326658c4fb2bc8d2c` (cleared spelling, which is the live digest). The same corner restores every other negative batch digest used by v6 and by the subgap rows. 584 negative rows are therefore recorded under the live cleared digest. Positive digests were already the live digests. Energies were unchanged.

The historical propagator object checked by the threshold closure is git blob `31818fdd6b8cc9c4d8b9d33f41e0f241b2dd1227`, not HEAD. The closure checked 67 source hashes and 64 archive digests at `eebbe0c3b0bfe76be0bae4a64d3b52ba8da7cfb9`.

## Budget

Completing family `(14, 1)` removes the remaining signed-row gap inside this one positive-angular family only. `physical_upstream_budget_component` stays null. Quadrature error, numerical field error, value-assembly arithmetic, and every other source family stay out. This record does not close the physical local gate.

## Files and checks

The aggregate owner is `scripts/derive_nsc_ks_source_operator_majorant_v7.py`. The result is `results/development/nsc-ks-source-operator-majorant-v7.json`, sha256 `c9a98578bf53ee02903403f4b63550ddb1d67ecf870d656f7b1f2ee74de44c12`. Accepted v6 remains sha256 `15615d1ee0206ed0fb6b8c3a10df8d860c581b2ebc206bef2242492f0154d517`. v4 and v5 bytes were not rewritten.

From the repository root, after the validation environment is ready:

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_ks_source_operator_majorant_v7.py --record
```

`--record` checks the witness tree, payload hashes, threshold closure, signed-zero spelling, disjoint census, mutated preparation, negative epsilon, and duplicate rows, then writes the file once. It does not call the 54-row proof replay. `--check` repeats those record checks and does not write when the bytes match. The validation environment is Python 3.12.14 with the recorded numerical pins. Its test runner, installed without upgrading those pins, is pytest 8.4.2 from `results/development/nsc-threshold18-import-v1/requirements.txt`. `tests/test_nsc_ks_source_operator_majorant_v7.py` passed, and `--check` reproduced sha256 `c9a98578bf53ee02903403f4b63550ddb1d67ecf870d656f7b1f2ee74de44c12`.


## Retained method chain and historical exclusions

The archived 2026-09-28 continuation recorded the local gate OPEN with
Stages2–5 incomplete. Stage2 lacked the physical source-column bound and
changed-history UV remainder; its family14 IR example was
$E=0.0013248831260437577$. Later v7 completes only the signed-row insertion
stated above. The scoped earlier exclusions remain: homogeneous
frequency-diagonal KS without an interface had shift residual about
`-0.0135927`; smooth compact histories with exact frozen incoming $C_0$
obeyed the strict bound $P_K<-0.022178625$, not an asserted exact value.
Compensator, first-order and jet-rigidity results retain that frozen-class
domain. Principal-matrix, flux and conditional local surface existence are
reusable geometry, not existence of the coupled prepared state. Certified
background regions and group32 `[24,32]` remain their original action
certificates. No double counting or freely chosen compensation is allowed.
These historical exclusions neither prove non-existence of the switched
state law nor authorize another campaign.

### Continuous preparation-error insertion

The v1 method behind the aggregate compares exact canonical propagators
over the backward cone. The carrier energy is a phase; characteristic row
growth depends on the off-diagonal potential. If $B$ bounds its integral,
and $M,M_z$ bound radius-difference and axial-coupling integrals, set

$$f=\sqrt2,\qquad d=\sqrt2e^BM,\qquad
d_{z,0}=\sqrt2e^BM_z,\qquad A=2fd+d^2,\quad J_0=(f+d)d_{z,0}.$$

The physical derivative of the difference retains both $|E|d$ and
$d_{z,0}$. For per-energy unweighted covariance errors $\epsilon_E$ and
original quadrature weights,

$$S_0=\sum_E\frac{w_E}{2\pi}\epsilon_E,\qquad
S_1=\sum_E\frac{w_E}{2\pi}|E|\epsilon_E,$$

$$\epsilon_\beta\le\mu(AS_1+J_0S_0),\qquad
\epsilon_N\le\mu\left[
\sqrt{m^2+\ell^2/r_\Sigma^2}\,AS_0+
\frac{AS_1+J_0S_0}{a_\Sigma}\right].$$

Both energy signs have their own covariance comparisons; weights and
signed-family multiplicity enter once, without a three-column factor.
The history-minus-reference difference is formed before the trace norm.
Zero deformation cancels this insertion but leaves baseline uncertainty.
Numerical field errors at the saved preparation are separate.

The physical $\rho=1$ ledger distinguished102 low/subgap classifier
panels and29356 rows (64 low,38 subgap) from threshold coverage:27828
rows below V5 joined thresholds,1528 at/above, and59 straddling panels.
The62 middle panels/20016 action-covered rows were not automatically
physical column bounds. Thus the source-accuracy universe was all164
panels/49372 positive rows until a proved bridge; production panel splits
must handle threshold straddling. Group12/32 order24 action successors
are reusable, but the old aggregate near `4.998e-9` is not the present
missing low/subgap quantity. Groups0,13,23 have exact ell-zero radius
response with baseline retained. Missing physical bounds stay null.

Numerical continuation from $\rho=1$ to the exact archived
$\rho_{\rm up}=1159676904047903/1125899906842624$ obeys

$$\|e(\rho_{\rm up})\|_2\le\|e(1)\|_2+
\int_1^{\rho_{\rm up}}\|Y'-GY\|_2\,d\rho$$

for the anti-Hermitian homogeneous generator. The DOP853/Bernstein
residual and order8 background Taylor remainder bound only continuation.
The physical initial $\rho=1$ error remains separate and must later enter
the correlated reference/difference propagation, rather than a substituted
full-reference stress estimate.

### Exterior, horizon and subgap source enclosures

In the original exterior coordinate $r$ (areal radius $\sqrt{1+r^2}$),
$A=1-3[(1+r^2)\arctan(1/r)-r]$, $V_1=\ell\sqrt A/\sqrt{1+r^2}$,
$V_2=m\sqrt A$. The subgap decaying ratio $z=e^{i\theta}$ satisfies

$$\theta'=\frac2A(V_1\cos\theta+V_2\sin\theta-E),\qquad
\theta(\infty)=\arcsin(E/m),\quad0\le E<m.$$

A whole-half-line polynomial defect enclosure $\varepsilon$ and positive
phase-tube lower $\lambda$ give outer error $\varepsilon/\lambda$,
including the actual floating initializer's rounding. The metric tail uses
the positive majorant of its convergent inverse-radius series; samples or
formal cancellation cannot replace it.

Radial transport uses the actual DOP853 anchored polynomial in
$y=\log(r-\hat h)$ and normalized defect
$\delta=\theta_x-hf_y$, with original $h<0$ entering once. If
$\lambda\le f_{y,\theta}$ on a valid metric/phase tube, put $k=-h\lambda$:

$$|e(1)|\le e^{-k}|e(0)|+\varepsilon\varphi(k),\qquad
\varphi(k)=\int_0^1e^{-kt}dt,
\quad\varphi(0)=1.$$

The signed successor permits negative $k$ with finite growth when both
endpoint errors stay inside the tube; intervals containing zero use the
integrated exponential or monotone endpoint bounds without dividing by
zero. This unblocked group14/low16_1 row15,
$E=0.24867511687395624$, whose843 cells included23 expanding cells.
Its pilot had inner phase upper `1.814e-9`, frame reflection tail
`6.460e-22`, Bloch upper `2.027e-10`, and positive/negative covariance
uppers `1.047e-10`. The earlier row0,
$E=0.0013248831260437577$, had906 transport cells and about `9.04e-12`
phase error after four defect subdivisions. These are scoped row proofs,
not N/beta gate bounds. Metric tails are composed into the actual cell
coordinate, and unresolved metric, tube or endpoint tests raise.

The full-profile horizon frame uses
$q=\pi/2+\arctan\rho$, its enclosed simple root $q_h$, and

$$W(q)=3(\pi-q)+\tfrac32\sin(2q)-\sin^2q,\qquad
H(u)=-W(q_h+u)/u,\quad H(0)=2\kappa_{\rm geometry}.$$

For $t=\sqrt u$, $F_s=u^{\lambda_s}\sum c_{s,n}t^n$,
$\lambda_s=isE/H(0)$, and
$K(t)=iE\sigma_3/H(t^2)+tP(t^2)/\sqrt{H(t^2)}$,
$P=-m\csc(q_h+u)\sigma_1+\ell\sigma_2$, the recurrence is

$$\left(\tfrac n2I+\lambda_sI-K_0\right)c_{s,n}
=\sum_{j=1}^nK_jc_{s,n-j}.$$

The inverse has induced1-norm at most $2/n$. Analytic disk bounds
majorize the full series by $(1-x/R)^{-2M}$; a positive-term ratio bound
encloses its entire omitted tail. Continue with $t=-i\sqrt\delta$ inside,
keeping the interior oscillatory factors and exact radial Jacobian
$u\,r(\rho)^2$. Geometric slope never replaces source occupation kappa.
The frame maps an enclosed exterior phase to

$$R_{\rm direct}=
\frac{-F_{21}+F_{11}e^{i\theta}}{F_{22}-F_{12}e^{i\theta}},$$

with denominator separated from zero. Its row0 order16 disk $R=0.1$
tail was below `2.50e-22`. Matching to the original producer needs
$D(c)$ with $c=E\log r_h^2/(2\kappa_{\rm geometry})-
\tfrac12\operatorname{atan2}(mr_h,\ell)$ and interior
$D_\pi=\operatorname{diag}(e^{-i\pi/4},e^{i\pi/4})$.
The sewing identity
$D(c)S(e^{2ic}R,T)=S(R,T)\operatorname{diag}(e^{ic},e^{ic},e^{-ic})$
preserves the full three-port covariance; altering only the reflection or
discarding relative coherence would not.

For complete homogeneous subgap covariance $Q=(I+n\cdot\sigma)/2$,
the actual interior equation is $n_y=2h\times n$, a norm-preserving
real generator. Continuous anchored-cell defects therefore add without
growth or a second width factor. The direct-frame horizon covariance is
$\left(\begin{smallmatrix}1-f&c_h\overline R\\c_hR&f\end{smallmatrix}\right)$,
$c_h^2=f(1-f)$; source kappa stays fixed. The saved covariance comparison
retains trace and Bloch discrepancies:

$$\|Q_{\rm true}-Q_{\rm num}\|_{\rm op}\le
\tfrac12\left(|1-t_{\rm num}|+\|n_{\rm end}-n_{\rm num}\|
+\epsilon_{\rm Bloch}\right).$$

The actual PG recovery/restriction maps differ by a common per-energy
phase, which cancels in this homogeneous covariance. At the exact upstream
slice the negative partner satisfies
$Q_-=I-\sigma_3\overline{Q_+}\sigma_3$, or
$(n_x,n_y,n_z)\mapsto(n_x,-n_y,-n_z)$. Compare that ball separately with
the actual negative columns to retain Gram/complement errors; copying the
positive scalar error can underbound them. This identity is not assumed
pointwise after inhomogeneous evolution. Saved row0 and row15 witness
replays validate phase/Bloch traces without repeating either solver, keep
covariance errors unweighted, and supply no complete upstream component.

### Vacuum and coherent-occupation remainder methods

For the original three-port affine source and $P=\operatorname{diag}(0,1,0)$,
the coherent horizon block squares to $fI$, hence

$$\|C_{\rm source}-P\|_{\rm op}=\max(\sqrt f,n_{\rm in}).$$

Normalized sewing and homogeneous unitary transport do not increase that
distance. The branch $n_{\rm in}=0$ is used only when an entire interval
lies strictly below the original incoming mass threshold; crossing or
equality keeps the nonzero allowance. The group14 row at
$E=2.071557525252559$ had analytic occupation upper below `1.383e-12`,
rather than the old generic `1.072e-6`. This tightens an estimate, not
the stored state; binary occupation/coherence arithmetic remains separate.

For the horizon vacuum Bloch branch, write $\delta=t^2$, $y=\log\delta$,
$L=\sqrt H(mr-i\ell)$, $w\sim t\sum P_jE^{-j}$ and
$n_z\sim\sum Z_jE^{-j}$, with $P_0=0,Z_0=1$. The recurrence is

$$P_{j+1}=\frac{iH}{2}(P_j/2+\delta\partial_\delta P_j)+LZ_j,
\qquad
Z_n=-\tfrac12\sum_{j=1}^{n-1}
(\delta\overline{P_j}P_{n-j}+Z_jZ_{n-j}).$$

Its exact transverse defect is $-2itP_{M+1}/(HE^M)$ and longitudinal
defect is zero through the retained polynomial. The skew generator gives

$$\|n-n_M\|_2\le C_M^{\rm prep}/E^M,\qquad
C_M^{\rm prep}=4\int_0^{\sqrt{\delta_{\rm up}}}
|P_{M+1}(t^2)|/H(t^2)\,dt.$$

Whole-cell directed integrals, not order/refinement differences, enclose
this constant; covariance error is half the Bloch error. The order8,
64-cell value at E=160 was below `8.46e-15`. The order16,128-cell pilot
compared384 positive rows at E≥32 and384 actual negative partners with
their stored columns, adding occupation/coherence and Gram errors. This
preparation constant is not the changed-history $C_M$ of the UV gate.

At moderate energy the small correction $e=n-n_M$ obeys
$e_y=A_Ee-r_M$. Starting its numerical curve at zero at $y_0=-18$
does not zero its true initial uncertainty:
$\|e(y_0)\|\le C_M^{\rm prep}(y_0)/E^M$. Each saved correction cell
has normalized defect $p_x-2h\,h_E\times p+hr_M$; initial, full-cell
and exact-endpoint bridge errors all enter the comparison. The first
group14/mid24_1 control near E=16.01925 compared both original signed
covariances. The middle-coverage adapter treated $16\le E<32$ as a
finite48-positive-row campaign requiring explicit row/CPU budgets;
partial coverage stayed OPEN, and replay solved nothing.

One upstream phase-fixed vacuum spinor can be chosen as
$u_0=\sqrt{(1+n_z)/2}>0$, $u_1=(n_x+in_y)/(2u_0)$.
Hold that common per-energy phase fixed during subsequent evolution. For
major lower $a_{\min}$, transverse bound $W$, and Bloch remainder
$C_nx^M$ at $x=1/E$,

$$|u-u(n_M)|\le
\left[\frac3{4a_{\min}}+\frac{W}{8a_{\min}^3}\right]C_nx^M.$$

An enclosing chart Taylor remainder $C_Tx^M$ gives their sum $C_u$.
The truncated Bloch polynomial need not be normalized; chart separation
includes its remainder. The A4 adapter bounds all discarded higher degrees
and supplies homogeneous envelope Sobolev initial data
$(\sqrt L\,C,0,0)$. Physical carrier derivatives remain present. The
negative occupied vacuum representative is $(v,a)$ for positive $(a,v)$,
not merely its sigma3 conjugate. Original thermal/coherent source matching
and propagated $L_0A_4$ bounds stay separate.

### Radius Fourier moments and actual Weyl defect

The p=16 profile geometry gives nonlinear
$q=1/(r_{\rm ref}+\delta r)-1/r_{\rm ref}$ spatial and Fourier-l1
bounds through order9, with original amplitude0.001, U=0 and period205/512.
For $A_j(f)=\sum|2\pi k/L|^j|c_k|+T_j$, the Wiener inverse uses
the positive Neumann margin $r_{\rm ref,min}-A_0(\delta r)$:

$$A_0(1/r)\le\frac1{r_{\rm ref,min}-d_0},\qquad
A_0(q)\le\frac{d_0}{r_{\rm ref,min}(r_{\rm ref,min}-d_0)}.$$

Positive-order moments follow the exact reciprocal Leibniz recurrence;
a pointwise radius lower alone does not prove the Wiener bound.
For this control $r_{\rm ref,min}=7/5$,
$d_0\le2.439484626404717\times10^{-5}$ and the margin lower is
`1.399975605153736`. Original amplitude is already inside W and is never
multiplied twice. Norm uppers and denominator lowers retain opposite
rounding directions.

The omega-only transfer helper uses
$q_4=\sum_{p=1}^4(-\chi s)^pW^p/r_{\rm ref}^{p+1}$ and bounds its
true reciprocal remainder by

$$A_j(q-q_4)\le j![t^j]
\frac{\eta(t)^5}{r_{\rm ref,min}(1-\eta(t))},\qquad
\eta_j=A_j(\delta r)/r_{\rm ref,min}.$$

It contracts retained mode balls, omitted profile tails and this remainder
against the nonnegative kernels, for $w=|\omega|$, $s\in\{0,1\}$,

$$f_{\rm small}(n,s)=\frac{(w/2)^n}{(5-s)\max(K,w)^{5-s}},$$
$$f_{\rm shift}(s)=1_{w>K}\frac{w^{s+1}-K^{s+1}}{s+1},\qquad
f_{\rm taylor}(n,p,s)=1_{w>K}(w/2)^n\int_K^wk^{s-p}dk.$$

These exclude $dk/(2\pi)$, the two k-signs and source/spin multiplicity;
callers apply them once. Both omega signs already enter the inventory.
K=1024 gave tighter nonnegative sums than K=256 on this control, but K is
a canonical momentum split, not an incoming source-energy cutoff.

For the actual auxiliary defect,

$$D=C_g-P,\quad R_\rho=P_\rho-[G_\rho,P],\quad
D(\rho)=U(\rho,u)D(u)U(\rho,u)^\dagger
-\int_u^\rho U(\rho,s)R_\rho(s)U(\rho,s)^\dagger ds.$$

Here $G_\rho=iH_\rho$ and $R_\rho=-R_T/a$; the kinetic and potential
Weyl commutators are exact shifted products:

$$[kC,p]_{\#,n}=k[C,p_n]+\frac{\omega_n}{2}\{C,p_n\},$$
$$[V,p]_{\#,n}=\sum_l
\left[V_lp_{n-l}(k-\omega_l/2)-p_{n-l}(k+\omega_l/2)V_l\right].$$

Taylor remainder bounds require derivative uppers over both closed shifted
momentum segments. For one mode, shift $v=\omega_l/2$ and order M,
the Frobenius upper is

$$\|V_l\|\frac{|v|^{M+1}}{(M+1)!}
\left(\sup_t\|\partial_k^{M+1}p(k-tv)\|_F
+\sup_t\|\partial_k^{M+1}p(k+tv)\|_F\right).$$

Formal order alone cannot replace this exact remainder. Fourth axial
derivatives of the J=4 defect need total mixed geometry order9. In the
small-transfer domain $|k|\ge K$, $|\omega|\le|k|$, the regular mu
interval must cover $[0,2/K]$. Its bound is
$\|\partial_z^qR_{T,\rm small}\|_F
\le C_{5,q}|k|^{-5}+C_{6,q}|k|^{-6}$, with

$$C_{6,q}=\sum_{j=0}^4\sum_{h=0}^q
\frac{128\binom qh|\ell|A_{5-j+h}(1/r_g-1/r_{\rm ref})}{2^{5-j}(5-j)!}
\sup\|\partial_z^{q-h}D[j+1,5-j]B_j\|_F.$$

The128 factor includes both segments and their minimum momentum $|k|/2$.
The homogeneous reference potential contributes zero to these positive
transfer moments. Nuclear Fourier moments can use

$$\sum_l\|a_l(k)\|_1\le\sup_z\|a(z,k)\|_1
+\frac{L^2}{12}\sup_z\|a_{zz}(z,k)\|_1,$$

then integrate in canonical k with its measure and signs once. An optional
bounded low-k auxiliary
$Q_{\rm aux}=P_0+\chi_K(P_1+P_2+P_3+P_4)$, with
$\chi_K=1-\mathrm{plateau}(|k|;K/4,K/2)$, leaves high-k small-transfer
identities unchanged. If used, its low-k residual and finite reconstruction
correction remain required. One-box small-transfer and integrated-moment
controls supply no whole-time, full-support, initial covariance or source
tail certificate; their positive norm estimates are not physical divergence.

### Cutoff edge, angular pairing and linear controls

For unchanged intrinsic incoming data the owned ordered kernel limit is

$$\Delta G_B^{\rm owned}=\Delta G_B^{\rm raw,\Lambda}
+\Delta G_B^{\rm edge}+\mathcal T_{B,\Lambda},$$
$$\Delta G_N^{\rm edge}=\frac\mu{2\pi a_\Sigma}(f'_+-f'_-),\qquad
\Delta G_\beta^{\rm edge}=-\frac\mu{2\pi}(f'_++f'_-).$$

The raw term contracts the actual retained source before coincidence;
the edge is from the original phase coefficient in the owned order of
limits, with original multiplicity and $de/(2\pi)$ once. It introduces
no new source, action or compensating force. The finite recurrence
$\Pi_sL_0A_j=0$,
$\Pi_{-s}A_{j+1}=a^2\Pi_{-s}L_0A_j/(2i)$ has defect
$e^{-M}L_0A_M$ and a conditional H2 error $C_Me^{-M}$.
Restoring the carrier costs one power; M≥3 makes the differentiated
remainder integrable, but does not evaluate its constant at cutoffs160/320.

Global unitarity alone does not remove the order-of-limits issue. In the
exact negative-Fourier circle control with
$B_{n+1,n}=1/[2(n+1/2)]$ and $U=e^{itB}$, the complete covariance is
unchanged, while the truncated separated kernel derivative is

$$\dot C_{N,\eta}(z)=
\frac{\sin z}{2\pi(N+1/2)}e^{-i(N+1/2)\eta},\qquad
-i\partial_\eta\dot C_{N,\eta}(z)=
-\frac{\sin z}{2\pi}e^{-i(N+1/2)\eta}.$$

The kernel tends uniformly to zero but its differentiated coincidence
does not. This control identifies a limit obligation and licenses no
particular NSC correction.

For equal-weight original angular pairs, real geometry gives
$A_j(-\ell)=(-1)^jD_s\overline{A_j(\ell)}$,
$D_s=sS_3$. Complete raw current coefficients, including density and
interference, have parity $(-1)^{p+1}$ at $e^{-p}$.
The $e^{-2}$ term cancels in the angular pair; retaining M=4 gives a
paired integrand $O(e^{-3})$ and tail $O(\Lambda^{-2})$. This strengthens
the older valid $O(\Lambda^{-1})$ estimate without supplying a numerical
constant or asserting that the first allowed cubic coefficient is nonzero.
Thermal/coherent and retarded-variation remainders remain required.

For linear fixed transfer, $Q_2=-\ell a^2\hat w\sigma_1/(4r^2)$ and
the projector relation forces
$(Q_3)_{\rm diag}=m\ell a^3\hat w\sigma_3/(4r^2)$.
With $H_0=-m\sigma_1+\ell\sigma_2/r$, $H_1=-\sigma_3/a$,
$\operatorname{tr}(H_1Q_2)=
\operatorname{tr}(H_0Q_2+H_1Q_3)=
\operatorname{tr}Q_2=\operatorname{tr}Q_3=0$.
Dropping forced $Q_3$ would leave a false N coefficient.
This leaves a contracted linear $O(E^{-3})$ vacuum remainder, with
constants not uniform in axial transfer.

The separate large-transfer linear estimate retains both source/output
measures $dE_i\,dE_o/(2\pi)^2$ and symmetric momentum
$-(E_o+E_i)/2$. At $\sigma=.03$, $r^2\ge2$,
$I_w\le9/40000$, $I_U\le27/1600000000$. Put
$D_p=I_w\|w^{(p)}\|_1+I_U\|U^{(p)}\|_1$ for p>3, and

$$A_p=\frac{2^{p+1}}{(p-1)(p-2)},\qquad
B_p=\frac{2^p}{p-3}\left(\frac4{p-1}+\frac1{p-2}\right).$$

On $|E_i|\ge E_0$, $|\omega|\ge|E_i|/2$, both energy signs and
transfer tails give

$$|\delta G_N|_{\cal T}\le\frac{\mu|\ell|D_p}{(2\pi)^2}
\left[4(m+|\ell|/r)A_pE_0^{2-p}+
\frac2aB_pE_0^{3-p}\right],$$
$$|\delta G_\beta|_{\cal T}\le
\frac{2\mu|\ell|D_p}{(2\pi)^2}B_pE_0^{3-p}.$$

Here $B_4=88/3$, $B_6=112/5$; W4,1 gives current power E0^-1,
whereas E0^-3 needs W6,1. These bounds exclude the complementary
frequency region and finite-history covariance remainder.

### Numerical-control protocols retained from retired producers

The September25 repair corrected the true CF4 Euclidean adjoint, including
Fourier signs, spatial-replication sums and both Y-to-A/Y-to-D couplings.
Primal-block DOP853 control avoids tangent-count dilution. Version2 caches
bind actual source arrays, complete history, channel ledger, settings,
nodes and implementation; historical version1 caches are not current inputs.
Reconstruction uses $U_g=U_{\rm ref}+(U_g-U_{\rm ref})$ with
`reference_mode=direct-original-energies`: reference propagation uses original
energies, and only the difference is interpolated. The complete-family
zero-history changes became exactly `(0,0)`, without subtracting a measured
stress drift. Degree32→48 movements were
`(2.52055321e-15,2.42102345e-15)` for32_-1 and
`(5.59956541e-17,4.15486601e-17)` for14_1. These are refinement indicators.
The exact degree48 difference-interpolation component was
`(4.831544944e-26,3.784029101e-26)`; node solves, direct-reference error,
source, quadrature, field, UV and arithmetic remain separate. Degree32's
conservative `(2.213119788e-9,1.733298516e-9)` exceeded its allocation.

The fine protocol captured one original14_1 value-only history,
amplitude0.001, U=0, N=2048, exact numerical period
$205/512=0.400390625$, `rtol=5e-14`, `atol=5e-19` and
`max_step=1/8192`. Its streaming chunks contain at most four accepted
segments with bitwise endpoint continuity. Capture is data availability,
not proof. The geometric table uses p=16, M=65536, K=8192 on that exact
period and unchanged physical support; full residual validation covers
all246 saved steps, including the cell122 pilot. Each continuous cell
includes polynomial, time, radius and profile remainders; no partial sum
is promoted to complete coverage.

At a causally unaffected point $z_0$ outside both padded support copies,
with the owned distance bound $(\rho_{\rm up}-1)/(4/5)^2$, the reference
row error obeys
$\|A_{\rm ref,true}-A_{\rm ref,num}\|_{\rm row}
\le e_0+\|D_{\rm num}(z_0)\|_{\rm row}$.
The numerical difference is retained, never forced to zero. Current and
reference errors both enter coherent matter contraction. For the literal
stored covariance use
$\gamma\ge\|C\|_2\le\sqrt{\|C\|_1\|C\|_\infty}$ in interval arithmetic,
and weighted field bounds f,z,e,ez give

$$D=\gamma(2fe+e^2),\qquad J=\gamma(ze+fe_z+ee_z),$$
$$\epsilon_N\le\mu\left[
\sqrt{m^2+\ell^2/r_{\min}^2}\,D+J/a_{\min}\right],
\qquad\epsilon_\beta\le\mu J.$$

The actual axial Fourier derivative is $2\pi k/L-E$. Ball evaluation
also compares the mathematical reconstruction to stored binary fields,
contraction and intrinsic geometry. The finite14_1 numerical-control
target is `1e-11`; even a pass neither certifies omitted source spectrum
nor solves either incoming constraint.

Retired adapters retain these distinct roles:

| Control | Preserved operation and scope |
|---|---|
| Evolved-state and constraint join | Fixed coherent source, explicit $\delta C_{\rm src}=0$, actual F/Fz tangents; saved-array constraint replay propagates no fields. Field derivative tolerance `3e-8`, flux algebra `3e-11`; full continuum/source errors remain open. |
| Same-kernel E/k/eta bridge | Saved source energies, independent canonical-k quadrature on `[-4,4]` with24 nodes, identity tolerance `3e-11`; source/reference tails and ordered coincidence limits remain unbounded. |
| Exact source phase | One bounded numerical-phase diagnostic with30-CPU-second cap; array replay, no fitted phase or full-source/continuum certificate. |
| Local reference accuracy | Same four-energy source on the approved I, h=.0005, existing951/3801 buffer inputs and90-CPU-second cap; partial reference indicator, with no nonzero history or stress-drift subtraction. |
| KS/PG envelope comparison | Same PG arrays, fibers and upstream columns; at most5 solves/120 CPU seconds. The successor's two512-node checks resolve spatial indicators, not error bounds; no PG/source regeneration. |
| Source-control comparison v2 | Sum positive/negative matter changes within each family before taking componentwise maxima, then triangle-sum selected-family movements. Saved-array replay performs no field evolution and does not certify full source. |
| Retained response sum | Validate original rows, panels, quadrature measure, signed-family identities and implementation before summing; partial/full finite control status stays OPEN without source, interpolation, UV or between-node bounds. |
| KS spacetime variation | Weak metric/CTP direction controls on saved modes retain `dE/(2*pi)` once; no nonlinear U[g], full transmitting endpoint family, complete renormalized variation or new stress is supplied. |
| Low-energy vacuum enclosure | One directed projector witness or propagation-free replay; no source replacement, column normalization or inferred `3e-11` physical accuracy. |
| Middle and paired-tail replay | Finite primitives and the unchanged Riccati16/ad4 source approximation, with coherence and original signed counting; no radial/field generators. Measured truncation/quadrature agreement is not a rigorous complete-source remainder. |
| Tail-quadrature control | Original group22, cutoff160, Riccati16/reference4; compare certified infinite integral with archived order13 values and add separately owned vacuum/thermal bounds. One group and unchanged source values do not certify all32 group tails or low modes. |

Historical method records and retired producer code are available in the
[archived evidence DOI](https://doi.org/10.5281/zenodo.23090015).
The retained v7 aggregate above remains the current source-insertion owner:
1168/1168 signed rows in family(14,1), physical upstream component null,
quadrature, field, value arithmetic and other source families excluded,
physical local gate OPEN.
