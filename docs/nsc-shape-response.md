# Quantum stress difference on the same smooth varying-radius geometry

The [calculation](../scripts/check_nsc_shape_response.py) and
[record](../results/nsc-4-shape-response.json) compute the Dirac vacuum energy
and local stress **difference** between periodic and antiperiodic axial spin
structures on the actual smooth radius profile. They vary all lapse, radial,
and sphere-radius metric degrees of freedom. No cylinder pressure or local
density approximation is substituted into a neck.

This is a measurable mathematical advance toward the proposed self-sourcing
geometry: changing spin structure lowers the computed neck null source on the
same varying geometry. It does **not** determine the absolute antiperiodic
stress or match the complete quantum stress to the geometric Einstein tensor.

## The four-dimensional operator and its metric variations

Use signature \((+---)\), natural units \(\hbar=c=1\), and

\[
ds^2=N(x)^2dt^2-q(x)^2dx^2-r(x)^2d\Omega_2^2,
\qquad x\sim x+L.
\]

The static Hilbert measure is \(qr^2dx\,d\Omega\). After angular separation,
the normalization \(u=r\sqrt q\,\psi\) changes it to \(dx\), and the diagonal
tetrad's radial spin connection obeys

\[
r\sqrt q\,\frac Nq\left(\partial_x+\frac{r'}r+
\frac{N'}{2N}\right)\frac{u}{r\sqrt q}
=c\,u'+\frac{c'}2u,\qquad c=N/q.
\]

The implementation checks this identity and the measure removal symbolically.
The radial Hamiltonian is

\[
H_\kappa=-i\sigma_2(c\partial_x+c'/2)+v\sigma_1,
\qquad v=N\kappa/r.
\]

The node-to-edge block uses
\(c_e=(c_i+c_j)/2\), \(j=i+1\),

\[
A_{e,i}=-\frac{\sqrt{c_ec_i}}h+\frac{v_i+v_j}4,
\qquad
A_{e,j}=\frac{\sqrt{c_ec_j}}h+\frac{v_i+v_j}4.
\]

Its seam receives phase \(+1\) for periodic and \(-1\) for antiperiodic spin
structure. At constant \(c=1,v\), an independent Fourier calculation gives
singular values squared
\(4h^{-2}\sin^2(ph/2)+v^2\cos^2(ph/2)\), and agrees with the matrix SVD.

For one complex four-component massless Dirac field,
\(E_\eta=-4\sum_{\kappa\ge1}\kappa\sum_j s_j(A_{\kappa,\eta})\).
Both sphere signs and the particle/antiparticle determinant factor are included,
with the same counting as the independently derived
[cylinder benchmark](nsc-covariant-measure.md). Each angular sector's periodic
and antiperiodic quantities are subtracted before summing to reduce cancellation
of extensive ultraviolet terms.

The polar-factor Hellmann–Feynman identity
\(\delta\sum s_j=\operatorname{ReTr}[(UV^\dagger)^\dagger\delta A]\)
gives derivatives with respect to **every** \(N_i,q_i,r_i\), retaining the
variation of the edge coefficient. Only after differentiation are \(N=q=1\)
set for the reported smooth profile. Independent finite changes in each metric
field, both smooth and localized at the throat node, check these derivatives
on a nonconstant-lapse, nonconstant-\(q\) background.

Let \(g_N,g_q,g_r\) denote derivatives of
\(\Delta E=E_{\rm P}-E_{\rm AP}\). The local stress difference is

\[
\Delta\rho=\frac{g_N}{4\pi qr^2h},\quad
\Delta p_x=-\frac{g_q}{4\pi Nr^2h},\quad
\Delta p_\perp=-\frac{g_r}{8\pi Nqrh}.
\]

The discretization preserves the pointwise Weyl identity
\(Ng_N+qg_q+rg_r=0\), hence
\(\Delta\rho-\Delta p_x-2\Delta p_\perp=0\). This is the relative stress:
identical local anomalies and geometric counterterms cancel between spin
structures. It is not an assertion that either absolute quantum trace vanishes.
The global identity \(\sum_iN_i g_{N_i}=\Delta E\) independently checks the
lapse normalization. A centered spatial derivative tests

\[
(\Delta p_x)' +\frac{N'}N(\Delta\rho+\Delta p_x)
+2\frac{r'}r(\Delta p_x-\Delta p_\perp)=0.
\]

## Actual smooth-profile results and convergence

The radius is exactly the profile already used in
[the geometric Einstein-residual calculation](nsc-smooth-geometry.md):

\[
r(x)=\sqrt{1+\left(\frac{\sin(kx)}k\right)^2},\qquad
k=\frac{\pi}{2R},\quad L=2R.
\]

The sphere radius and all quantum stress components are evaluated on this
varying geometry. At 512 axial points and \(\kappa_{\max}=8\):

| \(R\) | \(\Delta E\) | Throat \(\Delta(\rho+p_x)\) | \(\int\Delta(\rho+p_x)dx\) | RMS conservation residual |
|---:|---:|---:|---:|---:|
| 2 | 0.1586544847 | 0.0186082412 | 0.0409127521 | \(4.44\times10^{-7}\) |
| 4 | 0.0194597880 | 0.0017988910 | 0.0041368012 | \(1.53\times10^{-7}\) |

All values use the declared throat-length unit; stress has inverse-length to
the fourth-power units. The record retains the **entire** finest-grid local
profile, all angular sectors, full-grid trace and conservation residuals,
unweighted proper axial null integrals, and sphere-weighted integrals. The
independently calculated geometric Einstein components accompany each point.

Axial grids of 64, 128, 256, and 512 points show changes consistent with second
order convergence. The conservation residual decreases approximately fourfold
per grid doubling. Local trace residuals are at floating-point roundoff. At
256 points, angular cutoffs 2, 4, 6, 8, and 10 are varied separately. The high
angular differences eventually reach subtractive SVD roundoff; those tiny
signed values are not interpreted as a proved tail sign or rigorous error
bound. The stencil rejects \(\max h\kappa q/r\ge2\).

Constant-radius controls at \(L/a=1,2,4\), with sufficiently many angular
sectors for each ratio, recover the full-angular Bessel energy **and all three
stress components** at second order. The radial and angular truncations are
numerical regulators being removed. In particular, a cutoff on this spatial
Hamiltonian is not relabeled as the physical covariant finite \(\Lambda\) for
a nonconstant lapse.

## What this says about self-sourcing

The geometric calculation requires
\(8\pi G(\rho+p_x)_{\rm required}|_{x=0}=-2\), and a negative integrated
null source, for an Einstein-equation realization of these imposed profiles.
The positive differences in the table mean that changing from P to AP lowers
the computed null source by those amounts. This is a state/domain-induced
change in the direction needed at the neck. It does not give the absolute
\(\langle T_{\mu\nu}\rangle_{\rm AP}\), and the geometric and quantum columns
have different dimensional coefficients: no \(G\) or subtraction coefficient
has been fitted to equate them.

Both calculations here use a **compact axial \(S^1\)**. This is distinct from
the unwrapped infinite periodic throat array and from compactifying the
carrier's transverse \(y\) direction. The actual domain and spin structure of
the proposed parent/child solution still require derivation.

The smallest remaining absolute-stress calculation is to vary the adopted
**four-dimensional covariant proper-time determinant and its complete
compensator** in a specified AP state on this same smooth metric, retaining
lapse, radial, and sphere variations. The Dirac state in the present control
is already specified; its absolute ultraviolet normalization and the finite
local terms of the common action are not fixed by the spin difference.
Unlike on the constant-radius product, the varying-radius geometry does not
force all local counterterms to have zero axial-null projection. Their
coefficients must come from the adopted ultraviolet/measure prescription and
common action, not from demanding that the chosen neck pass. No additional
scalar field is needed for this next calculation, and none has been added here.

Only the full absolute tensor, its constraints and Ward identities, and the
same action's metric equations can establish self-sourcing. This result earns
the narrower claim that the same Dirac operator now supplies a computed local
quantum stress **change on the actual varying geometry**, with independent
normalization, derivative, and conservation checks.

Reproduce with `python3 scripts/check_nsc_shape_response.py --check`.
New records use exclusive `--output PATH`. Every recorded field is checked:
floating values use `atol=2e-7, rtol=1e-6` to accommodate subtractive SVD/BLAS
roundoff; integers, scope booleans, equations, keys and source/input hashes are
exact. These reproduction tolerances are separate from measured grid changes,
the finite-variation tests, and the reported conservation residuals.
