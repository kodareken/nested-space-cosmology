# The next stationary equation: the covariant measure and vacuum response

The [calculation](../scripts/check_nsc_covariant_measure.py) and
[record](../results/nsc-4-covariant-measure.json) make two missing parts of the
common action explicit: its finite Weyl-invariant part and the meaning of its
scale-field measure. They also provide a full-angular, four-dimensional vacuum
response against which the next metric variation can be checked. No coefficient
is chosen to produce a stationary root.

## The anomaly leaves a physically relevant shape equation undetermined

The current finite regulated prescription fixes its matrix calculus. It does
not yet specify the continuum covariant determinant, quantum state, or complete
compensator. The primary spectral-anomaly derivation distinguishes cutoff from
normalization, separates invariant and anomalous factors, and states that the
factorization is not unique. [Andrianov, Kurkov and Lizzi, §§3–4](https://arxiv.org/html/1106.3263v1)

For the smooth product metric

\[
ds_E^2=d\tau^2+(L/2\pi)^2dx^2+a^2d\Omega_2^2,
\qquad x\sim x+2\pi,
\]

the implementation constructs Christoffel symbols and contracts the Riemann
tensor directly. It obtains

\[
\mathcal R=2/a^2,\quad R_{\mu\nu}^2=2/a^4,\quad
R_{\mu\nu\rho\sigma}^2=4/a^4,\quad
C^2=4/(3a^4),\quad E_4=0.
\]

Consider the difference between two otherwise unspecified completions,
\(\Delta\Gamma=\alpha\int\sqrt g\,C^2\). In four dimensions this term is
invariant under local Weyl transformations: the volume and squared-Weyl-tensor
weights cancel pointwise. Nevertheless, per unit Euclidean time,

\[
\Delta E=\frac{16\pi\alpha L}{3a^2},\qquad
\partial_L\Delta E=\frac{16\pi\alpha}{3a^2}.
\]

At fixed spatial volume \(V=4\pi a^2L\), the difference becomes
\(64\pi^2\alpha L^2/(3V)\), with nonzero spacing derivative. Thus preserving
volume does not eliminate this shape ambiguity. A hypothetical nondegenerate
root of a reference energy obeys

\[
\left.\frac{dL_\star}{d\alpha}\right|_{\alpha=0}
=-\frac{16\pi}{3a^2 E_{0,LL}(L_0)}
\]

at fixed \(a\). No numerical root or value of \(\alpha\) is selected here.
This demonstrates what anomaly matching alone leaves unresolved. A fully
specified ultraviolet functional could fix this finite invariant part.
Properly matched observables remain unchanged when renormalization scheme and
couplings are transformed consistently; this is not a scheme-dependent physical
prediction. The missing invariant functional is explicit in the primary
anomaly-induced-action construction. [Fabris, Pelinson and Shapiro, §2](https://arxiv.org/html/gr-qc/9810032v2)

In a four-dimensional torsionless geometric sector, a useful finite local basis
is \(M^4,M^2\mathcal R,C^2,\mathcal R^2,E_4,\Box\mathcal R\). The Euler and
total-derivative terms do not independently vary in the bulk on this closed
fixed-topology geometry. Conformal assumptions constrain the permitted terms;
the \(C^2\) ambiguity survives even exact Weyl invariance. Dynamical link,
gauge, and internal fields require their own covariant terms. The basis lists
what the same measure must determine, not independently weighted additions to
the already induced action. Covariant conserved Dirac stress and its
renormalization are developed by [Dappiaggi, Hack and Pinamonti](https://arxiv.org/abs/0904.0612).

## The finite flat scale average diverges

Using the unchanged [regulated operator](../src/recursive_horizons/nsc_regulated.py),

\[
\Gamma(D)=\tfrac12\operatorname{Tr}E_1(D^2/\Lambda^2)
+N[\log(M/\Lambda)+\gamma_E/2],\quad D_\varphi=e^{-\varphi}D,
\]

every nonzero fixed eigenvalue contributes an \(E_1\) term tending to zero as
\(\varphi\to-\infty\). Therefore \(e^{-\Gamma(D_\varphi)}\) tends to a positive
constant, and its flat integral over the real \(\varphi\) line diverges.

For \(D=\operatorname{diag}(1,-1)\), \(M=\Lambda=2\), the integrals over
\([-5,0],[-10,0],[-20,0]\) are respectively 2.615559345854003,
5.422856763688429, and 11.037451599357281. Their asymptotic slope is
\(e^{-\gamma_E}=0.5614594835668851\).

This is a finite fixed-\(N\) result with a flat measure, not a continuum no-go
theorem. A gauge volume or physical dilaton measure must be specified before the
displayed averaging can define a stationary action. Finite integration bounds
would be additional inputs. If common Weyl rescaling is instead an exact gauge
symmetry, with \(\delta g_{\mu\nu}=2\sigma g_{\mu\nu}\) and
\(\delta\varphi=-\sigma\), its Ward identity is

\[
2g_{\mu\nu}\frac{\delta\Gamma}{\delta g_{\mu\nu}}
-\frac{\delta\Gamma}{\delta\varphi}=0.
\]

The common-scale equation is then dependent. A relative parent–child scale may
remain physical, but requires its own linked metric/field variation.

## A complete angular and frequency control for the next stress calculation

Use natural units \(\hbar=c=1\), one complex four-component massless Dirac
field, and the static zero-temperature ground state on
\(\mathbb R_t\times S^1_L\times S^2_a\). There is no inserted fermion mass.
The two declared spatial spin structures are periodic (\(\eta=0\)) and
antiperiodic (\(\eta=1/2\)); a thermal time condition is not being substituted
for the spatial condition. This is the UV-finite continuum one-loop
compactification contribution. No finite physical spectral cutoff is imposed
in these Bessel/Fermi reference calculations. A separate calculation below
retains the proper-time cutoff and quantifies the difference from this limit.
The earlier scale-average control also retains its explicitly fixed regulator.

The sphere eigenvalues are \(\pm k/a\), with \(2k\) eigenspinors per sign.
This is the standard spinor-sphere spectrum. [Camporesi and Higuchi](https://arxiv.org/abs/gr-qc/9505009)
Tensoring the sphere with the two-component \((t,x)\) Dirac factor produces
spatial energies \(\pm\sqrt{p^2+k^2/a^2}\), each with multiplicity \(4k\).
Consequently \(-\operatorname{Tr}|H|/2=-4k\sqrt{p^2+k^2/a^2}\) per circle
momentum. Explicit Hamiltonian blocks and the rank-four spatial Weyl law check
this determinant counting independently.

Removing the zero-winding/decompactified-cylinder term and integrating the time
frequency gives

\[
E_\eta=\frac8{\pi a}\sum_{k\ge1}k^2\sum_{n\ge1}
\frac{\cos(2\pi n\eta)}n K_1(nkL/a).
\]

The difference \(\Delta E=E_0-E_{1/2}\) compares the same local metric, so local
counterterms cancel. An independently integrated Gaussian representation is

\[
\Delta E=\frac{4L}{\pi}\int_0^\infty\frac{dt}{t^2}
\left[\sum_{k\ge1}k e^{-tk^2/a^2}\right]
\left[\sum_{n\ge1,\ n\text{ odd}}e^{-n^2L^2/(4t)}\right].
\]

This formula includes both frequency and circle Poisson factors. The Bessel
derivative is checked against centered differences of this independent heat
integral. The pressure difference is
\(\Delta p_x=-\partial_L\Delta E/(4\pi a^2)\), at fixed \(a\). Energy has units
of inverse length and pressure inverse length to the fourth power.

For \(L/a=4\), the calculation gives
\(a\Delta E=0.06685499963800426\) and
\(a^2\partial_L\Delta E=-0.0798428968284928\).
Angular and winding truncations are varied independently; quadrature errors
are distinguished from truncation differences. No rigorous infinite-tail bound
is asserted. Positive odd windings and \(K_1'(z)<0\) prove
\(\Delta E>0\) and \(\partial_L\Delta E<0\) for every \(L,a>0\).

These are vacuum/domain response controls. They do not determine the full
absolute vacuum stress tensor, the physical spin structure, or a gravitational spacing. The
Casimir method, Weyl anomaly, and sphere spectrum are established physics; the
contribution here is their reproducible application to the missing closure
equations and normalization checks.

## A negative radial-null source from the same quantum field

The finite compactification energy determines its complete relative diagonal
stress on this constant-radius product:

\[
\rho_\eta=\frac{E_\eta}{4\pi a^2L},\qquad
p_{x,\eta}=-\frac{E_{\eta,L}}{4\pi a^2},\qquad
p_{\perp,\eta}=-\frac{E_{\eta,a}}{8\pi aL}.
\]

Here \(p_\perp\) is the pressure in either sphere direction, not their sum.
The radius derivative is calculated from the spectrum, checked against an
independently differentiated integral, and checked against a finite change in
radius. For an orthonormal radial null vector \(k^a=(1,1,0,0)\), the projection
is \(T_{ab}k^ak^b=\rho+p_x\).

There is a useful exact sign result. The antiperiodic energy also has the
independent Fermi-log representation

\[
E_{\rm AP}=-\frac8\pi\sum_{k\ge1}k\int_0^\infty dp\,
\log\!\left(1+e^{-L\omega_k}\right),\qquad
\omega_k=\sqrt{p^2+k^2/a^2}.
\]

The logarithm is positive, and

\[
E_{{\rm AP},L}=\frac8\pi\sum_{k\ge1}k\int_0^\infty dp\,
\frac{\omega_k}{e^{L\omega_k}+1}>0.
\]

Consequently \(E_{\rm AP}<0\) and \(\rho_{\rm AP}+p_{x,\rm AP}<0\) for every
positive \(L,a\). The positive periodic Bessel series instead has
\(E_{\rm P}>0\), \(E_{{\rm P},L}<0\), and positive radial-null stress. The
Fermi-log quadrature independently agrees with the Bessel energy and both metric
derivatives at \(L/a=1,2,4\).

| \(L/a\) | \(a^4(\rho+p_x)_{\rm P}\) | \(a^4(\rho+p_x)_{\rm AP}\) |
|---:|---:|---:|
| 1 | 1.6994660879849266 | -1.507047697114369 |
| 2 | 0.09620919543527895 | -0.08845021288528555 |
| 4 | 0.003879491274996693 | -0.003804242533368134 |

Unlike the separate density and pressures, this radial-null projection is
unambiguous even before the local geometric coefficients are fixed. Every
locally covariant metric counterterm tensor on the constant-radius product
respects boosts in its locally flat \((t,x)\) factor. Its two-dimensional block
is therefore proportional to the Minkowski metric and has zero radial-null
contraction. This is a statement about local tensors; the compactification and
its vacuum response break the global boost symmetry. The decompactified
reference ground state is itself boost invariant in \((t,x)\), so its null
projection also vanishes. The calculated compactification projection is thus
the absolute radial-null projection of this specified compact ground state.

The earlier ambiguous \(\alpha C^2\) term demonstrates the distinction:

\[
\rho_{\alpha C^2}=-p_{x,\alpha C^2}
=p_{\perp,\alpha C^2}=\frac{4\alpha}{3a^4},
\qquad (\rho+p_x)_{\alpha C^2}=0.
\]

For the finite massless compactification contribution,
\(LE_L+aE_a=-E\), and the independently varied stress obeys
\(-\rho+p_x+2p_\perp=0\). This does not discard the Dirac trace anomaly: the
identical local anomaly cancels when the decompactified cylinder is subtracted.

The compact direction here is axial \(x\), and the negative null projection
points along that direction. This topology differs both from the **unwrapped**
periodic throat array and from compactifying the carrier's **transverse** \(y\)
direction. They cannot be identified by calling each setup periodic.

This supplies a constructive physical connection. The ordinary quantum Dirac
field in a declared compact state produces a negative null source of a type
that a flaring throat can require in an Einstein-equation realization. No negative
kinetic-energy field or fitted stress is introduced. The cylinder is not a
flaring throat and is not a self-consistent gravitational solution. Its stress
must not be inserted into a varying-radius neck as though it were exact there;
the next calculation must derive that geometry's actual domain, topology,
state-dependent response, and all metric equations with the same action. The result is a benchmark for
the mechanism, not a discovery of Casimir stress or a completed nesting model.

## The negative axial-null sign survives the declared proper-time cutoff

The continuum result connects to the project's regulator by retaining its
lower proper-time limit after the same frequency integration and
compactification subtraction:

\[
E_\eta(L,a;\Lambda)=\frac{2L}{\pi}\int_{\Lambda^{-2}}^\infty
\frac{dt}{t^2}K_s(t)S_{0,\eta}\!\left(\frac{L^2}{4t}\right),\qquad
K_s(t)=\sum_{k\ge1}k e^{-tk^2/a^2},
\]

where \(S_{j,\eta}(b)=\sum_{n\ge1}n^j\cos(2\pi n\eta)e^{-bn^2}\).
The zero-winding term is removed at the same cutoff. This is the specified
proper-time modulus's compactification response, without an independently
weighted action or an asserted complete anomaly/compensator.

Holding dimensional \(\Lambda\) fixed while varying \(L\), the term from the
explicit prefactor in \(E_L\) cancels the density in the null projection:

\[
\rho+p_x=\frac{L^2}{4\pi^2a^2}\int_{\Lambda^{-2}}^\infty
\frac{dt}{t^3}K_s(t)S_{2,\eta}\!\left(\frac{L^2}{4t}\right).
\]

For AP spin structure, set \(q=e^{-b}\). Then
\(S_{2,\rm AP}=q\partial_q\theta_4(0,q)/2<0\): the positive product
\(\theta_4(0,q)=\prod_{m\ge1}(1-q^{2m})(1-q^{2m-1})^2\)
strictly decreases on \(0<q<1\). Its logarithmic derivative is an absolutely
convergent sum of negative terms on compact subintervals. The product identity
is [DLMF 20.5.4](https://dlmf.nist.gov/20.5#E4).
All other factors in the null integral are positive. Consequently the AP axial
null projection is negative for every \(\Lambda,L,a>0\) in this prescription.
Periodic \(S_2\) has positive terms and gives the opposite sign.

The implementation uses direct winding sums for \(b\ge1\) and differentiated
Poisson-dual theta identities for \(b<1\). The latter avoids alternating
cancellation when proper time is large. For example, at \(b=0.02\) it resolves
\(S_{2,\rm AP}=-2.030214246791232\times10^{-49}\), agreeing with an
independent 100-digit theta derivative to relative error below \(10^{-12}\).
For much smaller \(b\), floating-point underflow can give zero; the exact sign
does not rely on such floating-point values.

| \(L/a\) | \(\Lambda a=1\) | \(\Lambda a=2\) | \(\Lambda a=4\) | Continuum |
|---:|---:|---:|---:|---:|
| 1 | -0.00000612338 | -0.06447565 | -1.12415939 | -1.50704770 |
| 2 | -0.00313775 | -0.06511492 | -0.08844858 | -0.08845021 |
| 4 | -0.00250995 | -0.00380414 | -0.00380424 | -0.00380424 |

Entries are \(a^4(\rho+p_x)_{\rm AP}\). The record contains full precision,
both spin structures, energy errors, pressure errors, and independent finite
differences of the regulated energy. Auxiliary \(\Lambda a=8,16\) calculations
check approach to the continuum; all reported \(\Lambda a=16\) null responses
agree with the Bessel reference within \(10^{-11}\). Angular truncations are
checked separately. The sign survives a low cutoff even when the magnitude is
strongly suppressed. The continuum proof using \(E_{{\rm AP},L}>0\) does not
extend to every finite cutoff: this derivative can become negative, while the
theta identity still fixes \(\rho+p_x<0\).

The metric derivatives keep \(\Lambda\) fixed. A joint unit change instead
rescales it inversely with lengths, giving
\(E+LE_L+aE_a-\Lambda E_\Lambda=0\). The finite-cutoff compactification trace
therefore need not vanish; the calculation checks this cutoff term explicitly.
Neither the surviving sign nor these numbers identify the compact axial source
with the original throat's radial source or determine the unknown full
covariant action. They connect the same Dirac response to the accepted
proper-time prescription on a precisely declared product geometry.

## The physical equation this unlocks

Next use a smooth periodic static cell
\(ds^2=-N(x)^2dt^2+q(x)^2dx^2+a(x)^2d\Omega_2^2\). Retain \(N,q,a\) through
variation, and compute the common action's metric and link equations before
fixing coordinate or Weyl gauge:

\[
\frac{\delta\Gamma_{\rm one}}{\delta N}=
\frac{\delta\Gamma_{\rm one}}{\delta q}=
\frac{\delta\Gamma_{\rm one}}{\delta a}=
\frac{\delta\Gamma_{\rm one}}{\delta\Phi}=0.
\]

Setting the lapse to one before variation loses the energy constraint. Check
the independently obtained stress against
\(p_x'+(N'/N)(\rho+p_x)+2(a'/a)(p_x-p_\perp)=0\).
The present calculation supplies only the constant-radius reference response,
not the varying-radius stress in those equations. It includes the complete
angular tower; one \(|\kappa|=1\) reduction cannot supply four-dimensional
stress. Equal-cell periodicity sets \(\Omega=1\); unequal coupled cells and
their metric matching are required to determine inheritance. This ultrastatic
four-dimensional control is separate from the trapped five-dimensional carrier.

Reproduce with `python3 scripts/check_nsc_covariant_measure.py --check`.
Every scientific field and every source hash is compared, with numeric
`rtol=atol=1e-8`; exact expressions, keys and scope booleans are exact. New
records use exclusive `--output PATH` creation and cannot silently replace an
existing result. The runner rejects simultaneous `--output` and `--check`.
