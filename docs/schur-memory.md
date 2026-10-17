# Schur memory and the local effective Lambda target

For $D=\begin{pmatrix}D_p&\Phi\\\Phi^\dagger&D_c\end{pmatrix}$, eliminating
the exterior gives $\Sigma(z)=\Phi(z-D_c)^{-1}\Phi^\dagger$ and the interior
resolvent $(z-D_p-\Sigma)^{-1}$. The corresponding causal memory carries an
exterior initial-state drive and cross correlations. The intended physical
target is a local effective $\Lambda$ response from this unresolved exterior
structure. The finite identities, spectral controls and geometric candidates
below supply its retained dictionary; the physical identification of
$\Sigma$ with the actual exterior remains open.

## Step 7: the current candidate's first vacuum projection

The current conformal-scalar cap–wall candidate gives the exact smooth-bulk
result

$$
\boxed{\Delta\Lambda_{\rm classical}\,L_{\rm ref}^2=0.}
$$

This is the additional trace projection of its classical scalar stress. The
cap already contains the supplied term $\Lambda_{0,-}\ell^2=3(H\ell)^2=3$;
the exterior contains $\Lambda_{0,+}=0$. Those baseline values were inputs
to the geometry. The zero above holds for any positive reference length,
so it requires no choice of physical ruler or observed cosmological value.
For displaying the baseline, $L_{\rm ref}=\ell$ uses the existing reference
units; the P001 value $\ell=10000\,\mathrm m$ remains a fiducial scale.

Reuse the [owned action and trace identity](nsc-atlas-cap-wall-response.md),
with $\Lambda_0=\kappa_4 V_0$, rather than extending the cap computation:

$$
T^\phi_{\mu\nu}=\nabla_\mu\phi\nabla_\nu\phi
-\tfrac12g_{\mu\nu}(\nabla\phi)^2
+\tfrac16\left[G_{\mu\nu}\phi^2+
(g_{\mu\nu}\Box-\nabla_\mu\nabla_\nu)\phi^2\right],
$$
$$
T^\phi{}^\mu{}_\mu=\phi(\Box-\mathcal R_J/6)\phi=0,
\qquad \mathcal R_J=4\Lambda_0,
\qquad
\Delta\Lambda_{\rm classical}L_{\rm ref}^2
=\tfrac14(\mathcal R_J-4\Lambda_0)L_{\rm ref}^2=0.
$$

This is exact on shell, including nonzero classical pulses. The pulse still
has nonzero stress, mass transport and wall response: a zero trace does not
make its stress tensor zero. In the decomposition
$T^\phi_{\mu\nu}=-(\Delta\Lambda/\kappa_4)g_{\mu\nu}
+T^{\rm TF}_{\mu\nu}$ with $g^{\mu\nu}T^{\rm TF}_{\mu\nu}=0$, the
vacuum projection is zero. This calculation concerns the unaveraged smooth
bulk. It does not assign a four-dimensional vacuum density to the material
wall, perform nonlinear metric averaging, or include quantum stress.

The [executable calculation](../src/recursive_horizons/nsc_atlas_effective_vacuum.py)
contracts the full improved stress and checks the result independently in
the Einstein frame. With $q=\Psi/\sqrt6$, $F=\operatorname{sech}^2q$ and
$\widehat V=\Lambda_0\cosh^4q$, use
$\mathcal R_E=6(\nabla_Eq)^2+4\widehat V$ and
$\Box_Eq=\widehat V_{,q}/6$. Then

$$
\mathcal R_J=F\left[\mathcal R_E+3\Box_E\log F
-\tfrac32(\nabla_E\log F)^2\right]=4\Lambda_0.
$$

The apparent field-dependent Einstein-frame potential therefore supplies
no extra Jordan-frame vacuum curvature. Seven
[checks](../tests/test_nsc_atlas_effective_vacuum.py) also reconstruct this
curvature from the actual nonlinear bulk operator's metric and field jets,
in both regions at zero and nonzero amplitude. The exact result has no
quadrature error; these checks do not estimate errors in an evolved solution.
Reproduce with
`python scripts/run_python.py -m recursive_horizons.nsc_atlas_effective_vacuum`
and `python scripts/run_python.py -m pytest tests/test_nsc_atlas_effective_vacuum.py -q`.
No endpoint replay, new spatial grid or observational fit enters this calculation.

For Step 7 this resolves the classical scalar route to an additional smooth
vacuum term in the current candidate. A nonzero contribution needs an
explicit additional mechanism, such as a specified quantum effective action
or a defined geometric averaging operation. The existing finite Schur
matrices have not been derived from this candidate's exterior. Its physical
response operator, state and metric variation must be supplied before those
matrices can yield an effective stress. In particular, the zero-frequency
coefficient of a scalar self-energy multiplies a field quadratic term;
it is not by itself a metric volume term or a cosmological constant.
The older two-sheet and Gauss–Bonnet results below retain their own actions.

## The transition and vacuum-form coefficient

The frozen positive-$F$ DEV3 transition uses radius $32/25$, expansion
magnitude $225/64$, thickness $1/50$ and endpoint $F=4$. Its coefficients obey

$$
c_{\rm transition}=-\frac{18}{1015},\qquad
v_{\rm room}=-6c_{\rm transition}=\frac{108}{1015},\qquad
\alpha_{\rm GB}/L_*^2=-\frac1{8c_{\rm transition}}=\frac{1015}{144},\qquad
\alpha_{\rm GB}v_{\rm room,physical}=\frac34.
$$

The explicit two-sheet carrier is

$$
ds_5^2=e^{-36y^2/1015}
\left[A(\rho)dt^2-\frac{d\rho^2}{A(\rho)}-(\rho^2+1)d\Omega_2^2-dy^2\right].
$$

Both areal branches are needed: one radius can have distinct parent and
child metric values, joining at the throat. Pure Einstein gravity requires
a negative null source for this carrier; a canonical positive-kinetic scalar
alone cannot provide it. The transition-fixed Gauss–Bonnet coefficient
cancels that local null projection geometrically. On the room, its full
component identity is

$$
[G^A{}_B+\alpha_{\rm GB}H^A{}_B]_{y=0}
=\operatorname{diag}\left(\frac{108}{1015},\frac{108}{1015},
\frac{108}{1015},\frac{108}{1015},
-\frac{R_4+(1015/144)\,GB_4}{2}\right).
$$

Thus the four tangential components have vacuum form and the normal component
contains the room curvature action. The same $\alpha_{\rm GB}=1015/144$ lies
at the five-dimensional EGB critical point: with
$\Lambda_0=-108/1015$, $\alpha_{\rm GB}\Lambda_0=-3/4$ and the vacuum
discriminant vanishes. Kinetic and characteristic conclusions must therefore
use the actual two-sheet background. Bulk-source microphysics and full
stability remain open.

The expanding-child Ricci sign selects positive $\Lambda_{\rm effective}$
in the retained doubled-FLRW convention $R=-2\Lambda_{\rm effective}$.
It selects $\Phi^2/\Lambda^2=1006/1015$. The retarded boundary pole then has
positive weights and subluminal group velocity, and its propagating band
ends at $p^2+\Phi^2=\Lambda^2$. The squared bandwidth is $9/1015$, half the
transition-curvature magnitude. This boundary excitation is not an electron
identification or the full graviton propagator.

The imported RS2 black-universe equation uses
$G^\nu{}_\mu=-\Lambda_4\delta^\nu{}_\mu-\kappa_4^2T^\nu{}_\mu
-\kappa_5^4\Pi^\nu{}_\mu-E^\nu{}_\mu$ and
$\Lambda_4=(\Lambda_5+\kappa_5^4\lambda^2/6)/2$.
The outside-kernel variation projects to $-E^\nu{}_\mu$ in that convention.
It matches the black-universe exterior response dictionary; the bulk metric
and its healthy response spectrum still have to determine that projection.

## Spectral and boundary controls

Room geometry forces originate in the regulated fermionic spectrum and its
scale anomaly. The compact allocation separates one massless Dirac field,
two fields at each $m_n=n\pi/\ell$ for $n\geq1$, and the warped complement.
It supplies a Euclidean subtraction, not a real-time state. With
$J=\ell_{\rm proper}/\ell$, the compact light limit is

$$
\lim_{u\downarrow0}[h_\Lambda(u)-E_1(u/\nu^2)]
=2\log\frac{\Lambda J}{\nu}
+2\sum_{n\geq1}E_1\!\left[\left(\frac{n\pi}{\ell J\Lambda}\right)^2\right].
$$

The light singular slope is $1/J^2$ and its canonical mass stays zero.
The imported spectral bulk/boundary constraint is $c_K=2c_R$, leaving one
coefficient for $c_R\int R+c_K\int K$. This universal ratio differs from
the particular transition coefficient $-18/1015$.

For the exponential heat profile,
$h(z)=\int_0^1e^{-\alpha(1-\alpha)z}d\alpha$ and
$F_2(z)=(1+z/4)h(z)-2$. Its nearest complex zeros lie beyond the unit cutoff
disk. Canceling the stationary vacuum background gives
$K_2(z)=F_2(z)+1$, $K_2(0)=0$, $K'_2(0)=1/12$ and a positive massless
propagator residue $12$, with no additional pole in that disk. These flat
and Euclidean controls leave the curved retarded positivity test open.
Reciprocal room symmetry requires $B(\chi)=B(-\chi)$: for
$B=b_0+b_2\chi^2+b_4\chi^4$, the stationary linear scalar mixing is zero,
while $B''(0)=2b_2$ retains the nonlinear transition channel.

The first-order curvature EFT change of inverse metric is

$$
\delta g^{\mu\nu}=-\frac{2c_W}{A_4}
\left[R^{\mu\nu}-\frac{R}{6}g^{\mu\nu}
+\frac{T^{\mu\nu}-Tg^{\mu\nu}/3}{2A_4}\right]
+\frac{c_R}{A_4}\left(R-\frac{T}{2A_4}\right)g^{\mu\nu}.
$$

Here old inverse metric equals new inverse metric plus this shift. Boundary
data, observables and the quantum Jacobian must follow the same perturbative
map; it is not a fundamental fourth-order theory.

## Finite state and energy controls

The coordinate generator is gravity plus $M\operatorname{Tr}(CH)$ with
$M=4\kappa$ counted once. Nodal force partials must be divided by node spacing
before conversion to stress densities. With $N=rL$, $q=rQ$ and
$V=4\pi qr^2$, the unit normal is $N^{-1}(\partial_t-\beta\partial_x)$,
whose dual is $Ndt$. Coordinate energy
$e=L F_L+\beta F_\beta=NV\rho-\beta qVj$ differs from normal-observer energy.
The latter obeys

$$
\partial_t(V\rho)+\partial_x[V(Nj/q-\beta\rho)]
=-NV(p_rK_r+2p_\perp K_\perp)-V(j/q)N_x.
$$

Independent metric variations control this identity and the chart-work
ledger. Mode-projection leakage is extra exchange; the controls establish
neither a continuum limit nor renewal.

The frozen six-mode incoherent imbalance witness uses the Cesàro projection
onto the commutant, with finite-time remainder
$\|C_T-\bar C\|_F\leq2\|Y\|_F/(\gamma T)$. Its retained mean populations are
$(1.3837895948404104,1.0058995585808204,0.6103108465787691)$;
the channel current is $0.0012511140479781876$. Reprojection moves the
response by $0.004920298239642384$, whereas the frozen shortcut gives
$0.00030421076749148856$. These finite controls do not close the physical
local gate. The smooth-motif Weyl-disk tail also controls a projected
resolvent only; it is not identified with a dark-energy density or a cosmic
fraction. Original calculations are at
[Zenodo](https://doi.org/10.5281/zenodo.23090015).

The retained coordinate Hamiltonian and normal-observer stress API lives in
[nested.energy_exchange](../src/recursive_horizons/nested/energy_exchange.py).
It keeps the generator gravity energy plus `M Tr(C H)` with `M=4 kappa`
counted once, the finite regional controls, and the separate normal-observer
shell energy. Renewal and a continuum limit remain unclaimed.
