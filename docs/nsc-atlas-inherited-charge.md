# Internal rotation, horizon transport and the inherited charge state

The [smooth-horizon calculation](nsc-atlas-smooth-horizon.md) supplied a
regular geometry and an offset-current comparison, but its two equal-mass
field components had been introduced as a comparison. Here we derive a
candidate pair from the higher-dimensional real field already used in the
[volume reduction](nsc-atlas-scale-selection.md), then calculate a
conditional inherited quantum state and its charge distribution.
The [compact record](../results/summaries/nsc-atlas-inherited-charge-v1.json)
contains the numbers, assumptions and reproduction bindings.

**Decision and scope.** Three routes were compared: identifying the two
unequal homogeneous frequencies as charge components; transporting an
unspecified complex probe; and restoring internal harmonics omitted by the
volume truncation. The first is not a conserved charge pair, and the second
does not derive its field. The third supplies a concrete pair from geometry.
We complete the internal Einstein space as a round S12, preserving the
action's curvature coefficient. Roundness is an explicit shape assumption,
not a selection theorem. The parent is the existing stationary late vacuum;
this is a quadratic perturbation calculation, not a complete 16-dimensional
formation solution or a consistent finite-mode nonlinear truncation.

The likely failures are a hidden new shape/state assumption, confusing
internal charge with baryons, and treating a thermal suppression as an
asymmetry. They are tested separately below. The budget is exact algebra
and finite spectral arithmetic; no field grid, scan or observation fit.

## 1. The field pair comes from a real higher-dimensional field

The existing metric ansatz and potential are

$$
ds_{16}^2=e^{-12s}ds_{4E}^2+e^{2s}\gamma_{ab}dy^ady^b,
\qquad V(s)=7e^{-12s}-\frac{15}{2}e^{-14s}+\frac{25}{32}e^{-24s}.
$$

The Einstein-Hilbert curvature contribution is
-R_gamma exp(-14s)/(2 kappa4), with kappa4=1. Thus R_gamma=15.
For a round twelve-sphere of radius a0,

$$
a_0^2=\frac{12\cdot11}{15}=\frac{44}{5},\qquad
\lambda_l=\frac{l(l+11)}{a_0^2},\qquad \lambda_1=\frac{15}{11}.
$$

Let n_i be unit embedding coordinates in R13. The thirteen normalized real
l=1 harmonics are Y_i=sqrt(13/V_gamma)n_i, so integral Y_i Y_j=delta_ij.
Their normalized complex combinations obey

$$
Y_+=\frac{Y_1+iY_2}{\sqrt{2}},\qquad
\delta\Phi=xY_1+yY_2=\psi Y_+^*+\psi^*Y_+,
\qquad \psi=\frac{x+iy}{\sqrt{2}}.
$$

The same real higher-dimensional field therefore contains a complex
generated-parent mode. Rotation in the (n1,n2) plane gives charges +/-1.
Choosing six commuting rotation planes organizes the thirteen real modes
as six complex pairs plus one neutral component. Under a single selected
plane, only its pair is charged and the other eleven components are neutral.
The choice of plane is not selected by the round geometry.

This uses the Kaluza-Klein relation between internal isometries, harmonics
and lower-dimensional charges; [Duff, sections 2 and 4](https://arxiv.org/pdf/hep-th/9410046)
is the borrowed construction. Canonical reduction here gives

$$
S_2=-\int\sqrt{-g}\left(|D\psi|^2+M_1^2(s)|\psi|^2\right)d^4x,
\qquad M_1^2=\omega_0^2e^{-12s}+\frac{15}{11}e^{-14s}.
$$

The external kinetic warp factors cancel. The internal gradient term carries
exp(-14s), while the existing scalar/top-form mass carries exp(-12s).
Gauge fields arise from mixed metric components and must be included when
these charged modes backreact. This is an internal rotational charge,
not yet electric charge or baryon number. Standard field quantization is
imported for the state calculation; quantum probabilities are not derived
from geometry alone here.

At the retained (12,2) minimum and omega0=20:

| Quantity | Value in existing model units |
| --- | --- |
| s_star | -.035646782553 |
| Lambda_star | .220983024213 |
| H_star | .271405615400 |
| Homogeneous mass squared | 613.528021388795 |
| First harmonic mass-squared increment | 2.246132580192 |
| Charged pair mass squared | 615.774153968986 |
| Charged pair mass | 24.814797076925 |

No SI mass or observational input enters these numbers. The homogeneous
volume and matter modes remain unequal; the new degeneracy comes from
internal rotational symmetry, not from setting those frequencies equal.

## 2. The flux does not select a sign on this round geometry

Flipping only n2 would reverse internal orientation. Instead take the
proper rotation P that flips n2 and n13: det(P)=+1. It conjugates the
rotation generator, P J12 P^-1=-J12, and exchanges Y_+ with Y_+^*.
It preserves the round metric, homogeneous scalar, and oriented top-form
flux. Consequently, the supplied nonzero flux does not break this
charge-exchange symmetry. The ancestor count j=2 supplies no internal
angular momentum or chemical potential in the current action/state.

Without a membrane current, q is continuous through the smooth horizon.
For a neutral exterior with constant s, both charge components obey the
same Klein-Gordon operator. In advanced EF coordinates, put
psi=exp(-i omega v)u(r)Y_lm/r. Then

$$
f u''+(f'-2i\omega)u'
-\left(\frac{f'}r+M_1^2+\frac{l(l+1)}{r^2}\right)u=0.
$$

At a simple root r_h, the indicial powers are p=0 and
p=2i omega/f'(r_h). The regular branch has a free complex amplitude.
Horizon regularity does not specify its charge direction or occupation.
This is an explicit transfer equation, not a phase imposed at a wall.

## 3. An invariant preparation and its complete local charge distribution

On the stationary de Sitter parent, impose the regular Euclidean Gaussian
preparation with internal rotation invariance. The continued S4 operator
has eigenvalues

$$
H_*^2 L(L+3)+M_1^2>0,\qquad L=0,1,2,\ldots.
$$

Its inverse is unique with those conditions. The two real components have
equal covariance and zero cross covariance; their complex Gaussian has
uniform phase. There is no coherent ellipse offset to extract. This is a
specified invariant preparation rule, not a derivation of the universe's
entire wavefunction from its action. De Sitter invariance alone is weaker
than Euclidean/Hadamard regularity: see [Allen](https://doi.org/10.1103/PhysRevD.32.3136).
A de Sitter Euclidean state must also not be silently installed as a
global equilibrium state on a generic two-horizon Kottler spacetime.

In the static patch this state has temperature T=H_star/(2pi)=.043195545274
with respect to the central clock, following
[Gibbons--Hawking](https://doi.org/10.1103/PhysRevD.15.2738). At each fixed
stationary frequency omega, define z=exp(-2pi omega/H_star). The two charge
occupations have equal geometric distributions. For their difference Q,

$$
P(Q=k)=\frac{1-z}{1+z}z^{|k|},\qquad
\langle Q\rangle=0,\qquad
\operatorname{Var}(Q)=\frac{2z}{(1-z)^2}.
$$

Thus zero expectation does not mean every local measurement has zero charge.
For the illustrative central-clock gap omega=M1, the log detailed-balance
ratio is -574.475838182, the ratio is 3.223394725e-250, and the single-mode
charge RMS is 2.539052865e-125. This frequency is a declared comparison gap,
not a selected static mode or the whole population. A regional fluctuation
requires its mode weights, volume and photon count before it can become an
eta distribution. None of those is fitted to 6.1e-10.

## 4. Inheritance preserves the state symmetry

A neutral two-mode Bogoliubov map has

$$
a_{\rm out}=\alpha a_{\rm in}+\beta^*b_{\rm in}^{\dagger},\qquad
b_{\rm out}=\alpha b_{\rm in}+\beta^*a_{\rm in}^{\dagger},
\qquad |\alpha|^2-|\beta|^2=1.
$$

Direct expansion, including the pair cross terms, gives the operator identity

$$
N_{+,\rm out}-N_{-,\rm out}=N_{+,\rm in}-N_{-,\rm in}.
$$

More generally, a charge-exchange-equivariant transport channel maps an
invariant density operator to another invariant density operator. Induction
therefore gives zero regional expected charge at every nesting depth when
the internal symmetry preserves the region. An open region still has the
boundary flux from the previous chapter; this is not an assumption that its
charge is separately conserved. An inherited nonsymmetric state can carry
nonzero charge, but neutral transport does not select it from an invariant one.

Smooth charge-conjugate coherent test states illustrate the remaining
freedom: their external quadratic stresses agree and their charges have
opposite signs. Mixtures vary mean charge as (2p-1)Q0 with unchanged external
stress. This is a fixed-background statement. In a backreacting solution,
the gauge and mixed-metric constraints also track charge; identical external
stress is not identical complete higher-dimensional initial data.

## 5. Result against the requested imbalance

| Readout | Calculated result | Target | Miss |
| --- | --- | --- | --- |
| Mean geometric charge asymmetry of invariant preparation | 0 | 6.1e-10 comparison | -100% |
| eta if this charge is identified with baryons and photon density is positive | 0 | 6.1e-10 | -100% |
| Physical baryon/photon eta from the atlas | Not selected | 6.1e-10 | Undefined |

The new contribution is the geometrically derived charged pair, its horizon
operator, and an explicit inherited state/distribution. The requested nonzero
number has not been obtained. The missing connection is now specific:
an inherited internal current or charge-asymmetric interaction/preparation,
with its baryon/photon interpretation. A round invariant vacuum cannot supply
it merely through the position or motion of its horizon.

Eight focused checks verify the warp factors, harmonic charge and proper
rotation, covariant radial equation, positive Euclidean kernel, exact
Bogoliubov charge balance and local thermal variance. Symbolic identities
and independent operator reductions are used rather than a parameter scan.
This closes the invariant harmonic route; a later asymmetric preparation
must identify its source instead of assigning a phase to match the target.

## 6. An added spinor candidate and its actual parent modes

The new [source calculation](../results/summaries/nsc-atlas-valve-source-gate-v1.json)
asks which operator could supply the [valve's](nsc-atlas-valve-formation.md)
missing turnaround stress. The preceding real scalar and its harmonics
remain bosons. This candidate adds a complex massless bulk Dirac field;
fermionic quantization/statistics are imported, not derived from that scalar
or from the offset generator. No Weyl/Majorana projection, occupation or
baryon identification is silently imposed.

For the existing ds16²=exp(-12s)g4+exp(2s)gamma12 and round R_gamma=15,
the graded Dirac operator is

$$
\not\nabla_{16}=e^{6s}(\not\nabla_4-3\not ds)\otimes1
+e^{-s}\Gamma_{4,*}\otimes\not\nabla_\gamma.
$$

The rescaling Psi=exp(3s) sum(psi_n tensor eta_n) gives canonical external
kinetic terms. The internal term retains a four-dimensional chirality
factor; its convention-dependent mass phase must be matched in every
vertex. Its positive squared mass is unambiguous:

$$
\lambda_{n,\pm}=\pm\frac{n+6}{a_0},\quad
d_n=64\binom{n+11}{11}\quad\hbox{per sign},\qquad
m_n^2=\lambda_n^2e^{-14s}.
$$

[Camporesi--Higuchi, sections3.1 and5.2](https://arxiv.org/abs/gr-qc/9505009)
supplies the round spectrum/representations. Applying it to the existing
a0²=44/5 gives

$$
m_{F0}^2=\frac{45}{11}e^{-14s}
=3\,\Delta M_{B,1}^2,
\qquad \frac{m_{F0}}{H_*}=9.56443878784.
$$

The logarithmic radion derivative is -7/sqrt84, while the actual mass
derivative at the minimum is-1.98260737801. Each sign's64-dimensional
ground representation has32 weights+1/2 and32 weights-1/2 under J12.
Both Dirac signs give128 four-dimensional Dirac mode fields. Spectral
sign, gauge-charge sign and particle/antiparticle excitation are distinct.
The proper charge-exchange rotation used above also pairs these weights,
so the added field does not select a charge bias or eta.

The exact old gauge normalization gives the new relations

$$
\frac{g_{\rm KK}^2}{\kappa m_{F0}^2}=\frac{13}{36},\qquad
\frac{e_{\min}^2}{\kappa m_{F0}^2}=\frac{13}{144},
\quad e_{\min}=g_{\rm KK}/2.
$$

These are leading geometric internal-U(1) relations, not electromagnetic
or baryon identifications. This e_min is not the membrane four-form charge.
The absolute coupling still contains g_grav=kappa_nat/Lunit², which is
unselected. No SI mass, cosmological ruler or measured fraction is fixed.

Positive internal curvature also imposes a robust obstruction: the twistor
norm and Lichnerowicz identity give lambda²>=d R_min/[4(d-1)]. For d=12
and R_min=15 this is45/11. Round S12 saturates it; no compact untwisted spin
shape with that positive-curvature lower bound has a zero mode. Saturation
is a mathematical variational fact, not an inherited shape-selection law.
A gauge twist/boundary or different source changes the hypotheses.

### Spatial modes on the actual initial shared parent

The [spatial calculation](../results/summaries/nsc-atlas-valve-dirac-modes-v1.json)
uses the actual degree40 epsilon1 constrained metric, not flat space or a
homogeneous average. Set Omega=Psi40²/H so h=Omega² gamma3. Left quaternion
multiplication Z=(-X4,X3,-X2,X1) puts old X1 into radial Z4 while preserving
the existing right frame. Conformal covariance gives

$$
D_h(\Omega^{-1}f)=\Omega^{-2}D_\gamma f,
\quad D_\gamma f=\lambda_\Sigma\Omega f,
\quad \int\Omega|f|^2dV_\gamma=1.
$$

The exact round operator is D_gamma=i sigma_i D_i+3/2. Its square identity
(D_gamma-3/2)²+2(D_gamma-3/2)=-Delta supplies the independent round check.
All angular sectors are kept. Coupling orbital spin to J=k+1/2 reduces the
problem to small radial blocks and exact one-dimensional Gegenbauer
integrals. No physical angular averaging or fitted Gram correction is used.

At cutoff24, the first positive doublet has lambda_Sigma=.679256139175;
the negative partner agrees in magnitude. With the candidate mass above,
E=sqrt(lambda_Sigma²+m_F0²)=2.68324181601. Full physical residuals are
1.90e-8 and7.79e-8; the much smaller projected Ritz residual is not used
as the full-operator error. A compact-inverse norm estimate isolates the
first positive cluster for the defined finite metric; quadrature/rounding
and the approximation to the Einstein initial geometry remain explicit.

For the normalized k=0 profile, f=(F I-iG sigma.n)xi/sqrt(4pi), so
rho_h=Omega^(-2)(F²+G²)/(4pi). The actual overlap is

$$
I_4=\int\rho_h^2dV_h=.00456033592910,
\qquad V_{\rm eff}=I_4^{-1}=219.282091396.
$$

The actual parent volume is232.733148481. This effective volume is a
calculated mode overlap, not a free cell size or a physical occupation.
The round-spectrum, old-frame representation, angular-density and full
residual controls are separate. The calculation selects no incoming
fermion state or backreaction.

### Source-state gates: thermal and correlated cases

Minimal four-dimensional EC has L_contact=-3kappa A_mu A^mu/16 in the
declared -+++ convention with {gamma,gamma}=-2eta. The
[corrected action](https://arxiv.org/html/2008.11686v2), equation2 at zero
nonminimal couplings, supplies this sign. The current two real scalars
have no such spin current. A four-dimensional axial-only comparison is
also not automatically the full16D torsion reduction.

In the explicitly vacuum-Wick-ordered, free, unpolarized, mu=0 thermal
benchmark, the literal four-component contraction gives A²=-m²T⁴/36
at T>>m. Its pressure correction is+kappa m²T⁴/192; total free plus
interaction NEC is positive. Its contact/free NEC ratio is
.09143871896 g_grav for this gap, in the controlled regime
g_grav T_hat²<<1 and m_hat<<T_hat. The half-Wick diagram convention in
[Diakonov--Tumanov--Vladimirov, equations64--76](https://arxiv.org/html/1104.2432v2)
has an explicit extra1/2; mixing conventions changes that coefficient.
None of these calculations chooses a temperature or fixes g_grav.

For the full complex massless16D Dirac field, retain both internal signs
and all transformed current vertices. The coincident addition theorem is
P_n,plus(y,y)=P_n,minus(y,y)=d_n I64/(64V_gamma). In the further specified
unpolarized nonrotating rest-frame ideal thermal state, equal plus/minus
occupations cancel the coincident thermal covariance difference. Its
fully vacuum-Wick-ordered full torsion contact is zero. This does not
declare the gas empty, the actual evolving state thermal, or an unspecified
renormalized vacuum stress zero. Derivatives of the two-point function
still carry energy. Exact matching, rather than a partial mode/vertex
truncation, is required.

The original MMP massless-support solution is also not automatically
available. Even if a magnetic RN sector were supplied, its normalization
r_e²=pi Q²G/e_min² and the above charge/mass relation imply

$$
m_{F0}r_e=\sqrt{18/13}\,|Q|\geq1.17669681.
$$

Thus this untwisted KK field has no parametrically massless throat limit.
This is a compatibility gate for [MMP's massless field construction](https://arxiv.org/abs/1807.04726),
not an exclusion of every massive or twisted-bundle source. Magnetic
Maxwell flux, a return loop and its state are additional missing data;
four-form flux is not their replacement.

### A connected four-point witness and an actual geometric kernel

Free energy/charge and quadratic covariance do not determine the contact.
An exact16-state CAR construction compares a one-particle/one-antiparticle
singlet with an isotropic triplet mixture. Both have N=2,Q=0, covariance
diag(1/2), zero anomalous covariance, equal free energy, and rotation/C
invariance. In one rest orbital,

$$
V^2\langle:A^2:\rangle_{S}=-10,\qquad
V^2\langle:A^2:\rangle_{T}=2.
$$

Their quartic difference12/V² survives common constant/quadratic vacuum
subtractions. The matched contact Hamiltonian favors the singlet by
9kappa/(4V) in that fixed sector. This is not the full interacting vacuum
or a derived equation of state; no curved-space mode was inferred from
that rest-orbital example.

The actual parent mode above removes that arbitrary volume and includes
its spatial kinetic energy. In a positive-mass basis its pointwise field
is sqrt(rho_h)*(Q_tau tensor U_spin)*psi_rest. Here Q mixes particle
components by E=sqrt(m²+lambda_Sigma²), while U_spin is a local SU2 rotation.
Spatial axial currents rotate; the temporal current becomes
(m/E)A0+(lambda_Sigma/E)(N-2)rho_h. Both states lie at N=2, so their full
normal-ordering difference is [8+4(m/E)²]rho_h². Consequently the bare
externally indexed torsion channel has the actual geometric kernel

$$
\Delta H_A=\frac{3\kappa}{4}[2+(m/E)^2]I_4,
\qquad \frac{\Delta(H_A L_{\rm unit})}{g_{\rm grav}}
=.0100415734715.
$$

The exact rational16-state mass-mixing/spin-rotation control checks this
formula; the spatial spectrum and I4 use the actual shared parent.
The singlet and triplet contact energies in the stated finite-mode Wick
convention are-.00844104g and+.00160053g in ELunit units. Their difference
is not an arbitrarily assigned offset. Other internal/mixed torsion
currents, inherited occupation, vacuum counterterms and backreaction have
not been matched. Neither an isolated negative contact energy nor this
splitting is a total radial-pressure or NEC calculation.

This is the specific next input: the inherited covariance AND connected
four-point function on the real parent/receiver mode space, with the
complete stress tensor. A free Gaussian transport cannot select an
independent connected quartic moment. All the full-goal matter, eta,
vacuum/scale and regular-child gates remain open until their physical
connections are supplied and tested.

## 7. The complete bare ground-tower vertex

The [source record](../results/summaries/nsc-atlas-valve-source-gate-v1.json)
now gives every Killing-sign block of the bare minimal16D EC coupling,
not only the external axial component used by the example above. Use
{gamma4,gamma4}=-2eta, gamma5²=1 and Gamma16^i=i gamma5 tensor gamma12^i.
The general-dimensional three-form operator is borrowed from
[Castillo-Felisola and collaborators, equations13 and18](https://arxiv.org/abs/1405.0397);
its Clifford/adjoint convention is converted together. In our convention
L16_contact=-kappa16 J_ABC J^ABC/32, and the independent four-dimensional
identity J_munuomega²=6A² calibrates the external coefficient.

In ambient Spin(13) Killing trivialization, eta_plus(v,n)=v/sqrt(V_gamma)
and eta_minus(v,n)=Gamma(n)v/sqrt(V_gamma). Let Cq be the sum of products
of ambient gamma_[I1...Iq] over ALL ordered indices; C0=I tensor I.
For two bilinear Killing-sign pairs (a,b),(c,d), the internal projection is:

| Internal index count r | Both pairs diagonal | Both pairs off-diagonal |
| --- | --- | --- |
| 0 | C0 | C1/13 |
| 1 | -ac C2/13 | -12ac C1/13 |
| 2 | 11 C2/13 | C3/13 |
| 3 | -ac C4/13 | -10ac C3/13 |

A diagonal/off-diagonal mixed block vanishes by odd n parity. These are
ordered sums: converting C2 or C4 to increasing indices contributes2!
or4!, not an optional degeneracy factor. The full vertex is

$$
\mathcal L_4=-\frac{\kappa_4}{32}
\sum_{r=0}^3\binom3r(-1)^r
B_r^{av,bw}B_r^{cx,dz}K_r^{ab,cd}[vw;xz],
$$

where B_r=barpsi_a gamma_[mu1...mu(3-r)]gamma5^r psi_b.
The warp/field factors cancel: exp(-12s) from volume times exp(12s)
from four canonical spinors. Normalized Killing profiles give
kappa16/V_gamma=kappa4. No new contact parameter is introduced.
The positive-mass rotations psi_a=exp(-ia pi gamma5/4)chi_a must be
carried into EVERY left/right vertex, especially cross-sign terms.

Pointwise intrinsic gamma products and exact signed-axis sphere averages
check all128 sign/channel blocks in internal dimensions2 and4. A direct
six-dimensional Lorentz three-form contraction independently checks the
binomial, ordered-index and i factors; both mass rotations are checked.
The action projection is now explicit. It is still not a consistent
nonlinear finite-KK truncation: higher modes can be sourced. Occupation,
renormalized stress and backreaction have not been inferred from the
coefficient table. Isolated contact energy is not its complete radial NEC.

## 8. A conditional state response driven by the same parent

The [state-onset calculation](../results/summaries/nsc-atlas-valve-dirac-onset-v1.json)
uses the actual initial parent, its just-computed Dirac doublet, and the
same Einstein--scalar evolution. Its extra preparation is explicit:
initially empty selected low modes in the instantaneous initial quadratic
Hamiltonian basis. The UV completion and full interacting state are not
selected by this prescription. Fermionic statistics and unitary quantum
evolution are imported. No temperature, damping, incident waveform or
observational number is adjusted to produce the response.

At K=Pi=0, Hdot=0. Write q0=Psi²/H and rr0=u q0-(1-u²)q0,u. If
d=(ka-kb)/(1-u²), the actual initial metric acceleration gives

$$
Z_u=\dot d,\quad Z(0)=0,\qquad
\ddot\Omega=-\dot k_bq_0-rr_0 Z.
$$

This is obtained by matching radial and angular metrics after an
isotropic-coordinate change, not by replacing the inhomogeneous metric
with a homogeneous scale factor. Hellmann--Feynman and the geometric
radion coupling then give

$$
\ddot\lambda_\Sigma=-\lambda_\Sigma
\int\sqrt{1-u^2}\,\ddot\Omega(F^2+G^2)du,
\qquad
\overline{\ddot m}=-\frac\kappa2m\mu_\chi^2A^2\langle u^2\rangle_h.
$$

The physical mode has <u²>_h=.37046961549. The same initial parent gives
lambda_ddot=-.105656578937 and projected mass_ddot=-.849754399712.
The off-diagonal same-mode positive/negative-energy Hamiltonian Hessian is

$$
M_{pn}=\frac mE\ddot\lambda_\Sigma
-\frac{\lambda_\Sigma}E\overline{\ddot m}
=-.102215098653+.215113259398=.112898160745.
$$

The extrinsic trace connection cancels under canonical half-density
normalization. The radial-coordinate change has its own shift/unitary
connection; its same-lambda pair projection vanishes, as does the relevant
commutator diagonal. This control is needed before calling the result a
physical-slice Hamiltonian coefficient.

The fixed-initial-reference propagator component is

$$
(P_+^0 U(\tau)P_-^0)_{pn}
=-i\frac{\tau^3}{6}M_{pn}+O(\tau^4)
=-i\,.0188163601242\,\tau^3+O(\tau^4).
$$

An instantaneous FINAL particle projector is a different operator:
its leading basis-rotation term is tau² M_pn/(4E), coefficient.0105188209344.
The cubic fixed-reference amplitude is not sold as that particle count,
a baryon excess, a time-origin theorem or the total inherited source.
The isotropic doublet pairs in the spin-singlet channel, with zero charge.

Independent finite evolution at proper times.02,.01,.005,.0025 reconstructs
the actual spatial metric and its Dirac spectrum. The inferred M_pn errors
decrease2.92e-5 ->7.57e-6 ->1.91e-6 ->4.79e-7. The coordinate map uses
J(u)=integral0^u(c/b)dv, newv=(u+tanhJ)/(1+u tanhJ), and
Omega_new=b(coshJ+u sinhJ), avoiding endpoint singularities without
changing any field. Two interpolation/integration orders are reported.
Closed dS gives lambda_ddot=-H²lambda and mass_ddot=0; the four-dimensional
massless control removes this same-mode pump.

This computes a state-response component from a declared initial state
and the shared parent dynamics. It does not choose the full initial
quantum state, turn a probe into a complete semiclassical solution, or
close the child/source gate. The remaining calculation is the complete
renormalized radial stress from the transported state and vertex, including
kinetic stress, all required modes and interaction-generated correlations.

## 9. The physical uplift does not hide a turn at the admitted cut

The declared reconstruction is
ds16²=e^(-12s) g4+e^(2s) gamma12. Its full closed codimension-two section
is S2 times S12. The two warp factors cancel in its area:

$$
\mathcal A_{16}=(4\pi e^{-12s}R^2)(e^{12s}V_\gamma)
=4\pi V_\gamma R^2,\qquad
\theta_{16,\pm}=e^{6s}\theta_{4,\pm}.
$$

The fixed internal shape, complete fiber and declared block-diagonal
ansatz are retained. The null normals are the corresponding future
orthonormal external null directions, with no internal component.
Thus the full lifted section
remains future trapped for any finite smooth roll. This statement concerns
these actual compact sections; it does not assert conformal invariance
of every apparent-horizon construction or exclude a later transition.

Even the external physical radius Rhat=e^(-6s)R is still contracting at
the current cut. Its physical-time derivative is Rdot-6R sdot. Put
t=e^(-2s); the potential is
V=t6[7-(15/2)t+(25/32)t6]. The bracket has global minimum
7-(25/4)(8/5)^(1/5)=.133996604337>0. Positivity also follows exactly from
(28/25)^5>8/5. Hence the canonical density bounds
rho>=42 sdot²/kappa. On the retained action and source,

$$
\frac{d\widehat R}{d\widehat\tau}
\le R_\tau+R\sqrt{\frac{6\kappa\rho}{7}}.
$$

The [admitted restart](../results/summaries/nsc-atlas-valve-restart-v1.json)
has R=.149532673785, R_tau=-2.84970627770 and rho=.697153540142.
The right side is -2.73411450537. Equivalently, the contraction rate
19.05741538 exceeds the largest allowed roll correction .773020166144.
Changing to the physical uplift therefore does not supply an expanding
receiving sphere at this cut. No additional evolution or grid is used.

## 10. The action determines stress and radion exchange together

Write the complete matched bare fermion action as
$\mathcal L_D=\mathcal K-\sum_a m_a(S)B_a+\mathcal C_4$, with
$B_a=\bar\psi_a\psi_a$ and $\mathcal C_4$ the full quartic vertex of Section7.
Here $\mathcal K=(i/2)\sum_a\bar\psi_a\gamma^\mu
\overleftrightarrow\nabla_\mu\psi_a$.
Vary before using the nonlinear Dirac equation. In the stated Clifford
convention the kinetic tensor and full stress are

$$
K_{\mu\nu}=-\frac i4\sum_a\bar\psi_a
\left(\gamma_\mu\overleftrightarrow\nabla_\nu+
\gamma_\nu\overleftrightarrow\nabla_\mu\right)\psi_a,
\qquad T^D_{\mu\nu}=K_{\mu\nu}+g_{\mu\nu}\mathcal L_D.
$$

At fixed canonical spinors the algebraic contact varies through the volume:
$T^{\rm direct}_{\mu\nu}=g_{\mu\nu}\mathcal C_4$, hence
$T^{\rm direct}_{kk}=0$. A negative direct contact energy does not specify
negative radial NEC. Quartic homogeneity and the interacting Dirac equation give

$$
\mathcal K=\sum_a m_aB_a-2\mathcal C_4,\qquad
\mathcal L_{D,\rm on}=-\mathcal C_4,\qquad
T^D_{\mu\nu}=K_{\mu\nu}-g_{\mu\nu}\mathcal C_4.
$$

Thus $T^D_{kk}=K_{kk}$ contains the **interacting** kinetic response. Substituting
free probe derivatives into this on-shell expression omits that response.
The bare reduction has $\partial_S\mathcal C_4=0$: its warp factors cancel.
With $m_{a,S}=-7m_a/\sqrt{84}$ in the retained units, the bare local Ward
identity and radion equation are

$$
\nabla^\mu T^D_{\mu\nu}
=-\sum_a m_{a,S}B_a\,\partial_\nu S,
\qquad
\Box S-W_{,S}=\sum_a m_{a,S}B_a.
$$

The bosonic stress carries the opposite exchange. A renormalized source
must vary its counterterms too: define
J_S=sum m_a,S<B_a>_ren-(1/sqrt(-g))delta S_ct/delta S, including any
derivative terms. Then div<T_D+T_ct>_ren=-J_S dS and
Box S-W_,S=J_S. A subtraction of a number from the energy is not this
matched procedure. The free-field conserved point-splitting construction
in [Dappiaggi--Hack--Pinamonti, section4](https://arxiv.org/html/0904.0612v1)
is a method input, not a computed interacting source on this parent.
At first contact order,
the direct quartic expectation and the retarded correction to quadratic
stress must both be retained. In the current units free fermion backreaction
starts at order $g_{\rm grav}$; these two interaction contributions start
together at order $g_{\rm grav}^2$. Its physical value remains unselected.

An exact homogeneous sanity check is $n=N/a^3$,
$\rho=m(S)n-cn^2$, $p=-cn^2$ for constant $c$. Direct differentiation gives
$\dot\rho+3H(\rho+p)=m_{,S}n\dot S$, closing the opposite radion work.
This checks the exchange algebra; it supplies no spin-fluid replacement,
occupation rule or turnaround stress for the actual parent.

## 11. Gaussian coherence is a state input carried by inheritance

The [Gaussian-source record](../results/summaries/nsc-atlas-valve-gaussian-source-v1.json)
separates two questions. The balanced charge-conjugation, parity and chirality
benchmark has zero coincident covariance difference and its specified
fully vacuum-Wick-ordered contact cancels. Derivative stress can still be
nonzero. This cancellation does **not** classify every Gaussian state.

For a finite CAR covariance $0\leq K\leq I$, let $\Gamma_*$ denote the full
sixteen-dimensional chirality grading. Hermitian $A$ with $A^2=I$ and
$\{A,\Gamma_*\}=0$ gives the admissible Gaussian family

$$
K_\epsilon=\frac12(I+\epsilon A),\qquad |\epsilon|\leq1,\qquad
[K_\epsilon,\Gamma_*]=\epsilon A\Gamma_*.
$$

It can carry cross-chirality coherence and a nonzero Wick contact while
remaining Gaussian. An independent connected four-point function is therefore
one possible additional datum, not the only way to leave the cancelling
benchmark. Choosing the Gaussian coherence is itself a state choice.

For the free massless bulk field, the complete graded one-particle evolution
intertwines chirality. More generally an admitted spinor handoff satisfies
$\Gamma_*^{\rm out}U=\varepsilon U\Gamma_*^{\rm in}$, where
$\varepsilon=1$, or $-1$ for an orientation-reversing spin lift. Consequently

$$
K_{\rm out}=UK_{\rm in}U^\dagger,\qquad
\Xi_{\rm out}=\varepsilon U\Xi_{\rm in}U^\dagger,
\qquad \Xi=[K,\Gamma_*].
$$

Zero coherence is preserved by this full free transport; a nonzero coherence
is carried rather than selected. A chirality-preserving projection also
preserves zero. A projection that mixes gradings or an open channel can
change reduced data, but needs its actual state, projection and discarded
correlations. The operator inheritance law
$\mathcal T^*\mathbb D=\mathbb D$ does not choose those data or $K_{\rm in}$.
The known vertex and pump now constrain a source calculation; a conserved
renormalized radial source and its backreaction remain to be obtained.

## 12. The source operator couples different spatial modes

The [full initial operator record](../results/summaries/nsc-atlas-valve-dirac-source-v1.json)
removes the earlier same-mode closure. Retain the actual proper-time chart
and write the spatial metric as a² dchi²+R² dOmega2². A polar half-density
spinor has the ordinary dchi inner product. In each total-angular-momentum
block J=k+1/2, with multiplicity2k+2, its self-adjoint spatial operator is

$$
D_\Sigma=\epsilon\left(A\partial_\chi+\frac12\partial_\chi A\right)
-\frac{k+1}{R}\sigma_1,\qquad A=1/a,\qquad \epsilon=-i\sigma_2.
$$

The polar radial map is O=[[cos(chi/2),-sin(chi/2)],
[-sin(chi/2),-cos(chi/2)]], multiplied by sqrt(q)sin(chi) on the initial
slice q=Psi40²/H. It is a basis reflection. The existing right-frame
F-iG sigma.n convention and prior spectral/onset values are unchanged.
No time-dependent isotropic-coordinate shift is omitted, because this
calculation stays in the proper chart. With zero initial extrinsic curvature,

$$
\ddot A=\frac{\dot K_a}{q},\qquad
\frac{d^2}{d\tau^2}\left(\frac1R\right)=\frac{\dot K_b}{R},\qquad
\ddot m=-\frac{\kappa}{2}m\mu_\chi^2\epsilon_0^2 A_\chi^2 u^2.
$$

Here epsilon0 A_chi is the already declared scalar amplitude, not the
radial matrix epsilon. The parent ADM equations supply K_a,dot and K_b,dot.
Projection gives all entries of H_ddot=tau1 tensor D_ddot+tau3 tensor m_ddot,
including unequal modes. For initial spatial eigenvalue lambda_i set
E_i=sqrt(m²+lambda_i²), U_i=sqrt((E_i+m)/(2E_i)),
V_i=lambda_i/sqrt(2E_i(E_i+m)). The positive/negative energy-pair Hessian is

$$
\mathcal M_{ij}=(U_iU_j-V_iV_j)\ddot D_{ij}
-(U_iV_j+V_iU_j)\ddot m_{ij}.
$$

At k=0 the existing cutoffs8/16/24 give the finest diagonal
M_00=.112898160745309, agreeing with the old independent isotropic
Hellmann--Feynman calculation to1.2e-15. The squared ground-column norm is
.0448473832573;71.5791786% lies in other modes. This is a squared
initial-reference amplitude-coefficient fraction, not an energy fraction
or an observed particle count. Its final adjacent spread is3.9e-11.
The strongest other spatial eigenvalue is1.31422697688, with
|M_i0|=.154109235008. No observation selected these numbers.

The independent polar reconstruction recovers the full original round
Dirac projection to3.62e-12, the metric Gram matrix to1.78e-15, and its
Hermitian Hessian to8.88e-16. Round de Sitter controls recover the complete
analytic diagonal and no unequal-mode mixing. A hypothetical four-dimensional
m=0 control has zero same-mode coefficient yet nonzero unequal-mode
couplings on the inhomogeneous parent; it is not a physical massless
sector of the declared internal compactification.

This fixes a genuine missing part of the initial generator. A whole-field
occupied state still needs its full covariance/UV completion, stress
renormalization and actual transport. The finite matrix does not establish
a negative NEC, a regular child, a chosen g or a cosmological observable.
