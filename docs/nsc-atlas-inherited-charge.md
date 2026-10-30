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
four-dimensional mode. Rotation in the (n1,n2) plane gives charges +/-1.
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
