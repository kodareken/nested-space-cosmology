# A geometric interaction behind the cosmohedron seed

This continues the [assessment](nsc-cosmohedra-assessment.md) and the
[sourced gravitational body](nsc-atlas-gravitational-response.md). The aim is
one interaction with frequencies and coupling calculated on the actual nested
geometry, followed by a test of the two-vertex wavefunction energy kernel.

## Preparation and chosen route

The original supporting wall scalar has a nonzero background; its fluctuations
mix with polar gravity. Treating them as independent oscillators would not be a
valid normal-mode calculation. An alternative is the full cubic Einstein action
for three gravitational modes. The bounded route chosen here is an explicit
probe extension: add one canonical surface scalar chi with zero background on
each existing wall. Its background and linear stress vanish, preserving the
sourced static solution. Its linear positive harmonic modes and its cubic
coupling to axial metric fluctuations can then be derived directly.

Inputs fixed before evaluating the interaction:

- The existing two-wall body, M=1, outer wall 2.5M, inner wall .25M,
  intermediate mass .1M and the same solved vacuum/source data.
- A finite conditional perturbation domain through r=10M, with regular center
  and Dirichlet axial master at the outer worldtube. This is a stated quantum
  preparation/domain, not an additional physical spacetime wall or the earlier
  driven-worldtube scattering state.
- The first ell=2,m=0 axial mode, normalized through its quadratic action.
- On each wall the real scalar harmonics ell=1,m=1 cosine and ell=2,m=1 sine.
  These have different frequencies, avoiding the equal-frequency current
  selection zero. The scalar zero mode is outside this retained sector.
- The positive-frequency Gaussian reference state uses the one shared static
  time. Quantum theory is imported explicitly; it is not derived from fluid
  averaging. The calculation is a projected tree-exchange contribution,
  not a nonperturbative vacuum or the complete Einstein four-point functional.

Executed finite route:

1. Derive the scalar-mode normalization, angular selection integral and cubic
   metric/current interaction from the action, including the Legendre transform.
2. Compute the physical axial mode using regular-cap data and the existing
   signed interface transfer; obtain a norm independently by the eigenvalue
   Wronskian identity and spatial quadrature.
3. Calculate both wall vertices. Test the actual momentum-type interaction in
   the perturbative Schrodinger equation and identify its four-point pole kernel.
4. Compare the same outer geometry with inner radius .30M, re-solving only the
   inner stiff junction. It must independently pass radial stability.
5. Test scale dependence with G fixed, and distinguish geometric homogeneity
   from inheritance of the dimensionless quantum coupling G/M^2.

Allow 180 CPU seconds and 400000 RHS/quadrature evaluations, with up to four
useful workers. Preserve compact mode values, vertices, convergence indicators
and source hashes, without mode histories. No old scalar or mass-closure grids
are reopened. Complete only after canonical normalization, all required cubic
conditions and a nonzero numerically resolved vertex pass independent checks.

Likely failures are an omitted harmonic/clock factor, a pure-gauge or
equal-frequency coupling, and confusing a propagating exchange contribution
with the full four-point coefficient. Normalization is checked against the
published gravitational action/flux convention; direct angular integration
checks the selection rule; constraint/contact terms are identified separately.
A later cubic-action or direct quantum-mode calculation disagreeing with the
derived vertex would falsify this finite result.

## Completed action and normalization

The extension is explicit:

$$
S_\chi=-\frac12\sum_s\int_{\Sigma_s}\sqrt{-\gamma}\,
\gamma^{ab}\partial_a\chi_s\partial_b\chi_s\,d^3y.
\tag{C1}
$$

Its zero background has zero stress and zero first stress variation. The
supporting theta_s and the Israel data are unchanged. Theta and chi may be
viewed as two real canonical components with only theta rolling. Vacuum
loops and semiclassical backreaction are outside this tree calculation.

With wall lapse N_s=A sqrt(f), real unit harmonics and the common time,

$$
\chi_s=\frac{\sqrt{N_s}}{R_s}\sum_i Q_{si}Y_i,\qquad
L_{2s}=\frac12\sum_i(\dot Q_{si}^2-\omega_{si}^2Q_{si}^2),\qquad
\omega_{s\ell}=\frac{N_s}{R_s}\sqrt{\ell(\ell+1)}.
\tag{C2}
$$

Only positive harmonics are retained; the shift zero mode is excluded.
For delta gamma_tA=h0 X_A, use X^A=epsilon^AB D_B Y20 and
Y20=sqrt(5/(16pi))(3cos(theta)^2-1). This X is minus Martel--Poisson's
vector harmonic. Use the correspondingly sign-reversed master coordinate;
the two vertices reverse together and their product stays the same.

Direct metric variation and the Legendre transform give

$$
L_{3s}=-\frac{h_0}{R_s^2}\sum_{ij}C_{ij}\dot Q_{si}Q_{sj},\qquad
C_{ij}=\int Y_i X^A D_A Y_j\,d\Omega=-C_{ji}.
\tag{C3}
$$

For Y1=sqrt(3/(4pi)) sin(theta) cos(phi) and
Y2=sqrt(15/(4pi)) sin(theta) cos(theta) sin(phi), C12=3/(2sqrt(pi)).
The polar integral is integral mu^2(1-mu^2)dmu=4/15 and the azimuthal
integral is pi. There is no freely supplied cubic coupling.

The actual axial operator is L=-d_x^2+U with dx/dr=1/(Af) and the existing
shell jump [e_x]=q_s e. Impose e(10M)=0 and regularity at the center.
This conditional mode domain differs from the earlier open characteristic
experiment; it does not add a physical spacetime wall or select a global
astrophysical vacuum. A trial Rayleigh quotient brackets the first mode.

The continuous phase, amplitude and eigenvalue sensitivity give its root
and norm. In M=1 units, e=exp(B)sin(theta), p=exp(B)cos(theta),
lambda=omega_g^2, and

$$
\theta_x=\cos^2\theta+(\lambda-U)\sin^2\theta,\qquad
B_x=(1+U-\lambda)\sin\theta\cos\theta,
\tag{C4}
$$

$$
\eta_x=2(\lambda-U-1)\sin\theta\cos\theta\,\eta+\sin^2\theta,
\qquad\eta=\partial_\lambda\theta.
\tag{C5}
$$

The existing real interface shear transports all three quantities. At the
first Dirichlet root theta_end=pi, the Wronskian identity gives
integral e^2 dx=p_end^2 eta_end. Independent spatial quadrature checks the
unit norm. The implementation restores k=1/M in the phase coordinates for
the scale control; it does not reuse the M=1 formulas without their units.

[Martel--Poisson](https://arxiv.org/html/gr-qc/0502028), the master relation
and flux normalization (HTML equations 66/84), fix

$$
S_g=\frac{\alpha_\ell}{2}\int dt\,dx
(\dot\Psi_C^2-\Psi_{C,x}^2-U\Psi_C^2)
-\frac{\alpha_\ell}{2}\sum_s\int dt\,q_s\Psi_C(x_s)^2,
\tag{C6}
$$

$$
\alpha_\ell=\frac{(\ell-1)\ell(\ell+1)(\ell+2)}{64\pi G},\quad
\Psi_C=\frac{e(r)Q_g}{\sqrt{\alpha_\ell}},\quad
h_0=\frac{Af\Psi_C+r\Psi_{C,x}}2.
\tag{C7}
$$

The real tensor-harmonic norm is 12 for ell=2, independently giving
alpha_2=3/(8piG). The interface contributes a potential, not a kinetic
norm. The retained Hamiltonian interaction is therefore

$$
H_3=Q_g\sum_s c_s(P_{s1}Q_{s2}-P_{s2}Q_{s1}),\quad
c_s=\frac{C_{12}}{2R_s^2\sqrt{\alpha_2}}
(Af e+R_s e_x)_{\Sigma_s}.
\tag{C8}
$$

Both sides of each junction give the same h0 and c_s. The vacuum pair
numerator is n_s=c_s(omega_s1-omega_s2), so degenerate frequencies give
zero pair creation. The selected unequal harmonics pass this check.

## Actual-mode numbers and the quantum graph

For the .25M child and G=1 used to tabulate shapes,

| Quantity in exterior M units | Inner wall | Outer wall |
| --- | --- | --- |
| First scalar frequency | 1.36230331533 | .252982212813 |
| Second scalar frequency | 2.35957855748 | .438178046004 |
| Current coupling c_s | .001441720226 | .028346822842 |
| Pair numerator n_s | -.001437791887 | -.005249713475 |

The actual mode has M omega_g=.392804111132, below the Rayleigh upper
bound .437281067846, with no interior node. Its normalization and shape
determine the vertex; the prior assessment's arbitrary g_L/g_R is removed.

The positive Gaussian reference uses this one global time. It is a
perturbative preparation, not a nonperturbative ground-state claim for a
cubic truncation. With Psi=exp(-S), the first term is

$$
S_{3s}=\frac{i n_s}{\omega_g+\omega_{s1}+\omega_{s2}}
Q_g Q_{s1}Q_{s2}.
\tag{C9}
$$

The mixed coefficient in log(Psi) for the two distinct wall pairs is

$$
\Gamma_{1234}^{(g_1)}=-n_{\rm in}n_{\rm out}K_4,\qquad
K_4=\frac1{E_tE_LE_R},
\tag{C10}
$$

with Et=sum omega_i, EL=omega1+omega2+omega_g and
ER=omega3+omega4+omega_g. Exact substitution in the Schrodinger equation
checks this current-vertex result, including its sign and imaginary cubic
phase. Its denominator matches [Cosmohedra, equation 49](https://arxiv.org/html/2412.19881v1);
the derivative numerator is calculated from our selected material action.

This is one retained propagating axial-mode tree contribution. Other radial
modes, even parity, nonpropagating constraint/contact terms and higher
orders are not summed. It is not the full Einstein chi^4 functional or a
complete positive-geometry construction for gravity. Mode energies replace
plane-wave momentum magnitudes; no literal polygon/spacetime identification
is asserted. The selected pole and its basis are specified, so the isolated
value is not advertised as the full measurable four-point function.

Within this projected Hamiltonian, its leading normalized connected moment
is C_selected=Gamma/(2sqrt(omega1 omega2 omega3 omega4)). The same scope
restrictions apply to that channel contribution.

## Interior contrast and inheritance

The .30M control keeps the exterior, outer wall and parent parameters fixed.
The inner junction is re-solved under the same material law and has positive
radial potential curvature. These are separately prepared static states;
no instantaneous change or causal signal is claimed.

| Quantity | Child radius .25M | Child radius .30M |
| --- | --- | --- |
| M omega_g | .392804111132 | .392804157226 |
| Gamma_selected M^4/G | -3.83479234069e-7 | -6.32140173389e-7 |
| C_selected M^2/G | -3.21208794678e-7 | -4.91767335107e-7 |

The selected coefficient changes 64.8434%; this one linear frequency changes
1.17345e-7 relatively. Inner scalar frequencies also change; no universal
degeneracy of all linear observations is claimed. The contrast is resolved
above the declared numerical comparison indicators.

At fixed dimensionless geometry, normalized e scales as M^(-1/2) and
frequencies as M^(-1). The actual action gives

$$
c_s\sim\sqrt G\,M^{-5/2},\quad n_s\sim\sqrt G\,M^{-7/2},\quad
\Gamma_{\rm selected}\sim GM^{-4},\quad C_{\rm selected}\sim G/M^2.
\tag{C11}
$$

The M=2 control preserves Gamma M^4/G and C M^2/G. Thus geometric shape
inherits, while quantum strength epsilon_q=G/M^2 runs at fixed G. Restored
units make this the squared Planck length over squared geometric mass length.
A physical T*D=D map must carry this parameter and state; identical quantum
interactions do not follow from a geometric rescaling alone. The tree/probe
regime uses epsilon_q small. Tabulating shapes with G=1 does not assert a
controlled Planck-scale cubic model.

## Evidence and continuation

The [compact record](../results/summaries/nsc-atlas-cosmohedron-interaction-v1.json)
binds the action, conditional state/domain, modes, vertices, three solver
tolerances and M=2 control. It uses 100597 RHS calls and about 1.6 mode-solver
CPU seconds in four workers. Independent cap/spatial normalization and
direct Regge--Wheeler propagation check the phase method. Eight new controls
pass. An [additional Fock-space check](../tests/test_nsc_cosmohedron_fock_exchange.py)
independently recovers the mixed moment and kernel. In that calculation,
the vacuum-to-second-order path and the two first-order states combine to
cancel the apparent extra internal-frequency denominators. Grok independently
reviewed the action and mode normalization; the coefficient algebra is checked
by the two explicit quantum derivations above.

Next: include the remaining constrained axial terms and/or a controlled pole
sum before claiming a full cross-wall four-point prediction, then extend to
a time-dependent parent. This finite vertex/kernel result is complete at its
stated scope. No black/white transition or Lambda value follows from the
polygon terminology.

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_cosmohedron_interaction.py
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_cosmohedron_interaction.py tests/test_nsc_cosmohedron_fock_exchange.py
```
