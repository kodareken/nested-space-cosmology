# Charged monopole shell for the later RN experiment

This note prepares the charged experiment that follows the classical
Reissner–Nordström run. It does not evolve that run, does not change the
neutral producer, and does not accept the historical charged-neck or
\(V_{\rm full}\) numbers as its result. The classical parent remains the
large-\(r\) exterior and the child a distinct collar on its own areal
radius near the outer horizon. Both use the one ingoing flat-spatial chart

\[
ds^2=N^2dt^2-(dr+\beta\,dt)^2-r^2d\Omega^2,
\]

with \(\beta>0\), so the unit normal to a \(t=\)constant slice,

\[
u=N^{-1}(\partial_t-\beta\partial_r),
\]

moves toward decreasing \(r\). The outer boundary stays nonperiodic and
absorbing. The excision between \(r_-\) and \(r_+\) stays pure outflow: this
shell adds no boundary source there. \(N\) and \(\beta\) are the coupled
fields supplied by that chart. They are not replaced by \(N=1\).

## Common action, units, and the declared observer

Signature \(+---\), with \(c=\hbar=1\). The active classical branch uses the
native bulk term of the [RN reference](nsc-rn-reference.md),

\[
S=\int\sqrt{|g|}\,\bigl[-A R_L-C_F F_{\mu\nu}F^{\mu\nu}\bigr]+S_{\rm Dirac}.
\]

The laboratory Riemann convention of that owner is

\[
R^a{}_{bcd}=\partial_c\Gamma^a{}_{db}-\partial_d\Gamma^a{}_{cb}
+\Gamma^a{}_{ce}\Gamma^e{}_{db}-\Gamma^a{}_{de}\Gamma^e{}_{cb},
\qquad R_{bd}=R^a{}_{bad}.
\]

\(R_L\) is the scalar built from this tensor. Older matching notes write
\(A R\) with the opposite scalar \(R_E=-R_L\) recorded in the
[compact matching](nsc-compact-matching.md), so \(A R_E=-A R_L\). Those
labels are the same bulk term only after that sign. This note keeps
\(-A R_L\).

\(A>0\) has dimension length\(^{-2}\), \(C_F>0\) is dimensionless, and
\(A_\mu\) has dimension length\(^{-1}\). The reference identities are

\[
G_N=\frac1{16\pi A},\qquad g_4^2=\frac1{4C_F}.
\]

On this branch the cosmological coefficient is \(V_4=0\). A general declared
\(V\) belongs only to a separate future branch. It is not an input of the
asymptotically flat prototype. The Maxwell term and the Dirac term each
contribute once. The minimal coupling sits inside \(S_{\rm Dirac}\) and is
not a second force. No four-dimensional Weyl counterterm is added. The
neutral prototype stays the pending production path.

The volume density of the chart is \(\sqrt{|g|}=N r^2\sin\theta\). The
declared frequency observer is the exterior Killing observer of
`nsc_rn_observables`: modes \(e^{-i\omega t}\) in the gauge \(A_t(\infty)=0\),
accepted when the exterior Killing frequency satisfies \(\omega>0\). The
Eulerian normal \(u\) measures the slice energy and the electric field. It
does not replace that frequency test.

One four-dimensional massless Dirac field is the charged matter. The paired
compact copies already supply that one field; they are not counted again.
Its charge sign \(s=\pm1\) enters only through

\[
D_\mu=\nabla_\mu-isA_\mu.
\]

Both signs are later runs. This wave selects neither.

## Monopole harmonics and the lowest-Landau multiplicity

The fixed magnetic background is the integer monopole bundle already stated
for this geometry,

\[
F=\frac{q_{\rm mag}}{2}\sin\theta\,d\theta\wedge d\phi,
\qquad
\frac1{2\pi}\int_{S^2}F=q_{\rm mag}\in\mathbb Z.
\]

Flux is inherited and not refit. On the unit sphere the orbital generators
of a unit-charge section obey \(\mathbf L\cdot\hat r=-q_{\rm mag}/2\). The
total angular momentum is \(\mathbf J=\mathbf L+\boldsymbol\sigma/2\), so

\[
\mathbf J\cdot\hat r=\frac{\sigma_r-q_{\rm mag}}{2}.
\]

For \(q_{\rm mag}>0\) the eigenvalue \(\sigma_r=+1\) gives
\(|\mathbf J\cdot\hat r|=(q_{\rm mag}-1)/2\), and \(\sigma_r=-1\) gives the
larger value \((q_{\rm mag}+1)/2\). The lowest multiplet therefore has a
single chirality and

\[
j=\frac{|q_{\rm mag}|-1}{2},\qquad
2j+1=|q_{\rm mag}|,\qquad
m=-j,\ldots,j.
\]

The Dirac angular eigenvalue on that multiplet is

\[
\kappa=\pm\sqrt{\Bigl(j+\tfrac12\Bigr)^2-\Bigl(\frac{q_{\rm mag}}{2}\Bigr)^2}=0.
\]

The radial coefficient \(\kappa/r\) vanishes with \(\kappa\). It is an angular
eigenvalue divided by the areal radius, not a four-dimensional rest mass.
The field stays massless. The opposite bundle orientation,
\(q_{\rm mag}<0\), exchanges the chirality and leaves the multiplicity
\(|q_{\rm mag}|\).

The higher multiplets are \(j=(|q_{\rm mag}|-1)/2+n\) with \(n\ge1\),

\[
\kappa=\pm\sqrt{n(n+|q_{\rm mag}|)},
\]

with \(|q_{\rm mag}|+2n\) states on each sign and \(2(|q_{\rm mag}|+2n)\)
states together. The coordinate section of the lowest multiplet, for
\(q_{\rm mag}>0\) and \(s=+1\), is the holomorphic monopole harmonic

\[
\eta_-\propto e^{im\phi}
\Bigl(\sin\frac\theta2\Bigr)^{j-m}
\Bigl(\cos\frac\theta2\Bigr)^{j+m}.
\]

Each \(m\) channel is normalized in the round measure. The canonical radial
unknown of this shell has shape \((2,\text{points},|q_{\rm mag}|)\): the
first axis is the \(\sigma_2\)-diagonal characteristic basis, and the rank is
exactly this multiplicity. There is no extra weight of \(4\) and no second
copy.

The neutral closed shell is a different eigenspace. Its labels are
\(\kappa=\pm1,\pm2,\ldots\), and the shell \(\kappa=1\) has multiplicity
\(4\times1=4\), entered once. That integer equals \(|q_{\rm mag}|\) only when
\(|q_{\rm mag}|=4\). The operators still differ: the neutral shell has
\(|\kappa|=1\), while the lowest Landau shell has \(\kappa=0\).

The canonical radial column is the rescaled orthonormal spinor
\(\phi=r\psi\), the same radial identification the neutral shell uses.
A normalized packet in the weighted radial Gram,

\[
\sum_a w_a\,\phi^\dagger\phi=1,
\]

with \(w_a\) the radial summation-by-parts weights, has canonical charge
\(s\) in that Gram. The angular coefficient \(\kappa/r\) is identically zero,
and the charge row is not. Vanishing \(\kappa\) is therefore an empty angular
barrier with a nonzero Gauss source. The later occupation of each channel is
the actual Gaussian weight in this same Gram. Electric current of the
classical neutral run remains absent; the current below belongs only to this
later shell.

## Electric field, Hamiltonian shift, and Gauss law

The Eulerian coframe has \(e^{\hat0}=N\,dt\) and
\(e^{\hat1}=dr+\beta\,dt\). Its dual satisfies
\(E_{\hat0}{}^t=1/N\) and \(E_{\hat0}{}^r=-\beta/N\), while
\(E_{\hat1}=\partial_r\). Therefore the coordinate gamma matrix on this chart is

\[
\gamma^t=\frac1N\gamma^{\hat0}.
\]

The shift sits in \(E_{\hat0}{}^r\), not in \(\gamma^t\). For the declared
coupling \(D_\mu=\nabla_\mu-isA_\mu\) and \(\bar\psi=\psi^\dagger\gamma^{\hat0}\),
with \((\gamma^{\hat0})^2=1\),

\[
\bar\psi\gamma^t\psi=\frac{\psi^\dagger\psi}{N}.
\]

The Noether current that enters \(S_{\rm Dirac}\) is
\(J^\mu=s\bar\psi\gamma^\mu\psi\), so

\[
J^t=\frac{s\,\psi^\dagger\psi}{N}.
\]

Here \(\psi\) is the Eulerian orthonormal spinor. The canonical column is
\(\phi=r\psi\). The lapse in the volume density and the lapse in \(J^t\)
cancel in the scalar density,

\[
N r^2 J^t=s\,\phi^\dagger\phi,
\]

before any Gauss integral is imposed.

The Maxwell components are fixed from the same chart, with
\(A=A_t(r)\,dt\). The inverse of the \((t,r)\) block gives

\[
F_{rt}=\partial_r A_t,\qquad
F^{rt}=-\frac{\partial_r A_t}{N^2},\qquad
F_{tr}=-F_{rt},\qquad
F^{tr}=-F^{rt},
\]

and the contraction

\[
F_{\mu\nu}F^{\mu\nu}=-2\mathcal E^2,
\qquad
\mathcal E=\frac1N\partial_r A_t.
\]

Variation of \(-C_F F_{\mu\nu}F^{\mu\nu}\) defines the Gauss charge by that
field strength,

\[
Q_E=16\pi C_F\, r^2\mathcal E.
\]

The Legendre transform of the same coupling on \(\phi\) puts \(-s A_t\) in
the Hamiltonian. These three relations are the sign convention. Only then
is the conservation law stated:

\[
\partial_r Q_E=\int_{S^2} N r^2 J^t\,d\Omega
=\int_{S^2} s\,\phi^\dagger\phi\,d\Omega.
\]

Positive canonical charge makes \(Q_E\) increase outward, so \(\mathcal E>0\)
and \(\partial_r A_t>0\). The gauge \(A_t(\infty)=0\) then forces \(A_t<0\).
The Coulomb tail at unit lapse is

\[
A_t=-\frac{Q_E}{16\pi C_F r},\qquad
\Phi:=-A_t=\frac{Q_E}{16\pi C_F r},
\]

and the Hamiltonian shift is \(-s A_t=s\Phi\). Like charges repel. On a
derived lapse the same sign is the inward integral

\[
A_t(r)=-\int_r^{\infty}\frac{N(r')Q_E(r')}{16\pi C_F r'^2}\,dr'.
\]

The radial operator sees \(\omega-s\Phi=\omega+s A_t\). At infinity this is
the exterior Killing frequency. Acceptance remains \(\omega>0\) in the gauge
\(A_t(\infty)=0\).

The orthonormal electric energy density is \(2C_F\mathcal E^2\). Outside a
sphere of radius \(R\) at unit lapse,

\[
U(R)=\frac{Q_E^2}{32\pi C_F R}.
\]

At \(C_F=1/4\) this is \(Q_E^2/(8\pi R)\), and \(Q_E\Phi(R)=2U(R)\). The
Einstein equation uses the Maxwell stress of \(-C_F F^2\) once.

The canonical integral of \(s\phi^\dagger\phi\) is the Gauss charge \(Q_E\).
The [RN reference](nsc-rn-reference.md) records the executable magnetic piece
of the fixed-charge branch,

\[
P_{\rm mag}^2=\frac{C_F q_{\rm mag}^2}{4A},
\]

where \(q_{\rm mag}\) is that executable flux and not the prose label
\(2\pi q_{\rm mag}\). The same Maxwell term converts an electric Gauss charge
by

\[
P_{\rm electric}^2=\frac{Q_E^2}{256\pi^2 C_F A}.
\]

The charge square in the Reissner–Nordström function is the sum

\[
P_{\rm tot}^2=P_{\rm mag}^2+P_{\rm electric}^2.
\]

With mass parameter \(M\),

\[
f=1-\frac{2M}{r}+\frac{P_{\rm tot}^2}{r^2}.
\]

Extremality is the ratio \(M/\sqrt{P_{\rm tot}^2}=1\). A subextremal hole has
ratio greater than 1, outer root

\[
r_+=M+\sqrt{M^2-P_{\rm tot}^2},
\]

and an open interval between the roots. The total extremal radius of that
solution is \(r_{m,\rm tot}=\sqrt{P_{\rm tot}^2}\). It is not the magnetic
radius \(\sqrt{P_{\rm mag}^2}\). The Misner–Sharp correction uses
\(P_{\rm tot}^2/(2r)\) and is not \(Q_E\).

On the vacuum chart \(N=1\), so the background potential derived from this
Gauss charge is \(\Phi(r)=Q_E/(16\pi C_F r)\), and

\[
\Phi_H=\frac{Q_E}{16\pi C_F r_+}.
\]

\(\Phi_H\) is not an independent datum. \(V_4\) remains 0. No negative vacuum
source is added to open a frequency window.

Parent and child keep separate areal nodes. For a field that satisfies the
null energy condition, the horizon inequality of Natário, Queimada and
Vicente uses this potential through \(\omega\ge s\Phi_H\). A classical Dirac
test field can leave that hypothesis inside \(\omega-s\Phi_H<0\), as stated
by Tóth. That comparison is initial data for a packet. It is not a new force
and not a cosmic-censorship claim. The horizonless MMP return stays in the
[MMP source study](nsc-rn-mmp-source-study.md).

## Pure-magnetic control and the dyonic window

Experiment 2 keeps the same action and fixes the executable magnetic flux at
\(q_{\rm mag}=1\). The default is one subextremal initial slice, not a pair of
charging and infall episodes. A later charging episode is admissible only
when its output state is preserved as the next initial slice. The smaller
default is the dyonic slice itself.

The pure-magnetic arm sets the black-hole electric Gauss charge to zero.
Then \(P_{\rm electric}=0\), \(P_{\rm tot}^2=P_{\rm mag}^2\) and \(\Phi_H=0\).
Every accepted exterior frequency \(\omega>0\) has

\[
\omega-s\Phi_H=\omega>0.
\]

That arm cannot enter the electric window \(\omega-s\Phi_H<0\). The charged
packet on it is the capture control.

The dyonic arm declares a nonzero black-hole electric Gauss charge of either
sign, derives \(P_{\rm electric}\), \(P_{\rm tot}\), \(r_+\) and \(\Phi_H\) as
above, and requires \(M/\sqrt{P_{\rm tot}^2}>1\). Both fermion signs are run.
For the sign with \(s\Phi_H>0\) the positive frequencies split into

\[
B_<=\{\omega:0<\omega<s\Phi_H\},
\qquad
B_>=\{\omega:\omega>s\Phi_H\}.
\]

Only \(B_<\) has \(\omega-s\Phi_H<0\). For the opposite sign, \(s\Phi_H<0\), so
every \(\omega>0\) stays in \(\omega-s\Phi_H>0\) on the same dyonic geometry.

The closed comparison uses one packet family at the same \(M\) and the same
\(q_{\rm mag}=1\). It distinguishes the electric window by where a response
sits:

- a response on the pure-magnetic arm is outside \(\omega-s\Phi_H<0\);
- a response in \(B_>\), or at the opposite fermion sign, is also outside
  that window;
- a response that is present in \(B_<\) and absent from those three controls
  is the electric window’s contribution.

No horizon flux is assigned a value here. The test reports the band, not a
censorship verdict.

Low bands are not the neutral frequency \(\omega r_{m,\rm tot}=1\). Whenever
\(s\Phi_H r_{m,\rm tot}<1\), that neutral frequency lies in \(B_>\). The
budget is finite. The parent cap for a solution is \(32\,r_{m,\rm tot}\) of
that solution, not \(32\) magnetic radii. A packet centred at \(\omega\) is
given the budget width \(\sigma_\omega=1/\omega\) and must satisfy

\[
4\sigma_\omega\le 32\,r_{m,\rm tot}-r_+.
\]

If \(B_<\) has no frequency that meets this inequality, the band is outside
the budget. The parent is not enlarged in that case, and \(\omega\) is not
sent to zero.

## Minimal variables the charged experiment adds

The classical state already carries \((N,\beta)\), the neutral \(\kappa=1\)
shell of multiplicity 4, the fixed flux, and the radial weights. The later
experiment adds only

| Variable | Role |
|---|---|
| \(s=\pm1\) | Both signs. The electric window opens only for \(s\Phi_H>0\). |
| \(q_{\rm mag}=1\) | Executable magnetic flux fixed on both arms. Rank of the lowest shell is 1. |
| \(Q_E^{\rm BH}\) | Declared black-hole Gauss charge. Zero on the pure-magnetic control; nonzero on the dyonic arm. |
| \(P_{\rm tot}^2\) | \(P_{\rm mag}^2+P_{\rm electric}^2\). Extremality ratio \(M/\sqrt{P_{\rm tot}^2}>1\). |
| \(\Phi_H\) | Derived, \(Q_E^{\rm BH}/(16\pi C_F r_+)\). Bands \(B_<\) and \(B_>\) use it. |
| \(A_t(r)\) | Rebuilt from the consistent \(Q_E\), \(N\) and \(C_F\), with \(A_t(\infty)=0\). |
| \(J^t\) | \(s\psi^\dagger\psi/N\). Then \(\int N r^2 J^t\,d\Omega\) is the canonical charge. |
| \(\omega>0\) | Exterior Killing frequency. The window test is the sign of \(\omega-s\Phi_H\). |

\(A\), \(C_F\), \(N\) and \(\beta\) are read. \(V_4\) stays \(0\) and is not a
state variable. No production trajectory and no scientific record are part
of this note.

## Manufactured check

Source-free algebra under root `scripts/lab.py`, with every BLAS and OpenMP
thread count set to 1. At \(N=1.7\) and \(\beta=0.35\), with
\(\partial_r A_t=0.3\), the contraction gives \(F_{rt}=0.3\),
\(F^{rt}=-0.10380622837370243\), \(\mathcal E=0.17647058823529413\) and
\(F_{\mu\nu}F^{\mu\nu}=-2\mathcal E^2\). The positive-charge Coulomb sign is
\(A_t<0\) and the Hamiltonian coefficient \(-s A_t\) equals \(s\Phi\).

The lapse cancellation uses five radii
\((1.0,1.4,1.8,2.3,2.9)\) and lapses
\((0.55,1.35,2.20,0.80,1.65)\), none equal to 1, with positive trapezoid
weights, \(\phi=r\psi\) and \(J^t=s\psi^\dagger\psi/N\). The weighted canonical
charge and the weighted Gauss integrand are both
\(2.6752110396318973\). Replacing \(J^t\) by \(s\psi^\dagger\psi\) leaves a
maximum pointwise gap \(3.888\) and a ratio equal to the lapse. With
\(C_F=0.5\), \(A=0.125\) and the sample \(Q_E=17.74075851438942\), the
geometric conversion gives \(q_{\rm geom}=1.4117647058823528\) and
\(r_Q^2=1.9930795847750862\), distinct from the Misner–Sharp length
\(0.49826989619377154\). That sample is pure electric, so its square is
\(P_{\rm electric}^2\) and not a dyonic total.

The window algebra uses \(A=0.125\), \(C_F=0.5\), \(q_{\rm mag}=1\),
\(Q_E=4\) and mass ratio 2. These are manufactured samples, not selected
initial data and not a measured flux. They give \(P_{\rm mag}^2=1\),
\(P_{\rm electric}^2=0.10132118364233778\), \(P_{\rm tot}^2=1.1013211836423378\),
\(r_{m,\rm mag}=1\), \(r_{m,\rm tot}=1.0494385087475768\) and
\(\Phi_H=0.04063643378570953\). The same mass with the electric charge removed
has \(\Phi_H=0\). The combination \(\omega-s\Phi_H\) at
\(\omega=\Phi_H/2\) is \(-0.020318216892854766\) on the dyonic arm and
\(+0.020318216892854766\) on the pure-magnetic arm. The neutral frequency
\(1/r_{m,\rm tot}=0.9528905139886873\) lies above \(\Phi_H\). The budget width
\(1/\omega=49.21691727543603\) exceeds the room
\(32 r_{m,\rm tot}-r_+=29.665474445857186\), so this low band is outside the
cap. The maximum absolute residual of the run is \(6.66\times10^{-16}\).
User time is 0.04 s. No production trajectory was run.
