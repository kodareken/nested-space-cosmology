# Parent-matched Dirac covariance and the child source budget

The standard fermionic characteristic-state prescription, combined with
the **complex reflection amplitude of the actual exterior Dirac operator**,
gives a definite massless source candidate on the expanding child.
In the same canonical heat convention as the previous reference, the
computed neck value is approximately

\[
\rho+p_\parallel=-0.05280.
\]

Both radial null contractions are negative after including the flux.
No scalar source, gravitational coefficient or target stress was fitted.
The full source budget does **not** close: the density has the opposite
sign to the positive-Einstein neck requirement, the pressures are unequal,
and the state carries outward parent power. The remaining common-action
contributions and backreaction must satisfy those equations together.

## What is imported, and what is calculated

Use the existing [DHP applicability and regulator audit](nsc-massless-reference.md).
DHP supplies the renormalized Dirac stress construction and its local
ambiguities. [Gérard–Häfner–Wrochna, Theorem 1.1](https://arxiv.org/html/2008.10995)
supplies a rigorous massless fermionic Unruh construction on its stated
Kerr–Kruskal domain: affine-frequency horizon data and incoming-vacuum
data at past null infinity. Neither theorem automatically proves the
characteristic trace and Hadamard hypotheses on the NSC black-universe
extension. Here that prescription defines a **boundary-state candidate**.
Global extension, spin/current trace maps and microlocal regularity remain
explicit conditions of its physical use.

Reuse the recorded [geometry and horizon](nsc-horizon-source.md), including
surface gravity 0.23832579963401956. Hawking temperature and the general
state construction are prior results. New work evaluates the actual
complex Dirac reflection, binds its phase to the real interior spin frame,
and transports the resulting occupied covariance through the full angular
source calculation. The expensive auxiliary reference is read from its
authenticated record; its generator is not run again.

The field here is one untwisted massless four-dimensional Dirac field.
The spatial domain is the unwrapped longitudinal line times S2. The
angular channels have |kappa|>=1 and multiplicity 4*kappa, with the
two-component Dirac trace inside each channel. The angular reduction's
effective mass is not a four-dimensional rest mass or a new coupling Phi.

## The phase-resolved state interface

In the exterior current basis,

\[
H_{\rm ext}=-i\sigma_3\partial_x+V_\kappa\sigma_1,
\quad V_\kappa=\kappa\sqrt A/r,\quad dx/d\rho=1/A.
\]

The up mode is outgoing plus R_up times incoming at the horizon and
outgoing at infinity. Its conserved current gives |R|^2+|T|^2=1.
Integrate the first-order ratio equation backwards from its outgoing
large-radius series. A uniform Frobenius frame removes the finite-collar
mixing before extracting R; a plane-wave extraction at a finite collar
would introduce a spurious phase error.

For positive omega, set f=(1+exp(2*pi*omega/kappa_h))^-1 and s=sqrt(f(1-f)).
The occupied affine-horizon covariance in a real spin frame is

\[
C_H=\begin{pmatrix}f&-is\\is&1-f\end{pmatrix},\qquad C_H^2=C_H.
\]

Compress the horizon inputs and empty incoming parent input with
M=((0,1,0),(R,0,T)). Since MM^dagger=I, this preserves the CAR bound.
For the interior k_z=-omega and H=kappa*sigma1+(k_z/b)*sigma3,

\[
\boxed{C_{\rm in}=\begin{pmatrix}
1-f&i sR^*\\-i sR&f|R|^2
\end{pmatrix}.}
\]

Its eigenvalues are 0 and 1-f|T|^2. The charge-related momentum channel
has \(C(+\omega)=I-\sigma_3 C(-\omega)^*\sigma_3\). These two channels give equal
diagonal stress and a nonzero momentum flux. The off-diagonal entries
are state correlations, not additional terms in the operator or a mass.

The factor i is necessary in this real frame. At zero frequency the known
decaying Rindler solution has R=-i; the interior covariance then becomes
(I-sigma1)/2, the negative-energy massive two-dimensional projector.
Omitting i would give the wrong energy projector. A development version
with that omitted phase was corrected before any accepted result.
Nonzero-frequency Rindler reflection is checked against the known
Gamma-function solution; it is not a new scattering theorem.

Under x -> x+c, R -> exp(2*i*omega*c)R. Opposite phases in the two
interior frame columns cancel this change. Independent covariance
transport checks the phase invariance and the collar limit. The optional
thermal incoming control adds f|T|^2 to C22 and cancels net power. It is
not identified as a proved global Hartle–Hawking state on this extension.

## Absolute source in the inherited canonical convention

Use q=arctan(-1/rho) on the child, continued to its horizon branch,
r=csc(q), b=sqrt(W), and

\[
\bar g=g/r^2=dq^2/W-Wdz^2-d\Omega^2,\qquad
W=3(\pi-q)+3\sin q\cos q-\sin^2q.
\]

Forward proper time corresponds to decreasing q. Unitary evolution from
the horizon transports the above covariance; no instantaneous vacuum is
reset later. For its rank-one occupied vector v, let beta and alpha be
the instantaneous positive- and negative-energy amplitudes and lambda
its norm. The mode energy above the instantaneous negative projector is
2*omega_inst*|beta|^2+omega_inst*(1-lambda). The second term retains the
missing occupation associated with outward power; dropping it would
violate the normalization of the energy source.

The unchanged [angular-stress owner](nsc-angular-stress.md) supplies the
orders 0,2,4 subtraction, joint angular finite part, static-cylinder
normalization and full conformal anomaly. Units are hbar=c=L_throat=1
and canonical normalization mu=1. No complete physical vacuum, Einstein
or curvature-squared coefficient is set to zero by this allocation.

The convergent angular remainder has a leading kappa^-3 term from the
sixth adiabatic order. Its integrated coefficient is compared with the
resolved channels; coefficient*zeta(3,N+1) estimates the omitted tail.
This is a numerical tail approximation, not a physical counterterm or
a rigorous bound. N=16 and N=32 corrected tensors agree within 4e-6
over the stated five radii. Separate frequency, quadrature and time
checks are recorded. Above omega=4 only exponentially small reflection
coherences are omitted numerically, with an explicit finite-N bound;
the theoretical covariance retains them at every frequency.

At the neck the N=32 corrected candidate is approximately

| Component | Value in the declared units |
|---|---:|
| rho | -0.0033069 |
| p_parallel | -0.0494946 |
| p_sphere | -0.042996 |
| rho+p_parallel | -0.0528015 |
| T_hatT_hatz | -0.0000030486 |
| outward parent Killing power | 0.00014222068 |

For null vectors (1,+/-1,0,0), the contractions are
rho+p_parallel +/- 2*T_hatT_hatz. Thus rho+p_parallel alone is the
average of the two contractions when flux is present. Parent Killing
power is not yet a cosmological density transfer Q: deposition, proper
volume, clock matching and a solved geometry are still needed.

The auxiliary reference gave rho+p_parallel=1.2128458668. Its difference
from this candidate is a state correction of about -1.26564736 in the
same convention. The state was specified by horizon/parent data before
this sign was calculated. Neither source is asserted to be the complete
finite-cutoff quantum stress of NSC.

## Independent source equations and next decision

At the imposed neck, the recorded geometry's two-derivative equations
with a_EH=A_EH*L_throat^2>0 require

\[
(\rho,p_\parallel,p_\perp,T_{\hat T\hat z})_{\rm required}
=\big(2a_{\rm EH},\ 2a_{\rm EH}(1-3\pi),\
2a_{\rm EH}(1-3\pi),\ 0\big).
\]

Their null sum is 4*a_EH*(1-3*pi/2), as already recorded. Matching only
that equation by choosing a_EH would conceal the density, anisotropy and
flux residuals. Leave a_EH symbolic. Remaining local/nonlocal terms
must supply the componentwise difference of this requirement and the
computed source, or the geometry must evolve under the same action.

The next source owner is the **same-state finite-cutoff causal conversion
and complementary functional**, together with the domain applicability
conditions above. Reuse the established conversion
Gamma5=Gamma_light,ren+[Gamma_H,nu+C_nu,mu] from the
[normalization map](nsc-massless-reference.md). Adding this canonical
stress to an unchanged cutoff complement would double-count terms.
Massive compact modes, gauge/boundary interactions and recursive state
correlations retain their own derived contributions. None may be
assigned a preferred value merely to cancel the residual.

## Reproduction

The immutable [result](../results/development/unruh-state.json) records
source and input hashes, exact interface checks, all stored source and
control rows, tail coefficients and the open source budget.

```sh
python -B scripts/check_nsc_unruh_state.py --check
```

This runs only the new parent-source calculation and checks every
published field with the recorded tolerances. It authenticates and reuses
earlier evidence. Raw reflection grids are working integration data;
they are not silently treated as published all-field comparisons.
