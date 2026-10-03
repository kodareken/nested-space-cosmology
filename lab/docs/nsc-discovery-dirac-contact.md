# Dirac first-order contact, mean and connected

The curvature-EFT owner supplies one off-shell bulk contact and one
inverse-metric map. This note applies that contact to a Dirac-shaped
stress and separates what a finite Gaussian state actually determines.
The scalar/collective specialization in `spherical_contact` is not used
as a Dirac stress. No physical value of \(c_W\), \(c_R\) or \(A\) is
chosen. The renormalized four-fermion expectation stays open, and that
opening does not block the parent canonical evolution.

Owners are
[the calculation](../src/recursive_horizons/nsc_discovery_dirac_contact.py),
[CLI](../scripts/derive_nsc_discovery_dirac_contact.py) and
[tests](../tests/test_nsc_discovery_dirac_contact.py).
The contact and map are
[nsc-curvature-eft.md](nsc-curvature-eft.md). The force-to-stress frame is
[nsc-adm-source-constraints.md](nsc-adm-source-constraints.md). The raw
one-block trace is
[nsc-conformal-adm-source.md](nsc-conformal-adm-source.md). Wick
contractions use the existing finite Gaussian influence.

## Unconstrained Dirac components

In the orthonormal frame already used to contract the curvature-EFT
stress, a Dirac source has four independent entries

\[
T_{\hat0\hat0}=\rho,\quad
T_{\hat1\hat1}=p_r,\quad
T_{\hat2\hat2}=T_{\hat3\hat3}=p_\perp,\quad
T_{\hat0\hat1}=j.
\]

The owned contraction gives the trace \(\rho-p_r-2p_\perp\). Nothing in
that matrix forces the collective relations
\(p_r=\mathrm{kin}-\rho_{\mathrm{EM}}-V\) and
\(p_\perp=\rho_{\mathrm{EM}}-V\). The bulk contact remains

\[
\Delta\mathcal L_m=
-\frac{c_W}{2A^2}\left(T_{\mu\nu}T^{\mu\nu}-\frac{T^2}{3}\right)
-\frac{c_R}{4A^2}T^2.
\]

Evaluated on the Dirac-shaped matrix and factored, this is

\[
\begin{aligned}
\Delta\mathcal L_m={}&-\frac{1}{12A^2}\Big(
12 c_R p_\perp^2+12 c_R p_\perp p_r-12 c_R p_\perp\rho
+3 c_R p_r^2-6 c_R p_r\rho+3 c_R\rho^2\\
&-12 c_W j^2+4 c_W p_\perp^2-8 c_W p_\perp p_r+8 c_W p_\perp\rho
+4 c_W p_r^2+4 c_W p_r\rho+4 c_W\rho^2\Big).
\end{aligned}
\]

The same identity that removes \(-c_W C^2-c_R R^2\) was repeated on this
Dirac-shaped stress against a generic symmetric Einstein residual. The
gap is zero. Euler, box-R and the boundary variation of \(S_0\) stay
outside `stress_contact` and are not added a second time. The metric
Jacobian is not set to one, and no ghost initial datum is added.

A pure radial current \(\rho=p_r=p_\perp=0\), \(j=1\), at the probe
\(A=1,c_W=1,c_R=-2\), has contact \(1\). The collective specialization
at vanishing kinetic scalar, charge and vacuum has contact \(0\). The
two expressions are different.

## What the raw trace is

On the existing direct-ordering block, nodal partials of
\(\operatorname{Tr}(CH)\) are varied in lapse \(N\), conformal \(Q\),
sphere radius \(r\) and shift \(\beta\). The grid is the 12-point
`smooth_metric(..., general=True)` with
\(\beta=0.08\cos(2\pi x/L)\) and the fixed Gaussian covariance of that
block. Finite differences at step \(10^{-6}\) match the owned nodal
derivatives with maximum absolute gaps

| Direction | Maximum gap |
|---|---:|
| lapse | \(1.72\times10^{-10}\) |
| \(Q\) | \(1.70\times10^{-10}\) |
| \(r\) | \(1.13\times10^{-10}\) |
| shift | \(1.33\times10^{-10}\) |

The conformal chain \(F_q=F_Q/r\) and the mass-radius identity are the
owned residuals, at most about \(10^{-16}\) and \(10^{-17}\). Multiplicity
\(M=4\kappa\) at \(\kappa=1\) multiplies the mean once: the coupled energy
is exactly four times the one-block energy. It is not applied to the noise.

The owned orthonormal map, with force densities equal to nodal partials
divided by the grid spacing, produces \((\rho,j,p_r,p_\perp)\). On this
raw trace the combination \(NF_N+qF_q+rF_r\) is the massless Ward
identity of the direct ordering, so

\[
\rho-p_r-2p_\perp
=\frac{NF_N+qF_q+rF_r}{4\pi Nqr^2}
\]

is at most about \(10^{-17}\), and the gap between the projected trace and
that Ward expression is below \(10^{-12}\). This is an ordering identity
of \(\operatorname{Tr}(CH)\). It is not a renormalized Dirac trace and not
a conformal anomaly. The Weyl structure
\(T_{\mu\nu}T^{\mu\nu}-T^2/3\) built from the same raw mean does not
vanish: its maximum absolute value on the grid is \(1.53\times10^{-3}\).
Those scalars are not multiplied by a chosen \(c_W\) or \(c_R\).

The partials are not identified with \(F_A=-\delta\Gamma/\delta A\).
That identification needs the vacuum-branch variation
\(\partial_A B\), \(B=\Gamma_{\mathrm{heat}}-\Gamma_{\mathrm{canonical\,sea}}\),
counted once, plus the subtraction of the composite. Neither is evaluated.

## Mean product and Wick connected piece

For number-conserving bilinears the second moment splits as

\[
\langle AB\rangle=\langle A\rangle\langle B\rangle
+\operatorname{Tr}\bigl[C V(I-C)W\bigr].
\]

The product of means is not part of the connected trace. A one-mode
number state with occupation \(0.4\) and weights \(1.7\) and \(-0.6\) has
full moment \(-0.408\), product of means \(-0.1632\) and connected piece
\(-0.2448\), matching \(n(1-n)\) and the trace formula with gap zero.
An independent two-mode Fock density matches
\(\operatorname{Tr}[CV(I-C)W]\) with gap about \(10^{-17}\). The
connected piece there is \(-0.060196-0.01815\,i\), not the product of
the means.

On the geometric block, the representative kernels \(\partial H\) have
real connected variances \(8.380\), \(19.146\), \(0.6024\) and \(7.516\)
for lapse, \(Q\), radius and shift. Their mean squares are the separate
numbers \(0.1568\), \(0.3408\), \(3.43\times10^{-4}\) and \(0.1091\).
These are variances of the one-block generators, not a renormalized
contact.

## Angular copies and an optional sign pair

For an explicit product of \(M=4\) identical one-mode states, the
connected noise of the sum is \(4\) times one copy (\(2.7744/0.6936\)).
The same observable scaled by \(M\) inside one state has noise \(M^2=16\)
(\(11.0976/0.6936\)). The default is the product sum. Noise is not set
to \(M^2\) unless that scaled state is the state being used. The saved
one-block covariance is neither declaration.

An optional completion of the opposite angular sign is the block-diagonal
state \(C\oplus UCU^\dagger\) with \(U=\sigma_2\otimes I\). It sends
\(H(\kappa)\) to \(H(-\kappa)\) with conjugation gap zero on this grid.
It is not unique: \(U\exp(iH/\|H\|)\) conjugates the same Hamiltonian,
with gap about \(1.8\times10^{-15}\), and changes the completed state by
Frobenius distance \(0.271\). The Pauli map \(\sigma_1\otimes I\) does
not implement the sign pair; its gap is \(6.65\). The completed energy
is twice the one-block energy. That factor is the sign pair only. It is
not multiplied again by \(M=4\kappa\), and the \(2|\kappa|\) spherical
degeneracy is not constructed. On the lapse kernel, the connected noise
of this explicit product is twice the one-block noise, not sixteen times
it. The four-dimensional angular vertex and the angular current remain
absent.

## The quartic is not a Gaussian mean force

On a three-mode Fock space the generator \((\psi^\dagger V\psi)^2\),
evolved for time \(0.35\), leaves the Gaussian family. The exact state
differs from the Gaussian rebuilt from its one-body covariance by
\(0.02267\), and from the initial state by \(0.09960\). Replacing the
quartic by its scalar expectation times the identity does not move the
state (distance about \(10^{-17}\)). The mean-field bilinear
\(\langle\psi^\dagger V\psi\rangle\,(\psi^\dagger V\psi)\) stays Gaussian
and misses the exact covariance by \(0.125\). Neither replacement is
installed. The negative spectral projector of the benchmark Hamiltonian
is distance \(2.573\) from the fixed Gaussian and is not used as the
physical vacuum.

## Open relation

The missing relation is specific. The contact expectation needs

\[
F_A=\frac{1}{\Delta x}\partial_A\operatorname{Tr}(CH)+\partial_A B,
\]

with \(B\) the vacuum branch above, the heat-kernel subtraction of
\(:T_{\mu\nu}T^{\mu\nu}:\), and angular vertices that are operators rather
than the scalar \(M\). Those ingredients are not in this calculation, so
the renormalized number is not returned. The parent mean force remains
\(M\operatorname{Tr}(CH)\). This benchmark does not replace that force
and is not a prerequisite for the parent evolution.

The executable benchmark is the polynomial, the off-shell gap, the four
geometric variations, the Ward trace identity, the Fock split of mean and
connected pieces, the \(M\) versus \(M^2\) noise, the optional sign pair,
and the three-mode quartic witness. No evolution, ultraviolet campaign or
production record is run.

## Commands

The default command prints the benchmark and writes nothing:

```sh
python scripts/lab.py scripts/derive_nsc_discovery_dirac_contact.py
python scripts/lab.py -m pytest tests/test_nsc_discovery_dirac_contact.py -q
```

An exclusive JSON+NPZ record is a later root step, after every declared
source byte matches one full frozen commit. The two files together must
stay at most 64 MiB. The benchmark payload is about 33 KiB.

```sh
python scripts/lab.py scripts/derive_nsc_discovery_dirac_contact.py --write --science-commit FULL_FROZEN_COMMIT
python scripts/lab.py scripts/derive_nsc_discovery_dirac_contact.py --check
```
