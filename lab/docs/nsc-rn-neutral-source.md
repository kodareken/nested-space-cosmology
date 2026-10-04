# Neutral RN scattering packet and common-action source

This owner prepares one neutral massless Dirac packet on a static
exterior Reissner–Nordström background and evaluates that packet's
stress. The geometry is the ingoing flat-spatial chart

\[
ds^2=N^2dt^2-(q\,dx+\beta\,dt)^2-r^2d\Omega^2,
\]

signature \(+---\), with \(c=\hbar=1\). The areal section used for the
modes is \(x=r\) and \(q=1\). That value of \(q\) is not the coordinate
gauge \(q=1/r\). Metric forces are taken before that identification.

The background is classical Einstein–Maxwell. With
\(f=1-2M/r+Q^2/r^2\),

\[
\beta=\sqrt{2M/r-Q^2/r^2},\qquad N=\sqrt{f+\beta^2}.
\]

For this flat slicing the second expression equals \(1\). The lapse array
is that root; the Dirac operator reads the array. \(r_m=\sqrt{P^2}\) is
the magnetic radius, not the mass. The Maxwell field is held fixed at the
classical monopole. The integer-\(q=1\) chart rescales \(C_F\) so that
\(P^2\) matches the audited flux-4 geometry; it is not that quantized
charge sector. The Dirac shell carries no electric current. Extremal
holes are outside this regular branch.
There is no Weyl variation and no second copy of the Maxwell stress.

## Characteristic Dirac operator

The canonical half-density is \(u=r\sqrt{q}\,\psi\). In the Pauli frame
fixed by the tetrad reduction,

\[
H=-i\big(C\partial_x+\partial_x C\big)/2+( \kappa N/r)\,\sigma_1,
\qquad
C=(N/q)\,\sigma_2-(\beta/q)\,I.
\]

\(\kappa/r\) is the angular eigenvalue, not a four-dimensional mass.
The unitary \(U\) with columns \((1,i)/\sqrt2\) and \((1,-i)/\sqrt2\)
puts \(\sigma_2\) on \(\mathrm{diag}(1,-1)\). Speeds in that basis are

\[
a_+=(N-\beta)/q,\qquad a_-=-(N+\beta)/q.
\]

Outside the outer horizon \(a_+>0\) and \(a_-<0\). On the future horizon
\(a_+\to0\). A mode of Killing frequency \(E>0\) is integrated from the
regular root

\[
\chi_+=\frac{\kappa N/r}{iE-a_+'/2}\,\chi_-,
\]

then sampled on a nonperiodic summation-by-parts grid. This is not a
positive spectral projection of an absorbing matrix, and it does not
multiply a waveform by such a projector. Negative-frequency contamination
of the resulting span is zero because no negative root is admitted.
The old periodic black-universe box and its \(\kappa=14.5\) calibration
are not this background.

## Packet

The default center, width and frequency are \(12r_m\), \(2r_m\) and
\(\omega r_m=1\), with \(r_m=\sqrt{P^2}\). Spectral samples lie on the
positive-frequency window set by

\[
\sigma=|a_-(r_{\mathrm{center}})|/\mathrm{width}.
\]

Each mode is weighted-normalized. One constant phase makes its ingoing
characteristic real and positive at the packet center; the radial transport
phase is kept, so the Gaussian spectral sum interferes there. Its
filling \(\nu\in(0,1)\) is fixed only after the measured Killing energy
\(E_{\mathrm{packet}}\) of the normalized spinor is known:

\[
\nu=\frac{\varepsilon M}{G_N\,4\,E_{\mathrm{packet}}},\qquad
\varepsilon=G_N\frac{4\nu E_{\mathrm{packet}}}{M}\in\{10^{-3},10^{-2}\}.
\]

The factor \(4\) is the isotropic \(\kappa=1\) shell, applied once to the
stress and the energy. It is not inserted into the per-channel
anticommutator. That anticommutator is the weighted Gram of the mode
columns. On a finite exterior interval those modes need not be orthogonal:
the Gram is stored, not replaced by the identity. The occupation matrix on those columns is the rank-one packet
matrix, and its diagonal is the Gaussian spectral population. The angular
covariance is \(\nu I_4\), with trace \(4\nu\).

`prepare_radial_packet` returns the production column `phi` of shape
`(2, points, 1)` with occupation \([\nu]\). The scattering modes stay beside
that column. The per-channel CAR is \(G^{1/2}CG^{1/2}\) and is exactly
rank one for this packet. Multiplicity 4 multiplies the stress once and
does not enter the CAR. Constructors do not write evidence files.

## Source

For a per-channel matrix \(C_{jk}\),

\[
I_\pm=\sum_{jk}C_{jk}\,
\frac{\chi_{\pm j}^{*}\partial\chi_{\pm k}
-(\partial\chi_{\pm j})^{*}\chi_{\pm k}}{2i},
\qquad
B=\sum_{jk}C_{jk}\,\chi_j^{\dagger}\sigma_2\chi_k.
\]

One channel then has the flat-measure derivatives

\[
\begin{aligned}
F_N&=(I_+-I_-)/q+(\kappa/r)B,\\
F_\beta&=-(I_++I_-)/q,\\
F_q&=-(a_+I_++a_-I_-)/q,\\
F_r&=-(\kappa N/r^2)B.
\end{aligned}
\]

The shell multiplies these by \(4\) once. Coordinate energy density is
\(N F_N+\beta F_\beta\). Angular-integrated densities divide by \(q r^2\)
only. Physical orthonormal densities divide by \(4\pi q r^2\):

\[
\rho=\frac{F_N}{4\pi q r^2},\qquad
j=-\frac{F_\beta}{4\pi q r^2},\qquad
p_r=\frac{F_N-(\kappa/r)B_{\mathrm{shell}}}{4\pi q r^2},\qquad
p_\Omega=-\frac{F_r}{2 N q r\,4\pi}.
\]

The massless projection identity is \(\rho-p_r-2p_\Omega=0\). The radial
Noether current is \(a_+|\chi_+|^2+a_-|\chi_-|^2\), constant on a
stationary mode. Electric current is the zero array. No force with
respect to magnetic flux is formed.

The same evaluator accepts a general positive \((N,\beta,q,r)\) and the
caller’s SBP derivative. It does not allocate a dense spatial covariance.

## What this slice does not do

It does not evolve the packet, choose a parent or child window, impose
excision between the horizons, or run a campaign. Charged Landau modes
and the electric Gauss constraint are a later experiment. Horizon-thermal
and MMP covariances are not copied in.
