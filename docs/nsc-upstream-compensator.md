# One common upstream metric compensator

Reuse the [fixed-interval rank certificate](nsc-upstream-compensation-rank.md),
the frozen unitary radial generator and the exact first-order Fourier pair
formula. The new experiment asks whether one actual metric compensator,
calibrated at one diagonal pair, also matches other pairs and channels.
The three unknowns are raw lapse, axial scale and radius coefficients of
one geometric bump. They are not action couplings or independently supplied
matrix forces.

The bump is `b(rho)=profile(200*rho−204)`, centered at51/50 with halfwidth1/200.
It vanishes on an open neighborhood of rho1, so all its incoming geometric
jets are exactly zero. Retain the original radius pulse, axial bump w(z),
signed energies, source law and scales. The total raw tangent is
`(deltaN,deltaa,deltar)=(cN*b*w,ca*b*w,s*chi(s)*w+cr*b*w)` with betaKS=0.
The same real triple acts on every frequency and angular/compact channel.

For Ebar=(Eo+Ei)/2 the three Weyl vertices are

\[
V_N=-m\sigma_1+\ell\sigma_2/r-\bar E\sigma_3/a,\qquad
V_a=\bar E\sigma_3/a^2,\qquad V_r=-\ell\sigma_2/r^2.
\]

The fixed chart gives forcing `(i/a0)*w_hat(Eo−Ei)*b*Vj*Ui`; the reference
clock factor a0 is not varied. Ebar retains the axial Weyl derivative terms.
Ui starts at I2 as an operator fundamental, not as a quantum state. Four
direction responses (original radius and N/a/r) are integrated together,
and `K_oi=B_oi*Ui†`. Operator CAR requires `K_oi+K_io†=0`.

At14_1 and the positive selected energy0.5562120090641313, check that each
diagonal K is traceless anti-Hermitian, then solve its real3x3 Pauli system
directly. Freeze the fine coefficients before any validation channel.
Exactly six short solves are allowed: two14_1 tolerances(2e-10/2e-12,
2e-13/2e-15), followed by one fine solve each for14_-1,13_1,0_1,24_1.
Maximum radial step is1/2000; total CPU limit60 seconds. Stop at these cases,
PASS or OPEN, with no adaptive experiment extension.

## Source covariance and verification

All canonical stationary columns f come from authenticated801-node history
fields and their original complete source fibers. No reference resampling,
horizon or scattering solve occurs. Form `C_i=f_i*Csrc_i*f_i†` without
normalization. In this unchanged-upstream slab `deltaA_oi=K_oi*f_i`, hence
`deltaChat_oi=K_oi*C_i+C_o*K_io†`, including open fibers. This is a local
factorization, not permission to drop source-complement terms when upstream
preparation changes. No energy weight or angular multiplicity is inserted.

The checks use3e-11 for local vertex algebra, unitarity, CAR, selected
cancellation, fixed-fine-coefficient coarse consistency and agreement of
the original14_1 radius response with the frozen radial artifact. Report
every uncompensated and compensated pair, including channels whose original
radius response vanishes. Background preparation accuracy remains inherited;
numerical zero at selected pairs is not a rigorous matching certificate.

The union support retains the original retarded past/inflow inequalities.
A conservative strictly positive range for the overall family amplitude eta
keeps rawN,a,r positive and qPG²>0. Its rational construction uses|eta*cN|<.01,
|eta*ca|<.008 and|eta|*(.03+|cr|)<.5. It is a local chart-validity range,
not physical initial data or a chosen evolution interval.

Likely failure modes are a vertex/frame mismatch, inaccurate narrow-bump
transport, and failure of the frozen triple at uncalibrated channels. Each
is reported separately; no source adjustment or channel-dependent refit
is allowed. The [record](../results/development/nsc-upstream-compensator.json)
and saved arrays retain the actual coefficients, residual pairs, all
provenance and CPU accounting. Global C0 matching and the common action
constraints remain OPEN whenever residuals remain, even if solver checks
pass. No state correction, normalization, new interaction or metric timestep
is introduced.

```sh
python3 scripts/derive_nsc_upstream_compensator.py --run
python3 scripts/derive_nsc_upstream_compensator.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_upstream_compensator.py
```
