# Short-shell radial Fourier response

Reuse the owned canonical KS frame and the authenticated group14_1 columns.
The new connection is an independent check of the
[spacetime Fourier pairs](nsc-retarded-fourier-response.md), without another
global PDE, horizon, scattering or source-state calculation. Its equations are

\[
A_i'=G_iA_i,\qquad
B_{oi}'=G_oB_{oi}-\frac{i\ell s\chi(s)\widehat w(E_o-E_i)}{a r^2}\sigma_2A_i,
\quad G_E=\frac{i}{a}\left(-m\sigma_1+\frac\ell r\sigma_2-\frac E a\sigma_3\right).
\]

Here a=sqrt(beta²−1), s=T(rho)−T(1), and the Fourier convention is
`w_hat(omega)=integral exp(i*omega*z)*w(z) dz`. The inherited frame/clock
identity cancels the radial phase from the forcing. B_oi(1) is the axial
Fourier amplitude itself. No2pi, energy quadrature weight, occupation or
angular multiplicity enters this ODE or the first-order pair contraction.
The source covariance in that contraction remains the original full matrix.

Take the first saved3201-grid node at or above1.03, outside the exact pulse
support, recording its actual coordinate. Convert the saved characteristic
columns with MODE_TO_CURRENT, the owned trace inverse/r and exp(i*E*S).
The unchanged upstream preparation sets B=0. Keep four signed energies and
all three source columns. Use the unchanged normal cutoff and axial bump.

The bounded decision is two short DOP853 solves with rtol/atol
(2e-10,2e-12) and(2e-13,2e-15), using the192-node scalar Gauss rule.
Compare that Fourier rule with96 nodes as a separate quadrature indicator.
Total CPU budget is30 seconds, including preparation. Require local
PG/KS generator and forcing-map residuals, homogeneous endpoint mismatch,
radial refinement and CAR-pair residuals below3e-11. Compare every mode and
covariance pair against the saved fine spacetime result at1% relative.
Stop after the two solves, PASS or OPEN, without expanding the frequencies.

Likely failure modes are a frame/sign inconsistency, remaining spacetime
discretization error, and a radial tolerance/CPU failure. Each is reported
at its named check without fitting the pulse or hiding an unresolved term.

The [record](../results/development/nsc-retarded-radial-response.json) stores
actual comparisons and authenticated source inputs. It also reports the
diagonal K=B f† anti-Hermitian residual and independently K=B_active/f_active,
without normalizing columns. Positive-energy covariance pairs are separated
mathematically into the Pv=diag(0,1,0) contribution and the retained remainder;
the actual full covariance response is unchanged. No physical remainder
bound is inferred from this numerical separation.

This is conditional numerical agreement. The low-energy physical input was
prepared using a finite horizon working frame and a subgap reflection phase;
its physical preparation error remains OPEN. Tolerance differences are not
rigorous bounds. No fullC0 exclusion, source replacement, action modification,
stationarity or metric-solution claim follows. Source fibers and all Gamma
terms are retained. Replay contracts saved arrays only.

```sh
python3 scripts/derive_nsc_retarded_radial_response.py --run
python3 scripts/derive_nsc_retarded_radial_response.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_retarded_radial_response.py
```
