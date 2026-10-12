# Paired incoming source tail in the retained mode approximation

The [finite incoming source](nsc-incoming-spectral-source.md) already fixes
the common KS surface, state and fourth-order subtraction. This calculation
adds its missing energy tail for all 32 retained non-LLL groups. It reuses the
same order-16 Riccati mode approximation and adiabatic Bloch expressions;
it runs no horizon, scattering, packet or metric evolution.

For each group, the actual angular signs are summed before expansion in
$x=1/E$. The existing multiplicity is
$n_{\mathrm{copies}}d/(4\pi^2r_1^2a_1n_{\mathrm{signs}})$, with
$r_1^2=2$ and $a_1^2=3\pi/2-4$. The archived split is 160, except for
groups 10–12 and 31–32, whose split is 320. These are the existing integration
endpoints, not new physical cutoffs.

## Endpoint identity before numerical differentiation

Write $L=\lambda/r$, $A=-a^2$. The existing Riccati recurrence begins with

$$
c_1=L+im,\qquad c_2=-mA_\rho/2+i(-a^2L_\rho+A_\rho L/2).
$$

The common chart gives $\dot a=A_\rho/2$ and $\dot L=-aL_\rho$.
The vacuum Bloch vector and the zeroth reference therefore have equal leading
transverse coefficients and equal longitudinal coefficients through $x^2$.
The higher adiabatic terms have transverse order at least $x^2$ and
longitudinal order at least $x^3$. Each subtracted stress kernel is consequently
$x^2g(x)$; its coefficients of $x^0$ and $x^1$ vanish structurally.

The implemented exact coefficient check further compares the transverse
$x^2$ and longitudinal $x^3$ terms with the already-owned first adiabatic
recursion $b_1=-h\times\dot b_0/(2|h|^2)$. Their three symbolic residuals are
zero separately for each angular sign. This establishes $g(0)=0$.
The owner differentiates this regular function and prepends the two exact
structural zeros. It uses `chop=False`; no small coefficient is thresholded
away. The tiny first derivative returned by differentiating the unfactored
floating expression is retained as an arithmetic diagnostic, not interpreted
as a nonintegrable physical tail.

## Tail integration and verification

For the resulting coefficients $a_n$, the tail of the specified approximation
is integrated as

$$
\int_{E_{\max}}^\infty a_n E^{-n}\,dE
=\frac{a_n}{(n-1)E_{\max}^{n-1}}.
$$

Orders7,9,11 and13 are retained. An independent evaluation of the same
arbitrary-precision expressions integrates $g(x)$ over
$0<x<1/E_{\max}$ with16 and24 Gaussian nodes, at50 and60 decimal digits.
The record separates order differences, quadrature/precision differences and
agreement with already stored resolved source bands. These are numerical
controls for the defined approximation, not a rigorous bound on the exact
Dirac mode remainder.

The immutable coefficient artifact also contains the raw unfactored
derivative diagnostic. `--check` authenticates all inputs, reconstructs every
recorded field from that artifact and the existing source panels, and does
not rerun coefficient generation or old scientific generators.

## Thermal bound and remaining physical scope

The same horizon covariance and incoming occupation obey

$$
\|C_{\mathrm{src}}-C_{\mathrm{vac}}\|_1
\le 2e^{-\pi E/\kappa_h}
+e^{-2\pi E/(\Omega\kappa_h)}.
$$

The first term retains the horizon-pair coherence. Current-normalized
compression cannot increase this trace norm. The existing normal-energy,
pressure and current vertices then give elementary exponential tail bounds;
the JSON stores their natural logarithms and decimal strings so numerical
underflow is not reported as an exact zero. These bounds concern the thermal
state correction, not the vacuum mode approximation.

The physical remainder of the order-16 Riccati approximation remains separate
and OPEN. So do independent subgap refinement, retained angular/compact
truncation, complete source matching and extended stationarity. This result
does not assign new constraints, null stresses, a selected geometry or a
metric timestep. The LLL allocation, scales, seeds, earlier records, PDF and
public repository remain unchanged.

```sh
python3 scripts/derive_nsc_incoming_source_tail.py --prepare
python3 scripts/derive_nsc_incoming_source_tail.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_source_tail.py
```
