# Closed coefficient for the incoming third-normal lapse response

For `u=delta(r_TTT)`, the already declared local/reference action gives

$$
c_u=c_{u,\mathrm{ref}}+c_{u,\mathrm{local}}
\simeq 0.1012760934543-0.0438768158680
=0.0573992775863>0.
$$

This coefficient belongs to the [conditional constraint plane](nsc-incoming-constraint-germ.md).
It is fixed by the existing action, not a fitted gravitational or source
parameter. A directed interval evaluation establishes its positive sign
without evaluating the physical matter source or choosing an initial geometry.

## Reuse and missing connection

Reuse the existing homogeneous adiabatic Bloch recurrence, retained channel
weights and local action densities. The small `u=0.01` contour probe resolves
an independent lapse direction, but dividing its small difference amplifies
roundoff in the coefficient. The missing connection is a closed coefficient
with controlled input interpretation. Integrate the existing fourth-order
reference coefficient with the standard beta-function moments and extract
the homogeneous local mixed density derivative. Stop after verifying these
identities, their directed enclosure and agreement with the owned functions;
no reference, radial, physical mode or metric evolution is needed.

Risks are a factor-two/factorial error, hidden change of the locked coefficient
interpretation, and attributing probe cancellation error to new physics.
Full-line normalization, independent density checks and explicit exact/float
input hulls keep these separate.

## Reference contribution

Let `h=(-m,L,p)`, `L=ell/r`, `p=k/a`, `M^2=m^2+L^2`, `H_a=a_T/a`
and `H_r=r_T/r`. The reused Bloch recurrence and normalization give

$$
\delta(h\cdot b_4)
=-\frac{(h\times\dot h)\cdot(h\times\delta h''')}{16\omega^7},
\quad
(h\times\dot h)\cdot(h\times\delta h''')
=\frac{L^2u}{r}[m^2H_r+p^2(H_r-H_a)].
$$

Using the ordinary full-line moments `I0=16/(15 M^6)` and `I2=4/(15 M^4)`
with `dk=a dp` and the existing `d/(2 pi)` measure yields

$$
c_{u,\mathrm{ref},g}
=-\frac{d a L^2}{120\pi r}
\left[\frac{4m^2H_r}{M^6}+\frac{H_r-H_a}{M^4}\right].
$$

`d=copy_count*degeneracy` already includes both angular signs. No additional
frequency folding is inserted. The LLL reference response is exactly zero;
massive zero-angular channels also contribute zero without a singular formula.

## Local contribution

For homogeneous raw KS fields, the owned curvature densities have no genuine
`N_TT` dependence. Thus `c_u,local=-partial_N_T partial_r_TT L`. The existing
second-derivative densities give

$$
c_{u,\mathrm{local}}=ar\left[
-\frac{(h_q+\log r)(H_a-H_r)}{30\pi}
-\frac{H_a+H_r}{60\pi}
+\frac{32\pi C_W}{3}(H_a-H_r)\right].
$$

The three terms collect barred radial curvature/WZ Weyl, WZ BoxR and compact
Weyl. Euler, Einstein, cylinder, gauge and LLL-normal-order terms give zero
to this mixed coefficient. The harmonic is the existing
`h_q=EulerGamma/2-H4/4`, not a new normalization. A symmetric four-corner
stencil of the **original density owners** is exact at this derivative grade
and independently checks the formula; its floating arithmetic remains visible.

The [record](../results/development/nsc-incoming-lapse-coefficient.json)
encloses exact geometric/harmonic expressions and their stored evaluations.
The unchanged binary `C_Weyl` ledger value is the declared coefficient used
here; no unprovided error interval for a different exact coefficient is assumed.
The old small-probe coefficient remains recorded as a numerical control.
Source accuracy, the second response coefficient's rigorous enclosure,
physical initial data, a surface solution and extended stationarity remain OPEN.

```sh
python3 scripts/derive_nsc_incoming_lapse_coefficient.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_lapse_coefficient.py
```
