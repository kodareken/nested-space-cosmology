# GMF-1B-PF1: focusing and thermal-spin-fluid preflight

GMF-1B-PF1 is a negative-result preflight for two *specified* routes toward a
regular spherical transition. It does not construct a global spacetime. It
does not reject full Einstein--Cartan--Dirac theory, EC-1's homogeneous
background algebra, noncanonical/modified gravity, or topology-changing
physics.

Reproduce the deterministic JSON record with:

```bash
python3 scripts/reproduce_gmf1b_preflight.py
```

## Canonical scalar null-focusing gate

For Einstein gravity, a canonical scalar, and a cosmological constant,

\[
T_{\mu\nu}k^\mu k^\nu=(k^\mu\nabla_\mu\Phi)^2\equiv T_{kk}\ge0.
\]

The cosmological constant does not contribute to a metric-null contraction.
For an affine, future-directed, twist-free null congruence,

\[
\frac{d\theta}{d\lambda}
=-\frac{\theta^2}{2}-\sigma_{\mu\nu}\sigma^{\mu\nu}
-8\pi G T_{kk}\le0.
\]

Consequently a negative expansion cannot rise to zero along that same future
generator before a caustic or endpoint, under precisely those assumptions.
This is not a no-go for Einstein--Cartan matter, noncanonical fields, modified
gravity, non-affine parameterizations without the corresponding term, twistful
congruences, or topology-changing/endpoint physics.

Accordingly, the generated record explicitly says that no **general Recursive
Horizons no-go theorem** has been proved. It also says that this focusing gate
does **not** support a regular bridge made from a canonical scalar under the
listed assumptions. Those are boundary statements about this route, not
evidence for an alternative global geometry.

## Naive thermal Einstein--Cartan closure gate

For the averaged ultrarelativistic thermal closure used to motivate EC-1,

\[
z=\frac{\alpha h_n^2T^2}{h_\star},
\qquad
\frac{\rho}{h_\star T^4}=1-z,
\qquad
\frac{p}{h_\star T^4}=\frac13-z,
\]

and its barotropic characteristic derivative is

\[
c_s^2=\frac{dp/dT}{d\rho/dT}
=\frac{\frac43-6z}{4-6z}.
\]

GMF-1B-PF1 classifies the exact intervals:

| Range | Result |
| --- | --- |
| `z < 2/9` | causal perfect-fluid branch |
| `z = 2/9` | zero-sound-speed degeneracy |
| `2/9 < z < 2/3` | gradient instability |
| `z = 2/3` | singular enthalpy and density-temperature derivative |
| `z > 2/3` | superluminal relative to the metric cone |

At an EC-1 lower turning point, with

\[
\dot a^2=-1+\frac A{a^2}-\frac B{a^4},
\qquad A,B>0,\quad A^2>4B,
\]

the stable roots obey

\[
x_{\max}=\frac{A+\sqrt{A^2-4B}}2,
\qquad x_{\min}=\frac B{x_{\max}},
\]

and the thermal mapping is

\[
z_{\rm bounce}=\frac{B}{Ax_{\min}}=\frac{x_{\max}}A
=\frac12\left(1+\sqrt{1-\frac{4B}{A^2}}\right).
\]

Therefore `1/2 < z_bounce < 1`: the EC-1 bounce never lies in
`0 < c_s^2 <= 1` for this naive thermal averaged perfect-fluid closure. A
continuous low-density branch must pass through the unstable or singular
regime before it reaches that bounce.

For very small `B/A^2`, binary64 may display `z_bounce` as `1.0`; the record
also carries the positive stable complement `one_minus_z_bounce=B/(A*x_max)`
and the strict interval conclusion is algebraic, not inferred from rounded
float comparison.

This rejects **only** that material-law closure as a healthy GMF-1B medium.
It does not reject full Einstein--Cartan--Dirac dynamics or the exact EC-1
background identities. It derives no external exchange `Q`, late dark energy,
or variable local `c`.

Primary sources: [A. Raychaudhuri, “Relativistic Cosmology. I,” *Physical
Review* 98, 1123 (1955)](https://doi.org/10.1103/PhysRev.98.1123); [N.
Popławski, *Gravitational Collapse with Torsion and Universe in a Black Hole*
(2023)](https://arxiv.org/abs/2307.12190); and [A. H. Ziaie, P. V. Moniz, A.
Ranjbar, and H. R. Sepangi, *Einstein--Cartan gravitational collapse of a
homogeneous Weyssenhoff fluid* (2014)](https://doi.org/10.1140/epjc/s10052-014-3154-2).
