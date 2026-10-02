# Constrained flat curvature sector of the current local action

The [symbolic owner](../src/recursive_horizons/nsc_discovery_curvature_sector.py)
derives the flat-patch equations from the actual first-order density,
canonical momenta and lapse/shift constraints. It assesses the current
local Einstein–Weyl truncation **when treated as exact dynamics**. It does not
change a coefficient, source, force, trajectory or older record. The calculation
is not periodic whole-domain initial data or an all-source instability proof.

## Exact flat patch and constrained equations

The analytic domain has $x>0$, $r=x$, $Q=L=1/x$, beta zero, chi zero,
zero canonical momenta, zero field covariance and zero magnetic flux. The
physical metric is $ds^2=dt^2-dx^2-x^2d\Omega^2$. Its Euler equations and
both constraints vanish. In the owned Riemann convention, $R_h=2$ and
the four-dimensional curvature is zero.

Set $g=8\pi A$, $f=2\alpha$, $a=\delta r$, $q=\delta\log Q$,
$c=\delta\chi$ and $\Box=\partial_t^2-\partial_x^2$. The owned canonical
momentum linearization is

$$
\delta p_\chi=2f q_t,\qquad
\delta p_r=-6g a_t-2gxq_t,\qquad
\delta p_Q=2fxc_t-2gx^2a_t.
$$

The radial and auxiliary equations give

$$
\Box a=\frac c{6x},\qquad
\Box q=-\frac{c+4q}{2x^2}.
$$

Define $b=a_x-a/x-q$ and $T_c=c_x/x+c/x^2$. The constraints are

$$
\frac{\delta C}{2gx^2}=b_x+\frac bx+
\frac fg\left(\frac c{x^3}-\frac{c_x}{x^2}-\frac{c_{xx}}x\right),
\qquad
\frac{\delta D}{2gx}=\partial_t\left(b-\frac fgT_c\right).
$$

Their joint solution is $b=(f/g)T_c+K/x$, with K constant in both
coordinates. The remaining metric equation then becomes

$$
\Box c+\frac{2c_x}{x}+\frac{4c}{x^2}-m_*^2c
=-\frac{2gK}{fx},\qquad
m_*^2=\frac g{6f}=\frac{2\pi A}{3\alpha}=-\frac A{2C_W}.
$$

The static Einstein/Coulomb component is $c_E=12K/x$, with
$K=-\delta M$ for the linear Schwarzschild mass. Removing it and writing
$\psi=(c-c_E)/x$ gives

$$
\boxed{\psi_{tt}-\psi_{xx}+\left(\frac6{x^2}-m_*^2\right)\psi=0.}
$$

The $6/x^2$ radial tensor/tidal potential is retained. Chi does not obey a
bare scalar massive-wave equation at finite radius. The scalar auxiliary's
background is identically zero, so its linear perturbation is invariant
under an infinitesimal coordinate transformation. The propagating chi
component is nongauge. The symbolic owner also checks that its wave equation,
the reconstructed constrained q and the forced radial equation satisfy the
auxiliary equation together; the pole is not inferred from an off-constraint
force.

## Proper-flat symbol and physical sign

At large R, retain the locally flat spherical jets $r_x=1$,
$Q_x=-1/R^2$ and their derivatives before taking the limit. Do not replace
the patch by a constant-radius cylinder. In variables
$[q,a/R,c/R^2,p_Q/R^3,p_r/R,p_\chi]$, the six-variable symbol has

$$
\det(\lambda I-J)=(\lambda^2+k^2)^2
 (\lambda^2+k^2-m_*^2).
$$

Its extra eigenvector satisfies the limiting constraints, for example
$q=-3f\bar c/g$, $a/R=f\bar c/g$, $p_Q/R^3=p_r/R=0$ and
$p_\chi=-6f^2\lambda\bar c/g$, where $\bar c=c/R^2$.
The physical metric has $N=rQ=1$, so proper time equals t and radial
proper distance equals x. Consequently

$$
\omega_{\rm proper}^2=k_{\rm proper}^2-m_*^2
$$

in the large-radius limit. Current locked coefficients give
$m_*^2=26.6666666667$, growth rate $m_*=5.163977795$ for long waves,
and proper growth time $m_*^{-1}=0.1936491673$ in the declared local-action
units. A finite-R approximation additionally contains $6/R^2$. These proper
flat scales are not clock measurements from a curved production trajectory.

The [owned sign convention](nsc-spherical-action.md) is signature +---,
$R^a{}_{bcd}=\partial_c\Gamma^a{}_{db}-\partial_d\Gamma^a{}_{cb}
+\Gamma^a{}_{ce}\Gamma^e{}_{db}-\Gamma^a{}_{de}\Gamma^e{}_{cb}$,
and $R_{bd}=R^a{}_{bad}$. Eliminating chi returns
$\alpha\sqrt{|h|}(R_h-2)^2$. With the owned
$C^2=(R_h-2)^2/(3r^4)$ identity, the four-dimensional local action is

$$
S=\int\sqrt{|g|}\,[-A R-C_W C^2],\qquad
\gamma_W=-C_W.
$$

Thus the current negative C_W means a **positive actual Lorentzian Weyl
coefficient**. The flat spin-2 operator factors, up to its overall
normalization, as

$$
\Box\left(\Box+\frac A{2C_W}\right)=\Box(\Box-m_*^2).
$$

The extra pole's physical squared mass is $\mu_2^2=A/(2C_W)<0$.
Its residue has the opposite sign from the Einstein pole. A negative mass
square gives a tachyonic branch; the opposite residue is the separate ghost
property of this local fourth-order truncation. The declared action and
variational implementation agree on the sign. Changing a Riemann label
alone does not change the physical pole of a consistently transformed action.

## Existing perturbative EFT branch and limits

The separate [curvature EFT owner](nsc-curvature-eft.md) retains the curvature
correction only through its declared first order, with its source interaction,
physical-metric map and boundary/state requirements. That branch uses the
leading theory's initial data rather than adding the extra modes of a resummed
fourth-order truncation. This assessment does not rederive or implement it.

Whether the extra pole belongs within the valid spectral/EFT range requires
the matched Lorentzian action and controlled nonlocal remainder. Neither the
current spectral CAR positivity nor positive initial lapse-source density
validates a healthy physical mass sign for the curvature pole. This result
does not prove full nonlocal mode health, all-source instability, failure of
bounded regimes, or an NSC verdict. Existing production records remain in
their original finite domains.

## Reproduction and immutable record

Default calculation is read-only and symbolic. No initial solve, time march,
campaign or old evidence rewrite occurs. Source and coefficient hashes are
captured before calculation and checked afterward. Explicit creation requires
a frozen producing commit matching all source blobs. The writer uses an
exclusive create and refuses unfrozen or unresolved records. The checker
authenticates producing inputs and repeats the symbolic identities without
writing.

```sh
.venv/validation/bin/python scripts/lab.py \
  scripts/assess_nsc_discovery_curvature_sector.py
.venv/validation/bin/python scripts/lab.py \
  scripts/assess_nsc_discovery_curvature_sector.py \
  --producer-commit HEAD --write
.venv/validation/bin/python scripts/lab.py \
  scripts/assess_nsc_discovery_curvature_sector.py --check
```

The [tests](../tests/test_nsc_discovery_curvature_sector.py) check the owned
constrained identities, pole signs, deliberate inconsistent-action detection,
and creation-only replay in temporary files. Root owns source freeze and
production record generation.
