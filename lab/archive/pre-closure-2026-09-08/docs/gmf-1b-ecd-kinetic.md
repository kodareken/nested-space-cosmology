# GMF-1B-ECD-KIN1: contact tetrad variation and free-Dirac control

KIN1 is a corrective gate between the local interaction identities and the
still-open full ECD initial-data problem. It retracts the earlier fixed-current
NEC claim and establishes the physically valid contact-sector result before
the programme proceeds to ECD-ID1.

## Corrected action bookkeeping

INT1 already fixes the full action after the independent connection has been
eliminated everywhere:

```text
canonical (+---):  L_4 = -(3 kappa/16) J^I J_I
VC (-+++):         L_4 = +(3 kappa/16) A^I A_I
J_+^2 = -A_-^2
```

For the exact `F=3/4`, `G=1/2`, `C=1`, `kappa=1` fixture,
`J^2=+9/4`, `A^2=-9/4`, and both representations give `L_4=-27/64`.
Varying this quartic interaction with respect to a spinor produces the
Hehl--Datta coefficient `3 kappa/8`; the factor of two comes from varying the
quadratic current product. It is not a second action contribution.

The original KIN1 added a raw Dirac-contortion insertion to a quantity already
identified as the reduced contact interaction. Those are
convention-dependent intermediate pieces of the unreduced connection action.
They cannot be added again after the connection has been eliminated. The
correct artifact therefore imports the full reduced sign from INT1 and records
that no separate connection pieces were double counted.

## Why the fixed-current NEC term was false

The axial current entering the reduced action is an internal Lorentz vector,

```text
J^I = bar(psi) Gamma^I Gamma5 psi,
J^2 = eta_IJ J^I J^J.
```

With spinor components fixed, `J^I` and `J^2` do not depend on the tetrad. In
coordinate components,

```text
J_mu = e^I_mu J_I,
J^2 = g^munu J_mu J_nu.
```

The two coordinate factors vary together. For every coframe perturbation,

```text
delta(J^2)
  = J_mu J_nu delta(g^munu)
  + 2 g^munu J_mu delta(J_nu)
  = 0.
```

The executable uses a generic non-diagonal coframe and a generic internal
current, checks all `16` independent `delta e^I_mu` directions, and finds a
nonzero fixed-coordinate-current metric partial (`4.5493986797...`) canceled by
the current response to a maximum residual below `2e-15`.

The earlier expression

```text
T_ab = c J^2 g_ab - 2 c J_a J_b
```

is valid only for a separate coordinate covector artificially held fixed. It is
not the tetrad variation of the Dirac axial current and its claimed negative
null contraction is retracted.

## Physical contact stress

Because the internal invariant is tetrad independent, the algebraic contact
sector varies only through the determinant:

```text
T_ab^(4) = L_4 g_ab.
```

For the VC fixture this gives

```text
rho = 27/64,
p = -27/64,
w = -1,
rho + p = 0,
T_ab^(4) k^a k^b = 0
```

for every metric-null `k`. KIN1 checks two independent null directions. The
contact term is therefore NEC-saturating; it supplies neither the earlier
local defocusing claim nor the homogeneous `rho+p<0` bounce. No isotropic
average of the retracted tensor is retained.

This result agrees with the torsion-eliminated action and stress displayed in
[Choudhury, Maity, and Lahiri (2024), Eqs. 7--10](https://link.springer.com/article/10.1140/epjc/s10052-024-13618-4).
The separation between an anisotropic intermediate Dirac expression and the
metric-proportional total torsion contribution is also explicit in
[Cabral, Lobo, and Rubiera-Garcia (2019), Eqs. 20--24](https://arxiv.org/abs/1902.02222).

## Free-Dirac control

The intended next calculation begins with the torsion-free Dirac sector. KIN1
now includes only the smallest safe control: in the canonical `(+---)` gamma
representation,

```text
p^a = (1,0,0,1),
u   = (1,0,1,0),
gamma^a p_a u = 0,
j^a = bar(u) gamma^a u = 2 p^a,
T_ab = (p_a j_b + p_b j_a)/2 = 2 p_a p_b.
```

The parallel null probe has contraction `0`; the opposite null probe has
contraction `8`. This verifies one positive-frequency massless plane wave. It
does **not** prove the null energy condition for arbitrary classical Dirac
fields or for the nonlinear spherical two-spinor configuration.

## What remains ECD-ID1

KIN1 completes the tetrad variation of the algebraic contact sector only. The
following remain open:

- the full torsion-free Dirac kinetic/spin-connection tetrad source;
- complete spin-weighted harmonic and SO(3) source closure;
- a regular-centre series and finite initial stress, curvature, and torsion;
- the Einstein constraints and a constraint-compatible exact vacuum buffer;
- finite Misner--Sharp mass, conservation, evolution, and the actual
  focusing/defocusing history.

Accordingly the result keeps `full_ecd_effective_stress_verified=false`,
`general_classical_dirac_nec_proven=false`,
`metric_null_defocusing_derived=false`, and `bounce_constructed=false`.

## Reproduction

```bash
python3 scripts/reproduce_gmf1b_ecd_kinetic.py \
  --output results/gmf-1b-ecd-kinetic.json
```

The exact fixture and scope flags are audited by `scripts/check_repo.py`.
