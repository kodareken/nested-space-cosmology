# GMF-1B-ECD-TOR1: spin tensor, torsion, and effective-contact cross-check

TOR1 is the algebraic bridge from SYM1's axial current to the
torsion-eliminated interaction fixed by INT1. It derives the spin tensor and
torsion invariant directly, then cross-checks the full reduced contact scalar
without assigning that scalar to one intermediate connection term.

## What it derives

For the Ventrella--Choptuik left-chiral spherical pair in the VC `(-+++)`
representation, the totally antisymmetric Dirac spin tensor is

```text
S^{abc}
  = -(1/2) sum bar(psi) gamma^{[a} gamma^b gamma^{c]} psi
  =  (1/2) eps^{abcd} A_d.
```

The executable evaluates both sides independently. For `F=3/4`, `G=1/2`, and
normalization `C=1`, the nonzero control components include
`S^{012}=5/16` and `S^{123}=13/16`, with residual at machine precision.

The algebraic Cartan relation then gives

```text
T^{abc} = kappa S^{abc},
T_abc T^abc = -(3/2) kappa^2 A^2.
```

For `A^2=-9/4` and `kappa=1`, the torsion invariant is `27/8`.

## Full reduced interaction, not an isolated piece

After the independent connection is solved and substituted into the complete
Einstein--Cartan--Dirac action, the VC interaction is

```text
L_4 = (3 kappa/16) A^2 = -27/64.
```

TOR1 cross-checks this value against INT1. It no longer claims that the
Einstein--Hilbert torsion-squared term alone equals the full Hehl--Datta
interaction. The gravitational and Dirac connection terms are
convention-dependent intermediate contributions of the unreduced action; the
full coefficient is meaningful only after the connection has been eliminated
everywhere. Adding an intermediate term to the already-reduced contact scalar
would double count it.

## Contact-sector tetrad stress

The internal scalar `A^I A_I` is tetrad independent at fixed spinor
components. Therefore the complete variation of the algebraic contact sector
is

```text
T_ab^(4) = L_4 g_ab.
```

In the flat VC fixture,

```text
rho = 27/64,
p = -27/64,
w = -1,
rho + p = 0,
rho + 3p = -27/32.
```

This is a positive cosmological-constant-shaped local contact contribution,
whose magnitude depends on the unconstrained spinor amplitude. Its null
contraction is zero, so it does not establish metric-null defocusing or a
bounce.

The action and stress form are independently displayed in
[Choudhury, Maity, and Lahiri (2024), Eqs. 7--10](https://link.springer.com/article/10.1140/epjc/s10052-024-13618-4)
and in the total source obtained by
[Cabral, Lobo, and Rubiera-Garcia (2019), Eqs. 20--24](https://arxiv.org/abs/1902.02222).

## Scope

TOR1 does not derive the torsion-free Dirac kinetic stress, full ECD tetrad
source, harmonic closure, regular centre, constraints, vacuum exterior,
defocusing history, or global spacetime. Those remain ECD-ID1. It also does not
derive the homogeneous `w=+1` spin-fluid closure assumed in EC-1 or fix a
late-time dark-energy value.

## Reproduction

```bash
python3 scripts/reproduce_gmf1b_ecd_torsion.py \
  --output results/gmf-1b-ecd-torsion.json
```

The exact rational fixture is audited by `scripts/check_repo.py`.
