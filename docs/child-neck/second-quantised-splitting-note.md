# Does the torsion splitting survive second quantisation? The vector-axial term

7 Oct 2026. Script `second_quantised_splitting.py` (sympy mean field plus numbers), data `second_quantised_splitting.json`.
Sources read: Freidel, Minic, Takeuchi, "Quantum gravity, torsion, parity violation and all that", hep-th/0507253 (Eqs. 21-23, 37-48).
Follows `torsion-splitting-and-inherited-spin-note.md`. Working assumption: the splitting is needed; this checks which term can supply it with field operators.

## Effective interaction

FMT Eq. 23: `S_int = -(3/2) pi G (gamma^2/(gamma^2+1)) Int [ A.A + (2 alpha/gamma) A.V - alpha^2 V.V ]`,
`V^mu = psibar gamma^mu psi`, `A^mu = psibar gamma^mu gamma5 psi`, `alpha` the non-minimal coupling (zero for minimal coupling), `gamma` the Immirzi parameter.
Only the `A.A` term survives at `alpha = 0`; it is the paper-of-record contact term used for the bounce. The vector-axial term is parity violating.

For a gas at rest the energy density is `H = kappa' [ h^2 - s^2 + (2 alpha/gamma) h n_V - alpha^2 n_V^2 ]`, `kappa' = (3 kappa/16) gamma^2/(gamma^2+1)`:
`h = A^0 = N_R - N_L` is the net helicity density of all quanta (particles and antiparticles), `n_V = V^0` the net number, `s` the spin density.

## Result of the mean-field check (sympy)

Single-quantum shifts for the four types (helicity plus or minus, particle or antiparticle):

- `A.A` term: splitting between a particle and an antiparticle of the same helicity is exactly `0`. The spin-spin contact term is C even and shifts both equally.
- Vector-axial term: splitting `4 kappa' (alpha/gamma) h`. It needs a net helicity density `h` (a P-odd background), not a spin polarisation.
- `V.V` term: splitting `-4 kappa' alpha^2 n_V`. It needs net number already present and lowers the matter levels relative to antimatter (an amplifier).

So the particle-antiparticle splitting exists in second quantisation if and only if `alpha != 0`, with a source `h` or `n_V`.

## Strength at the bounce

Level shift (fermion `+delta`, antifermion `-delta`): `delta = 2 kappa' (alpha/gamma) h = (3 kappa/8)(beta/2) h`, with `beta = 2 gamma alpha/(gamma^2+1)` (FMT's parity-violating parameter).
This has the same form as the classical Hehl-Datta shift with the polarisation `p` replaced by `(beta/2) chi_h`, `chi_h = h/n`. With the Fermi-Dirac result of the previous note, `eta = 34.5 (beta/2) chi_h`.

| target `eta` | 6e-10 | 9e-8 |
|---|---|---|
| `chi_h` needed for `beta = 1` | `3.5e-11` | `5.2e-9` |

Experimental bound on `beta` (FMT, atomic parity violation): `beta < ~1e30`, practically none.

## Where `h` can come from

Helicity (chirality) is a number-type quantity. If it is inherited from the parent with no regeneration, `chi_h <= N0 / N_f = M*/M` (same bound as for `eta` earlier):

| parent mass | 3 solar | 10 solar | 30 solar |
|---|---|---|---|
| `chi_h` bound | `3.0e-20` | `8.9e-21` | `3.0e-21` |
| `beta` needed for `eta = 6e-10` | `1.2e9` | `3.9e9` | `1.2e10` |
| `beta` needed for `eta = 9e-8` | `1.7e11` | `5.8e11` | `1.7e12` |

For `beta = 1` the parent would weigh `5.1e21 kg` (`eta = 6e-10`) or `3.4e19 kg` (`eta = 9e-8`).
By contrast the classical route of the previous note takes its polarisation from angular momentum, which scales with the area (`p_max = 1.26 a*`, any mass), so it had ten orders to spare.
The two routes differ in what they need: spin angular momentum (plentiful) versus a P-odd helicity density (number bounded).

## What is not accounted for

- The source of a net helicity density in the child. Candidates, none computed: weak interactions (parity violating, left-handed couplings) acting on the produced fermions out of equilibrium; a pseudoscalar built from the geometry (the parent's angular momentum, an axial vector, dotted into the stretch direction, a polar vector).
- The sign: a helicity density flips sign under parity, and the splitting is C odd, so the combination is CP even while the helicity background is CP odd. Turning a space inside out is an orientation reversal; whether the two sides of a throat then carry opposite `h` is a question for the throat geometry, not answered here.
- The amplifier `V.V` over `A.A` for spin-averaged gas is `6 (gamma^2/(gamma^2+1)) alpha^2 eta^2`, about `1e-18` at `eta = 1e-9`; it matters only for `eta` of order `0.4/alpha`.

## Limits

Mean field at the bounce, thermal equilibrium, `h` held fixed, the effective action of FMT taken as given. Massless fermions (helicity equals chirality); masses and Yukawa couplings mix chirality and would regenerate or erase `h`. Species sum not done; one Dirac species.
