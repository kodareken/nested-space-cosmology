# Binding the NSC charged sector to the imported fermionic throat

Maldacena, Milekhin and Popov already establish the central existence result
needed here: charged massless-fermion Casimir energy can support a
semiclassical Einstein–Maxwell throat. This calculation does not repeat their
gravity solution. It asks whether the current NSC operator, domain, state and
coefficients restrict to that established system.

## Exact coefficient dictionary

Write the retained Lorentzian local action as

\[
S=\int\sqrt{-g}\,(A R-CF_{\mu\nu}F^{\mu\nu}-V).
\]

The MMP normalization is

\[
G_N=\frac1{16\pi A},\qquad
g^2=\frac1{4C},\qquad
r_e^2=\frac{q^2C}{4A},\qquad
\ell=\frac{16r_e^3}{|q|G_N}.
\]

The first three identities follow exactly by matching the Einstein and
Maxwell actions. The fourth is the imported self-sourced throat relation.
Magnetic flux \(q\) remains an integer sector rather than a continuous fit.

## Field and state map

The existing opposite-parity compact domain supplies one anomaly-free
four-dimensional Dirac zero field and two massive Dirac fields at every
nonzero compact level. Magnetic flux gives \(|q|\) complex lowest-Landau
two-dimensional channels. Two opposite parity assignments remain admissible.

The existing closed-cell potential has its stationary minimum at effective
antiperiodic phase \(\alpha=1/2\) modulo one, with positive phase Hessian.
This supplies the state-selection match on a closed return path. MMP's magnetic
field-line return path is the imported reference domain. It is not identified
with the unwrapped NSC radial line.

## Result of the first embedding gate

At the existing matching point \(\nu=1\), the retained partial Dirac
coefficients give

\[
\frac{\Xi_{\rm partial}}{q^2}
=\frac{V_D C_D}{8A_D^2}
=3.4459938566757384.
\]

The charged RN–de Sitter seed bound recorded by the existing coefficient
calculation is \(\Xi/q^2\leq1/4\). The retained partial contribution therefore
does not embed that seed. More importantly, the repository has not fixed the
complete common-action vacuum coefficient \(V_{\rm full}\), and the physical
closed recursive return state is also not derived. The full embedding cannot
be decided from the partial coefficients.

This identifies the first unresolved owner without rerunning a throat:

\[
\boxed{V_{\rm full}\ \text{from the same normalized common functional}.}
\]

No value is assigned to it. The recursive geometry solve remains downstream
of this gate.

## Reproduction

The compact record authenticates the charged-sector, compact-matching,
compact-Casimir and vacuum-charge inputs. It evaluates only the dictionary and
gate:

```sh
python3 scripts/derive_nsc_mmp_embedding.py --check
```

The MMP Einstein equations, compact spectrum, Casimir potential and historical
NSC generators are not rerun.
