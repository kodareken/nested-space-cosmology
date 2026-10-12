# Finite Information-Map Demonstrator

## Status and boundary

This document specifies a small, executable mathematical demonstrator for the
information-map requirement in Recursive Horizons.  It is **not** a
black-hole-to-cosmology solution, a quantum theory of gravity, a factorization
of the physical gravitational Hilbert space, or evidence that a child universe
exists.  Its value is narrower: it makes the minimum information-theoretic
requirements explicit and lets them fail under deterministic tests.

The code uses a stated semiclassical code-subspace approximation.  It assumes a
fixed effective background, a chosen relational dressing/reference structure,
and a finite number of low-energy modes with negligible backreaction.  Only in
that approximation do we write effective observable algebras

```text
A_in = B(H_in^code)
A_R  = B(H_R) tensor I_C
A_C  = I_R tensor B(H_C).
```

An actual covariant completion must identify the gravitational constraints,
edge/boundary data, dressing, slicing, and an error tolerance for this
approximation.  Gauge constraints and gravitational dressing are precisely why
this tensor product cannot be promoted to an exact claim about gravity.  See
Donnelly and Freidel, [Local subsystems in gauge theory and gravity](https://arxiv.org/abs/1601.04744).

## The finite map

Let

```text
H_in^code = span{|i> : i = 0, ..., d-1}.
```

The effective radiation/exterior sector `R` and child sector `C` each have
dimension `d + 1`.  Their last basis element, `|blank>`, is orthogonal to all
code states.  For `p` in `[0, 1]`, the demonstrator defines

```text
V_p |i> = sqrt(1-p) |i>_R |blank>_C
        + sqrt(p)   |blank>_R |i>_C.
```

The blank states make the two terms orthogonal, so

```text
V_p^dagger V_p = I_d.
```

Consequently `Phi_p(rho) = V_p rho V_p^dagger` is a completely positive,
trace-preserving channel with one Kraus operator.  This is a finite
Stinespring-isometry example, following the general framework of W. F.
Stinespring, [Positive Functions on C*-Algebras](https://doi.org/10.1090/S0002-9939-1955-0069403-4), *Proceedings of the AMS* **6**, 211–216 (1955).

It is best understood as a coherent *location/erasure encoding*: `p=0` leaves
the code in `R`, `p=1` transfers it to `C`, and an interior value encodes it in
the joint output.  It does not derive `p`, describe a real emission history,
or show that a parent transition has this form.

Tracing out one output gives exact flagged-erasure channels:

```text
rho_R = (1-p) rho + p |blank><blank|
rho_C = p rho + (1-p) |blank><blank|.
```

Here `rho` occupies the `d`-dimensional code block.  The inverse `V_p^dagger`
recovers a state only when the full joint output remains in `range(V_p)`; it is
not recovery from either marginal alone.

## No cloning and tested invariants

For a pure input `|psi>`, the fidelity of the `R` code block with `|psi>` is
`1-p`, while that of the `C` code block is `p`.  For an interior `p`, neither
marginal is a deterministic complete copy.  The test suite also records the
standard no-cloning overlap contrast: a linear isometry must preserve the
overlap of `|0>` and `|+>` (`1/sqrt(2)`), whereas a hypothetical perfect cloning
map would produce its square (`1/2`).  See Wootters and Zurek,
[A single quantum cannot be cloned](https://doi.org/10.1038/299802a0), *Nature*
**299**, 802–803 (1982).

The implementation verifies, for dimensions `d=2` and `d=3`, endpoint and
interior probabilities, pure and mixed inputs:

- `V^dagger V = I`;
- trace preservation, Hermiticity, and positive semidefiniteness of the output;
- `V^dagger (V rho V^dagger) V = rho` on the full joint code range;
- exact partial-trace erasure blocks;
- pure-state preservation of the global fine-grained entropy;
- analytic marginal-entropy formulas; and
- input validation for illegal dimensions, types, probabilities, state shapes,
  and nonfinite amplitudes or entropy inputs.

These are mathematical checks.  They do not show that an actual covariant
transition preserves information.

## Entropy ledger

All logarithms below are natural and all entropies are in nats.

| Quantity | Definition | What the toy shows | What it does not show |
| --- | --- | --- | --- |
| Fine-grained entropy | `S_vN(rho) = -Tr rho ln rho` | An isometry preserves it: `S_vN(V rho V^dagger)=S_vN(rho)`. A pure input produces a pure joint `R tensor C` output. | That a gravitational transition is unitary. |
| Marginal/code entropy | `S(R)=h_2(p)+(1-p)S_vN(rho)`, `S(C)=h_2(p)+pS_vN(rho)` | Accessible effective sectors can have entropy because of correlations/erasure. | A thermodynamic horizon entropy. |
| Generalized entropy | `S_gen[Sigma]=A(Sigma)/(4 G hbar)+S_vN^bulk+renormalization terms` | Nothing; no surface, counterterms, or gravitational state is implemented. | A bare area-plus-code entropy identity. |
| Horizon entropy | semiclassical `S_BH=A_BH/(4 l_P^2)` and the corresponding de Sitter expression | Nothing; the toy has no horizon geometry. | The logarithm of this toy's code dimension or its marginal entropy. |

The generalized-entropy comparison requires a specified surface, observable
algebra, and renormalization prescription.  For context see Engelhardt and
Wall, [Quantum Extremal Surfaces](https://arxiv.org/abs/1408.3203).  Island/Page
curve results in controlled settings, such as Penington,
[Entanglement Wedge Reconstruction and the Information Paradox](https://arxiv.org/abs/1905.08255), motivate careful code-subspace language but do not establish this parent-to-child process.

## Logical relation to H-SAT

There is no implication in either direction:

```text
V^dagger V = I  does not imply  S_BH,parent = S_dS,child.
H-SAT            does not imply  an isometry, sector factorization, or unitarity.
```

The map preserves fine-grained state information within its assumed code
subspace.  It does not set a horizon area, `Lambda`, `p`, or a parent mass.
H-SAT remains the repository's separate optional coarse-entropy boundary
postulate.  A future covariant microphysical derivation could connect the two,
but the connection must be derived rather than imported from unitarity.

The same separation applies to [DSF-1 and DST-1](domain-settings.md). The
chosen map parameter `p` is not a cosmological constant, a domain setting, a
scalar-transfer rate, or evidence of inherited physical parameters. A future
theory would have to derive both its effective settings and its information
channel from the same covariant transition; neither present demonstrator
supplies that relation.

## Explicit failure criteria

The finite demonstrator fails if any tested map has an output smaller than its
input, violates `V^dagger V=I`, changes trace, produces a non-Hermitian or
non-positive output, fails the full-joint code-range round trip, accepts
nonfinite state data, or is claimed to give two identity marginal channels on
all code states.

The scientific transition proposal still fails—even if every finite test
passes—if it cannot specify approximate dressed algebras and their error,
cannot derive this or another channel from a covariant regular background, or
compares generalized/fine-grained/coarse horizon entropies without a common
physical definition.  Those are the next physical requirements, not optional
interpretations of this demonstrator.
