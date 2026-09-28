# Same-surface incoming constraint assembly

The lapse and shift coefficients of the declared action are assembled as

$$
\mathcal E_B[j]=\mathcal E_B^{\rm matter}[j_0,C_0]
+\mathcal E_B^{\rm local}[j]
+\sum_\alpha\frac{d_\alpha}{2\pi}\int_{\mathbb R}dk\,
\operatorname{tr}\{(P_\alpha[j]-P_\alpha[j_0])V_B\},
\qquad B=N,\beta.
$$

Here `j` belongs to the existing incoming normal-jet domain at rho=1. The
intrinsic metric, canonical frame, Hamiltonian vertices and physical C0 are
unchanged. The change of the positive subtraction action is added with the
displayed sign. Each angular sign receives its existing share of the group
multiplicity. The ordinary insertion is the bulk coefficient under the
complete unweighted phase-space trace and compact variations, using the
[paired reference-band exchange](nsc-reference-band-bulk.md). Finite-cutoff
and physical endpoint forces do not follow from this bulk identity.

The baseline matter term reuses the authenticated finite non-LLL panels,
replaces the group13 subgap contribution, adds the existing approximate
paired tail and adds physical LLL state once. Group13's exact closed-channel
current is zero; its raw floating zero residual remains separately recorded.
The local action already contains all light geometry and the locked compact
terms. Adding the finite source's entire light allocation would count the
light geometry twice and is rejected by an independent allocation check.

## Current decision and numerical domain

Reuse the spatial order-four projector, local Euler owner, finite panels,
group13 refinement and paired source tails. The missing connection is the
integrated reference change when normal data vary. The smallest experiment
uses one labeled `delta(dT dz r)=0.01` control with unchanged other jets.
This is a response probe, not a physical IV choice. It determines whether
the complete reference/local response needed by the joint solve is now
available and numerically resolved; it does not select a constraint root.

Both signs of k are integrated using a Gauss rule under
`k=a sqrt(m^2+(lambda/r)^2) s tan(theta)`. Numerical map scales and node
counts change while the full momentum domain stays fixed. Leading P0 is
exactly unchanged and checked before subtracting orders one through four.
The LLL uses its existing fixed chiral charts, with zero reference change;
its c=4 geometric allocation remains in the local owner. A content-addressed
artifact retains every new nodal insertion and quadrature weight. Replay
contracts those arrays and authenticates inputs without repeating reference
or local scientific generation.

Stop once reference-response refinement and the assembly identities meet
`3e-11`, or report the concrete numerical obstruction. Risks are cancellation
at large momentum, finite endpoints outside the bulk identity, and existing
physical mode errors not bounded by quadrature changes. The latter remains
OPEN even if this assembly passes. Baseline base/refined panel differences
already exceed the separate stationarity tolerance after action projection.

The [record](../results/development/nsc-incoming-joint-constraints.json)
reports the two raw bulk action coefficients per dT dz and their opposite
forces. Its `certified_constraint_residual` and `full_source_error_bound`
remain null. No optimizer, physical history, endpoint completion, new term,
scale refit, or metric timestep is used. The next work is physical source
error closure and constraint-compatible normal data, followed by extended
stationarity; the fixed compact-pulse exclusion remains in force.

```sh
python3 scripts/derive_nsc_incoming_joint_constraints.py --prepare
python3 scripts/derive_nsc_incoming_joint_constraints.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_reference_response.py
```
