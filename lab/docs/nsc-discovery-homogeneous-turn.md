# Homogeneous source-selected areal turn

This consumer restricts the unchanged e6 leading four-geometry action to the
homogeneous symmetry sector authenticated in PREPv2. It does not introduce a
new force, vacuum interpretation, or solver framework. The three occupied AP
sectors have positive momenta `pi/8, 3pi/8, 5pi/8`; their negative-momentum
partners are `sigma1 u_j`, with the same occupation. Momentum parity is not an
antimatter label. Six standing columns are recovered through fixed recorded
frame coefficients, rather than replacing the source by newly chosen columns.
The actual inactive-modal tail and cross-momentum covariance are reported;
exact closure refers to the symmetry idealization, with recorded roundoff in
the retained source. The coherent source is outside this reduction.

Set `a=-8pi A`, `Z=3a`, `mag=2pi C_F flux²`, carrier period `ell=8` and
`M=4`. The lapse is `Q`; `ell` is not a lapse. The momenta `p=p_Q,v=p_r` are
nodal densities decoded from canonical momenta. The 16 real state components
are `(Q,r,p,v,Re u_1..u_3,Im u_1..u_3)`, and
`H_j=p_j sigma2+kappa Q sigma1`, `i udot_j=H_j u_j`.
Multiplicity applies once outside the one-particle Hamiltonian:
`E=2M sum c_j u_j† H_j u_j`, `S=2M kappa sum c_j u_j† sigma1 u_j`.

The unchanged homogeneous geometry equations are

```text
Qdot = Q v/(2ar) - 3Q² p/(2ar²)
rdot = Q p/(2ar)
pdot = -p v/(2ar) + 3Q p²/(2ar²) + 2Q(-ar²-mag) - S/ell
vdot = Q p v/(2ar²) - 3Q² p²/(2ar³) - 2ar Q².
```

The constrained total energy is
`ell[Qpv/(2ar)-3Q²p²/(4ar²)-Q²(-ar²-mag)]+E=0`.
The field work identity is `Edot=S Qdot`; geometry receives its opposite.
An auxiliary work integral checks this identity without changing the physical
16-state equations. At the downward `p_Q=0` event,

```text
r_turn² = mag/(8pi A) + E_turn/(8pi A ell Q_turn²)
rdd_turn = (2E_turn-Q_turn S_turn)/(2a r_turn ell).
```

Thus a submagnetic radius requires negative field energy in this declared
finite source convention. The radius is selected by the actual source
history. The saved T3 full-band geometry and energy provide a comparison;
the saved minimum sample at T2.6 is not an independently rooted event.
Finite-CAR state activity, physical filled-sea Gamma matching and the
strong-curvature EFT domain remain explicit limits. No universal bounce,
regeneration, holding, or antimatter claim follows.

The fixed source perturbation is `alpha=.002`, direction `(1,0,-1)` per
momentum pair: held occupations `.752,.5,.248`, repeated twice, trace three.
The initial radius is recomputed from the same constraint, with fixed `Q0`
and zero initial momenta. Its analytic tangent includes `delta r0`; the
spinors are unchanged. The baseline full tangent predicts
`rstar(alpha)=rstar(0)+alpha delta r(tstar)`, since the first event-time term
in radius vanishes at `rdot=0`. The time derivative is
`delta tstar=-delta p_Q/pdot_Q`. Centered event solves at `.0001,.00005`
check this derivative; the held `.002` trajectory is never a forecast input.

Default CLI operation only authenticates and previews. After freezing both
producer modules and tests in one commit, run these stages in order:

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_homogeneous_turn.py --prepare
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_homogeneous_turn.py --predict
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_homogeneous_turn.py --measure
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_homogeneous_turn.py --check
```

`--output` may select a new child of
`lab/results/development/nsc-discovery-homogeneous-turn-v1`. JSON/NPZ records
are creation-only, hash-bound and limited to 64 MiB. Prepare and predict seal
both producer closures before the independent NF32 full-field backend reads
the lock. That backend independently folds the original AP source, derives
its own constraint radius and locates two held events at matched step caps;
it does not use the predicted radius/time to select its initial state or
root. The aggregate production allowance is 30 CPU seconds. Read-only checks
authenticate original producer bytes through their recorded Git commit,
without rewriting historical hashes. These comparisons are numerical
indicators rather than continuous error certificates.
