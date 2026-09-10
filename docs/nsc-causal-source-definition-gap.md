# The missing definition for a causal initial-data calculation

The next initial-data solve is blocked by the full geometric/UV source
definition. The repository has a calculated Euclidean measure, a real
geometric phase for complete periodic histories, canonical state evolution,
and the required geometric vertices. It does not yet define the full causal
geometric response and its state on the actual parent–child domain.

This is a limitation of the present realization, not a claim that nested
spaces or every possible nonlocal theory are inconsistent. A solver cannot
supply the missing dynamics by choosing its desired residual.

## What the current owners actually supply

- [VacuumMatchedDiracAction](../src/recursive_horizons/nsc_vacuum_matched_ctp.py)
  supplies a full-frequency quadratic radius vertex and a phase for an
  entire finite-Fourier history. Its continued_force_kernel is explicitly
  a bare geometric vertex, not a retarded correlator.
- [The ADM owner](../src/recursive_horizons/nsc_adm_source.py) supplies the
  covariant operator and metric vertices. It does not evaluate the complete
  nonlocal lapse and shift sources on a physical history.
- The canonical Dirac state, the local warp contribution and the
  [parent flux condition](nsc-parent-backreaction-gate.md) remain valid in
  their recorded scopes. They do not determine the missing source at an
  initial slice.

The initial constraints need F_N and F_beta from the complete action,
with an admissible geometric state/history prescription. Knowing their
Ward identities is different from knowing their values.

## Why the real branch phase does not complete that definition

To quadratic order, write the real geometric phase as

\[
\mathcal B[J]=\frac12\int dt\,ds\,J(t)K(t,s)J(s),
\qquad K(t,s)=K(s,t).
\]

With J_Delta=J_plus-J_minus and J_c=(J_plus+J_minus)/2,

\[
\mathcal B[J_+]-\mathcal B[J_-]
=\int dt\,ds\,J_\Delta(t)K(t,s)J_c(s).
\]

Its mean source depends on the complete history through K. If a symmetric
kernel were also retarded on ordinary distributional histories, symmetry
would force its support onto t=s. A translation-invariant tempered kernel
with that support is a finite combination of delta-function derivatives,
and its frequency dependence is polynomial.

The implemented function is not such a kernel. This follows directly from
its already derived finite proper-time expression; no frequency scan is
needed. Let a=Lambda^-2, E=max_i |e_i| and choose a nonzero vertex element
w=|V_ij|^2. For omega^2>4 E^2, restrict the positive part of its integral to
s in [a/2,a] and alpha in [1/3,2/3]. The continued force vertex satisfies

\[
\mathcal B_2(-\omega^2)\ \ge\
\frac{w\,a^{3/2}}{24\sqrt{2\pi}}
(\omega^2-4E^2)\,
e^{-aE^2+a\omega^2/9}.
\]

The contact term is nonnegative, and all the other bubble contributions
are nonnegative in this range. The bound establishes super-polynomial
growth. The finite-Fourier calculation remains well defined, but it is not
a definition on arbitrary smooth compactly supported histories and cannot
simply be used as a retarded initial-value source.

This argument does not require a bare geometric action to be an independent
bath influence functional. It identifies what the phase-only product has
not established: causal metric evolution and admissible initial data.
A constraint reduction might remove a nonphysical channel, or a complete
microscopic sector might supply additional contour and state terms. Neither
has been constructed for the full source.

## The missing contour data are not fixed by a one-branch determinant

The underdetermination can be stated directly. For a real antisymmetric
kernel A(t,s)=-A(s,t), a mixed contribution

\[
\Delta\Gamma=\int dt\,ds\,J_\Delta(t)A(t,s)J_c(s)
\]

vanishes for a single nonzero branch, because J A J=0. It also vanishes on
identical histories, but it changes the mixed response on the physical
branch. Thus one-branch action values and equal-history normalization do
not determine the real-time response. Quantum positivity additionally
requires the corresponding fluctuation/state construction; a guessed
mixed kernel or noise term is not a solution.

A full microscopic action and state could fix these data. The recursive
organizing identity alone does not specify them.

## What would unblock the calculation

The required input is a complete causal metric/UV realization, or a
specified whole-history quantum boundary formulation with its physical
state and observable rule. It must determine:

1. the full metric source and mixed real-time response, including the
   lapse and shift constraints;
2. its initial state/memory or admissible boundary data;
3. the consistent covariant normalization and relation to the retained
   Dirac sector and cutoff profile.

The existing finite Hamiltonian and Dirac covariance do not specify this
additional geometric/UV dynamics. Adopting a different microscopic
realization would require an explicit new physical definition and changed
claims. A local truncation, an invented bath, or a prescribed mass history
would not complete the requested theory.

The previous three goal checkpoints already left this same condition
unresolved: the ADM source construction supplied vertices; the local warp
calculation supplied a summand with an unbounded remainder; and the flux
gate supplied a necessary evolution condition. This inspection confirms
that the missing full causal source still has no owner that an initial-data
solver can call.

No old scientific generator or new stress-test campaign was run for this
decision. The completed results remain available. The unified physical
goal is not achieved, and further equivalent rearrangements or static
checks do not remove this blocking definition.
