# Fourth-order time evolution of the same prepared state

The measured midpoint time change in the nonzero local response is about
1.7e-7, above its1e-11 indicator target. This activates the approved change
of numerical time method. The Dirac operator, source, geometry, restriction,
subtraction and retarded tangent equation are unchanged.

Reuse [Blanes--Moan, equation43](https://personales.upv.es/~serblaza/2006APNUM.pdf)
on the same augmented source/field/tangent matrix M(t). For each step h,

$$
c_{1,2}=\tfrac12\mp\tfrac{\sqrt3}{6},\qquad
a_{1,2}=\frac{3\mp2\sqrt3}{12},
$$
$$
X_{n+1}=e^{h(a_1M_1+a_2M_2)}e^{h(a_2M_1+a_1M_2)}X_n,
\qquad M_j=M(t_n+c_jh).
$$

Both exponentials include the harmonic source block. Since a1+a2=1/2,
their total source phase is the same exact exp(-i E h). The tangent blocks
carry actual dL at the two Gauss nodes, so the derivative is that of the
complete numerical map. F_z and dF_z still come from the actual node PDE.
Fourth order is a method property, not a continuum or full-source bound.

`prepare_cf4_incoming` returns `CF4Incoming`, a specialized exact-phase
prepared-state result. It retains the earlier state, axial derivative,
weight and preparation interfaces, and additionally binds both evaluated
Gauss-stage metrics/directions. The immutable midpoint owner and its
records are not changed. The subclass explicitly declares CF4 and validates
its stage bindings rather than pretending to have used midpoint evolution.

Three independent synthetic controls compare the constant case to a dense
matrix exponential, compare the nonzero full retarded/PDE tangent to
centered changes, and measure convergence against adaptive DOP853 on the
same small finite operator. None uses authenticated physical source runs.
Earlier-past support remains the caller's retarded-family contract; the
initial-slice check alone is not a theorem about the entire continuum past.

The physical local gate, full spectrum, continuum and between-node errors
remain OPEN. No metric timestep or source/coupling refit occurs.

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_cf4_prepared_state.py
```
