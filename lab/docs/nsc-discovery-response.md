# Stage-5 first wave: full-state geometry response

This note specifies the retarded response of the existing conformal
Galerkin / nested parent–child system and records the directional
primitive that this wave actually evaluates. It does not publish a
result and it does not evolve a history.

The executable owner is
[nsc_discovery_response.py](../src/recursive_horizons/nsc_discovery_response.py).
Checks live in
[test_nsc_discovery_response.py](../tests/test_nsc_discovery_response.py).

## What is differentiated

The rate is the existing map in
`nsc_spherical_galerkin_coupling.compose_fine_hamiltonian`, then the
constant NestedPair encode

\[
a=W^{T}g,\qquad \pi=\Delta x_{g}W^{T}p .
\]

A tangent may move every geometry coefficient \((Q,r,\chi)\), every
momentum, both spinor columns, and the six source occupations. The
representative column rate remains \(-iH\Phi\). The factor \(M=4\kappa\)
is already inside the nodal forces and is not applied again.

In the conformal chart \(L=Q\) and \(\beta=0\), and \(\dot p_Q\) receives
both \(F_Q/\Delta x_q\) and \(F_L/\Delta x_q\). The primitive reuses the
quadrature derivative and the antiperiodic momentum stored on the fine
system. It does not build another Fourier matrix.

The Jacobian-vector product is analytic. A centred difference is only a
check. If the owned projected radius Jacobian were singular, the radius
tangent would be returned as a named missing primitive rather than
replaced by a difference.

## Initial source tangent

The initial radius remains the root of the owned projected Hamilton
residual at the constant ratio \(Q=b_0/a_0\), with \(\chi\) and the
momenta held at zero. For a source tangent \((\delta\Phi,\delta c)\) the
preparation solves one linear system

\[
J\,\delta r=-\mathrm{pull}(\delta\rho),\qquad \rho=F_L/\Delta x_q,
\]

with \(J\) the existing `_projected_radius_jacobian`. The base radius is
not updated. No Newton loop is run.

## Clock and readouts

The child clock rate is the owned periodic sample

\[
\dot\tau=(rQ)(x=2).
\]

In the conformal chart this is \(rL\). On a declared coordinate interval
of length \(\Delta t\), with the rate held at the slice,

\[
\delta\tau=\delta\dot\tau\,\Delta t.
\]

That product is not an integral along an evolved trajectory. Matched
proper time uses

\[
\delta O\big|_{\tau}=\delta O\big|_{t}-\dot O\,\frac{\delta\tau}{\dot\tau}.
\]

The primary readout is the child regional content: the integral over
\(I_C=(1,3)\) of the periodic interpolant of the one-body density of the
full column Gaussian. The secondary readout is the child proper mean of
\(r\), with weight \(rQ\), the same quotient as the nested metrics.

## Gaussian state

The quasifree many-body state is the full one-particle covariance

\[
C=\Phi\,\mathrm{diag}(c)\,\Phi^{\dagger}
\]

on the column space, including the unoccupied complement. Complement
occupation 0 is the support of this finite \(C\). It is not a vacuum
identification. Connected force and noise use the existing influence
identity

\[
\operatorname{Tr}\big(C A(I-C)B\big)
\]

and are not rewritten. For one orthonormal column of occupation \(n\in(0,1)\)
the diagonal value is \(n(1-n)\), which is nonzero. The induced spectral
force remains `SpectralInducedSource` inside `CausalCommonFunctional`.
This wave does not add a second copy.

## Future effective stress

The CTP future effective stress is specified as the sum of three
contributions:

1. the complementary-branch force, matter plus the single existing induced coefficient;
2. the initial cross, the child-observer block of \(C\) against the rest of the one-particle space;
3. the retarded change of the full state.

Only the initial cross is evaluated, and only at the current slice. The
complementary future branch and the retarded state transport are not
implemented. Their slots stay empty. The sum is not returned as a number,
and a finite difference is not written into the empty slots. Equal-time
\(F_L\) and \(F_Q\) remain inside the rate; they are not a substitute for
the future complementary force.

## The 2×2 kernel

The child observer kernel is the \(2\times 2\) block

\[
K_C=V_C^{\dagger}CV_C
\]

on the fixed middle pair, indices \((2,3)\). It is a modal covariance.
It is not the spatial exterior of \(I_C\), and a nonzero entry is not
spatial leakage. Spatial content is the collar integral above.

## Frozen checks and blockers

Checks:

- centred two-sided differences of selected full-state directions;
- one small manufactured positive-chart state;
- the saved v5 initial arrays `nf256`, read in the conformal chart, with no time step;
- the connected bilinear on the full column Gaussian;
- one linear solve with the owned projected radius Jacobian;
- the matched-\(\tau\) formula on a manufactured clock.

Blockers, left explicit:

- `ctp_future_effective_stress`;
- `retarded_delta_state`;
- no evolution campaign in this wave;
- the induced coefficient is not rederived;
- the \(2\times 2\) kernel is not spatial outside.

## Reproduction

```sh
python3 scripts/lab.py -m pytest tests/test_nsc_discovery_response.py -q
```

The command reads the saved initial state and compares derivatives. It
does not write a record, evolve a trajectory, or touch a manuscript.
