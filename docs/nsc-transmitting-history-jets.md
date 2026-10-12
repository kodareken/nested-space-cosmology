# Nonlinear metric jets of the transmitting Gaussian history

The metric derivative now acts on the evolved field of a supplied nonzero
history. The Jacobian is evaluated at the actual metric, and the Gaussian
contraction retains the same initial source covariance. No new functional is
introduced and no sampled response matrix is exponentiated as a closed bulk.

## Differentiate the same Dirac history

In the canonical PG characteristic frame, write $\partial_\tau\chi=L[g]\chi$,
including the already owned transmitting exterior data. For a raw KS metric
direction, $Y_B=\delta_B\chi_g$ satisfies

$$
\partial_\tau Y_B=L[g]Y_B+\delta_B L[g]\chi_g,
\qquad Y_B(\tau_i)=0.
$$

The control directions vanish with all time jets at the initial surface,
so the prepared state and its source phases are fixed. The tangent and state
are propagated together, rather than driving the tangent with the unevolved
reference field. Lapse, shift, radial scale and radius are independent
directions; their KS-to-PG Jacobian is evaluated at the current metric.

The incoming traces are those of the same global reference modes. A direct
spatial evolution is used for both reference and varied fields, so its finite
current identity and derivative use one discretization. Drift of the discrete
reference from the analytic reference is reported as a numerical diagnostic.
It is not subtracted as a stress or added as a compensating force.

## Couple the field jets to the existing CTP differential

Form the source-mode kernel **after** evolving the spatial fields:

$$
\mathcal A_B(E_o,E_i)
=\int d\rho\,\chi_{g,E_o}^{\dagger}(\tau_f,\rho)
Y_{B,E_i}(\tau_f,\rho)
=\bigl(U[g]^\dagger\delta_B U[g]\bigr)(E_o,E_i).
$$

This includes the propagation through unsampled intermediate energies.
At the physical equal-history point, move both branches to the current-history
frame, $W_\pm=U[g]^\dagger U_\pm$. Then

$$
W_\pm=I,\qquad \delta W_\pm=\pm\tfrac12\mathcal A_B,
\qquad
\delta_B\Gamma_G=-i\,\mathrm{Tr}(C_0\mathcal A_B).
$$

The existing `ctp_first_variation` contracts this first derivative on the
inherited energy quadrature. The weights multiply the kernel as
$\sqrt{w_ow_i}/(2\pi)$; $C_0(E)$ remains the source-fiber multiplication
operator. Its two horizon components and incoming occupation are retained.
Using $I$ in this **relative first variation** is not a claim that the sampled
frequencies provide a closed finite representation of $U[g]$.

The computed source values are finite spectral-sample controls of the bare
Gaussian action derivative, per reduced radial family. They are not the
complete energy integral, a renormalized four-dimensional stress, or a
stationarity residual of the full action.

## Numerical owner and verification

The second-order spatial scheme showed a source-variation refinement error
above the fixed tolerance. The new jet owner therefore imports the standard
diagonal-norm SBP(4,2) coefficient table, as tabulated in
[SummationByPartsOperators.jl](https://github.com/ranocha/SummationByPartsOperators.jl/blob/21414744283d75d6172cd1763dd5becff5988909/src/SBP_coefficients/MattssonNordstr%C3%B6m2004.jl#L36-L74).
The continuum Dirac operator, state and coefficients are unchanged. Earlier
second-order records retain their bytes and scope.

The focused checks compare analytic jets with independent variations of the
same nonzero metric history, and assess spatial and temporal refinement of
the actual Gaussian contraction. Norm/flux tangents, the anti-Hermitian overlap
identity and the Gaussian trace identity are separate algebraic checks. A small
algebraic residual is not substituted for the refinement test.

The old eight diagnostic KS end-node coordinates are not assigned a new
physical meaning. These jets differentiate the declared response directions
on the fixed PG Cauchy domain, with raw KS metric components. A physical
endpoint variation family and the complete spectral/reference allocation
must still be bound before declaring full B1 or evaluating B2. Smooth-seam
geometric cancellation remains imported; the $93.54264532195464$ time-node
diagnostic is not a boundary-force target. No selected metric, new absolute
stress, new null contractions or coupled metric evolution is produced here.
