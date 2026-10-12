# Nonlinear transmitting mode histories on the common PG slice

The mode propagator now evaluates a supplied metric history in the existing
PG field representation. Its numerical unknown is a **spatial field**, not
a closed matrix made from the few frequencies used to inspect its kernel.
The metric is an input to this conditional propagator; no geometry is evolved
by an Einstein equation in this calculation.

## Whole-field evolution relative to the known preparation

The [prepared modes](nsc-pg-retained-covariance.md) give $\Phi_E(\rho)$ and
the unchanged horizon/incoming source covariance. Use the already owned
canonical Dirac operator $H[g]$ on the same transmitting domain and define

$$
w_E(\tau)=(U[g](\tau,\tau_i)-U_0(\tau,\tau_i))\Phi_E,
\qquad w_E(\tau_i)=0.
$$

Applying the imported variation-of-constants identity to that operator gives

$$
i\partial_\tau w_E
=H[g](\tau)w_E+
(H[g](\tau)-H_0)e^{-iE(\tau-\tau_i)}\Phi_E.
$$

The metric dependence is evaluated nonlinearly through the existing KS–PG
map and canonical lapse, shift, radial-scale and radius coefficients. No
instantaneous $B$ or another matter interaction replaces the Dirac operator.
The updated mode field is $e^{-iE(\tau-\tau_i)}\Phi_E+w_E$.

The actual whole-field covariance would follow by integrating these updated
fields with the same $C_H\oplus n_{\rm in}$ and the complete spectral measure.
The present computation does **not** substitute a sampled frequency matrix
for that integration or assign a new absolute stress.

## Computational domain and time prescription

The perturbation uses the existing compact spatial collar. Its finite
propagation domain lies inside the trapped chart. Both principal velocities
$\pm N/q-\beta$ point toward decreasing $\rho$. The **difference** field has
zero initial data and zero incoming difference at the right computational
edge. The left edge is outflow, with its norm flux retained. Neither edge is
a reflecting physical boundary or a new compact carrier.

The spatial split/SAT convention is imported from `nsc_lorentzian`; its
norm matrix and flux identity are checked for the new metric coefficients.
The coefficient history is sampled at time midpoints. A sparse exponential
integrates each forced step. Auxiliary exponential amplitudes inside this
variation-of-constants algorithm are deterministic solver variables, not new
particles, a reservoir or terms in the action.

The response control uses a smooth time envelope with all derivatives zero
at both ends. Consequently the initial metric and its time jets agree with
the reference on which the state was prepared. The interval and amplitudes
are labeled numerical controls, not a selected physical duration or a
source-generated trajectory. The grid buffer and traces of the difference
field are reported; a depleted sampled causal buffer is rejected.

## What the checks establish

The nonlinear response is compared under separate spatial and temporal
refinements. Its small-amplitude limit is compared with the already verified
[spacetime Duhamel derivative](nsc-ks-spacetime-variation.md). This distinguishes
the new conditional propagation from merely applying the old linear kernel.
Horizon, angular-spectrum and seed generators are not replayed. The same
radial mode equation is sampled on the newly required propagation grid from
the authenticated horizon/infinity columns.

All channel labels, scales and initial source occupations retain their
recorded meanings. Real-energy mode samples inspect this propagator. Complex
subgap contour data are not relabeled as canonical frequency states; their
analytic covariance representation must be retained when propagating that
part of the complete state.

The full physical endpoint/history binding, complete evolved covariance
integral and renormalized source still have to be composed. In particular,
the finite projected response is not itself a unitary branch matrix for
`EndpointBranchJets`: it omits the unsampled source energies. The full
Gaussian action includes that complement once. No additional
$\Gamma_{\rm rest}$ is introduced, and the old Weyl time-node diagnostic is
not treated as an interface force. Extended stationarity remains OPEN and
metric stepping remains closed.
