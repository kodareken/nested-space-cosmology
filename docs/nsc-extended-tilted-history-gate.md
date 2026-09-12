# The three tilted-history components compose, but stationarity remains open

The general-KS reference, node-wise local history and transmitting tilted
interface now have separate executable records.  Their composition gives a
single residual vector without selecting either diagnostic interface witness.

The fourth-order reference passes its Bloch, normalization, parity-momentum
and four-force gates.  The local bulk action passes its node-gradient and
Euler/box boundary-support controls.  Its free Weyl endpoint variation remains

$$
\max|\delta S_W|_{\partial}=93.54264532195464
$$

against a $3\times10^{-11}$ stationarity tolerance.  The tilted interface
family passes weighted current conservation, unitarity, CAR and full-rank
checks, but its weighted-isometry constraint leaves an admissible family of
dimension $447488$.  The current boundary data do not select one $V_c$.

Therefore the extended action derivative is not yet a single-valued function
of the history.  Running an optimizer would choose an interface kernel or
endpoint completion implicitly and would not test the declared action.  The
extended gate consequently claims neither existence nor non-existence, and no
finite stress, null sign or updated constraint is emitted.

The two open residuals have one common, already declared owner:
`ModeResolvedTransmittingBoundaryHistoryAction`.  It is the existing
$\Gamma_{\rm rest}$ in

$$
\pi_p^{ab}+\pi_c^{ab}
+\frac{\delta\Gamma_{\rm rest}}{\delta h_{ab}}=0.
$$

It must derive the mode-resolved scattering/embedding kernel $V_c$ from the
link/boundary dynamics and supply the same interface metric variation that
completes the Weyl endpoint momentum.  This is an implementation of an existing
term, not a new counterflow.

The killed homogeneous no-interface non-existence result remains a passing
regression.  Coupled evolution remains closed.

Reproduce the composition with

```sh
python3 scripts/derive_nsc_extended_tilted_history_gate.py --check
```
