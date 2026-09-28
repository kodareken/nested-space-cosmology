# The transmitting CTP action has an executable first differential

Step F applies the existing finite Gaussian action to the serialized NSC
channel covariances. It reuses `nsc_influence.influence`; no state evolution
or source integration is repeated. With the initial covariance held fixed,

$$
\Gamma_G=-i\log\det Q,\qquad
Q=I-C_0+C_0V_-^\dagger V_+,
$$

Here $V_\pm$ are full canonical branch unitaries within a retained channel,
not the rectangular parent-from-child response block also called $V$ in
`GaussianBoundaryState.reconstruct`.

the imported matrix differential gives

$$
\delta\Gamma_G=-i\,\mathrm{Tr}
\big[Q^{-1}C_0(\delta V_-^\dagger V_++V_-^\dagger\delta V_+)\big].
$$

`ctp_first_variation` evaluates this expression, checks the unitary tangent
conditions and retains the linear-solve residual. The [record](../results/development/nsc-transmitting-ctp-variation.json)
compares it with a finite difference of the existing action for all 33
retained channels. This is an application to the NSC payload, not a new
determinant identity.

## Which variation could select the kernel?

A common change of both branches keeps $V_-^\dagger V_+=I$. Its action
variation vanishes for every $V$. Minimizing this equal-history value cannot
select a member of the 447,488-dimensional admissible family.

The source uses a difference-history variation. The explicit convention is

$$
V_+(s)=e^{isX/2}V,\qquad V_-(s)=e^{-isX/2}V,
$$

which at $s=0$ gives

$$
\frac{d\Gamma_G}{ds}=\mathrm{Tr}(VC_0V^\dagger X).
$$

For the normalized identity generator $X=I/\sqrt{d_c}$ this derivative is
$\mathrm{Tr}(C_{c,0})/\sqrt{d_c}$ and is independent of $V_c$. The new
record evaluates it in every channel and checks its sign against the existing
action. The dimensionless finite-difference step is a numerical derivative
parameter; it is not a metric-history duration.

This is a bare finite Gaussian differential, before the reference and
same-action local/boundary allocation. It is not a renormalized metric source
or a non-existence condition on the full action. In particular, declaring
every unrestricted unitary direction to be an independent physical equation
would change the variational problem. Its physical directions must instead
come from the actual $U[g,B,\mathrm{embedding}]$ and the varied fields.

## Result

**Conditional differential PASS; physical selection OPEN.** Neither common
branch flatness nor this source differential selects $V_c$. The earlier
state-changing mixing diagnostic remains an authenticated input and is not
rerun or declared resolved. The missing ingredient is the transmitting
link/embedding functional and its admissible variations, together with the
remaining derivative of the same boundary action.

The code can now evaluate the Gaussian contribution once those derivatives
are supplied. [Step G](nsc-transmitting-metric-pullback.md) connects it to
the existing endpoint object, and [Step H](nsc-extended-action-completion-gate.md)
retains the corresponding OPEN decision.

```sh
python3 scripts/derive_nsc_transmitting_action_completion.py --check
```
