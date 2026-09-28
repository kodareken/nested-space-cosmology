# Contracting the same CTP differential into KS endpoint coordinates

`EndpointBranchJets` extends the [existing endpoint object](nsc-weyl-endpoint-match.md)
without changing its eight Weyl coefficients, basis identifier or exact
extraction residual. For each retained channel, the required arrays are

$$
J^\pm_{A,e}=\frac{\partial V_\pm}{\partial g^A_{\Delta,e}},
\qquad A\in\{N,\beta,q_{\rm ADM},r\},\quad e\in\{0,L\}.
$$

Their shape is `(4, 2, d_c, d_c)`. At the physical equal-history limit,
`endpoint_ctp_pullback` contracts these arrays with the first differential
of the existing Gaussian action and returns its contribution in exactly the
old `EndpointVariation` basis. Numerical tests exercise the chain rule and
endpoint signs using explicitly algebraic test inputs. Those test arrays
are not assigned to the NSC geometry.

The complete matching still requires the same-action boundary remainder and
the geometric mapping of its induced metric and normal derivatives:

$$
b_{A,e}=\big[\mathrm{pullback}
(\delta\Gamma_{\rm rest}/\delta h)\big]_{A,e},
\qquad R_{\partial,A,e}=w_{A,e}+b_{A,e}.
$$

This displayed formula also corrects the missing backslash before `left` in
the authenticated v0.19.0 endpoint note. That note and its certificates retain
their original bytes; this is the current rendering of the same expression.

Renaming `q` to `q_ADM` in the causal interface is a known coordinate alias.
It does not supply $J^\pm$, the transmitting hypersurface embedding, or
the additional boundary-jet data of the higher-derivative action. The
existing bulk ADM vertices on their declared spatial domain do not specify
these tilted, frequency-mixing derivatives.

**OPEN:** the physical branch jets, remaining boundary derivative and
geometric pullback are absent from the locked inputs. The
[record](../results/development/nsc-transmitting-metric-pullback.json)
therefore keeps each as `null`. The physical two-sided mismatch is also
`null`, with target tolerance $3\times10^{-11}$. The diagnostic local maximum
$93.54264532195464$ is preserved; it is neither cancelled nor replaced by
the test covectors. Independent first boundary-jet slots remain explicitly
unavailable.

The Gaussian contribution is counted once. The conditional differential does
not add a new action, refit a coefficient, or claim the remaining local and
reference derivatives are zero.
