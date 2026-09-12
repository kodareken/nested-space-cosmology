# The Cauchy propagator is unitary once a history is supplied

The mode-state payload now makes blockwise Dirac propagation executable.  For
any declared homogeneous ADM history, the canonical half-density Hamiltonian
is

$$
H_{jn}=N\left[-m_j\sigma_1+\frac{\lambda_n}{r}\sigma_2
 +\frac{k}{a_\parallel}\sigma_3\right]-\beta k I,
$$

and its time-ordered exponential gives a unitary map.  The implementation uses
exact $2\times2$ exponentials at every midpoint and pushes every physical
covariance forward as

$$
C_L=U C_0U^\dagger.
$$

## Endpoint data do not select the history

The source-selected record fixes endpoint values but does not fix the
intervening lapse, shift, radius, axial scale or embedding of the tilted
Landau surface.  To test whether this omission matters, the new owner evaluates
two smooth diagnostic Kantowski--Sachs histories.  They share every recorded
endpoint, including the retained axial scale and $H_\parallel$, the selected
radius, and the Landau lapse/shift representation.  Their durations differ.

Both maps are unitary, both preserve the full covariance spectrum and both use
the same magnetic and Dirac channel data.  Nevertheless their resulting
$U$ matrices differ.  This is the direct finite witness that

$$
(v,\eta,r_\star,a_{\parallel,0},H_{\parallel,0})
\quad\text{does not determine}\quad U_{L0}.
$$

The two histories are controls, not candidate physical metrics.  Choosing one
would add an unowned duration/profile prescription.  A pointwise spin boost is
not used as the isometry.

## Consequence for the finite stress

Because endpoint-identical histories produce different target covariances,
there is no unique finite $T_{\mu\nu}(r_\star,\Sigma_L)$ yet.  The same
fourth-order reference and local induced allocation also require the selected
metric history, so neither can be evaluated before this choice.  The old
tensor and the linearized density diagnostic are not substituted.

This is a conditional API pass and a physical-selection fail.  The missing
owner is a joint history/state boundary-value solution or a same-action
history-selection equation.  The present sequencing is circular: the history
would normally be produced by coupled metric evolution, while metric stepping
was held until the Cauchy map existed.

The record and its content-addressed control payload are reproduced by

```sh
python3 scripts/derive_nsc_landau_cauchy_isometry.py --check
```

The command loads the existing mode-state artifact, computes both blockwise
maps and checks every payload array.  It performs no source, MMP or metric-field
solve.
