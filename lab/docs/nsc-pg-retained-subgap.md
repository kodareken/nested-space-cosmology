# The compact subgap covariance panels

All **38 declared signed families in the 20 compact groups** now have their
$[1,m_j]$ covariance panels, with $m_j=j\pi/2$. Both angular signs are
evaluated where the compact/angular representation requires them. The
zero-angular sectors are not doubled.

The [record](../results/development/nsc-pg-retained-subgap.json) preserves the
positive-energy Gram and centered-covariance matrices, their physical signed
combinations, the quadrature receipts and the locked action/scale ledger.
The source remains the existing horizon covariance and inherited occupation.

## A stable representation of the same state

The [subgap contour identity](nsc-pg-threshold-projection.md) removes the
rapid real-energy winding near the mass threshold. For large angular labels,
directly adding its horizon-basis terms is poorly conditioned: in group 22,
terms of order $1.9\times10^9$ cancel to an answer of order one.

The same source can instead be expressed through the bounded whole-line
retarded packet resolvent. With $w$ the interior horizon partner and $a$ the
exterior-horizon source column continued into the interior,

$$
K(E)=(f-\tfrac12)G_E-2(f-\tfrac12)ww^\dagger
+is(wa^\dagger-aw^\dagger),
\qquad
G_E=\frac{R^+(E)-R^+(E)^\dagger}{i}.
$$

This is a rearrangement of the original horizon covariance. It preserves
its coherence and introduces no additional field, counterflow or subtraction.
The [packet-resolvent owner](../src/recursive_horizons/nsc_pg_retarded_packets.py)
uses the same global ingoing/outgoing boundary data and retains the causal
zero for exterior response to trapped-region sources. Its first four packet
entries agree with the existing trapped resolvent control.

The [adaptive owner](../src/recursive_horizons/nsc_pg_adaptive_threshold.py)
integrates this expression on the same pole-free frequency contour. Each
family reports its own estimated error and integration status. Those estimates
are numerical error diagnostics, not a new cosmological parameter or a proof
of the complete high-energy remainder.

## What has closed

Every declared family has a successful integration status. The largest sum
of adaptive error estimates is $3.14\times10^{-11}$. Occupied and complementary
positive-energy matrices are positive within numerical roundoff; their largest
negative excursions are below $1.7\times10^{-13}$. Signed assembly is checked
against the stored opposite-angular family, preserving the physical pairing.

This completes the compact **subgap panels**. The remaining real-energy
windows and high-energy tails are required for each full covariance. The
[first complete massive covariance](nsc-pg-group13-covariance.md) is already
available for group 13. Full retained C1b, physical stress, transmitting
endpoint jets and metric evolution remain separate gates. PDF v0.25.0 and
the earlier records are unchanged.

```sh
python3 scripts/check_nsc_pg_retained_subgap.py --check

# New/cold evaluation of one declared family, only when needed:
python3 scripts/derive_nsc_pg_retarded_threshold.py --channel 14 --sign -1 --adaptive
```
