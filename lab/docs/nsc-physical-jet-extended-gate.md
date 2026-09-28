# J→M: partial Dirac jets computed, physical transmission OPEN

The [new gate](../results/development/nsc-physical-jet-extended-gate.json)
composes the Hamiltonian-to-endpoint derivative check with the existing
reference, local, interface and Gaussian-action records.

| Step | Executed result | Status |
|---|---|---|
| J | Known bulk Hamiltonian vertices propagated into the existing endpoint-jet object; new numerical chain-rule checks | Bulk component PASS / physical jets OPEN |
| K | Same endpoint object and remainder allocation retained | OPEN: physical boundary derivative unevaluated |
| L | Physical selector prerequisites checked | OPEN: no transmitting jets or complete boundary variation |
| M | Extended dependency gate | **OPEN** |

The new computation verifies the derivative path from the owned Hamiltonian
to the owned CTP differential. It does not recreate the previous norm-family
or determinant-differential experiments. Its one frozen control is a
derivative regression only; it is not an attempted stationary history in the
already excluded homogeneous class.

To evaluate the physical jets, the action/domain must specify the
transmitting $B[g,\mathrm{embedding}]$ and Cauchy restriction/normal map.
To form the physical boundary mismatch, the same action must supply its
remaining boundary derivative in the existing metric and normal-jet basis.
Neither function is determined by the supplied diagnostic vertices.

No optimizer is invoked and no extended existence or non-existence is
claimed. The old homogeneous non-existence certificate and Gaussian
differential checks are imported unchanged. All ledger coefficients and
scales remain locked; finite stress, nulls and updated constraints remain
unevaluated. Metric stepping and coupled evolution remain closed.
