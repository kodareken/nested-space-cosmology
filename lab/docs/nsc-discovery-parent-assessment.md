# Parent episode assessment

Read-only measurement of the completed magnetic parent episodes at producer commit `8c47d0c17fa5962f1d2210f55db3a96c6457a82e`. The consumer is [assess_nsc_discovery_parent_episode.py](../scripts/assess_nsc_discovery_parent_episode.py). It authenticates each snapshot NPZ SHA-256 and the dtype-and-shape array hash, and it reads producer bytes from that commit. It does not import the live parent-episode module, resign momenta, or rebuild the source.

```sh
python scripts/lab.py scripts/assess_nsc_discovery_parent_episode.py
```

Stdout is the record. `--output PATH` creates one new file, mode `0444`, at most 64 MiB, outside the episode directories. Two identical reads on this tree gave content SHA-256 `79dbc6a7bca6010c160a92784dd5d46f57e56ba29b79e6e9d77eb5b79f992bdf`.

## Finding

On the plus-balanced branch the child contracts and depletes, and the NF256 state agrees with NF128 on that matter while the late curvature does not settle. From \(T=0\) to \(T=3\) at \(n_f=128\), the child proper length falls from \(1.0449\) to \(0.0162013156\) and the minimum lapse \(N=rQ\) falls from \(0.99984\) to \(0.006481\). The child probability fraction falls from \(0.37292\) to \(0.059278\). The same fraction, length, and centre proper time at \(T=3\) move only from \(0.0592781809\) to \(0.0590581862\), \(0.0162013156\) to \(0.0162016671\), and \(1.47232277\) to \(1.47221932\) when \(n_f\) goes from 128 to 256. The projected \(|R_4|\) maximum falls from \(4741.23\) to \(982.58\), and the weighted Euler residual \(-8\pi A r^3 Q^2 R_4\) falls from \(0.25698\) to \(0.04995\). That drop is unresolved: \(982\) is not a continuum value, and the raw same-point Hamiltonian \(|R_4|\) stays below \(4\times 10^{-4}\). It is a diagnostic and is not the metric.

The parent-heavy scalar series, both signs, passes a child-fraction trough near \(0.174\) at \(T=1.25\) and a maximum near \(0.251\) at \(T=2.25\), then is lower at the \(T=3\) snapshot (\(0.106\) plus, \(0.109\) minus). Those intermediate times have no full state, so they have no column ancestry. Ancestry at \(T=0,1,3\) is the evolved field times the original column weights; the two columns sum to the actual child probability, and they are not particles. The rise through \(T=2.25\) is not an autonomous renewal. Across every completed case the Hamiltonian total stays inside \(6\times 10^{-13}\) while the field energy is about \(0.018\), and the Gram distance and CAR eigenvalues are unchanged through \(10^{-13}\) and \(10^{-15}\). The saved normal-event ledger is a trapezoid indicator, not a bound. Child normal-energy rates are present on the series and at the snapshots; no separate integrated child-energy stock was stored.

No \(n_f=512\) episode is in these directories. Of the NF256 confirmation list, only plus balanced reached \(T=3\).

## Unmeasured

Column ancestry exists only on the full snapshots \(T=0,1,3\). The other five NF256 populations and every NF512 case are absent. Raw jets, the ledger quadrature, a continuum limit, and a physical \(R_4\) are outside this measurement.
