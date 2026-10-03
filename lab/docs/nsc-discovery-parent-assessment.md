# Parent episode assessment

Read-only measurement of the completed magnetic parent episodes at producer commit `8c47d0c17fa5962f1d2210f55db3a96c6457a82e`. The consumer is [assess_nsc_discovery_parent_episode.py](../scripts/assess_nsc_discovery_parent_episode.py). It authenticates each snapshot NPZ SHA-256 and the dtype-and-shape array hash, and it reads producer bytes from that commit. It does not import the live parent-episode module, resign momenta, or rebuild the source.

```sh
python scripts/lab.py scripts/assess_nsc_discovery_parent_episode.py
```

Stdout is now assessment v2. `--output PATH` creates one new file, mode `0444`, at most 64 MiB, outside the episode directories. Two identical reads gave content SHA-256 `9febeebd6d71e91ec965ca0c6a9d99cbff972789b8bb44453fb7ac69620e730c`. The sealed v1 file remains byte-unchanged, content digest `79dbc6a7bca6010c160a92784dd5d46f57e56ba29b79e6e9d77eb5b79f992bdf`, and its campaigns and comparison match this successor (`measurement_sha256` `1f18663db29612d0e2412d4b6997c2c6942914173361e9fd3ab7e2748e4fe5f7`).

## Finding

On the plus-balanced branch the child contracts and depletes, and the NF256 state agrees with NF128 on that matter while the late curvature does not settle. From \(T=0\) to \(T=3\) at \(n_f=128\), the child proper length falls from \(1.0449\) to \(0.0162013156\) and the minimum lapse \(N=rQ\) falls from \(0.99984\) to \(0.006481\). The child probability fraction falls from \(0.37292\) to \(0.059278\). The same fraction, length, and centre proper time at \(T=3\) move only from \(0.0592781809\) to \(0.0590581862\), \(0.0162013156\) to \(0.0162016671\), and \(1.47232277\) to \(1.47221932\) when \(n_f\) goes from 128 to 256. The projected \(|R_4|\) maximum falls from \(4741.23\) to \(982.58\), and the weighted Euler residual \(-8\pi A r^3 Q^2 R_4\) falls from \(0.25698\) to \(0.04995\). That drop is unresolved: \(982\) is not a continuum value, and the raw same-point Hamiltonian \(|R_4|\) stays below \(4\times 10^{-4}\). It is a diagnostic and is not the metric.

The parent-heavy scalar series, both signs, passes a child-fraction trough near \(0.174\) at \(T=1.25\) and a maximum near \(0.251\) at \(T=2.25\), then is lower at the \(T=3\) snapshot (\(0.106\) plus, \(0.109\) minus). Those intermediate times have no full state, so they have no column ancestry. Ancestry at \(T=0,1,3\) is the evolved field times the original column weights; the two columns sum to the actual child probability, and they are not particles. The rise through \(T=2.25\) is not an autonomous renewal. Across every completed case the Hamiltonian total stays inside \(6\times 10^{-13}\) while the field energy is about \(0.018\), and the Gram distance and CAR eigenvalues are unchanged through \(10^{-13}\) and \(10^{-15}\). The saved normal-event ledger is a trapezoid indicator, not a bound. Child normal-energy rates are present on the series and at the snapshots; no separate integrated child-energy stock was stored.

No \(n_f=512\) episode is in these directories. Of the NF256 confirmation list, only plus balanced reached \(T=3\).

## Unmeasured

Column ancestry exists only on the full snapshots \(T=0,1,3\). The other five NF256 populations and every NF512 case are absent. Raw jets, the ledger quadrature, a continuum limit, and a physical \(R_4\) are outside this measurement.

## Responsive controls

```sh
python scripts/lab.py scripts/assess_nsc_discovery_parent_episode.py --responsive-controls
```

This path is separate from the default successor. Two reads gave content SHA-256 `05eb4e4a45d0199fb95066d69f1e4fbc030c22e2673dfc52d8d4ed63232600e1`. Every quoted length, fraction, clock, and curvature below was recomputed from the saved arrays. The strong \(T=3\) checkpoints at \(T=1.25\), \(1.5\), and \(2.25\) still carry the \(T=1\) observation header; that header is not the measurement.

The strong balanced run is sign-dependent, and it is not a source-force-only experiment. Its occupations are \(173.16\) times the magnetic plus-balanced handoff. \(Q\), \(r\), and \(\Phi\) match that handoff, while canonical \(\pi_r\) differs by \(0.2797\). The authoritative parent-record \(k\) is \(0.002621213673327129\) on the weak handoff and \(0.45388971484879426\) on the strong handoff. The weak episode `common_k` of about \(1\times 10^{-8}\) is a legacy label and is not that comparison. No completed episode changes only the weights or only \(k\).

From the arrays, plus goes from child fraction \(0.32617\) and proper length \(0.31386\) at \(T=1\) to \(0.059670\) and \(0.00053105\) at \(T=3\). Its centre proper time reaches only \(0.70642\), and the projected \(|R_4|\) maximum is \(5.30275\times 10^{7}\) against a raw same-point reference of \(1.036\). Minus goes from fraction \(0.34039\) and length \(1.76097\) to \(0.090496\) and \(0.14933\), with centre proper time \(3.05079\) and projected \(|R_4|\) \(18.798\). Those projected values are not an instability and they set no new threshold.

The frozen plus parent-heavy control has the same canonical handoff as the magnetic plus parent-heavy state, gap \(0\). \(Q\) and \(r\) do not move through \(T=3\), so the child length stays \(1.04491\). The child fraction is \(0.21242\), \(0.21889\), \(0.33464\), and \(0.29834\) at \(T=0\), \(1.25\), \(2.25\), and \(3\). Projected and raw \(|R_4|\) agree at \(12.784\). The coupled parent-heavy depletion is not this fixed-geometry motion.

The empty-source control keeps the balanced weights, deletes \(\Phi\), and re-solves the momenta (\(\pi_r\) gap \(0.01816\)). The child fraction stays \(0\). By \(T=3\) the child length is \(0.016335\) and the projected \(|R_4|\) maximum is \(4768.1\), close to the coupled balanced \(0.016201\) and \(4741\). That metric contraction is not carried by the source force. Hamiltonian totals and CAR eigenvalues stay at the cancellation residual in every one of these controls. The ledger quadrature remains an indicator, not a bound.
