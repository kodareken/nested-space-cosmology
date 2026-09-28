# Joint linear gain inside the measured quadratic term

Best measured point remains `ad759424…`, residual
`(2.61952127e-09, 7.29271610e-08)`, prediction minus measurement
`(5.79501185e-16, -6.83770152e-16)`. Its step norm
`2.58599619e-10` is the largest norm at which a measured step stayed
linear to about `6e-16`.

On that Jacobian, no further step was evolved. Rank 56, 64 and 70,
damped to this norm, raise both residuals. Rank 74, 80, 88 and 94 have
raw norms `0.3002`, `0.7905`, `12.812` and `13.764`. Damped to the same
norm, their N gains are `2.597e-18`, `8.547e-19`, `5.291e-20` and
`4.922e-20`. Modes 74–94 have N Lipschitz constant `1.34117058e-7`.
Against the measured coefficient `9.10350293e7` from `baf90351…`, that
tail meets its quadratic term at norm `1.47324672e-15`, where the N
gain is `1.97587516e-22`.

In the first 48 equilibrated modes the min-L1 Euclidean step at the
same norm can hold beta and predict an N gain `4.60865858e-12`.
Requiring beta to fall by `1e-10` leaves an N gain
`4.49339876e-12`. The quadratic term of the measured coefficient at
this norm is `6.08785497e-12`, above both gains. The joint slope
`0.0173758909` meets that coefficient at norm `1.90870383e-10`, where
the N gain is `3.31654294e-12` and the residual remains
`2.61620472e-09`. Half of that crossover, where a pure quadratic model
would net the most, nets `8.29135736e-13`. None of these gains reaches
`3e-11`.

The coefficient belongs to the rejected five-times step. It was not
remeasured on these unevolved directions. The search does not exclude
the declared `(w,U)` class: the joint linear gain is positive, and the
continuous cokernel of compact `(w,U)` geometry is already recorded as
trivial. Named stall: `n64_primal_norm_joint_gain_inside_quadratic`.
UV, field, low/subgap and full between-node remain `None`.

```sh
python3 scripts/derive_nsc_ks_n64_primal_norm_joint_quadratic.py --write
python3 scripts/derive_nsc_ks_n64_primal_norm_joint_quadratic.py --check
```
