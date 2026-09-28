# N nonlinearity of the second beta-primary step

Parent `6534e55b…` measured `(2.61985695e-09, 7.30697217e-08)`.
The 48-mode step of norm `9.99927046e-08`, identity `b11573aa…`, predicted
`(2.59876829e-09, 7.00697217e-08)` and measured
`(2.78238122e-09, 7.00743079e-08)`.

Beta followed the Jacobian. N rose, and the N prediction error
`1.83612937e-10` is larger than the predicted N gain `2.10886620e-11`.
The step is not accepted. Best measured point remains `6534e55b…`.

```sh
python3 scripts/derive_nsc_ks_n64_primal_norm_beta_nonlinearity.py --write
python3 scripts/derive_nsc_ks_n64_primal_norm_beta_nonlinearity.py --check
```
