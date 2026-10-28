# Fixed-parent likelihood scoreboard and SI route ladder

Checkpoint: 485b3b9. Charge, flux, curvature and source normalization are unchanged. P004/P005 remain frozen.

| Comparison | NSC parent | Flat LambdaCDM | Fitted quantities in each |
| --- | ---: | ---: | --- |
| Pantheon+ calibration, 1322 light curves (fixed checkpoint result) | 1179.38727 | 1178.20879 | 1 shape + 1 magnitude intercept |
| DESI DR2, all 13 entries, one shared ruler; chi-square at joint optimum | 18.84303 | 10.78167 | 1 ruler added; 3 total in joint fit |
| SN contribution at that joint optimum | 1184.38115 | 1180.39513 | same 3 joint quantities |
| Joint SN+DESI minimum chi-square | 1203.22418 | 1191.17680 | 3 total; 1335 data values |
| DESI posterior-predictive deviance, SN shape posterior integrated | 27.24410 | 16.46889 | 1 common ruler optimised; not a chi-square variable |
| Conditional Planck distance diagnostic at the joint shape, unknown sound ratio profiled | 1408.82 | 3.04 | 1 extra diagnostic nuisance; not a complete CMB fit |

Delta joint chi-square (NSC minus flat) = **12.04739**. Delta posterior-predictive BAO deviance = 10.77521.

Joint best-fit NSC (Omega_Lambda, Omega_m, Omega_k) = (0.76310716, 0.32320464, -0.08631180), A_o=1.46560035, r_d h=101.91822 Mpc.
Joint best-fit flat control = (0.69630912, 0.30369088, 0), r_d h=101.05339 Mpc.

The posterior r_d h 16/50/84-percent values are [101.22911223704735, 101.90797159479726, 102.58829132982252] Mpc for NSC and [100.36031248417756, 101.0333715586062, 101.7080508133851] Mpc for the control. The complete one-dimensional SN posterior is propagated, with stated priors and a common flat positive-beta prior. No source coupling is fitted.

NSC 3*Omega_Lambda: best fit 2.289321, posterior 2.288667 +/- 0.021121; against 2.054 +/- .022 the standardized difference including posterior width is 7.695.
Flat control: best fit 2.088927, posterior 2.087833 +/- 0.023993; standardized difference 1.039.

The original six AP scores (7.15813 and 8.89401) remain unchanged. Requiring one common ruler is a stronger new comparison. The new joint chi-square does not overwrite the fixed SN calibration score. Omega_m in NSC denotes the scalar pressureless comparison, not a derived baryon/CDM split.

The CMB row is conditional on a supplied z_star=1089.92, the radiation-free scalar/dust continuation, and the curvature-relaxed compression reconstructed from Chen et al. Table 1. It profiles an unknown r_s(z_star)/r_d rather than inventing one. It is not added to the primary SN+DESI total. A full physical CMB likelihood is **undefined** until a compatible photon/baryon/recombination state is supplied.
The conditional geometric R values are 1.55147562 and 1.75179174; the common compression mean is 1.7429. The NSC unaveraged scalar enclosure is [1.5514009663585506, 1.5515502848626992], giving an optimistic phase-adjusted diagnostic lower bound 1407.72. This does not constitute a general test against nesting.

## Route ladder

1. **Flux/charge quantization: q/e=3; screened vacuum levels 9/32, 1/8, 1/32, 0.** q/e=3 and the charge steps leave the fundamental SI charge and reference length free. **No inherited SI selection.**
2. **Wall EOS: frozen stationary radii 1.846303 and 1.155844; V_second=-1.28718 and -11.32694.** The fixed membranes are evolving initial data, not an inherited stable equilibrium selecting metres; the stable surface-scalar branch has a different action/state. **No inherited SI selection.**
3. **Recursion: vacuum-length ratios 1.5 and 2; flux fixed-point residual -1/2.** The supplied map fixes no root length, level or full M/sigma/ell fixed point; its nonzero charge shift has no flux fixed point. **No inherited SI selection.**
4. **Scalar clock: Lambda*c^2/omega^2=2.8125e-19; calibrated matter-vacuum equality redshift .331586.** The epoch and omega/H ratios are dimensionless; no measured SI scalar clock is inherited. **No inherited SI selection.**
5. **Formal compact period: Lambda*lP^2=1/n; illustrative n=1 gives 3.8280731e69 m^-2 and L_ref=5.6795331e-35 m.** The formal period gives Lambda*lP^2=1/n, but neither an admissible parent cycle nor its winding is selected; n=1 gives a Planck-scale value. **No inherited SI selection.**
6. **Transition weights: pure-flux action coefficients .6956584, 3.1684104, 70.9872067 multiply L_unit^2/kappa_nat.** The action exponent, determinant/state and parent preparation are not a scale-selection equation; isolated vacuum bounces are not this interacting parent. **No inherited SI selection.**
7. **Occupation: an extra one-quantum assumption would give 9.17759e61 m^-2.** An inherited occupation/entropy count and its quantum state are absent; setting one quantum changes the prepared classical source. **No inherited SI selection.**
8. **Anomaly/running: prior added classical trace is zero; a quantum stationary scale still needs its action/state/subtractions.** The compatible renormalized quantum action, state and subtraction conditions do not fix b or an absolute vacuum; old spectral closures are different actions. **No inherited SI selection.**
9. **Induced/Casimir/compactification: alternative stable body gives Lambda*M^2=.1144906.** Boundary state, compactification modulus or source mass in SI remains supplied; the prior branches do not select the fixed parent ruler. **No inherited SI selection.**
10. **Infinite hierarchy: ell_n=r^n ell_0.** An inherited ratio or infinite operator composition leaves ell_0 and observer level unselected. **No inherited SI selection.**
11. **Legacy spectral fixed point: threshold 3*pi/2, alternative coefficient 18*pi, supplied length_anchor.** The legacy coefficient dictionary retains a length_anchor and has no solved physical fixed point; its action is not the fixed scalar/four-form parent. **No inherited SI selection.**

**Surviving SI-fixed route: none in this evaluated ladder.** This does not exclude a future completed quantum/global preparation, and does not reject nesting.

The fixed nonzero coefficient is **Lambda*L_ref^2=12.348252398929**. Matching the target by inversion would require L_ref=3.3651918430e+26 m; that is a target-derived diagnostic, not a prediction. The compact-period illustration is 3.51070535e+121 times the target and leaves its winding unselected. The requested Lambda_obs derivation remains **unresolved**.

Sources, input ownership and reproduction:
- [Shared-ruler method and posterior](nsc-atlas-shared-ruler.md); [record](../results/summaries/nsc-atlas-shared-ruler-v1.json).
- [Scale derivations and one-line dispositions](nsc-atlas-scale-selection.md); [record](../results/summaries/nsc-atlas-SI-scale-ladder-v1.json).
- [Conditional CMB record](../results/summaries/nsc-atlas-cmb-geometric-diagnostic-v1.json); [Planck distance compression](https://arxiv.org/html/1808.05724v1).
- [Stationary wall control](../results/summaries/nsc-atlas-wall-scale-stationary-v1.json).
