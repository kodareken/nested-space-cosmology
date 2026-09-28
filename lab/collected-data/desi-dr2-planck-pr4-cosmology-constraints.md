# DESI DR2 cosmological constraints and Planck directional-isotropy limits

Fetched from primary web sources on 2026-08-15.

Original search topics:

- `DESI DR2 2025 BAO full-shape best fit Omega_m Omega_Lambda H0 constraints arXiv 2503.14738 cosmological parameters`
- `Planck 2018 limits direction-dependent spectral index n_s hemispherical asymmetry tilt constraints isotropy paper PR4`

## Executive summary

Two terminology distinctions are essential:

1. [arXiv:2503.14738](https://arxiv.org/abs/2503.14738) is the DESI DR2 **BAO cosmology** analysis, not a DR2 full-shape analysis.
2. “Planck 2018” refers to PR3, whereas Planck PR4/NPIPE is a later reprocessing. Results from the Planck 2018 papers and later PR4 analyses must therefore be quoted separately.

The main DESI DR2 BAO+CMB flat-ΛCDM result is

\[
\boxed{\Omega_m=0.3027\pm0.0036,\quad
\Omega_\Lambda\simeq0.6973\pm0.0036,\quad
H_0=68.17\pm0.28\ \mathrm{km\,s^{-1}\,Mpc^{-1}}.}
\]

Here \(\Omega_\Lambda\) is derived from flatness, not directly tabulated.

For directional variation of the scalar spectral index, the most direct Planck PR4 patch-based analysis reports an \(n_s\)-dipole probability-to-exceed of 5.8% with temperature plus polarization and 38.8% with temperature alone. It does not claim a direction-dependent tilt and does not publish a conventional PR4 95% upper limit on hemispheric \(\Delta n_s\).

## 1. DESI DR2 BAO constraints

The principal source is:

- DESI Collaboration, [“DESI DR2 Results II: Measurements of Baryon Acoustic Oscillations and Cosmological Constraints,” arXiv:2503.14738](https://arxiv.org/abs/2503.14738), published in *Physical Review D* 112, 083515 (2025).

For flat \(\Lambda\)CDM, Table V and the associated equations report the following marginalized central values and 68% credible intervals:

| Dataset | \(\Omega_m\) | \(\Omega_\Lambda\), derived | \(H_0\) [km s\(^{-1}\) Mpc\(^{-1}\)] |
|---|---:|---:|---:|
| DESI DR2 BAO | \(0.2975\pm0.0086\) | \(\simeq0.7025\pm0.0086\) | Not separately determined |
| DESI DR2 BAO + BBN | \(0.2977\pm0.0086\) | \(\simeq0.7023\pm0.0086\) | \(68.51\pm0.58\) |
| DESI DR2 BAO + BBN + \(\theta_*\) | \(0.2967\pm0.0045\) | \(\simeq0.7033\pm0.0045\) | \(68.45\pm0.47\) |
| DESI DR2 BAO + CMB | \(0.3027\pm0.0036\) | \(\simeq0.6973\pm0.0036\) | \(68.17\pm0.28\) |

The BAO-only result is more precisely

\[
\Omega_m=0.2975\pm0.0086,\qquad
h\,r_d=(101.54\pm0.73)\ \mathrm{Mpc},
\]

with correlation coefficient \(r=-0.92\). BAO alone measures \(H_0r_d\), not \(H_0\) independently. An early-Universe calibration such as BBN or CMB information is needed to obtain \(H_0\). The [DESI DR2 paper’s ΛCDM section](https://arxiv.org/html/2503.14738#S6) contains the full discussion.

### 1.1 Reported values versus “best fit”

The entries above are marginalized posterior constraints, not the complete maximum-posterior best-fit parameter vector. The paper’s discussion informally rounds the DESI and DESI+CMB best-fitting matter densities to approximately 0.297 and 0.303, respectively, while Table V gives the posterior summaries above.

### 1.2 Derivation of \(\Omega_\Lambda\)

The \(\Omega_\Lambda\) entries in the table are derived rather than directly reported. For flat \(\Lambda\)CDM,

\[
\Omega_\Lambda=1-\Omega_m-\Omega_r.
\]

At the precision shown, \(1-\Omega_m\) is sufficient; the present radiation fraction changes only the fourth decimal place.

### 1.3 Non-flat comparison

When curvature is allowed, DESI+CMB gives

\[
\Omega_m=0.3034\pm0.0037,\qquad
H_0=68.50\pm0.33,\qquad
\Omega_K=0.0023\pm0.0011.
\]

This is not a significant preference for nonzero curvature. The corresponding derived dark-energy density is approximately

\[
\Omega_\Lambda\simeq1-\Omega_m-\Omega_K-\Omega_r\simeq0.6943.
\]

## 2. What “full shape” means in the 2025 DESI material

The relevant 2025 companion paper is:

- DESI Collaboration, [“Constraints on Neutrino Physics from DESI DR2 BAO and DR1 Full Shape,” arXiv:2503.14744](https://arxiv.org/abs/2503.14744).

It analyzes or compares:

- DESI DR2 BAO;
- DESI DR1 full shape plus DR1 BAO;
- several CMB priors and likelihoods.

The authors explicitly state that they did **not** present a joint DR1-full-shape + DR2-BAO analysis in 2025 because the required DR1/DR2 cross-correlation from mock catalogues was not yet available. Therefore, there is no official 2025 combined “DR2 BAO + full shape” \(\Omega_m,H_0\) entry in arXiv:2503.14738 or arXiv:2503.14744. The [limitation is explained in the companion paper](https://arxiv.org/html/2503.14744#S2.SS2).

For the companion paper’s baseline model in which \(\sum m_\nu\) is allowed to vary, DESI DR2 BAO+CMB gives

\[
\Omega_m=0.3009\pm0.0037,\qquad
H_0=68.36\pm0.29,
\]

and

\[
\sum m_\nu<0.0642\ \mathrm{eV}\qquad(95\%).
\]

These numbers differ slightly from arXiv:2503.14738 because neutrino mass is a free parameter.

## 3. Later joint DESI DR1 full-shape + DR2 BAO result

A genuine joint analysis appeared in February 2026:

- D. Forero-Sánchez et al., [“Cosmological constraints from a joint DESI DR1 Full-Shape and DR2 BAO,” arXiv:2602.18761](https://arxiv.org/abs/2602.18761), submitted to JCAP.

It estimates the cross-covariance using mocks. For flat \(\Lambda\)CDM with a BBN prior and a wide \(n_s\) prior, it reports

\[
\Omega_m=0.3035\pm0.0085,\qquad
h=0.6876\pm0.0059,\qquad
\sigma_8=0.822\pm0.034.
\]

Therefore,

\[
H_0=68.76\pm0.59\ \mathrm{km\,s^{-1}\,Mpc^{-1}},
\qquad
\Omega_\Lambda\simeq0.6965\pm0.0085.
\]

This is the closest match to a literal “DESI DR1 full shape + DR2 BAO” constraint, but it is a 2026 preprint rather than arXiv:2503.14738.

## 4. Planck constraints on a direction-dependent \(n_s\)

Three related but non-equivalent quantities must be separated:

1. the ordinary, sky-averaged scalar tilt \(n_s\);
2. a spatial gradient or dipole in \(n_s\);
3. a phenomenological hemispherical modulation of CMB power.

They are not interchangeable.

### 4.1 Ordinary isotropic scalar tilt

The Planck 2018/PR3 inflation analysis finds

\[
n_s=0.9649\pm0.0042
\]

at 68% confidence from temperature, polarization, and lensing, with no evidence for scale dependence of \(n_s\):

- Planck Collaboration, [“Planck 2018 results. X. Constraints on inflation,” arXiv:1807.06211](https://arxiv.org/abs/1807.06211).

Using the later PR4/NPIPE maps and the HiLLiPoP/LoLLiPoP likelihoods, the baseline PR4 TTTEEE analysis gives approximately

\[
n_s=0.9681\pm0.0039.
\]

This remains a sky-averaged isotropic value. The PR4 cosmological-parameter source is:

- M. Tristram et al., [“Cosmological parameters derived from the final (PR4) Planck data release,” arXiv:2309.10034](https://arxiv.org/abs/2309.10034), *Astronomy & Astrophysics* 682, A37 (2024).

### 4.2 Planck 2018 \(n_s\)-gradient model

Planck 2018 tested a model in which the scalar tilt has a linear spatial gradient across the observable volume. Its asymmetry spectrum is parameterized as

\[
C_\ell^{\mathrm{lo}}
=
-\frac{\Delta n_s}{2}\frac{dC_\ell}{dn_s},
\]

where \(\Delta n_s\) is the dipolar modulation amplitude.

The actual Planck 2018 polarization data gave a **57% p-value relative to statistically isotropic polarization simulations** for this \(n_s\)-gradient model. Polarization therefore did not support or confirm a physical direction-dependent tilt. The same analysis returned p-values of 43% and 30% for its power-law and tanh modulation models. The Planck Collaboration concluded that polarization did not support the proposed physical scale-dependent dipolar-modulation models. See section 10.1 of the [Planck 2018 inflation paper](https://arxiv.org/pdf/1807.06211).

A historical temperature-only fit based on Planck 2015 SMICA data had found the model-dependent best-fitting parameters

\[
\Delta n_s=0.014,\qquad
k_*=0.10\ \mathrm{Mpc}^{-1},\qquad
(l,b)=(214^\circ,-17^\circ).
\]

Source:

- D. Contreras et al., [“Testing physical models for dipolar asymmetry with CMB polarization,” arXiv:1704.03143](https://arxiv.org/abs/1704.03143), *Physical Review D* 96, 123522 (2017).

The value 0.014 must not be described as a Planck PR4 limit:

- it was a best-fitting modulation amplitude, not a 95% upper bound;
- it used Planck 2015 temperature data;
- the later Planck 2018 polarization test did not confirm it.

The Planck 2018 paper does not tabulate a conventional marginalized 95% upper limit on \(|\Delta n_s|\).

### 4.3 Direct PR4 test of a dipole in fitted \(n_s\)

The paper most directly matching “direction-dependent \(n_s\), statistical isotropy, PR4” is:

- C. Gimeno-Amo et al., [“Exploring Statistical Isotropy in Planck Data Release 4: Angular Clustering and Cosmological Parameter Variations Across the Sky,” arXiv:2504.05597](https://arxiv.org/abs/2504.05597).

It fits cosmological parameters independently in 12 disjoint sky patches and tests the amplitude of a dipole in each parameter against 600 PR4 end-to-end simulations. Its \(n_s\) results are:

| Analysis | \(n_s\)-dipole probability-to-exceed |
|---|---:|
| Temperature only | 38.8% |
| Temperature + \(E\)-mode polarization | 5.8% |
| Debiased temperature + polarization | 6.0% |
| More restrictive multipole cuts | 10.2% or 12.8% |

The study therefore classifies the \(n_s\) variation as consistent with \(\Lambda\)CDM and does not claim a direction-dependent spectral tilt. Its only potentially anomalous parameter dipole is \(A_s\), with a probability-to-exceed of 0.83%, corresponding to 5 of 600 simulations exceeding the observed amplitude. The [numerical probability-to-exceed table is available in the full paper](https://arxiv.org/html/2504.05597#A1.T3).

The plotted main-case \(n_s\) dipole amplitude is visually about \(1.7\times10^{-2}\), but the paper does not tabulate that value or provide a conventional 95% upper limit. The robust, explicitly reported statistic is the 5.8% probability-to-exceed.

### 4.4 PR4 hemispherical power asymmetry without interpreting it as \(n_s\)

A separate PR4 study measures a local-variance dipole rather than fitting a hemispheric \(n_s\):

- C. Gimeno-Amo et al., [“Hemispherical Power Asymmetry in intensity and polarization for Planck PR4 data,” arXiv:2306.14880](https://arxiv.org/abs/2306.14880), JCAP 12 (2023) 029.

It finds the following.

#### Temperature

- No simulation among 600 had a larger dipole for \(4^\circ\), \(6^\circ\), or \(8^\circ\) discs, hence \(p<0.17\%\).
- Inferred modulation amplitude: approximately 7%.
- Direction at full resolution: approximately \((l,b)=(205^\circ,-20^\circ)\).

#### PR4 SEVEM \(E\)-mode polarization

- \(p=2.8\%\), with a split-dependent range of 1.8%–3.8%.
- Direction: \((l,b)=(234^\circ,-14^\circ)\pm5^\circ\).
- Model-dependent modulation amplitude: approximately 9%, with an inferred 68% range of roughly 6%–13%.
- Temperature–polarization axis-alignment p-value: 4.5%.

The polarization result changes materially with mask choice, reaching \(p=7.4\%\) for one sky fraction. The authors conclude that Planck polarization is not sensitive enough for a robust detection. The [PR4 numerical results and qualifications are given in the paper’s results section](https://arxiv.org/html/2306.14880#S4.SS3).

A 7%–9% hemispherical **power modulation is not a measurement of \(\Delta n_s\)**. A tilt gradient is only one possible scale-dependent physical explanation for such a modulation.

## 5. Compact conclusions

### DESI

DESI DR2 BAO alone, flat \(\Lambda\)CDM:

\[
\boxed{\Omega_m=0.2975\pm0.0086,\qquad
\Omega_\Lambda\simeq0.7025\pm0.0086}
\]

but BAO alone does not independently determine \(H_0\).

DESI DR2 BAO+CMB, flat \(\Lambda\)CDM:

\[
\boxed{\Omega_m=0.3027\pm0.0036,\quad
\Omega_\Lambda\simeq0.6973\pm0.0036,\quad
H_0=68.17\pm0.28.}
\]

Later 2026 joint DESI DR1 full-shape + DR2 BAO:

\[
\boxed{\Omega_m=0.3035\pm0.0085,\quad
\Omega_\Lambda\simeq0.6965\pm0.0085,\quad
H_0=68.76\pm0.59.}
\]

### Planck

Planck 2018 ordinary isotropic tilt:

\[
\boxed{n_s=0.9649\pm0.0042.}
\]

Direct PR4 dipole test of fitted \(n_s\):

\[
\boxed{\mathrm{PTE}=5.8\%\ \text{with }T+E,\qquad
\mathrm{PTE}=38.8\%\ \text{with }T\text{ alone}.}
\]

The authors regard these values as consistent with statistical isotropy. No standard PR4 95% upper limit on a hemispheric \(\Delta n_s\) was reported.

## Primary-source links

1. [DESI DR2 Results II — arXiv:2503.14738](https://arxiv.org/abs/2503.14738)
2. [DESI DR2 BAO and DR1 Full Shape neutrino analysis — arXiv:2503.14744](https://arxiv.org/abs/2503.14744)
3. [Joint DESI DR1 Full-Shape and DR2 BAO — arXiv:2602.18761](https://arxiv.org/abs/2602.18761)
4. [Planck 2018 results X: Constraints on inflation — arXiv:1807.06211](https://arxiv.org/abs/1807.06211)
5. [Planck 2018 results VII: Isotropy and statistics — arXiv:1906.02552](https://arxiv.org/abs/1906.02552)
6. [Final Planck PR4 cosmological parameters — arXiv:2309.10034](https://arxiv.org/abs/2309.10034)
7. [PR4 hemispherical power asymmetry — arXiv:2306.14880](https://arxiv.org/abs/2306.14880)
8. [PR4 cosmological-parameter variations across the sky — arXiv:2504.05597](https://arxiv.org/abs/2504.05597)
9. [Physical models for dipolar asymmetry — arXiv:1704.03143](https://arxiv.org/abs/1704.03143)
