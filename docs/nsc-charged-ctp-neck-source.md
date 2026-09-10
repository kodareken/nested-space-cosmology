# Charged angular CTP source on the horizon-penetrating neck

This is the direct implementation of the locked source gate.  It keeps

$$
q=4,\qquad
\Omega=3.973074368754331,\qquad
\zeta=15.785319939652627,\qquad
V_{\rm full}=0
$$

fixed.  No flux sector or scale is searched.

## Operator and covariance

For one charged four-dimensional Dirac zero field, magnetic flux gives the
positive angular levels

$$
\lambda_n=\sqrt{n(n+|q|)},
\qquad d_n=2(|q|+2n),\qquad n\ge1,
$$

plus the exact massless LLL of degeneracy $|q|$.  Each massive level uses the
actual exterior reflection, the horizon-regular spin frame and canonical
unitary transport to the stored Bronnikov neck.  The covariance combines the
parent affine occupation with the inherited child occupation

$$
f_p(\omega)=\frac1{1+e^{2\pi\omega/\kappa_p}},
\qquad
f_c(\omega)=\frac1{1+e^{2\pi\omega/(\Omega\kappa_p)}}.
$$

Its finite eigenvalues remain in $[0,1]$ to roundoff.  The calculation uses
the Lorentzian covariance and does not continue a Euclidean spatial heat
operator through the horizon.  The existing opposite-parity copies and their
anomaly cancellation are unchanged.

The massive angular sum uses the previously adopted fourth-order adiabatic
subtraction.  The locked compact complement contributes its local Weyl term;
its vacuum and gauge terms have zero radial-null contraction, while Euler and
total-derivative terms have no interior bulk source.

## Retained null result

At 12 charged angular levels, 12 quadrature points per frequency interval and
5,000 transport steps, the retained four-component source is

$$
(\rho,T_{01},p_\parallel,p_\perp)
=(0.04681379270,\ 0.00122387016,\ -0.10607899355,\ 0.03025113955),
$$

with

$$
T_{++}=-0.05681746054,
\qquad
T_{--}=-0.06171294116.
$$

Both radial null components are negative.  Angular refinement to 16 levels
and a separate frequency/time refinement retain both signs with margins of
order $10^{-2}$.  Scattering-current defects are below $3\times10^{-9}$ and
the covariance/evolution checks remain at their stated tolerances.

The channel decomposition is informative: the scale-inherited LLL is positive
in both null directions, while the massive charged angular tower and the
locked compact Weyl response are negative.  Their retained sum passes the
sign test without changing $q$, $\Omega$, $\zeta$ or the vacuum law.

## Hard gate status

The full requested gate is **not yet accepted**.  The first positive compact
Dirac levels lie in the physical cutoff window, but their nonlocal Lorentzian
CTP covariance and four ADM variations on the horizon domain have not been
evaluated; only their locked local Wilsonian response is present.  That exact
channel—not a new scale, compensator, MMP solve or metric ansatz—is the sole
remaining owner before the sign result can be promoted to the complete
charged angular+compact source.

The focused command is

```sh
python3 scripts/derive_nsc_charged_ctp_neck.py --check
```

It runs only this new charged source calculation and its two declared
refinements.  It does not invoke any historical generator.
