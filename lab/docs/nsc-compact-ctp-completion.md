# Positive compact levels complete the charged CTP neck gate

The preceding charged angular record left one named channel: the nonlocal CTP
covariance and four metric projections of the first positive compact Dirac
levels.  This calculation adds exactly that channel while keeping

$$
q=4,\qquad \Omega=3.973074368754331,\qquad
\zeta=15.785319939652627,\qquad V_{\rm full}=0
$$

fixed.

## Canonical positive compact fields

The established free compact reduction gives two four-dimensional Dirac
copies at every

$$
m_j=\frac{j\pi}{2L_\star}.
$$

At the locked cutoff, $j=1,2$ are the positive levels inside the physical
window.  Their magnetic angular modes use the same $q=4$ spectrum, actual
massive exterior reflection, horizon frame and inherited Gaussian occupation
as the zero compact level.  The child threshold scales with $\Omega$, so no
incoming reservoir is inserted.

The CTP covariance is subtracted with the negative-energy superadiabatic
projector through fourth derivative order.  In the $m_j\to0$ limit its energy
and axial-pressure projections reproduce the existing $E_2+E_4$ and
$P_2+P_4$ subtractions exactly.  This removes the invalid instantaneous-vacuum
tail.  The already locked Wilsonian compact term restores the local response
once, and the new calculation supplies only its nonlocal state/geometry part.

## Completed source and gate

The completed child-frame tensor is the sum of the authenticated charged
angular/local record and the new positive-level contribution.  At the selected
compact resolution it is

$$
(\rho,T_{01},p_\parallel,p_\perp)
=(0.10074836289,\ 0.00122387016,\ -0.34978809527,\ 0.12958962232),
$$

with

$$
\boxed{T_{++}=-0.24659199207,\qquad
T_{--}=-0.25148747269.}
$$

Its two radial null components remain negative under compact-angular,
frequency-extent and frequency/time refinement.

This passes the requested neck gate in the declared free Gaussian one-loop
realization.  It introduces no metric ansatz, compensator, scalar, flux-sector
change or scale refit.  The already calculated interacting compact vertex is
not enabled because its torsion stiffness and interacting state are not fixed;
that is a different realization rather than a hidden adjustment of this one.

Passing the local gate permits the homogeneous projection to be recorded.  It
does not yet define $H(z)$: conservation, proper-volume/clock conversion and
component deposition still have to be applied to this same tensor.

The local homogeneous projection gives

$$
p_{\rm iso}=-0.03020295021,
\qquad w_{\rm iso,neck}=-0.29978601484,
$$

with nonzero anisotropic stress.  These are neck data, not a fitted cosmic
equation of state.

The focused reproducer is

```sh
python3 scripts/derive_nsc_compact_ctp_completion.py --check
```

It reads the previous charged angular result and does not rerun it, MMP, or any
historical generator.
