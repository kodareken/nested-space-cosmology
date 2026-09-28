# Proper-volume and child-clock projection of the completed neck tensor

This step takes the completed charged angular plus compact CTP tensor as its
input.  It does not calculate another source.  The scale branch remains

$$
q=4,\qquad \Omega=3.973074368754331,\qquad
\zeta=15.785319939652627,\qquad V_{\rm full}=0.
$$

## Imported conservation ledger

On the stored Bronnikov child branch,

$$
ds^2=dT^2-a_\parallel(T)^2dz^2-r(T)^2d\Omega_2^2,
\qquad dT=-\frac{d\rho}{\sqrt{-A}},
$$

and the proper volume per coordinate interval $dz$ is

$$
V=4\pi a_\parallel r^2,
\qquad \frac{\dot V}{V}=H_\parallel+2H_\perp.
$$

The standard homogeneous conservation equations are reused from the existing
child-state owner:

$$
\dot\rho+H_\parallel(\rho+p_\parallel)
 +2H_\perp(\rho+p_\perp)=Q,
$$

$$
\dot j+2(H_\parallel+H_\perp)j=S_j.
$$

No GR conservation identity or Friedmann equation is rederived here.

## Neck deposition ledger

At the minimum sphere, $\rho_{\rm coordinate}=0$,

$$
a_\parallel=\sqrt{3\pi/2-1}=1.92675607703,
\quad H_\parallel=1.55702116929,
\quad H_\perp=0,
$$

$$
V_0=24.21233094747,\qquad \dot V_0=12\pi.
$$

Insert the already completed child-frame tensor

$$
(\rho,j,p_\parallel,p_\perp)
=(0.10074836289,\ 0.00122387016,\ -0.34978809527,\ 0.12958962232).
$$

The accepted realization is free and block diagonal in compact level, so its
regular inter-block energy transfer is $Q=0$.  Covariant conservation then
fixes the first child-time density derivative:

$$
\boxed{\dot\rho_0=0.38776013531.}
$$

Its two contributions are

$$
-\frac{\dot V}{V}\rho=-0.15686733379,
\qquad
-H_\parallel p_\parallel-2H_\perp p_\perp=0.54462746910.
$$

Thus the proper energy of a unit-$dz$ cell is

$$
E_0=V_0\rho_0=2.43935270463,
\qquad
\dot E_0=13.18670052493.
$$

This increase is the work of the negative directional pressure during volume
expansion.  It is not an added fluid source.

## What the parent power becomes in child variables

Inside the horizon, the stationary parent coordinate is the homogeneous
spatial coordinate $z$.  With future child time defined by
$dT=-d\rho/a_\parallel$, the recorded outward parent Killing power is the
negative of the child spatial-Killing momentum charge:

$$
P_z=a_\parallel Vj,
\qquad P_{\rm parent}=-P_z=-4\pi a_\parallel^2r^2j.
$$

Numerically,

$$
P_z=+0.057095079694875,
\qquad P_{\rm parent}=-0.057095079694875,
$$

which agrees with the stored parent Killing power to
$4.4\times10^{-16}$.  The momentum equation gives $\dot P_z=0$ for the
homogeneous free bulk.  The parent power is therefore not a continuing child
energy-deposition rate $Q$ in this chart.

## Spectrum allocation and the exact scope boundary

The recorded functional allocation reconstructs the tensor from two pieces:

| Allocation | $\rho$ | $p_\parallel$ | $p_\perp$ |
|---|---:|---:|---:|
| Compact-zero charged angular sector | 0.00629636705 | -0.03213569135 | -0.02697922438 |
| Positive compact levels plus locked Wilsonian complement | 0.09445199584 | -0.31765240392 | 0.15656884669 |

Their sum reproduces every completed tensor component exactly at stored
precision.  No interaction vertex is enabled between these free compact
blocks, so the derived regular inter-block $Q$ is zero.  This allocation is a
spectral/renormalization split; it does not identify either row with visible
matter, dark matter, or dark energy.

Define the anisotropic volume scale $a_V=(V/V_0)^{1/3}$.  The result supplies
the local Cauchy jet

$$
\rho(T)=0.10074836289+0.38776013531T+O(T^2),
$$

$$
\rho(a_V)=0.10074836289
 +0.74711919715\log a_V+O((\log a_V)^2).
$$

A single neck tensor does not determine $\dot p_\parallel$ or
$\dot p_\perp$.  Consequently it does not determine a global
$(\rho(a),p(a))$ history.  The next owner is metric backreaction together
with evolution of this same CTP state and tensor on the locked child branch;
only that evolution can extend the local jet without inserting an equation of
state.

The focused reproducer is

```sh
python3 scripts/derive_nsc_background_projection.py --check
```

It authenticates the completed source and existing geometry/ADM records.  It
does not run a source, angular, compact, MMP, metric, or historical generator.
