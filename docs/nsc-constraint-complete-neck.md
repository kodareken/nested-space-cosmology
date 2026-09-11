# A source-selected neck closes both initial constraints

The stored unit-radius Bronnikov neck fails the initial ADM gate for the
completed CTP tensor.  This follow-up uses the permitted
**B: geometry-from-source** route.  It keeps

$$
q=4,\quad \Omega=3.973074368754331,\quad
\zeta=15.785319939652627,\quad
A=0.04501936182826115,\quad V_{\rm full}=0
$$

and every state occupation fixed.

## Why source-completion A is not selected

No already evaluated channel owns both recorded residuals.  The homogeneous
magnetic Maxwell contribution is diagonal and cannot supply the required
counterflow.  The relational vacuum is fixed at zero.  The compact interaction
has an undetermined torsion stiffness and interacting state, while the
unevaluated nonlocal remainder has no recorded finite tensor that could be
assigned to the gap.  None of these terms is inserted or tuned.

## The tensor selects its own normal

The radial stress block is type-I because

$$
(\rho+p_\parallel)^2-4T_{01}^2
=0.0620147968724349>0.
$$

It therefore has a unique future, subluminal Landau normal.  The boost from
the old child frame is fixed by

$$
T_{01}(1+v^2)+v(\rho+p_\parallel)=0,
$$

giving

$$
v=0.00491447570673408,
\qquad \eta=0.00491451527223027.
$$

In this source-selected frame,

$$
(\rho_L,T_{01,L},p_{\parallel,L},p_{\perp,L})
=(0.100754377566935,\ 0,\ -0.349782080588978,\ 0.129589622317672).
$$

The stress trace and product of the two radial null contractions are invariant
to machine precision.  Both null components become

$$
T_{++,L}=T_{--,L}=-0.249027703022043.
$$

The boost changes the normal used for the ADM split; it does not change the
state occupations or erase the previously recorded parent charge.

## The density selects the areal radius

Retain the existing Kantowski--Sachs metric class and a temporal minimum,
$H_\perp=0$.  The Hamiltonian constraint then fixes, rather than fits,

$$
\boxed{
r_\star=\sqrt{\frac{2A}{\rho_L}}
=0.945328394434129.
}
$$

The axial scale and $H_\parallel$ are retained from the seed because the two
constraints do not determine either at $H_\perp=0$.  The replacement initial
geometry is therefore

$$
ds^2=dT_L^2-a_{\parallel,0}^2dz_L^2-r_\star^2d\Omega_2^2,
$$

$$
a_{\parallel,0}=1.92675607703328,
\quad H_{\parallel,0}=1.55702116929053,
\quad H_{\perp,0}=0.
$$

It is a source-selected Kantowski--Sachs minimum, not the stored unit-radius
Bronnikov profile and not a new metric ansatz.  The original MMP unit-radius
condition is not imposed on this replacement geometry.

The imported constraints now give

$$
\boxed{
\mathcal C_H=2.22\times10^{-16},
\qquad \mathcal C_M=0,
}
$$

below the declared $2\times10^{-12}$ tolerance.  No change to
$A,q,\Omega,\zeta,V_{\rm full}$ or $H_\perp$ was used.

This passes the initial assembly gate.  The next calculation may now evolve
the same CTP covariance and this source-selected metric together, evaluating
the stress again at every step.

The focused reproducer is

```sh
python3 scripts/derive_nsc_constraint_complete_neck.py --check
```

It reads the recorded tensor and source-owner diagnostics.  It does not rerun
the angular, compact, gauge, MMP, metric or historical generators.
