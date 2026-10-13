# Bronnikov child-frame stress ledger

This retained calculation projects a supplied angularly homogeneous neck
tensor onto the fixed expanding Bronnikov child. It supplies no state,
metric equation solution, or tensor extension away from that slice.
Use dimensionless throat units and finite $\rho\le0$ with

$$
A=1+3\rho+3(1+\rho^2)(\arctan\rho-\pi/2),\quad
a_\parallel=\sqrt{-A},\quad r=\sqrt{1+\rho^2},\quad
dT=-d\rho/a_\parallel.
$$

The proper-volume element per coordinate longitudinal length is
$V=4\pi a_\parallel r^2$. Set
$H_\parallel=A'/(2a_\parallel)$,
$H_\perp=-a_\parallel\rho/(1+\rho^2)$ and
$\theta=H_\parallel+2H_\perp$, so $\dot V=V\theta$.
For orthonormal tensor components $(\varepsilon,T^{01},p_\parallel,p_\perp)$
and declared proper-volume transfer densities $Q,P$, the homogeneous ledger is

$$
\dot\varepsilon=-\theta\varepsilon-H_\parallel p_\parallel
-2H_\perp p_\perp+Q,\qquad
\dot T^{01}=-2(H_\parallel+H_\perp)T^{01}+P.
$$

Thus $\dot{(V\varepsilon)}=-V(H_\parallel p_\parallel+2H_\perp p_\perp-Q)$.
The spatial Killing charge is $a_\parallel VT^{01}$, with rate
$a_\parallel VP$; oriented parent Killing power is its negative.
The locked free block-diagonal Gaussian case has $Q=P=0$. Conservation
residuals check this local projection, while agreement with a stored parent
power remains a separate input check. The projection does not determine a
physical source, a cosmological matter fraction or a coupled child evolution.
