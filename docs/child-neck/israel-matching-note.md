# Exploratory Israel matching of the closed child to the parent's Schwarzschild exterior

7 Oct 2026. Script `israel_matching.py`, data `israel_matching.json`, figure `israel_matching.png`. Needs the restored `child_*` scripts and `junction_neck.py`.
Units G = c = 1 in formulas; the example is a 10 solar-mass parent (r_g = 29.5 km, G M_p/c^3 = 49.3 microseconds) and the child with N_need = 1.82e78 fermions (bare mass 7.9e40 solar masses, a_b = 7.1e-10 m, radiation radius r_rad = 4.08e17 m).
The torsion contact term is treated as an effective fluid in the metric (the Einstein equations hold with the effective stress), so the junction conditions are the standard GR ones. Spin delta-function terms at the neck are not included.

## 0. Correction to the earlier neck note (own error)

The earlier note gave the neck at the bounce as `R_b = (3 M_p/(4 pi eps_b))^(1/3) = 3.57e-23 m`, independent of N. That used the bare density. The Misner-Sharp mass is geometric, so it contains the negative torsion piece:
**10 October correction:** `a_b` below is the density-cancellation scale,
where `eps = alpha*n^2`, used as a high-density approximation in
[child_expansion_ladder.py](../../scripts/child_neck/child_expansion_ladder.py).
It is not the exact real turning point of the closed model. The full equation
is `adot^2 = (kappa/3)(eps-alpha*n^2)a^2 - 1`. At a real bounce `adot=0`,
so `1+adot^2=1`, not zero. Exact density cancellation instead gives
`adot^2=-1`, outside the real closed-FLRW trajectory. The former claim of
zero geometric child mass at a real bounce is withdrawn. In the stated
mass identity its value at an actual turning point is
`m = a*sin(chi)^3/2` in geometric units.

The following near-cancellation neck onset and its tiny `x^2-1` estimates
remain historical approximation outputs, not established exact-bounce
predictions. Retaining a small curvature correction may be essential in that
limit. A consistent near-turning-point treatment is missing; no replacement
trajectory or corrected numerical result is claimed here. The standard angular
junction identities in section 1 can still be reused in their own domains.

With the geometric mass `R_j^3 = a_b^3 (M_p/M_b) x^6/(x^2-1)`, `x = a/a_b`. A neck with mass M_p first exists at `x^2 - 1 = M_p/M_b = 1.3e-40`, where it sits at the equator (`R_j = a_b`). Its smallest value is `R_min = 1.89 a_b (M_p/M_b)^(1/3) = 6.7e-23 m` at `a = 1.22 a_b` (the old number times 1.89).
The larger consequence: this mass-matched surface is spacelike (`|a chi'| >> 1`) everywhere before `t_exit/27` (`|a chi'| = 6.8e45` at `x^2 - 1 = 1e-39`, `6.1e15` at `1e-3`, `6.8e12` at `x^2 = 2`). So the earlier "neck grows as t^(2/3) from the bounce" holds only as a timelike surface from `t_exit/27` on.
Unchanged and rechecked: the exit at `t = G M_p/c^3` (49 microseconds), `a_exit = 1.10e11 m`, and the hidden-child limits `N_max = 3.5e58` (radiation) and `9.7e57` (all matter).

## 1. Formulas (checked)

Surface `chi = chi(t)` in the closed FRW child, areal radius `R = a sin(chi)`, proper time `tau`. On both sides `K^theta_theta = beta/R` with `beta = n . grad R` and `beta^2 = 1 - 2m/R + Rdot^2`. Child side `2m = a sin^3(chi)(1 + adot^2)`; exterior `m = M_p`.
Child side: `K^tau_tau = (a chi'' + 2 adot chi' - a^2 adot chi'^3)/(1 - a^2 chi'^2)^(3/2)`. Exterior: `K^tau_tau = (Rddot + M_p/R^2)/beta_+`.
Shell content, normal pointing from the child to the exterior, `[X] = X_+ - X_-`: `sigma = -[beta]/(4 pi R)`, `P = ([K^tau_tau] + [K^theta_theta])/(8 pi)`.
Checks that pass: `beta^2 - Rdot^2 = 1 - 2m/R` symbolically for arbitrary `a(t)`, `chi(t)`; the Codazzi (conservation) identity `d sigma/d tau + (sigma + P) Adot/A = -(eps + p) a chi' gamma^2` on a radiation child with an arbitrary timelike neck and arbitrary `M`, equal to 25 digits (this fixes the sign conventions and tests both `K^tau_tau` expressions); Oppenheimer-Snyder (comoving cap, dust) gives `sigma = P = 0` to 1e-16 and 1e-14.

## 2. When no shell is needed

`sigma = P = 0` means `K_ab` is continuous. The Gauss and Codazzi equations then give `[T_nn] = [T_n tau] = 0`. With a vacuum exterior and a perfect-fluid child this requires the neck to be comoving with the child fluid and `p = 0` at the neck (Oppenheimer-Snyder).
The shell-free case is the cap child (the child is only the ball inside the neck, `beta_- > 0`) filled with dust: it holds at most M_p, and at turnaround it is the bag-of-gold equality `M_turn = a_max/2 = M_p` (`sigma = 0` there; `P` was not checked).
A closed child whose bare mass exceeds M_p, glued to the exterior through a small neck (`beta_- < 0`, `beta_+ > 0`), has `sigma = -(|beta_+| + |beta_-|)/(4 pi R) < 0`. The negative layer is the Einstein-Rosen flare-out condition, and no choice of neck path removes it (except at an instant of time symmetry with `beta = 0` on both sides).
Since `beta_-^2 - beta_+^2 = 2(M_p - m_-)/R` (`m_-` the child's mass inside the neck sphere), every path obeys `|E| >= sqrt(2R|M_p - m_-|)` with `E = 4 pi R^2 sigma`, with equality when `beta_+ = 0`.

## 3. The two neck paths for the radiation child (N_need, 10 solar masses)

**Mass-matched path** (`m_- = M_p`). Timelike only from `a = 3^(-3/2) a_exit = 0.192 a_exit`, `t = t_exit/27 = 1.82 microseconds`, `R = r_g/9` (analytic; the numbers agree to 1e-14). For an all-matter child the onset is at `a = 0.0836 a_exit`, `t = 0.0121 t_exit`, `R = 0.0582 r_g`. Analytically independent of N and M_p for a radiation child; numerically unchanged for eta up to 5.8e-6.
`E = 4 pi R^2 sigma` in units of M_p, and `w = P/sigma`:

| a / a_exit | R / r_g | E / M_p | w |
|---|---|---|---|
| 0.30 | 0.201 | -3.20 | -0.112 |
| 0.50 | 0.397 | -3.44 | -0.346 |
| 1.0 (exit) | 1.00 | -5.66 | -0.4375 |
| 3.0 | 4.33 | -18.9 | -0.484 |
| 10 | 21.5 | -87.7 | -0.497 |

At exit `|beta| = sqrt 2`, `E = -4 sqrt 2 M_p`. The cap topology on the same path has `sigma = 0` and a surface tension: `P A = -0.53 M_p` at exit, tending to `-M_p/2`.

**Thread path** (`beta_+ = 0`: the exterior side follows the radial geodesic of the white-hole region, `R = (r_g/2)(1 - cos eta)`, `tau = (r_g/2)(eta - sin eta)`, `Rdot = cot(eta/2)`). This is the parent's sphere radius read in reverse, so the mirror in the first note is exact for this path. The child time `t(eta)` follows from the timelike condition (`|a chi'| < 1` throughout; within 2e-7 of null for `R < 1e-7 r_g`).
`E = -sqrt(2R(M_p - m_-))` (checked to 5e-11), the smallest demand of any path; `w = -1/2` for `R << r_g`.

| R / r_g | E / M_p | m_- / M_p |
|---|---|---|
| 2.9e-15 | -1.1e-7 | 2.9e-15 |
| 8.6e-8 | -5.9e-4 | 8.6e-8 |
| 4.6e-4 | -0.043 | 4.6e-4 |
| 0.386 | -1.10 | 0.22 |
| 1 | -1.91 | 0.090 |

The path ends at `R = r_g` with `Rdot = 0` (the bifurcation sphere, the throat of the Einstein-Rosen bridge) after parent proper time `pi G M_p/c^3` (155 microseconds); the child clock reads `3.34 G M_p/c^3` (164 microseconds).

## 4. Sizing the negative layer

Energy held by the neck volume at the bounce density `eps_b = 20.2 eps_P` against the thread-path demand `sqrt(R r_g)`: the volume holds enough for `R >= 1.9e-28 m` (`1.2e7` Planck lengths). As a thin layer at `eps_b`, the needed thickness is `1.8e-70 m` at `R = 1 m` and `3.3e-77 m` at `R = r_g`.
The effective density `eps - alpha n^2` is zero at `n_b` and negative above it (`-0.53 eps_b` at `1.5 n_b`, `-1.48 eps_b` at `2 n_b`, `-16 eps_b` at `5 n_b`), so `eps_b` is the right order for the scale of a layer at `n > n_b`. That the child's matter can supply such a layer at the neck is not shown here.

## 5. Not accounted for

- What the layer is made of (torsion above `n_b`, or something else), and its equation of state (`w = -1/2` on the thread path).
- When the thread meets the child (here the neck starts at `R = 0` at the child's bounce).
- Spin delta-function terms in the junction at the neck.
- The child's recollapse and a second neck; the bag-of-gold case with `P` checked.
